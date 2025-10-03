"""
Improved N-glycan builder with strict biosynthetic rules.

This module provides an improved build_glycan function that follows
strict N-glycan biosynthesis rules for layer-based construction.
"""

import networkx as nx
from itertools import combinations, product


def build_nglycan_strict(hexnac_count, hex_count, fuc_count, sia_count):
    """
    Build N-glycan structure following strict biosynthetic layering rules.

    Rules:
    - Layer 0: 1 HexNAc (reducing end)
    - Layer 1: 1 HexNAc
    - Layer 2: 1 Hex (core mannose)
    - Layer 3: 2 Hex (branching mannoses)

    High Mannose (HexNAc == 2):
    - Layer 4: Up to 4 Hex, connected to layer 3 Hex
    - Layer 5+: Additional Hex, connected to layer 4 Hex (Hex-Hex chains)

    Complex/Hybrid (HexNAc > 2):
    - Layer 4: Additional HexNAc (max 4), connected to layer 3 Hex
    - Layer 5: Additional Hex, connected to layer 4 HexNAc
    - Layer 6+: Sia, connected to Hex in layer 5 or later

    General:
    - Fuc: Can attach to any HexNAc EXCEPT layer 1 HexNAc
    - Additional HexNAc (>6): Added to Hex in layer 5
    - Extended structures: Hex chains can extend from layer 5+

    Parameters:
        hexnac_count: Total number of HexNAc residues
        hex_count: Total number of Hex residues
        fuc_count: Total number of Fuc residues
        sia_count: Total number of Sia residues

    Returns:
        NetworkX DiGraph with nodes having 'type' and 'layer' attributes

    Raises:
        ValueError: If core structure cannot be built (need minimum 2 HexNAc, 3 Hex)
    """
    # Validate minimum requirements for N-glycan core
    if hexnac_count < 2:
        raise ValueError(f"N-glycan requires minimum 2 HexNAc, got {hexnac_count}")
    if hex_count < 3:
        raise ValueError(f"N-glycan requires minimum 3 Hex, got {hex_count}")

    G = nx.DiGraph()
    node_id = 0

    # Track nodes by layer for easy access
    nodes_by_layer = {}

    def add_node(node_type, layer, parent=None):
        """Helper to add node and track by layer."""
        nonlocal node_id
        node_id += 1
        G.add_node(node_id, type=node_type, label=node_type, layer=layer)
        if parent is not None:
            G.add_edge(parent, node_id)

        if layer not in nodes_by_layer:
            nodes_by_layer[layer] = []
        nodes_by_layer[layer].append(node_id)
        return node_id

    # === CORE STRUCTURE (Layers 0-3) ===
    # Layer 0: Reducing end HexNAc
    n0_hexnac = add_node("HexNAc", 0)

    # Layer 1: Second HexNAc (chitobiose core)
    n1_hexnac = add_node("HexNAc", 1, parent=n0_hexnac)

    # Layer 2: Core mannose
    n2_hex = add_node("Hex", 2, parent=n1_hexnac)

    # Layer 3: Two branching mannoses (3-arm and 6-arm)
    n3_hex_3arm = add_node("Hex", 3, parent=n2_hex)  # 3-arm
    n3_hex_6arm = add_node("Hex", 3, parent=n2_hex)  # 6-arm

    # Track what we've used
    hexnac_used = 2
    hex_used = 3
    fuc_used = 0
    sia_used = 0

    # === LAYER 4: High mannose vs Complex type branching ===
    layer4_hexnac_nodes = []
    layer4_hex_nodes = []
    layer3_hex_nodes = [n3_hex_3arm, n3_hex_6arm]

    # For high mannose (HexNAc == 2), add Hex to layer 4 instead of HexNAc
    if hexnac_count == 2:
        # High mannose: Add up to 4 Hex in layer 4
        for i in range(min(4, hex_count - hex_used)):
            parent_hex = layer3_hex_nodes[i % 2]  # Alternate between 3-arm and 6-arm
            hex_node = add_node("Hex", 4, parent=parent_hex)
            layer4_hex_nodes.append(hex_node)
            hex_used += 1
    else:
        # Complex/Hybrid: Add HexNAc to layer 4 (max 4)
        for i in range(min(4, hexnac_count - hexnac_used)):
            parent_hex = layer3_hex_nodes[i % 2]  # Alternate
            hexnac_node = add_node("HexNAc", 4, parent=parent_hex)
            layer4_hexnac_nodes.append(hexnac_node)
            hexnac_used += 1

    # === LAYER 5+: Additional Hex ===
    layer5_hex_nodes = []

    # For complex glycans: Add Hex connected to layer 4 HexNAc
    for hexnac_node in layer4_hexnac_nodes:
        if hex_used >= hex_count:
            break
        hex_node = add_node("Hex", 5, parent=hexnac_node)
        layer5_hex_nodes.append(hex_node)
        hex_used += 1

    # For high mannose: Add Hex connected to layer 4 Hex
    for hex_parent in layer4_hex_nodes:
        if hex_used >= hex_count:
            break
        hex_node = add_node("Hex", 5, parent=hex_parent)
        layer5_hex_nodes.append(hex_node)
        hex_used += 1

    # Add any remaining Hex to available layer 5+ positions
    # This handles extended structures
    available_hex_parents = list(layer5_hex_nodes)  # Can extend from layer 5 Hex
    while hex_used < hex_count and available_hex_parents:
        parent_hex = available_hex_parents.pop(0)
        parent_layer = G.nodes[parent_hex]['layer']
        new_hex = add_node("Hex", parent_layer + 1, parent=parent_hex)
        hex_used += 1
        # New Hex can also be a parent for more Hex if needed
        available_hex_parents.append(new_hex)

    # If still need more Hex, try adding to layer 3 terminals (if capacity allows)
    if hex_used < hex_count:
        for parent_hex in layer3_hex_nodes:
            if hex_used >= hex_count:
                break
            # Check if this Hex has capacity (max 2 children per node)
            if G.out_degree(parent_hex) < 2:
                new_hex = add_node("Hex", 4, parent=parent_hex)
                hex_used += 1
                layer4_hex_nodes.append(new_hex)
                # This new hex can also be extended
                available_hex_parents.append(new_hex)

    # === ADDITIONAL HexNAc (>6): Add to layer 5 Hex ===
    if hexnac_used < hexnac_count and layer5_hex_nodes:
        for hex_node in layer5_hex_nodes:
            if hexnac_used >= hexnac_count:
                break
            # Add HexNAc to layer 5 Hex
            hexnac_node = add_node("HexNAc", 6, parent=hex_node)
            hexnac_used += 1

            # If this HexNAc needs Hex, add it
            if hex_used < hex_count:
                new_hex = add_node("Hex", 7, parent=hexnac_node)
                hex_used += 1

    # === SIALIC ACID: Add to Hex in layer 5 or later ===
    # Get all Hex nodes in layer 5+
    hex_for_sia = []
    for layer in sorted(nodes_by_layer.keys()):
        if layer >= 5:
            for node in nodes_by_layer[layer]:
                if G.nodes[node]['type'] == 'Hex':
                    hex_for_sia.append(node)

    for hex_node in hex_for_sia:
        if sia_used >= sia_count:
            break
        hex_layer = G.nodes[hex_node]['layer']
        sia_node = add_node("Neu5Ac", hex_layer + 1, parent=hex_node)
        sia_used += 1

    # If not enough Hex in layer 5+, warn
    if sia_used < sia_count:
        print(f"Warning: Could only add {sia_used}/{sia_count} Sia. "
              f"Need more Hex in layer 5+ for Sia attachment.")

    # === FUCOSE: Add to HexNAc (except layer 1) ===
    # Get all HexNAc nodes except layer 1
    hexnac_for_fuc = []
    for node in G.nodes():
        if G.nodes[node]['type'] == 'HexNAc' and G.nodes[node]['layer'] != 1:
            # Check if this HexNAc can accept more children (max 2 outgoing edges)
            if G.out_degree(node) < 2:
                hexnac_for_fuc.append(node)

    for hexnac_node in hexnac_for_fuc:
        if fuc_used >= fuc_count:
            break
        hexnac_layer = G.nodes[hexnac_node]['layer']
        fuc_node = add_node("Fuc", hexnac_layer + 1, parent=hexnac_node)
        fuc_used += 1

    # If not enough capacity, warn
    if fuc_used < fuc_count:
        print(f"Warning: Could only add {fuc_used}/{fuc_count} Fuc. "
              f"Not enough available HexNAc with capacity.")

    # Final validation
    if hexnac_used < hexnac_count:
        print(f"Warning: Used only {hexnac_used}/{hexnac_count} HexNAc")
    if hex_used < hex_count:
        print(f"Warning: Used only {hex_used}/{hex_count} Hex")

    return G


def validate_nglycan_structure(G, expected_counts):
    """
    Validate that an N-glycan structure follows the biosynthetic rules.

    Parameters:
        G: NetworkX DiGraph
        expected_counts: dict with keys 'HexNAc', 'Hex', 'Fuc', 'Neu5Ac'

    Returns:
        tuple: (is_valid, list of error messages)
    """
    errors = []

    # Check that each Hex has maximum 2 children
    for node in G.nodes():
        if G.nodes[node]['type'] == 'Hex':
            if G.out_degree(node) > 2:
                errors.append(f"Hex node {node} has {G.out_degree(node)} children (max 2 allowed)")

    # Check core structure (layers 0-3)
    if 0 not in {G.nodes[n]['layer'] for n in G.nodes() if G.nodes[n]['type'] == 'HexNAc'}:
        errors.append("Missing HexNAc in layer 0 (reducing end)")

    layer1_hexnac = [n for n in G.nodes()
                     if G.nodes[n]['type'] == 'HexNAc' and G.nodes[n]['layer'] == 1]
    if len(layer1_hexnac) != 1:
        errors.append(f"Layer 1 should have exactly 1 HexNAc, found {len(layer1_hexnac)}")

    layer2_hex = [n for n in G.nodes()
                  if G.nodes[n]['type'] == 'Hex' and G.nodes[n]['layer'] == 2]
    if len(layer2_hex) != 1:
        errors.append(f"Layer 2 should have exactly 1 Hex, found {len(layer2_hex)}")

    layer3_hex = [n for n in G.nodes()
                  if G.nodes[n]['type'] == 'Hex' and G.nodes[n]['layer'] == 3]
    if len(layer3_hex) != 2:
        errors.append(f"Layer 3 should have exactly 2 Hex, found {len(layer3_hex)}")

    # Check Fuc is not on layer 1 HexNAc
    for hexnac in layer1_hexnac:
        for child in G.successors(hexnac):
            if G.nodes[child]['type'] == 'Fuc':
                errors.append("Fuc cannot be attached to layer 1 HexNAc")

    # Check Sia is only on Hex in layer 5+
    for node in G.nodes():
        if G.nodes[node]['type'] == 'Neu5Ac':
            parents = list(G.predecessors(node))
            if parents:
                parent = parents[0]
                if G.nodes[parent]['type'] != 'Hex':
                    errors.append(f"Sia (node {node}) attached to {G.nodes[parent]['type']}, must be Hex")
                elif G.nodes[parent]['layer'] < 5:
                    errors.append(f"Sia (node {node}) attached to Hex in layer {G.nodes[parent]['layer']}, must be layer 5+")

    # Check counts
    actual_counts = {}
    for node in G.nodes():
        node_type = G.nodes[node]['type']
        actual_counts[node_type] = actual_counts.get(node_type, 0) + 1

    for sugar_type in ['HexNAc', 'Hex', 'Fuc', 'Neu5Ac']:
        expected = expected_counts.get(sugar_type, 0)
        actual = actual_counts.get(sugar_type, 0)
        if actual != expected:
            errors.append(f"{sugar_type}: expected {expected}, found {actual}")

    return (len(errors) == 0, errors)


def generate_all_isomers(hexnac_count, hex_count, fuc_count, sia_count):
    """
    Generate all possible structural isomers for a given glycan composition.

    This function explores all combinatorial possibilities for:
    - Antenna placement (which layer 3 Hex gets which antennae)
    - Fucose attachment sites (core vs. antenna fucosylation)
    - Sialic acid placement (which antenna termini)
    - High mannose branching patterns

    Parameters:
        hexnac_count: Total number of HexNAc residues
        hex_count: Total number of Hex residues
        fuc_count: Total number of Fuc residues
        sia_count: Total number of Sia residues

    Returns:
        list of NetworkX DiGraph objects, each representing a unique isomer
    """
    isomers = []

    # Validate minimum requirements
    if hexnac_count < 2 or hex_count < 3:
        return isomers

    # Build core structure (same for all isomers)
    def build_core():
        G = nx.DiGraph()
        node_id = [0]  # Use list to allow modification in nested function

        def add_node(node_type, layer, parent=None):
            node_id[0] += 1
            nid = node_id[0]
            G.add_node(nid, type=node_type, label=node_type, layer=layer)
            if parent is not None:
                G.add_edge(parent, nid)
            return nid

        # Core layers 0-3
        n0 = add_node("HexNAc", 0)
        n1 = add_node("HexNAc", 1, parent=n0)
        n2 = add_node("Hex", 2, parent=n1)
        n3_a = add_node("Hex", 3, parent=n2)  # 3-arm
        n3_b = add_node("Hex", 3, parent=n2)  # 6-arm

        return G, add_node, n0, n1, n3_a, n3_b

    # Determine if high mannose or complex type
    is_high_mannose = (hexnac_count == 2)

    if is_high_mannose:
        # High mannose: generate all branching patterns for Hex
        isomers.extend(_generate_high_mannose_isomers(hexnac_count, hex_count, fuc_count, sia_count))
    else:
        # Complex/Hybrid: generate all antenna and fucose combinations
        isomers.extend(_generate_complex_isomers(hexnac_count, hex_count, fuc_count, sia_count))

    return isomers


def _generate_high_mannose_isomers(hexnac_count, hex_count, fuc_count, sia_count):
    """Generate all high mannose isomers (HexNAc=2)."""
    isomers = []

    # For high mannose, main variation is in how Hex distribute between 3-arm and 6-arm
    # Layer 4 can have up to 4 Hex total, distributed between two arms
    remaining_hex = hex_count - 3  # After core

    # Generate all ways to distribute Hex in layer 4 between two arms
    max_layer4 = min(4, remaining_hex)

    for n_layer4 in range(0, max_layer4 + 1):
        # Try all distributions between 3-arm and 6-arm
        for n_3arm in range(min(2, n_layer4) + 1):
            n_6arm = n_layer4 - n_3arm
            if n_6arm > 2:  # Max 2 per arm
                continue

            # Build this isomer
            G = nx.DiGraph()
            node_id = [0]

            def add_node(node_type, layer, parent=None):
                node_id[0] += 1
                nid = node_id[0]
                G.add_node(nid, type=node_type, label=node_type, layer=layer)
                if parent is not None:
                    G.add_edge(parent, nid)
                return nid

            # Core
            n0 = add_node("HexNAc", 0)
            n1 = add_node("HexNAc", 1, parent=n0)
            n2 = add_node("Hex", 2, parent=n1)
            n3_a = add_node("Hex", 3, parent=n2)
            n3_b = add_node("Hex", 3, parent=n2)

            hex_used = 3
            layer4_hex = []

            # Add Hex to 3-arm (max 2 children per Hex)
            if n_3arm > 2:
                continue
            for _ in range(n_3arm):
                h = add_node("Hex", 4, parent=n3_a)
                layer4_hex.append(h)
                hex_used += 1

            # Add Hex to 6-arm (max 2 children per Hex)
            if n_6arm > 2:
                continue
            for _ in range(n_6arm):
                h = add_node("Hex", 4, parent=n3_b)
                layer4_hex.append(h)
                hex_used += 1

            # Extend with layer 5+ Hex (chain extensions)
            # Each Hex can only have max 2 children
            remaining = hex_count - hex_used
            if remaining > 0 and layer4_hex:
                # Build chains from layer 4 Hex, respecting max 2 children per node
                extension_queue = list(layer4_hex)
                layer_offset = 1

                while hex_used < hex_count and extension_queue:
                    parent = extension_queue.pop(0)
                    # Check if parent already has 2 children
                    if G.out_degree(parent) >= 2:
                        continue

                    parent_layer = G.nodes[parent]['layer']
                    h = add_node("Hex", parent_layer + 1, parent=parent)
                    hex_used += 1

                    # This new Hex can also be extended
                    if hex_used < hex_count:
                        extension_queue.append(h)

            # Add fucose if any (to HexNAc except layer 1)
            if fuc_count > 0:
                add_node("Fuc", 1, parent=n0)  # Core fucose

            if hex_used == hex_count:
                isomers.append(G.copy())

    return isomers


def _generate_complex_isomers(hexnac_count, hex_count, fuc_count, sia_count):
    """Generate all complex/hybrid type isomers."""
    isomers = []

    # Calculate antenna count (HexNAc in layer 4)
    n_antennae = hexnac_count - 2  # Minus core

    # Max direct antennae on layer 3 Hex is 4 (2 per arm)
    # Additional HexNAc go to layer 6 as extended antennae (HexNAc-Hex-HexNAc)
    max_direct_antennae = min(4, n_antennae)

    # For extended structures (>4 HexNAc), we need enough Hex
    # Each direct antenna needs 1 Hex, extended antennae don't need extra Hex
    # So minimum Hex needed = 3 (core) + number of direct antennae

    # Check if this is a hybrid glycan
    # Hybrid: more Hex than expected for pure complex type with max direct antennae
    hex_needed_for_max_complex = 3 + max_direct_antennae
    is_hybrid = hex_count > hex_needed_for_max_complex

    # Generate all antenna distributions
    for n_ant in range(1, max_direct_antennae + 1):
        # n_ant is the number of direct antennae in layer 4
        # Additional HexNAc will be added as extended antennae

        # Distribute between 3-arm and 6-arm
        for n_3arm in range(min(2, n_ant) + 1):
            n_6arm = n_ant - n_3arm
            if n_6arm > 2:  # Max 2 per arm
                continue
            if n_3arm == 0 and n_6arm == 0:
                continue

            # For hybrid glycans, also try different mannose arm configurations
            if is_hybrid:
                # Generate hybrid variants with mannose arm on either 3-arm or 6-arm
                for hybrid_config in _generate_hybrid_configs(n_3arm, n_6arm, hex_count, n_antennae):
                    for fuc_variant in _generate_fucose_variants(n_ant, fuc_count):
                        G = _build_hybrid_glycan(hexnac_count, hex_count, fuc_count, sia_count,
                                                hybrid_config, fuc_variant)
                        if G is not None:
                            isomers.append(G)
            else:
                # Pure complex type
                for fuc_variant in _generate_fucose_variants(n_ant, fuc_count):
                    G = _build_complex_glycan(hexnac_count, hex_count, fuc_count, sia_count,
                                             n_3arm, n_6arm, fuc_variant)
                    if G is not None:
                        isomers.append(G)

    return isomers


def _generate_hybrid_configs(n_3arm_ant, n_6arm_ant, hex_count, n_antennae):
    """
    Generate hybrid glycan configurations.

    Hybrid glycans have one arm with HexNAc antennae and one with mannose (Hex) extensions.
    Returns configs: (n_3arm_hexnac, n_6arm_hexnac, mannose_arm, n_mannose)
    where mannose_arm is '3' or '6' indicating which arm has mannose extension.
    """
    configs = []

    # Calculate extra Hex (beyond core + terminal Hex for antennae)
    core_hex = 3
    terminal_hex = n_antennae  # One Hex per antenna
    extra_hex = hex_count - core_hex - terminal_hex

    if extra_hex <= 0:
        return configs

    # Try mannose arm on 3-arm (6-arm has antennae)
    if n_6arm_ant > 0:
        configs.append({
            'n_3arm_hexnac': 0,
            'n_6arm_hexnac': n_6arm_ant,
            'mannose_arm': '3',
            'n_mannose': extra_hex
        })

    # Try mannose arm on 6-arm (3-arm has antennae)
    if n_3arm_ant > 0:
        configs.append({
            'n_3arm_hexnac': n_3arm_ant,
            'n_6arm_hexnac': 0,
            'mannose_arm': '6',
            'n_mannose': extra_hex
        })

    return configs


def _generate_fucose_variants(n_antennae, fuc_count):
    """Generate all fucosylation patterns (core vs. antenna)."""
    variants = []

    # Fucose can attach to: layer 0 HexNAc (core) or antenna HexNAc (layer 4)
    # Total positions: 1 (core) + n_antennae
    total_positions = 1 + n_antennae

    if fuc_count == 0:
        return [set()]

    # Generate all combinations
    for fuc_positions in combinations(range(total_positions), min(fuc_count, total_positions)):
        variants.append(set(fuc_positions))

    return variants


def _build_hybrid_glycan(hexnac_count, hex_count, fuc_count, sia_count,
                        hybrid_config, fuc_variant):
    """Build a hybrid glycan isomer (one mannose arm, one complex arm)."""
    G = nx.DiGraph()
    node_id = [0]

    def add_node(node_type, layer, parent=None):
        node_id[0] += 1
        nid = node_id[0]
        G.add_node(nid, type=node_type, label=node_type, layer=layer)
        if parent is not None:
            G.add_edge(parent, nid)
        return nid

    # Core
    n0 = add_node("HexNAc", 0)
    n1 = add_node("HexNAc", 1, parent=n0)
    n2 = add_node("Hex", 2, parent=n1)
    n3_a = add_node("Hex", 3, parent=n2)
    n3_b = add_node("Hex", 3, parent=n2)

    hexnac_used = 2
    hex_used = 3
    fuc_used = 0
    sia_used = 0

    # Add core fucose if specified
    if 0 in fuc_variant:
        add_node("Fuc", 1, parent=n0)
        fuc_used += 1

    # Determine which arm gets mannose and which gets antennae
    mannose_arm = hybrid_config['mannose_arm']
    n_mannose = hybrid_config['n_mannose']
    n_3arm_hexnac = hybrid_config['n_3arm_hexnac']
    n_6arm_hexnac = hybrid_config['n_6arm_hexnac']

    # Build mannose arm (Hex extensions)
    mannose_parent = n3_a if mannose_arm == '3' else n3_b
    mannose_nodes = []

    # Add Hex to mannose arm (max 2 direct children)
    for i in range(min(2, n_mannose)):
        h = add_node("Hex", 4, parent=mannose_parent)
        mannose_nodes.append(h)
        hex_used += 1

    # Extend with more Hex if needed (chain extensions)
    remaining_mannose = n_mannose - len(mannose_nodes)
    if remaining_mannose > 0 and mannose_nodes:
        extension_queue = list(mannose_nodes)
        while remaining_mannose > 0 and extension_queue:
            parent = extension_queue.pop(0)
            if G.out_degree(parent) >= 2:
                continue
            parent_layer = G.nodes[parent]['layer']
            h = add_node("Hex", parent_layer + 1, parent=parent)
            hex_used += 1
            remaining_mannose -= 1
            extension_queue.append(h)

    # Build antenna arm (HexNAc + Hex + Sia)
    antenna_hexnac = []
    antenna_idx = 1  # 0 is core

    complex_parent = n3_b if mannose_arm == '3' else n3_a

    # Add antennae to complex arm
    n_ant = n_3arm_hexnac if mannose_arm == '6' else n_6arm_hexnac
    for _ in range(n_ant):
        hn = add_node("HexNAc", 4, parent=complex_parent)
        antenna_hexnac.append(hn)
        hexnac_used += 1

        # Add fucose if specified
        if antenna_idx in fuc_variant:
            add_node("Fuc", 5, parent=hn)
            fuc_used += 1
        antenna_idx += 1

    # Add terminal Hex to antennae
    antenna_hex = []
    for hn in antenna_hexnac:
        if hex_used >= hex_count:
            break
        h = add_node("Hex", 5, parent=hn)
        antenna_hex.append(h)
        hex_used += 1

    # Add sialic acids to terminal Hex
    for h in antenna_hex[:sia_count]:
        if sia_used >= sia_count:
            break
        add_node("Neu5Ac", 6, parent=h)
        sia_used += 1

    # Validate counts
    if hexnac_used == hexnac_count and hex_used == hex_count and fuc_used == fuc_count and sia_used == sia_count:
        return G
    return None


def _build_complex_glycan(hexnac_count, hex_count, fuc_count, sia_count,
                         n_3arm, n_6arm, fuc_variant):
    """Build a specific complex glycan isomer."""
    G = nx.DiGraph()
    node_id = [0]

    def add_node(node_type, layer, parent=None):
        node_id[0] += 1
        nid = node_id[0]
        G.add_node(nid, type=node_type, label=node_type, layer=layer)
        if parent is not None:
            G.add_edge(parent, nid)
        return nid

    # Core
    n0 = add_node("HexNAc", 0)
    n1 = add_node("HexNAc", 1, parent=n0)
    n2 = add_node("Hex", 2, parent=n1)
    n3_a = add_node("Hex", 3, parent=n2)
    n3_b = add_node("Hex", 3, parent=n2)

    hexnac_used = 2
    hex_used = 3
    fuc_used = 0
    sia_used = 0

    # Add core fucose if specified
    if 0 in fuc_variant:
        add_node("Fuc", 1, parent=n0)
        fuc_used += 1

    # Add antennae (max 2 per arm)
    if n_3arm > 2 or n_6arm > 2:
        return None

    antenna_hexnac = []
    antenna_idx = 1  # Start from 1 (0 is core)

    for _ in range(n_3arm):
        hn = add_node("HexNAc", 4, parent=n3_a)
        antenna_hexnac.append(hn)
        hexnac_used += 1

        # Add fucose if specified
        if antenna_idx in fuc_variant:
            add_node("Fuc", 5, parent=hn)
            fuc_used += 1
        antenna_idx += 1

    for _ in range(n_6arm):
        hn = add_node("HexNAc", 4, parent=n3_b)
        antenna_hexnac.append(hn)
        hexnac_used += 1

        # Add fucose if specified
        if antenna_idx in fuc_variant:
            add_node("Fuc", 5, parent=hn)
            fuc_used += 1
        antenna_idx += 1

    # Add terminal Hex to antennae
    antenna_hex = []
    for hn in antenna_hexnac:
        if hex_used >= hex_count:
            break
        h = add_node("Hex", 5, parent=hn)
        antenna_hex.append(h)
        hex_used += 1

    # Add extended antennae (HexNAc-Hex-HexNAc for LacdiNAc structures)
    # This handles cases where we have extra HexNAc beyond the 4 direct antennae
    if hexnac_used < hexnac_count and antenna_hex:
        for h in antenna_hex:
            if hexnac_used >= hexnac_count:
                break
            # Add HexNAc to terminal Hex
            hn = add_node("HexNAc", 6, parent=h)
            hexnac_used += 1

    # Add sialic acids to terminal Hex
    for h in antenna_hex[:sia_count]:
        if sia_used >= sia_count:
            break
        add_node("Neu5Ac", 6, parent=h)
        sia_used += 1

    # Check if we used the right amounts
    if hexnac_used == hexnac_count and hex_used == hex_count and fuc_used == fuc_count and sia_used == sia_count:
        return G
    return None


if __name__ == "__main__":
    # Test cases
    print("Testing N-glycan builder with strict rules\n")

    # Test 1: Basic complex-type glycan
    print("Test 1: Complex-type N-glycan (HexNAc=4, Hex=5, Fuc=1, Sia=2)")
    G1 = build_nglycan_strict(hexnac_count=4, hex_count=5, fuc_count=1, sia_count=2)
    valid, errors = validate_nglycan_structure(G1, {'HexNAc': 4, 'Hex': 5, 'Fuc': 1, 'Neu5Ac': 2})
    print(f"Valid: {valid}")
    if errors:
        for err in errors:
            print(f"  - {err}")
    print(f"Nodes: {G1.number_of_nodes()}, Edges: {G1.number_of_edges()}\n")

    # Test 2: High-mannose type
    print("Test 2: High-mannose N-glycan (HexNAc=2, Hex=9, Fuc=0, Sia=0)")
    G2 = build_nglycan_strict(hexnac_count=2, hex_count=9, fuc_count=0, sia_count=0)
    valid, errors = validate_nglycan_structure(G2, {'HexNAc': 2, 'Hex': 9, 'Fuc': 0, 'Neu5Ac': 0})
    print(f"Valid: {valid}")
    if errors:
        for err in errors:
            print(f"  - {err}")
    print(f"Nodes: {G2.number_of_nodes()}, Edges: {G2.number_of_edges()}\n")

    # Test 3: Highly branched
    print("Test 3: Highly branched N-glycan (HexNAc=6, Hex=7, Fuc=2, Sia=4)")
    G3 = build_nglycan_strict(hexnac_count=6, hex_count=7, fuc_count=2, sia_count=4)
    valid, errors = validate_nglycan_structure(G3, {'HexNAc': 6, 'Hex': 7, 'Fuc': 2, 'Neu5Ac': 4})
    print(f"Valid: {valid}")
    if errors:
        for err in errors:
            print(f"  - {err}")
    print(f"Nodes: {G3.number_of_nodes()}, Edges: {G3.number_of_edges()}\n")
