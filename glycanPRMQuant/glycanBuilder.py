import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.path import Path
import itertools
import pandas as pd

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

def build_glycan(composition, glycan_type="N"):
    """
    Build a glycan structure following custom rules:
    1. Extra HexNAc are added to the core first and all placed at layer 4
    2. Extra Hex are added to the core second
    3. Extra Fuc are added to the core third
    4. Extra Neu5Ac are added to the core fourth
    5. If there is Neu5Ac, there must be a Hex attached after HexNAc and the Neu5Ac must be on the Hex
       forming a HexNAc-Hex-Neu5Ac linkage
    6. If there is a Fuc, it must be on the HexNAc
    7. For high mannose structures (HexNAc=2, Hex>=5), one branch remains linear while the other branches
    8. Each hexose can only branch twice maximum
    9. Only one of the terminal mannose on the core can have 2 outgoing connections
    10. Any mannose at or after layer 4 can only have 1 outgoing connection
    
    Parameters:
        composition: A tuple of (total_hexnac, total_hex, total_fuc, total_neu5ac)
        glycan_type: The type of glycan to build. Default is "N" for N-glycan.
    
    Returns:
        A networkx graph representing the glycan structure
    """
    total_hexnac, total_hex, total_fuc, total_neu5ac = composition
    
    # Define the core structure
    core_hexnac, core_hex = 2, 3
    extra_hexnac = max(0, total_hexnac - core_hexnac)
    extra_hex = max(0, total_hex - core_hex)
    
    # Check if this is a high mannose structure
    is_high_mannose = (total_hexnac == 2 and total_hex >= 5)
    
    # Build Core Scaffold
    G = nx.DiGraph()
    
    # Add core HexNAc (chitobiose core)
    G.add_node(1, type="HexNAc", label="HexNAc", layer=0)
    G.add_node(2, type="HexNAc", label="HexNAc", layer=1)
    G.add_edge(1, 2)
    
    # Add core Hex (mannose branches)
    G.add_node(3, type="Hex", label="Hex", layer=2)  # Central mannose
    G.add_edge(2, 3)
    G.add_node(4, type="Hex", label="Hex", layer=3)  # Left mannose branch
    G.add_node(5, type="Hex", label="Hex", layer=3)  # Right mannose branch
    G.add_edge(3, 4)
    G.add_edge(3, 5)
    
    current_node_id = 5
    
    # Helper function to get maximum allowed connections for a node
    def get_max_connections(graph, node):
        """Get the maximum number of outgoing connections allowed for a node"""
        node_type = graph.nodes[node]["type"]
        node_layer = graph.nodes[node]["layer"]
        
        if (node_type == "Hex") & is_high_mannose:
            # For terminal mannose in core (nodes 4 and 5)
            if node == 4:  # Left terminal mannose - can have only 1 connection
                return 1
            elif node == 5:  # Right terminal mannose - can have up to 2 connections
                return 2
            elif node_layer >= 4:  # Any mannose at or after layer 4 can only have 1 connection
                return 1
            else:  # Other mannose can have up to 2 connections
                return 2
        elif node_type == "Hex":
            # For non-terminal mannose (nodes 3, 4, 5)
            if node == 4:
                return 2
            elif node == 5:
                return 2
            elif node_layer >= 4:
                return 1
            else:
                return 2
        elif (node_type == "HexNAc") & (node_layer == 1):
            return 1  # Layer 1 HexNAc can only have 1 connection
        elif node_type == "HexNAc":
            return 2
        else:
            return 0  # Terminal residues (Fuc, Neu5Ac) cannot have connections
    
    # Helper function to get available branch points
    def get_available_branch_points(graph, sugar_type=None, min_layer=None, max_layer=None):
        """
        The function `get_available_branch_points` filters and returns nodes in a graph based on
        specified criteria such as sugar type, minimum and maximum layers, and available connections.
        
        :param graph: The `graph` parameter in the `get_available_branch_points` function is expected to
        be a graph data structure, such as a networkx graph object. This graph represents a network or a
        system where nodes are connected by edges, and each node may have attributes like "type" and
        "layer"
        :param sugar_type: The `sugar_type` parameter in the `get_available_branch_points` function is
        used to filter nodes based on their type. If a `sugar_type` is specified, only nodes with that
        specific type will be considered as available branch points. If `sugar_type` is `None`,
        :param min_layer: The `min_layer` parameter in the `get_available_branch_points` function is
        used to filter nodes based on their layer in the graph. If a `min_layer` value is provided, only
        nodes with a layer greater than or equal to the specified `min_layer` value will be considered
        as available
        :param max_layer: The `max_layer` parameter in the `get_available_branch_points` function is
        used to filter nodes based on their layer in the graph. If specified, only nodes with a layer
        less than or equal to the `max_layer` value will be considered as available branch points. Nodes
        with a layer greater
        :return: The function `get_available_branch_points` returns a list of nodes in a graph that can
        accept more branches based on certain criteria such as sugar type, minimum layer, and maximum
        layer.
        """
        available_nodes = []
        for node in graph.nodes():
            # Skip the check if node doesn't exist anymore (safeguard)
            if node not in graph:
                continue
                
            # Get current out-degree and maximum allowed connections
            out_degree = graph.out_degree(node)
            max_connections = get_max_connections(graph, node)
            
            # Check if node can accept more branches
            if out_degree < max_connections:
                node_type = graph.nodes[node]["type"]
                node_layer = graph.nodes[node]["layer"]
                
                # Apply type filter if specified
                if sugar_type and node_type != sugar_type:
                    continue
                
                # Apply layer filters if specified
                if min_layer is not None and node_layer < min_layer:
                    continue
                if max_layer is not None and node_layer > max_layer:
                    continue
                
                available_nodes.append(node)
        
        return available_nodes
    
    # For high mannose structures:
    if is_high_mannose:
        # Left branch (node 4) remains linear
        # Right branch (node 5) will branch further
        
        # We've used 3 Hex already (nodes 3, 4, 5)
        remaining_hex = extra_hex
        
        # First extend the linear branch (left branch, node 4)
        current_parent = 4
        for i in range(min(2, remaining_hex)):  # Add up to 2 Hex in linear branch
            current_node_id += 1
            G.add_node(current_node_id, type="Hex", label="Hex", layer=G.nodes[current_parent]["layer"] + 1)
            G.add_edge(current_parent, current_node_id)
            current_parent = current_node_id
            remaining_hex -= 1
        
        # Then build the branching side (right branch, node 5)
        if remaining_hex > 0:
            # First add one Hex in a linear fashion
            current_node_id += 1
            next_branch_point = current_node_id
            G.add_node(next_branch_point, type="Hex", label="Hex", layer=G.nodes[5]["layer"] + 1)
            G.add_edge(5, next_branch_point)
            remaining_hex -= 1
            
            # Keep track of available branch points (only for Hex)
            available_branch_points = get_available_branch_points(G, "Hex")
            
            # Add remaining Hex while respecting the branching constraints
            while remaining_hex > 0 and available_branch_points:
                parent = available_branch_points.pop(0)  # Get the first available branch point
                parent_layer = G.nodes[parent]["layer"]
                
                # Add a new Hex
                current_node_id += 1
                new_node_layer = parent_layer + 1
                G.add_node(current_node_id, type="Hex", label="Hex", layer=new_node_layer)
                G.add_edge(parent, current_node_id)
                remaining_hex -= 1
                
                # Update available branch points
                available_branch_points = get_available_branch_points(G, "Hex")
                
                # If we can't add all Hex, that's okay - we respect the branching constraint
                if not available_branch_points and remaining_hex > 0:
                    break
    else:
        # Standard structure processing
        # Define potential branch points for extra HexNAc (terminal mannose only)
        # Only terminal mannose residues can be branch points for HexNAc
        branch_nodes = get_available_branch_points(G, "Hex", min_layer=3, max_layer=3)
        
        # 1. Attach extra HexNAc - ALL at layer 4
        hexnac_nodes = []
        core_terminals = [4, 5]
        terminal_rotation = itertools.cycle(core_terminals)
        for _ in range(extra_hexnac):
            # Find the next eligible parent
            parent = None
            for _ in range(len(core_terminals)):
                candidate = next(terminal_rotation)
                if candidate in get_available_branch_points(G, "Hex", min_layer=3, max_layer=3):
                    parent = candidate
                    break
            if not parent:
                print("Warning: No available parent for extra HexNAc.")
                break
            current_node_id += 1
            G.add_node(current_node_id, type="HexNAc", label="HexNAc", layer=4)
            G.add_edge(parent, current_node_id)
            hexnac_nodes.append(current_node_id)
        # 2. Attach extra Hex

        required_hex_for_neu5ac = total_neu5ac
        # Rule 5: If Neu5Ac is present, ensure Hex is attached to HexNAc
        if total_neu5ac > 0:
            if extra_hex < required_hex_for_neu5ac:
                print(f'Warning: Increasing Hex Count to {required_hex_for_neu5ac} to accommodate Neu5Ac.')
                extra_hex = required_hex_for_neu5ac
            elif extra_hex > required_hex_for_neu5ac:
                pass
        


        # First attach Hex to extra HexNAc nodes if needed for Neu5Ac
        hex_to_neu5ac = []
        hex_nodes = []
        for i in range(min(extra_hex, len(hexnac_nodes))):
            parent = hexnac_nodes[i]
            current_node_id += 1
            G.add_node(current_node_id, type="Hex", label="Hex", layer=5)
            G.add_edge(parent, current_node_id)
            if i < total_neu5ac:
                hex_to_neu5ac.append(current_node_id)
            hex_nodes.append(current_node_id)
            extra_hex -= 1
        
        # Then attach remaining Hex to other branch points, respecting branching constraint
        if extra_hex > 0:
            # Get all available branch points (respecting custom branching rules)
            branch_points = get_available_branch_points(G)

            # Split branch points into two groups: core terminals (4,5) and others
            core_terminals = [n for n in branch_points if n in {4, 5}]
            other_points = [n for n in branch_points if n not in {4, 5}]

            core_terminals.sort(key = lambda x: (G.out_degree(x), x))
            
            priority_nodes = core_terminals + other_points
            
            while extra_hex > 0 and priority_nodes:
                parent = priority_nodes.pop(0)
                parent_layer = G.nodes[parent]["layer"]
                
                current_node_id += 1
                G.add_node(current_node_id, type="Hex", label="Hex", layer=parent_layer + 1)
                G.add_edge(parent, current_node_id)
                hex_nodes.append(current_node_id)
                extra_hex -= 1
                
                # Update priorities after each addition
                branch_points = get_available_branch_points(G)
                core_terminals = sorted([n for n in branch_points if n in {4, 5}],
                                        key = lambda x: (G.out_degree(x), x))
                other_points = [n for n in branch_points if n not in {4, 5}]
                priority_nodes = core_terminals + other_points
                # If we can't add all Hex, that's okay
                if not priority_nodes and extra_hex > 0:
                    break
    
    # 3. Attach Fuc (to HexNAc nodes only)
    # Rule 6: Fuc must be on HexNAc
    # HexNAc can only have 1 Fuc
    fuc_added = 0
    while fuc_added < total_fuc:
        # Get all HexNAc nodes that can generally take more branches [1, 2]
        # get_available_branch_points considers the max_connections rule (up to 2 for HexNac)
        all_generally_available_hexnac = get_available_branch_points(G, "HexNAc")

        # Filter this list to find HexNAc nodes that DO NOT already have a Fuc child
        suitable_hexnac = []
        for hexnac_node in all_generally_available_hexnac:
            has_fuc_child = False
            # Check outgoing edges (children) of potential parent node
            for _, child_node in G.edges(hexnac_node):
                if G.nodes[child_node]["type"] == "Fuc":
                    has_fuc_child = True
                    break
            # If no FUc child, add to suitable list
            if not has_fuc_child:
                suitable_hexnac.append(hexnac_node)
        
        # If no suitable parents are foudn that can accept a Fuc
        # We cannot add more Fuc
        if not suitable_hexnac:
            # Add a warning if not all requested Fuc could be added
            if fuc_added < total_fuc:
                print(f"Warning: Could not add all requested Fuc. Only {fuc_added} added.")
            break

        # Pick a suitable parent node
        parent = suitable_hexnac.pop(0)
        parent_layer = G.nodes[parent]["layer"]

        current_node_id += 1
        G.add_node(current_node_id, type="Fuc", label="Fuc", layer=parent_layer + 1)

        # Add the edge from the parent to the new Fuc node
        G.add_edge(parent, current_node_id)

        # Increment the count of added Fuc
        fuc_added += 1
    


    # 4. Attach Neu5Ac (to Hex nodes only)
    # Rule 5: Neu5Ac must be on Hex in HexNAc-Hex-Neu5Ac linkage
    if not is_high_mannose:  # Only for non-high-mannose structures
        available_hex = get_available_branch_points(G, "Hex")
        
        # Prioritize hex_to_neu5ac nodes first
        priority_hex = [n for n in hex_to_neu5ac if n in available_hex]
        other_hex = [n for n in available_hex if n not in hex_to_neu5ac]
        prioritized_hex = priority_hex + other_hex
        
        neu5ac_added = 0
        
        # Fixed: Iterate while we have both Neu5Ac to add and available Hex
        while neu5ac_added < total_neu5ac and prioritized_hex:
            parent = prioritized_hex.pop(0)  # Get the next available Hex node
            parent_layer = G.nodes[parent]["layer"]
            current_node_id += 1
            G.add_node(current_node_id, type="Neu5Ac", label="Neu5Ac", layer=parent_layer + 1)
            G.add_edge(parent, current_node_id)
            neu5ac_added += 1
            
            # Update available Hex nodes after each addition
            available_hex = get_available_branch_points(G, "Hex")
            priority_hex = [n for n in hex_to_neu5ac if n in available_hex]
            other_hex = [n for n in available_hex if n not in hex_to_neu5ac]
            prioritized_hex = priority_hex + other_hex
    
    if is_high_mannose & (total_neu5ac > 0):
        print("Warning: Neu5Ac cannot be added to high mannose structures.")
    
    return G

def visualize_glycan_snfg(G, composition):
    """
    Visualize the glycan structure using Symbol Nomenclature for Glycans (SNFG).
    
    Parameters:
        G: NetworkX DiGraph representing the glycan structure
    
    Returns:
        Matplotlib figure object with the visualization
    """
    # Define color scheme for glycan components
    colors = {
        "HexNAc": "blue",
        "Hex": "green",
        "Fuc": "red",
        "Neu5Ac": "purple"
    }
    
    # Create position layout based on layers
    pos = {}
    
    # start from teh root node (layer 0)
    root_node = [n for n in G.nodes() if G.in_degree(n) == 0][0]
    max_layer = max(data["layer"] for _, data in G.nodes(data=True))

    def hierarchical_layout(G, root, base_width = 0.8, vert_gap = 0.4, xcenter = 0.5):
        pos = {}

        def _hierarchical_layout_recursive(node, depth, parent_x, width):
            nonlocal pos
            y = depth * vert_gap
            num_children = len(list(G.successors(node)))

            if num_children == 0:
                x = parent_x
                pos[node] = (x, y)
                return [(x, y)]
            
            else:
                # Position children recursively
                child_width = width * num_children * (0.6 ** depth)
                child_x_start = parent_x - child_width / 2
                x_step = width / num_children
                child_positions = []
                
                for i, child in enumerate(G.successors(node)):
                    child_x = child_x_start + x_step * (i + 0.5)
                    child_pos = _hierarchical_layout_recursive(child, depth + 1, child_x,
                                                               child_width)
                    child_positions.extend(child_pos)

                avg_x = sum(x for x,_ in child_positions) / len(child_positions)
                pos[node] = (avg_x, y)
                return [(avg_x, y)] + child_positions
            
        _hierarchical_layout_recursive(root, 0, xcenter, base_width)
        return pos
    pos = hierarchical_layout(G, root_node)
    
    # Create figure and axis
    fig, ax = plt.subplots(figsize=(12, 10))

    ax.set_aspect('equal', adjustable='box')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    
    # Draw edges first (so they're behind nodes)
    for u, v in G.edges():
        ax.annotate("",
                   xy=pos[v], xycoords='data',
                   xytext=pos[u], textcoords='data',
                   arrowprops=dict(arrowstyle="-|>", color="black", 
                                  shrinkA=12, shrinkB=12,
                                  connectionstyle="arc3,rad=0.0"))
    
    # Node size
    node_size = 0.05
    
    # Draw nodes with appropriate shapes and colors
    for node, data in G.nodes(data=True):
        x, y = pos[node]
        sugar_type = data["type"]
        color = colors[sugar_type]
        
        if sugar_type == "HexNAc":
            # Blue square for HexNAc
            square = plt.Rectangle((x - node_size/2, y - node_size/2), 
                                  node_size, node_size, 
                                  facecolor=color, edgecolor='black')
            ax.add_patch(square)
        
        elif sugar_type == "Hex":
            # Green circle for Hex
            circle = plt.Circle((x, y), node_size/2, 
                               facecolor=color, edgecolor='black')
            ax.add_patch(circle)
        
        elif sugar_type == "Fuc":
            # Red triangle for Fucose
            triangle_vertices = [
                (x, y + node_size/2),
                (x - node_size/2, y - node_size/2),
                (x + node_size/2, y - node_size/2)
            ]
            triangle = plt.Polygon(triangle_vertices, 
                                  facecolor=color, edgecolor='black')
            ax.add_patch(triangle)
        
        elif sugar_type == "Neu5Ac":
            # Purple diamond for Sialic acid
            diamond_vertices = [
                (x, y + node_size/2),
                (x + node_size/2, y),
                (x, y - node_size/2),
                (x - node_size/2, y)
            ]
            diamond = plt.Polygon(diamond_vertices, 
                                 facecolor=color, edgecolor='black')
            ax.add_patch(diamond)
        
        # Add node ID for debugging (uncomment if needed)
        # ax.text(x, y, str(node), fontsize=8, ha='center', va='center', color='white')
    
    # Create legend
    # legend_elements = [
    #     mpatches.Rectangle((0,0), 1, 1, facecolor=colors["HexNAc"], edgecolor='black', label='HexNAc'),
    #     mpatches.Circle((0.5,0.5), 0.5, facecolor=colors["Hex"], edgecolor='black', label='Hex'),
    #     mpatches.Polygon([(0.5,1), (0,0), (1,0)], facecolor=colors["Fuc"], edgecolor='black', label='Fucose'),
    #     mpatches.Polygon([(0.5,1), (1,0.5), (0.5,0), (0,0.5)], facecolor=colors["Neu5Ac"], edgecolor='black', label='Sialic acid')
    # ]
    # ax.legend(handles=legend_elements, loc='upper right')
    
    # Set axis limits with some padding
    ax.set_xlim(-0.1, 1.1)
    ax.set_ylim(-0.1, (max_layer + 0.1) * 0.5)
    
    # Show layer numbers on y-axis
    ax.set_yticks([])
    
    # Set title to glycan name
    ax.set_title(f"{get_glycan_name(composition)}", fontsize=16)
    
    # Remove x ticks
    ax.set_xticks([])
    
    plt.tight_layout()
    return fig

def get_glycan_name(composition):
    """Generate a standardized name for the glycan based on composition"""
    total_hexnac, total_hex, total_fuc, total_neu5ac = composition
    
    # Check if this is a high mannose structure
    if total_hexnac == 2 and total_hex >= 5 and total_fuc == 0 and total_neu5ac == 0:
        return f"High-Mannose-N{total_hexnac}{total_hex}{total_fuc}{total_neu5ac}"
    
    return f"N{total_hexnac}{total_hex}{total_fuc}{total_neu5ac}"

def analyze_branching(G):
    """Analyze the branching pattern of the glycan structure"""
    node_branching = {}
    
    # Analyze branching for each node type and layer
    for node, out_degree in G.out_degree():
        node_type = G.nodes[node]["type"]
        node_layer = G.nodes[node]["layer"]
        
        # Create a key that combines type and layer
        key = f"{node_type}-L{node_layer}"
        
        if key not in node_branching:
            node_branching[key] = []
        
        node_branching[key].append(out_degree)
    
    # Create a summary of branching
    branching_summary = {}
    for key, degrees in node_branching.items():
        branching_summary[key] = {
            "count": len(degrees),
            "average_branches": sum(degrees) / len(degrees) if degrees else 0,
            "max_branches": max(degrees) if degrees else 0
        }
    
    # Special analysis for terminal mannose (nodes 4 and 5)
    if 4 in G.nodes and 5 in G.nodes:
        branching_summary["terminal_mannose"] = {
            "left_mannose_branches": G.out_degree(4),
            "right_mannose_branches": G.out_degree(5)
        }
    
    return branching_summary

def save_glycan_image(composition, output_dir="./"):
    """Build a glycan structure and save its visualization as an image"""
    import os
    G = build_glycan(composition)
    fig = visualize_glycan_snfg(G)
    
    # Create name for the glycan
    name = get_glycan_name(composition)
    
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Save figure
    filename = os.path.join(output_dir, f"glycan_{name}.png")
    fig.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    # Count nodes by type
    node_types = {}
    for _, data in G.nodes(data=True):
        node_type = data["type"]
        if node_type not in node_types:
            node_types[node_type] = 0
        node_types[node_type] += 1
    
    # Get the actual number of each monosaccharide in the final structure
    actual_composition = (
        node_types.get("HexNAc", 0),
        node_types.get("Hex", 0),
        node_types.get("Fuc", 0),
        node_types.get("Neu5Ac", 0)
    )
    
    # Analyze branching patterns
    branching_analysis = analyze_branching(G)
    
    # Return comprehensive information about the glycan
    info = {
        "filename": filename,
        "requested_composition": composition,
        "actual_composition": actual_composition,
        "branching_analysis": branching_analysis
    }
    
    return info

def calculate_mass(composition, permethylation = False):
    """
    Calculate the mass of the glycan based on its composition.
    
    Parameters:
        composition: A tuple of (total_hexnac, total_hex, total_fuc, total_neu5ac)
    
    Returns:
        The mass of the glycan.
    """
    if permethylation:
        # Adjust the masses for permethylation
        mass_hexnac = 245.1263
        mass_hex = 204.0998
        mass_fuc = 174.0892
        mass_neu5ac = 361.1737
    else:
        mass_hexnac = 203.0794  
        mass_hex = 162.0528    
        mass_fuc = 146.079 
        mass_neu5ac = 291.0954  
        
    total_mass = (composition[0] * mass_hexnac +
                  composition[1] * mass_hex +
                  composition[2] * mass_fuc +
                  composition[3] * mass_neu5ac)
    if permethylation:
        total_mass = total_mass + 15.0235 + 47.0497
    else:
        total_mass = total_mass + 20.0262

    return total_mass

def get_fragment_composition(fragment_graph):
    """
    Get the composition of a fragment as a string representation.
    
    Parameters:
        fragment_graph: NetworkX DiGraph representing a glycan fragment
    
    Returns:
        String representation of the fragment composition
    """
    # Count nodes by type
    node_types = {}
    for _, data in fragment_graph.nodes(data=True):
        node_type = data["type"]
        if node_type not in node_types:
            node_types[node_type] = 0
        node_types[node_type] += 1
    
    # Create a composition string
    hexnac = node_types.get("HexNAc", 0)
    hex = node_types.get("Hex", 0)
    fuc = node_types.get("Fuc", 0)
    neu5ac = node_types.get("Neu5Ac", 0)
    
    return f"HexNAc{hexnac}Hex{hex}Fuc{fuc}Neu5Ac{neu5ac}"

def visualize_fragments(fragments, output_dir="./fragments/"):
    """
    Visualize all fragments and save them as images.
    
    Parameters:
        fragments: Dictionary of fragment graphs
        output_dir: Directory to save the fragment images
    
    Returns:
        Dictionary mapping fragment IDs to their image filenames
    """
    import os
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Dictionary to store filenames
    filenames = {}
    
    # Define color scheme for glycan components
    colors = {
        "HexNAc": "blue",
        "Hex": "green",
        "Fuc": "red",
        "Neu5Ac": "purple"
    }
    
    for fragment_id, G in fragments.items():
        # Create position layout based on layers
        pos = {}
        layers = {}
        
        # Group nodes by layer
        for node, data in G.nodes(data=True):
            layer = data["layer"]
            if layer not in layers:
                layers[layer] = []
            layers[layer].append(node)
        
        # Position nodes by layer
        max_layer = max(layers.keys()) if layers else 0
        for layer, nodes in layers.items():
            y = max_layer - layer  # Reverse y-axis to have root at top
            x_step = 1.0 / (len(nodes) + 1) if len(nodes) > 0 else 0.5
            for i, node in enumerate(sorted(nodes)):
                pos[node] = ((i + 1) * x_step, y)
        
        # Create figure and axis
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Draw edges first (so they're behind nodes)
        for u, v in G.edges():
            ax.annotate("",
                       xy=pos[v], xycoords='data',
                       xytext=pos[u], textcoords='data',
                       arrowprops=dict(arrowstyle="-|>", color="black", 
                                      shrinkA=12, shrinkB=12,
                                      connectionstyle="arc3,rad=0.0"))
        
        # Node size
        node_size = 0.04
        
        # Draw nodes with appropriate shapes and colors
        for node, data in G.nodes(data=True):
            x, y = pos[node]
            sugar_type = data["type"]
            color = colors[sugar_type]
            
            if sugar_type == "HexNAc":
                # Blue square for HexNAc
                square = plt.Rectangle((x - node_size/2, y - node_size/2), 
                                      node_size, node_size, 
                                      facecolor=color, edgecolor='black')
                ax.add_patch(square)
            
            elif sugar_type == "Hex":
                # Green circle for Hex
                circle = plt.Circle((x, y), node_size/2, 
                                   facecolor=color, edgecolor='black')
                ax.add_patch(circle)
            
            elif sugar_type == "Fuc":
                # Red triangle for Fucose
                triangle_vertices = [
                    (x, y + node_size/2),
                    (x - node_size/2, y - node_size/2),
                    (x + node_size/2, y - node_size/2)
                ]
                triangle = plt.Polygon(triangle_vertices, 
                                      facecolor=color, edgecolor='black')
                ax.add_patch(triangle)
            
            elif sugar_type == "Neu5Ac":
                # Purple diamond for Sialic acid
                diamond_vertices = [
                    (x, y + node_size/2),
                    (x + node_size/2, y),
                    (x, y - node_size/2),
                    (x - node_size/2, y)
                ]
                diamond = plt.Polygon(diamond_vertices, 
                                     facecolor=color, edgecolor='black')
                ax.add_patch(diamond)
            
            # Add node ID for reference
            ax.text(x, y, str(node), fontsize=8, ha='center', va='center', color='white')
        
        # Create legend
        legend_elements = [
            mpatches.Rectangle((0,0), 1, 1, facecolor=colors["HexNAc"], edgecolor='black', label='HexNAc'),
            mpatches.Circle((0.5,0.5), 0.5, facecolor=colors["Hex"], edgecolor='black', label='Hex'),
            mpatches.Polygon([(0.5,1), (0,0), (1,0)], facecolor=colors["Fuc"], edgecolor='black', label='Fucose'),
            mpatches.Polygon([(0.5,1), (1,0.5), (0.5,0), (0,0.5)], facecolor=colors["Neu5Ac"], edgecolor='black', label='Sialic acid')
        ]
        ax.legend(handles=legend_elements, loc='upper right')
        
        # Set axis limits with some padding
        ax.set_xlim(-0.1, 1.1)
        ax.set_ylim(-0.1, max_layer + 0.1)
        
        # Show layer numbers on y-axis
        yticks = list(range(max_layer + 1))
        ax.set_yticks(yticks)
        ax.set_yticklabels([f"Layer {max_layer-y}" for y in yticks])
        
        # Set title
        ax.set_title(f"Glycan Fragment: {fragment_id}")
        
        # Remove x ticks
        ax.set_xticks([])
        
        plt.tight_layout()
        
        # Save the figure
        filename = os.path.join(output_dir, f"{fragment_id}.png")
        fig.savefig(filename, dpi=300, bbox_inches='tight')
        plt.close(fig)
        
        filenames[fragment_id] = filename
    
    return filenames

import networkx as nx # Make sure nx is imported

def fragment_glycans(G, permethylation=False):
    """
    Generates standard 'b' and 'y' ions resulting ONLY from single glycosidic bond cleavages.

    This function iterates through each bond, calculates the resulting 'b' ion
    (non-reducing side fragment) and 'y' ion (reducing side fragment),
    and stores their information. It does not produce internal fragments or
    fragments resulting from multiple simultaneous cuts.

    Parameters:
        G: NetworkX DiGraph representing the glycan structure.
           Assumes nodes have 'type' and 'layer' attributes.
           Assumes edges represent glycosidic bonds pointing from parent to child
           (e.g., towards the non-reducing end).
        permethylation: Boolean indicating if masses should be calculated for
                       permethylated fragments using calculate_fragment_mass.

    Returns:
        Dictionary of fragment information keyed by unique fragment ID
        (e.g., "b-HexNAc1Hex1-1", "y-HexNAc2Hex2Fuc1-1").
        Includes 'b' and 'y' fragments only, each resulting from 1 cut.
    """
    original_graph = G  # Use the original graph directly for analysis
    fragments = {}
    unique_fragments_tracker = {} # To handle structural isomers giving same composition

    # Get the root node (reducing end, node with in-degree 0)
    try:
        # Ensure only one root node exists or handle appropriately if needed
        root_nodes = [n for n in original_graph.nodes() if original_graph.in_degree(n) == 0]
        if len(root_nodes) != 1:
            # If multiple roots, maybe pick the one with the lowest ID or raise error
            # Sticking to the assumption of a single root for N-glycans usually node 1
             if 1 in root_nodes:
                 root_node = 1
             else:
                 raise ValueError(f"Expected exactly one root node (in-degree 0) or node 1 as root, found {len(root_nodes)}")
        else:
            root_node = root_nodes[0]

    except IndexError:
        raise ValueError("Could not find root node (node with in-degree 0) in the graph.")

    # Iterate through each edge representing a single glycosidic bond to cleave
    for u, v in original_graph.edges():
        # --- Generate 'b' fragment (non-reducing side) ---
        # Nodes in the 'b' fragment are the child 'v' and all its descendants
        b_fragment_nodes = nx.descendants(original_graph, v)
        b_fragment_nodes.add(v) # Add the node v itself

        # Create the subgraph for the 'b' fragment
        b_fragment_graph = original_graph.subgraph(b_fragment_nodes).copy()

        # Calculate properties for 'b' fragment
        b_composition_str = get_fragment_composition(b_fragment_graph)
        b_mass = calculate_fragment_mass(b_fragment_graph, "b", permethylation)
        b_nodes_tuple = tuple(sorted(list(b_fragment_graph.nodes())))
        b_unique_key = ("b", b_composition_str, b_nodes_tuple)

        # Add the fragment if it's unique based on type, composition, and structure
        if b_unique_key not in unique_fragments_tracker:
            base_id = f"b-{b_composition_str}"
            count = 1
            fragment_id = f"{base_id}-{count}"
            while fragment_id in fragments:
                count += 1
                fragment_id = f"{base_id}-{count}"

            fragments[fragment_id] = {
                "graph": b_fragment_graph,
                "mass": b_mass,
                "composition": b_composition_str,
                "type": "b",
                "cuts": 1, # Always 1 cut in this function
                "cut_edges": {(u, v)}
            }
            unique_fragments_tracker[b_unique_key] = fragment_id

        # --- Generate 'y' fragment (reducing side) ---
        y_fragment_nodes = set(original_graph.nodes()) - b_fragment_nodes
        if not y_fragment_nodes:
             continue

        y_fragment_graph = original_graph.subgraph(y_fragment_nodes).copy()

        y_composition_str = get_fragment_composition(y_fragment_graph)
        y_mass = calculate_fragment_mass(y_fragment_graph, "y", permethylation)
        y_nodes_tuple = tuple(sorted(list(y_fragment_graph.nodes())))
        y_unique_key = ("y", y_composition_str, y_nodes_tuple)

        if y_unique_key not in unique_fragments_tracker:
            base_id = f"y-{y_composition_str}"
            count = 1
            fragment_id = f"{base_id}-{count}"
            while fragment_id in fragments:
                count += 1
                fragment_id = f"{base_id}-{count}"

            fragments[fragment_id] = {
                "graph": y_fragment_graph,
                "mass": y_mass,
                "composition": y_composition_str,
                "type": "y",
                "cuts": 1, # Always 1 cut in this function
                "cut_edges": {(u, v)}
            }
            unique_fragments_tracker[y_unique_key] = fragment_id

    return fragments

def get_path_nodes(G, source, target):
    """Helper to get nodes in the shortest path. Returns None if no path."""
    try:
        return nx.shortest_path(G, source=source, target=target)
    except nx.NetworkXNoPath: return None
    except nx.NodeNotFound: return None

def are_edges_on_different_branches(G, edge1, edge2, root_node):
    """Checks if two edges originate from different branches relative to the root."""
    u1, _ = edge1; u2, _ = edge2
    if u1 == u2: return G.out_degree(u1) > 1
    path1_nodes = get_path_nodes(G, root_node, u1)
    path2_nodes = get_path_nodes(G, root_node, u2)
    if path1_nodes is None or path2_nodes is None: return False
    last_common_node = None; len1, len2 = len(path1_nodes), len(path2_nodes); idx = 0
    while idx < len1 and idx < len2 and path1_nodes[idx] == path2_nodes[idx]:
        last_common_node = path1_nodes[idx]; idx += 1
    if last_common_node is None: return False
    if last_common_node == u1 or last_common_node == u2: return False
    if G.out_degree(last_common_node) > 1:
        node_after_common1 = path1_nodes[idx] if idx < len1 else None
        node_after_common2 = path2_nodes[idx] if idx < len2 else None
        return node_after_common1 != node_after_common2
    else: return False
# --- End helper function definitions ---

def calculate_fragment_mass(fragment_graph, fragment_type, permethylation=False):
    """
    Calculate the mass of a fragment based on fragment type.
    
    Parameters:
        fragment_graph: NetworkX DiGraph representing a glycan fragment
        fragment_type: String indicating the type of fragment ("b", "y", or "by")
        permethylation: Boolean indicating whether the fragment is permethylated
    
    Returns:
        Calculated mass of the fragment
    """
    # Count nodes by type
    node_types = {}
    for _, data in fragment_graph.nodes(data=True):
        node_type = data["type"]
        if node_type not in node_types:
            node_types[node_type] = 0
        node_types[node_type] += 1
    
    # Get composition tuple (HexNAc, Hex, Fuc, Neu5Ac)
    composition = (
        node_types.get("HexNAc", 0),
        node_types.get("Hex", 0),
        node_types.get("Fuc", 0),
        node_types.get("Neu5Ac", 0)
    )
    
    # Use the existing calculate_mass function
    mass = calculate_mass(composition, permethylation)
    
    # Apply adjustments based on fragment type
    if fragment_type == "b":
        # Subtract 62.0732 for the reducing end for b fragments
        mass -= 62.0732
        # Add 14.0156 for b-ion
        mass += 14.0156
    elif fragment_type == "y":
        # Subtract 14.0156 for y-ion
        mass -= 14.0156
    elif fragment_type == "by":
        # For true "by" fragments, only subtract reducing end mass
        mass -= 62.0732
        # No ion-type adjustment for mixed fragmentation
    
    return mass

# # Main script for DataFrame creation and Excel export
# import pandas as pd

# # Build a glycan
# composition = (4, 5, 0, 2)  # Example: HexNAc=4, Hex=5, Fuc=0, Neu5Ac=2
# G = build_glycan(composition)

# # Fragment the glycan
# fragments = fragment_glycans(G, permethylation=True)

# # Create lists to store the data
# fragment_ids = []
# masses = []
# compositions = []
# types = []
# cuts = []

# # Extract data from fragments dictionary
# for fragment_id, info in fragments.items():
#     fragment_ids.append(fragment_id)
#     masses.append(round(info['mass'], 4))
#     compositions.append(info['composition'])
#     types.append(info['type'])
#     cuts.append(info['cuts'])

# # Create a DataFrame
# fragments_df = pd.DataFrame({
#     'Fragment ID': fragment_ids,
#     'Mass (Da)': masses,
#     'Composition': compositions,
#     'Type': types,
#     'Cuts': cuts
# })

# # Sort by mass
# fragments_df = fragments_df.sort_values(by='Mass (Da)')

# # Save to Excel
# excel_file = 'glycan_fragments.xlsx'
# fragments_df.to_excel(excel_file, index=False)

# print(f"Fragment data saved to {excel_file}")

# # Display the DataFrame
# print(fragments_df)

def fragment_glycans_yy(G, permethylation=False):
    """
    Improved function to fragment glycans focusing on y-like ions from branches after layer 2.
    
    Parameters:
        G: NetworkX DiGraph representing the glycan structure.
        permethylation: Boolean indicating if masses should be calculated for permethylated fragments.
    
    Returns:
        Dictionary of fragment graphs keyed by unique fragment ID.
        This version returns the graphs directly to be compatible with visualize_fragments.
    """
    fragments = {}
    unique_fragments_tracker = {}
    fragment_info = {}  # Store additional information separately
    
    # Get the root node (reducing end, node with in-degree 0)
    root_nodes = [n for n in G.nodes() if G.in_degree(n) == 0]
    if len(root_nodes) != 1:
        if 1 in G.nodes and G.in_degree(1) == 0:
            root_node = 1
        else:
            raise ValueError(f"Expected exactly one root node, found {len(root_nodes)}")
    else:
        root_node = root_nodes[0]
    
    # Get the first mannose node (layer 2)
    layer2_nodes = [n for n in G.nodes() if G.nodes[n].get('layer') == 2]
    if not layer2_nodes:
        print("No nodes found at layer 2.")
        return fragments
    
    # Get all edges from branches after layer 2
    branch_edges = []
    for u, v in G.edges():
        if G.nodes[u].get('layer') >= 2:  # Edge parent is at layer 2 or higher
            branch_edges.append((u, v))
    
    # Process all pairs of edges for double fragmentation
    for i in range(len(branch_edges)):
        for j in range(i+1, len(branch_edges)):
            edge1 = branch_edges[i]
            edge2 = branch_edges[j]
            
            # Check if these edges are from different branches
            is_different_branches = are_edges_on_different_branches(G, edge1, edge2, root_node)
            
            # Only process if edges are from different branches
            if is_different_branches:
                # Create a temporary graph with these edges removed
                G_temp = G.copy()
                G_temp.remove_edges_from([edge1, edge2])
                
                # Find the connected component containing the root
                root_component = None
                try:
                    for component in nx.weakly_connected_components(G_temp):
                        if root_node in component:
                            root_component = component
                            break
                    else:
                        continue  # Root not found in any component
                except Exception as e:
                    print(f"Error finding components: {e}")
                    continue
                
                # Create the y-fragment from the root component
                y_fragment_graph = G.subgraph(root_component).copy()
                if len(y_fragment_graph.nodes()) <= 1:
                    continue  # Skip if fragment is too small
                
                # Get composition and check if this is a unique fragment
                y_composition_str = get_fragment_composition(y_fragment_graph)
                y_nodes_tuple = tuple(sorted(list(y_fragment_graph.nodes())))
                y_unique_key = ("y", y_composition_str, y_nodes_tuple)
                
                if y_unique_key not in unique_fragments_tracker:
                    # Calculate mass with proper adjustment for multi-branch cut
                    node_types = {}
                    for _, data in y_fragment_graph.nodes(data=True):
                        node_type = data["type"]
                        node_types[node_type] = node_types.get(node_type, 0) + 1
                    
                    composition_tuple = (
                        node_types.get("HexNAc", 0),
                        node_types.get("Hex", 0),
                        node_types.get("Fuc", 0),
                        node_types.get("Neu5Ac", 0)
                    )
                    
                    raw_mass = calculate_mass(composition_tuple, permethylation)
                    
                    # Double cut in different branches gets double adjustment
                    y_ion_2_cut_multi_branch_adjustment = -2 * 14.0156
                    final_mass = raw_mass + y_ion_2_cut_multi_branch_adjustment
                    
                    # Generate a unique ID for this fragment
                    layer1 = G.nodes[edge1[1]].get('layer')
                    layer2 = G.nodes[edge2[1]].get('layer')
                    base_id = f"yL{layer1}+L{layer2}-{y_composition_str}"
                    count = 1
                    fragment_id = f"{base_id}-{count}"
                    while fragment_id in fragments:
                        count += 1
                        fragment_id = f"{base_id}-{count}"
                    
                    # Store the fragment graph directly in fragments (for visualize_fragments compatibility)
                    fragments[fragment_id] = y_fragment_graph
                    
                    # Store additional info in a separate dictionary
                    fragment_info[fragment_id] = {
                        "mass": final_mass,
                        "composition": y_composition_str,
                        "type": "y",
                        "cuts": 2,
                        "multi_branch_cut": True,
                        "effective_cuts": 2,
                        "cut_edges": {edge1, edge2},
                        "cut_layer": (layer1, layer2)
                    }
                    
                    unique_fragments_tracker[y_unique_key] = fragment_id
    
    print(f"Total unique y-like fragments from multiple branches: {len(fragments)}")
    return fragments, fragment_info  # Return both the graphs and their info

# Helper function to verify branch points and branches
def identify_branches(G, root_node):
    """
    Helper function to identify branch points and branches in a glycan structure.
    
    Parameters:
        G: NetworkX DiGraph representing the glycan structure.
        root_node: The root node of the glycan.
        
    Returns:
        Dictionary with branch points and their branches.
    """
    branch_points = {}
    
    for node in G.nodes():
        if G.out_degree(node) > 1:
            # This is a branch point
            branch_points[node] = {
                'layer': G.nodes[node].get('layer'),
                'type': G.nodes[node].get('type'),
                'branches': []
            }
            
            # Get the direct children of this branch point
            for _, child in G.edges(node):
                # For each child, find all descendants
                descendants = nx.descendants(G, child)
                descendants.add(child)  # Include the child itself
                branch_points[node]['branches'].append({
                    'root': child,
                    'nodes': descendants
                })
    
    return branch_points

# Function to visualize branches for debugging
def visualize_branches(G, branch_points, root_node):
    """
    Create a visualization showing branch points and different branches in different colors.
    
    Parameters:
        G: NetworkX DiGraph representing the glycan structure.
        branch_points: Dictionary with branch point information.
        root_node: The root node of the glycan.
        
    Returns:
        Matplotlib figure with visualization.
    """
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    import random
    
    # Create position layout based on layers
    pos = {}
    layers = {}
    
    # Group nodes by layer
    for node, data in G.nodes(data=True):
        layer = data["layer"]
        if layer not in layers:
            layers[layer] = []
        layers[layer].append(node)
    
    # Position nodes by layer
    max_layer = max(layers.keys()) if layers else 0
    for layer, nodes in layers.items():
        y = max_layer - layer  # Reverse y-axis to have root at top
        x_step = 1.0 / (len(nodes) + 1)
        for i, node in enumerate(sorted(nodes)):
            pos[node] = ((i + 1) * x_step, y)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Assign colors to branches
    colors = list(mcolors.TABLEAU_COLORS)
    random.shuffle(colors)
    
    # Draw edges with branch coloring
    for u, v in G.edges():
        edge_color = 'black'  # Default color
        
        # Check if this edge is part of a branch
        for bp, bp_data in branch_points.items():
            for i, branch in enumerate(bp_data['branches']):
                if v in branch['nodes']:
                    edge_color = colors[i % len(colors)]
                    break
            if edge_color != 'black':
                break
        
        ax.annotate("",
                   xy=pos[v], xycoords='data',
                   xytext=pos[u], textcoords='data',
                   arrowprops=dict(arrowstyle="-|>", color=edge_color, 
                                  shrinkA=12, shrinkB=12,
                                  connectionstyle="arc3,rad=0.0"))
    
    # Draw nodes
    node_size = 0.08
    colors = {
        "HexNAc": "blue",
        "Hex": "green",
        "Fuc": "red",
        "Neu5Ac": "purple"
    }
    
    for node, data in G.nodes(data=True):
        x, y = pos[node]
        sugar_type = data["type"]
        color = colors[sugar_type]
        
        if sugar_type == "HexNAc":
            square = plt.Rectangle((x - node_size/2, y - node_size/2), 
                                  node_size, node_size, 
                                  facecolor=color, edgecolor='black')
            ax.add_patch(square)
        elif sugar_type == "Hex":
            circle = plt.Circle((x, y), node_size/2, 
                               facecolor=color, edgecolor='black')
            ax.add_patch(circle)
        elif sugar_type == "Fuc":
            triangle_vertices = [
                (x, y + node_size/2),
                (x - node_size/2, y - node_size/2),
                (x + node_size/2, y - node_size/2)
            ]
            triangle = plt.Polygon(triangle_vertices, 
                                  facecolor=color, edgecolor='black')
            ax.add_patch(triangle)
        elif sugar_type == "Neu5Ac":
            diamond_vertices = [
                (x, y + node_size/2),
                (x + node_size/2, y),
                (x, y - node_size/2),
                (x - node_size/2, y)
            ]
            diamond = plt.Polygon(diamond_vertices, 
                                 facecolor=color, edgecolor='black')
            ax.add_patch(diamond)
        
        # Add node ID
        ax.text(x, y, str(node), fontsize=8, ha='center', va='center', color='white')
        
        # Highlight branch points
        if node in branch_points:
            plt.scatter(x, y, s=300, facecolors='none', edgecolors='yellow', linewidths=2)
    
    # Set axis limits with some padding
    ax.set_xlim(-0.1, 1.1)
    ax.set_ylim(-0.1, max_layer + 0.1)
    
    # Show layer numbers on y-axis
    yticks = list(range(max_layer + 1))
    ax.set_yticks(yticks)
    ax.set_yticklabels([f"Layer {max_layer-y}" for y in yticks])
    
    # Set title
    ax.set_title("Glycan Structure with Branch Points Highlighted")
    
    # Remove x ticks
    ax.set_xticks([])
    
    plt.tight_layout()
    return fig

# Test function to validate the improved fragmentation
def test_branch_detection(composition):
    """
    Test function to validate branch detection and fragmentation.
    
    Parameters:
        composition: Tuple (HexNAc, Hex, Fuc, Neu5Ac)
        
    Returns:
        Fragments and branch information.
    """
    # Build the glycan
    G = build_glycan(composition)
    
    # Get the root node
    root_nodes = [n for n in G.nodes() if G.in_degree(n) == 0]
    root_node = root_nodes[0] if root_nodes else 1
    
    # Identify branches
    branch_points = identify_branches(G, root_node)
    
    # Test the improved fragmentation
    fragments, fragment_info = fragment_glycans_yy(G, permethylation=False)
    
    # Print branch information
    print(f"Glycan with composition {composition}:")
    print(f"Found {len(branch_points)} branch points:")
    for bp, data in branch_points.items():
        print(f"  Node {bp} (layer {data['layer']}, {data['type']}) has {len(data['branches'])} branches")
        for i, branch in enumerate(data['branches']):
            print(f"    Branch {i+1}: {len(branch['nodes'])} nodes starting with {branch['root']}")
    
    print(f"\nFound {len(fragments)} y-like fragments from multiple branches")
    
    return {
        'glycan': G,
        'root_node': root_node,
        'branch_points': branch_points,
        'fragments': fragments,
        'fragment_info': fragment_info
    }
