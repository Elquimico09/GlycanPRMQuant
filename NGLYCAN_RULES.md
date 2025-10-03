# N-Glycan Biosynthetic Rules - Layer-by-Layer Construction

## Core Structure (Always Required)

```
         Layer 6+: [Sia] (on Hex from layer 5+)
                    |
         Layer 5:  [Hex] (on HexNAc from layer 4)
                    |
         Layer 4:  [HexNAc] (max 4, on Hex from layer 3)
                  /   \
         Layer 3: [Hex]  [Hex] (3-arm and 6-arm)
                  \   /
         Layer 2:  [Hex] (core mannose)
                    |
         Layer 1:  [HexNAc] (cannot have Fuc)
                    |
         Layer 0:  [HexNAc] (reducing end)
```

## Detailed Rules

### Layer 0 (Reducing End)
- **Always**: 1 HexNAc
- **Function**: Reducing end of N-glycan
- **Can attach Fuc**: ✅ Yes

### Layer 1 (Chitobiose Core)
- **Always**: 1 HexNAc
- **Function**: Forms chitobiose core with layer 0
- **Can attach Fuc**: ❌ NO (Special restriction!)

### Layer 2 (Central Mannose)
- **Always**: 1 Hex (Mannose)
- **Function**: Central branching point
- **Children**: Exactly 2 (layer 3 mannoses)

### Layer 3 (Branching Mannoses)
- **Always**: 2 Hex (Mannoses)
- **Names**: 3-arm and 6-arm
- **Function**: Primary branch points for complex glycans
- **Can have**: Up to 2 children each

### Layer 4 (Antenna HexNAc)
- **Content**: Additional HexNAc beyond core (max 4 in this layer)
- **Attachment**: To layer 3 Hex (alternates between 3-arm and 6-arm)
- **Can attach Fuc**: ✅ Yes
- **Typical**: 0-4 HexNAc in complex/hybrid glycans

### Layer 5 (Antenna Hex/Galactose)
- **Content**: Hex residues (typically Galactose)
- **Attachment**: MUST attach to layer 4 HexNAc first
- **Can attach Sia**: ✅ Yes (primary site)
- **Function**: Forms LacNAc antennae

### Layer 6+ (Terminal Modifications)
- **Content**: Sialic acid (Neu5Ac/Neu5Gc)
- **Attachment**: To Hex in layer 5 or later ONLY
- **Restriction**: Cannot attach directly to layer 3 or 4

## Special Rules

### Fucosylation (Fuc)
- ✅ **Can attach to**: Any HexNAc EXCEPT layer 1
- ✅ **Common sites**:
  - Layer 0 HexNAc (core fucosylation)
  - Layer 4 HexNAc (antenna fucosylation)
- ❌ **Cannot attach to**: Layer 1 HexNAc (strict rule)
- **Max per HexNAc**: 1 Fuc per HexNAc node

### Sialylation (Neu5Ac/Neu5Gc)
- ✅ **Can attach to**: Hex in layer 5 or later
- ❌ **Cannot attach to**: Hex in layers 2, 3, or 4
- **Reasoning**: Ensures proper LacNAc formation before sialylation

### Extended Structures (HexNAc > 6)
- **When**: More than 6 total HexNAc (2 core + 4 layer 4)
- **Action**: Add additional HexNAc to layer 5 Hex nodes
- **New layer**: Creates layer 6 HexNAc
- **Can extend**: Additional Hex can attach to these layer 6 HexNAc

## Example Structures

### Example 1: Simple Complex (HexNAc=4, Hex=5, Fuc=1, Sia=2)
```
Composition: 4-5-1-2

                    [Sia]     [Sia]      Layer 6
                      |         |
              [Hex]-------[Hex]          Layer 5
                |           |
            [HexNAc]    [HexNAc]         Layer 4
                 \       /
              [Hex]   [Hex]              Layer 3
                  \ /
                 [Hex]                   Layer 2
                   |
               [HexNAc]                  Layer 1 (no Fuc!)
                   |
            [HexNAc]--[Fuc]              Layer 0 (core fucosylation)
```

### Example 2: High-Mannose (HexNAc=2, Hex=9, Fuc=0, Sia=0)
```
Composition: 2-9-0-0

                [Hex]   [Hex]            Layer 6
                  |       |
              [Hex]   [Hex]              Layer 5
                  |       |
              [Hex]   [Hex]              Layer 4
                  \   /
              [Hex] [Hex]                Layer 3
                  \ /
                 [Hex]                   Layer 2
                   |
               [HexNAc]                  Layer 1
                   |
               [HexNAc]                  Layer 0
```

### Example 3: Tri-antennary (HexNAc=5, Hex=6, Fuc=1, Sia=3)
```
Composition: 5-6-1-3

      [Sia]        [Sia]      [Sia]      Layer 7
        |            |          |
      [Hex]        [Hex]      [Hex]      Layer 6
        |            |          |
    [HexNAc]     [HexNAc]  [HexNAc]      Layer 4/5
           \        |      /
           [Hex]  [Hex]                  Layer 3
               \  /
              [Hex]                      Layer 2
                |
            [HexNAc]                     Layer 1 (no Fuc!)
                |
         [HexNAc]--[Fuc]                 Layer 0 (core fucosylation)
```

### Example 4: Extended Structure (HexNAc=7, Hex=8, Fuc=0, Sia=2)
```
Composition: 7-8-0-2
(Demonstrates HexNAc >6 rule)

                    [Sia]     [Sia]                Layer 8
                      |         |
              [Hex]-------[Hex]                    Layer 7
                |           |
            [HexNAc]    [HexNAc]                   Layer 6 (extra HexNAc!)
                |           |
              [Hex]-------[Hex]                    Layer 5
                |           |           |
            [HexNAc]    [HexNAc]    [HexNAc]       Layer 4
                 \          |       /
              [Hex]-------[Hex]                    Layer 3
                      \ /
                     [Hex]                         Layer 2
                       |
                   [HexNAc]                        Layer 1
                       |
                   [HexNAc]                        Layer 0
```

## Algorithm Flow

```python
# Pseudocode for build_nglycan_strict()

1. Validate: Check HexNAc >= 2, Hex >= 3

2. Build Core (Layers 0-3):
   - Add HexNAc to layer 0 (reducing end)
   - Add HexNAc to layer 1 (chitobiose)
   - Add Hex to layer 2 (central mannose)
   - Add 2 Hex to layer 3 (3-arm and 6-arm)

3. Add Layer 4 HexNAc (up to 4):
   - Alternate between 3-arm and 6-arm parents
   - Track layer 4 HexNAc nodes for next step

4. Add Layer 5 Hex:
   - First, attach to all layer 4 HexNAc nodes
   - Then, add any remaining Hex to available positions
   - Track layer 5 Hex nodes for Sia

5. Handle Extended HexNAc (>6):
   - Attach to layer 5 Hex nodes
   - Creates layer 6 HexNAc
   - Can add additional Hex to these if needed

6. Add Sialic Acid:
   - Only to Hex nodes in layer 5 or later
   - Priority: layer 5 Hex first

7. Add Fucose:
   - To any HexNAc EXCEPT layer 1
   - Typical: layer 0 (core) or layer 4 (antenna)
   - Max 1 per HexNAc (respect out-degree limit)

8. Validate final structure
```

## Key Differences from Old Builder

| Feature | Old Builder | New Strict Builder |
|---------|-------------|-------------------|
| Layer 4 composition | Mixed HexNAc/Hex | Strictly HexNAc (max 4) |
| Layer 5 composition | Variable | Strictly Hex first |
| Sia placement | Complex priority | Simple: layer 5+ Hex only |
| Fuc on layer 1 HexNAc | Allowed with check | Forbidden (hard constraint) |
| Extended structures | Limited | Full support (HexNAc >6) |
| Validation | None | Built-in validator |

## Usage in Code

```python
from glycanPRMQuant.glycanBuilder_improved import (
    build_nglycan_strict,
    validate_nglycan_structure
)

# Build a glycan
G = build_nglycan_strict(
    hexnac_count=4,
    hex_count=5,
    fuc_count=1,
    sia_count=2
)

# Validate it follows rules
valid, errors = validate_nglycan_structure(
    G,
    {'HexNAc': 4, 'Hex': 5, 'Fuc': 1, 'Neu5Ac': 2}
)

if not valid:
    print("Validation errors:")
    for error in errors:
        print(f"  - {error}")

# Access node information
for node in G.nodes():
    node_data = G.nodes[node]
    print(f"Node {node}: {node_data['type']} at layer {node_data['layer']}")
```

## Common Glycan Compositions

| Composition | Type | Description |
|-------------|------|-------------|
| 2-3-0-0 | Core | Minimal N-glycan core |
| 2-5-0-0 | High-mannose | Man5 |
| 2-9-0-0 | High-mannose | Man9 |
| 3-3-0-0 | Hybrid | One antenna |
| 4-5-0-0 | Complex | Biantennary, asialylated |
| 4-5-1-0 | Complex | Core-fucosylated biantennary |
| 4-5-0-2 | Complex | Disialylated biantennary |
| 4-5-1-2 | Complex | Core-fucosylated, disialylated |
| 5-6-1-3 | Complex | Triantennary, trisialylated |
| 6-7-1-4 | Complex | Tetraantennary, tetrasialylated |
