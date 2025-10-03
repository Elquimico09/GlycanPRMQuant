# N-Glycan Builder GUI - Quick Start Guide

## Installation

```bash
# Navigate to package directory
cd /Users/vishalsandilya/Documents/glycanPRMQuant

# Install the package with GUI support
pip install -e .[gui]

# Or install GUI dependencies manually
pip install PyQt5 matplotlib networkx
```

## Launching the GUI

### Option 1: Direct Python Script (Recommended)
```bash
python launch_glycan_builder.py
```

### Option 2: Demo Script (Interactive)
```bash
python demo_glycan_visualization.py
# Then select option 1 for GUI mode, or 2 for matplotlib-only demo
```

### Option 3: Python Import
```python
from glycanPRMQuant.glycanBuilderGUI import launch_gui
launch_gui()
```

### Option 4: Command Line (After pip install)
```bash
glycan-builder-gui
```

## Using the GUI

### Quick Start

1. **Adjust Residue Counts** (left panel):
   - HexNAc: 2-20 (minimum 2 for N-glycan core)
   - Hex: 3-20 (minimum 3 for N-glycan core)
   - Fuc: 0-10
   - Sia: 0-10

2. **Click "Build Glycan"** button

3. **View Structure** (right panel):
   - Blue squares = HexNAc
   - Green circles = Hex in core (layers ≤3)
   - Yellow circles = Hex in antenna (layers >3)
   - Red triangles = Fuc
   - Magenta diamonds = Sia

4. **Check Structure Info** (bottom left):
   - Composition notation
   - Validation status
   - Layer breakdown
   - Warning messages

### Using Presets

Click any preset button for instant visualization:

- **Core (2-3-0-0)**: Minimal N-glycan structure
- **Man5 (2-5-0-0)**: High-mannose with 5 mannoses
- **Man9 (2-9-0-0)**: High-mannose with 9 mannoses
- **Complex (4-5-1-2)**: Biantennary complex glycan
- **Triantennary (5-6-1-3)**: Triantennary complex glycan

### Understanding the Visualization

#### Node Colors and Shapes
```
◼️ Blue Square      = HexNAc (N-acetylhexosamine)
⬤ Green Circle     = Hex (core mannose, layers 0-3)
⬤ Yellow Circle    = Hex (antenna galactose, layers 4+)
▲ Red Triangle     = Fuc (fucose)
◆ Magenta Diamond  = Sia (sialic acid)
```

#### Layer Labels (Right Side)
Shows which layer each residue belongs to:
- Layer 0: Reducing end HexNAc
- Layer 1: Chitobiose HexNAc
- Layer 2: Core mannose
- Layer 3: Branching mannoses
- Layer 4+: Antenna structures

#### Node Labels (White Numbers)
Unique ID for each node in the graph

### Example Workflow

**Building a Complex Biantennary Glycan (4-5-1-2):**

1. Set values:
   ```
   HexNAc: 4
   Hex:    5
   Fuc:    1
   Sia:    2
   ```

2. Click "Build Glycan"

3. Expected structure:
   ```
   Layer 0: HexNAc (1) ← Fuc (12) attached here (core fucosylation)
   Layer 1: HexNAc (2)
   Layer 2: Hex (3)
   Layer 3: Hex (4), Hex (5)
   Layer 4: HexNAc (6), HexNAc (7)
   Layer 5: Hex (8), Hex (9)
   Layer 6: Neu5Ac (10), Neu5Ac (11)
   ```

4. Check info panel:
   ```
   ✓ Valid
   Total residues: 12
   Nodes: 12
   Edges: 11
   ```

## Troubleshooting

### GUI Won't Start

**Error: "ModuleNotFoundError: No module named 'PyQt5'"**
```bash
pip install PyQt5
```

**Error: "ModuleNotFoundError: No module named 'glycanPRMQuant'"**
```bash
# Make sure you're in the package directory
cd /Users/vishalsandilya/Documents/glycanPRMQuant

# Install in development mode
pip install -e .
```

### Validation Warnings

**"Could not add all requested Fuc"**
- Fuc attaches to HexNAc (except layer 1)
- Each HexNAc can have max 2 children
- If all HexNAc already have 2 children, no more Fuc can attach
- Solution: Increase HexNAc count

**"Could only add X/Y Sia"**
- Sia must attach to Hex in layer 5 or later
- Not enough Hex in antenna positions
- Solution: Increase Hex count to create more antenna positions

**"Used only X/Y HexNAc"**
- Structure cannot accommodate all requested HexNAc
- Layer 4 has max 4 HexNAc
- Additional HexNAc need Hex in layer 5 to attach to
- Solution: Reduce HexNAc or increase Hex

### Display Issues

**Structure appears cut off:**
- Resize the window
- Structure automatically adjusts to canvas size

**Shapes overlap:**
- Normal for complex structures
- Use node labels to identify residues
- Check layer info panel for details

## Alternative: Matplotlib-Only Mode

If PyQt5 causes issues, use matplotlib-only visualization:

```python
from glycanPRMQuant.glycanBuilderGUI import simple_visualize

# Build and display glycan
simple_visualize(hexnac=4, hex_val=5, fuc=1, sia=2)
```

Or run the demo:
```bash
python demo_glycan_visualization.py
# Select option 2 for matplotlib demo
```

## Programmatic Usage

### Build Without GUI
```python
from glycanPRMQuant.glycanBuilder_improved import build_nglycan_strict

# Build structure
G = build_nglycan_strict(
    hexnac_count=4,
    hex_count=5,
    fuc_count=1,
    sia_count=2
)

# Access node information
for node in G.nodes():
    node_type = G.nodes[node]['type']
    layer = G.nodes[node]['layer']
    print(f"Node {node}: {node_type} at layer {layer}")

# Access edges (bonds)
for u, v in G.edges():
    print(f"Bond: {u} → {v}")
```

### Validate Structure
```python
from glycanPRMQuant.glycanBuilder_improved import validate_nglycan_structure

valid, errors = validate_nglycan_structure(
    G,
    {'HexNAc': 4, 'Hex': 5, 'Fuc': 1, 'Neu5Ac': 2}
)

if valid:
    print("✓ Structure is valid")
else:
    print("Validation errors:")
    for error in errors:
        print(f"  - {error}")
```

### Custom Visualization
```python
import matplotlib.pyplot as plt
from glycanPRMQuant.glycanBuilderGUI import GlycanCanvas

# Build glycan
G = build_nglycan_strict(5, 6, 1, 3)

# Create custom canvas
fig, ax = plt.subplots(figsize=(14, 12))
canvas = GlycanCanvas()
canvas.axes = ax
canvas.plot_glycan(G, "My Custom Glycan")
plt.savefig('my_glycan.png', dpi=300, bbox_inches='tight')
plt.show()
```

## Common Glycan Compositions

| Composition | Name | Description |
|-------------|------|-------------|
| 2-3-0-0 | Core | Minimal N-glycan |
| 2-5-0-0 | Man5 | High-mannose |
| 2-9-0-0 | Man9 | High-mannose |
| 3-3-0-0 | Hybrid | One antenna |
| 4-5-0-0 | Complex | Biantennary |
| 4-5-1-0 | Complex-F | Core-fucosylated |
| 4-5-0-2 | Complex-S2 | Disialylated |
| 4-5-1-2 | Complex-F-S2 | Core-fuc, disialylated |
| 5-6-1-3 | Triantennary | Three antennae |
| 6-7-1-4 | Tetraantennary | Four antennae |

## Keyboard Shortcuts

- **Tab**: Navigate between input fields
- **Enter**: Build glycan (when focused on input)
- **↑/↓ Arrows**: Adjust spinner values
- **Ctrl+Q**: Quit application

## Tips and Tricks

1. **Start with presets**: Use preset buttons to understand typical structures
2. **Check layer info**: Use the info panel to verify layer assignments
3. **Incremental building**: Start with core, gradually increase complexity
4. **Validation first**: Check warnings before adjusting composition
5. **Save screenshots**: Use system screenshot tools to save visualizations

## File Locations

- **GUI Code**: `glycanPRMQuant/glycanBuilderGUI.py`
- **Builder Logic**: `glycanPRMQuant/glycanBuilder_improved.py`
- **Launcher**: `launch_glycan_builder.py`
- **Demo**: `demo_glycan_visualization.py`
- **Rules**: `NGLYCAN_RULES.md`

## Getting Help

1. Check biosynthetic rules: `NGLYCAN_RULES.md`
2. Review improvements: `IMPROVEMENTS.md`
3. Read GUI README: `GUI_README.md`
4. Check structure info panel for validation messages
5. Use matplotlib-only mode if GUI issues persist

## What's Next?

After building structures, you can:
- Export NetworkX graph for further analysis
- Use in fragmentation studies
- Integrate with MS/MS matching
- Build glycan databases
- Perform structural analysis

## Support

For issues:
- Check error messages in structure info panel
- Review validation warnings
- Consult `NGLYCAN_RULES.md` for biosynthetic rules
- Use matplotlib-only fallback if needed

Contact: vishal.sandilya@ttu.edu
