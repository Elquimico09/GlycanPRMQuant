# N-Glycan Builder GUI

Interactive graphical interface for building and visualizing N-glycan structures following strict biosynthetic rules.

![N-Glycan Builder](https://via.placeholder.com/800x600.png?text=N-Glycan+Builder+GUI)

## Features

- **Interactive Input**: Adjust HexNAc, Hex, Fuc, and Sia counts with spin boxes
- **Real-time Visualization**: See structures update as you build them
- **Preset Glycans**: Quick access to common N-glycan structures
- **Structure Validation**: Automatic validation of biosynthetic rules
- **Layer Information**: Visual layer labels and detailed breakdown
- **Custom Color Scheme**:
  - **HexNAc**: Blue squares ◼️
  - **Hex (core, layers ≤3)**: Green circles ⬤
  - **Hex (antenna, layers >3)**: Yellow circles ⬤
  - **Fuc**: Red triangles ▲
  - **Sia**: Magenta diamonds ◆

## Installation

### Standard Installation
```bash
# Install the package
pip install -e .

# Install GUI dependencies
pip install -e .[gui]
```

### Manual Installation
```bash
# Install required packages
pip install PyQt5 matplotlib networkx
```

## Usage

### Method 1: Command Line Launcher
```bash
python launch_glycan_builder.py
```

### Method 2: Console Script (after installation)
```bash
glycan-builder-gui
```

### Method 3: Python Import
```python
from glycanPRMQuant.glycanBuilderGUI import launch_gui
launch_gui()
```

### Method 4: Matplotlib-Only (No PyQt5)
```python
from glycanPRMQuant.glycanBuilderGUI import simple_visualize

# Visualize a glycan: HexNAc=4, Hex=5, Fuc=1, Sia=2
simple_visualize(4, 5, 1, 2)
```

## GUI Components

### Input Panel (Left Side)

#### Residue Count Spinners
- **HexNAc**: Total N-acetylhexosamine residues (minimum 2)
- **Hex**: Total hexose residues (minimum 3)
- **Fuc**: Total fucose residues (0-10)
- **Sia**: Total sialic acid residues (0-10)

#### Build Button
Click to construct and visualize the glycan structure

#### Preset Buttons
Quick access to common N-glycan compositions:
- **Core (2-3-0-0)**: Minimal N-glycan core
- **Man5 (2-5-0-0)**: High-mannose Man5
- **Man9 (2-9-0-0)**: High-mannose Man9
- **Complex (4-5-1-2)**: Typical complex-type biantennary
- **Triantennary (5-6-1-3)**: Triantennary complex glycan

#### Structure Info Panel
Displays:
- Composition notation (e.g., 4-5-1-2)
- Total residue count
- Number of nodes and edges
- Validation status
- Layer-by-layer breakdown
- Warning messages (if any)

### Visualization Panel (Right Side)

Shows the glycan structure with:
- **Nodes**: Colored shapes representing sugar residues
- **Edges**: Black lines representing glycosidic bonds
- **Node Labels**: Numerical IDs for each residue
- **Layer Labels**: Right-side layer annotations
- **Legend**: Color/shape key for residue types

## Biosynthetic Rules

The builder follows strict N-glycan biosynthetic layering:

```
Layer 0: 1 HexNAc (reducing end) - blue square
Layer 1: 1 HexNAc (chitobiose) - blue square (no Fuc allowed!)
Layer 2: 1 Hex (core mannose) - green circle
Layer 3: 2 Hex (branching mannoses) - green circles
Layer 4: Additional HexNAc (max 4) - blue squares
Layer 5: Hex on layer 4 HexNAc - yellow circles
Layer 6+: Sia on layer 5+ Hex - magenta diamonds
```

### Special Rules
1. **Fuc** can attach to any HexNAc **except** layer 1
2. **Sia** can only attach to Hex in layer 5 or later
3. Core structure (layers 0-3) is always built first
4. Maximum 4 HexNAc in layer 4, balanced between arms

## Example Workflows

### Building a Simple Complex Glycan

1. Set composition:
   - HexNAc: 4
   - Hex: 5
   - Fuc: 1
   - Sia: 2

2. Click "Build Glycan"

3. View structure:
   - Blue squares: 4 HexNAc (layers 0, 1, 4, 4)
   - Green circles: 3 Hex (layers 2, 3, 3)
   - Yellow circles: 2 Hex (layer 5)
   - Red triangle: 1 Fuc (typically on layer 0)
   - Magenta diamonds: 2 Sia (layer 6)

### Using Presets

1. Click "Complex (4-5-1-2)" preset button
2. Structure automatically builds and displays
3. Review structure info panel for details
4. Modify counts if desired and rebuild

### Validation Warnings

If you see warnings:
- **"Could only add X/Y Sia"**: Need more Hex in layer 5+ for Sia attachment
- **"Could only add X/Y Fuc"**: Not enough HexNAc with available capacity
- **"Used only X/Y residues"**: Structure cannot accommodate all requested residues

## Color Scheme Reference

| Residue Type | Shape | Color | Layer Restrictions |
|--------------|-------|-------|-------------------|
| HexNAc | Square ◼️ | Blue | All layers |
| Hex (core) | Circle ⬤ | Green | Layers 0-3 |
| Hex (antenna) | Circle ⬤ | Yellow | Layers 4+ |
| Fuc | Triangle ▲ | Red | Any (except on layer 1 HexNAc) |
| Sia | Diamond ◆ | Magenta | Layer 6+ (on layer 5+ Hex) |

## Keyboard Shortcuts

- **Up/Down Arrows**: Adjust focused spinner
- **Enter**: Build glycan
- **Tab**: Navigate between inputs

## Troubleshooting

### GUI Won't Launch
```bash
# Check PyQt5 installation
python -c "import PyQt5; print('PyQt5 OK')"

# If error, reinstall
pip install --upgrade PyQt5
```

### "Cannot build glycan" Error
- Check minimum requirements: HexNAc ≥ 2, Hex ≥ 3
- Reduce counts if exceeding capacity
- Check error message for specific issue

### Visualization Issues
- Resize window if structure is cut off
- Use matplotlib-only fallback if PyQt5 issues persist
- Update matplotlib: `pip install --upgrade matplotlib`

### Import Errors
```bash
# Install all dependencies
pip install PyQt5 matplotlib networkx pandas numpy scipy
```

## Advanced Usage

### Programmatic Access
```python
from glycanPRMQuant.glycanBuilder_improved import build_nglycan_strict
from glycanPRMQuant.glycanBuilderGUI import GlycanCanvas
import matplotlib.pyplot as plt

# Build structure
G = build_nglycan_strict(hexnac_count=6, hex_count=7, fuc_count=2, sia_count=4)

# Create canvas
fig, ax = plt.subplots(figsize=(12, 10))
canvas = GlycanCanvas()
canvas.axes = ax
canvas.plot_glycan(G, "Custom Glycan: 6-7-2-4")
plt.show()
```

### Export Structure
```python
import networkx as nx

# Build glycan
G = build_nglycan_strict(4, 5, 1, 2)

# Export to various formats
nx.write_gml(G, "glycan.gml")  # GraphML
nx.write_graphml(G, "glycan.graphml")  # GraphML
nx.write_edgelist(G, "glycan.edges")  # Edge list

# Get adjacency matrix
adj = nx.to_numpy_array(G)
```

## Requirements

- Python ≥ 3.6
- PyQt5 ≥ 5.15.0 (for GUI)
- matplotlib ≥ 3.0
- networkx ≥ 2.0

## Related Files

- `glycanBuilder_improved.py`: Core builder logic
- `NGLYCAN_RULES.md`: Detailed biosynthetic rules
- `IMPROVEMENTS.md`: Package improvements summary

## Future Enhancements

Planned features:
- [ ] Export to image (PNG, SVG, PDF)
- [ ] Save/load glycan structures
- [ ] IUPAC nomenclature display
- [ ] Fragment prediction overlay
- [ ] 3D structure visualization
- [ ] Batch processing mode
- [ ] Glycan database browser

## Support

For issues or questions:
- Check `NGLYCAN_RULES.md` for biosynthetic rules
- Review structure info panel for validation messages
- Use matplotlib-only fallback if GUI issues persist

## License

MIT License - See package LICENSE file

## Citation

If you use this tool in your research, please cite:
```
glycanPRMQuant: A package for glycan PRM quantification
Author: Vishal Sandilya
Email: vishal.sandilya@ttu.edu
```
