"""
N-Glycan Builder GUI

Interactive GUI for building and visualizing N-glycan structures
following strict biosynthetic rules with SNFG-like notation.
"""

import sys
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.patches as mpatches

try:
    from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout,
                                  QHBoxLayout, QLabel, QSpinBox, QPushButton,
                                  QGroupBox, QMessageBox, QTextEdit)
    from PyQt5.QtCore import Qt
    PYQT_AVAILABLE = True
except ImportError:
    PYQT_AVAILABLE = False
    print("Warning: PyQt5 not available. GUI will use matplotlib only.")

from .glycanBuilder_improved import build_nglycan_strict, validate_nglycan_structure, generate_all_isomers


class GlycanCanvas(FigureCanvas):
    """Canvas for drawing glycan structures."""

    def __init__(self, parent=None, width=10, height=8, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        super(GlycanCanvas, self).__init__(self.fig)
        self.setParent(parent)

    def plot_glycan(self, G, title="N-Glycan Structure"):
        """Plot glycan structure with custom color scheme."""
        self.axes.clear()

        if G is None or G.number_of_nodes() == 0:
            self.axes.text(0.5, 0.5, 'Enter composition and click "Build Glycan"',
                          ha='center', va='center', fontsize=14)
            self.axes.set_xlim(0, 1)
            self.axes.set_ylim(0, 1)
            self.axes.axis('off')
            self.draw()
            return

        # Calculate positions
        pos = self._calculate_positions(G)

        # Draw edges first (so they appear behind nodes)
        for u, v in G.edges():
            x1, y1 = pos[u]
            x2, y2 = pos[v]
            self.axes.plot([x1, x2], [y1, y2], 'k-', linewidth=2, zorder=1)

        # Draw nodes with custom shapes and colors
        node_size = 0.15
        for node in G.nodes():
            x, y = pos[node]
            node_type = G.nodes[node]['type']
            layer = G.nodes[node]['layer']

            # Draw shape based on type
            if node_type == "HexNAc":
                # Blue square
                self._draw_square(x, y, node_size, 'blue', node)
            elif node_type == "Hex":
                # Green circle if: layer <= 3 OR has incoming connection from another Hex
                # Otherwise yellow circle
                color = 'yellow'
                if layer <= 3:
                    color = 'green'
                else:
                    # Check if any parent is a Hex
                    for parent in G.predecessors(node):
                        if G.nodes[parent]['type'] == 'Hex':
                            color = 'green'
                            break
                self._draw_circle(x, y, node_size, color, node)
            elif node_type == "Fuc":
                # Red triangle
                self._draw_triangle(x, y, node_size, 'red', node)
            elif node_type == "Neu5Ac":
                # Magenta diamond
                self._draw_diamond(x, y, node_size, 'magenta', node)

        # Add layer labels on the right
        layers = sorted(set(G.nodes[n]['layer'] for n in G.nodes()))
        for layer in layers:
            # Find average y position for this layer
            layer_nodes = [n for n in G.nodes() if G.nodes[n]['layer'] == layer]
            if layer_nodes:
                y_positions = [pos[n][1] for n in layer_nodes]
                avg_y = sum(y_positions) / len(y_positions)
                self.axes.text(1.1, avg_y, f'Layer {layer}',
                              fontsize=10, va='center', style='italic')

        # Set axis properties
        self.axes.set_aspect('equal')
        self.axes.axis('off')
        self.axes.set_title(title, fontsize=16, fontweight='bold', pad=20)

        # Add legend
        legend_elements = [
            mpatches.Rectangle((0, 0), 1, 1, facecolor='blue', edgecolor='black', label='HexNAc'),
            mpatches.Circle((0.5, 0.5), 0.5, facecolor='green', edgecolor='black', label='Hex (core)'),
            mpatches.Circle((0.5, 0.5), 0.5, facecolor='yellow', edgecolor='black', label='Hex (antenna)'),
            mpatches.RegularPolygon((0.5, 0.5), 3, radius=0.5, facecolor='red', edgecolor='black', label='Fuc'),
            mpatches.RegularPolygon((0.5, 0.5), 4, radius=0.5, facecolor='magenta', edgecolor='black', label='Sia')
        ]
        self.axes.legend(handles=legend_elements, loc='upper left', fontsize=10)

        self.draw()

    def _calculate_positions(self, G):
        """Calculate node positions using hierarchical layout."""
        pos = {}
        root_node = [n for n in G.nodes() if G.in_degree(n) == 0][0]

        # Recursive layout
        def layout_tree(node, x, y, width, visited):
            if node in visited:
                return
            visited.add(node)
            pos[node] = (x, y)

            children = list(G.successors(node))
            if not children:
                return

            # Calculate positions for children
            child_width = width / len(children)
            start_x = x - width / 2 + child_width / 2

            for i, child in enumerate(children):
                child_x = start_x + i * child_width
                child_y = y - 0.25  # Move down for next layer
                layout_tree(child, child_x, child_y, child_width * 0.8, visited)

        # Start layout from root
        layout_tree(root_node, 0.5, 1.0, 1.0, set())

        return pos

    def _draw_square(self, x, y, size, color, label):
        """Draw a square (HexNAc)."""
        square = mpatches.Rectangle((x - size/2, y - size/2), size, size,
                                    facecolor=color, edgecolor='black',
                                    linewidth=2, zorder=2)
        self.axes.add_patch(square)
        self.axes.text(x, y, str(label), ha='center', va='center',
                      fontsize=8, color='white', fontweight='bold', zorder=3)

    def _draw_circle(self, x, y, size, color, label):
        """Draw a circle (Hex)."""
        circle = mpatches.Circle((x, y), size/2, facecolor=color,
                                edgecolor='black', linewidth=2, zorder=2)
        self.axes.add_patch(circle)
        self.axes.text(x, y, str(label), ha='center', va='center',
                      fontsize=8, color='white', fontweight='bold', zorder=3)

    def _draw_triangle(self, x, y, size, color, label):
        """Draw a triangle (Fuc)."""
        triangle = mpatches.RegularPolygon((x, y), 3, radius=size/2,
                                          facecolor=color, edgecolor='black',
                                          linewidth=2, zorder=2)
        self.axes.add_patch(triangle)
        self.axes.text(x, y, str(label), ha='center', va='center',
                      fontsize=8, color='white', fontweight='bold', zorder=3)

    def _draw_diamond(self, x, y, size, color, label):
        """Draw a diamond (Sialic acid)."""
        diamond = mpatches.RegularPolygon((x, y), 4, radius=size/2,
                                         orientation=0.785398,  # 45 degrees
                                         facecolor=color, edgecolor='black',
                                         linewidth=2, zorder=2)
        self.axes.add_patch(diamond)
        self.axes.text(x, y, str(label), ha='center', va='center',
                      fontsize=8, color='white', fontweight='bold', zorder=3)


class GlycanBuilderGUI(QMainWindow):
    """Main GUI window for N-glycan builder."""

    def __init__(self):
        super().__init__()
        self.glycan = None
        self.isomers = []
        self.current_isomer_idx = 0
        self.init_ui()

    def init_ui(self):
        """Initialize the user interface."""
        self.setWindowTitle('N-Glycan Builder')
        self.setGeometry(100, 100, 1200, 800)

        # Central widget
        central_widget = QWidget()
        self.setCentralWidget(central_widget)

        # Main layout
        main_layout = QHBoxLayout(central_widget)

        # Left panel - Controls
        left_panel = self._create_control_panel()
        main_layout.addWidget(left_panel, 1)

        # Right panel - Visualization
        right_panel = self._create_visualization_panel()
        main_layout.addWidget(right_panel, 3)

    def _create_control_panel(self):
        """Create the control panel with input fields."""
        panel = QWidget()
        layout = QVBoxLayout(panel)

        # Title
        title = QLabel('N-Glycan Composition')
        title.setStyleSheet('font-size: 16px; font-weight: bold; padding: 10px;')
        layout.addWidget(title)

        # Input group
        input_group = QGroupBox('Residue Counts')
        input_layout = QVBoxLayout()

        # HexNAc input
        hexnac_layout = QHBoxLayout()
        hexnac_layout.addWidget(QLabel('HexNAc:'))
        self.hexnac_spin = QSpinBox()
        self.hexnac_spin.setRange(2, 20)
        self.hexnac_spin.setValue(4)
        self.hexnac_spin.setStyleSheet('font-size: 14px;')
        hexnac_layout.addWidget(self.hexnac_spin)
        input_layout.addLayout(hexnac_layout)

        # Hex input
        hex_layout = QHBoxLayout()
        hex_layout.addWidget(QLabel('Hex:'))
        self.hex_spin = QSpinBox()
        self.hex_spin.setRange(3, 20)
        self.hex_spin.setValue(5)
        self.hex_spin.setStyleSheet('font-size: 14px;')
        hex_layout.addWidget(self.hex_spin)
        input_layout.addLayout(hex_layout)

        # Fuc input
        fuc_layout = QHBoxLayout()
        fuc_layout.addWidget(QLabel('Fuc:'))
        self.fuc_spin = QSpinBox()
        self.fuc_spin.setRange(0, 10)
        self.fuc_spin.setValue(1)
        self.fuc_spin.setStyleSheet('font-size: 14px;')
        fuc_layout.addWidget(self.fuc_spin)
        input_layout.addLayout(fuc_layout)

        # Sia input
        sia_layout = QHBoxLayout()
        sia_layout.addWidget(QLabel('Sia:'))
        self.sia_spin = QSpinBox()
        self.sia_spin.setRange(0, 10)
        self.sia_spin.setValue(2)
        self.sia_spin.setStyleSheet('font-size: 14px;')
        sia_layout.addWidget(self.sia_spin)
        input_layout.addLayout(sia_layout)

        input_group.setLayout(input_layout)
        layout.addWidget(input_group)

        # Build button
        self.build_button = QPushButton('Build Glycan')
        self.build_button.setStyleSheet('''
            QPushButton {
                background-color: #4CAF50;
                color: white;
                font-size: 16px;
                font-weight: bold;
                padding: 10px;
                border-radius: 5px;
            }
            QPushButton:hover {
                background-color: #45a049;
            }
        ''')
        self.build_button.clicked.connect(self.build_glycan)
        layout.addWidget(self.build_button)

        # Preset buttons
        preset_group = QGroupBox('Common Glycans')
        preset_layout = QVBoxLayout()

        presets = [
            ('Core (2-3-0-0)', 2, 3, 0, 0),
            ('Man5 (2-5-0-0)', 2, 5, 0, 0),
            ('Man9 (2-9-0-0)', 2, 9, 0, 0),
            ('Complex (4-5-1-2)', 4, 5, 1, 2),
            ('Triantennary (5-6-1-3)', 5, 6, 1, 3),
        ]

        for name, hexnac, hex_val, fuc, sia in presets:
            btn = QPushButton(name)
            btn.clicked.connect(lambda checked, h1=hexnac, h2=hex_val, f=fuc, s=sia:
                              self.load_preset(h1, h2, f, s))
            preset_layout.addWidget(btn)

        preset_group.setLayout(preset_layout)
        layout.addWidget(preset_group)

        # Isomer navigation
        isomer_group = QGroupBox('Isomer Navigation')
        isomer_layout = QVBoxLayout()

        self.isomer_label = QLabel('No isomers')
        self.isomer_label.setStyleSheet('font-weight: bold; padding: 5px;')
        isomer_layout.addWidget(self.isomer_label)

        nav_layout = QHBoxLayout()
        self.prev_button = QPushButton('◀ Previous')
        self.prev_button.clicked.connect(self.prev_isomer)
        self.prev_button.setEnabled(False)
        nav_layout.addWidget(self.prev_button)

        self.next_button = QPushButton('Next ▶')
        self.next_button.clicked.connect(self.next_isomer)
        self.next_button.setEnabled(False)
        nav_layout.addWidget(self.next_button)

        isomer_layout.addLayout(nav_layout)
        isomer_group.setLayout(isomer_layout)
        layout.addWidget(isomer_group)

        # Info text
        self.info_text = QTextEdit()
        self.info_text.setReadOnly(True)
        self.info_text.setMaximumHeight(200)
        self.info_text.setPlaceholderText('Structure information will appear here...')
        layout.addWidget(QLabel('Structure Info:'))
        layout.addWidget(self.info_text)

        layout.addStretch()

        return panel

    def _create_visualization_panel(self):
        """Create the visualization panel."""
        panel = QWidget()
        layout = QVBoxLayout(panel)

        # Canvas for glycan structure
        self.canvas = GlycanCanvas(panel, width=10, height=8)
        layout.addWidget(self.canvas)

        return panel

    def load_preset(self, hexnac, hex_val, fuc, sia):
        """Load a preset glycan composition."""
        self.hexnac_spin.setValue(hexnac)
        self.hex_spin.setValue(hex_val)
        self.fuc_spin.setValue(fuc)
        self.sia_spin.setValue(sia)
        self.build_glycan()

    def build_glycan(self):
        """Build and display all glycan isomers."""
        hexnac = self.hexnac_spin.value()
        hex_val = self.hex_spin.value()
        fuc = self.fuc_spin.value()
        sia = self.sia_spin.value()

        try:
            # Generate all isomers
            self.isomers = generate_all_isomers(hexnac, hex_val, fuc, sia)
            self.current_isomer_idx = 0

            if len(self.isomers) == 0:
                raise ValueError("No valid isomers could be generated")

            # Update navigation
            self.isomer_label.setText(f'Isomer 1 of {len(self.isomers)}')
            self.prev_button.setEnabled(len(self.isomers) > 1)
            self.next_button.setEnabled(len(self.isomers) > 1)

            # Display first isomer
            self.display_current_isomer()

        except ValueError as e:
            QMessageBox.critical(self, 'Error', f'Cannot build glycan:\n{str(e)}')
            self.info_text.setPlainText(f'Error: {str(e)}')

    def prev_isomer(self):
        """Show previous isomer."""
        if self.isomers and self.current_isomer_idx > 0:
            self.current_isomer_idx -= 1
            self.isomer_label.setText(f'Isomer {self.current_isomer_idx + 1} of {len(self.isomers)}')
            self.display_current_isomer()

    def next_isomer(self):
        """Show next isomer."""
        if self.isomers and self.current_isomer_idx < len(self.isomers) - 1:
            self.current_isomer_idx += 1
            self.isomer_label.setText(f'Isomer {self.current_isomer_idx + 1} of {len(self.isomers)}')
            self.display_current_isomer()

    def display_current_isomer(self):
        """Display the currently selected isomer."""
        if not self.isomers:
            return

        self.glycan = self.isomers[self.current_isomer_idx]
        hexnac = self.hexnac_spin.value()
        hex_val = self.hex_spin.value()
        fuc = self.fuc_spin.value()
        sia = self.sia_spin.value()

        # Validate
        expected = {'HexNAc': hexnac, 'Hex': hex_val, 'Fuc': fuc, 'Neu5Ac': sia}
        valid, errors = validate_nglycan_structure(self.glycan, expected)

        # Update visualization
        title = f"N-Glycan: {hexnac}-{hex_val}-{fuc}-{sia} (Isomer {self.current_isomer_idx + 1}/{len(self.isomers)})"
        self.canvas.plot_glycan(self.glycan, title)

        # Update info text
        info = self._generate_structure_info(hexnac, hex_val, fuc, sia, valid, errors)
        self.info_text.setHtml(info)

    def _generate_structure_info(self, hexnac, hex_val, fuc, sia, valid, errors):
        """Generate HTML info about the structure."""
        html = f'''
        <h3>Composition: {hexnac}-{hex_val}-{fuc}-{sia}</h3>
        <p><b>Total residues:</b> {hexnac + hex_val + fuc + sia}</p>
        <p><b>Nodes:</b> {self.glycan.number_of_nodes()}</p>
        <p><b>Edges:</b> {self.glycan.number_of_edges()}</p>
        <p><b>Validation:</b> <span style="color: {'green' if valid else 'red'};">
        {'✓ Valid' if valid else '✗ Warnings'}</span></p>
        '''

        if errors:
            html += '<p><b>Warnings:</b></p><ul>'
            for err in errors:
                html += f'<li>{err}</li>'
            html += '</ul>'

        # Isomer-specific features
        html += '<p><b>Isomer Features:</b></p><ul>'

        # Determine glycan type
        layer3_nodes = [n for n in self.glycan.nodes() if self.glycan.nodes[n]['layer'] == 3]
        is_hybrid = False
        glycan_type = 'Unknown'

        if hexnac == 2:
            glycan_type = 'High Mannose'
        else:
            # Check if hybrid (one arm has Hex, one has HexNAc in layer 4)
            arm_types = []
            for n3 in layer3_nodes:
                children = [c for c in self.glycan.successors(n3) if self.glycan.nodes[c]['layer'] == 4]
                child_types = set([self.glycan.nodes[c]['type'] for c in children])
                arm_types.append(child_types)

            has_hex_arm = any('Hex' in types for types in arm_types)
            has_hexnac_arm = any('HexNAc' in types for types in arm_types)

            if has_hex_arm and has_hexnac_arm:
                glycan_type = 'Hybrid'
                is_hybrid = True
            elif has_hexnac_arm:
                glycan_type = 'Complex'

        html += f'<li><b>Type:</b> {glycan_type}</li>'

        # Check fucosylation type
        fuc_types = []
        for node in self.glycan.nodes():
            if self.glycan.nodes[node]['type'] == 'Fuc':
                parents = list(self.glycan.predecessors(node))
                if parents:
                    p = parents[0]
                    p_layer = self.glycan.nodes[p]['layer']
                    if p_layer == 0:
                        fuc_types.append('Core fucose')
                    else:
                        fuc_types.append(f'Antenna fucose (layer {p_layer})')

        if fuc_types:
            for ft in fuc_types:
                html += f'<li>{ft}</li>'
        else:
            html += '<li>No fucosylation</li>'

        # Check antenna distribution
        layer4_nodes = [n for n in self.glycan.nodes() if self.glycan.nodes[n]['layer'] == 4]
        layer3_nodes = [n for n in self.glycan.nodes() if self.glycan.nodes[n]['layer'] == 3]

        if len(layer4_nodes) > 0:
            # Count antennae per arm
            n3_arms = {}
            for n3 in layer3_nodes:
                children = [c for c in self.glycan.successors(n3) if self.glycan.nodes[c]['layer'] == 4]
                n3_arms[n3] = len(children)

            arm_counts = list(n3_arms.values())
            html += f'<li>Antenna distribution: {"-".join(map(str, arm_counts))}</li>'

        html += '</ul>'

        # Layer breakdown
        html += '<p><b>Layer Breakdown:</b></p><ul>'
        layers = {}
        for node in self.glycan.nodes():
            layer = self.glycan.nodes[node]['layer']
            node_type = self.glycan.nodes[node]['type']
            if layer not in layers:
                layers[layer] = []
            layers[layer].append(node_type)

        for layer in sorted(layers.keys()):
            types = ', '.join(layers[layer])
            html += f'<li>Layer {layer}: {types}</li>'
        html += '</ul>'

        return html


def launch_gui():
    """Launch the N-glycan builder GUI."""
    if not PYQT_AVAILABLE:
        print("ERROR: PyQt5 is required for the GUI.")
        print("Install with: pip install PyQt5")
        return

    app = QApplication(sys.argv)
    app.setStyle('Fusion')  # Modern look
    window = GlycanBuilderGUI()
    window.show()
    sys.exit(app.exec_())


# Simple matplotlib-only version (fallback)
def simple_visualize(hexnac, hex_val, fuc, sia):
    """Simple visualization without PyQt5 (matplotlib only)."""
    try:
        G = build_nglycan_strict(hexnac, hex_val, fuc, sia)

        fig, ax = plt.subplots(figsize=(12, 10))
        canvas = GlycanCanvas(None, width=12, height=10)
        canvas.axes = ax
        canvas.plot_glycan(G, f"N-Glycan: {hexnac}-{hex_val}-{fuc}-{sia}")
        plt.tight_layout()
        plt.show()

    except Exception as e:
        print(f"Error: {e}")


if __name__ == '__main__':
    if PYQT_AVAILABLE:
        launch_gui()
    else:
        # Fallback to simple matplotlib visualization
        print("PyQt5 not available, using simple matplotlib visualization")
        print("Example: Complex glycan (4-5-1-2)")
        simple_visualize(4, 5, 1, 2)
