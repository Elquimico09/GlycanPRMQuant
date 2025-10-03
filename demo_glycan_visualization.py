#!/usr/bin/env python
"""
Demo script for N-Glycan visualization

This script demonstrates both GUI and matplotlib-only visualization modes.
"""

import sys
import os

# Add package to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

def demo_matplotlib_only():
    """Demonstrate matplotlib-only visualization (no PyQt5 required)."""
    print("=" * 60)
    print("N-Glycan Builder - Matplotlib Demo")
    print("=" * 60)

    from glycanPRMQuant.glycanBuilder_improved import build_nglycan_strict
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    # Build a few example glycans
    examples = [
        ("Core N-Glycan", 2, 3, 0, 0),
        ("Complex Biantennary", 4, 5, 1, 2),
        ("High-Mannose Man9", 2, 9, 0, 0),
        ("Triantennary", 5, 6, 1, 3),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(16, 14))
    axes = axes.flatten()

    for idx, (name, hexnac, hex_val, fuc, sia) in enumerate(examples):
        ax = axes[idx]

        print(f"\nBuilding: {name} ({hexnac}-{hex_val}-{fuc}-{sia})")

        try:
            # Build glycan
            G = build_nglycan_strict(hexnac, hex_val, fuc, sia)

            # Calculate positions
            pos = calculate_positions(G)

            # Draw edges
            for u, v in G.edges():
                x1, y1 = pos[u]
                x2, y2 = pos[v]
                ax.plot([x1, x2], [y1, y2], 'k-', linewidth=2, zorder=1)

            # Draw nodes
            node_size = 0.12
            for node in G.nodes():
                x, y = pos[node]
                node_type = G.nodes[node]['type']
                layer = G.nodes[node]['layer']

                if node_type == "HexNAc":
                    # Blue square
                    square = mpatches.Rectangle((x - node_size/2, y - node_size/2),
                                                node_size, node_size,
                                                facecolor='blue', edgecolor='black',
                                                linewidth=2, zorder=2)
                    ax.add_patch(square)
                elif node_type == "Hex":
                    # Green (layer ≤3) or Yellow (layer >3) circle
                    color = 'green' if layer <= 3 else 'yellow'
                    circle = mpatches.Circle((x, y), node_size/2,
                                            facecolor=color, edgecolor='black',
                                            linewidth=2, zorder=2)
                    ax.add_patch(circle)
                elif node_type == "Fuc":
                    # Red triangle
                    triangle = mpatches.RegularPolygon((x, y), 3, radius=node_size/2,
                                                      facecolor='red', edgecolor='black',
                                                      linewidth=2, zorder=2)
                    ax.add_patch(triangle)
                elif node_type == "Neu5Ac":
                    # Magenta diamond
                    diamond = mpatches.RegularPolygon((x, y), 4, radius=node_size/2,
                                                     orientation=0.785398,
                                                     facecolor='magenta', edgecolor='black',
                                                     linewidth=2, zorder=2)
                    ax.add_patch(diamond)

                # Add node label
                ax.text(x, y, str(node), ha='center', va='center',
                       fontsize=8, color='white', fontweight='bold', zorder=3)

            ax.set_aspect('equal')
            ax.axis('off')
            ax.set_title(f"{name}\n({hexnac}-{hex_val}-{fuc}-{sia})",
                        fontsize=14, fontweight='bold')

            print(f"  ✓ Built successfully: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

        except Exception as e:
            print(f"  ✗ Error: {e}")
            ax.text(0.5, 0.5, f"Error:\n{str(e)}", ha='center', va='center',
                   transform=ax.transAxes, fontsize=10, color='red')
            ax.axis('off')

    # Add overall legend
    legend_elements = [
        mpatches.Rectangle((0, 0), 1, 1, facecolor='blue', edgecolor='black', label='HexNAc'),
        mpatches.Circle((0.5, 0.5), 0.5, facecolor='green', edgecolor='black', label='Hex (core)'),
        mpatches.Circle((0.5, 0.5), 0.5, facecolor='yellow', edgecolor='black', label='Hex (antenna)'),
        mpatches.RegularPolygon((0.5, 0.5), 3, 0.5, facecolor='red', edgecolor='black', label='Fuc'),
        mpatches.RegularPolygon((0.5, 0.5), 4, 0.5, facecolor='magenta', edgecolor='black', label='Sia')
    ]
    fig.legend(handles=legend_elements, loc='lower center', ncol=5,
              fontsize=12, frameon=True, bbox_to_anchor=(0.5, -0.02))

    plt.suptitle('N-Glycan Structures - Color/Shape Legend',
                fontsize=16, fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0.03, 1, 0.96])

    print("\n" + "=" * 60)
    print("Displaying visualization...")
    print("Close the window to exit.")
    print("=" * 60)

    plt.show()


def calculate_positions(G):
    """Calculate hierarchical positions for glycan nodes."""
    pos = {}
    root_node = [n for n in G.nodes() if G.in_degree(n) == 0][0]

    def layout_tree(node, x, y, width, visited):
        if node in visited:
            return
        visited.add(node)
        pos[node] = (x, y)

        children = list(G.successors(node))
        if not children:
            return

        child_width = width / len(children)
        start_x = x - width / 2 + child_width / 2

        for i, child in enumerate(children):
            child_x = start_x + i * child_width
            child_y = y - 0.25
            layout_tree(child, child_x, child_y, child_width * 0.8, visited)

    layout_tree(root_node, 0.5, 1.0, 1.0, set())
    return pos


def demo_gui():
    """Launch the full GUI (requires PyQt5)."""
    print("=" * 60)
    print("N-Glycan Builder - GUI Mode")
    print("=" * 60)

    try:
        from glycanPRMQuant.glycanBuilderGUI import launch_gui
        print("Launching GUI...")
        launch_gui()
    except ImportError as e:
        print(f"Error: {e}")
        print("\nPyQt5 is required for GUI mode.")
        print("Install with: pip install PyQt5")
        print("\nFalling back to matplotlib demo...")
        demo_matplotlib_only()


def main():
    """Main entry point."""
    print("\nN-Glycan Builder Demo")
    print("Select mode:")
    print("  1. GUI Mode (requires PyQt5)")
    print("  2. Matplotlib Demo (no PyQt5 required)")
    print("  3. Exit")

    # Check if we can import PyQt5
    try:
        import PyQt5
        pyqt_available = True
    except ImportError:
        pyqt_available = False
        print("\nNote: PyQt5 not available. GUI mode will not work.")

    choice = input("\nEnter choice (1-3): ").strip()

    if choice == "1":
        if not pyqt_available:
            print("PyQt5 not available. Falling back to matplotlib demo...")
            demo_matplotlib_only()
        else:
            demo_gui()
    elif choice == "2":
        demo_matplotlib_only()
    elif choice == "3":
        print("Exiting...")
        sys.exit(0)
    else:
        print("Invalid choice. Using matplotlib demo...")
        demo_matplotlib_only()


if __name__ == '__main__':
    main()
