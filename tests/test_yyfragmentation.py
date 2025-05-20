from glycanPRMQuant.glycanBuilder import fragment_glycans_yy, fragment_glycans, build_glycan, visualize_glycan_snfg, get_glycan_name, calculate_mass, visualize_fragments
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Example usage
composition = (2, 9, 0, 0) 

# 2. Build the glycan graph from the composition
G = build_glycan(composition)
G.nodes(data =True)

# Visualize the glycan structure
fig = visualize_glycan_snfg(G, composition=composition)
plt.show()

# 3. Generate yy fragments
fragments, fragment_info = fragment_glycans_yy(G, permethylation=True)

# Generate regular fragments
fragments_single = fragment_glycans(G, permethylation=True)

for fragment_id in fragments:
    print(f"Fragment {fragment_id} has mass {fragment_info[fragment_id]['mass']}")

for fragment_id in fragments_single:
    print(f"Fragment {fragment_id} has mass {fragments_single[fragment_id]['mass']}")