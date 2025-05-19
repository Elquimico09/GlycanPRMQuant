from glycanPRMQuant.glycanBuilder import fragment_glycans_yy, fragment_glycans, build_glycan, visualize_glycan_snfg, get_glycan_name, calculate_mass
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

glycan_test = (2, 9, 0, 0)          # HexNAc, Hex, Fuc, Neu5Ac
G = build_glycan(glycan_test)

fig = visualize_glycan_snfg(G, glycan_test)  # returns a Figure
#plt.show()
mass = calculate_mass(glycan_test, permethylation=True)
print("Mass of the glycan:", mass)

frags = fragment_glycans_yy(G, permethylation=True)      # G is your full glycan graph

print(f"Total fragments: {len(frags)}")
print("First 5 IDs:", list(frags)[:5])

from pprint import pprint

frags_all = fragment_glycans_yy(G)      # outer list/tuple
frags = frags_all[0]                    # the actual dict of fragments

# now this works:
for fid, info in frags.items():
    print(fid, info["mass"])
