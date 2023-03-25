from scipy.cluster.hierarchy import linkage, dendrogram
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
# Perform hierarchical clustering using UPGMA method
from scipy.cluster.hierarchy import dendrogram as original_dendrogram
import heatmap as ht
ht.main('Delftia', 'Delftia_matrix.csv','Delftia_matrix.pdf')

pathway_matrix = pd.read_csv('Delftia_matrix.csv', index_col=0)
linked = linkage(pathway_matrix, method='average', metric='euclidean')


# def dendrogram(*args, **kwargs):
#     show_distances = kwargs.pop('show_distances', False)
#     ddata = original_dendrogram(*args, **kwargs)

#     if show_distances:
#         icoord = np.array(ddata['icoord'])
#         dcoord = np.array(ddata['dcoord'])
#         dists = np.array(ddata['ivl'])
#         y_dist = np.array(ddata['leaves'])

#         for (i, d, y) in zip(icoord, dcoord, y_dist):
#             x = i[1]
#             y = d[1]
#             plt.text(x, y, f'{pathway_matrix.index[y]} ({y:.2f})', va='center', ha='left')

#     return ddata

# Visualize the dendrogram
# plt.figure(figsize=(10, 10))
# dendrogram(linked, labels=pathway_matrix.index, orientation='right', leaf_font_size=10, show_distances=True)
# plt.title("Bacterial Phylogenetic Tree using UPGMA Method")
# plt.xlabel("Euclidean Distance")
# plt.tight_layout()

# # Save the figure as a PDF file
# plt.savefig('bacterial_phylogenetic_tree_with_distances.pdf', format='pdf', bbox_inches='tight')


# Determine the number of bacteria
num_bacteria = len(pathway_matrix.index)

# Calculate the figure height based on the number of bacteria
figure_height = max(7, num_bacteria * 0.2) # 0.2 is the scaling factor; adjust it as needed

# Visualize the dendrogram
plt.figure(figsize=(10, figure_height))
dendrogram(linked, labels=pathway_matrix.index, orientation='right', leaf_font_size=10)
plt.title("Bacterial Phylogenetic Tree using UPGMA Method")
plt.xlabel("Euclidean Distance")
plt.savefig('Delftia_phylogenetic_tree.pdf', format='pdf', bbox_inches='tight')

# Show the figure
# plt.show()
