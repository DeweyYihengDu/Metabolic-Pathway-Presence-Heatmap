from scipy.cluster.hierarchy import linkage, dendrogram
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import dendrogram as original_dendrogram
import heatmap as ht

x='human'

ht.main(x, x+'_matrix.csv',x+'_matrix.pdf')

pathway_matrix = pd.read_csv(x+'_matrix.csv', index_col=0)
linked = linkage(pathway_matrix, method='average', metric='euclidean')

# Determine the number of bacteria
num_bacteria = len(pathway_matrix.index)

# Calculate the figure height based on the number of bacteria
figure_height = max(7, num_bacteria * 0.3) # 0.2 is the scaling factor; adjust it as needed

# Visualize the dendrogram
plt.figure(figsize=(10, figure_height))
dendrogram(linked, labels=pathway_matrix.index, orientation='right', leaf_font_size=10)
plt.title("Bacterial Phylogenetic Tree using UPGMA Method")
plt.xlabel("Euclidean Distance")
plt.savefig(x+'_phylogenetic_tree.pdf', format='pdf', bbox_inches='tight')

# Show the figure
plt.show()
