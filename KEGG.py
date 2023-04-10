import requests
import re

# KEGG API base URL
KEGG_API_BASE = "http://rest.kegg.jp/"


def get_species_list(genus):
    organism_url = KEGG_API_BASE + "list/organism"
    response = requests.get(organism_url)
    lines = response.text.split("\n")
    species_list = []

    for line in lines:
        if line and genus in line:
            parts = line.split("\t")
            species_code = parts[1]
            species_name = parts[2]
            species_list.append((species_code, species_name))

    return species_list


def get_metabolic_pathways(species_code):
    # species_code = species_code
    pathways_url = KEGG_API_BASE + f"list/pathway/{species_code}"
    response = requests.get(pathways_url)
    lines = response.text.split("\n")

    pathways = []
    for line in lines:
        if line:
            parts = line.split("\t")
            pathway_id=parts[0]
            # pathway_id = parts[0].split(":")[1]
            pathways.append(pathway_id)

    return pathways


a=get_species_list("Planctomycetales")
b=[]
for c, d in a:
    b.append(get_metabolic_pathways(c))


import numpy as np
import pandas as pd

def pathways_to_matrix(pathway_list):
    # Flatten the list and extract unique prefixes
    flattened_list = [item for sublist in pathway_list for item in sublist]
    prefixes = sorted(set([item[:item.find("0")] for item in flattened_list]))

    # Extract unique pathway IDs
    pathway_ids = sorted(set([item[item.find("0"):] for item in flattened_list]))

    # Create a zero-filled matrix
    matrix = np.zeros((len(prefixes), len(pathway_ids)))

    # Fill the matrix with 1s for existing pathways
    for pathway in flattened_list:
        prefix = pathway[:pathway.find("0")]
        pathway_id = pathway[pathway.find("0"):]
        row_idx = prefixes.index(prefix)
        col_idx = pathway_ids.index(pathway_id)
        matrix[row_idx][col_idx] = 1

    df = pd.DataFrame(matrix, index=prefixes, columns=pathway_ids)
    return df



pathway_matrix = pathways_to_matrix(b)
print(pathway_matrix)

pathway_matrix.to_csv('matrix.csv', index=True, header=True, index_label='row_name', float_format='%.0f')



import seaborn as sns
import matplotlib.pyplot as plt

def plot_heatmap(matrix):
    plt.figure(figsize=(20, 10))
    sns.heatmap(matrix, cmap='coolwarm', annot=True, fmt=".0f", linewidths=.5)
    plt.xlabel("Pathway ID")
    plt.ylabel("Species Prefix")
    plt.title("Metabolic Pathway Presence Heatmap")
    plt.show()

plot_heatmap(pathway_matrix)


