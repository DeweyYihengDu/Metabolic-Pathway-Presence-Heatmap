# Metabolic-Pathway-Presence-Heatmap

This project generates a heatmap for the presence of metabolic pathways across different species within a given genus. The heatmap is visualized as a PDF file, and the underlying data is saved as a CSV file. The data is fetched from the [KEGG (Kyoto Encyclopedia of Genes and Genomes) API](https://www.kegg.jp/kegg/rest/keggapi.html).

## Prerequisites
To run this project, you will need Python 3.6 or later and the following Python packages:

+ pandas
+ numpy
+ seaborn
+ matplotlib
+ requests
You can install these packages using pip:

```
pip install pandas numpy seaborn matplotlib requests
```

## Usage
Run the script from the command line by providing the genus name, output CSV filename, and output PDF filename:

```
python heatmap.py Vibrio matrix.csv heatmap.pdf
```

Replace `Vibrio` with the genus name you want to analyze, and `matrix.csv` and `heatmap.pdf` with your desired output file names.

Alternatively, you can run the script from an IDE such as VSCode. If you do not provide command-line arguments, the script will use default parameters: "Vibrio" for the genus name, "matrix.csv" for the output matrix CSV, and "heatmap.pdf" for the output heatmap PDF.

## Output
The script generates two output files:

1. A CSV file containing a matrix representation of the presence of metabolic pathways across different species within the specified genus.
2. A PDF file with a heatmap visualization of the data.

## Citation

If you use this project in your research, please cite the following paper:

```
[1]Y.-H. Du and J.-H. Mu, “Metabolic-Pathway-Presence-Heatmap (MPPH):Constructing phylogenetic trees based on metabolic pathways.,” Jun. 2023, doi: https://doi.org/10.1101/2023.06.27.546232.
```