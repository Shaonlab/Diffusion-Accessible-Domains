# Diffusion-Accessible-Domains: Bioinformatics analysis

This repository contains all the files required to replicate the figure 4(c) results (for `H3K27me3` histone- modification) from Singh and Chakrabarti (`https://www.biorxiv.org/content/10.1101/2022.08.31.505992v1.full`). 

## Installation

See this on how to install Conda -> `https://conda.io/projects/conda/en/latest/user-guide/install/index.html`

To create an environment using the .yml file in this repository:
```bash
conda env create --file=environment.yml
```

## Usage
After installing `conda`, follow the steps below to run `script.r`:
```bash
# Clone this repository to your local computer
git clone https://github.com/Shaonlab/Diffusion-Accessible-Domains.git

# Create a conda environment using the provided yml file
conda env create --file=/path/to/environment.yml

# Activate the conda environment
conda activate hic_analysis

# Run the analysis script (-a and -b options indicate the median D(i1,i2) values for H3K27me3, without taking into account the HiC information)
Rscript script.r -e /path/to/repearly_1Mb.bed -l /path/to/replate_1Mb.bed  -c /path/to/HiC_matrices/ -o /path/to/output/directory -m H3K27me3 -a 0.434 -b 0.617

# to know more about the parameters passed to the R-script
Rscript script.r --help
```
## Downloading the HiC data
The HiC data is sourced from the original HiC paper -> `https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2858594/`

To download the processed HiC data required to run this script: 

-> search for `GSE18199` on `GEO` (`https://www.ncbi.nlm.nih.gov/geo/`) 

-> download `GSE18199_binned_heatmaps.zip.gz` from the supplementary file section

-> decompress the `.gz` file using tools like `7-Zip` (`https://www.7-zip.org/`)
