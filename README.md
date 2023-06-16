## ## Citation
if you use `DeepTarget`, please consider adding the following
citation (bioRxiv preprint, [direct link
here](https://www.biorxiv.org/content/10.1101/2022.10.17.512424v1)

## Installation
### Github installation
library(devtools)
install_github("CBIIT-CGBB/DeepTarget")

### Bioconductor installation

In order to install the package from Bioconductor, make sure
`BiocManager` is installed and then install the package
`DeepTarget`:

    if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

    BiocManager::install("DeepTarget")
## Usage

### Vignettes

`DeepTarget` performs deep characterization of cancer drugs mechanism of action by integrating large-scale genetic and drug screens from public datasets and aims to find drug's primary target(s), predict whether the drug specifically targets the wild-type or mutated target forms, and predict the secondary target(s) that mediate its response when the primary target is not expressed.

### Data
 For the purpose of demonstration, we provided a OntargetM object that contains a small subset of data. For full datasets used in the manuscript, please contact us at tinh.nguyen@nih.gov. For the detail of the links for downloading the data, please do ??DeepTarget.

### Brief examples for interesting drugs from the paper.

Please refer to the file named V_analyses.Rmd in the directory vignettes


