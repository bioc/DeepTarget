## Installation

### Bioconductor installation

In order to install the package from Bioconductor, make sure
`BiocManager` is installed and then install the package
`DeepTarget`:

    if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

    BiocManager::install("DeepTarget")
## Usage

### Vignettes

`DeepTarget` performs deep characterization of cancer drugs mechanism of action by integrating large-scale genetic and drug screens.

### Data
 For the purpose of demonstration, we provided a OntargetM object that contains a small subset of data. For full datasets used in the manuscript, please contact us at tinh.nguyen@nih.gov. For the detail of the links for downloading the data, please do ??DeepTarget.

### Brief example for 2 interesting drugs. Please note that each drug has two assays.
Please refer the script and user manual in the inst folder


