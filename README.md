
# DeepTarget
`DeepTarget` performs deep characterization of cancer drugs mechanism of action by integrating large-scale genetic and drug screens from the Cancer Dependency Map project run by the Broad Institute (depmap.org)  and aims to find drug's primary target(s), predict whether the drug specifically targets the wild-type or mutated target forms, and predict the secondary target(s) that mediate its response when the primary target is not expressed.
# Depmap license
The user of this package is required to agree to the terms and conditions of DepMap portal (https://depmap.org/portal/). 
Some of these terms and conditions are problematic for U.S. Federal Government employees, and they should consult their technology transfer office/legal office before agreeing to such terms and conditions.
## Installation
### GitHub installation
You can install the development version of DeepTarget from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("CBIIT-CGBB/DeepTarget")
```

### Bioconductor installation
In order to install the package from Bioconductor, make sure
`BiocManager` is installed and then install the package
`DeepTarget`:
``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DeepTarget")
```
## Example

Please refer to the file named DeepTarget_Vignette.Rmd in the directory vignettes for a demonstration of how the package can be used.

## Data
For the purpose of demonstration, we include an OntargetM object that contains a subset of data from the Cancer Dependency Map (depmap.org). The full data set it too large for package memory constraints. For details and links to download the full data set, please use the  `??OntargetM` command.

## Citation
If you use `DeepTarget`, please consider adding the following
citation from bioRxiv preprint: (https://www.biorxiv.org/content/10.1101/2022.10.17.512424v1)
@article {Sinha2022.10.17.512424,
	author = {Sanju Sinha and Neelam Sinha and Eytan Ruppin},
	title = {Deep characterization of cancer drugs mechanism of action by integrating large-scale genetic and drug screens},
	elocation-id = {2022.10.17.512424},
	year = {2022},
	doi = {10.1101/2022.10.17.512424},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2022/10/19/2022.10.17.512424},
	eprint = {https://www.biorxiv.org/content/early/2022/10/19/2022.10.17.512424.full.pdf},
	journal = {bioRxiv}
}
