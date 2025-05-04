# DeepIP (R version) 
DeepIP is an R package to train and test a deep learning based model for predicting internal priming artifacts from DNA sequences. For complete and detailed information about **DeepIP**, please refer to the [DeepIP repository](https://github.com/APAexplorer/DeepIP)

About
====================
DeepIp is a component of [PolyAseqTrap](https://github.com/APAexplorer/PolyAseqTrap), a deep learning model designed to remove internal priming artifacts in polyA site identification. Inspired by the **DeepPASTA** model ([Arefeen, et al., 2019](https://doi.org/10.1093/bioinformatics/btz283)) that predicts polyA sites from DNA sequences, we designed a deep learning model called DeepIP to predict internal priming artifacts from A-rich polyA sites. DeepIP utilizes both convolutional neural network (CNN) and recurrent neural network (RNN). CNN extracts features from sequences, and RNN is used to combine the extracted feature effects for predicting internal priming artifacts. The corresponding DeepIP scripts can also be found in the [PolyAseqTrap](https://github.com/APAexplorer/PolyAseqTrap) GitHub repository. The schematic diagram of DeepIP is shown below.

<img src="https://github.com/APAexplorer/PolyAseqTrap/blob/main/img/DeepIP_schema.png" alt="schema" style="width:90%;"/>

Installing DeepIP (R version)
=============
Mandatory 
---------

* R (>=3.5.3). [4.3.1 ](https://www.r-project.org/) is recommended.

Required R Packages
---------
* [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html), [IRanges](https://bioconductor.org/packages/release/bioc/html/IRanges.html), [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html),[GenomicFeatures](https://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html), [seqinr](https://cran.r-project.org/web/packages/seqinr/index.html), [tidyr](https://cran.r-project.org/web/packages/tidyr/index.html), [movAPA](https://github.com/BMILAB/movAPA), [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)
  
Installation
---------
* Install the R package using the following commands on the R console:
```
install.packages("devtools")
library(devtools)
install_github("APAexplorer/DeepIP_R")
library(DeepIP)
browseVignettes(DeepIP)

##or you can download ZIP, and then unzip
devtools::install_local("your_path_of_DeepIP-master.zip", build_vignettes = TRUE)
```

Application examples
=============
Prepare data for DeepIP model training and testing
---------

This documentation explains how to generate training and test FASTA files for deep learning models. For demonstration purposes, we used **Arabidopsis** as the model species.

All output files for this demo can be found [here](https://github.com/APAexplorer/DeepIP_R/tree/main/demo_data_ath).

For complete details, please refer to the vignette ([PDF](https://github.com/APAexplorer/DeepIP_R/blob/main/doc/DeepIP_data_preparation.pdf) or HTML).

**Note**: We currently provide pre-trained models for **human**, **mouse**, and **Arabidopsis** species (see **[Training Model](https://github.com/APAexplorer/DeepIP/tree/main/training_model)**).





  
