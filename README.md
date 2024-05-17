SLIDE :Significant Latent Factor Interaction Discovery and Exploration across biological domains
============================================
 
 
 Here, we present SLIDE, a novel data-distribution-free approach to analyze high-dimensional multiomic datasets and uncover latent factors that drive the outcome of interest.
 - SLIDE makes no assumptions regarding the distribution of the underlying data as it significantly builds on a unique latent-factor regression framework.
 - It takes into account an extremely large search space of relationships to converge on a very small subset of biologically relevant and actionable latent factors.
 - Critically, SLIDE incorporates both linear and nonlinear relationships, including complex hierarchical structures.
 - The discovery of SLIDE is also coupled to rigorous false discovery rate (FDR) control via our unique analytical framework that creatively adapts ultramodern methods for FDR control

The link for the paper:[website](https://www.nature.com/articles/s41592-024-02175-z)





## Contents

- [Installation Guide](#installation-guide)
- [Demo](#demo)
- [Results](#results)
- [Hardware Requirements](#Hardware-Requirements)



## Installation Procedure
# SLIDE
First install `SLIDE` in R by running the following code from a  R command terminal:



```library(devtools)```   
```devtools::install_github("jishnu-lab/SLIDE")```
Requires R ≥ 4.2


## Running the software

For running the software please refer to `vignette` folder. 




## Hardware Requirements
For minimal performance of the ```SLIDE```,  a computer with about 8 GB of RAM is required. For optimal performance, we recommend a computer with the following specs:

RAM: 16+ GB
CPU: 4+ cores, 3.3+ GHz/core






   
