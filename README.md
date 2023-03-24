# PAH-biomarker
Analysis pipeline for the PAH biomarker study by Sonnweber T et al. 

## Summary

A complete R pipeline for search of novel, clinically applicable signatures of pulmonary artherial hypertension mortality risk. 

<br>

<p align = "center"> 
<img src = "https://user-images.githubusercontent.com/80723424/227527668-fc2188c9-2c21-4298-893c-7e81c1f8d564.png" width = "80%">
</p>

<br>

## Terms of use

Please cite the repository, and the peer-reviewed publication of the analysis results when available. The raw data files will be made upon request to the study authors, [Dr. Thomas Sonnweber](mailto:thomas.sonnweber@tirol-kliniken.at) and [Prof. Judith Löffler-Ragg](mailto:judith.loeffler@i-med.ac.at).

## Usage

The following development packages are required to run the pipeline:

```r

devtools::install_github('PiotrTymoszuk/soucer') ## script sourcing
devtools::install_github('PiotrTymoszuk/ExDA') ## exploratory data analysis and staristical hypothesis testing
devtools::install_github('PiotrTymoszuk/clustTools') ## factor analysis and unsupervised clustering
devtools::install_github('PiotrTymoszuk/coxExtensions') ## fit statistics and quality control for Cox models
devtools::install_github('PiotrTymoszuk/figur') ## management of figures and tables in Rmd documents
devtools::install_github('PiotrTymoszuk/trafo') ## handling of tabular data

```

Source 'exec.R' to launch the entire pipeline:

```r

source('exec.R')

```

## Contact

The repository maintainer is [Piotr Tymoszuk](mailto:piotr.s.tymoszuk@gmail.com). Data requests should be addressed to [Dr. Thomas Sonnweber](mailto:thomas.sonnweber@tirol-kliniken.at) and [Prof. Judith Löffler-Ragg](mailto:judith.loeffler@i-med.ac.at).
