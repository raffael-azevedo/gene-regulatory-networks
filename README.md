---
output:
  pdf_document: default
  html_document: default
---
# Gene Regulatory Networks Generation using RTN R/Bioconductor package.

This is an step-by-step on doing a gene regulatory analysis. Data used for this analysis are publicly available in Gene Expression Omnibus or Array Express. 

Links to data specified on the scripts: 
* [GSE4607](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE4607)
* [GSE26378](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE26378)
* [GSE13904](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE13904)
* [GSE13205](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE13205)
* [GSE21942](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE21942)
* [GSE48780](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48780)

For package installing, there is two sources: The Comprehensive R Archive Network - [CRAN](https://cran.r-project.org/) - and [Bioconductor](https://bioconductor.org). To install CRAN packages, on the R console, just type:

```{r}
install.packages("package")
```

For Bioconductor packages, the first installation requires the package BiocManager (from CRAN). After installing `BiocManager`, you can install Bioconductor packages using the following command on R console:

```{r}
BiocManager::install("package")
```
The parameter `package` can be substituted by an character vector (`c()`) containing all packages to install. 

R Packages used in this analysis:
1. From CRAN
* BiocManager
* data.table
* classInt 
* RColorBrewer 
* ggplot2 
* dplyr
* Honorable Mention: tidyverse

2. From Bioconductor
* affy
* limma
* biomaRt
* Fletcher2013b 
* RTN
* RedeR 
* enrichR

Info: `data.table`, `dplyr`, and  `ggplot2` are part of the `tidyverse`. You are encouraged to install the entire plethora of packages from the Tidyverse and learn to work with them. They make R data analysis easier than it is with `base`.

More details will be added here.