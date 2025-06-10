
<!-- README.md is generated from README.Rmd. Please edit that file -->

![Genotools.](inst/rstudio/templates/genotools.png)

# genotools

The goal of the **genotools** package is to facilitate use of MoBa
genetic data that has undergone quality control (QC) in the
MoBaPsychGen\_v1 pipeline \[link to documentation\]. Currently, it
mostly includes functionality for making and using polygenic scores.

Please contact [Laurie
Hannigan](mailto:laurie.hannigan@bristol.ac.uk;laurie.hannigan@lds.no)
with bugs, feedback, or development ideas.

 

## Installation

**genotools** is built and runs entirely within the
[TSD](https://www.uio.no/english/services/it/research/sensitive-data/)
environment, in which MoBa data are accessed for analyses. As such, you
can’t install the package directly from github. Instead, you should
download the binary for the latest working version from
[here](https://osf.io/6g8bj/files/), import to a sensible location in
your project in TSD and install in R as follows, amending the path
appropriately:

You can install the latest working version as follows…

``` r
install.packages("//ess01/P471/data/durable/common/software/genotools_0.3.0.zip", 
                 repos=NULL,
                 type = "binary")
```

At present, any missing dependencies need to be installed manually from
the TSD CRAN copy, i.e.,

``` r
install.packages('dplyr',
                 repos = "file://tsd-evs/shared/R/cran")
```

You can also install from source if needed.

 

Find out how to use genotools by running `vignette("genotools")`.
