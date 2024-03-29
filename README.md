
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Last-changedate](https://img.shields.io/badge/last%20change-2019--12--12-brightgreen.svg)](https://github.com/jmcascalheira/ScaledPiecesVB/commits/master)
[![minimal R
version](https://img.shields.io/badge/R%3E%3D-3.2.4-brightgreen.svg)](https://cran.r-project.org/)
[![Licence](https://img.shields.io/github/license/mashape/apistatus.svg)](http://choosealicense.com/licenses/mit/)
[![ORCiD](https://img.shields.io/badge/ORCiD-/0000--0003--0321--8892-green.svg)](http://orcid.org/0000-0003-0321-8892)

# LGMTechnology

This repository contains the data and code for the paper:

> Cascalheira, J. (2019). *Territoriality and the organization of
> technology during the Last Glacial Maximum in southwestern Europe*.
> PLOS ONE 14(12): e0225828
> <https://doi.org/10.1371/journal.pone.0225828>

The pre-print is online here:

> Cascalheira, J. (2019). *Territoriality and the organization of
> technology during the Last Glacial Maximum in southwestern Europe*.
> Accessed 12 Dec 2019. Online at
> <https://doi.org/10.31235/osf.io/fgxtu>

### How to cite

Please cite this compendium as:

> Cascalheira, J. (2019). *Compendium of R code and data for
> Territoriality and the organization of technology during the Last
> Glacial Maximum in southwestern Europe*. Accessed 12 Dec 2019. Online
> at <https://doi.org/10.17605/OSF.IO/YD2VE>

### How to download or install

You can download the compendium as a zip from from this URL:
<https://github.com/jmcascalheira/LGMtechnology/archive/master.zip>

This repository is organized as an R package using rrtools by Ben
Marwick, which can be installed from github with:

``` r
# install.packages("devtools")
devtools::install_github("benmarwick/rrtools")
```

There are no actual R functions in this package - all the R code is in
the Rmd file. We simply used the R package structure to help manage
dependencies, to take advantage of continuous integration for automated
code testing, and so I didn’t have to think too much about how to
organise the files.

To download the package source as you see it on GitHub, for offline
browsing, use this line at the shell prompt (assuming you have Git
installed on your computer):

``` r
git clone https://github.com/jmcascalheira/LGMtechnology.git
```

Once the download is complete, open the `lgmtechnolgoy.Rproj` in RStudio
to begin working with the package and compendium files.

The package has a number of dependencies on other R packages, and
programs outside of R. Installing these can be time-consuming and
complicated, so we’ve included a packrat directory, which contains the
source code for all the packages we depend on. If all works well, these
will be installed on your computer when you open `lgmtechnology.Rproj`
in RStudio.

### Licenses

**Text and figures :**
[CC-BY-4.0](http://creativecommons.org/licenses/by/4.0/)

**Code :** See the [DESCRIPTION](DESCRIPTION) file

**Data :** [CC-0](http://creativecommons.org/publicdomain/zero/1.0/)
attribution requested in reuse

### Contributions

I welcome contributions from everyone. Before you get started, please
see my [contributor guidelines](CONTRIBUTING.md). Please note that this
project is released with a [Contributor Code of Conduct](CONDUCT.md). By
participating in this project you agree to abide by its terms.

### Dependencies:

See the colophon section of the docx file in `analysis/paper` for a full
list of the packages that this project depends on.

### Contact:

João Cascalheira, Post-Doc Researcher, ICArEHB University of Algarve,
Campus de Gambelas 8005-139 Faro PORTUGAL e. <jmcascalheira@ualg.pt> w.
<http://www.icarehb.com/cascalheira>
