Package "rutilstimflutre"
=========================

This directory contains the "rutilstimflutre" package for the R programming language.
This package contains Timothee Flutre's personal R code which can nevertheless be useful to others.
The development is funded by the Institut National de la Recherche Agronomique (INRA).

The copyright is owned by the INRA.
The code is available under the GNU Public License (version 3 and later).
See the COPYING file for usage permissions.
It is also available under the GNU Affero General Public License (version 3 and later), especially important when the software is run over a network.

The content of this directory is versioned using git, the central repository being hosted on [GitHub](https://github.com/timflutre/rutilstimflutre).
You can report issues on GitHub directly, but remember to copy-paste the output of sessionInfo().

I have invested a lot of time and effort in creating this package, please cite it when using it for data analysis:
```
R> citation("rutilstimflutre")
```
See also citation() for citing R itself.

**For users**, the easiest is to directly install the package from GitHub:
```
R> library(devtools); install_github("timflutre/rutilstimflutre")
```

In the case where the package is officially released as a "bundled package"
(e.g. as on CRAN), you can install it directly from inside R via:
```
R> install.packages("path/to/rutilstimflutre_<version>.tar.gz")
```

**For advanced users**, when retrieving the bundled package (that is, as a tar.gz), you can install it from the command-line via:
```
$ R CMD INSTALL rutilstimflutre_<version>.tar.gz
```

Once the package is installed, start using it:
```
R> library(rutilstimflutre)
R> help(package=rutilstimflutre)
```

**For developpers**, when editing the content of this repo, increment the version of the package in DESCRIPTION and execute the following commands:
```
$ Rscript -e 'library(devtools); devtools::document()'
$ R CMD build rutilstimflutre
$ R CMD check rutilstimflutre_<version>.tar.gz
```
It may be necessary to `export _R_CHECK_FORCE_SUGGESTS_=0` before `R CMD check`.

When developping, here are a few useful tips:
```
R> detach("package:rutilstimflutre", unload=TRUE)
R> testthat::test_file("tests/testthat/test_quantgen.R")
```
