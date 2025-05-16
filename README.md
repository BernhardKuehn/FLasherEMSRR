# FLasherEMSRR 

## This is a fork of the flr-package 'FLasher'! Currently this 'Readme' copies only what was written there. Will be updated if full functionality of the package is ready!


[![R-CMD-check](https://github.com/BernhardKuehn/FLasherEMSRR/workflows/R-CMD-check/badge.svg)](https://github.com/BernhardKuehn/FLasherEMSRR/actions)
[![License](https://flr-project.org/img/eupl_1.1.svg)](https://joinup.ec.europa.eu/licence/european-union-public-licence-version-11-or-later-eupl)
[![Codecov test coverage](https://codecov.io/gh/BernhardKuehn/FLasherEMSRR/branch/master/graph/badge.svg)](https://codecov.io/gh/BernhardKuehn/FLasherEMSRR?branch=master)

## Overview

Projection of future population and fishery dynamics is carried out for a given set of management targets. A system of equations is solved, using Automatic Differentation (AD), for the levels of effort by fishery/fleet that will result in the required abundances, catches or fishing mortalities. 

This fork explicitly allows for the incorporation of ANY (also non-parametric, flexible ones) stock-recruitment relationship (given that it is properly defined and can be called with a 'predict' function), allowing the integration of environmentally-driven stock recruitment relationships (EMSRRs) in the projection. 

## Installation
To install this package, start R and install it directly from this github repository by using:

```
remotes::install_github("BernhardKuehn/FLasherEMSRR")
```

**WARNING**: FLasher requires a 64 bit installation of R. Installation from source in R for Windows should be carried out using `--no-multiarch` for a 64 bit-only installation if both 32 and 64 bit R are available.

```
library(devtools)
options(devtools.install.args = "--no-multiarch")   
install_github("BernhardKuehn/FLasherEMSRR")
```

## Documentation
- [Forecasting on the Medium Term for advice using FLasher](https://flr-project.org/doc/Forecasting_on_the_Medium_Term_for_advice_using_FLasher.html)
- [Help pages](http://flr-project.org/FLasher)

## License
Copyright (c) 2016-2022 European Union. Released under the [EUPL v1.2](https://eupl.eu/1.2/en/).

## Contact
For info on the main-functionality of the 'FLasher' package:
- Author: Finlay Scott and Iago Mosqueira (EC-JRC).
- Maintainer: Iago Mosqueira <iago.mosqueira@wur.nl>
  
For the small extension allowing the incorporation of flexible SRRs:
- Contact: Bernhard Kuehn 

