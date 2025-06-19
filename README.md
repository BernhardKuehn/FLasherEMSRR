# FLasherEMSRR  <img src="hexsticker/hexsticker.png" align="right" width="170" />

[![R-CMD-check](https://github.com/BernhardKuehn/FLasherEMSRR/workflows/R-CMD-check/badge.svg)](https://github.com/BernhardKuehn/FLasherEMSRR/actions)
[![License](https://flr-project.org/img/eupl_1.1.svg)](https://joinup.ec.europa.eu/licence/european-union-public-licence-version-11-or-later-eupl)
[![Lifecycle:Experimental](https://img.shields.io/badge/Lifecycle-Experimental-339999)](<Redirect-URL>)

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

## Usage

> [!CAUTION]
> IMPORTANT: When using this package, also in combination with the [mseEMSRR-package](https://github.com/BernhardKuehn/mseEMSRR), make sure NOT to load the original 'FLasher' package within the same session, as this will cause errors due to overlapping function and class definitions.

To incorporate an explicit environmentally-mediated stock recruitment relationship one needs basically two things: 
- a fitted model object & an accompanying 'predict' function working with SSB and an environmental covariate
- an environmental covariate for the projection period
  
The way it is currently implemented is putting the covariate in the 'params' slot in the SR-object and having the 'predict' call working with the arguments 'ssb' and 'params'. A simple example showing some workable code in a vignette is still under construction and will be available in the upcoming months...

## Documentation
- [Forecasting on the Medium Term for advice using FLasher](https://flr-project.org/doc/Forecasting_on_the_Medium_Term_for_advice_using_FLasher.html)
- [Help pages](http://flr-project.org/FLasher)

## License
Copyright (c) 2016-2022 European Union. Released under the [EUPL v1.2](https://eupl.eu/1.2/en/).

## Contact
For info on the main-functionality of the 'FLasher' package:
- Author(s) of 'FLasher': Finlay Scott and Iago Mosqueira (EC-JRC).
- Maintainer of 'FlasherEMSRR': Bernhard Kuehn  
For the small extension allowing the incorporation of flexible SRRs:
- Contact: Bernhard Kuehn 

