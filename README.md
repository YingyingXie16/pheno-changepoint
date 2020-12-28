## Introduction
Timing of life cycle events (i.e. phenology) of vegetation is an important bio-indicator of climate change, and has substantial ecological functions and impacts. Identification of phenological transitions from satellite observations still has relatively large uncertainties due to limitation of raw data time series and methodology. I developed a new method using linear change points to estimate phenological transitions for deciduous forests in the central and eastern United States. This method increases estimation accuracy for senescence date in autumn, and provides quantification of estimation uncertainty. Other than four phenological dates, the method also estimates change rates of greenup, green-down in summer, and senesence in autumn, and the EVI values on four phenological transitio dates.

Here I'm sharing the R code used for change point estimation of phenological transitions using EVI time series from satellite imagery. The study is published in 2020 at the scientific journal, Remote Sensing of Environment.

## Code demostration
The code (demo_changepoint.r) is used to estimate phenological metrics using change point estimation method for 20,000 MODIS pixels sampled in the central and eastern US during 2000 and 2018. 

EVI values were calculated from MODIS products: twice-daily surface reflectance data from TERRA and
AQUA (MOD09GA and MYD09GA, Version 006), stored in netCDF files. Each file stores EVI value for all pixels in the study area at one time point.

Note: Initial estimation values were provided from logistic curve fitting results based on the same EVI time series data for each pixel.

### Reference
Xie Y, Wilson AM. (2020) Change Point Estimation of Deciduous Forest Land Surface Phenology. Remote Sensing and Environment, 240. doi: 10.1016/j.rse.2020.111698 
