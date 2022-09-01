# flu_hosp
Fitting and forecasting of daily influenza hospital admission data using a flexible time-dependent transmission term

##########################################################

MAIN CODE 

flu_hosp_forecast.R

S[I]2RH fitting and forecasting code - using daily hospitalization data 

COMPILING THE FORTRAN CODE

The R code requires the compilation of fortran code

After downloading the repository please go to the gsrc subdirectory

%cd gsrc

and use the provided python script to compile the code

%./compile.py

If successul three dynamic (.so) libraries will be created:

detsirh.so, stosirh.so and mcmc.so

DATA

The code assumes a certain format for the data and an example is provided in the
data sub-directory

If you choose to change the location of the data file you will need to update the 'data_path'
variable in the flu_hosp_forecast.R code 

DETAILS ABOUT THE MODEL

The number of values that the time-dependent transmission term is set by the parameter 'nb'
Reasonable values are 1-4
For more on our model for R(t) and the fitting procedure see our publications:

Consistent pattern of epidemic slowing across many geographies led to longer, flatter initial waves of the COVID-19 pandemic
Michal Ben-Nun, Pete Riley, James Turtle, Steven Riley

Research Article | published 15 Aug 2022 PLOS Computational Biology

https://doi.org/10.1371/journal.pcbi.1010375

Forecasting national and regional influenza-like illness for the USA
Michal Ben-Nun, Pete Riley, James Turtle, David P. Bacon, Steven Riley

Research Article | published 23 May 2019 PLOS Computational Biology

https://doi.org/10.1371/journal.pcbi.1007013

Accurate influenza forecasts using type-specific incidence data for small geographic units
James Turtle, Pete Riley, Michal Ben-Nun, Steven Riley

Research Article | published 29 Jul 2021 PLOS Computational Biology

https://doi.org/10.1371/journal.pcbi.1009230

HELP and SUPPORT

For questions and/or help please email us: mbennun@predsci.com


##########################################################

