#Load required libraries
pack <- c("dplyr", "foreign", "fixest", "here")

lapply(pack, FUN = library, character.only = TRUE)

#here

here::i_am("Code/WOASS_RUN.R")

#Load the dataset with foreign package
df <- read.dta(here("Data", "DG_CORRECTED_SAMPLE_2012.dta"))


# Run the WAOSS estimator

source("WAOSS_estimator.R")

#Run the estimator on the subsample restricting year in 1997 and 2002 as they do

#weighted (slightly different answer, possibly due to how weights handled in lm?)
result1 <- did_continuous_nostayers(Y_XX = "y",ID_XX = "fips", 
                                  T_XX =  "year",D_XX =  "dd89", 
                                   weight_XX = "fland", data = df[df$year> 1992,])

#unweighted (same answer as them)

result2 <- did_continuous_nostayers(Y_XX = "y",ID_XX = "fips", 
                                   T_XX =  "year",D_XX =  "dd89", 
                                   data = df[df$year> 1992,])
