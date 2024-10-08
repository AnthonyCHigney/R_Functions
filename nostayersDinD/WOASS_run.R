#Load required libraries
pack <- c("dplyr", "foreign", "fixest")

lapply(pack, FUN = library, character.only = TRUE)

#load the function 

source("WAOSS_covary_saturated_estimator.R")



#Load the dataset with foreign package
#data is from de Chaisemartin, Cl´ement, Xavier D’Haultfoeuille, and Gonzalo Vazquez-Bare.
#2024. “Difference-in-Differences Estimators with Continuous Treatments and NoStayers.”
#and orginally from:
#Deschênes, Olivier, and Michael Greenstone. 2012. “The economic impacts of climate
#change: evidence from agricultural output and random fluctuations in weather: reply.”
#American Economic Review, 102(7): 3761–3773.

df <- read.dta( file.path("DG_CORRECTED_SAMPLE_2012.dta"))




# Sort the data by ID and Time
data <- df[df$year>1992,] %>% arrange(fips, year)

# Generate deltaD_XX representing the difference between two periods of D_XX
data <- data %>% 
  group_by(fips) %>% 
  mutate(deltaD_XX = dd89 - dplyr::lag(dd89)) %>%
  ungroup()



# Generate deltaY_XX representing the difference between two periods of Y_XX
data <- data %>% 
  group_by(fips) %>% 
  mutate(deltaY_XX = y - dplyr::lag(y)) %>%
  ungroup()


#get D1 values
data <- data %>%
  group_by(fips) %>%
  mutate(d1 = first(dd89))


#we have all data on the rows for T2 now so can drop T1
data <- data[data$year == max(data$year),]



#-------------------------------------------------------------------------------
#fir model for reg without covariates


reg1 <- feols(data = data, fml = deltaY_XX ~ d1 + deltaD_XX +
                d1:deltaD_XX + I(deltaD_XX^2))



WAOSS1 <- did_continuous_nostayers_covary_saturated(model = reg1, data = data,
                                          deltaY_XX = "deltaY_XX", deltaD_XX = "deltaD_XX")


#Fit model with reg that has covariates
reg2 <- feols(data = data, fml = deltaY_XX ~ d1 + deltaD_XX + prcp  + d1:prcp + #this is the trend without treatment change
                d1:deltaD_XX + I(deltaD_XX^2)  + deltaD_XX:prcp + deltaD_XX:prcp:d1 + I(deltaD_XX^2):prcp )



WAOSS2 <- did_continuous_nostayers_covary_saturated(model = reg2, data = data, 
                                                deltaY_XX = "deltaY_XX", deltaD_XX = "deltaD_XX")
