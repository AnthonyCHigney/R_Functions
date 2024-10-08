
#This is a function to run the diff-in-diff estimator "Difference-in-Differences Estimators 
#with Continuous Treatments and no Stayers, Clément de Chaisemartin, Xavier D'Haultfoeuille and 
#Gonzalo Vazquez-Bare, (2024)>>) where from period 1 to 2 all the units change treatment status (NO STAYERS)."
#It also includes a slight extension of covariates.
#these covariates must be pre-treatment (i.e. period 1) 
#The model is a saturated parametric model
#It is assumed there are only two periods

#DEPENDECIES
#packages dplyr and fixest

#ARGUMENTS

#MODEL -
#model = a fitted fixest model


#DATA - 
#The data must be arranged such that there is only one row per panel ID
#So for example individual one is only one row. They have the change in Y as the Y variable in the feols model
#they have a D1 variable (level of treatment in period 1). and they have a change in treatment variable and the 
#period 1 X variables

#*****IMPORTANT**********************
#*DO not use the | symbol in your fixest model. I.e. the fixed effects separator or pipe.
#this will mess up the function
#so use it without | pipe, like this:
#feols(data = df, y ~ d1 + x1 + d1:x1 + deltaD_XX + I(deltaD_XX^2) +  d1:deltaD_XX +   d1:deltaD_X:x1
#deltaD_XX:x1 + I(deltaD_XX^2):x1)
#just takes a little longer to compute but still fast


#deltaY_XX -
#what variable name is the change in your Y variable

#deltaD_XX -
#What variable name is the change in your treatment variable

#The function will output a point estimate (theta-hat). This is the estimated WAOSS 
#It will output a standard error with the influence function 
#and it outputs 95% confidence interval upper and lower bounds











did_continuous_nostayers_covary_saturated <- function(model, data, deltaY_XX, deltaD_XX) {
  

  #make the model matrix
  fitted_trend_matrix <- model.matrix(model)
  
  #adjust model matrix so that all variables which include deltaD_XX are not included in counterfactual

  variable_name <- deltaD_XX
  
  
  # Find columns that contain the unneeded variable name string
  cols_to_zero <- grepl(variable_name, colnames(fitted_trend_matrix))
  
  #Columns that contain the deltaD_XX variable, as these are not used in counterfactual g function, are set to zero
  fitted_trend_matrix[, cols_to_zero] <- 0

  
  #now get residuals for variance calcs later on
  data$eps_XX <- model$residuals
  

  # Generate signdeltaD_XX variable (sign of change in treatment)
  data$signdeltaD_XX <- sign(data[[deltaD_XX]])

  
  #------------------------------------------------------------
  #COMPUTE POINT ESTIMATE
  #now we make the numerator of the WAOSS point estimate
  #first we get the g(d1,0) fitted counterfactual trend for each observation
  gfunction_counterfactual <- fitted_trend_matrix %*% as.matrix(model$coeftable[1])
  
  #now change in Y minus our fitted counterfactual trend
  data$num_theta_XX <- data[[deltaY_XX]] - gfunction_counterfactual
  
  #multiplied by sign of change in treatment for each row
  data$num_theta_XX <- data$signdeltaD_XX*data$num_theta_XX
  
  
  #now sum - numerator done
  num_theta <- sum(data$num_theta_XX, na.rm = TRUE)
  #denominator is just sum of absolute change in treatment
  denom_theta <- sum((abs(data[[deltaD_XX]])), na.rm = TRUE)
  
  

  #compute theta finally 
  theta_hat <- num_theta/denom_theta
  #Point estimate done
  
  
  
  
  #------------------------------------------------------------
  #COMPUTE VARIANCE
  #now we compute the standard errors as per the paper with an influence function
  
  #this function is made up of
  #our already calculated denominator of theta (1/E[absolute change in treatment])
  #the sample size
  #Our already calculated num_theta - sign of change in treatment*(change in Y - counterfactual)
  # expected value of [sign of change in treatment multiplied by derivative of counterfactual estimated function with respect to parameters
  #influence function
  #theta point estimate multiplied by absolute change in treatment
  
  
  #get sample size
  N_XX <- nrow(data)
  
  #get the absolute change in treatment
  data$absdeltaD_XX = abs(data[[deltaD_XX]])
  
  
  #the derivative of g function counterfactual with respect to the parameters is just our "trend_matrix"
  #we now just get expected value of each (signed) derivative by summing and dividing by N 
  #and multiply by sign before that
  
  
  expected_value_deriv <- as.matrix(colSums(fitted_trend_matrix*data$signdeltaD_XX)/N_XX)
  
  
  #now just need influence function
  #first make inverse XX
  #X from model
  X <- model.matrix(model)
  #crossprod
  Xcross_XX <- crossprod(X)
  
  #invert matrix
  Q_XX <- chol2inv(chol(Xcross_XX))
 
  #X as we are doing (XX)^-1 X'
  u_XX <- model.matrix(model)
  
  #matrix IF = Q*ut*predicted
  Qut_XX <- t(Q_XX %*% t(u_XX))
  
  
  
  
  #make df for the loop
  Qut_XX_df <- as.data.frame(Qut_XX)
  
  #-----------------------------------------------------------
  
  #empty matrix for IF_XX to fill
  IF_XX <-  matrix(0, nrow = nrow(data), ncol(fitted_trend_matrix))
  
  # Loop through each row of Qut_XX to get influence function
  for (i in 1:ncol(fitted_trend_matrix)) {
    #Calculate IF_i_XX
    IF_XX[,i] <- (Qut_XX_df[,i] * data$eps_XX * N_XX)
  }
  
  # Store IF_XX in a data frame
  IF_XX_df <- as.data.frame(IF_XX) 
  
  
  #Generate Expected value of derivative * influence function
  Expderiv_IF_XX <- t(expected_value_deriv) %*% t(as.matrix(IF_XX_df))
  E_w_denom_theta_XX <- denom_theta / N_XX
  
  #Generate Phi_XX from paper
  Phi_XX <- (data$num_theta_XX - t(Expderiv_IF_XX) - 
               theta_hat * data$absdeltaD_XX) / E_w_denom_theta_XX
  
  
  
  
  
  sum_phiXX = sum(Phi_XX, na.rm = TRUE)
  #expected value
  E_w_Phi_XX = sum_phiXX/N_XX
  
  
  
  
  
  
  #Now we have everything to calculate variance as standard
  #demeaning Phi_XX, taking the square 
  Phi_XX <- Phi_XX - E_w_Phi_XX
  Phi_XX <- Phi_XX * Phi_XX
  sum_Phi_XX <- sum(Phi_XX, na.rm = TRUE)
  
  #Calculate variance of theta_XX
  var_theta_XX <- mean(Phi_XX) / N_XX
  # Calculate standard error of theta_XX
  sd_theta_XX <- sqrt(var_theta_XX)
  # Calculate confidence intervals for theta_XX
  LB_XX <- theta_hat - 1.96 * sd_theta_XX
  UB_XX <- theta_hat + 1.96 * sd_theta_XX
  
  #return values
  return(list(theta_hat = theta_hat,
              se =sd_theta_XX, 
              LB = LB_XX, UB = UB_XX))
  
  
  
  
  
}





