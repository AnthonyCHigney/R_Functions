
#This is the R conversion of the stata ado file from the replication package
#of "Difference-in-Differences Estimators with Continuous Treatments and no Stayers, Clément de Chaisemartin, Xavier D'Haultfoeuille and 
#Gonzalo Vazquez-Bare, (2023)>>) where from period 1 to 2 all the units change treatment status (NO STAYERS)."




# Define the function for estimating WAOSS
did_continuous_nostayers <- function(Y_XX, ID_XX, T_XX, D_XX, weight_XX = NULL, data = NULL) {

  
  
 
  # If data is provided, extract variable names from the dataset
  if (!is.null(data)) {
    
    data = data.frame(Y_XX = data[[substituteDirect(Y_XX)]] ,# Extract values for Y_XX from the dataset
    ID_XX = data[[substituteDirect(ID_XX)]],  # Extract values for ID_XX from the dataset
    T_XX = data[[substituteDirect(T_XX)]] , # Extract values for T_XX from the dataset
    D_XX = data[[substituteDirect(D_XX)]],    # Extract values for D_XX from the dataset
    weight_XX = if(!is.null(weight_XX)) {  data[[substituteDirect(weight_XX)]]
    }
    else {
                       rep(1, length(Y_XX))
                        } )
    #equal weights if no weight specified
  }



  
  
  # Sort the data by ID_XX and T_XX
  data <- data %>% arrange(ID_XX, T_XX)
  
  # Generate deltaD_XX representing the difference between two periods of D_XX
  data <- data %>% 
    group_by(ID_XX) %>% 
    mutate(deltaD_XX = D_XX - dplyr::lag(D_XX)) %>%
    ungroup()
  

  
  # Generate deltaY_XX representing the difference between two periods of Y_XX
  data <- data %>% 
    group_by(ID_XX) %>% 
    mutate(deltaY_XX = Y_XX - dplyr::lag(Y_XX)) %>%
    ungroup()

  
#get D1_XX values
data <- data %>%
    group_by(ID_XX) %>%
    mutate(D1_XX = first(D_XX))

  

#we have all data on the rows for T2 now so can drop T1
data  <- data[data$T_XX == max(data$T_XX),]

#make scalars for Sample size and weight sum
W_XX = sum(data$weight_XX, na.rm = TRUE)
if (is.null(data)) {
  N_XX <- length(data$Y_XX)
} else {
  N_XX <- nrow(data)
}


#change these if you want to change the weighting system (analytical vs freq e.g.)

data$weights = data$weight_XX 


#NOW ESTIMATE LAMBDA PARAMETERS
#These are the fitted parameteres for our modelled counterfactual
#it uses quasi stayers

regLambda <- lm(data = data, formula  = deltaY_XX ~ D1_XX + deltaD_XX +
               D1_XX:deltaD_XX + I(deltaD_XX^2), weights = weights)


 
#now get residuals for variance calcs later on
data$eps_XX <- regLambda$residuals


#Now we calculate theta, the average of slopes of the treated

#first get sign of chnage in treatment deltaD_XX
data$signdeltaD_XX <- sign(data$deltaD_XX)

#get the individual values before averaging
data$num_theta_XX <- data$signdeltaD_XX*(data$deltaY_XX - regLambda$coefficients[["(Intercept)"]] - 
                        regLambda$coefficients[["D1_XX"]]*data$D1_XX)


#if weights used then we weight here otherwise it's just multiply by 1

#now sum

num_theta <- sum(data$num_theta_XX*data$weights, na.rm = TRUE)
  
denom_theta <- sum((abs(data$deltaD_XX)*data$weights), na.rm = TRUE)

#compute theta finally 
#------------------------------------------------------------

theta_hat <- num_theta/denom_theta
  
#------------------------------------------------------------
  
  
#==================================================================
#Now we compute the variance



#generate sqrt(w)*u
data$X2_XX <- sqrt(data$weights) * data$D1_XX * data$deltaD_XX
data$X3_XX <- sqrt(data$weights) * data$deltaD_XX * data$deltaD_XX
data$sqrt_weight_XX <- sqrt(data$weights)
data$D1_wXX <- sqrt(data$weights) * data$D1_XX
data$deltaD_wXX <- sqrt(data$weights) * data$deltaD_XX

# Create a matrix 'Xcross_XX' without a constant term
data_matrix <- cbind(data$sqrt_weight_XX, data$D1_wXX, 
                     data$deltaD_wXX, data$X2_XX, data$X3_XX)

Xcross_XX <- crossprod(data_matrix)






#invert
Q_XX <- solve(Xcross_XX)

#Generate w*u
data$X2_wXX <- data$weights * data$D1_XX * data$deltaD_XX
data$X3_wXX <- data$weights * data$deltaD_XX * data$deltaD_XX
data$D1_w2XX <- data$weights* data$D1_XX
data$deltaD_w2XX <- data$weights * data$deltaD_XX

#Create a matrix 'u_XX'
u_XX <- cbind(data$weights, data$D1_w2XX, 
              data$deltaD_w2XX, data$X2_wXX, 
              data$X3_wXX)

#Transpose 'u_XX'
ut_XX <- t(u_XX)


#matrix IF = Q*ut*predicted (IF = influence function)

Qut_XX <- Q_XX %*% ut_XX

Qut_XX <- t(Qut_XX)


Qut_XX_df <- as.data.frame(Qut_XX)


IF_XX <-  matrix(0, nrow = 2342, ncol = 5)

# Loop through each row of Qut_XX
for (i in 1:5) {
  #Calculate IF_i_XX
 IF_XX[,i] <- (Qut_XX_df[,i] * data$eps_XX * N_XX)
}

# Store IF_XX in a data frame
IF_XX_df <- as.data.frame(IF_XX)


#Compute the partial derivative : sgn(deltaD)*deriv(g)(D1,0)

data$derivg_alph0_XX <- 1
data$derivg_alph1_XX <- data$D1_XX
data$deriv0_XX <- data$signdeltaD_XX * data$derivg_alph0_XX * data$weights
sum_deriv0_XX <- sum(data$deriv0_XX, na.rm = TRUE)
exp_deriv0_XX <- sum_deriv0_XX / W_XX


data$deriv1_XX <- data$signdeltaD_XX * data$derivg_alph1_XX * data$weights
sum_deriv1_XX <- sum(data$deriv1_XX, na.rm = TRUE)
exp_deriv1_XX <- sum_deriv1_XX / W_XX



#Generate Expderiv_IF_XX
Expderiv_IF_XX <- exp_deriv0_XX * IF_XX[,1] +
  exp_deriv1_XX * IF_XX[,2]
E_w_denom_theta_XX <- denom_theta / W_XX

#Generate Phi_XX
data$absdeltaD_XX = abs(data$deltaD_XX)

Phi_XX <- (data$num_theta_XX - Expderiv_IF_XX - 
                  theta_hat * data$absdeltaD_XX) / E_w_denom_theta_XX

Phi_XX1 <- (data$num_theta_XX - Expderiv_IF_XX - 
             theta_hat * data$absdeltaD_XX) / E_w_denom_theta_XX


Phi_wXX <- data$weights*Phi_XX

sum_phiwXX = sum(Phi_wXX, na.rm = TRUE)
E_w_Phi_XX = sum_phiwXX/W_XX




#demeaning Phi_XX, taking the square and weighting it //weight option
Phi_XX <- Phi_XX - E_w_Phi_XX
Phi_XX <- Phi_XX * Phi_XX
Phi_XX <- Phi_XX*data$weights
sum_Phi_XX <- sum(Phi_XX, na.rm = TRUE)

#Calculate variance of theta_XX
var_theta_XX <- mean(Phi_XX) / W_XX
#Calculate standard deviation of theta_XX
sd_theta_XX <- sqrt(var_theta_XX)
#Calculate confidence intervals for theta_XX
LB_XX <- theta_hat - 1.96 * sd_theta_XX
UB_XX <- theta_hat + 1.96 * sd_theta_XX



#output results
return(list(theta_hat = theta_hat,
          se =sd_theta_XX, 
          LB = LB_XX, UB = UB_XX))




  
  
#end of function
}