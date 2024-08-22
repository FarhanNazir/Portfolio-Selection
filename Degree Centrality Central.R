library('xdcclarge')
library('rmgarch')
library('tibble')
library('dplyr')
library('readxl')
library('quadprog')
library('igraph') 


#Function
#1. Calculate the logarithmic rate of return
rate_of_return = function(data){
  data$Close <- as.Date(data$Close)
  n <- nrow(data)
  u <- ncol(data)
  allrate <- signif(log(data[2:n, 2:u]/data[1:n-1, 2:u]),6)
  allmat <- as.matrix(allrate)
  return(allmat)
}

#2. Calculate the correlation matrix and Create MST
#2.1 Pearson Correlation

pearson_corr = function(allmat){
  allco <- signif(cor(allmat, method = "pearson"),6)
  alldm <- sqrt(2*(1-allco))
  allgraph <- graph.adjacency(alldm, mode="undirected", weighted = TRUE)
  aa <- simplify(allgraph)
  allMST <- minimum.spanning.tree(aa)
  return(allMST)
  
}



#2.2 DCC
dcc_corr = function(allmat){
  total_stock <- ncol(allmat)
  garch11.spec = ugarchspec(mean.model = list(armaOrder = c(0,0)), 
                            variance.model = list(garchOrder = c(1,1), 
                                                  model = "sGARCH"), 
                            distribution.model = "norm")
  mspec = multispec(replicate(garch11.spec, n = total_stock))
  fitlist = multifit(multispec = mspec, data = allmat,  solver = "hybrid" )
  ht <-sigma(fitlist)^2
  residuals <- residuals(fitlist)
  dcc<-dcc_estimation(ini.para=c(0.05,0.93) ,ht ,residuals)
  Rt<-signif(matrix(dcc$dcc_Rt,total_stock,total_stock),6)
  colnames(Rt) <-  colnames(allmat)[1:total_stock]
  rownames(Rt) <-  colnames(allmat)[1:total_stock]
  alldm <- sqrt(2*(1-Rt))
  allgraph <- graph.adjacency(alldm, mode="undirected", weighted = TRUE)
  aa <- simplify(allgraph)
  allMST <- minimum.spanning.tree(aa)
  return(allMST)
  
}


#3. Calculate the Degree Centrality Central
degree_central = function(allMST,n){
  dg<-degree(allMST, normalized= FALSE)
  stocks<- as.data.frame(dg)
  stocks<- stocks %>% rownames_to_column(var = "Name")
  stocks<- stocks[order(-stocks$dg),]
  stocks <- stocks %>% slice(1:n)
  stocks.name <- stocks$Name
  return(stocks.name)
}

degree_central_name = function(allMSt){
  dg<-degree(allMST, normalized= FALSE)
  stocks<- as.data.frame(dg)
  return(stocks)
}


#5. Create portfolio
#5.1 Equal weight portfolio 

equal_weight = function(return, port ){
  num_port <- length(port)
  equal_weight = rep(1,num_port)/num_port
  port.ret<-return[, port]
  
  
  # Calculate the portfolio return
  
  cov_matrix<-cov(port.ret)
  mean_returns <-colMeans(port.ret)
  
  
  # Get number of assets
  n <- length(mean_returns)
  
  # Define equal weight portfolio weights
  portfolio_weights <- rep(1/n, n)
  
  # Calculate portfolio return
  portfolio_return <- sum(mean_returns * portfolio_weights)
  
  # Calculate portfolio standard deviation
  portfolio_std_dev <- sqrt(t(portfolio_weights) %*% cov_matrix %*% portfolio_weights)
  
  # Calculate Sharpe ratio
  sharpe_ratio <- (portfolio_return) / portfolio_std_dev
  data <- list("sr" = sharpe_ratio , "std" = portfolio_std_dev, "r" = portfolio_return)
  return(data)
  
  
}



data_table = function(data1, data2, data3, data4, data5){
  df <- data.frame(r = c(data1$r,data2$r, data3$r,data4$r,data5$r ) , 
                   sd = c(data1$std,data2$std, data3$std,data4$std,data5$std ),
                   sr=c(data1$sr,data2$sr, data3$sr,data4$sr,data5$sr ))
  
  rownames(df) <- c("2016","2017","2018","2019","2020")
  
  print (df)
  
}


#5.3. Minimum Variance 
minimum_variance = function(return, stock){
  
  port <- stock
  port.ret<-return[, port]
  mean_returns <- colMeans(port.ret)
  cov_matrix<-cov(port.ret)
  
  # Define constraints
  n <- length(mean_returns)
  I <- diag(n)
  A <- matrix(1, ncol=n, nrow=1)
  b <- 1
  
  # Solve optimization problem
  sol <- solve.QP(Dmat = 2 * cov_matrix, dvec = -mean_returns, Amat = t(A), bvec = b, meq = 1)
  
  # Get portfolio weights
  portfolio_weights <- sol$solution
  
  # Calculate portfolio return
  portfolio_return <- sum(mean_returns * portfolio_weights)
  portfolio_std_dev <- sqrt(portfolio_weights %*% cov_matrix %*% portfolio_weights)
  sharpe_ratio <- (portfolio_return ) / portfolio_std_dev
  
  
  
  
  datalist <- list("r" = portfolio_return , "std" =  portfolio_std_dev, "sr" = sharpe_ratio )
  return(datalist)
}

# Calculate the stock return for the each year
Data_2016 <- read_excel("./Data/Stocks/data 2016.xlsx")
Data_2017 <- read_excel("./Data/Stocks/data 2017.xlsx")
Data_2018 <- read_excel("./Data/Stocks/data 2018.xlsx")
Data_2019 <- read_excel("./Data/Stocks/data 2019.xlsx")
Data_2020 <- read_excel("./Data/Stocks/data 2020.xlsx")
Data_all  <- read_excel("./Data/Stocks/data all.xlsx")


return16 <-rate_of_return(Data_2016)
return17 <-rate_of_return(Data_2017)
return18 <-rate_of_return(Data_2018)
return19 <-rate_of_return(Data_2019)
return20 <-rate_of_return(Data_2020)
returnAll <-rate_of_return(Data_all)

print(returnAll)

# 1.1 Calculate Correlation

calculate_dcc2016 <-  dcc_corr(return16)
calculate_dcc2017 <-  dcc_corr(return17)
calculate_dcc2018 <-  dcc_corr(return18)
calculate_dcc2019 <-  dcc_corr(return19)
calculate_dcc2020 <-  dcc_corr(return20)
calculate_dccAll <- dcc_corr(returnAll)

print(calculate_dcc2016)



calculate_pearson2016 <-  pearson_corr(return16)
calculate_pearson2017 <-  pearson_corr(return17)
calculate_pearson2018 <-  pearson_corr(return18)
calculate_pearson2019 <-  pearson_corr(return19)
calculate_pearson2020 <-  pearson_corr(return20)



#Degree Centrality Central

degree_central_name_2016 <- degree_central_name(calculate_dcc2016)


# DCC
# n = 70 
degree_central_70_dcc_2016 <-  degree_central(calculate_dcc2016, 70)
degree_central_70_dcc_2017 <-  degree_central(calculate_dcc2017, 70)
degree_central_70_dcc_2018 <-  degree_central(calculate_dcc2018, 70)
degree_central_70_dcc_2019 <-  degree_central(calculate_dcc2019, 70)
degree_central_70_dcc_2020 <-  degree_central(calculate_dcc2020, 70)


degree_central_70_equalweight_dcc_2016 <- equal_weight(return16, degree_central_70_dcc_2016)
degree_central_70_equalweight_dcc_2017 <- equal_weight(return17, degree_central_70_dcc_2017)
degree_central_70_equalweight_dcc_2018 <- equal_weight(return18, degree_central_70_dcc_2018)
degree_central_70_equalweight_dcc_2019 <- equal_weight(return19, degree_central_70_dcc_2019)
degree_central_70_equalweight_dcc_2020 <- equal_weight(return20, degree_central_70_dcc_2020)


degree_central_70_equalweight_dcc = data_table (degree_central_70_equalweight_dcc_2016,
                                                   degree_central_70_equalweight_dcc_2017,
                                                   degree_central_70_equalweight_dcc_2018,
                                                   degree_central_70_equalweight_dcc_2019,
                                                   degree_central_70_equalweight_dcc_2020)



degree_central_70_minvar_dcc_2016 <- minimum_variance(return16, degree_central_70_dcc_2016)
degree_central_70_minvar_dcc_2017 <- minimum_variance(return17, degree_central_70_dcc_2017)
degree_central_70_minvar_dcc_2018 <- minimum_variance(return18, degree_central_70_dcc_2018)
degree_central_70_minvar_dcc_2019 <- minimum_variance(return19, degree_central_70_dcc_2019)
degree_central_70_minvar_dcc_2020 <- minimum_variance(return20, degree_central_70_dcc_2020)


degree_central_70_minvar_dcc =  data_table (degree_central_70_minvar_dcc_2016,
                                               degree_central_70_minvar_dcc_2017,
                                               degree_central_70_minvar_dcc_2018,
                                               degree_central_70_minvar_dcc_2019,
                                               degree_central_70_minvar_dcc_2020)

# n = 80 

degree_central_80_dcc_2016 <-  degree_central(calculate_dcc2016, 80)
degree_central_80_dcc_2017 <-  degree_central(calculate_dcc2017, 80)
degree_central_80_dcc_2018 <-  degree_central(calculate_dcc2018, 80)
degree_central_80_dcc_2019 <-  degree_central(calculate_dcc2019, 80)
degree_central_80_dcc_2020 <-  degree_central(calculate_dcc2020, 80)


degree_central_80_equalweight_dcc_2016 <- equal_weight(return16, degree_central_80_dcc_2016)
degree_central_80_equalweight_dcc_2017 <- equal_weight(return17, degree_central_80_dcc_2017)
degree_central_80_equalweight_dcc_2018 <- equal_weight(return18, degree_central_80_dcc_2018)
degree_central_80_equalweight_dcc_2019 <- equal_weight(return19, degree_central_80_dcc_2019)
degree_central_80_equalweight_dcc_2020 <- equal_weight(return20, degree_central_80_dcc_2020)


degree_central_80_equalweight_dcc = data_table (degree_central_80_equalweight_dcc_2016,
                                                degree_central_80_equalweight_dcc_2017,
                                                degree_central_80_equalweight_dcc_2018,
                                                degree_central_80_equalweight_dcc_2019,
                                                degree_central_80_equalweight_dcc_2020)



degree_central_80_minvar_dcc_2016 <- minimum_variance(return16, degree_central_80_dcc_2016)
degree_central_80_minvar_dcc_2017 <- minimum_variance(return17, degree_central_80_dcc_2017)
degree_central_80_minvar_dcc_2018 <- minimum_variance(return18, degree_central_80_dcc_2018)
degree_central_80_minvar_dcc_2019 <- minimum_variance(return19, degree_central_80_dcc_2019)
degree_central_80_minvar_dcc_2020 <- minimum_variance(return20, degree_central_80_dcc_2020)


degree_central_80_minvar_dcc =  data_table (degree_central_80_minvar_dcc_2016,
                                            degree_central_80_minvar_dcc_2017,
                                            degree_central_80_minvar_dcc_2018,
                                            degree_central_80_minvar_dcc_2019,
                                            degree_central_80_minvar_dcc_2020)

# n = 90 

degree_central_90_dcc_2016 <-  degree_central(calculate_dcc2016, 90)
degree_central_90_dcc_2017 <-  degree_central(calculate_dcc2017, 90)
degree_central_90_dcc_2018 <-  degree_central(calculate_dcc2018, 90)
degree_central_90_dcc_2019 <-  degree_central(calculate_dcc2019, 90)
degree_central_90_dcc_2020 <-  degree_central(calculate_dcc2020, 90)


degree_central_90_equalweight_dcc_2016 <- equal_weight(return16, degree_central_90_dcc_2016)
degree_central_90_equalweight_dcc_2017 <- equal_weight(return17, degree_central_90_dcc_2017)
degree_central_90_equalweight_dcc_2018 <- equal_weight(return18, degree_central_90_dcc_2018)
degree_central_90_equalweight_dcc_2019 <- equal_weight(return19, degree_central_90_dcc_2019)
degree_central_90_equalweight_dcc_2020 <- equal_weight(return20, degree_central_90_dcc_2020)


degree_central_90_equalweight_dcc = data_table (degree_central_90_equalweight_dcc_2016,
                                                degree_central_90_equalweight_dcc_2017,
                                                degree_central_90_equalweight_dcc_2018,
                                                degree_central_90_equalweight_dcc_2019,
                                                degree_central_90_equalweight_dcc_2020)



degree_central_90_minvar_dcc_2016 <- minimum_variance(return16, degree_central_90_dcc_2016)
degree_central_90_minvar_dcc_2017 <- minimum_variance(return17, degree_central_90_dcc_2017)
degree_central_90_minvar_dcc_2018 <- minimum_variance(return18, degree_central_90_dcc_2018)
degree_central_90_minvar_dcc_2019 <- minimum_variance(return19, degree_central_90_dcc_2019)
degree_central_90_minvar_dcc_2020 <- minimum_variance(return20, degree_central_90_dcc_2020)


degree_central_90_minvar_dcc =  data_table (degree_central_90_minvar_dcc_2016,
                                            degree_central_90_minvar_dcc_2017,
                                            degree_central_90_minvar_dcc_2018,
                                            degree_central_90_minvar_dcc_2019,
                                            degree_central_90_minvar_dcc_2020)

# n = 100 

degree_central_100_dcc_2016 <-  degree_central(calculate_dcc2016, 100)
degree_central_100_dcc_2017 <-  degree_central(calculate_dcc2017, 100)
degree_central_100_dcc_2018 <-  degree_central(calculate_dcc2018, 100)
degree_central_100_dcc_2019 <-  degree_central(calculate_dcc2019, 100)
degree_central_100_dcc_2020 <-  degree_central(calculate_dcc2020, 100)


degree_central_100_equalweight_dcc_2016 <- equal_weight(return16, degree_central_100_dcc_2016)
degree_central_100_equalweight_dcc_2017 <- equal_weight(return17, degree_central_100_dcc_2017)
degree_central_100_equalweight_dcc_2018 <- equal_weight(return18, degree_central_100_dcc_2018)
degree_central_100_equalweight_dcc_2019 <- equal_weight(return19, degree_central_100_dcc_2019)
degree_central_100_equalweight_dcc_2020 <- equal_weight(return20, degree_central_100_dcc_2020)


degree_central_100_equalweight_dcc = data_table (degree_central_100_equalweight_dcc_2016,
                                                degree_central_100_equalweight_dcc_2017,
                                                degree_central_100_equalweight_dcc_2018,
                                                degree_central_100_equalweight_dcc_2019,
                                                degree_central_100_equalweight_dcc_2020)



degree_central_100_minvar_dcc_2016 <- minimum_variance(return16, degree_central_100_dcc_2016)
degree_central_100_minvar_dcc_2017 <- minimum_variance(return17, degree_central_100_dcc_2017)
degree_central_100_minvar_dcc_2018 <- minimum_variance(return18, degree_central_100_dcc_2018)
degree_central_100_minvar_dcc_2019 <- minimum_variance(return19, degree_central_100_dcc_2019)
degree_central_100_minvar_dcc_2020 <- minimum_variance(return20, degree_central_100_dcc_2020)


degree_central_100_minvar_dcc =  data_table (degree_central_100_minvar_dcc_2016,
                                            degree_central_100_minvar_dcc_2017,
                                            degree_central_100_minvar_dcc_2018,
                                            degree_central_100_minvar_dcc_2019,
                                            degree_central_100_minvar_dcc_2020)




# pearson
# n = 70 
degree_central_70_pearson_2016 <-  degree_central(calculate_pearson2016, 70)
degree_central_70_pearson_2017 <-  degree_central(calculate_pearson2017, 70)
degree_central_70_pearson_2018 <-  degree_central(calculate_pearson2018, 70)
degree_central_70_pearson_2019 <-  degree_central(calculate_pearson2019, 70)
degree_central_70_pearson_2020 <-  degree_central(calculate_pearson2020, 70)


degree_central_70_equalweight_pearson_2016 <- equal_weight(return16, degree_central_70_pearson_2016)
degree_central_70_equalweight_pearson_2017 <- equal_weight(return17, degree_central_70_pearson_2017)
degree_central_70_equalweight_pearson_2018 <- equal_weight(return18, degree_central_70_pearson_2018)
degree_central_70_equalweight_pearson_2019 <- equal_weight(return19, degree_central_70_pearson_2019)
degree_central_70_equalweight_pearson_2020 <- equal_weight(return20, degree_central_70_pearson_2020)


degree_central_70_equalweight_pearson = data_table (degree_central_70_equalweight_pearson_2016,
                                                 degree_central_70_equalweight_pearson_2017,
                                                 degree_central_70_equalweight_pearson_2018,
                                                 degree_central_70_equalweight_pearson_2019,
                                                 degree_central_70_equalweight_pearson_2020)



degree_central_70_minvar_pearson_2016 <- minimum_variance(return16, degree_central_70_pearson_2016)
degree_central_70_minvar_pearson_2017 <- minimum_variance(return17, degree_central_70_pearson_2017)
degree_central_70_minvar_pearson_2018 <- minimum_variance(return18, degree_central_70_pearson_2018)
degree_central_70_minvar_pearson_2019 <- minimum_variance(return19, degree_central_70_pearson_2019)
degree_central_70_minvar_pearson_2020 <- minimum_variance(return20, degree_central_70_pearson_2020)


degree_central_70_minvar_pearson =  data_table (degree_central_70_minvar_pearson_2016,
                                             degree_central_70_minvar_pearson_2017,
                                             degree_central_70_minvar_pearson_2018,
                                             degree_central_70_minvar_pearson_2019,
                                             degree_central_70_minvar_pearson_2020)

# n = 80 

degree_central_80_pearson_2016 <-  degree_central(calculate_pearson2016, 80)
degree_central_80_pearson_2017 <-  degree_central(calculate_pearson2017, 80)
degree_central_80_pearson_2018 <-  degree_central(calculate_pearson2018, 80)
degree_central_80_pearson_2019 <-  degree_central(calculate_pearson2019, 80)
degree_central_80_pearson_2020 <-  degree_central(calculate_pearson2020, 80)


degree_central_80_equalweight_pearson_2016 <- equal_weight(return16, degree_central_80_pearson_2016)
degree_central_80_equalweight_pearson_2017 <- equal_weight(return17, degree_central_80_pearson_2017)
degree_central_80_equalweight_pearson_2018 <- equal_weight(return18, degree_central_80_pearson_2018)
degree_central_80_equalweight_pearson_2019 <- equal_weight(return19, degree_central_80_pearson_2019)
degree_central_80_equalweight_pearson_2020 <- equal_weight(return20, degree_central_80_pearson_2020)


degree_central_80_equalweight_pearson = data_table (degree_central_80_equalweight_pearson_2016,
                                                    degree_central_80_equalweight_pearson_2017,
                                                    degree_central_80_equalweight_pearson_2018,
                                                    degree_central_80_equalweight_pearson_2019,
                                                    degree_central_80_equalweight_pearson_2020)



degree_central_80_minvar_pearson_2016 <- minimum_variance(return16, degree_central_80_pearson_2016)
degree_central_80_minvar_pearson_2017 <- minimum_variance(return17, degree_central_80_pearson_2017)
degree_central_80_minvar_pearson_2018 <- minimum_variance(return18, degree_central_80_pearson_2018)
degree_central_80_minvar_pearson_2019 <- minimum_variance(return19, degree_central_80_pearson_2019)
degree_central_80_minvar_pearson_2020 <- minimum_variance(return20, degree_central_80_pearson_2020)


degree_central_80_minvar_pearson =  data_table (degree_central_80_minvar_pearson_2016,
                                                degree_central_80_minvar_pearson_2017,
                                                degree_central_80_minvar_pearson_2018,
                                                degree_central_80_minvar_pearson_2019,
                                                degree_central_80_minvar_pearson_2020)

# n = 90 

degree_central_90_pearson_2016 <-  degree_central(calculate_pearson2016, 90)
degree_central_90_pearson_2017 <-  degree_central(calculate_pearson2017, 90)
degree_central_90_pearson_2018 <-  degree_central(calculate_pearson2018, 90)
degree_central_90_pearson_2019 <-  degree_central(calculate_pearson2019, 90)
degree_central_90_pearson_2020 <-  degree_central(calculate_pearson2020, 90)


degree_central_90_equalweight_pearson_2016 <- equal_weight(return16, degree_central_90_pearson_2016)
degree_central_90_equalweight_pearson_2017 <- equal_weight(return17, degree_central_90_pearson_2017)
degree_central_90_equalweight_pearson_2018 <- equal_weight(return18, degree_central_90_pearson_2018)
degree_central_90_equalweight_pearson_2019 <- equal_weight(return19, degree_central_90_pearson_2019)
degree_central_90_equalweight_pearson_2020 <- equal_weight(return20, degree_central_90_pearson_2020)


degree_central_90_equalweight_pearson = data_table (degree_central_90_equalweight_pearson_2016,
                                                    degree_central_90_equalweight_pearson_2017,
                                                    degree_central_90_equalweight_pearson_2018,
                                                    degree_central_90_equalweight_pearson_2019,
                                                    degree_central_90_equalweight_pearson_2020)



degree_central_90_minvar_pearson_2016 <- minimum_variance(return16, degree_central_90_pearson_2016)
degree_central_90_minvar_pearson_2017 <- minimum_variance(return17, degree_central_90_pearson_2017)
degree_central_90_minvar_pearson_2018 <- minimum_variance(return18, degree_central_90_pearson_2018)
degree_central_90_minvar_pearson_2019 <- minimum_variance(return19, degree_central_90_pearson_2019)
degree_central_90_minvar_pearson_2020 <- minimum_variance(return20, degree_central_90_pearson_2020)


degree_central_90_minvar_pearson =  data_table (degree_central_90_minvar_pearson_2016,
                                                degree_central_90_minvar_pearson_2017,
                                                degree_central_90_minvar_pearson_2018,
                                                degree_central_90_minvar_pearson_2019,
                                                degree_central_90_minvar_pearson_2020)

# n = 100 

degree_central_100_pearson_2016 <-  degree_central(calculate_pearson2016, 100)
degree_central_100_pearson_2017 <-  degree_central(calculate_pearson2017, 100)
degree_central_100_pearson_2018 <-  degree_central(calculate_pearson2018, 100)
degree_central_100_pearson_2019 <-  degree_central(calculate_pearson2019, 100)
degree_central_100_pearson_2020 <-  degree_central(calculate_pearson2020, 100)


degree_central_100_equalweight_pearson_2016 <- equal_weight(return16, degree_central_100_pearson_2016)
degree_central_100_equalweight_pearson_2017 <- equal_weight(return17, degree_central_100_pearson_2017)
degree_central_100_equalweight_pearson_2018 <- equal_weight(return18, degree_central_100_pearson_2018)
degree_central_100_equalweight_pearson_2019 <- equal_weight(return19, degree_central_100_pearson_2019)
degree_central_100_equalweight_pearson_2020 <- equal_weight(return20, degree_central_100_pearson_2020)


degree_central_100_equalweight_pearson = data_table (degree_central_100_equalweight_pearson_2016,
                                                    degree_central_100_equalweight_pearson_2017,
                                                    degree_central_100_equalweight_pearson_2018,
                                                    degree_central_100_equalweight_pearson_2019,
                                                    degree_central_100_equalweight_pearson_2020)



degree_central_100_minvar_pearson_2016 <- minimum_variance(return16, degree_central_100_pearson_2016)
degree_central_100_minvar_pearson_2017 <- minimum_variance(return17, degree_central_100_pearson_2017)
degree_central_100_minvar_pearson_2018 <- minimum_variance(return18, degree_central_100_pearson_2018)
degree_central_100_minvar_pearson_2019 <- minimum_variance(return19, degree_central_100_pearson_2019)
degree_central_100_minvar_pearson_2020 <- minimum_variance(return20, degree_central_100_pearson_2020)


degree_central_100_minvar_pearson =  data_table (degree_central_100_minvar_pearson_2016,
                                                degree_central_100_minvar_pearson_2017,
                                                degree_central_100_minvar_pearson_2018,
                                                degree_central_100_minvar_pearson_2019,
                                                degree_central_100_minvar_pearson_2020)




