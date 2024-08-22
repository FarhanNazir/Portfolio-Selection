library('xdcclarge')
library('rmgarch')
library('tibble')
library('dplyr')


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


#3. Calculate the Betweeness Centrality Central
betweenness_central = function(allMST , n){
  dg<-betweenness(allMST, normalized= FALSE)
  stocks<- as.data.frame(dg)
  stocks<- stocks %>% rownames_to_column(var = "Name")
  stocks<- stocks[order(-stocks$dg),]
  stocks <- stocks %>% slice(1:n)
  stocks.name <- stocks$Name
  return(stocks.name)
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



# Calculate the stock return for the each year
Data_2016 <- read_excel("./Data/Stocks/data 2016.xlsx")
Data_2017 <- read_excel("./Data/Stocks/data 2017.xlsx")
Data_2018 <- read_excel("./Data/Stocks/data 2018.xlsx")
Data_2019 <- read_excel("./Data/Stocks/data 2019.xlsx")
Data_2020 <- read_excel("./Data/Stocks/data 2020.xlsx")

return16 <-rate_of_return(Data_2016)
return17 <-rate_of_return(Data_2017)
return18 <-rate_of_return(Data_2018)
return19 <-rate_of_return(Data_2019)
return20 <-rate_of_return(Data_2020)


# 1.1 Calculate Centrality

calculate_dcc2016 <-  dcc_corr(return16)
calculate_dcc2017 <-  dcc_corr(return17)
calculate_dcc2018 <-  dcc_corr(return18)
calculate_dcc2019 <-  dcc_corr(return19)
calculate_dcc2020 <-  dcc_corr(return20)


calculate_pearson2016 <-  pearson_corr(return16)
calculate_pearson2017 <-  pearson_corr(return17)
calculate_pearson2018 <-  pearson_corr(return18)
calculate_pearson2019 <-  pearson_corr(return19)
calculate_pearson2020 <-  pearson_corr(return20)



#Degree Centrality Central

# DCC
# n = 70 
betweenness_central_70_dcc_2016 <-  betweenness_central(calculate_dcc2016, 70)
betweenness_central_70_dcc_2017 <-  betweenness_central(calculate_dcc2017, 70)
betweenness_central_70_dcc_2018 <-  betweenness_central(calculate_dcc2018, 70)
betweenness_central_70_dcc_2019 <-  betweenness_central(calculate_dcc2019, 70)
betweenness_central_70_dcc_2020 <-  betweenness_central(calculate_dcc2020, 70)


betweenness_central_70_equalweight_dcc_2016 <- equal_weight(return16, betweenness_central_70_dcc_2016)
betweenness_central_70_equalweight_dcc_2017 <- equal_weight(return17, betweenness_central_70_dcc_2017)
betweenness_central_70_equalweight_dcc_2018 <- equal_weight(return18, betweenness_central_70_dcc_2018)
betweenness_central_70_equalweight_dcc_2019 <- equal_weight(return19, betweenness_central_70_dcc_2019)
betweenness_central_70_equalweight_dcc_2020 <- equal_weight(return20, betweenness_central_70_dcc_2020)


betweenness_central_70_equalweight_dcc = data_table (betweenness_central_70_equalweight_dcc_2016,
                                                        betweenness_central_70_equalweight_dcc_2017,
                                                        betweenness_central_70_equalweight_dcc_2018,
                                                        betweenness_central_70_equalweight_dcc_2019,
                                                        betweenness_central_70_equalweight_dcc_2020)



betweenness_central_70_minvar_dcc_2016 <- minimum_variance(return16, betweenness_central_70_dcc_2016)
betweenness_central_70_minvar_dcc_2017 <- minimum_variance(return17, betweenness_central_70_dcc_2017)
betweenness_central_70_minvar_dcc_2018 <- minimum_variance(return18, betweenness_central_70_dcc_2018)
betweenness_central_70_minvar_dcc_2019 <- minimum_variance(return19, betweenness_central_70_dcc_2019)
betweenness_central_70_minvar_dcc_2020 <- minimum_variance(return20, betweenness_central_70_dcc_2020)


betweenness_central_70_minvar_dcc =  data_table (betweenness_central_70_minvar_dcc_2016,
                                                    betweenness_central_70_minvar_dcc_2017,
                                                    betweenness_central_70_minvar_dcc_2018,
                                                    betweenness_central_70_minvar_dcc_2019,
                                                    betweenness_central_70_minvar_dcc_2020)


# n = 80 

betweenness_central_80_dcc_2016 <-  betweenness_central(calculate_dcc2016, 80)
betweenness_central_80_dcc_2017 <-  betweenness_central(calculate_dcc2017, 80)
betweenness_central_80_dcc_2018 <-  betweenness_central(calculate_dcc2018, 80)
betweenness_central_80_dcc_2019 <-  betweenness_central(calculate_dcc2019, 80)
betweenness_central_80_dcc_2020 <-  betweenness_central(calculate_dcc2020, 80)


betweenness_central_80_equalweight_dcc_2016 <- equal_weight(return16, betweenness_central_80_dcc_2016)
betweenness_central_80_equalweight_dcc_2017 <- equal_weight(return17, betweenness_central_80_dcc_2017)
betweenness_central_80_equalweight_dcc_2018 <- equal_weight(return18, betweenness_central_80_dcc_2018)
betweenness_central_80_equalweight_dcc_2019 <- equal_weight(return19, betweenness_central_80_dcc_2019)
betweenness_central_80_equalweight_dcc_2020 <- equal_weight(return20, betweenness_central_80_dcc_2020)


betweenness_central_80_equalweight_dcc = data_table (betweenness_central_80_equalweight_dcc_2016,
                                                     betweenness_central_80_equalweight_dcc_2017,
                                                     betweenness_central_80_equalweight_dcc_2018,
                                                     betweenness_central_80_equalweight_dcc_2019,
                                                     betweenness_central_80_equalweight_dcc_2020)



betweenness_central_80_minvar_dcc_2016 <- minimum_variance(return16, betweenness_central_80_dcc_2016)
betweenness_central_80_minvar_dcc_2017 <- minimum_variance(return17, betweenness_central_80_dcc_2017)
betweenness_central_80_minvar_dcc_2018 <- minimum_variance(return18, betweenness_central_80_dcc_2018)
betweenness_central_80_minvar_dcc_2019 <- minimum_variance(return19, betweenness_central_80_dcc_2019)
betweenness_central_80_minvar_dcc_2020 <- minimum_variance(return20, betweenness_central_80_dcc_2020)


betweenness_central_80_minvar_dcc =  data_table (betweenness_central_80_minvar_dcc_2016,
                                                 betweenness_central_80_minvar_dcc_2017,
                                                 betweenness_central_80_minvar_dcc_2018,
                                                 betweenness_central_80_minvar_dcc_2019,
                                                 betweenness_central_80_minvar_dcc_2020)

# n = 90 

betweenness_central_90_dcc_2016 <-  betweenness_central(calculate_dcc2016, 90)
betweenness_central_90_dcc_2017 <-  betweenness_central(calculate_dcc2017, 90)
betweenness_central_90_dcc_2018 <-  betweenness_central(calculate_dcc2018, 90)
betweenness_central_90_dcc_2019 <-  betweenness_central(calculate_dcc2019, 90)
betweenness_central_90_dcc_2020 <-  betweenness_central(calculate_dcc2020, 90)


betweenness_central_90_equalweight_dcc_2016 <- equal_weight(return16, betweenness_central_90_dcc_2016)
betweenness_central_90_equalweight_dcc_2017 <- equal_weight(return17, betweenness_central_90_dcc_2017)
betweenness_central_90_equalweight_dcc_2018 <- equal_weight(return18, betweenness_central_90_dcc_2018)
betweenness_central_90_equalweight_dcc_2019 <- equal_weight(return19, betweenness_central_90_dcc_2019)
betweenness_central_90_equalweight_dcc_2020 <- equal_weight(return20, betweenness_central_90_dcc_2020)


betweenness_central_90_equalweight_dcc = data_table (betweenness_central_90_equalweight_dcc_2016,
                                                     betweenness_central_90_equalweight_dcc_2017,
                                                     betweenness_central_90_equalweight_dcc_2018,
                                                     betweenness_central_90_equalweight_dcc_2019,
                                                     betweenness_central_90_equalweight_dcc_2020)



betweenness_central_90_minvar_dcc_2016 <- minimum_variance(return16, betweenness_central_90_dcc_2016)
betweenness_central_90_minvar_dcc_2017 <- minimum_variance(return17, betweenness_central_90_dcc_2017)
betweenness_central_90_minvar_dcc_2018 <- minimum_variance(return18, betweenness_central_90_dcc_2018)
betweenness_central_90_minvar_dcc_2019 <- minimum_variance(return19, betweenness_central_90_dcc_2019)
betweenness_central_90_minvar_dcc_2020 <- minimum_variance(return20, betweenness_central_90_dcc_2020)


betweenness_central_90_minvar_dcc =  data_table (betweenness_central_90_minvar_dcc_2016,
                                                 betweenness_central_90_minvar_dcc_2017,
                                                 betweenness_central_90_minvar_dcc_2018,
                                                 betweenness_central_90_minvar_dcc_2019,
                                                 betweenness_central_90_minvar_dcc_2020)
# n = 100 

betweenness_central_100_dcc_2016 <-  betweenness_central(calculate_dcc2016, 100)
betweenness_central_100_dcc_2017 <-  betweenness_central(calculate_dcc2017, 100)
betweenness_central_100_dcc_2018 <-  betweenness_central(calculate_dcc2018, 100)
betweenness_central_100_dcc_2019 <-  betweenness_central(calculate_dcc2019, 100)
betweenness_central_100_dcc_2020 <-  betweenness_central(calculate_dcc2020, 100)

betweenness_central_100_equalweight_dcc_2016 <- equal_weight(return16, betweenness_central_100_dcc_2016)
betweenness_central_100_equalweight_dcc_2017 <- equal_weight(return17, betweenness_central_100_dcc_2017)
betweenness_central_100_equalweight_dcc_2018 <- equal_weight(return18, betweenness_central_100_dcc_2018)
betweenness_central_100_equalweight_dcc_2019 <- equal_weight(return19, betweenness_central_100_dcc_2019)
betweenness_central_100_equalweight_dcc_2020 <- equal_weight(return20, betweenness_central_100_dcc_2020)


betweenness_central_100_equalweight_dcc = data_table (betweenness_central_100_equalweight_dcc_2016,
                                                     betweenness_central_100_equalweight_dcc_2017,
                                                     betweenness_central_100_equalweight_dcc_2018,
                                                     betweenness_central_100_equalweight_dcc_2019,
                                                     betweenness_central_100_equalweight_dcc_2020)



betweenness_central_100_minvar_dcc_2016 <- minimum_variance(return16, betweenness_central_100_dcc_2016)
betweenness_central_100_minvar_dcc_2017 <- minimum_variance(return17, betweenness_central_100_dcc_2017)
betweenness_central_100_minvar_dcc_2018 <- minimum_variance(return18, betweenness_central_100_dcc_2018)
betweenness_central_100_minvar_dcc_2019 <- minimum_variance(return19, betweenness_central_100_dcc_2019)
betweenness_central_100_minvar_dcc_2020 <- minimum_variance(return20, betweenness_central_100_dcc_2020)


betweenness_central_100_minvar_dcc =  data_table (betweenness_central_100_minvar_dcc_2016,
                                                 betweenness_central_100_minvar_dcc_2017,
                                                 betweenness_central_100_minvar_dcc_2018,
                                                 betweenness_central_100_minvar_dcc_2019,
                                                 betweenness_central_100_minvar_dcc_2020)




# pearson
# n = 70 
betweenness_central_70_pearson_2016 <-  betweenness_central(calculate_pearson2016, 70)
betweenness_central_70_pearson_2017 <-  betweenness_central(calculate_pearson2017, 70)
betweenness_central_70_pearson_2018 <-  betweenness_central(calculate_pearson2018, 70)
betweenness_central_70_pearson_2019 <-  betweenness_central(calculate_pearson2019, 70)
betweenness_central_70_pearson_2020 <-  betweenness_central(calculate_pearson2020, 70)


betweenness_central_70_equalweight_pearson_2016 <- equal_weight(return16, betweenness_central_70_pearson_2016)
betweenness_central_70_equalweight_pearson_2017 <- equal_weight(return17, betweenness_central_70_pearson_2017)
betweenness_central_70_equalweight_pearson_2018 <- equal_weight(return18, betweenness_central_70_pearson_2018)
betweenness_central_70_equalweight_pearson_2019 <- equal_weight(return19, betweenness_central_70_pearson_2019)
betweenness_central_70_equalweight_pearson_2020 <- equal_weight(return20, betweenness_central_70_pearson_2020)


betweenness_central_70_equalweight_pearson = data_table (betweenness_central_70_equalweight_pearson_2016,
                                                      betweenness_central_70_equalweight_pearson_2017,
                                                      betweenness_central_70_equalweight_pearson_2018,
                                                      betweenness_central_70_equalweight_pearson_2019,
                                                      betweenness_central_70_equalweight_pearson_2020)



betweenness_central_70_minvar_pearson_2016 <- minimum_variance(return16, betweenness_central_70_pearson_2016)
betweenness_central_70_minvar_pearson_2017 <- minimum_variance(return17, betweenness_central_70_pearson_2017)
betweenness_central_70_minvar_pearson_2018 <- minimum_variance(return18, betweenness_central_70_pearson_2018)
betweenness_central_70_minvar_pearson_2019 <- minimum_variance(return19, betweenness_central_70_pearson_2019)
betweenness_central_70_minvar_pearson_2020 <- minimum_variance(return20, betweenness_central_70_pearson_2020)


betweenness_central_70_minvar_pearson =  data_table (betweenness_central_70_minvar_pearson_2016,
                                                  betweenness_central_70_minvar_pearson_2017,
                                                  betweenness_central_70_minvar_pearson_2018,
                                                  betweenness_central_70_minvar_pearson_2019,
                                                  betweenness_central_70_minvar_pearson_2020)


# n = 80 

betweenness_central_80_pearson_2016 <-  betweenness_central(calculate_pearson2016, 80)
betweenness_central_80_pearson_2017 <-  betweenness_central(calculate_pearson2017, 80)
betweenness_central_80_pearson_2018 <-  betweenness_central(calculate_pearson2018, 80)
betweenness_central_80_pearson_2019 <-  betweenness_central(calculate_pearson2019, 80)
betweenness_central_80_pearson_2020 <-  betweenness_central(calculate_pearson2020, 80)


betweenness_central_80_equalweight_pearson_2016 <- equal_weight(return16, betweenness_central_80_pearson_2016)
betweenness_central_80_equalweight_pearson_2017 <- equal_weight(return17, betweenness_central_80_pearson_2017)
betweenness_central_80_equalweight_pearson_2018 <- equal_weight(return18, betweenness_central_80_pearson_2018)
betweenness_central_80_equalweight_pearson_2019 <- equal_weight(return19, betweenness_central_80_pearson_2019)
betweenness_central_80_equalweight_pearson_2020 <- equal_weight(return20, betweenness_central_80_pearson_2020)


betweenness_central_80_equalweight_pearson = data_table (betweenness_central_80_equalweight_pearson_2016,
                                                         betweenness_central_80_equalweight_pearson_2017,
                                                         betweenness_central_80_equalweight_pearson_2018,
                                                         betweenness_central_80_equalweight_pearson_2019,
                                                         betweenness_central_80_equalweight_pearson_2020)



betweenness_central_80_minvar_pearson_2016 <- minimum_variance(return16, betweenness_central_80_pearson_2016)
betweenness_central_80_minvar_pearson_2017 <- minimum_variance(return17, betweenness_central_80_pearson_2017)
betweenness_central_80_minvar_pearson_2018 <- minimum_variance(return18, betweenness_central_80_pearson_2018)
betweenness_central_80_minvar_pearson_2019 <- minimum_variance(return19, betweenness_central_80_pearson_2019)
betweenness_central_80_minvar_pearson_2020 <- minimum_variance(return20, betweenness_central_80_pearson_2020)


betweenness_central_80_minvar_pearson =  data_table (betweenness_central_80_minvar_pearson_2016,
                                                     betweenness_central_80_minvar_pearson_2017,
                                                     betweenness_central_80_minvar_pearson_2018,
                                                     betweenness_central_80_minvar_pearson_2019,
                                                     betweenness_central_80_minvar_pearson_2020)


# n = 90 

betweenness_central_90_pearson_2016 <-  betweenness_central(calculate_pearson2016, 90)
betweenness_central_90_pearson_2017 <-  betweenness_central(calculate_pearson2017, 90)
betweenness_central_90_pearson_2018 <-  betweenness_central(calculate_pearson2018, 90)
betweenness_central_90_pearson_2019 <-  betweenness_central(calculate_pearson2019, 90)
betweenness_central_90_pearson_2020 <-  betweenness_central(calculate_pearson2020, 90)


betweenness_central_90_equalweight_pearson_2016 <- equal_weight(return16, betweenness_central_90_pearson_2016)
betweenness_central_90_equalweight_pearson_2017 <- equal_weight(return17, betweenness_central_90_pearson_2017)
betweenness_central_90_equalweight_pearson_2018 <- equal_weight(return18, betweenness_central_90_pearson_2018)
betweenness_central_90_equalweight_pearson_2019 <- equal_weight(return19, betweenness_central_90_pearson_2019)
betweenness_central_90_equalweight_pearson_2020 <- equal_weight(return20, betweenness_central_90_pearson_2020)


betweenness_central_90_equalweight_pearson = data_table (betweenness_central_90_equalweight_pearson_2016,
                                                         betweenness_central_90_equalweight_pearson_2017,
                                                         betweenness_central_90_equalweight_pearson_2018,
                                                         betweenness_central_90_equalweight_pearson_2019,
                                                         betweenness_central_90_equalweight_pearson_2020)



betweenness_central_90_minvar_pearson_2016 <- minimum_variance(return16, betweenness_central_90_pearson_2016)
betweenness_central_90_minvar_pearson_2017 <- minimum_variance(return17, betweenness_central_90_pearson_2017)
betweenness_central_90_minvar_pearson_2018 <- minimum_variance(return18, betweenness_central_90_pearson_2018)
betweenness_central_90_minvar_pearson_2019 <- minimum_variance(return19, betweenness_central_90_pearson_2019)
betweenness_central_90_minvar_pearson_2020 <- minimum_variance(return20, betweenness_central_90_pearson_2020)


betweenness_central_90_minvar_pearson =  data_table (betweenness_central_90_minvar_pearson_2016,
                                                     betweenness_central_90_minvar_pearson_2017,
                                                     betweenness_central_90_minvar_pearson_2018,
                                                     betweenness_central_90_minvar_pearson_2019,
                                                     betweenness_central_90_minvar_pearson_2020)


# n = 100 

betweenness_central_100_pearson_2016 <-  betweenness_central(calculate_pearson2016, 100)
betweenness_central_100_pearson_2017 <-  betweenness_central(calculate_pearson2017, 100)
betweenness_central_100_pearson_2018 <-  betweenness_central(calculate_pearson2018, 100)
betweenness_central_100_pearson_2019 <-  betweenness_central(calculate_pearson2019, 100)
betweenness_central_100_pearson_2020 <-  betweenness_central(calculate_pearson2020, 100)


betweenness_central_100_equalweight_pearson_2016 <- equal_weight(return16, betweenness_central_100_pearson_2016)
betweenness_central_100_equalweight_pearson_2017 <- equal_weight(return17, betweenness_central_100_pearson_2017)
betweenness_central_100_equalweight_pearson_2018 <- equal_weight(return18, betweenness_central_100_pearson_2018)
betweenness_central_100_equalweight_pearson_2019 <- equal_weight(return19, betweenness_central_100_pearson_2019)
betweenness_central_100_equalweight_pearson_2020 <- equal_weight(return20, betweenness_central_100_pearson_2020)


betweenness_central_100_equalweight_pearson = data_table (betweenness_central_100_equalweight_pearson_2016,
                                                         betweenness_central_100_equalweight_pearson_2017,
                                                         betweenness_central_100_equalweight_pearson_2018,
                                                         betweenness_central_100_equalweight_pearson_2019,
                                                         betweenness_central_100_equalweight_pearson_2020)



betweenness_central_100_minvar_pearson_2016 <- minimum_variance(return16, betweenness_central_100_pearson_2016)
betweenness_central_100_minvar_pearson_2017 <- minimum_variance(return17, betweenness_central_100_pearson_2017)
betweenness_central_100_minvar_pearson_2018 <- minimum_variance(return18, betweenness_central_100_pearson_2018)
betweenness_central_100_minvar_pearson_2019 <- minimum_variance(return19, betweenness_central_100_pearson_2019)
betweenness_central_100_minvar_pearson_2020 <- minimum_variance(return20, betweenness_central_100_pearson_2020)


betweenness_central_100_minvar_pearson =  data_table (betweenness_central_100_minvar_pearson_2016,
                                                     betweenness_central_100_minvar_pearson_2017,
                                                     betweenness_central_100_minvar_pearson_2018,
                                                     betweenness_central_100_minvar_pearson_2019,
                                                     betweenness_central_100_minvar_pearson_2020)




