
#install.packages("R.matlab")
library(R.matlab)
library(dplyr)

home_dir = "C:\\Users\\unnav\\Dropbox\\Education\\OSU\\Ongoing_Research\\"
project_dir = paste0(home_dir, "2nd Year Paper\\political-polarization")

data_dir = paste0(project_dir,"\\d")
graph_dir = paste0(project_dir,"\\v")

setwd(data_dir)

agrid = readMat("agrid.mat")
amu = readMat("amu.mat") %>% unlist()

## want to pull in all the shares for all combos of taxes and utilities
shares_grid = matrix(1:5,ncol = 5)

numlist = c("0.25", "0.10", "0.40", "0.05", "0.09", "0.12", "0.14", "0.15")

for (num in numlist) {
  u_results = readMat(paste0(data_dir, "\\personal_return",num,"\\results.mat"))
  
  shares_grid = rbind(shares_grid, c(u_results$share2A,
                                     u_results$share2B,
                                     as.numeric(num),
                                     0.086, 0.181))
}

shares_grid = shares_grid[2:nrow(shares_grid),] #getting rid of initialization row

# now tax levels

taxlist = c("0.00000.9750", "0.00000.9500", "0.00000.9000","0.00001.0000")

for (num in taxlist) {
  u_results = readMat(paste0(data_dir, "\\tax_scheme_",num,"\\results_pval10.mat"))
  
  t1 = substr(num, start = 1, stop = 6) %>% as.numeric()
  t2 = substr(num, start = 7, stop = 12) %>% as.numeric()
  
  shares_grid = rbind(shares_grid, c(u_results$share2A,
                                     u_results$share2B,
                                     0.1,
                                     t1, t2))
}


taxlist = c("0.00000.9500", "0.00000.9562","0.00000.9625")

for (num in taxlist) {
  u_results = readMat(paste0(data_dir, "\\tax_scheme_",num,"\\results_pval05.mat"))
  
  t1 = substr(num, start = 1, stop = 6) %>% as.numeric()
  t2 = substr(num, start = 7, stop = 12) %>% as.numeric()
  
  shares_grid = rbind(shares_grid, c(u_results$share2A,
                                     u_results$share2B,
                                     0.05,
                                     t1, t2))
}

shares_grid = as.data.frame(shares_grid)
colnames(shares_grid) = c("share2A", "share2B", "ubonus", "t1", "t2")

## getting laffer curves

results_u1 =  readMat(paste0(data_dir, "\\personal_return0.25\\results.mat"))
results_u2 =  readMat(paste0(data_dir, "\\personal_return0.05\\results.mat"))

results_t1 = readMat(paste0(data_dir, "\\tax_scheme_0.00000.9750\\results_pval10.mat"))
results_t2 = readMat(paste0(data_dir, "\\tax_scheme_0.00000.9000\\results_pval10.mat"))
results_t3 = readMat(paste0(data_dir, "\\tax_scheme_0.00000.9500\\results_pval05.mat"))
results_t4 = readMat(paste0(data_dir, "\\tax_scheme_0.00000.9500\\results_pval10.mat"))
results_t5 = readMat(paste0(data_dir, "\\tax_scheme_0.00001.0000\\results_pval10.mat"))
results_t6 = readMat(paste0(data_dir, "\\tax_scheme0.9250\\resultspval05.mat"))
results_t7 = readMat(paste0(data_dir, "\\tax_scheme0.9125\\resultspval05.mat"))

results_t5B = readMat(paste0(data_dir, "\\tax_scheme_Beqm_0.00000.9900\\results_pval10.mat"))

lorenz_A10 = getLorenz(results_t5)
gini_A10 = getGini(lorenz_A10)
quartiles_A10 = getQuartiles(lorenz_A10)

lorenz_B10 = getLorenz(results_t5B)
gini_B10 = getGini(lorenz_B10)
quartiles_B10 = getQuartiles(lorenz_B10)

lorenz_u05 = getLorenz(results_u2)
giniu05 = getGini(lorenz_u05)
quartiles = getQuartiles(lorenz_u05)

lorenz_u25 = getLorenz(results_u1)
giniu25 = getGini(lorenz_u25)
quartiles = getQuartiles(lorenz_u25)

lorenz_tB = getLorenz(results_t1)

#### functions --------

getQuartiles = function(lorenz) {
  quantile(lorenz$cdf, c(.10, 0.25, .5, .75, .9))
}


getGini = function(lorenz) {
  dc=lorenz$cdf[2:length(amu)] - lorenz$cdf[1:length(amu)-1]
  db = lorenz$wealth[2:length(amu)] - lorenz$wealth[1:length(amu)-1]
  
  B = sum(dc*db)
  A = .5-B
  GiniC1 = 2.0*A
}

getLorenz = function(results) {
  lorenz = t(results$adistrA) %>%
    rowSums()
  
  lorenz = lorenz %>% as.data.frame()
  colnames(lorenz) = c("pdf")
  
  lorenz = lorenz %>%
    mutate(cdf = cumsum(pdf))
  
  lorenz$amu = amu
  
  lorenz = lorenz %>%
    rowwise() %>%
    mutate(here = pdf*amu) %>%
    ungroup() %>%
    mutate(wealth = cumsum(here))
}

graphVotes = function(results, colfunc) {
  VOTESa = results$VOTESa
  VOTESb = results$VOTESb
  
  filled.contour(VOTESa, col = colfunc(2))
}

colfunc <- colorRampPalette(c("#C75DAB", "#009B9E"))