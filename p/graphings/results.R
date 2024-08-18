
#install.packages("R.matlab")
library(R.matlab)
library(dplyr)

home_dir = "C:\\Users\\unnav\\Dropbox\\Education\\OSU\\Ongoing_Research\\"
project_dir = paste0(home_dir, "2nd Year Paper\\political-polarization")

data_dir = paste0(project_dir,"\\d")
graph_dir = paste0(project_dir,"\\v")

setwd(data_dir)

## want to pull in all the shares for all combos of taxes and utilities
shares_grid = matrix(1:5,ncol = 5)

numlist = c("0.25", "0.10", "0.40", "0.05", "0.09")

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

shares_grid = as.data.frame(shares_grid)
colnames(shares_grid) = c("share2A", "share2B", "ubonus", "t1", "t2")

