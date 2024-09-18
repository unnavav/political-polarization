# DATA VIZ with VAAS
#   June 2024
#   
#   I use this program to construct the visualizations for the second year
#   field paper. This program pulls data from what is output from the matlab
#   programs, as well as the data on presidential elections.

homedir = "C:/Users/unnav/Dropbox/Education/OSU/Ongoing_Research/2nd Year Paper/political-polarization/"
setwd(paste0(homedir, "d"))

library(tidyverse)
library(readxl)
library(data.table)
library(R.matlab)
library(DescTools)

library(ggplot2)
library(extrafont)
library(showtext)
library(RColorBrewer)
library(viridis)
library(hrbrthemes)

library(ggridges)

## Presidential Data ----

pres = read_xlsx("elections_results.xlsx")
pres = pres[order(pres$year),]

pres$diff = abs(pres$`candidate 1` - pres$`candidate 2`)*100

ggplot(pres, aes(x=year, y=diff)) + 
  geom_point(
    color="black",
    shape=21,
    alpha=0.5,
    size=2,
    stroke = 2
  ) + 
  theme_classic() + 
  ylim(-10, 30) + 
  geom_smooth(method = lm, color = "#ff6361", fill="darkgoldenrod1", se=TRUE, size=2)  + 
  geom_hline(yintercept = 5.3, linetype = "dashed", color = "azure4", size = 1) + 
  geom_hline(yintercept = 0, color = "black") + 
  ylab("Percent (%)") + 
  xlab("Year") + 
  theme(text = element_text(size = 18)) 


ggsave(
  filename = paste0(homedir,"v/presidents.png"),
  plot = last_plot(),
  scale = 1,
  units = "in",
  height = 5,
  width = 8,
  dpi = 400,
  bg = NULL
)

# decision rules ----

pinks = c("#06365c", "#7d6eae", "#ffa5ec")

ochres = c("#5c0000", "#bf6801", "#ffe005")


agrid = readMat("agrid.mat") %>% unlist() %>% as.vector()
lgrid = readMat("lgrid.mat") %>%  data.table::transpose() %>%
  unlist() %>% as.vector()

lgrid = lgrid %>% round(digits = 2)

setwd(paste0(homedir, "d\\personal_return0.05"))
results = readMat("results.mat")


GaA = results$Ga[,,1]
GaB = results$Ga[,,2]
GbA = results$Gb[,,1]
GbB = results$Gb[,,2]

VOTESa = results$VOTESa
VOTESb = results$VOTESb

acondA = readMat('acondA.mat') %>% as.data.frame()
acondB = readMat('acondB.mat') %>% as.data.frame()

setwd(paste0(homedir, "d\\personal_return0.25"))
results = readMat("results.mat")


## Distributions ----
colgrid = lgrid
acondA2 = cbind(colgrid, acondA)
acondA2 = acondA2[c(2,4,7),]
colnames(acondA2) = c('Productivity', as.character(round(agrid,digits=2)))
ldat = gather(acondA2, assets,probability, 2:251)
ldat = ldat %>% 
  group_by(Productivity) %>%
  mutate(Productivity = round(Productivity, digits = 2),
         Productivity = as.factor(Productivity), 
         assets = as.numeric(assets))

colnames(ldat) = c("Productivity", "Savings", "Probability")

ldat_filt = ldat %>% filter(Savings < 50)

ggplot(ldat_filt, aes(x = Savings, y = Probability, fill = Productivity)) + 
  geom_area(alpha = 0.6, colour = "black") + 
  # scale_y_continuous(trans="log10") 
  labs(x="Asset Grid", y="Probability",
       title="Household Asset Distribution by Productivity") + 
  theme_ipsum_rc() + 
  scale_fill_manual(values = rev(pinks))  
  

ggsave(
  filename = paste0(homedir,"v/acondA_u05.png"),
  plot = last_plot(),
  scale = 1,
  units = "in",
  height = 5,
  width = 8,
  dpi = 400,
  bg = NULL
)

acondB = cbind(colgrid, acondB)
acondB = acondB[c(2,4,7),]
colnames(acondB) = c('Productivity', as.character(round(agrid,digits=2)))
ldatB = gather(acondB, assets,probability, 2:251)
ldatB = ldatB %>% 
  group_by(Productivity) %>%
  mutate(Productivity = round(Productivity, digits = 2),
         Productivity = as.factor(Productivity), 
         assets = as.numeric(assets))

colnames(ldatB) = c("Productivity", "Savings", "Probability")
# Create ridge plot
ggplot(ldatB, aes(x = Savings, y = Probability, fill = Productivity)) + 
  geom_area(alpha = 0.6, colour = "black") + 
  labs(x="Asset Grid", y="Probability",
       title="Household Asset Distribution by Productivity") + 
  theme_ipsum_rc() + 
  scale_color_viridis(discrete = TRUE)

ggsave(
  filename = paste0(homedir,"v/acondB_u05.png"),
  plot = last_plot(),
  scale = 1,
  units = "in",
  height = 5,
  width = 8,
  dpi = 400,
  bg = NULL
)

dist_diff = ldatB$probability-ldat$Probability

ldat$dist_diff = dist_diff

ggplot(ldat, aes(x = Savings, y = dist_diff, fill = Productivity)) + 
  geom_area(alpha = 0.6, colour = "black") + 
  labs(x="Asset Grid", y="Probability",
       title="Household Asset Distribution by Productivity") + 
  theme_ipsum_rc() + 
  scale_color_viridis(discrete = TRUE)

## decision rules -----

getDecRuleA(GaA, GaB, "u5")
getDecRuleB(GbA, GbB, "u5")

setwd(paste0(homedir, "d\\personal_return0.25"))
results = readMat("results.mat")

GaA = results$Ga[,,1]
GaB = results$Ga[,,2]
GbA = results$Gb[,,1]
GbB = results$Gb[,,2]

getDecRuleA(GaA, GaB, "u25")
getDecRuleB(GbA, GbB, "u25")

setwd(paste0(homedir, "d\\tax_scheme_0.00000.9500"))
results = readMat("results_pval10.mat")

GaA = results$Ga[,,1]
GaB = results$Ga[,,2]
GbA = results$Gb[,,1]
GbB = results$Gb[,,2]

getDecRuleA(GaA, GaB, lgrid, "t95")
getDecRuleB(GbA, GbB, lgrid, "t95")

dat = GbB - GbA
dat = cbind(lgrid, dat %>% as.data.frame())
dat = dat[c(1,4,7),]
colnames(dat) = c('Productivity', as.character(agrid))
ldat = gather(dat, assets,savings, 2:251)
colnames(ldat) = c("Productivity", "Savings", "Savings Choice")
ldat = ldat %>% 
  group_by(Productivity) %>%
  mutate(Productivity = as.factor(Productivity), 
         `Savings Choice` = as.numeric(`Savings Choice`))

ldat2 = lapply(ldat, DescTools::Winsorize)
# Create ridge plot
ggplot(ldat, aes(x = Savings, y = `Savings Choice`, 
                 group = Productivity, color = Productivity)) + 
  geom_line(size=1) + 
  scale_color_manual(values = rev(ochres)) + 
  labs(x="Asset Grid", y="Difference in Savings Choice",
       title="Differences in Asset Choices by Productivity",
       subtitle = "Household B") + 
  theme_ipsum() + 
  # scale_x_discrete(n.breaks=10)
  scale_x_discrete(breaks=c(seq(1,250,50), 350),
                     labels= as.character(agrid[c(seq(1,250,50), 50)]))

ggsave(
  filename = paste0(homedir,"v/Gbdiff.png"),
  plot = last_plot(),
  scale = 1,
  units = "in",
  height = 5,
  width = 8,
  dpi = 400,
  bg = NULL
)

dat = GaAdat[,c(1,2:251)] - GaBdat[,c(1,2:251)]
colnames(dat) = c('Productivity', as.character(agrid))
ldat = gather(dat, assets,savings, 2:251)
colnames(ldat) = c("Productivity", "Savings", "Savings Choice")
ldat = ldat %>% 
  group_by(Productivity) %>%
  mutate(Productivity = as.factor(Productivity), 
         `Savings Choice` = as.numeric(`Savings Choice`))

# Create ridge plot
ggplot(ldat, aes(x = Savings, y = `Savings Choice`, 
                 group = Productivity, color = Productivity)) + 
  geom_line(size=1) + 
  scale_color_manual(values = rev(ochres)) + 
  labs(x="Asset Grid", y="Difference in Savings Choice",
       title="Differences in Asset Choices by Productivity",
       subtitle = "Household A") + 
  theme_ipsum() + 
  # scale_x_discrete(n.breaks=10)
  scale_x_discrete(breaks=c(seq(1,250,50), 350),
                   labels= as.character(agrid[c(seq(1,250,50), 50)]))

ggsave(
  filename = paste0(homedir,"v/Gadiff.png"),
  plot = last_plot(),
  scale = 1,
  units = "in",
  height = 5,
  width = 8,
  dpi = 400,
  bg = NULL
)

### functions ----

getDecRuleA = function(GaA, GaB, colgrid, fileext) {
  
  GaAdat = cbind(colgrid, GaA)
  GaAdat = GaAdat[c(1,4,7),] %>% as.data.frame()
  colnames(GaAdat) = c('Productivity', as.character(agrid))
  ldat = gather(GaAdat, assets,probability, 2:251)
  ldat = ldat %>% 
    group_by(Productivity) %>%
    mutate(Productivity = round(Productivity, digits = 2),
           Productivity = as.factor(Productivity), 
           assets = as.numeric(assets))
  
  colnames(ldat) = c("Productivity", "Savings", "Savings Choice")
  # Create ridge plot
  ggplot(ldat, aes(x = Savings, y = `Savings Choice`, group = Productivity,
                   color = Productivity)) + 
    geom_line(size = 1) + 
    scale_color_viridis(discrete = TRUE) +
    labs(x="Asset Grid", y="Savomgs Choice",
         title="Asset Choices by Productivity",
         subtitle = "Household A under Party A") + 
    theme_ipsum()
  
  ggsave(
    filename = paste0(homedir,"v/GaA_",fileext,".png"),
    plot = last_plot(),
    scale = 1,
    units = "in",
    height = 5,
    width = 5,
    dpi = 400,
    bg = NULL
  )
  
  GaBdat = cbind(colgrid, GaB)
  GaBdat = GaBdat[c(1,4,7),] %>% as.data.frame()
  colnames(GaBdat) = c('Productivity', as.character(agrid))
  ldat = gather(GaBdat, assets,probability, 2:251)
  ldat = ldat %>% 
    group_by(Productivity) %>%
    mutate(Productivity = round(Productivity, digits = 2),
           Productivity = as.factor(Productivity), 
           assets = as.numeric(assets))
  
  colnames(ldat) = c("Productivity", "Savings", "Savings Choice")
  # Create ridge plot
  ggplot(ldat, aes(x = Savings, y = `Savings Choice`, group = Productivity,
                   color = Productivity)) + 
    geom_line(size = 1) + 
    scale_color_manual(values = rev(pinks)) +
    labs(x="Asset Grid", y="Savomgs Choice",
         title="Asset Choices by Productivity",
         subtitle = "Household A under Party B") + 
    theme_ipsum()
  
  ggsave(
    filename = paste0(homedir,"v/GaB_",fileext,".png"),
    plot = last_plot(),
    scale = 1,
    units = "in",
    height = 5,
    width = 5,
    dpi = 400,
    bg = NULL
  )
  
}

getDecRuleB = function(GbA, GbB, colgrid, fileext) {
  
  GaAdat = cbind(colgrid, GbA)
  GaAdat = GaAdat[c(1,4,7),] %>% as.data.frame()
  colnames(GaAdat) = c('Productivity', as.character(agrid))
  ldat = gather(GaAdat, assets,probability, 2:251)
  ldat = ldat %>% 
    group_by(Productivity) %>%
    mutate(Productivity = round(Productivity, digits = 2),
           Productivity = as.factor(Productivity), 
           assets = as.numeric(assets))
  
  colnames(ldat) = c("Productivity", "Savings", "Savings Choice")
  # Create ridge plot
  ggplot(ldat, aes(x = Savings, y = `Savings Choice`, group = Productivity,
                   color = Productivity)) + 
    geom_line(size = 1) + 
    scale_color_viridis(discrete = TRUE) +
    labs(x="Asset Grid", y="Savomgs Choice",
         title="Asset Choices by Productivity",
         subtitle = "Household A under Party A") + 
    theme_ipsum()
  
  ggsave(
    filename = paste0(homedir,"v/GbA_",fileext,".png"),
    plot = last_plot(),
    scale = 1,
    units = "in",
    height = 5,
    width = 5,
    dpi = 400,
    bg = NULL
  )
  
  GaBdat = cbind(colgrid, GbB)
  GaBdat = GaBdat[c(1,4,7),] %>% as.data.frame()
  colnames(GaBdat) = c('Productivity', as.character(agrid))
  ldat = gather(GaBdat, assets,probability, 2:251)
  ldat = ldat %>% 
    group_by(Productivity) %>%
    mutate(Productivity = round(Productivity, digits = 2),
           Productivity = as.factor(Productivity), 
           assets = as.numeric(assets))
  
  colnames(ldat) = c("Productivity", "Savings", "Savings Choice")
  # Create ridge plot
  ggplot(ldat, aes(x = Savings, y = `Savings Choice`, group = Productivity,
                   color = Productivity)) + 
    geom_line(size = 1) + 
    scale_color_manual(values = rev(pinks)) +
    labs(x="Asset Grid", y="Savomgs Choice",
         title="Asset Choices by Productivity",
         subtitle = "Household A under Party B") + 
    theme_ipsum()
  
  ggsave(
    filename = paste0(homedir,"v/GbB_",fileext,".png"),
    plot = last_plot(),
    scale = 1,
    units = "in",
    height = 5,
    width = 5,
    dpi = 400,
    bg = NULL
  )
  
}


