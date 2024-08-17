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


setwd(paste0(homedir, "d/results-20240602"))
agrid = read.csv("agrid.csv", header = F)
lgrid = read.csv("lgrid.csv", header = F)
GaA = read.csv("GaA.csv", row.names = NULL)
GaB = read.csv("GaB.csv", row.names = NULL)
GbA = read.csv("GbA.csv", row.names = NULL)
GbB = read.csv("GbB.csv", row.names = NULL)

VOTESa = read.csv("VOTESa.csv", row.names = NULL)
VOTESb = read.csv("VOTESb.csv", row.names = NULL)
acondA = read.csv("acondA.csv", row.names = NULL) 
acondB = read.csv("acondB.csv", row.names = NULL)

## Distributions ----

acondA = acondA[c(2,4,7),]
colnames(acondA) = c('Productivity', as.character(agrid))
ldat = gather(acondA, assets,probability, 2:251)
ldat = ldat %>% 
  group_by(Productivity) %>%
  mutate(Productivity = as.factor(Productivity), 
         assets = as.numeric(assets))

colnames(ldat) = c("Productivity", "Savings", "Probability")
# Create ridge plot
ggplot(ldat, aes(x = Savings, y = Probability, fill = Productivity)) + 
  geom_area(alpha = 0.6, colour = "black") + 
  # scale_y_continuous(trans="log10") 
  labs(x="Asset Grid", y="Probability",
       title="Household Asset Distribution by Productivity") + 
  theme_ipsum_rc() + 
  scale_fill_manual(values = rev(pinks))

ggsave(
  filename = paste0(homedir,"v/acondA.png"),
  plot = last_plot(),
  scale = 1,
  units = "in",
  height = 5,
  width = 8,
  dpi = 400,
  bg = NULL
)

acondB = acondB[c(2,4,7),]
colnames(acondB) = c('Productivity', as.character(agrid))
ldat = gather(acondB, assets,probability, 2:251)
ldat = ldat %>% 
  group_by(Productivity) %>%
  mutate(Productivity = as.factor(Productivity), 
         assets = as.numeric(assets))

colnames(ldat) = c("Productivity", "Savings", "Probability")
# Create ridge plot
ggplot(ldat, aes(x = Savings, y = Probability, fill = Productivity)) + 
  geom_area(alpha = 0.6, colour = "black") + 
  labs(x="Asset Grid", y="Probability",
       title="Household Asset Distribution by Productivity") + 
  theme_ipsum_rc() + 
  scale_color_viridis(discrete = TRUE)

ggsave(
  filename = paste0(homedir,"v/acondB.png"),
  plot = last_plot(),
  scale = 1,
  units = "in",
  height = 5,
  width = 8,
  dpi = 400,
  bg = NULL
)

## decision rules -----

GaAdat = GaA[c(1,4,7),]
colnames(GaAdat) = c('Productivity', as.character(agrid))
ldat = gather(GaAdat, assets,probability, 2:251)
ldat = ldat %>% 
  group_by(Productivity) %>%
  mutate(Productivity = as.factor(Productivity), 
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
  filename = paste0(homedir,"v/GaA.png"),
  plot = last_plot(),
  scale = 1,
  units = "in",
  height = 5,
  width = 5,
  dpi = 400,
  bg = NULL
)

GaBdat = GaB[c(1,4,7),]
colnames(GaBdat) = c('Productivity', as.character(agrid))
ldat = gather(GaBdat, assets,probability, 2:251)
ldat = ldat %>% 
  group_by(Productivity) %>%
  mutate(Productivity = as.factor(Productivity), 
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
  filename = paste0(homedir,"v/GaB.png"),
  plot = last_plot(),
  scale = 1,
  units = "in",
  height = 5,
  width = 5,
  dpi = 400,
  bg = NULL
)

## gb

GbAdat = GbA[c(1,4,7),]
colnames(GbAdat) = c('Productivity', as.character(agrid))
ldat = gather(GbAdat, assets,probability, 2:251)
ldat = ldat %>% 
  group_by(Productivity) %>%
  mutate(Productivity = as.factor(Productivity), 
         assets = as.numeric(assets))

colnames(ldat) = c("Productivity", "Savings", "Savings Choice")
# Create ridge plot
ggplot(ldat, aes(x = Savings, y = `Savings Choice`, group = Productivity,
                 color = Productivity)) + 
  geom_line(size = 1) + 
  scale_color_viridis(discrete = TRUE) +
  labs(x="Asset Grid", y="Savomgs Choice",
       title="Asset Choices by Productivity",
       subtitle = "Household B under Party A") + 
  theme_ipsum()+
  scale_x_discrete(breaks=c(seq(1,10,2), 2), 
                   labels= as.character(agrid[c(seq(1,250,50), 50)]))
ggsave(
  filename = paste0(homedir,"v/GbA.png"),
  plot = last_plot(),
  scale = 1,
  units = "in",
  height = 5,
  width = 5,
  dpi = 400,
  bg = NULL
)

GbBdat = GbB[c(1,4,7),]
colnames(GbBdat) = c('Productivity', as.character(agrid))
ldat = gather(GbBdat, assets,probability, 2:251)
ldat = ldat %>% 
  mutate(Productivity = as.factor(Productivity), 
         assets = as.numeric(assets))

colnames(ldat) = c("Productivity", "Savings", "Savings Choice")
# Create ridge plot
ggplot(ldat, aes(x = Savings, y = `Savings Choice`, group = Productivity,
                 color = Productivity)) + 
  geom_line(size = 1) + 
  scale_color_manual(values = rev(pinks)) +
  labs(x="Asset Grid", y="Savings Choice",
       title="Asset Choices by Productivity",
       subtitle = "Household B under Party B") + 
  theme_ipsum() +
  scale_x_discrete(breaks=c(seq(1,10,2), 2), 
                   labels= as.character(agrid[c(seq(1,250,50), 50)]))

ggsave(
  filename = paste0(homedir,"v/GbB.png"),
  plot = last_plot(),
  scale = 1,
  units = "in",
  height = 5,
  width = 5,
  dpi = 400,
  bg = NULL
)

dat = GbA[c(1,4,7),c(1,2:251)] - GbB[c(1,4,7),c(1,2:251)]
dat[,1] =  GbA[c(1,4,7),1]
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

dat = GaA[c(1,4,7),c(1,2:251)] - GaB[c(1,4,7),c(1,2:251)]
dat[,1] =  GbA[c(1,4,7),1]
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



