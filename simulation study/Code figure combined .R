
library(plotly)
library(ggplot2)
library(DataCombine)
library(darthpack)
library(dampack)
library(reshape2)
library(dampack)
# code for the DARTH colors for the figures
DARTHgreen      <- '#009999'  
DARTHyellow     <- '#FDAD1E'  
DARTHblue       <- '#006699' 
DARTHlightgreen <- '#00adad'
DARTHgray       <- '#666666'

install.packages("remotes")
remotes::install_github("feralaes/dampack")

data_combined <- read.csv("data_combined.csv")

data_combined$value[data_combined$variable == "Relative memory: dynamic / traditional"] <- 
  1 / data_combined$value[data_combined$variable == "Relative memory: dynamic / traditional"]
levels(data_combined$variable)[levels(data_combined$variable) == "Relative time: dynamic / traditional"] <- 
  "Relative speedup: dynamic / traditional"
levels(data_combined$Approach)

data_combined$yint <- 1
data_combined$yint[data_combined$variable == "Time (in seconds)" | data_combined$variable == "Memory (in MB)"] <- NA
data_combined$variable <- factor(data_combined$variable,levels(data_combined$variable)[c(4,3,1,2)])
data_combined$Approach <- factor(data_combined$Approach,levels(data_combined$Approach)[c(2,3,1)])

ggplot(data_combined, aes(x = n_states, value), group = variable) +
  geom_point(aes(color = Approach, shape = Approach)) +
  geom_line(aes(color = Approach)) + 
  scale_y_continuous(breaks = number_ticks(7)) +
  scale_x_continuous(breaks = number_ticks(6)) +
  xlab("Number of health states") + 
  ylab(element_blank()) +  
  facet_wrap(~variable, scales ="free_y", nrow = 2, 
             labeller = as_labeller(c)) + 
  geom_hline(aes(yintercept = yint), linetype = "dashed") + 
  theme_bw(base_size = 16) + 
  theme(legend.position = "bottom", legend.title = element_blank()) +
  theme(strip.background = element_rect(fill = DARTHgreen)) +
  theme(strip.text = element_text(color = 'white')) +
  scale_color_manual(values = c(DARTHyellow, DARTHblue,DARTHgreen))
