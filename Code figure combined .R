
library(plotly)
library(ggplot2)
library(DataCombine)
library(darthpack)
library(dampack)
library(reshape2)

# code for the DARTH colors for the figures
DARTHgreen      <- '#009999'  
DARTHyellow     <- '#FDAD1E'  
DARTHblue       <- '#006699' 
DARTHlightgreen <- '#00adad'
DARTHgray       <- '#666666'

data_combined <- read.csv("data_combined.csv")

ggplot(data_combined, aes(x = n_states, value), group = variable) +
  geom_point(aes(color = Approach, shape = Approach)) +
  geom_line(aes(color = Approach)) + 
  scale_y_continuous(breaks = number_ticks(7)) +
  scale_x_continuous(breaks = number_ticks(6)) +
  xlab("Number of health states") + 
  ylab(element_blank()) +  
  facet_wrap(~ variable, scales ="free_y", nrow = 2, 
             labeller = as_labeller(c)) + 
  theme_bw(base_size = 16) + 
  theme(legend.position = "bottom", legend.title = element_blank()) +
  theme(strip.background = element_rect(fill = DARTHgreen)) +
  theme(strip.text = element_text(color = 'white')) +
  scale_color_manual(values = c(DARTHgreen, DARTHyellow, DARTHblue))

