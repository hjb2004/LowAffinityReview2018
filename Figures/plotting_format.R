# Load packages
library(grid)
library(hexbin)
library(scales)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(extrafont)
library(gridExtra)
library(Biostrings)
library(MultinomialMethods)

# ggplot theme
ggtheme = theme_bw() + theme(text=element_text(size=23, family="Helvetica"), axis.title.y=element_text(margin=margin(0,20,0,0)),
                axis.title.x=element_text(margin=margin(20,0,0,0)), aspect.ratio=1, axis.line=element_line(color="black", size=1), 
                axis.ticks=element_line(color="black", size=1), panel.border=element_blank(), legend.justification=c(1,0),
                legend.title.align=0.5, legend.position=c(.99,0.01))