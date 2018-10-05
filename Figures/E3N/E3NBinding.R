# Get path of current file
curr.fpath = dirname(sys.frame(1)$ofile)
# Source plot formatting
setwd(curr.fpath)
source("../plotting_format.R")
library(BSgenome.Dmelanogaster.UCSC.dm3)

# Define genomic regions from Crocker, 2015 paper
genome = Dmelanogaster
E3N    = getSeq(genome, "chrX", 4915195,4915486)
H7     = getSeq(genome, "chrX", 4933945,4935001)

# Load single-mode nucleotide model for Exd-UbxIVa
ubx = file.parser("UbxIVa.csv")
ubx.full = ubx[[2]][[20]]$NB
ubx.12 = ubx.full[(3*4+1):(15*4)]
ubx.8  = ubx.12[(1*4+1):(9*4)]

# Plot logos
logo(betas=ubx.full, display = FALSE, save.path = "UbxFull.pdf", ylim=c(-4, 4))
logo(betas=ubx.12, display = FALSE, save.path = "Ubx12.pdf", ylim=c(-4, 4))
logo(betas=ubx.8, display = FALSE, save.path = "Ubx8.pdf", ylim=c(-4, 4))

# Create new 'fits'
ubx[[1]]$k[1] = 8
ubx[[1]]$Di[1]= FALSE
ubx[[2]][[1]]$NB = ubx.8
ubx[[2]][[1]]$DB = NULL
ubx[[1]]$k[2] = 12
ubx[[1]]$Di[2]= FALSE
ubx[[2]][[2]]$NB = ubx.12
ubx[[2]][[2]]$DB = NULL

# Plot tracks normalized to THEORETICAL MAX. Start with full
plot.E3N = function(idx) {
  E3N.score = score.genome(genomicSequence = E3N, fits = ubx, index = idx)
  windows = ncol(E3N.score)
  # Build a data frame to plot the sequences
  E3N.score = data.frame(Group = as.character(c(rep(1, each=2*windows), rep(3:2, each=windows))), 
                         Position = rep(1:windows, 4), 
                         Affinity=rep(c(E3N.score[1,], -E3N.score[2,]), 2))
  # Plot E3N genomic regions
  ggplot(E3N.score, aes(xmin=Position, xmax=Position+(as.numeric(ubx[[1]]$k[idx])-1), ymin=0, ymax=Affinity, fill=Group)) +
    coord_fixed(ylim=c(-max(E3N.score$Affinity),max(E3N.score$Affinity))) +
    expand_limits(x=300) +
    geom_rect(color=NA) +
    scale_fill_manual(values=c("#FFFFFF", alpha(c("#FF0000", "#000000"), .5)), 
                      labels=c("", "Reverse", "Forward")) +
    guides(fill=guide_legend(reverse=TRUE)) + 
#    scale_y_continuous(breaks=c(-scale:scale)/1000, labels=formatC(round(abs(-scale:scale)), digits=1, format="f")) + 
    xlab("Position") +
    ylab(expression(paste("Relative Affinity", sep=""))) +
    ggtheme +
    theme(legend.title=element_blank(), legend.position=c(0.01,0.01), legend.justification=c(0,0)) 
}

plot.E3N(20)
ggsave(filename = "E3N_Full.pdf", width=6, height=6)

plot.E3N(1)
ggsave(filename = "E3N_8.pdf", width=6, height=6)

plot.E3N(2)
ggsave(filename = "E3N_12.pdf", width=6, height=6)