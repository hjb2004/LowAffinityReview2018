# Get path of current file
curr.fpath = dirname(sys.frame(1)$ofile)
# Source plot formatting
setwd(curr.fpath)
source("../plotting_format.R")

# Load Labial and UbxIVa first
lab = file.parser("Lab.csv")
ubx = file.parser("UbxIVa.csv")

# Plot Ubx/Lab monomer NUCLEOTIDE fits in the right orientation and remove extra bases (align)
logo(fits = lab, index = 9, rc = T, ylim = c(-5,5), save.path = "LabMonomer.pdf", display = F, r.del=1)
logo(fits = ubx, index = 8, rc = F, ylim = c(-5,5), save.path = "UbxMonomer.pdf", display = F, l.del=1)

# Plot differential logo. Begin by computing the differential betas
lab.betas = lab[[2]][[9]]$NB
lab.betas = rev(lab.betas)
ubx.betas = ubx[[2]][[8]]$NB

lab.betas = lab.betas[1:(12*4)]
ubx.betas = ubx.betas[5:(13*4)]

# 'Fix' previous fit outputs to enable scoring later on
lab[[2]][[9]]$NB = lab.betas
lab[[1]]$k[9] = 12
ubx[[2]][[8]]$NB = ubx.betas
ubx[[1]]$k[8] = 12

# Create an 'artificial' fit to override MultinomialMethods restrictions
lab[[1]]$k[1] = 12
lab[[1]]$Di[1]= FALSE
lab[[2]][[1]]$NB = lab.betas-ubx.betas
lab[[2]][[1]]$DB = NULL

logo(fits = lab, index = 1, ylim=c(-5,5), save.path = "LabUbxDiff.pdf", display = F)
# Find top and bottom sequences
lab.top.seq = max.seq(lab, 1)$BestSeq
lab[[2]][[1]]$NB = ubx.betas-lab.betas 
ubx.top.seq = max.seq(lab, 1)$BestSeq
logo(fits = lab, index = 1, ylim=c(-5,5), save.path = "LabUbxDiffFlipped.pdf", display = F)

# record affinities
df = data.frame(Protein="Lab", SeqType="Labial Preference for Lab/Ubx", Affinity = 1, stringsAsFactors = FALSE)
df$Affinity[1] = max(score.genome(DNAString(lab.top.seq), lab, 9))/max.seq(lab, 9)$MaxAffinity
df = rbind(df, c("Ubx", "Labial Preference for Lab/Ubx", max(score.genome(DNAString(lab.top.seq), ubx, 8))/max.seq(ubx, 8)$MaxAffinity),
           c("Lab", "UbxIVa Preference for Lab/Ubx", max(score.genome(DNAString(ubx.top.seq), lab, 9))/max.seq(lab, 9)$MaxAffinity),
           c("Ubx", "UbxIVa Preference for Lab/Ubx", max(score.genome(DNAString(ubx.top.seq), ubx, 8))/max.seq(ubx, 8)$MaxAffinity))
df$Affinity = as.numeric(df$Affinity)

# Bring Max into the picture. Need a 10bp logo, so fix Lab again
lab.betas = lab[[2]][[9]]$NB
lab.betas = lab.betas[1:(10*4)]
lab[[2]][[9]]$NB = lab.betas
lab[[1]]$k[9] = 10

# Load Max
max = file.parser("Max_NRLB_SELEX.csv")
max.betas = max[[2]][[27]]$NB
max.betas = max.betas[(2*4+1):(12*4)]
max[[2]][[27]]$NB = max.betas
max[[1]]$k[27] = 10

# Create an artificial fit to override MultinomialMethods
lab[[1]]$k[2] = 10
lab[[1]]$Di[2]= FALSE
lab[[2]][[2]]$NB = lab.betas-max.betas
lab[[2]][[1]]$DB = NULL

logo(fits = lab, index = 2, ylim = c(-7,7), save.path = "LabMaxDiff.pdf", display = F)
lab.max.top.seq = max.seq(lab, 2)$BestSeq
lab[[2]][[2]]$NB = max.betas-lab.betas
max.top.seq = max.seq(lab, 2)$BestSeq

# record affinities (normalized to the best sequence)
df = rbind(df, c("Lab", "Labial Preference for Lab/Max", max(score.genome(DNAString(lab.max.top.seq), lab, 9))/max.seq(lab, 9)$MaxAffinity),
           c("Max", "Labial Preference for Lab/Max", max(score.genome(DNAString(lab.max.top.seq), max, 27))/max.seq(max, 27)$MaxAffinity),
           c("Lab", "Max Preference for Lab/Max", max(score.genome(DNAString(max.top.seq), lab, 9))/max.seq(lab, 9)$MaxAffinity),
           c("Max", "Max Preference for Lab/Max", max(score.genome(DNAString(max.top.seq), max, 27))/max.seq(max, 27)$MaxAffinity))
df$Affinity = as.numeric(df$Affinity)

# Write output
write.table(x = df, file = "DifferentialAffinities.csv", append = FALSE, quote = FALSE, sep = ",", row.names = FALSE)