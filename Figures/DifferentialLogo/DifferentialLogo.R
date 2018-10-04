# Get path of current file
curr.fpath = dirname(sys.frame(1)$ofile)
# Source plot formatting
setwd(curr.fpath)
source("../plotting_format.R")

# Load Labial and UbxIVa first
lab = file.parser("Lab.csv")
ubx = file.parser("UbxIVa.csv")

# Plot Ubx/Lab monomer NUCLEOTIDE fits in the right orientation and remove extra bases (align)
logo(fits = lab, index = 9, rc = T, ylim = c(-4,4), save.path = "LabMonomer.pdf", display = F, r.del=1)
logo(fits = ubx, index = 8, rc = F, ylim = c(-4,4), save.path = "UbxMonomer.pdf", display = F, l.del=1)

# Plot differential logo. Begin by computing the differential betas
lab.betas = lab[[2]][[9]]$NB
lab.betas = rev(lab.betas)
ubx.betas = ubx[[2]][[8]]$NB

lab.betas = lab.betas[1:(12*4)]
ubx.betas = ubx.betas[5:(13*4)]
logo(betas=(lab.betas-ubx.betas), ylim=c(-4,4), save.path = "LabUbxDiff.pdf", display = F)