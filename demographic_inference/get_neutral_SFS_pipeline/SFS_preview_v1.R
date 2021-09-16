## Author: Paulina G. Nuñez Valencia, LCG-UNAM, Mexico: pnunez@lcg.unam.mx, modified by Sergio Nigenda for Vaquita data
## NOTES: R 4.0,v1.0, July 23 2020
#--------------------------------------------------------------------------------------------------
#                                            GOAL
# Join previews SFS projections for chromosomes
#--------------------------------------------------------------------------------------------------

#################################################################################################
#                                        MAIN PROGRAM                                           #
#################################################################################################


# Load the R packages
library(plyr)
library(readr)
library(ggplot2)
library(RColorBrewer)

#Define functions ------------------------------ 
readtables <- function(filenames, dataframe){
  c <- 1
  for(file in filenames){
    dataframe <- cbind(dataframe,read_csv(file,skip = 1,col_names = c("prj", paste0("SNPs_", c)) )[2]) 
    c <- c+1
  }
  return(dataframe)
}


#Define variables --------------------
dir <- "/u/project/rwayne/jarobins/vaquita/Phsi_neutral_sfs/SFS_neutral/Neutral/SNPs"
dir.create("/u/project/rwayne/jarobins/vaquita/Phsi_neutral_sfs/SFS_neutral/Neutral/Preview")
outputdir <- "/u/project/rwayne/jarobins/vaquita/Phsi_neutral_sfs/SFS_neutral/Neutral/Preview"
pop <- "Vaquita"

#Main---------------------

#Obtain tables for each pop that contain projection values of each chromosome file

preview.files <- list.files(dir,pattern= pop,full.names = TRUE)
preview.values <- read_csv(preview.files[1], col_names = c("Projection","snp"), skip = 1)[1]
preview <- readtables(preview.files,preview.values)
write.table(preview, paste0(outputdir,"/SFS_PreviewProjection_",pop,".txt"),quote = F, row.names = F)
preview$SNPs <- rowSums(preview[c(2:22)])
preview$pop <- pop
assign(paste0("preview_",pop), preview)

#Join populations
preview_allpops <- rbind(preview_Vaquita)
colors <- brewer.pal(3,"Dark2")

#Plot preview
p1 <- ggplot(preview_allpops,aes(x = Projection, y=SNPs), group = 1, color=pop) + geom_line() + geom_point() + 
  facet_wrap(~pop,scales="free_x") + theme_light() +  scale_x_continuous(breaks=seq(2,80,2)) +
  labs(title = "Site Frequency Spectrum Projection Preview", y="Number of segregating sites") + 
  guides(color = F) + scale_color_manual(values = colors)

ggsave(paste0(outputdir,"/easySFS_PreviewProjections_Vaquita.pdf"),p1,height=5,width=10)
