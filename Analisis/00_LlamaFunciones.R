code_dir <- "./Rfun/"
file_path <- list.files(code_dir,full.names = T)
for(f_path in file_path){source(f_path,encoding = "UTF-8")}

library(ggplot2)
library(ggpubr)
library(reshape2)
source(paste0('Rfun/','mi.tema.R')) #personalizacion ggplot2
source(paste0('Rfun/','read.admb.R'))
