
##R version 3.4.1  was used to run the following scripts

#####Libraries used

#flowCore_1.42.3
library("flowCore")
##flowViz_1.40.0
library("flowViz")
##ggplot2_2.2.1
library(ggplot2)
##lattice_0.20-35
library(lattice)
##gridExtra_2.3
library(gridExtra)
##dplyr_0.7.4
library(dplyr)
##plyr_1.8.4
library(plyr)
##reshape2_1.4.3
library(reshape2)
##lemon_0.4.1
library(lemon)


###########################################
#Flow cytometry data analysis             #
###########################################
##########################################
# to transform a .fcs file to a .csv file#
##########################################

#source("https://bioconductor.org/biocLite.R")
biocLite("flowCore")
biocLite("flowViz")
library("flowCore")
library("flowViz")
library(ggplot2)
library(lattice)
library(gridExtra)
library(dplyr)

#function to rename strains with their real number in the plate

num2coords = function (file_number) {
  ligne = (file_number-1) %/% 12 + 1
  colonne= (file_number-1) %% 12 +1 
  return(list (x= colonne, y = ligne))
}

coords2strain = function(list) {
  strain = ((list$x - 1)*8) + list$y
  return (strain)
}

####Transform FCS files in to CSV files

setwd("C:/Users/marsit/Dropbox/Ploidy-MA-CRL/")

##### Control 1 diploid
d = read.FCS("controls/2017-07-14_ploidy controls-v2-1.fcs")

df = as.data.frame(exprs(d)[,1:ncol(d)])
ech = rep("Diplo", nrow(df))
ctrl1 = cbind(ech, df)
write.csv(ctrl1,"controls/Control1.csv", row.names = FALSE, quote = F)
csv_ctrl1 = read.table("controls/Control1.csv", header = T, sep = ",")


##### Control 2 haploid
d = read.FCS("controls/2017-07-14_ploidy controls-v2-2.fcs")

df = as.data.frame(exprs(d)[,1:ncol(d)])
ech = rep("Haplo", nrow(df))
ctrl2 = cbind(ech, df)
write.csv(ctrl2,"controls/Control2.csv", row.names = FALSE, quote = F)
csv_ctrl2 = read.table("controls/Control2.csv", header = T, sep = ",")

# Tranform samples files from FCS to CSV files with the good strain number
########################################################################

setwd("C:/Users/marsit/Dropbox/Ploidy-MA-CRL/2018-06-22_MA-Ploidy-ABDE-tf")

file = dir(pattern="fcs")
file

for(i in 1:length(file)){

  d = read.FCS(file[i])
  df = as.data.frame(exprs(d)[,1:ncol(d)])
  
  file_number= strsplit (strsplit(file[i], ".", fixed = T)[[1]][1],"-", fixed = T) [[1]] %>% tail(n=1) %>% as.numeric()
  strain_number= coords2strain(num2coords(file_number))
  
  ech1 = cbind(strain_number, df)
  csv_name = paste0(strsplit(file[i], ".", fixed = T)[[1]][1], "_strain_",strain_number, ".csv")
  write.csv(ech1, csv_name, row.names = FALSE, quote = F)
}  

################################################################################################################################### 
##################################################################################
# Generate a table with only needed informations for fluorescence- cytometry data#
##################################################################################

library(ggplot2)
library(cowplot)
library(lattice)
library(gridExtra)
library(dplyr)
library(plyr)
library(reshape2)

setwd("C:/Users/marsit/Dropbox/Ploidy-MA-CRL")

#####################
# control 1 haploid #
#####################

csv_ctrl2 = read.table('C:/Users/marsit/Desktop/Dropbox/Ploidy-MA-CRL/2018-06-07_controls/Control3.csv', header = T, sep = ',')
head(csv_ctrl2)
csv_ctrl2$ech <- as.character(csv_ctrl2$ech)

#calculate fluorescence for a bin of 400 window

ctrl2_bin <- csv_ctrl2 %>% mutate(bin = GRN.B.HLog %/% 400)
sample = ctrl2_bin[1,1]

all_bins <- c(1:max(ctrl2_bin$bin)) %>% as.data.frame
colnames(all_bins) <- 'bin'

ctrl2_allbins <- full_join(ctrl2_bin, all_bins, by = 'bin') %>% mutate(ech = ifelse(is.na(ech), sample, ech)) %>% mutate(GRN.B.HLog = ifelse(is.na(ech), 0, GRN.B.HLog))

test0 <- ctrl2_allbins %>% group_by(ech, bin, GRN.B.HLog) %>% dplyr::summarise(count = length(GRN.B.HLog)) %>% as.data.frame

#head(test0)

#####################
# control 2 diploid #
#####################
csv_ctrl1 = read.table('C:/Users/marsit/Dropbox/Ploidy-MA-CRL/controls/Control8.csv', header = T, sep = ',')
csv_ctrl1$ech <- as.character(csv_ctrl1$ech)
#calculate fluorescence for a bin of 400 window
ctrl1_bin <- csv_ctrl1 %>% mutate(bin = GRN.B.HLog %/% 400)
sample = ctrl1_bin[1,1]

ctrl1_allbins <- full_join(ctrl1_bin, all_bins, by = 'bin') %>% mutate(ech = ifelse(is.na(ech), sample, ech)) %>% mutate(GRN.B.HLog = ifelse(is.na(ech), 0, GRN.B.HLog))

test0 <- rbind(test0,ctrl1_allbins %>% group_by(ech, bin, GRN.B.HLog) %>% dplyr::summarise(count = length(GRN.B.HLog)) %>% as.data.frame)
test_ctrl <- rbind(test0,ctrl1_allbins %>% group_by(ech, bin, GRN.B.HLog) %>% dplyr::summarise(count = length(GRN.B.HLog)) %>% as.data.frame)
head(test0)

###############
# Samples     #
###############

setwd("C:/Users/marsit/Dropbox/Ploidy-MA-CRL/2018-12-13_MA-Ploidy-CxB-11-bad")
dossier = strsplit(getwd(), split = "/" )[[1]] %>% tail(n=1)

#read csv files and calculate fluorescence for a bin of 400 
file_csv = dir(pattern='csv')

for(i in 1:length(file_csv)){
  #i=5
  print(i)
  csv_ech = read.table(file_csv[i], header = T, sep = ',')
  ech_bin <- csv_ech %>% mutate(bin = GRN.B.HLog %/% 400)
  file_name = strsplit(file_csv[i], '.', fixed = T)[[1]][1]
  strain_nb = strsplit(file_name,split = '_', fixed = T)[[1]] %>% tail(n=1)
  
  ech_allbins <- full_join(ech_bin, all_bins, by = 'bin') %>% mutate(ech = strain_nb) %>% mutate(GRN.B.HLog = ifelse(is.na(ech), 0, GRN.B.HLog))
  
  test0 <- rbind(test0,ech_allbins %>% group_by(ech, bin, GRN.B.HLog) %>% dplyr::summarise(count = length(GRN.B.HLog)) %>% as.data.frame)
  
}

test0$ech = as.factor(test0$ech)
head(test0)
tail(test0)
nrow(test0)


#################################################################################################################
# register table with fluorescence data only needed informations
write.table(test0, file="C:/Users/marsit/Dropbox/Ploidy-MA-CRL/2018-12-13_MA-Ploidy-CxB-11-bad/tab_heatmap_CxB-11-bad.txt", quote = FALSE, sep = "\t", row.names = FALSE )

##################################################################
#G1 ploidy heatmap for all crosses at different times-fig3-A     #
##################################################################

setwd("C:/Users/marsit/Desktop/Dropbox/Ploidy-MA-CRL")

#read tables with fluorescence data for each cross
table_PloidD_0 = read.table("C:/Users/marsit/Desktop/Dropbox/Ploidy-MA-CRL/MA_paper_fig3_data/tab_heatmap_D_t0.txt",header=T,sep="\t") %>% mutate(replicate="1", cross="L_Tini") 
table_PloidD_f = read.table("C:/Users/marsit/Desktop/Dropbox/Ploidy-MA-CRL/MA_paper_fig3_data/tab_heatmap_D_tf.txt",header=T,sep="\t") %>% mutate(replicate="1", cross="L_Tend") 
table_PloidB_0 = read.table("C:/Users/marsit/Desktop/Dropbox/Ploidy-MA-CRL/MA_paper_fig3_data/tab_heatmap_B_t0.txt",header=T,sep="\t") %>% mutate(replicate="1", cross="M_Tini")  
table_PloidB_f = read.table("C:/Users/marsit/Desktop/Dropbox/Ploidy-MA-CRL/MA_paper_fig3_data/tab_heatmap_B_tf.txt",header=T,sep="\t") %>% mutate(replicate="2", cross="M_Tend") 
table_PloidF_0 = read.table("C:/Users/marsit/Desktop/Dropbox/Ploidy-MA-CRL/MA_paper_fig3_data/tab_heatmap_F_t0.txt",header=T,sep="\t") %>% mutate(replicate="2", cross="H_Tini") 
table_PloidF_f = read.table("C:/Users/marsit/Desktop/Dropbox/Ploidy-MA-CRL/MA_paper_fig3_data/tab_heatmap_F_tf.txt",header=T,sep="\t") %>% mutate(replicate="2", cross="H_Tend")  
table_PloidI_0 = read.table("C:/Users/marsit/Desktop/Dropbox/Ploidy-MA-CRL/MA_paper_fig3_data/tab_heatmap_I_t0.txt",header=T,sep="\t") %>% mutate(replicate="1", cross="VL_Tini")
table_PloidI_f = read.table("C:/Users/marsit/Desktop/Dropbox/Ploidy-MA-CRL/MA_paper_fig3_data/tab_heatmap_I_tf.txt",header=T,sep="\t") %>% mutate(replicate="2", cross="VL_Tend")
table_Ploid_ctl = read.table("C:/Users/marsit/Desktop/Dropbox/Ploidy-MA-CRL/MA_paper_fig3_data/Heatmap_superposed/tab_heatmap_ctl.txt",header=T,sep="\t") 
table_Ploid_ctl = filter (table_Ploid_ctl, ech == "ctl3" | ech == "Haplo")

##############################################
#group all tables (all crosses, all timings) in one
Tab_P_heatmap= rbind (table_PloidD_0, table_PloidD_f, table_PloidB_0, table_PloidB_f, table_PloidF_0, table_PloidF_f, table_PloidI_0, table_PloidI_f) %>% mutate(cross = ifelse(ech == "ctl3" | ech == "Haplo", "c", cross))

######################################
## To keep only G1 peak fluorescence #
######################################

head(Tab_P_heatmap)
Tab_P_heatmap = select (Tab_P_heatmap, -replicate)

#separate fluorescence values corresonding to 2n peaks from 4n peaks
test2n = filter (Tab_P_heatmap, GRN.B.HLog < 7800)
test4n = filter (Tab_P_heatmap, GRN.B.HLog > 7800)

#calculate maximum count cell for 2n fluorescence values
test_max2n = ddply(test2n, .(ech, Hyb), subset, count==max(count))

### Select samples different from 2n 
test_max2n0 = filter(test_max2n, count<20) # ploidy= 3n & 4n
pas2n = test_max2n0 %>% select(ech, Hyb, Generations, cross)
testpas2n = left_join(pas2n, test4n)
test_max4n = ddply(testpas2n, .(ech, Hyb), subset,  count==max(count))

### Select diploid samples 2n 
test_max2ntrue= filter (test_max2n, count>20) # true diploids

###keep the 4n control
test_max4nC= ddply(test4n, .(ech, Hyb), subset, count==max(count))

max4nC = test_max4nC %>% select(ech, Hyb, Generations, cross, count, GRN.B.HLog, bin)
test_max4nctl = filter(max4nC, ech== "Haplo", Hyb== "t0")

test_maxtot2n = full_join(test_max2ntrue, test_max4nctl)

#calculate the maximum cell count for each sample
GRN.B.HLog_max2n = test_maxtot2n %>% group_by(ech, Hyb, Generations, cross) %>% dplyr::summarise(GRN.B.HLog_max2n = max(GRN.B.HLog)) %>% as.data.frame()
GRN.B.HLog_max4n = test_max4n %>% group_by(ech, Hyb, Generations, cross) %>% dplyr::summarise(GRN.B.HLog_max4n = max(GRN.B.HLog)) %>% as.data.frame()

head (test_maxtot2n)
nrow(test_maxtot2n)

Tab_P_heatmap_max2n = left_join(GRN.B.HLog_max2n %>% mutate(GRN.B.HLog = GRN.B.HLog_max2n), test_maxtot2n) %>% unique() %>% mutate(bin_max = bin)
Tab_P_heatmap_max4n = left_join(GRN.B.HLog_max4n %>% mutate(GRN.B.HLog = GRN.B.HLog_max4n), test_max4n) %>% unique() %>% mutate(bin_max = bin)
Tab_P_heatmap_max2n = select (Tab_P_heatmap_max2n, -GRN.B.HLog_max2n)
Tab_P_heatmap_max4n = select (Tab_P_heatmap_max4n, -GRN.B.HLog_max4n)

Tab_P_heatmap_max = rbind (Tab_P_heatmap_max2n , Tab_P_heatmap_max4n)

head(Tab_P_heatmap_max)
Tab_P_heatmap$ech = as.character(Tab_P_heatmap$ech)
str(Tab_P_heatmap)

#Keep fluorescence values 5 bins +/- around the maximum cell count
my_table = NULL
for(i in 1:nrow(Tab_P_heatmap_max)){
  #i=5
  print(i)
  #ligne = 1
  # les infos dont j'ai besoin
  ech = as.character(Tab_P_heatmap_max[i,1])
  Hyb = as.character(Tab_P_heatmap_max[i,2])
  cross = as.character(Tab_P_heatmap_max[i,4])
  bin_max = Tab_P_heatmap_max[i,8]
  my_bins = seq(bin_max - 5, bin_max + 5, by = 1)
  
  truc = Tab_P_heatmap %>% filter(bin %in% my_bins & ech == !!ech & Hyb == !!Hyb & cross == !!cross)
  my_table = rbind(my_table,truc)
}

#order the cross lines and timing
my_table$cross_f = factor(my_table$cross, levels=c("c", "VL_Tini", "VL_Tend", "L_Tini", "L_Tend", "M_Tini", "M_Tend", "H_Tini", "H_Tend"))

###########################
head(Tab_P_heatmap)

#generate a table for fluorescence values and a table for the maximum cell count
my_table_mean = my_table %>% group_by(ech, Hyb, Generations, cross, cross_f, bin) %>% dplyr::summarise(count_mean = mean(count)) %>% as.data.frame()
my_table_max = ddply(my_table_mean, .(ech, Hyb), subset,  count_mean==max(count_mean))

write.table(my_table_mean, file="C:/Users/marsit/Desktop/Dropbox/Ploidy-MA-CRL/fertility_exp2/Tab_heatmap_mean_count_bin_03_2019.txt", quote = FALSE, sep = "\t", row.names = FALSE )

head(Tab_P_heatmap)
summary(Tab_P_heatmap)

#### Figure 3A ##########

blue= rgb(0,0,(139/255))
red= rgb((238/255),0,0)
green= rgb(0,(139/255),0)
black= rgb(0,0,0,)

pdf("Heatmap_Ploid_fig3A.pdf", width=4, height= 6)
fig3A <- 
  ggplot(my_table_mean , aes(bin, ech)) + 
  geom_raster(aes(fill = count_mean)) + scale_fill_gradient(low = 'white', high = 'royalblue4', name="Cell count") +
  geom_point(data = my_table_max, mapping = aes(x=bin, y=ech), size=0.5) +
  ylab("Strain") +
  xlab("Fluorescence")+
  facet_grid(cross_f~., scales= "free", space = "free")+ theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), strip.text = element_text(colour = 'white', face="bold", size = 9))+
  xlim(9, 40)

g <- ggplot_gtable(ggplot_build(fig3A))
stripr <- which(grepl('strip-r', g$layout$name))
fills <- c("grey", red, red , blue ,blue, green, green, black ,black)

k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g)

dev.off()


##################################################################
#G1 ploidy points for different crosses at different times-fig3-B#
##################################################################

setwd("C:/Users/marsit/Dropbox/Ploidy-MA-CRL")

#read tables with fluorescence data for each cross generated for the previous script of Figure 3A (the same data). Here we represent the maximum cell counts normalized
table_fert_PloidA = read.table("C:/Users/marsit/Desktop/Dropbox/Ploidy-MA-CRL/MA_paper_fig3_data/tab_ploid_fert_A.txt",header=T,sep="\t") %>% mutate(replicate="1", cross="L")  %>% filter(Generations!="ctrl")
table_fert_PloidB = read.table("C:/Users/marsit/Desktop/Dropbox/Ploidy-MA-CRL/MA_paper_fig3_data/tab_ploid_fert_B.txt",header=T,sep="\t") %>% mutate(replicate="1", cross="M")  %>% filter(Generations!="ctrl") 
table_fert_PloidC = read.table("C:/Users/marsit/Desktop/Dropbox/Ploidy-MA-CRL/MA_paper_fig3_data/tab_ploid_fert_C.txt",header=T,sep="\t") %>% mutate(replicate="1", cross="H")  %>% filter(Generations!="ctrl")
table_fert_PloidD = read.table("C:/Users/marsit/Desktop/Dropbox/Ploidy-MA-CRL/MA_paper_fig3_data/tab_ploid_fert_D.txt",header=T,sep="\t") %>% mutate(replicate="2", cross="L")  %>% filter(Generations!="ctrl")
table_fert_PloidE = read.table("C:/Users/marsit/Desktop/Dropbox/Ploidy-MA-CRL/MA_paper_fig3_data/tab_ploid_fert_E.txt",header=T,sep="\t") %>% mutate(replicate="2", cross="M")  %>% filter(Generations!="ctrl")%>% filter(ech!="88")
table_fert_PloidF = read.table("C:/Users/marsit/Desktop/Dropbox/Ploidy-MA-CRL/MA_paper_fig3_data/tab_ploid_fert_F.txt",header=T,sep="\t") %>% mutate(replicate="2", cross="H")  %>% filter(Generations!="ctrl")
table_fert_PloidH = read.table("C:/Users/marsit/Desktop/Dropbox/Ploidy-MA-CRL/MA_paper_fig3_data/tab_ploid_fert_H.txt",header=T,sep="\t") %>% mutate(replicate="1", cross="VL")
table_fert_PloidI = read.table("C:/Users/marsit/Desktop/Dropbox/Ploidy-MA-CRL/MA_paper_fig3_data/tab_ploid_fert_I.txt",header=T,sep="\t") %>% mutate(replicate="2", cross="VL")

# bind all tables in one
Tab_fert_Ploid= rbind (table_fert_PloidA, table_fert_PloidB, table_fert_PloidC, table_fert_PloidD, table_fert_PloidE, table_fert_PloidF, table_fert_PloidH, table_fert_PloidI) %>% mutate(cross = ifelse(Generations == "ctrl", "c", cross))

head(Tab_fert_Ploid)

#reorder samples for the figure
Tab_fert_Ploid$cross_f = factor(Tab_fert_Ploid$cross, levels=c("c", "VL", "L", "M", "H"))

Tab_fert_Ploid$Generations_f = factor(Tab_fert_Ploid$Generations, levels=c("ctrl", "0", "~500", "~1000"))
Tab_fert_Ploid= filter(Tab_fert_Ploid, Generations!="NA")

head(Tab_fert_Ploid)

#######Figure 3B##########

blue_GC= rgb(0,0,(139/255))
red_GC= rgb((238/255),0,0)
green_GC= rgb(0,(139/255),0)
black_GC= rgb(0,0,0,)

pdf("plot_Ploid_point_fig3B.pdf", width=12, height= 4)

fig3B <- 
  ggplot(Tab_fert_Ploid[!is.na(Tab_fert_Ploid$Generations_f),], aes(x = Generations_f, y = as.numeric(normalisation), group=ech, group=Generations)) + 
  geom_line(alpha=2/10) + 
  geom_point(size=2.5, fill = "black") + 
  #geom_point(aes(color=fert), alpha = 0.85, size = 2.5) + 
  #scale_color_gradient2(midpoint=0.5, low="darkblue", mid="yellow", high="red", space ="Lab", name = "Fertility") + 
  scale_x_discrete(labels = c("ctrl" = "ctrl", "0"="Tini", "~500"= "Tmid", "~1000"="Tend", size=18 )) + 
  ylab("Ploidy")+ theme(axis.text=element_text(size=18))+
  xlab("")+
  facet_rep_grid(.~cross_f, scales= "free", space= "free")+ theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(strip.text = element_text(colour = 'white', face="bold", size = 18))+ theme(axis.title.y=element_text(size=18))+
  theme(axis.text.x=element_text(size=18, colour = "black"))+ theme(axis.text.y=element_text(size=18, colour = "black"))
ylim(1.5,4.5)

g <- ggplot_gtable(ggplot_build(fig3B))
stripr <- which(grepl('strip-t', g$layout$name))
fills <- c("grey", red_GC, blue_GC, green_GC, black_GC)

k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

grid.draw(g) 

dev.off()

#################################################################################################################
#################################################################################################################
##########################################################################
# GBS data analysis after generating the VCF merged files for figure 3 C #
##########################################################################

###########Analysis of VCF merge file to get the Allele frequency

library(dplyr)
library(reshape2)
library(ggplot2)
library(plyr)

#the path to the file
setwd("/Users/souhirmarsit/dropbox/Ploidy/GBSdata")

#read the VCf file for each cross (here exemple of L cross)
d = read.table("merge_A_10DP_sepf_3AD.txt", header=T) %>% mutate_all(as.character)

str(d)
names (d)

#####filter for heterozygous variants
d2 = mutate(d, "filtre_parent"=
              ifelse (LL11004.GT == MSH604.GT, "fail",
                      ifelse (LL11004.GT == "./." | MSH604.GT == "./.", "fail",
                              ifelse (LL11004.GT != "1/1" & LL11004.GT != "0/0" & LL11004.GT != "2/2", "fail", 
                                      ifelse (MSH604.GT != "1/1" & MSH604.GT != "0/0" & MSH604.GT != "2/2", "fail", "ok"))))) 
head(d2)
str(d2)

# change from character to facteur
d2$filtre_parent=as.factor(d2$filtre_parent)
summary(d2)

#filter d2 to keep the heterezygous sites between parents
d3= filter (d2, filtre_parent == "ok")
head(d3)
summary(d3)

#identify the parental alleles

dtest= select (d3, CHROM , POS, REF, ALT, starts_with ("LL11004"), starts_with ("MSH604"), starts_with ("A"))

dtestlong = melt (dtest, id.vars = c("CHROM", "POS", "REF", "ALT"))

dtestvar= mutate(dtestlong, variable = gsub("AD.", "AD_", variable), strain = gsub("(^.*)[.].*$", "\\1", variable), var = gsub (".*[.](.*)", "\\1", variable))
head(dtestvar)

dtest=dcast(dtestvar, CHROM + POS + REF + ALT+ strain ~ var)

dtestallele= mutate(dtest, GT.AD_1 = gsub ("(.*)[/].*", "\\1", GT), GT.AD_2 = gsub (".*[/](.*)", "\\1", GT), GT.AD_0 = "A" )
d_REF = dtestallele %>% filter(strain == "LL11004" | strain == "MSH604") %>%  select(-REF, -ALT, -starts_with("AD"), -GT, -GT.AD_2, -GT.AD_0) %>% dcast(CHROM + POS ~ strain)
head(d_REF)

dalBC= full_join(d_REF, dtestallele, by = c("CHROM", "POS"))

dalleleBC = mutate(dalBC, GT.AD_1=
                     ifelse (GT.AD_1 == LL11004, "C",
                             ifelse (GT.AD_1 == MSH604, "B", "MUT")),
                   GT.AD_2 = ifelse (GT.AD_2 == LL11004, "C",
                                     ifelse (GT.AD_2 == MSH604, "B", "MUT")))

head(dalleleBC)

AD_0 = dalleleBC %>% select(CHROM, POS, strain, DP, AD, GT.AD_0) %>% melt(id.vars = c("CHROM", "POS", "strain", "DP", "GT.AD_0")) %>% mutate(allele = GT.AD_0, AllelicDepth = value) %>% select(-GT.AD_0, -value)
head(AD_0)

AD_1 = dalleleBC %>% select(CHROM, POS, strain, DP, AD_1, GT.AD_1) %>% melt(id.vars = c("CHROM", "POS", "strain", "DP", "GT.AD_1")) %>% mutate(allele = GT.AD_1, AllelicDepth = value) %>% select(-GT.AD_1, -value)
head(AD_1)

AD_2 = dalleleBC %>% select(CHROM, POS, strain, DP, AD_2, GT.AD_2) %>% melt(id.vars = c("CHROM", "POS", "strain", "DP", "GT.AD_2")) %>% mutate(allele = GT.AD_2, AllelicDepth = value) %>% select(-GT.AD_2, -value) %>% mutate(AllelicDepth = ifelse(AllelicDepth < 0, "NA", AllelicDepth), allele = ifelse( AllelicDepth == "NA", "NA", allele))
head(AD_2)

MONSUPERTABLEAU = rbind(AD_0, AD_1, AD_2)

#register the table with read depth of identified parental alleles
write.table (MONSUPERTABLEAU, "/Users/souhirmarsit/dropbox/Ploidy/GBSdata/tableauA.txt", row.names = F, quote = F, sep= "\t")
head(MONSUPERTABLEAU)
MONSUPERTABLEAU$DP=as.integer(MONSUPERTABLEAU$DP)
MONSUPERTABLEAU$AllelicDepth=as.integer(MONSUPERTABLEAU$AllelicDepth)

#calculate allele frequency
tFA = mutate(MONSUPERTABLEAU, "FA"= (AllelicDepth/DP))
head (tFA)
# filter for SpC parent alleles
tFA_f = filter (tFA, allele == "C")

head(tFA_f)
summary(tFA_f)

#register the table with allele frequency of SpC parent alleles
write.table (MONSUPERTABLEAU, "/Users/souhirmarsit/dropbox/Ploidy/GBSdata/tableau-xA1-spC.txt", row.names = F, quote = F, sep= "\t")

##################################################################
#Allele frequency heatmap-fig3-C#
##################################################################

library(ggplot2)
library(cowplot)
library(lattice)
library(gridExtra)
library(dplyr)
library(plyr)
library(lemon)

setwd("C:/Users/marsit/Desktop/Dropbox/Ploidy/GBSdata-4n")
#read tables with allele frequency data for each cross generated in the script above
table_GBSA1 = read.table("C:/Users/marsit/Desktop/Dropbox/Ploidy-MA-CRL/MA_paper_fig3_data/tableau-xA1-spC.txt",header=T,sep="\t") %>% mutate(cross="A1")  
table_GBSA = read.table("C:/Users/marsit/Desktop/Dropbox/Ploidy-MA-CRL/MA_paper_fig3_data/tableau-xA-spC.txt",header=T,sep="\t") %>% mutate(cross="A",) 
table_GBSD1 = read.table("C:/Users/marsit/Desktop/Dropbox/Ploidy-MA-CRL/MA_paper_fig3_data/tableau-xD1-spC.txt",header=T,sep="\t") %>% mutate(cross="D1")
table_GBSD = read.table("C:/Users/marsit/Desktop/Dropbox/Ploidy-MA-CRL/MA_paper_fig3_data/tableau-xD-spC.txt",header=T,sep="\t") %>% mutate(cross="D")
table_ploidy = read.table("C:/Users/marsit/Desktop/Dropbox/Ploidy-MA-CRL/MA_paper_fig3_data/GBS_strain_ploidy.txt",header=T,sep="\t") 

#bind all tables in one
Tab_AF= rbind (table_GBSA1, table_GBSA, table_GBSD1, table_GBSD)
Tab_AF1 <- full_join(Tab_AF, table_ploidy , by = 'strain') 

#select a sample of strains to show (all 5 tetraploids in L crosses and a subset (6) of diploids and triploids)
Tab_AF2=filter(Tab_AF1, strain == "D40"|strain == "D61" |strain == "D84"|strain == "D85"|strain == "D89"|strain == "D93"|strain == "D75"|strain == "D42"|strain == "D36_P0"|strain == "D86"|strain == "A51_t0"|strain == "A87_t0"|strain == "A74"|strain == "A67"|strain == "A51_t35"|strain == "A87_P35"|strain == "D36_P35"|strain == "D45_P35"|strain == "LL11009")

#separate data according to the different levels of ploidy
Tab_par = filter(Tab_AF2, strain== "LL11009")

Tab_N = filter(Tab_AF2, Ploidy== 2)

Tab_3N= filter(Tab_AF2, Ploidy== 3)

Tab_4N= filter(Tab_AF2, Ploidy== 4)

#reorder the ploidy levels
Tab_AF2$Ploidy_f = factor(Tab_AF2$Ploidy, levels=c("1", "2", "3", "4"))

####### Figure 3C ############
pdf("Heatmap-GBS-fig3C.pdf", width=6, height= 5)  

fig3C <-
  g_ref <- ggplot(Tab_par,aes(x= interaction(POS, CHROM), y= strain)) +
  geom_raster(aes (fill=FA)) + scale_x_discrete(breaks = NULL)+ scale_y_discrete(breaks = NULL)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5)) +
  scale_fill_gradient2(low = "darkred", high = "midnightblue", mid = "white", midpoint = 0.5) + 
  labs(x=NULL) + labs(y=NULL) + theme(axis.ticks=element_blank()) + theme(axis.text.x=element_blank()) + theme(axis.text.y=element_blank())+ theme(legend.position= "none")+ theme(plot.margin= margin(0,0,0,2))

g_2n <- ggplot(Tab_N,aes(x= interaction(POS, CHROM), y= strain)) +
  geom_raster(aes (fill=FA)) + scale_x_discrete(breaks = NULL)+ scale_y_discrete(breaks = NULL)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5)) +
  scale_fill_gradient2(low = "darkred", high = "midnightblue", mid = "white", midpoint = 0.5) +
  labs(x=NULL) + labs(y=NULL) + theme(axis.ticks=element_blank()) + theme(axis.text.x=element_blank()) + theme(axis.text.y=element_blank())+ theme(legend.position= "none")+ theme(plot.margin= margin(0,0,0,2))

g_3n <- ggplot(Tab_3N,aes(x= interaction(POS, CHROM), y= strain)) +
  geom_raster(aes (fill=FA)) + scale_x_discrete(breaks = NULL)+ scale_y_discrete(breaks = NULL)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5)) +
  scale_fill_gradient2(low = "darkred", high = "midnightblue", mid = "white", midpoint = 0.5, name = "AF") +
  labs(x=NULL) + labs(y=NULL)+ theme(axis.ticks=element_blank()) + theme(axis.text.x=element_blank()) + theme(axis.text.y=element_blank())+ theme(plot.margin= margin(0,0,0,2))

g_4n <- ggplot(Tab_4N,aes(x= interaction(POS, CHROM), y= strain)) +
  geom_raster(aes (fill=FA)) + scale_x_discrete(breaks = NULL)+ scale_y_discrete(breaks = NULL)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5)) +
  scale_fill_gradient2(low = "darkred", high = "midnightblue", mid = "white", midpoint = 0.5) +
  labs(x=NULL)+ labs(y=NULL) + theme(axis.ticks=element_blank()) + theme(axis.text.x=element_blank()) + theme(axis.text.y=element_blank())+ theme(legend.position= "none")+ theme(plot.margin= margin(0,0,0,2))

g_3n_a <- g_3n + theme(legend.position= "none")

cowplot::plot_grid(g_ref, g_2n, g_3n_a, g_4n, nrow=4, rel_heights = c(50, 300, 300, 250)) %>% 
  cowplot::plot_grid(cowplot::get_legend(g_3n), rel_widths= c(1, 0.2))

dev.off()

#############################################################################################################
##############################################################################
#Allele frquency analysis- Whole genome sequencing -Figure3D and Figure S16  #
##############################################################################
###########Analyse VCF merge and generate tables with allele frequency of all variants corresponding to hybrid parent 2 for each cross (here exemple for L cross)
library(dplyr)
library(reshape2)
library(ggplot2)
library(plyr)

setwd("/Users/souhirmarsit/Desktop/MA_Sequences")

#Read VCF merged table for each cross (here exemle for L cross)
d = read.table("merge_SpB_A_vcf_sep.txt", header=T, sep='\t') %>% mutate_all(as.character)
d = select(d, CHROM, POS, REF, ALT, contains("A31"), contains("A51"), contains("A87"), contains("LL2011_004"), contains("MSH.604"))

str(d)
names (d)

#####filtrer for heterozygous variants between the two hybrid parents

d2 = mutate(d, "filtre_parent"=
              ifelse (LL2011_004.GT == MSH.604.GT, "fail",
                      ifelse (LL2011_004.GT == "./." | 	MSH.604.GT == "./.", "fail", "ok"))) 
head(d2)
str(d2)

d2$filtre_parent=as.factor(d2$filtre_parent)
summary(d2)

#filter d2 to keep only heterozygous variants

d3= filter (d2, filtre_parent == "ok")
head(d3)
summary(d3)

dtest= select (d3, CHROM , POS, REF, ALT, starts_with ("LL2011_004"), starts_with ("MSH.604"), starts_with ("A"))

#Identify parental alleles

dtestlong = melt (dtest, id.vars = c("CHROM", "POS", "REF", "ALT"))

dtestvar= mutate(dtestlong, variable = gsub("AD.", "AD_", variable), strain = gsub("(^.*)[.].*$", "\\1", variable), var = gsub (".*[.](.*)", "\\1", variable))
head(dtestvar)

dtest=dcast(dtestvar, CHROM + POS + REF + ALT+ strain ~ var)

dtestallele= mutate(dtest, GT.AD_1 = gsub ("(.*)[/].*", "\\1", GT), GT.AD_2 = gsub (".*[/](.*)", "\\1", GT), GT.AD_0 = "B")
d_REF = dtestallele %>% filter(strain == "LL2011_004") %>%  select(-REF, -ALT, -starts_with("AD"), -GT, -GT.AD_2, -GT.AD_0) %>% dcast(CHROM + POS ~ strain)
head(d_REF)

dalBC= full_join(d_REF, dtestallele, by = c("CHROM", "POS"))

dalleleBC = mutate(dalBC, GT.AD_1= "C", GT.AD_2 = "MUT")

head(dalleleBC)

AD_0 = dalleleBC %>% select(CHROM, POS, strain, DP, AD, GT.AD_0) %>% melt(id.vars = c("CHROM", "POS", "strain", "DP", "GT.AD_0")) %>% mutate(allele = GT.AD_0, AllelicDepth = value) %>% select(-GT.AD_0, -value)
head(AD_0)

AD_1 = dalleleBC %>% select(CHROM, POS, strain, DP, AD_1, GT.AD_1) %>% melt(id.vars = c("CHROM", "POS", "strain", "DP", "GT.AD_1")) %>% mutate(allele = GT.AD_1, AllelicDepth = value) %>% select(-GT.AD_1, -value)
head(AD_1)

AD_2 = dalleleBC %>% select(CHROM, POS, strain, DP, AD_2, GT.AD_2) %>% melt(id.vars = c("CHROM", "POS", "strain", "DP", "GT.AD_2")) %>% mutate(allele = GT.AD_2, AllelicDepth = value) %>% select(-GT.AD_2, -value) %>% mutate(AllelicDepth = ifelse(AllelicDepth < 0, "NA", AllelicDepth), allele = ifelse( AllelicDepth == "NA", "NA", allele))
head(AD_2)

MONSUPERTABLEAU = rbind(AD_0, AD_1, AD_2)

#Identification of contigs -> chromosome
ref_conversion <- read.table("SpB_tig_rearrangement.txt", header = F, col.names = c("Chr", "CHROM", "Orientation"))
head (ref_conversion)
SUPERTABLEAU <- left_join(MONSUPERTABLEAU, ref_conversion, by = "CHROM")
head (SUPERTABLEAU)
SUPERTABLEAU$Chr_ord <- factor(SUPERTABLEAU$Chr, levels = c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII", "chrIX", "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI"))

#Identification of strains name
Info_strain <- read.table("Table_info_A_LOH.txt", header = T, sep = "\t", fill = TRUE)
head (Info_strain)

SUPERTABLEAU$DP=as.integer(MONSUPERTABLEAU$DP)
SUPERTABLEAU$AllelicDepth=as.integer(MONSUPERTABLEAU$AllelicDepth)

#filter for variants with more than 20 reads 
tabDP20 = filter (SUPERTABLEAU, DP>20)

#Calculate allele frequency the ratio of allelic depth and the read depth of the locus   
tFA = mutate(tabDP20, "FA"= (AllelicDepth/DP))
head (tFA)

#Register the table of allele frequency of all variants 
write.table (tFA, "/Users/souhirmarsit/Desktop/MA_Sequences/LOH_4n/tableau_FA_A_31_51_87_all.txt", row.names = F, quote = F, sep= "\t")
#filter for the Hybrid parent 2 (SpC or SpA or S. cerevisiae)
tFA_f = filter (tFA, allele == "C")
tFA_f <- left_join(tFA_f, Info_strain, by = "strain")
head(tFA_f)
summary(tFA_f)
#Register the table of allele frequency corresponding to the hybrid parent2
write.table (tFA_f, "/Users/souhirmarsit/Desktop/MA_Sequences/LOH_4n/tableau_FA_A_31_51_87.txt", row.names = F, quote = F, sep= "\t")


####################figure3D script
library(dplyr)
library(reshape2)
library(ggplot2)
library(plyr)
library(cowplot)
library(lattice)
library(gridExtra)
library(lemon)
library(grid)

setwd("/Users/souhirmarsit/Desktop/MA_Sequences/LOH_4n")

I1 = read.table("name_strain_A.txt", header=T) 
I2 = read.table("name_strain_D.txt", header=T) 
I3 = read.table("name_strain_B.txt", header=T) 
I4 = read.table("name_strain_F.txt", header=T) 

d1 = read.table("table_FA_A_31_51_87.txt", header=T) 
d2 = read.table("table_FA_D_36_45_all.txt", header=T) 
d3 = read.table("table_FA_B_40_all.txt", header=T) 
d4 = read.table("table_FA_F_38_all.txt", header=T) 

d1 = select (d1, -strain_t, -time, -cross, -Ploidy, - Fertility)

d11 = full_join(d1, I1, by="strain")
d22 = full_join(d2, I2, by="strain")
d33 = full_join(d3, I3, by="strain")
d44 = full_join(d4, I4, by="strain")

d11 = filter (d11, variable_f != "NA")
d22 = filter (d22, variable_f != "NA")
d33 = filter (d33, variable_f != "NA")
d44 = filter (d44, variable_f != "NA")

dtot = rbind (d11, d22, d33, d44)
dtot$time_f = factor(dtot$time, levels=c('Tini', 'Tend'))

blue= rgb(0,0,(139/255))
red= rgb((238/255),0,0)
green= rgb(0,(139/255),0)
black= rgb(0,0,0,)

pdf(file = paste0("/Users/souhirmarsit/Desktop/MA_Sequences/LOH_4n/","LOH_4n_all_density.pdf"), height = 9, width = 6)

PlotA<- ggplot(dtot, aes(x = FA)) + geom_density(aes(linetype=time)) + scale_color_manual(values="black")+ 
  facet_grid(cross~., scales = "free") + 
  theme(axis.title.y = element_text(size=14), axis.text= element_text(size=16), strip.text.y= element_text(face="bold", color= "white", size=18))

g <- ggplot_gtable(ggplot_build(PlotA))
stripr <- which(grepl('strip-r', g$layout$name))
fills <- c(blue, blue, blue, blue, blue, green, green, black)

k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g)
dev.off()

############################################# Figure S16 ###############

d1_Tini = filter(d11, time == "Tini")
d1_Tend = filter(d11, time == "Tend")

d2_Tini = filter(d22, time == "Tini")
d2_Tend = filter(d22, time == "Tend")

d3_Tini = filter(d33, time == "Tini")
d3_Tend = filter(d33, time == "Tend")

d4_Tini = filter(d44, time == "Tini")
d4_Tend = filter(d44, time == "Tend")

d1_Tini$variable_ord = factor(d1_Tini$variable_f, levels=c("L1_87_Tini", "L1_51_Tini", "L1_31_Tini"))
d2_Tini$variable_ord = factor(d2_Tini$variable_f, levels=c("L2_45_Tini", "L2_36_Tini"))
d3_Tini$variable_ord = factor(d3_Tini$variable_f, levels=c("M1_49_Tini", "M1_40_Tini"))
d4_Tini$variable_ord = factor(d4_Tini$variable_f, levels=c("H2_38_Tini"))

d1_Tend$variable_ord = factor(d1_Tend$variable_f, levels=c("L1_87_Tend", "L1_51_Tend", "L1_31_Tend"))
d2_Tend$variable_ord = factor(d2_Tend$variable_f, levels=c("L2_45_Tend", "L2_36_Tend"))
d3_Tend$variable_ord = factor(d3_Tend$variable_f, levels=c("M1_49_Tend", "M1_40_Tend"))
d4_Tend$variable_ord = factor(d4_Tend$variable_f, levels=c("H2_38_Tend"))

pdf(file = paste0("/Users/souhirmarsit/Desktop/MA_Sequences/LOH_4n/","LOH_4n_density.pdf"), height = 10, width = 20)
jpeg("LOH_4n_nv.jpg", width = 1400, height = 350)
A_tini <- ggplot(d1_Tini, aes(x= interaction(POS, Chr), y= variable_ord)) + geom_raster(aes (fill=FA)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5)) + theme(axis.text.y = element_text( size = 18)) +
  scale_fill_gradient2(low = "darkred", high = "midnightblue", mid = "white", midpoint = 0.5, limits=c(0, 1), name= "Allele frequency")+ theme(axis.text.x=element_blank())+ labs(x="ChrIII Position") + labs(y=NULL) + labs(x=NULL)   
D_tini <- ggplot(d2_Tini, aes(x= interaction(POS, Chr), y= variable_ord)) + geom_raster(aes (fill=FA)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5)) + theme(axis.text.y = element_text( size = 18)) +
  scale_fill_gradient2(low = "darkred", high = "midnightblue", mid = "white", midpoint = 0.5, limits=c(0, 1))+ theme(axis.text.x=element_blank())+ theme(legend.position = "none") + labs(x="ChrIII Position") + labs(y=NULL) + labs(x=NULL) 
B_tini <- ggplot(d3_Tini, aes(x= interaction(POS, Chr), y= variable_ord)) + geom_raster(aes (fill=FA)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5)) + theme(axis.text.y = element_text( size = 18)) +
  scale_fill_gradient2(low = "darkred", high = "midnightblue", mid = "white", midpoint = 0.5, limits=c(0, 1))+ theme(axis.text.x=element_blank())+ theme(legend.position = "none") + labs(x="ChrIII Position") + labs(y=NULL) + labs(x=NULL) 
F_tini <- ggplot(d4_Tini, aes(x= interaction(POS, Chr), y= variable_ord)) + geom_raster(aes (fill=FA)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5)) + theme(axis.text.y = element_text( size = 18)) +
  scale_fill_gradient2(low = "darkred", high = "midnightblue", mid = "white", midpoint = 0.5, limits=c(0, 1))+ theme(axis.text.x=element_blank())+ theme(legend.position = "none") + labs(x="ChrIII Position") + labs(y=NULL) + labs(x=NULL) 

A_tend <- ggplot(d1_Tend, aes(x= interaction(POS, Chr), y= variable_ord)) + geom_raster(aes (fill=FA)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5)) + theme(axis.text.y = element_text( size = 18)) +
  scale_fill_gradient2(low = "darkred", high = "midnightblue", mid = "white", midpoint = 0.5, limits=c(0, 1))+ theme(axis.text.x=element_blank())+ theme(legend.position = "none") + labs(x="ChrIII Position") + labs(y=NULL) + labs(x=NULL) 
D_tend <- ggplot(d2_Tend, aes(x= interaction(POS, Chr), y= variable_ord)) + geom_raster(aes (fill=FA)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5)) + theme(axis.text.y = element_text( size = 18)) +
  scale_fill_gradient2(low = "darkred", high = "midnightblue", mid = "white", midpoint = 0.5, limits=c(0, 1))+ theme(axis.text.x=element_blank())+ theme(legend.position = "none") + labs(x="ChrIII Position") + labs(y=NULL) + labs(x=NULL) 
B_tend <- ggplot(d3_Tend, aes(x= interaction(POS, Chr), y= variable_ord)) + geom_raster(aes (fill=FA)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5)) + theme(axis.text.y = element_text( size = 18)) +
  scale_fill_gradient2(low = "darkred", high = "midnightblue", mid = "white", midpoint = 0.5, limits=c(0, 1))+ theme(axis.text.x=element_blank())+ theme(legend.position = "none") + labs(x="ChrIII Position") + labs(y=NULL) + labs(x=NULL) 
F_tend <- ggplot(d4_Tend, aes(x= interaction(POS, Chr), y= variable_ord)) + geom_raster(aes (fill=FA)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5)) + theme(axis.text.y = element_text( size = 18)) +
  scale_fill_gradient2(low = "darkred", high = "midnightblue", mid = "white", midpoint = 0.5, limits=c(0, 1))+ theme(axis.text.x=element_blank())+ theme(legend.position = "none") + labs(x="ChrIII Position") + labs(y=NULL) + labs(x=NULL) 

A_tini_a <- A_tini + theme(legend.position = "none")

cowplot::plot_grid(A_tini_a, A_tend, D_tini, D_tend, B_tini, B_tend, F_tini, F_tend, nrow=4, ncol=2, rel_heights= c(105, 75, 75, 42)) %>%
  cowplot::plot_grid(cowplot::get_legend(A_tini), rel_widths = c(1, 0.1))
dev.off()

###################################################################################
##################################################################
#G1 ploidy heatmap for SpC haploid parent stock Figure S9        #
##################################################################
library(ggplot2)
library(cowplot)
library(lattice)
library(gridExtra)
library(grid)
library(dplyr)
library(plyr)
library(lemon)
setwd("C:/Users/marsit/Desktop/Dropbox/Ploidy-MA-CRL/")   

#Read tables of fluorescence
table_Ploid_H = read.table("C:/Users/marsit/Desktop/Dropbox/Ploidy-MA-CRL/2018-11-28_MA-Ploidy-SpCparent3/tab_heatmap_SpCparent.txt",header=T,sep="\t") %>% mutate(cross="SpC rep1") 
table_Ploid_I = read.table("C:/Users/marsit/Desktop/Dropbox/Ploidy-MA-CRL/2019-01-16_MA-Ploidy-SpCparent4/tab_heatmap_SpCparent4.txt",header=T,sep="\t") %>% mutate(cross="SpC rep2") 
table_Ploid_ctl = read.table("C:/Users/marsit/Desktop/Dropbox/Ploidy-MA-CRL/tab_heatmap_ctrols.txt",header=T,sep="\t") %>% mutate(cross="ctl") 

##############################################

Tab_P_heatmap= rbind (table_Ploid_ctl, table_Ploid_H, table_Ploid_I)

Tab_P_heatmap= mutate(Tab_P_heatmap, cross = 
                        ifelse (ech == "Control1" | ech == "Control2" |ech == "Haplo", "c", 
                                ifelse(ech == "ctl3", "ctl", cross)))

###To take G1 peak and not G2 peak (1n peaks and 2n peaks)
######################################################
head(Tab_P_heatmap)
test1n = filter (Tab_P_heatmap, GRN.B.HLog < 4500)
test2n = filter (Tab_P_heatmap, GRN.B.HLog > 4500)
test_max1n = ddply(testn, .(ech, cross), subset, count==max(count))

###identify those without 2n or 3n G1 pic and remove them
test_max1n0 = filter(test_max1n, count<20) # = 2n or 3n or 4n

pas1n = test_max1n0 %>% select(ech, cross)
testpas1n = left_join(pas1n, test2n)
test_max2n = ddply(testpas1n, .(ech, cross), subset,  count==max(count))

test_max1ntrue= filter (test_max1n, count>20) # = 1n

###assemble 1n, 2n and 4n in G1 phase

###########to keep 4n controls
test_max2nC= ddply(test2n, .(ech, cross), subset, count==max(count))
max2nC = test_max2nC %>% select(ech, cross, count, GRN.B.HLog, bin)
test_max2nctl = filter(max2nC, ech== "Haplo", cross== "c")

test_maxtot1n = full_join(test_max1ntrue, test_max2nctl)

#calculate the max cell count
GRN.B.HLog_max1n = test_maxtot1n %>% group_by(ech, cross) %>% dplyr::summarise(GRN.B.HLog_max1n = max(GRN.B.HLog)) %>% as.data.frame()
GRN.B.HLog_max2n = test_max2n %>% group_by(ech, cross) %>% dplyr::summarise(GRN.B.HLog_max2n = max(GRN.B.HLog)) %>% as.data.frame()
head (test_maxtot1n)
nrow(test_maxtot1n)

Tab_P_heatmap_max1n = left_join(GRN.B.HLog_max1n %>% mutate(GRN.B.HLog = GRN.B.HLog_max1n), test_maxtot1n) %>% unique() %>% mutate(bin_max = bin)
Tab_P_heatmap_max2n = left_join(GRN.B.HLog_max2n %>% mutate(GRN.B.HLog = GRN.B.HLog_max2n), test_max2n) %>% unique() %>% mutate(bin_max = bin)
Tab_P_heatmap_max1n = select (Tab_P_heatmap_max1n, -GRN.B.HLog_max1n)
Tab_P_heatmap_max2n = select (Tab_P_heatmap_max2n, -GRN.B.HLog_max2n)
Tab_P_heatmap_max = rbind (Tab_P_heatmap_max1n , Tab_P_heatmap_max2n)
Tab_P_heatmap_max = filter (Tab_P_heatmap_max, cross== "c"|cross== "SpC rep1"| cross== "SpC rep2")
Tab_P_heatmap= filter (Tab_P_heatmap, cross== "c"|cross== "SpC rep1"| cross== "SpC rep2") 
Tab_P_heatmap = select (Tab_P_heatmap, -GRN.B.HLog)

head(Tab_P_heatmap_max)
Tab_P_heatmap$ech = as.character(Tab_P_heatmap$ech)
str(Tab_P_heatmap)
head(Tab_P_heatmap)

#Keep fluorescence values 5 bins around the max cell count
my_table = NULL
for(i in 1:nrow(Tab_P_heatmap_max)){
  print(i)
  ech = as.character(Tab_P_heatmap_max[i,1])
  cross = as.character(Tab_P_heatmap_max[i,2])
  bin_max = Tab_P_heatmap_max[i,6]
  my_bins = seq(bin_max -5, bin_max +5, by = 1)
  
  truc = Tab_P_heatmap %>% filter(bin %in% my_bins & ech == !!ech & cross == !!cross)
  my_table = rbind(my_table,truc)
}

#generate table for fluorescence values and a table for max cell count
my_table$cross = factor(my_table$cross, levels=c("c", "SpC rep1", "SpC rep2"))
my_table_mean = my_table %>% group_by(ech, cross, bin) %>% dplyr::summarise(count_mean = mean(count)) %>% as.data.frame()
my_table_max = ddply(my_table_mean, .(ech), subset,  count_mean==max(count_mean))
my_table_mean = mutate(my_table_mean, "fluo"= (bin*400))

head(Tab_P_heatmap)
summary(Tab_P_heatmap)

blue= rgb(0,0,(139/255))
red= rgb((238/255),0,0)
green= rgb(0,(139/255),0)
black= rgb(0,0,0,)
#### figure S9
pdf("Heatmap_figS9_revision.pdf", width=4, height= 7)

ggplot(my_table_mean , aes(fluo, ech)) + 
  geom_raster(aes(fill = count_mean)) + scale_fill_gradient(low = 'white', high = 'royalblue4', name='Cell count') +
  
  ylab("Strain") +
  xlab("Fluorescence(a.u.)")+
  facet_grid(cross~., scales= "free", space = "free")+ theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), strip.text = element_text(colour = 'black', face="bold.italic", size = 9))+
  xlim(2400, 16000)
dev.off()

pdf("Heatmap_Ploid_Point_fig_SpCf.pdf", width=4, height= 7)
plot_A 

grid.draw(g)
dev.off()

###################################################################################
########################################################################
#G1 ploidy heatmap for SpCxSpC, SpCxSpA, SpCxScer crosses Figure S11   #
########################################################################
library(ggplot2)
library(cowplot)
library(lattice)
library(gridExtra)
library(grid)
library(dplyr)
library(plyr)
library(lemon)
setwd("C:/Users/marsit/Desktop/Dropbox/Ploidy-MA-CRL/")   

#read tables with fluorescence data generated with the same sripts above for the other ploidy figures (Figure 3A, figure S9)
table_Ploid = read.table("C:/Users/marsit/Desktop/Dropbox/Ploidy-MA-CRL/Tab_fig_S11.txt",header=T,sep="\t")
table_Ploid = mutate (table_Ploid, "fluo" =(bin*400))

table_Ploid$cross_f = factor(table_Ploid$cross, levels=c('c','SpCxSpC','SpCxSpA','SpCxS.cer'))

pdf("Heatmap_figS11_revision_essai.pdf", width=4, height= 7)

ggplot(table_Ploid , aes(fluo, ech)) + 
  geom_raster(aes(fill = count_mean)) + scale_fill_gradient(low = 'white', high = 'royalblue4', name='Cell count') +
  
  ylab("Strain") +
  xlab("Fluorescence (a.u.)")+
  facet_grid(cross_f~., scales= "free", space = "free")+ theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), strip.text = element_text(colour = 'black', face="bold.italic", size = 9))+
  xlim(2400, 16000)
dev.off()


###################################################################################
########################################################################
#G1 ploidy heatmap for different SpCxSpB crosses  Figure S12           #
########################################################################
library(ggplot2)
library(cowplot)
library(lattice)
library(gridExtra)
library(grid)
library(dplyr)
library(plyr)
library(lemon)

setwd("C:/Users/marsit/Desktop/Dropbox/Ploidy-MA-CRL/")  
 
#read tables with fluorescence data generated with the same sripts described above for the other ploidy figures (Figure 3A, figure S9)
table_Ploid_1 = read.table("C:/Users/marsit/Desktop/Dropbox/Ploidy-MA-CRL/MA_paper_fig3_data/Tab_heatmap_FigS12_bin_L1.txt",header=T,sep="\t")
table_Ploid_2 = read.table("C:/Users/marsit/Desktop/Dropbox/Ploidy-MA-CRL/MA_paper_fig3_data/Tab_heatmap_FigS12_bin_L2.txt",header=T,sep="\t")
table_Ploid_3 = read.table("C:/Users/marsit/Desktop/Dropbox/Ploidy-MA-CRL/MA_paper_fig3_data/Tab_heatmap_FigS12_bin_L3.txt",header=T,sep="\t")
table_Ploid_4 = read.table("C:/Users/marsit/Desktop/Dropbox/Ploidy-MA-CRL/MA_paper_fig3_data/Tab_heatmap_FigS12_bin_L4.txt",header=T,sep="\t")
table_Ploid_5 = read.table("C:/Users/marsit/Desktop/Dropbox/Ploidy-MA-CRL/MA_paper_fig3_data/Tab_heatmap_FigS12_bin_L5.txt",header=T,sep="\t")
table_Ploid_6 = read.table("C:/Users/marsit/Desktop/Dropbox/Ploidy-MA-CRL/MA_paper_fig3_data/Tab_heatmap_FigS12_bin_L6.txt",header=T,sep="\t")

pdf("Heatmap_figS12.pdf", width=7, height= 10)

plot_1 <- 
  ggplot(table_Ploid_1 , aes(bin, ech)) + 
  geom_raster(aes(fill = Cell_count)) + scale_fill_gradient(low = 'white', high = 'royalblue4', name='Cell count') +
  
  ylab("Strain") + theme(plot.margin= margin(0,0,2,2)) +
  xlab("Fluorescence")+
  facet_grid(cross~., scales= "free", space = "free")+ theme_bw() + theme(legend.position= "none") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), strip.text = element_text(colour = 'black', face="bold", size = 11))+
  xlim(9, 40)

plot_2 <- 
  ggplot(table_Ploid_2 , aes(bin, ech)) + 
  geom_raster(aes(fill = Cell_count)) + scale_fill_gradient(low = 'white', high = 'royalblue4', name='Cell count') +
  
  labs(y=NULL) + theme(plot.margin= margin(0,0,2,2)) +
  xlab("Fluorescence")+
  facet_grid(cross~., scales= "free", space = "free")+ theme_bw() + theme(legend.position= "none") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), strip.text = element_text(colour = 'black', face="bold", size = 11))+
  xlim(9, 40)

plot_3 <- 
  ggplot(table_Ploid_3 , aes(bin, ech)) + 
  geom_raster(aes(fill = Cell_count)) + scale_fill_gradient(low = 'white', high = 'royalblue4', name='Cell count') +
  
  labs(y=NULL) + theme(plot.margin= margin(0,0,2,2)) +
  xlab("Fluorescence")+
  facet_grid(cross~., scales= "free", space = "free")+ theme_bw() + theme(legend.position= "none") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), strip.text = element_text(colour = 'black', face="bold", size = 11))+
  xlim(9, 40)

plot_4 <- 
  ggplot(table_Ploid_4 , aes(bin, ech)) + 
  geom_raster(aes(fill = Cell_count)) + scale_fill_gradient(low = 'white', high = 'royalblue4', name='Cell count') +
  
  ylab("Strain") + theme(plot.margin= margin(2,0,0,2)) +
  xlab("Fluorescence")+
  facet_grid(cross~., scales= "free", space = "free")+ theme_bw() + theme(legend.position= "none") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), strip.text = element_text(colour = 'black', face="bold", size = 11))+
  xlim(9, 40)

plot_5 <- 
  ggplot(table_Ploid_5 , aes(bin, ech)) + 
  geom_raster(aes(fill = Cell_count)) + scale_fill_gradient(low = 'white', high = 'royalblue4', name='Cell count') +
  
  labs(y=NULL) + theme(plot.margin= margin(2,0,0,2)) +
  xlab("Fluorescence")+
  facet_grid(cross~., scales= "free", space = "free")+ theme_bw() + theme(legend.position= "none") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), strip.text = element_text(colour = 'black', face="bold", size = 11))+
  xlim(9, 40)

plot_6 <- 
  ggplot(table_Ploid_6 , aes(bin, ech)) + 
  geom_raster(aes(fill = Cell_count)) + scale_fill_gradient(low = 'white', high = 'royalblue4', name='Cell count') +
  
  labs(y=NULL) + theme(plot.margin= margin(2,0,0,2)) +
  xlab("Fluorescence")+
  facet_grid(cross~., scales= "free", space = "free")+ theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), strip.text = element_text(colour = 'black', face="bold", size = 11))+
  xlim(9, 40)


plot_6_a <- plot_6 + theme(legend.position= "none")

cowplot::plot_grid(plot_1, plot_2, plot_3, plot_4, plot_5, plot_6_a, ncol=3, nrow=2) %>% 
  cowplot::plot_grid(cowplot::get_legend(plot_6), rel_widths= c(1, 0.2))

dev.off()

#############################################################################################################
##################################################################
#Allele frequency heatmap for H2 tetraploid hybrid-figS13        #
##################################################################

library(dplyr)
library(reshape2)
library(ggplot2)

setwd("C:/Users/marsit/Desktop/Dropbox/Ploidy/GBSdata-4n/")

#read vcf file merge for F cross, read mapped on S. cerevisiae (this exemple) or S. paradoxus reference genome
d = read.table("merge_Scer_10DP_sep.txt", header=T, sep="\t") %>% mutate_all(as.character)

str(d)
names (d)

#####filter to keep all heterozygous loci

d2 = mutate(d, "filtre_parent"=
              
              ifelse (F38_P0.Scer.GT == "./." | F38_P0.Scer.GT == "." | F38_P0.Scer.GT == "1/1" | F38_P0.Scer.GT == "0/0" | F38_P0.Scer.GT == "2/2", "fail", 
                      ifelse (F38_P35.Scer.GT == "./." | F38_P35.Scer.GT == "." | F38_P35.Scer.GT == "1/1" | F38_P35.Scer.GT == "0/0" | F38_P35.Scer.GT == "2/2", "fail",
                              ifelse (F50_P0.Scer.GT == "./." | F50_P0.Scer.GT == "." | F50_P0.Scer.GT == "1/1" | F50_P0.Scer.GT == "0/0" | F50_P0.Scer.GT == "2/2", "fail", 
                                      ifelse (F50_P35.Scer.GT == "./." | F50_P35.Scer.GT == "." | F50_P35.Scer.GT == "1/1" | F50_P35.Scer.GT == "0/0" | F50_P35.Scer.GT == "2/2", "fail","ok"))))) 
head(d2)
str(d2)

# change character to factor
d2$filtre_parent=as.factor(d2$filtre_parent)
summary(d2)

d3= filter (d2, filtre_parent == "ok")
head(d3)
summary(d3)

dtestlong = melt (d3, id.vars = c("CHROM", "POS", "REF", "ALT"))

dtestvar= mutate(dtestlong, variable = gsub("AD.", "AD_", variable), strain = gsub("(^.*)[.].*$", "\\1", variable), var = gsub (".*[.](.*)", "\\1", variable))
head(dtestvar)

dtest=dcast(dtestvar, CHROM + POS + REF + ALT+ strain ~ var)

dtestallele= mutate(dtest, GT.AD_1 = gsub ("(.*)[/].*", "\\1", GT), GT.AD_2 = gsub (".*[/](.*)", "\\1", GT), GT.AD_0 = "A" )

AD_0 = dtestallele %>% select(CHROM, POS, strain, DP, AD, GT.AD_0) %>% melt(id.vars = c("CHROM", "POS", "strain", "DP", "GT.AD_0")) %>% mutate(allele = GT.AD_0, AllelicDepth = value) %>% select(-GT.AD_0, -value)
head(AD_0)

AD_1 = dtestallele %>% select(CHROM, POS, strain, DP, AD_1, GT.AD_1) %>% melt(id.vars = c("CHROM", "POS", "strain", "DP", "GT.AD_1")) %>% mutate(allele = GT.AD_1, AllelicDepth = value) %>% select(-GT.AD_1, -value)
head(AD_1)

AD_2 = dtestallele %>% select(CHROM, POS, strain, DP, AD_2, GT.AD_2) %>% melt(id.vars = c("CHROM", "POS", "strain", "DP", "GT.AD_2")) %>% mutate(allele = GT.AD_2, AllelicDepth = value) %>% select(-GT.AD_2, -value) %>% mutate(AllelicDepth = ifelse(AllelicDepth < 0, "NA", AllelicDepth), allele = ifelse( AllelicDepth == "NA", "NA", allele))
head(AD_2)

MONSUPERTABLEAU = rbind(AD_0, AD_1, AD_2)
head(MONSUPERTABLEAU)

MONSUPERTABLEAU$DP=as.integer(MONSUPERTABLEAU$DP)
MONSUPERTABLEAU$AllelicDepth=as.integer(MONSUPERTABLEAU$AllelicDepth)

MONSUPERTABLEAU= filter(MONSUPERTABLEAU, DP >20)

#filter for all heterozygous alleles
AD1= filter(MONSUPERTABLEAU, variable == "AD_1")
tFA = mutate(AD1, "FA"= (AllelicDepth/DP))
head (tFA)
#register table with allele frequency of heterozygous loci
write.table(tFA, file="C:/Users/marsit/Desktop/Dropbox/Ploidy/GBS-FA-map-Scer.txt", quote = FALSE, sep = "\t", row.names = FALSE )

head(tFA_f)
summary(tFA_f)

##########################################script figure S13
library(ggplot2)
library(cowplot)
library(lattice)
library(gridExtra)
library(grid)
library(dplyr)
library(plyr)
library(lemon)
library(scales)
library(reshape2)

setwd("C:/Users/marsit/Desktop/Dropbox/Ploidy")

#Read tables of allele frequency for mapping on S. paradoxus or S. cerevisiae reference genome
table_GBS_Scer = read.table("C:/Users/marsit/Desktop/Dropbox/Ploidy/GBS-FA-map-Scer.txt",header=T,sep="\t")   
table_GBS_Spar = read.table("C:/Users/marsit/Desktop/Dropbox/Ploidy/GBS-FA-map-Spar.txt",header=T,sep="\t") 
Info_GBS_Scer = read.table("C:/Users/marsit/Desktop/Dropbox/Ploidy/Info-GBS-Scer.txt",header=T,sep="\t")  

#bind all tables in one
Tab_AF= rbind (table_GBS_Scer, table_GBS_Spar)
Tab_AF=full_join(Tab_AF, Info_GBS_Scer, by= "strain" )

#separate a table for S. cerevisiae reference and S. paradoxus reference
Tab_AF_Scer = filter (Tab_AF, ref=="scer")
Tab_AF_Spar = filter (Tab_AF, ref=="spar")

pdf("Heatmap-GBS-figS13.pdf", width=6, height= 4)  

g_scer <- ggplot(Tab_AF_Scer,aes(x= interaction(POS, CHROM), y= strain_t)) +
  geom_raster(aes (fill=FA)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5)) +
  scale_fill_gradient2(low = "darkred", high = "midnightblue", mid = "white", midpoint = 0.5) + 
  labs(x=NULL) + labs(y=NULL) + theme(axis.ticks=element_blank()) + theme(axis.text.x=element_blank()) + theme(axis.text.y=element_text(size=10))

g_spar <- ggplot(Tab_AF_Spar,aes(x= interaction(POS, CHROM), y= strain_t)) +
  geom_raster(aes (fill=FA)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5)) +
  scale_fill_gradient2(low = "darkred", high = "midnightblue", mid = "white", midpoint = (0.5)) +
  labs(x=NULL) + labs(y=NULL) + theme(axis.ticks=element_blank()) + theme(axis.text.x=element_blank())+ theme(axis.text.y=element_text(size=10))+ theme(legend.position= "none")

g_scer_a <-   g_scer + theme(legend.position= "none")

cowplot::plot_grid(g_scer_a, g_spar, nrow=2, rel_heights = c(300, 300)) %>% 
  cowplot::plot_grid(cowplot::get_legend(g_scer), rel_widths= c(1, 0.2))

dev.off()

#############################################################################################################
##################################################################
#Allele frequency heatmap for M1 tetraploid hybrids-figS14       #
##################################################################
library(dplyr)
library(reshape2)
library(ggplot2)

setwd("/Users/souhirmarsit/Dropbox/Ploidy/GBSdata_4n")

d = read.table("merge_B_10DP_sep.txt", header=T, sep="\t") %>% mutate_all(as.character)

str(d)
names (d)

#####filter to keep all heterozygous loci

d2 = mutate(d, "filtre_parent"=
              ifelse (B82_P0.GT == "./." | B82_P0.GT == "." | B82_P0.GT == "1/1" | B82_P0.GT == "0/0" | B82_P0.GT == "2/2", "fail", 
                      ifelse (B82_P35.GT == "./." | B82_P35.GT == "." | B82_P35.GT == "1/1" | B82_P35.GT == "0/0" | B82_P35.GT == "2/2", "fail", 
                              ifelse (B40_P0.GT == "./." | B40_P0.GT == "." | B40_P0.GT == "1/1" | B40_P0.GT == "0/0" | B40_P0.GT == "2/2", "fail", 
                                      ifelse (B40_P16.GT == "./." | B40_P16.GT == "." | B40_P16.GT == "1/1" | B40_P16.GT == "0/0" | B40_P16.GT == "2/2", "fail",
                                              ifelse (B40_P35.GT == "./." | B40_P35.GT == "." | B40_P35.GT == "1/1" | B40_P35.GT == "0/0" | B40_P35.GT == "2/2", "fail", 
                                                      ifelse (B49_P0.GT == "./." | B49_P0.GT == "." | B49_P0.GT == "1/1" | B49_P0.GT == "0/0" | B49_P0.GT == "2/2", "fail",
                                                              ifelse (B49_P35.GT == "./." | B49_P35.GT == "." | B49_P35.GT == "1/1" | B49_P35.GT == "0/0" | B49_P35.GT == "2/2", "fail",  
                                                                      ifelse (B23_P0.GT == "./." |B23_P0.GT == "." | B23_P0.GT == "1/1" | B23_P0.GT == "0/0" | B23_P0.GT == "2/2", "fail",  
                                                                              ifelse (B23_P35.GT == "./." | B23_P35.GT == "." | B23_P35.GT == "1/1" | B23_P35.GT == "0/0" | B23_P35.GT == "2/2", "fail", "ok")))))))))) 
head(d2)
str(d2)

# change character to factor
d2$filtre_parent=as.factor(d2$filtre_parent)
summary(d2)

d3= filter (d2, filtre_parent == "ok")
head(d3)
summary(d3)

dtestlong = melt (d3, id.vars = c("CHROM", "POS", "REF", "ALT"))

dtestvar= mutate(dtestlong, variable = gsub("AD.", "AD_", variable), strain = gsub("(^.*)[.].*$", "\\1", variable), var = gsub (".*[.](.*)", "\\1", variable))
head(dtestvar)

dtest=dcast(dtestvar, CHROM + POS + REF + ALT+ strain ~ var)

dtestallele= mutate(dtest, GT.AD_1 = gsub ("(.*)[/].*", "\\1", GT), GT.AD_2 = gsub (".*[/](.*)", "\\1", GT), GT.AD_0 = "A" )

AD_0 = dtestallele %>% select(CHROM, POS, strain, DP, AD, GT.AD_0) %>% melt(id.vars = c("CHROM", "POS", "strain", "DP", "GT.AD_0")) %>% mutate(allele = GT.AD_0, AllelicDepth = value) %>% select(-GT.AD_0, -value)
head(AD_0)

AD_1 = dtestallele %>% select(CHROM, POS, strain, DP, AD_1, GT.AD_1) %>% melt(id.vars = c("CHROM", "POS", "strain", "DP", "GT.AD_1")) %>% mutate(allele = GT.AD_1, AllelicDepth = value) %>% select(-GT.AD_1, -value)
head(AD_1)

AD_2 = dtestallele %>% select(CHROM, POS, strain, DP, AD_2, GT.AD_2) %>% melt(id.vars = c("CHROM", "POS", "strain", "DP", "GT.AD_2")) %>% mutate(allele = GT.AD_2, AllelicDepth = value) %>% select(-GT.AD_2, -value) %>% mutate(AllelicDepth = ifelse(AllelicDepth < 0, "NA", AllelicDepth), allele = ifelse( AllelicDepth == "NA", "NA", allele))
head(AD_2)

MONSUPERTABLEAU = rbind(AD_0, AD_1, AD_2)
head(MONSUPERTABLEAU)

MONSUPERTABLEAU$DP=as.integer(MONSUPERTABLEAU$DP)
MONSUPERTABLEAU$AllelicDepth=as.integer(MONSUPERTABLEAU$AllelicDepth)

#filter for read depth more than 20 reads
MONSUPERTABLEAU= filter(MONSUPERTABLEAU, DP >20)

AD1= filter(MONSUPERTABLEAU, variable == "AD_1")
tFA = mutate(AD1, "FA"= (AllelicDepth/DP))
head (tFA)
tFA = filter (tFA, strain != "B82_P0" )
tFA = filter (tFA, strain != "B82_P35")
tFA = filter (tFA, strain != "YPS644")
tFA = filter (tFA, strain != "MSH604" )

#Add chromosome names
setwd("/Users/souhirmarsit/Desktop/MA_Sequences/Chr3_cov_LOH")
chr_tab = read.table("CHr_name.txt", header=T) 
tFA = full_join (tFA, chr_tab, by=c("CHROM"))
tFA$chr_name_f = factor(tFA$chr_name, levels=c("I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI"))
#Add line names and cross identifications
setwd("/Users/souhirmarsit/Dropbox/Ploidy/GBSdata_4n")
info_tab = read.table("Info_GBS_map_B.txt", header=T) 
tFA = full_join (tFA, info_tab, by=c("strain"))

tFA$cross_f = factor(tFA$cross, levels=c('M1_23_Tend', 'M1_23_Tini', 'M1_40_Tend', 'M1_40_Tmid', 'M1_40_Tini', 'M1_49_Tend', 'M1_49_Tini'))
# heat map 

plot_A <- 
  pdf("Heatmap_GBS_4n_figS14.pdf", width=7, height= 4)
ggplot(tFA, aes(x= interaction(POS, chr_name_f), y= cross_f)) +
  geom_raster(aes (fill=FA)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5)) +
  scale_fill_gradient2(low = "darkred", high = "midnightblue", mid = "white", midpoint = 0.5)
plot_A 

dev.off()

#############################################################################################################
##################################################################
#Whole genome coverage of tetraploids -figS18                    #
##################################################################

library(ggplot2)
library(dplyr)
library(reshape2)

setwd("/Users/souhirmarsit/Desktop/MA_Sequences")
#read table with strain informations (name, cross etc)
Table_info="/Users/souhirmarsit/Desktop/MA_Sequences/Table_info_A.txt" 
Table_info <- read.table(Table_info, header = T, sep = "\t", fill = TRUE)

#read table of read depth for each cross (here exemple of L cross)
ref_cov_file="/Users/souhirmarsit/Desktop/MA_Sequences/SpB_cov_strains_A"  

ref_cov_wide <- read.table(ref_cov_file, header = T, sep = "\t", fill = TRUE)
ref_cov_long <- melt(ref_cov_wide, id.vars = c("contig", "pos"))

ref_cov_long = filter (ref_cov_long, variable == "A31_P1"| variable == "A31_P35"| variable == "A51_P1"| variable == "A51_P35"| variable == "A87_P1"| variable == "A87_P35")

#calculate coverage in a window of 10kb
depth_cov_10kb <- ref_cov_long_chr3 %>% group_by(contig) %>% mutate(bin = pos %/% 10000) %>% group_by(contig, variable, bin) %>% summarise(cov = sum(value) / 10000) %>% mutate(Contig = gsub("_pilon", "", contig)) %>% as.data.frame()

ref = strsplit(ref_cov_file, "/", fixed = T)[[1]] %>% tail(1) %>% strsplit(., "_", fixed = T) %>% unlist() %>% head(1)

if (ref != "Scer") {
  ref_conversion_file = paste0(ref, "_tig_rearrangement.txt")
  ref_conversion <- read.table(ref_conversion_file, header = F, col.names = c("Chr", "contig", "Orientation"))  
  depth_cov_10kb_conversion <- left_join(depth_cov_10kb, ref_conversion, by = "contig")
  depth_cov_10kb_conversion$Chr_ord <- factor(depth_cov_10kb_conversion$Chr, levels = c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII", "chrIX", "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI"))
  
} else {
  depth_cov_10kb_conversion <- depth_cov_10kb %>% mutate(Chr_ord = gsub("_pilon", "", contig), contig = Chr_ord)
}
#Calculate average coverage of the whole genome
summary_data <- depth_cov_10kb_conversion %>%  group_by(variable) %>% summarise(mean_cov = mean(cov), median_cov = median(cov)) %>% mutate(ref = ref, category = ifelse(variable == "MSH.604", "SpB_parent", ifelse(variable == "LL2011_004", "SpC_parent", "else"))) %>% as.data.frame()
summary_data = full_join(summary_data, Table_info, by = "variable")
summary_data= filter (summary_data,  variable == "A31_P1"| variable == "A31_P35"| variable == "A51_P1"| variable == "A51_P35"| variable == "A87_P1"| variable == "A87_P35")
#generate a table merge of read coverage and mean read coverage of whole genome
tab_essai = full_join(depth_cov_10kb_conversion, summary_data, by = "variable")
tab_essai = full_join(tab_essai, Table_info, by = "variable")
tab_cov_norm_4n= mutate (tab_essai, "pos"= (bin * 10000))

#Register table with read depth in a window of 10000kb For each cross generate a table
write.table (tab_cov_norm_4n, "/Users/souhirmarsit/Desktop/MA_Sequences/Chr3_cov_LOH/tab_cov_norm_4n_A.txt", row.names = F, quote = F, sep= "\t")

###################################script figure S18

library(ggplot2)
library(dplyr)
library(reshape2)
library(lattice)
library(gridExtra)
library(cowplot)
library(plyr)
library(lemon)
library(grid)

setwd("/Users/souhirmarsit/Desktop/MA_Sequences/Chr3_cov_LOH")
#read table with chromosome names
chr_tab = read.table("CHr_name.txt", header=T) 

#read table with correct names of crosses
I1 = read.table("nom_variable_A.txt", header=T) 
I2 = read.table("nom_variable_D.txt", header=T) 
I3 = read.table("nom_variable_B.txt", header=T) 
I4 = read.table("nom_variable_F.txt", header=T) 

#read table with coverages
d1 = read.table("tab_cov_norm_4n_A.txt", header=T) 
d2 = read.table("tab_cov_norm_4n_D.txt", header=T) 
d3 = read.table("tab_cov_norm_4n_B.txt", header=T) 
d4 = read.table("tab_cov_norm_4n_F.txt", header=T) 

#Join coverage and information tables
d1 = full_join(d1, I1, by="variable")
d2 = full_join(d2, I2, by="variable")
d3 = full_join(d3, I3, by="variable")
d4 = full_join(d4, I4, by="variable")
dall= rbind(d1, d2, d3, d4)
dall=full_join(dall, chr_tab, by="Chr")

#read table with whole genome mean coverages
ds1 = read.table("summary_data_A.txt", header=T) 
ds2 = read.table("summary_data_D.txt", header=T) 
ds3 = read.table("summary_data_B.txt", header=T) 
ds4 = read.table("summary_data_F.txt", header=T) 

#Join whole genome mean coverage 
ds1 = full_join(ds1, I1, by="variable")
ds2 = full_join(ds2, I2, by="variable")
ds3 = full_join(ds3, I3, by="variable")
ds4 = full_join(ds4, I4, by="variable")
dsall= rbind(ds1, ds2, ds3, ds4)


dtot = full_join(dall, dsall, by="variable")
write.table(dtot, file="/Users/souhirmarsit/Dropbox/Ploidy-MA-CRL/tab_fig_data_MA_paper/Tab_data_figS18.txt", quote = FALSE, sep = "\t", row.names = FALSE )

#Separate tables according to the timing

dall_Tini = filter(dall, time.x == "Tini")
dall_Tend = filter(dall, time.x == "Tend")

dsall_Tini = filter(dsall, time == "Tini")
dsall_Tend = filter(dsall, time == "Tend")


###############################################################################################################
blue= rgb(0,0,(139/255))
red= rgb((238/255),0,0)
green= rgb(0,(139/255),0)
black= rgb(0,0,0,)

pdf(file = paste0("/Users/souhirmarsit/Desktop/MA_Sequences/Chr3_cov_LOH/","Cov_4n_tini_all_10kb.pdf"), height = 55, width = 40)

PlotA<-
  dall_Plot_tini <- ggplot(dall_Tini, aes(x = pos, y = cov)) + geom_line() + geom_hline(data = dsall_Tini, aes(yintercept = median_cov), col = "red")+ 
  facet_grid(variable_f~chr_name, scales = "free") + ylim(0, 300)+ 
  theme(axis.title.y = element_text(size=70), axis.text= element_text(size=50), strip.text.y= element_text(face="bold", color= "white", size=73, margin = margin(0,1,0,1, "cm")), strip.text.x= element_text(face="bold", color= "black", size=80, margin = margin(1,0,1,0, "cm")))
+ labs(x=NULL) + labs(y='Read depth') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x=element_blank()) + theme(axis.ticks.x=element_blank())    

g <- ggplot_gtable(ggplot_build(PlotA))
stripr <- which(grepl('strip-r', g$layout$name))
fills <- c(blue, blue, blue, blue, blue, green, green, black)

k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g)
dev.off()

pdf(file = paste0("/Users/souhirmarsit/Desktop/MA_Sequences/Chr3_cov_LOH/","Cov_4n_tend_all_10kb_essai.pdf"), height = 55, width = 40)
PlotB<-
  dall_Plot_tend <- ggplot(dall_Tend, aes(x = pos, y = cov)) + geom_line() + geom_hline(data = dsall_Tend, aes(yintercept = median_cov), col = "red")+ 
  facet_grid(variable_f~chr_name, scales = "free") + ylim(0, 300)+ 
  theme(axis.title.y = element_text(size=70), axis.text= element_text(size=50), strip.text.y= element_text(face="bold", color= "white", size=73, margin = margin(0,1,0,1, "cm")), strip.text.x= element_text(face="bold", color= "black", size=80, margin = margin(1,0,1,0, "cm")))
+ labs(x=NULL) + labs(y='Read depth') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x=element_blank()) + theme(axis.ticks.x=element_blank())    

H <- ggplot_gtable(ggplot_build(PlotB))
stripr <- which(grepl('strip-r', H$layout$name))
fills <- c(blue, blue, blue, blue, blue, green, green, black)

l <- 1
for (z in stripr) {
  w <- which(grepl('rect', H$grobs[[z]]$grobs[[1]]$childrenOrder))
  H$grobs[[z]]$grobs[[1]]$children[[w]]$gp$fill <- fills[l]
  l <- l+1
}
grid.draw(H)
dev.off()

##########################################
#############################################################################################################
##################################################################
#Chr III coverage of tetraploids -figS19                         #
##################################################################
library(ggplot2)
library(dplyr)
library(reshape2)

setwd("/Users/souhirmarsit/Desktop/MA_Sequences")
#read table with strain informations (name, cross etc)
Table_info="/Users/souhirmarsit/Desktop/MA_Sequences/Table_info_A.txt" 
Table_info <- read.table(Table_info, header = T, sep = "\t", fill = TRUE)

#read table of read depth for each cross (here exemple of L cross)
ref_cov_file="/Users/souhirmarsit/Desktop/MA_Sequences/SpB_cov_strains_A"  

ref_cov_wide <- read.table(ref_cov_file, header = T, sep = "\t", fill = TRUE)
ref_cov_long <- melt(ref_cov_wide, id.vars = c("contig", "pos"))

#filter for read depth of chr III
ref_cov_long_chr3 = filter (ref_cov_long, contig == "utg584_pilon")
ref_cov_long_chr3 = filter (ref_cov_long_chr3, variable == "A31_P1"| variable == "A31_P35"| variable == "A51_P1"| variable == "A51_P35"| variable == "A87_P1"| variable == "A87_P35")
#head (ref_cov_wide)
#calculate coverage in a window of 100bp
depth_cov_10kb <- ref_cov_long_chr3 %>% group_by(contig) %>% mutate(bin = pos %/% 100) %>% group_by(contig, variable, bin) %>% summarise(cov = sum(value) / 100) %>% mutate(Contig = gsub("_pilon", "", contig)) %>% as.data.frame()

ref = strsplit(ref_cov_file, "/", fixed = T)[[1]] %>% tail(1) %>% strsplit(., "_", fixed = T) %>% unlist() %>% head(1)

if (ref != "Scer") {
  ref_conversion_file = paste0(ref, "_tig_rearrangement.txt")
  
  ref_conversion <- read.table(ref_conversion_file, header = F, col.names = c("Chr", "contig", "Orientation"))
  
  depth_cov_100bp_conversion <- left_join(depth_cov_10kb, ref_conversion, by = "contig")
  depth_cov_100bp_conversion$Chr_ord <- factor(depth_cov_10kb_conversion$Chr, levels = c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII", "chrIX", "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI"))
  
} else {
  depth_cov_100bp_conversion <- depth_cov_10kb %>% mutate(Chr_ord = gsub("_pilon", "", contig), contig = Chr_ord)
}
#Calculate average coverage of the whole chromosome
summary_data <- depth_cov_100bp_conversion %>%  group_by(variable) %>% summarise(mean_cov = mean(cov), median_cov = median(cov)) %>% mutate(ref = ref, category = ifelse(variable == "MSH.604", "SpB_parent", ifelse(variable == "LL2011_004", "SpC_parent", "else"))) %>% as.data.frame()
summary_data = full_join(summary_data, Table_info, by = "variable")
summary_data= filter (summary_data,  variable == "A31_P1"| variable == "A31_P35"| variable == "A51_P1"| variable == "A51_P35"| variable == "A87_P1"| variable == "A87_P35")
 #full join of coverage and average coverage of the whole chromosome
tab_essai = full_join(depth_cov_100bp_conversion, summary_data, by = "variable")
tab_essai = full_join(tab_essai, Table_info, by = "variable")
tab_cov_norm_4n_chr3= mutate (tab_essai, "pos"= (bin * 100))

#Register table with read depth of chr III in a window of 100bp. For each cross generate a table
write.table (tab_cov_norm_4n_chr3, "/Users/souhirmarsit/Desktop/MA_Sequences/Chr3_cov_LOH/tab_cov_norm_4n_chr3_A.txt", row.names = F, quote = F, sep= "\t")

##################################################################
#######################################script figS19#

library(ggplot2)
library(dplyr)
library(reshape2)
library(lattice)
library(gridExtra)
library(cowplot)
library(plyr)
library(lemon)

setwd("/Users/souhirmarsit/Desktop/MA_Sequences/Chr3_cov_LOH")

#Read tables of straines infomations (name, cross, time etc)
I1 = read.table("nom_variable_A.txt", header=T) 
I2 = read.table("nom_variable_D.txt", header=T) 
I3 = read.table("nom_variable_B.txt", header=T) 
I4 = read.table("nom_variable_F.txt", header=T) 

#Read tables of chr III read depth
d1 = read.table("tab_cov_norm_4n_Chr3_A.txt", header=T) 
d2 = read.table("tab_cov_norm_4n_Chr3_D.txt", header=T) 
d3 = read.table("tab_cov_norm_4n_Chr3_B.txt", header=T) 
d4 = read.table("tab_cov_norm_4n_Chr3_F.txt", header=T) 

#full join tables of chr III read depth and strains informations
d1 = full_join(d1, I1, by="variable")
d2 = full_join(d2, I2, by="variable")
d3 = full_join(d3, I3, by="variable")
d4 = full_join(d4, I4, by="variable")
dall= rbind(d1, d2, d3, d4)

#Read tables of average whole chr III read depth
ds1 = read.table("summary_data_A.txt", header=T) 
ds2 = read.table("summary_data_D.txt", header=T) 
ds3 = read.table("summary_data_B.txt", header=T) 
ds4 = read.table("summary_data_F.txt", header=T) 

ds1 = full_join(ds1, I1, by="variable")
ds2 = full_join(ds2, I2, by="variable")
ds3 = full_join(ds3, I3, by="variable")
ds4 = full_join(ds4, I4, by="variable")
dsall= rbind(ds1, ds2, ds3, ds4)

#Separate tables according to the timing

dall_Tini = filter(dall, time.x == "Tini")
dall_Tend = filter(dall, time.x == "Tend")

dsall_Tini = filter(dsall, time == "Tini")
dsall_Tend = filter(dsall, time == "Tend")


blue= rgb(0,0,(139/255))
red= rgb((238/255),0,0)
green= rgb(0,(139/255),0)
black= rgb(0,0,0,)

pdf(file = paste0("/Users/souhirmarsit/Desktop/MA_Sequences/Chr3_cov_LOH/","Cov_4n_tini_all_10kb.pdf"), height = 55, width = 40)

PlotA<-
  dall_Plot_tini <- ggplot(dall_Tini, aes(x = pos, y = cov)) + geom_line() + geom_hline(data = dsall_Tini, aes(yintercept = median_cov), col = "red")+ 
  facet_grid(variable_f~chr_name, scales = "free") + ylim(0, 300)+ 
  theme(axis.title.y = element_text(size=70), axis.text= element_text(size=50), strip.text.y= element_text(face="bold", color= "white", size=73, margin = margin(0,1,0,1, "cm")), strip.text.x= element_text(face="bold", color= "black", size=80, margin = margin(1,0,1,0, "cm")))
+ labs(x=NULL) + labs(y='Read depth') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x=element_blank()) + theme(axis.ticks.x=element_blank())    

g <- ggplot_gtable(ggplot_build(PlotA))
stripr <- which(grepl('strip-r', g$layout$name))
fills <- c(blue, blue, blue, blue, blue, green, green, black)

k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g)
dev.off()

pdf(file = paste0("/Users/souhirmarsit/Desktop/MA_Sequences/Chr3_cov_LOH/","Cov_4n_tend_all_10kb_essai.pdf"), height = 55, width = 40)
PlotB<-
  dall_Plot_tend <- ggplot(dall_Tend, aes(x = pos, y = cov)) + geom_line() + geom_hline(data = dsall_Tend, aes(yintercept = median_cov), col = "red")+ 
  facet_grid(variable_f~chr_name, scales = "free") + ylim(0, 300)+ 
  theme(axis.title.y = element_text(size=70), axis.text= element_text(size=50), strip.text.y= element_text(face="bold", color= "white", size=73, margin = margin(0,1,0,1, "cm")), strip.text.x= element_text(face="bold", color= "black", size=80, margin = margin(1,0,1,0, "cm")))
+ labs(x=NULL) + labs(y='Read depth') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x=element_blank()) + theme(axis.ticks.x=element_blank())    

H <- ggplot_gtable(ggplot_build(PlotB))
stripr <- which(grepl('strip-r', H$layout$name))
fills <- c(blue, blue, blue, blue, blue, green, green, black)

l <- 1
for (z in stripr) {
  w <- which(grepl('rect', H$grobs[[z]]$grobs[[1]]$childrenOrder))
  H$grobs[[z]]$grobs[[1]]$children[[w]]$gp$fill <- fills[l]
  l <- l+1
}
grid.draw(H)
dev.off()

#############################################################################################################
#########################################################################
#Whole genome allele frequency of chromosome III of tetraploids -figS20 #
#########################################################################

###########Analyse VCF merge and generate tables with allele frequency of all variants or corresponding to hybrid parent 2 for each cross (here exemple for L cross)
library(dplyr)
library(reshape2)
library(ggplot2)
library(plyr)

setwd("/Users/souhirmarsit/Desktop/MA_Sequences")

#Read VCF merged table for each cross (here exemle for L cross)
d = read.table("merge_SpB_A_vcf_sep.txt", header=T, sep='\t') %>% mutate_all(as.character)
d = select(d, CHROM, POS, REF, ALT, contains("A31"), contains("A51"), contains("A87"), contains("LL2011_004"), contains("MSH.604"))

#filter for allele frequency of chromosome III
d = filter (d, CHROM == "utg584_pilon")
str(d)
names (d)

#####filtrer for heterozygous variants between the two hybrid parents

d2 = mutate(d, "filtre_parent"=
              ifelse (LL2011_004.GT == MSH.604.GT, "fail",
                      ifelse (LL2011_004.GT == "./." | 	MSH.604.GT == "./.", "fail", "ok"))) 
head(d2)
str(d2)

d2$filtre_parent=as.factor(d2$filtre_parent)
summary(d2)

#filter d2 to keep only heterozygous variants

d3= filter (d2, filtre_parent == "ok")
head(d3)
summary(d3)

dtest= select (d3, CHROM , POS, REF, ALT, starts_with ("LL2011_004"), starts_with ("MSH.604"), starts_with ("A"))

#Identify parental alleles

dtestlong = melt (dtest, id.vars = c("CHROM", "POS", "REF", "ALT"))

dtestvar= mutate(dtestlong, variable = gsub("AD.", "AD_", variable), strain = gsub("(^.*)[.].*$", "\\1", variable), var = gsub (".*[.](.*)", "\\1", variable))
head(dtestvar)

dtest=dcast(dtestvar, CHROM + POS + REF + ALT+ strain ~ var)

dtestallele= mutate(dtest, GT.AD_1 = gsub ("(.*)[/].*", "\\1", GT), GT.AD_2 = gsub (".*[/](.*)", "\\1", GT), GT.AD_0 = "B")
d_REF = dtestallele %>% filter(strain == "LL2011_004") %>%  select(-REF, -ALT, -starts_with("AD"), -GT, -GT.AD_2, -GT.AD_0) %>% dcast(CHROM + POS ~ strain)
head(d_REF)

dalBC= full_join(d_REF, dtestallele, by = c("CHROM", "POS"))

dalleleBC = mutate(dalBC, GT.AD_1= "C", GT.AD_2 = "MUT")

head(dalleleBC)

AD_0 = dalleleBC %>% select(CHROM, POS, strain, DP, AD, GT.AD_0) %>% melt(id.vars = c("CHROM", "POS", "strain", "DP", "GT.AD_0")) %>% mutate(allele = GT.AD_0, AllelicDepth = value) %>% select(-GT.AD_0, -value)
head(AD_0)

AD_1 = dalleleBC %>% select(CHROM, POS, strain, DP, AD_1, GT.AD_1) %>% melt(id.vars = c("CHROM", "POS", "strain", "DP", "GT.AD_1")) %>% mutate(allele = GT.AD_1, AllelicDepth = value) %>% select(-GT.AD_1, -value)
head(AD_1)

AD_2 = dalleleBC %>% select(CHROM, POS, strain, DP, AD_2, GT.AD_2) %>% melt(id.vars = c("CHROM", "POS", "strain", "DP", "GT.AD_2")) %>% mutate(allele = GT.AD_2, AllelicDepth = value) %>% select(-GT.AD_2, -value) %>% mutate(AllelicDepth = ifelse(AllelicDepth < 0, "NA", AllelicDepth), allele = ifelse( AllelicDepth == "NA", "NA", allele))
head(AD_2)

MONSUPERTABLEAU = rbind(AD_0, AD_1, AD_2)

#Identification of contigs -> chromosome
ref_conversion <- read.table("SpB_tig_rearrangement.txt", header = F, col.names = c("Chr", "CHROM", "Orientation"))
head (ref_conversion)
SUPERTABLEAU <- left_join(MONSUPERTABLEAU, ref_conversion, by = "CHROM")
head (SUPERTABLEAU)
SUPERTABLEAU$Chr_ord <- factor(SUPERTABLEAU$Chr, levels = c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII", "chrIX", "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI"))

#Identification of strains name
Info_strain <- read.table("Table_info_A_LOH.txt", header = T, sep = "\t", fill = TRUE)
head (Info_strain)

SUPERTABLEAU$DP=as.integer(MONSUPERTABLEAU$DP)
SUPERTABLEAU$AllelicDepth=as.integer(MONSUPERTABLEAU$AllelicDepth)

#filter for variants with more than 20 reads 
tabDP10 = filter (SUPERTABLEAU, DP>20)

#Calculate allele frequency the ratio of allelic depth and the read depth of the locus   
tFA = mutate(tabDP10, "FA"= (AllelicDepth/DP))
head (tFA)

#filter for the allele frequency of Hybrid parent 2 (SpC or SpA or S. cerevisiae) of chr III
tFA_f = filter (tFA, allele == "C")
tFA_f <- left_join(tFA_f, Info_strain, by = "strain")
head(tFA_f)
summary(tFA_f)
#Register the table of allele frequency corresponding to the hybrid parent2
write.table (tFA_f, "/Users/souhirmarsit/Desktop/MA_Sequences/LOH_4n/tableau_LOH_A_chr3.txt", row.names = F, quote = F, sep= "\t")


##################Script figure S20

library(dplyr)
library(reshape2)
library(ggplot2)
library(plyr)

setwd("/Users/souhirmarsit/Desktop/MA_Sequences/Chr3_cov_LOH")

#read tables with strains informations (name, cross etc)
I1 = read.table("nom_strain_A.txt", header=T) 
I2 = read.table("nom_strain_D.txt", header=T) 
I3 = read.table("nom_strain_B.txt", header=T) 
I4 = read.table("nom_strain_F.txt", header=T) 

#read tables of allele frequency of chr III for each cross
d1 = read.table("tableau_LOH_A_chr3.txt", header=T) 
d2 = read.table("tableau_LOH_D_chr3.txt", header=T) 
d3 = read.table("tableau_LOH_B_chr3.txt", header=T) 
d4 = read.table("tableau_LOH_F_chr3.txt", header=T) 

d1 = full_join(d1, I1, by="strain")
d2 = full_join(d2, I2, by="strain")
d3 = full_join(d3, I3, by="strain")
d4 = full_join(d4, I4, by="strain")

str(d)
names (d)

d1 = filter (d1, FA>0.01)
d2 = filter (d2, FA>0.01)
d3 = filter (d3, FA>0.01)
d4 = filter (d4, FA>0.01)

dtot = rbind (d1, d2, d3, d4)
write.table(dtot, file="/Users/souhirmarsit/Dropbox/Ploidy-MA-CRL/tab_fig_data_MA_paper/Tab_data_figS17.txt", quote = FALSE, sep = "\t", row.names = FALSE )

#separate tables according to the timing Tini and Tend
d1_Tini = filter(d1, time == "Tini")
d1_Tend = filter(d1, time == "Tend")

d2_Tini = filter(d2, time == "Tini")
d2_Tend = filter(d2, time == "Tend")

d3_Tini = filter(d3, time == "Tini")
d3_Tend = filter(d3, time == "Tend")

d4_Tini = filter(d4, time == "Tini")
d4_Tend = filter(d4, time == "Tend")

#reorder the line names
d1_Tini$variable_ord = factor(d1_Tini$variable_f, levels=c("L1_87_Tini", "L1_51_Tini", "L1_31_Tini"))
d2_Tini$variable_ord = factor(d2_Tini$variable_f, levels=c("L2_45_Tini", "L2_36_Tini"))
d3_Tini$variable_ord = factor(d3_Tini$variable_f, levels=c("M1_49_Tini", "M1_40_Tini"))
d4_Tini$variable_ord = factor(d4_Tini$variable_f, levels=c("H2_38_Tini"))

d1_Tend$variable_ord = factor(d1_Tend$variable_f, levels=c("L1_87_Tend", "L1_51_Tend", "L1_31_Tend"))
d2_Tend$variable_ord = factor(d2_Tend$variable_f, levels=c("L2_45_Tend", "L2_36_Tend"))
d3_Tend$variable_ord = factor(d3_Tend$variable_f, levels=c("M1_49_Tend", "M1_40_Tend"))
d4_Tend$variable_ord = factor(d4_Tend$variable_f, levels=c("H2_38_Tend"))


pdf(file = paste0("/Users/souhirmarsit/Desktop/MA_Sequences/Chr3_cov_LOH/","LOH_4n_Chr3_verif.pdf"), height = 5, width = 15)

A_tini <- ggplot(d1_Tini, aes(x= interaction(POS, Chr), y= variable_ord)) + geom_raster(aes (fill=FA)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5)) + theme(axis.text.y = element_text( size = 12)) +
  scale_fill_gradient2(low = "darkred", high = "midnightblue", mid = "white", midpoint = 0.5, limits=c(0, 1), name= "Allele frequency")+ theme(axis.text.x=element_blank())+ labs(x="ChrIII Position") + labs(y=NULL) + labs(x=NULL)   
D_tini <- ggplot(d2_Tini, aes(x= interaction(POS, Chr), y= variable_ord)) + geom_raster(aes (fill=FA)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5)) + theme(axis.text.y = element_text( size = 12)) +
  scale_fill_gradient2(low = "darkred", high = "midnightblue", mid = "white", midpoint = 0.5, limits=c(0, 1))+ theme(axis.text.x=element_blank())+ theme(legend.position = "none") + labs(x="ChrIII Position") + labs(y=NULL) + labs(x=NULL) 
B_tini <- ggplot(d3_Tini, aes(x= interaction(POS, Chr), y= variable_ord)) + geom_raster(aes (fill=FA)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5)) + theme(axis.text.y = element_text( size = 12)) +
  scale_fill_gradient2(low = "darkred", high = "midnightblue", mid = "white", midpoint = 0.5, limits=c(0, 1))+ theme(axis.text.x=element_blank())+ theme(legend.position = "none") + labs(x="ChrIII Position") + labs(y=NULL) + labs(x=NULL) 
F_tini <- ggplot(d4_Tini, aes(x= interaction(POS, Chr), y= variable_ord)) + geom_raster(aes (fill=FA)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5)) + theme(axis.text.y = element_text( size = 12)) +
  scale_fill_gradient2(low = "darkred", high = "midnightblue", mid = "white", midpoint = 0.5, limits=c(0, 1))+ theme(axis.text.x=element_blank())+ theme(legend.position = "none") + labs(x="ChrIII Position") + labs(y=NULL) + labs(x=NULL) 

A_tend <- ggplot(d1_Tend, aes(x= interaction(POS, Chr), y= variable_ord)) + geom_raster(aes (fill=FA)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5)) + theme(axis.text.y = element_text( size = 12)) +
  scale_fill_gradient2(low = "darkred", high = "midnightblue", mid = "white", midpoint = 0.5, limits=c(0, 1))+ theme(axis.text.x=element_blank())+ theme(legend.position = "none") + labs(x="ChrIII Position") + labs(y=NULL) + labs(x=NULL) 
D_tend <- ggplot(d2_Tend, aes(x= interaction(POS, Chr), y= variable_ord)) + geom_raster(aes (fill=FA)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5)) + theme(axis.text.y = element_text( size = 12)) +
  scale_fill_gradient2(low = "darkred", high = "midnightblue", mid = "white", midpoint = 0.5, limits=c(0, 1))+ theme(axis.text.x=element_blank())+ theme(legend.position = "none") + labs(x="ChrIII Position") + labs(y=NULL) + labs(x=NULL) 
B_tend <- ggplot(d3_Tend, aes(x= interaction(POS, Chr), y= variable_ord)) + geom_raster(aes (fill=FA)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5)) + theme(axis.text.y = element_text( size = 12)) +
  scale_fill_gradient2(low = "darkred", high = "midnightblue", mid = "white", midpoint = 0.5, limits=c(0, 1))+ theme(axis.text.x=element_blank())+ theme(legend.position = "none") + labs(x="ChrIII Position") + labs(y=NULL) + labs(x=NULL) 
F_tend <- ggplot(d4_Tend, aes(x= interaction(POS, Chr), y= variable_ord)) + geom_raster(aes (fill=FA)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5)) + theme(axis.text.y = element_text( size = 12)) +
  scale_fill_gradient2(low = "darkred", high = "midnightblue", mid = "white", midpoint = 0.5, limits=c(0, 1))+ theme(axis.text.x=element_blank())+ theme(legend.position = "none") + labs(x="ChrIII Position") + labs(y=NULL) + labs(x=NULL) 

A_tini_a <- A_tini + theme(legend.position = "none")

cowplot::plot_grid(A_tini_a, A_tend, D_tini, D_tend, B_tini, B_tend, F_tini, F_tend, nrow=4, ncol=2, rel_heights= c(105, 75, 75, 42)) %>%
  cowplot::plot_grid(cowplot::get_legend(A_tini), rel_widths = c(1, 0.1))
dev.off()


#############################################################################################################
##################################################################
#copy number analysis of MAT locus -figS21A                      #
##################################################################

library(ggplot2)
library(dplyr)
library(reshape2)

setwd("/Users/souhirmarsit/desktop/MA_Sequences/Bam_files_Mat_locus/MAT_locus_paradoxus_data")
#read table of strains information (name, cross etc)
Table_info="/Users/souhirmarsit/desktop/MA_Sequences/Bam_files_Mat_locus/MAT_locus_paradoxus_data/info_strain_MAT_tot.txt" 
Table_info <- read.table(Table_info, header = T, sep = "\t", fill = TRUE)
#read table with coverage data 
ref_cov_file="/Users/souhirmarsit/Desktop/MA_Sequences/Bam_files_Mat_locus/MAT_locus_paradoxus_data/cov_MAT_locus_strains"  

ref_cov_wide <- read.table(ref_cov_file, header = T, sep = "\t", fill = TRUE)
ref_cov_long1 <- melt(ref_cov_wide, id.vars = c("contig", "pos"))
ref_cov_long = full_join(ref_cov_long1, Table_info, by = "variable")
ref_cov_long=filter (ref_cov_long, contig=="MAT"| contig=="V01313.1")

#keep positions that are different between MAT-a and MAT-alpha
ref_cov_long=filter (ref_cov_long, pos<3770)
ref_cov_long=filter (ref_cov_long, pos>2435)

#Mean of coverage for the MAT locus a and Alpha
summary_data <- ref_cov_long %>% group_by(variable) %>% summarise(mean_cov = mean(value), median_cov = median(value))  %>% as.data.frame()

summary_data_chr <- ref_cov_long %>% group_by(contig, variable) %>% summarise(mean_cov = mean(value), median_cov = median(value)) %>% as.data.frame()

summary_data_chr_MAT_A = filter (summary_data_chr, contig=="V01313.1")

summary_data_chr_MAT_alpha = filter (summary_data_chr, contig=="MAT")

tab_MAT_A_alpha = full_join(summary_data_chr_MAT_A, summary_data_chr_MAT_alpha, by = "variable")
#coverage of MAT-a /MAT-alpha
tab_MAT_CN = mutate(tab_MAT_A_alpha, "CN"=(median_cov.x/median_cov.y))

tab_MAT_CN = full_join(tab_MAT_CN, Table_info, by = "variable")

tab_MAT_CN1 <- tab_MAT_CN %>% group_by(strain_f) %>% summarise(mean_CN = mean(CN)) %>% as.data.frame()
tab_MAT_CN1 = full_join(tab_MAT_CN1, Table_info, by = "strain_f")
tab_MAT_CN1= select(tab_MAT_CN1, -R, -variable)
tab_MAT_CN1=unique(tab_MAT_CN1)

write.table(tab_MAT_CN1, file="/Users/souhirmarsit/Dropbox/Ploidy-MA-CRL/tab_fig_data_MA_paper/Tab_data_figS21A.txt", quote = FALSE, sep = "\t", row.names = FALSE )

###################################script figure S21A

library(lattice)
library(gridExtra)
library(cowplot)
#library(plyr)
library(lemon)

tab_MAT_CN1$strain_ord = factor(tab_MAT_CN1$cross, levels=c("L1_31", "L1_51", "L1_87", "L2_36", "L2_45", "M1_40", "M1_49", "H2_38", "SpC_parent1", "SpC_parent2", "SpA_parent1", "Scer_parent2", "SpB_parent1", "SpB_parent2"))
tab_MAT_CN1$time_ord = factor(tab_MAT_CN1$time, levels=c("Tini", "Tend", "Parent"))

pdf(file = paste0("/Users/souhirmarsit/Desktop/MA_Sequences/Bam_files_Mat_locus/MAT_locus_paradoxus_data/","MAT_locus_paradoxus_cov.pdf"), height = 8, width = 12)
ggplot(tab_MAT_CN1, aes(x = strain_ord, y = mean_CN, group=time_ord, fill = factor(time_ord)))  + 
  geom_bar(stat="identity", position="dodge", colour="black")+ theme(legend.title = element_blank())+ theme(legend.text = element_text(size=18))+ scale_fill_manual(values=c("lightcyan2","lightblue4", "grey"))+
  ylab("cov(MAT_a/MAT_alpha)")+ theme(axis.ticks=element_blank()) + theme(axis.text.x=element_text(size=18))+ labs(x=NULL) +
  theme(axis.text.y=element_text(size=18), axis.title.y=element_text(size=18)) + ylim(0, 2.1)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

#############################################################################################################
##################################################################
#Allele frquency analysis of MAT locus -figS21B                  #
##################################################################
###########Analyses of VCF merge for each cross (here exemple of L cross)
library(dplyr)
library(reshape2)
library(ggplot2)
library(plyr)

setwd("/Users/souhirmarsit/Desktop/MA_Sequences/Bam_files_Mat_locus/MAT_locus_paradoxus_data")

#read vcf file for each cross
d = read.table("merge_MAT_locus_A_sep.txt", header=T) %>% mutate_all(as.character)
str(d)
names (d)

#filter for heterozygous variants (alleles heterozygous between the two hybrid parents)
d2 = mutate(d, "filtre_parent"=
              ifelse (LL2011_004.GT == MSH.604.GT, "fail",
                      ifelse (LL2011_004.GT == "./." | 	MSH.604.GT == "./.", "fail", "ok"))) 
head(d2)
str(d2)

d2$filtre_parent=as.factor(d2$filtre_parent)
summary(d2)

#filter for only heterozygous variants

d3= filter (d2, filtre_parent == "ok")
head(d3)
summary(d3)

#Identify parental alleles
dtest= select (d3, CHROM , POS, REF, ALT, starts_with ("LL2011_004"), starts_with ("MSH.604"), starts_with ("A"))

dtestlong = melt (dtest, id.vars = c("CHROM", "POS", "REF", "ALT"))

dtestvar= mutate(dtestlong, variable = gsub("AD.", "AD_", variable), strain = gsub("(^.*)[.].*$", "\\1", variable), var = gsub (".*[.](.*)", "\\1", variable))
head(dtestvar)

dtest=dcast(dtestvar, CHROM + POS + REF + ALT+ strain ~ var)

dtestallele= mutate(dtest, GT.AD_1 = gsub ("(.*)[/].*", "\\1", GT), GT.AD_2 = gsub (".*[/](.*)", "\\1", GT), GT.AD_0 = "A")
d_REF_C = dtestallele %>% filter(strain == "LL2011_004") %>%  select(-REF, -ALT, -starts_with("AD"), -GT, -GT.AD_2, -GT.AD_0) %>% dcast(CHROM + POS ~ strain)
d_REF_C = filter(d_REF_C, LL2011_004 == 1 | LL2011_004 == 2)

d_REF_B = dtestallele %>% filter(strain == "MSH.604") %>%  select(-REF, -ALT, -starts_with("AD"), -GT, -GT.AD_2, -GT.AD_0) %>% dcast(CHROM + POS ~ strain)
d_REF_B = filter(d_REF_B, MSH.604 == 1 | MSH.604 == 2)

dalC= full_join(d_REF_C, dtestallele, by = c("CHROM", "POS"))
dalC= filter(dalC, LL2011_004 == 1 | LL2011_004 == 2)
dalB= full_join(d_REF_B, dtestallele, by = c("CHROM", "POS"))
dalB = filter(dalB, MSH.604 == 1 | MSH.604 == 2)

dalleleC = mutate(dalC, GT.AD_1= "C", GT.AD_2 = "MUT")%>%  select(-LL2011_004)
dalleleB = mutate(dalB, GT.AD_1= "B", GT.AD_2 = "MUT")%>%  select(-MSH.604)
dalleleBC= rbind (dalleleC, dalleleB)
head(dalleleBC)

AD_0 = dalleleBC %>% select(CHROM, POS, strain, DP, AD, GT.AD_0) %>% melt(id.vars = c("CHROM", "POS", "strain", "DP", "GT.AD_0")) %>% mutate(allele = GT.AD_0, AllelicDepth = value) %>% select(-GT.AD_0, -value)
head(AD_0)

AD_1 = dalleleBC %>% select(CHROM, POS, strain, DP, AD_1, GT.AD_1) %>% melt(id.vars = c("CHROM", "POS", "strain", "DP", "GT.AD_1")) %>% mutate(allele = GT.AD_1, AllelicDepth = value) %>% select(-GT.AD_1, -value)
head(AD_1)

AD_2 = dalleleBC %>% select(CHROM, POS, strain, DP, AD_2, GT.AD_2) %>% melt(id.vars = c("CHROM", "POS", "strain", "DP", "GT.AD_2")) %>% mutate(allele = GT.AD_2, AllelicDepth = value) %>% select(-GT.AD_2, -value) %>% mutate(AllelicDepth = ifelse(AllelicDepth < 0, "NA", AllelicDepth), allele = ifelse( AllelicDepth == "NA", "NA", allele))
head(AD_2)

MONSUPERTABLEAU = rbind(AD_0, AD_1, AD_2)

head(MONSUPERTABLEAU)
MONSUPERTABLEAU$DP=as.integer(MONSUPERTABLEAU$DP)
MONSUPERTABLEAU$AllelicDepth=as.integer(MONSUPERTABLEAU$AllelicDepth)

#filter for variants with more than 20 reads
tabDP20 = filter (MONSUPERTABLEAU, DP>20)
#calculate Allele frequency, the ratio of allelic depth and read depth of the locus
tFA = mutate(tabDP20, "FA"= (AllelicDepth/DP))
head (tFA)
#register the table of allele frequency for each cross
write.table (tFA, "/Users/souhirmarsit/Desktop/MA_Sequences/tableau_MAT_FA_A.txt", row.names = F, quote = F, sep= "\t")

##################################script figure S21B

library(dplyr)
library(reshape2)
library(ggplot2)

setwd("/Users/souhirmarsit/Desktop/MA_Sequences/Bam_files_Mat_locus/MAT_locus_paradoxus_data")

#read tables with allele frequency for each cross
dA = read.table("tableau_MAT_FA_A.txt", header=T)
dD = read.table("tableau_MAT_FA_D.txt", header=T)
dB = read.table("tableau_MAT_FA_B.txt", header=T)
dF = read.table("tableau_MAT_FA_F.txt", header=T) 
dall=rbind(dA, dD, dB, dF)

Table_info="/Users/souhirmarsit/desktop/MA_Sequences/Bam_files_Mat_locus/MAT_locus_paradoxus_data/info_strain_MAT_tot.txt" 
Table_info <- read.table(Table_info, header = T, sep = "\t", fill = TRUE)
Table_Mat="/Users/souhirmarsit/desktop/MA_Sequences/Bam_files_Mat_locus/MAT_locus_paradoxus_data/Matname.txt" 
Table_Mat <- read.table(Table_Mat, header = T, sep = "\t", fill = TRUE)

tFA_A= full_join(dA, Table_info, by="variable")
tFA_D= full_join(dD, Table_info, by="variable")
tFA_B= full_join(dB, Table_info, by="variable")
tFA_F= full_join(dF, Table_info, by="variable")

tFA_all= full_join(dall, Table_Mat, by="CHROM")
tFA_all= full_join(tFA_all, Table_info)
tFA_alll= select (tFA_all, CHROM, strain_f,allele , FA) 

#calculate the mean allele frequency for the two mating type sequences
tFA_all_mean <- tFA_alll %>% group_by(strain_f, CHROM, allele) %>% summarise(mean_FA = mean(FA))  %>% as.data.frame()
tFA_all_mean= full_join(tFA_all_mean, Table_info, by="strain_f") %>% filter (mean_FA > 0.01)
tFA_all_mean= full_join(tFA_all_mean, Table_Mat, by="CHROM")
tFA_all_mean= select(tFA_all_mean, -R, -variable)
tFA_all_mean=unique(tFA_all_mean)
tFA_all_mean1= filter  (tFA_all_mean, allele=="C"|  allele=="B")

write.table(tFA_all_mean1, file="/Users/souhirmarsit/Dropbox/Ploidy-MA-CRL/tab_fig_data_MA_paper/Tab_data_figS21B.txt", quote = FALSE, sep = "\t", row.names = FALSE )

#############################################
#############################################
#Separate data according to the time and the parental alleles
tFA_mean_tini= filter (tFA_all_mean, time=="Tini")
tFA_mean_tend= filter (tFA_all_mean, time=="Tend")
tFA_mean_tini_C= filter (tFA_mean_tini, allele=="C")
tFA_mean_tend_C= filter (tFA_mean_tend, allele=="C")
tFA_mean_tini_B= filter (tFA_mean_tini, allele=="B")
tFA_mean_tend_B= filter (tFA_mean_tend, allele=="B")

#############################################

tFA_mean_C= filter (tFA_all_mean, allele=="C")%>% filter (mean_FA < 0.99)
tFA_mean_B= filter (tFA_all_mean, allele=="B")%>% filter (mean_FA < 0.99)

tFA_mean_C$time_f = factor(tFA_mean_C$time, levels=c("Tini", "Tend"))
tFA_mean_B$time_f = factor(tFA_mean_B$time, levels=c("Tini", "Tend"))

tFA_mean_C$cross_ord = factor(tFA_mean_C$cross, levels=c("H2_38", "M1_49", "M1_40","L2_45", "L2_36", "L1_87", "L1_51", "L1_31"))
tFA_mean_B$cross_ord = factor(tFA_mean_B$cross, levels=c("H2_38", "M1_49", "M1_40","L2_45", "L2_36", "L1_87", "L1_51", "L1_31"))

pdf(file = paste0("/Users/souhirmarsit/Desktop/MA_Sequences/Bam_files_Mat_locus/MAT_locus_paradoxus_data/","AF_MAT_locus_4n.pdf"), height = 7, width = 10)
g_C <- ggplot(tFA_mean_C, aes(x= locusname, y= cross_ord)) +
  ggtitle("Allele parent2") +geom_raster(aes (fill=mean_FA)) + theme(plot.title = element_text(color="grey40", size=18, face="bold"), axis.text.x = element_text(angle = 90, hjust = 1, size = 18), axis.title.x = element_blank(), strip.text = element_text(colour = 'black', face="bold", size = 18)) +
  theme(axis.text.y = element_text(hjust = 1, size = 18), axis.title.y = element_blank())+ theme(legend.position = "none") +
  scale_fill_gradient2(low = "darkred", high = "midnightblue", mid = "white", midpoint = 0.5, limits=c(0.1, 1))+ facet_grid(.~time_f) 

g_B <- ggplot(tFA_mean_B, aes(x= locusname, y= cross_ord)) +
  ggtitle("Allele parent1 (SpB)") +geom_raster(aes (fill=mean_FA)) + theme(plot.title = element_text(color="grey40", size=18, face="bold"), axis.text.x = element_text(angle = 90, hjust = 1, size = 18), axis.title.x = element_blank(), strip.text = element_text(colour = 'black', face="bold", size = 18)) +
  theme(axis.text.y = element_text(hjust = 1, size = 18), axis.title.y = element_blank()) +
  scale_fill_gradient2(low = "darkred", high = "midnightblue", mid = "white", midpoint = 0.5, limits=c(0.1, 1), name= "Mean allele frequency")+ facet_grid(.~time_f) 

g_B_a <- g_B + theme(legend.position = "none")

cowplot::plot_grid(g_B_a, g_C, ncol=2, rel_heights= c(105, 75, 75, 42)) %>%
  cowplot::plot_grid(cowplot::get_legend(g_B), rel_widths = c(3, 0.8))

dev.off()

##################################################
                 #   THE END #
##################################################