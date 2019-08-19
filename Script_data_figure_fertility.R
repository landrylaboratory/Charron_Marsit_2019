##R version 3.5.2 was used to run the following script

#####Libraries for paper
##survival_2.43-3
library(survival)
##survminer_0.4.3 
library(survminer)
##beanplot_1.2
library(beanplot)
##jpeg_0.1-8 
library(jpeg)

####Colors for graphs
blue= rgb(0,0,(139/255))
red= rgb((238/255),0,0)
green= rgb(0,(139/255),0)
black= rgb(0,0,0)
orange = rgb((255/255),(128/255),(0/255))
purple = rgb((117/255),(112/255),(179/255))
grey = rgb((160/255),(160/255),(160/255))

####Load the data needed
viability= read.csv("Fertility_data.csv", header=T, sep=",")
ITC_data= read.csv("ITC_data.csv", header=T, sep=",")
Autodip_data= read.csv("Autodiploidization_data.csv", header=T, sep=",")
Bilan=read.csv("Survival_data.csv", header=T, sep=",")
generation= read.csv("Generations_data.csv", sep=",", header=T)

###FRS column 
####Calculate FRS scores for the evolution lines
T0_Tend_scores = viability$t_end-viability$t_ini
####append column 
viability= cbind(viability, T0_Tend_scores)
####calculate FRS scores for the ITC data
frs_ITC = ITC_data$P12-ITC_data$P0
#Append the column
ITC_data= cbind(ITC_data, frs_ITC)

##get FRS data for data source file
FRS_data = viability[c(1,2,3,10,12,31)]
colnames(FRS_data)[4]="Initial fertility"
colnames(FRS_data)[5]="Final fertility"
colnames(FRS_data)[6]="FRS"
write.table(FRS_data, "FRS_data.csv", sep=",", row.names=F)


####Subsets of lineages per cross type
L=subset(viability, viability$Cross_type == "BC" & viability$Ploidy_tend != "NA")
M=subset(viability, viability$Cross_type == "BA")
H=subset(viability, viability$Cross_type == "Bcer" & viability$Ploidy_tend != "NA")
VL=subset(viability, viability$Cross_type == "BB")
##subset per cross
L1=subset(viability, viability$Divergence == "L1" & viability$Ploidy_tend != "NA")
L2=subset(viability, viability$Divergence == "L2" & viability$Ploidy_tend != "NA")
M1=subset(viability, viability$Divergence == "M1" & viability$Ploidy_tend != "NA")
M2=subset(viability, viability$Divergence == "M2" & viability$Ploidy_tend != "NA")
H1=subset(viability, viability$Divergence == "H1" & viability$Ploidy_tend != "NA")
H2=subset(viability, viability$Divergence == "H2" & viability$Ploidy_tend != "NA")
VL1=subset(viability, viability$Divergence == "VL1" & viability$Ploidy_tend != "NA")
VL2=subset(viability, viability$Divergence == "VL2" & viability$Ploidy_tend != "NA")



####ITC data table reformated
ITC_table = data.frame(Characters = character())

for(i in 1:nrow(ITC_data)){

	line1 = cbind(ITC_data[i, c(1:4,8:12)], "ITC_0")
	colnames(line1)[c(4,5:10)] = c("Fertility", "sp_0", "sp_1", "sp_2", "sp_3", "sp_4","ITC_number")
	line2 = cbind(ITC_data[i, c(1:3,5,13:17)], "ITC_1")
	colnames(line2)[c(4,5:10)] = c("Fertility", "sp_0", "sp_1", "sp_2", "sp_3", "sp_4","ITC_number")
	line3 = cbind(ITC_data[i, c(1:3,6,18:22)], "ITC_6")
	colnames(line3)[c(4,5:10)] = c("Fertility", "sp_0", "sp_1", "sp_2", "sp_3", "sp_4","ITC_number")
	line4 = cbind(ITC_data[i, c(1:3,7,23:27)], "ITC_12")
	colnames(line4)[c(4,5:10)] = c("Fertility", "sp_0", "sp_1", "sp_2", "sp_3", "sp_4","ITC_number")
	

	ITC_table= rbind(ITC_table,line1,line2,line3,line4)

}

row.names(ITC_table) = 1:nrow(ITC_table)


#####Fertility data reformated

fertility = data.frame(Characters = character())

for(i in 1:nrow(viability)){

	line1 = cbind(viability[i, c(1:7,10,13:17,28)], "Tini")
	colnames(line1)[7:15] = c("Ploidy", "Fertility", "sp_0", "sp_1", "sp_2", "sp_3", "sp_4","Sporulation", "Time_point")
	line2 = cbind(viability[i, c(1:6,8,11,18:22,29)], "Tmid")
	colnames(line2)[7:15] = c("Ploidy", "Fertility", "sp_0", "sp_1", "sp_2", "sp_3", "sp_4","Sporulation", "Time_point")
	line3 = cbind(viability[i, c(1:6,9,12,23:27,30)], "Tend")
	colnames(line3)[7:15] = c("Ploidy", "Fertility", "sp_0", "sp_1", "sp_2", "sp_3", "sp_4","Sporulation", "Time_point")

	fertility= rbind(fertility,line1,line2,line3)

}

row.names(fertility) = 1:nrow(fertility)

###

test=aov(Fertility~Divergence*Time_point, data=fertility)



####FRS histograms data
frs_VL=hist(VL$T0_Tend_scores, xlim=c(-1,1), breaks=seq(-1,1, 0.1))
frs_L=hist(L$T0_Tend_scores, xlim=c(-1,1), breaks=seq(-1,1, 0.1))
frs_M=hist(M$T0_Tend_scores, xlim=c(-1,1), breaks=seq(-1,1, 0.1))
frs_H=hist(H$T0_Tend_scores, xlim=c(-1,1), breaks=seq(-1,1, 0.1))
frs_itc=hist(ITC_data$frs_ITC, xlim=c(-1,1), breaks=seq(-1,1, 0.1))

####Autodiploidization data
#average fertility for each category
avfert = apply(Autodip_data[,c(3:8)], 2, FUN=mean, na.rm=T)
#standard deviation for each category
sdfert = apply(Autodip_data[,c(3:8)], 2, FUN=sd, na.rm=T)
#n in each category
autodip_n = c(nrow(Autodip_data[complete.cases(Autodip_data[,3]),]), nrow(Autodip_data[complete.cases(Autodip_data[,4]),]),
nrow(Autodip_data[complete.cases(Autodip_data[,5]),]),nrow(Autodip_data[complete.cases(Autodip_data[,6]),]),
nrow(Autodip_data[complete.cases(Autodip_data[,7]),]),nrow(Autodip_data[complete.cases(Autodip_data[,8]),]))

####Survival curves data
passage_num=strsplit(as.character(Bilan[,4]), "P")
table_surv=cbind(Bilan[,1:4], as.numeric(unlist(passage_num)[seq(2,768*2,2)]))
colnames(table_surv)[5] = "Passage_lost"
table_surv=subset(table_surv, table_surv$Cross != "G")
status= c("")
for(i in 1:nrow(table_surv)){
	
	if (table_surv[i,5] == 35){
		status[i] = 0
	}else{
		status[i] = 1
}
}
status=as.numeric(status)



table_surv=cbind(table_surv, status)

Block = c(rep("ABC64",64), rep("ABC96", 32), rep("ABC64",64), rep("ABC96", 32),rep("ABC64",64), rep("ABC96", 32),rep("DEF64",64), rep("DEF96", 32),rep("DEF64",64), rep("DEF96", 32),rep("DEF64",64), rep("DEF96", 32), rep("H", 48), rep("I", 48))
Group = c(rep("LH",96),rep("VLM",96),rep("LH",96),rep("LH",96) ,rep("VLM",96),rep("LH",96),rep("VLM",96))
table_surv=cbind(table_surv, Block, Group)
table_surv=droplevels(table_surv)
survival_data = Surv(table_surv$Passage_lost, table_surv$status)
fit_type = survfit(survival_data~Cross_type, data=table_surv)
fit_cross = survfit(survival_data~Cross, data=table_surv)
fit_block = survfit(survival_data~Block, data=table_surv)
fit_group = survfit(survival_data~Group, data=table_surv)

diff_type=pairwise_survdiff(Surv(Passage_lost, status)~Cross_type, data=table_surv)
diff_cross=pairwise_survdiff(Surv(Passage_lost, status)~Cross, data=table_surv)
diff_block=pairwise_survdiff(Surv(Passage_lost, status)~Block, data=table_surv)
diff_group=pairwise_survdiff(Surv(Passage_lost, status)~Group, data=table_surv)


#Proportions of FRS scores and binomial test over and under 0
T0TeVLov=length(VL$T0_Tend_scores[(which(VL$T0_Tend_scores > 0))])
T0TeVLun=length(VL$T0_Tend_scores[(which(VL$T0_Tend_scores < 0))])
T0TeVL = binom.test(c(T0TeVLov, T0TeVLun))

T0TeLov=length(L$T0_Tend_scores[(which(L$T0_Tend_scores > 0))])
T0TeLun=length(L$T0_Tend_scores[(which(L$T0_Tend_scores < 0))])
T0TeL = binom.test(c(T0TeLov, T0TeLun))

T0TeMov=length(M$T0_Tend_scores[(which(M$T0_Tend_scores > 0))])
T0TeMun=length(M$T0_Tend_scores[(which(M$T0_Tend_scores < 0))])
T0TeM=binom.test(c(T0TeMov, T0TeMun))

T0TeHov=length(H$T0_Tend_scores[(which(H$T0_Tend_scores > 0))])
T0TeHun=length(H$T0_Tend_scores[(which(H$T0_Tend_scores < 0))])
T0TeH=binom.test(c(T0TeHov, T0TeHun))

T0TeITCov=length(ITC_data$frs_ITC[(which(ITC_data$frs_ITC > 0))])
T0TeITCun=length(ITC_data$frs_ITC[(which(ITC_data$frs_ITC < 0))])
T0TeITC=binom.test(c(T0TeITCov, T0TeITCun))


###fisher test between the timepoints 
#get a table with the viable/ non viable spores
##remove the NA for the lines that we couldn't get colonies from the glycerol
fertility_comp=fertility[c(as.numeric(which(complete.cases(fertility[,14])))),]

#create empty data frame
viable_spo_tab=data.frame()

#for each line in the fertility data frame calculate the number of viable/inviable spores
for (i in 1:nrow(fertility_comp)){
	
	if (fertility_comp[i,14] == "yes"){
	live_spores=sum(fertility_comp[i,c(9:13)]*c(0,1,2,3,4))
	dead_spores=96-live_spores
	viable_spo_line=cbind(fertility_comp[i,c(1:3,14,15)], live_spores, dead_spores)
	viable_spo_tab=rbind(viable_spo_tab, viable_spo_line)
	}else if (fertility_comp[i,14] == "no"){
	live_spores="NA"
	dead_spores="NA"
	viable_spo_line=cbind(fertility_comp[i,c(1:3,14,15)], live_spores, dead_spores)
	viable_spo_tab=rbind(viable_spo_tab, viable_spo_line)
}
}

###fisher test between the timepoints (Tini-Tend)

#get a table with the viable/ non viable spores
##remove the NA for the lines taht couldn't get colonies from the glycerol
fertility_comp=fertility[c(as.numeric(which(complete.cases(fertility[,14])))),]
###get a vector of line names for the loop to go through
lines=unique(viable_spo_tab[,1])
###create empty data frame to put the data in
pvalue_fisher=data.frame()

##loop through the line names
for (i in lines){

	sub=subset(viable_spo_tab, viable_spo_tab[,1] == i)
	
	Tini = subset(sub[,c(6,7)], sub[,5] == "Tini")
	Tmid = subset(sub[,c(6,7)], sub[,5] == "Tmid")
	Tend = subset(sub[,c(6,7)], sub[,5] == "Tend")
	
	if(dim(Tini)[1] == 1 & dim(Tmid)[1]== 1){
		I_M = matrix (nrow=2,ncol=2, c(as.numeric(Tini[1,1]), as.numeric(Tini[1,2]),as.numeric(Tmid[1,1]), as.numeric(Tmid[1,2])), byrow=T)
	}else{
		I_M=matrix (nrow=2,ncol=2, c(NA, NA, NA,NA), byrow=T)
	}
	if(anyNA(I_M) == "FALSE"){
		I_M_test = fisher.test(I_M)
		pval_IM = I_M_test$p.value		
	}else{
		pval_IM="NA"
	}
	
	if(dim(Tmid)[1] == 1 & dim(Tend)[1]== 1){
		M_E = matrix (nrow=2,ncol=2, c(as.numeric(Tmid[1,1]), as.numeric(Tmid[1,2]),as.numeric(Tend[1,1]), as.numeric(Tend[1,2])), byrow=T)
	}else{	
		M_E=matrix (nrow=2,ncol=2, c(NA, NA, NA,NA), byrow=T)
	}
	if(anyNA(M_E) == "FALSE"){
		M_E_test = fisher.test(M_E)
		pval_ME = M_E_test$p.value		
	}else{
		pval_ME="NA"
	}

	if(dim(Tini)[1] == 1 & dim(Tend)[1]== 1){
		I_E = matrix (nrow=2,ncol=2, c(as.numeric(Tini[1,1]), as.numeric(Tini[1,2]),as.numeric(Tend[1,1]), as.numeric(Tend[1,2])), byrow=T)
	}else{
		I_E=matrix (nrow=2,ncol=2, c(NA, NA, NA,NA), byrow=T)
	}
	if(anyNA(I_E) == "FALSE"){
		I_E_test = fisher.test(I_E)
		pval_IE = I_E_test$p.value		
	}else{
		pval_IE="NA"
	}
		
	test_line=cbind(i, pval_IM, pval_ME, pval_IE)
	pvalue_fisher=rbind(pvalue_fisher, test_line)
	
}

pvalue_fisher=cbind(as.numeric(as.character(pvalue_fisher[,2])), as.numeric(as.character(pvalue_fisher[,3])), as.numeric(as.character(pvalue_fisher[,4])))
row.names(pvalue_fisher)=lines
colnames(pvalue_fisher) = c("Tini-Tmid", "Tmid-Tend", "Tini-Tend")

P_fisher = as.numeric(c(pvalue_fisher[,1], pvalue_fisher[,2], pvalue_fisher[,3]))

P_fisher_corr = p.adjust(P_fisher, method="fdr")

pvalue_fisher_corr = cbind(pvalue_fisher[,1], P_fisher_corr[1:214], pvalue_fisher[,2], P_fisher_corr[215:428], pvalue_fisher[,3], P_fisher_corr[429:642])
colnames(pvalue_fisher_corr) = c("Tini-Tmid","Tini-Tmid_corr", "Tmid-Tend", "Tmid-Tend_corr", "Tini-Tend", "Tini-Tend_corr")

write.table(pvalue_fisher_corr, "pval_fertility.csv", sep=",", row.names=T)



####################################
############Main Figures############
####################################


####Figure 1 - Overview Figure######

###survival curves
#pdf("survival.pdf", height=8, width=8)
pdf("survival_BA_uncollapsed.pdf", height=8, width=8)
par(mar=c(18,5,1,1))
plot.new()
plot.window(xlim=c(0,35), ylim=c(.5,1))
title(xlab="Passage", ylab="Proportion of surviving lines (%)", cex.lab=1.5)
axis(1, cex.axis=1.5, lwd=2.5)
axis(2, las=1, cex.axis=1.5, lwd=2.5)
#lines(fit_type, lwd=2, col=c(green, red, blue, black), lty=c(1,1,1,1))
lines(fit_type[c(2:4)], lwd=2, col=c(red,blue,black), lty=c(1,1,1,1), conf.int=F)
lines(fit_cross[c(2,5)], lwd=2, col=c(green), lty=c(1,2), conf.int=F)
#lines(fit_type, lwd=2, col=c(rgb(0,(139/255),0,.01), rgb((238/255),0,0,.01), rgb(0,0,(139/255),.01), rgb(0,0,0,.01)), lty=c(1,1,1,1), conf.times=seq(0,36,.01), conf.cap=0)
lines(fit_type[c(2:4)], lwd=2, col=c(rgb((238/255),0,0,.01), rgb(0,0,(139/255),.01), rgb(0,0,0,.01)), lty=c(1,1,1,1), conf.times=seq(0,36,.01), conf.cap=0)
lines(fit_cross[c(2,5)], lwd=2, col=c(rgb(0,(139/255),0,.01)), lty=c(1,1,1,1), conf.times=seq(0,36,.01), conf.cap=0)

par(lty=0)
for(i in c(0,1,seq(4,31,3),35)){
	segments(i,.4,i,1, lty=1, lwd=6, col=rgb(160/255,160/255,160/255,0.2))
}
par(xpd=NA)
#legend(9, .39, ncol=4L, title="Log-rank test P-values", legend=c("", "      L", "      M", "      H", "   VL", c("< 0.01","< 0.05","< 0.01"), "    L", c("    -","< 0.01","   .667"), "    M", c("    -","    -","< 0.01")), bg="white", cex=1.3)
#legend(12, .35, ncol=5L, title="Log-rank test P-values", legend=c("","L","M1","M2","H", "VL","< 0.01","< 0.05", "n.s.","< 0.01","L"," -","n.s.","< 0.01","n.s.","M1"," -"," -","< 0.01","n.s.","M2"," -"," -"," -","< 0.01"), bg="white", cex=1.2, text.width=2)
legend(-1, .375, ncol=8L, title="          Log-rank test P-values", legend=c("","VL2","L1","L2","M1","M2","H1","H2","  VL1","  0.42","< 0.01","< 0.01","< 0.05","  0.24","< 0.01","< 0.01","  VL2","","< 0.01","< 0.01","< 0.05","  0.75","< 0.01","< 0.01","  L1","","","  0.63","  0.15","< 0.01","  0.81","  0.42","  L2","","","","  0.35","< 0.01","  0.75","  0.75","  M1","","","","","< 0.05","  0.35","  0.47","  M2","","","","","","< 0.01","< 0.01","  H1","","","","","","","  0.60"), bg="white", cex=1.4, text.width=1.7)
legend(-1,.75, title="Crosses", c("VL (96)", "L (192)", "M1 (96)","M2 (96)", "H (192)"), col=c(red, blue, green, green, black),lty=c(1,1,1,2,1), lwd=2, bg="white", cex=1.2)
par(xpd=0)

dev.off()
#mtext("A", side=3, line=0.5, adj=-.15, cex=3)

###non sporulating strains data

T0lostspo = subset(viability, viability$Sporulation_t0 == "no")
T0lostspo$Cross_type = factor(T0lostspo$Cross_type, levels=c("BB", "BC", "BA", "Bcer"))
Tmlostspo = subset(viability, viability$Sporulation_tmid == "no")
Tmlostspo$Cross_type = factor(Tmlostspo$Cross_type, levels=c("BB", "BC", "BA", "Bcer"))
Telostspo = subset(viability, viability$Sporulation_tend == "no")
Telostspo$Cross_type = factor(Telostspo$Cross_type, levels=c("BB", "BC", "BA", "Bcer"))

####GLM for no spo data

no_spo = viability[,c(2,30)]

div=c(rep(.02,72), rep(.04, 46), rep(.15, 48), rep(.004,48))

no_spo = cbind(no_spo, div)

colnames(no_spo) = c("Cross_type", "Sporulation_tend", "Gen_div")

no_spo=no_spo[c(-29,-130),]

test = glm(Sporulation_tend~Gen_div, data=no_spo, family=binomial)	
summary(test)

###fisher tests on the two lineages per cross

table(Telostspo$B_Strain, Telostspo$Cross_type)

no_spo_BB=fisher.test(matrix(nrow=2,ncol=2, c(1,23,2,22), byrow=T))
no_spo_BC=fisher.test(matrix(nrow=2,ncol=2, c(3,33,5,31), byrow=T))
no_spo_BA=fisher.test(matrix(nrow=2,ncol=2, c(1,23,7,17), byrow=T))
no_spo_Bcer=fisher.test(matrix(nrow=2,ncol=2, c(7,17,11,13), byrow=T))

fishpval_nospo=c(no_spo_BB$p.value,no_spo_BC$p.value, no_spo_BA$p.value, no_spo_Bcer$p.value)
p.adjust(fishpval_nospo, method="fdr")
#1.00000000000000000 0.94680209579077146 0.19169452456523581 0.74293762433441346


####barplots
#pdf("Figure_supp_NoSpo_roman_vertical.pdf", height=16, width=4)
#pdf("Figure_supp_NoSpo_roman_vertical.pdf", height=4, width=8)
pdf("Figure_supp_NoSpo.pdf", height=5, width=8)
#par(mfrow=c(2,1), lty=0, mar=c(5,5,4,4))
par(mfrow=c(1,3), lty=0, mar=c(5,5,4,4))

bp=barplot(table(T0lostspo$B_Strain, T0lostspo$Cross_type)/c(48,48,72,72,46,46,48,48), beside=F, ylim=c(0,.5), las=1, col=c("grey"), axes=F, names.arg=c(expression("VL"["div"]),expression("L"["div"]),expression("M"["div"]),expression("H"["div"])), cex.names=1.4)
axis(2, lwd=2.5, las=1, cex.axis=1.4)
title(main=expression("T"["ini"]),ylab="Proportion of lines", cex.lab=1.4, cex.main=2.5)
text(bp, .04,  labels=c("48","72","46","48"), pos=3, cex=1.4)
#legend("topleft", title=expression(paste(italic("SpB"),"parental strain")), bty="n", c("UWOPS_91_202", "MSH604"), fill=c("black", "grey"))
#mtext("i", side=3, line=1.2, adj=-.2, cex=2)
#mtext("i", side=3, line=1.2, adj=-.3, cex=2)

bp=barplot(table(Tmlostspo$B_Strain, Tmlostspo$Cross_type)/c(48,48,72,72,46,46,48,48), beside=F, ylim=c(0,.5), las=1, col=c("grey"), axes=F, names.arg=c(expression("VL"["div"]),expression("L"["div"]),expression("M"["div"]),expression("H"["div"])), cex.names=1.4)
axis(2, lwd=2.5, las=1, cex.axis=1.4)
title(main=expression("T"["mid"]), ylab="Proportion of lines", cex.lab=1.4, cex.main=2.5)
text(bp, summary(Tmlostspo$Cross_type)/c(48,72,46,48)+.04,  labels=c("48","72","46","48"), pos=3, cex=1.4)
#mtext("ii", side=3, line=1.2, adj=-.2, cex=2)
#mtext("ii", side=3, line=1.2, adj=-.3, cex=2)

bp=barplot(table(Telostspo$B_Strain, Telostspo$Cross_type)/c(48,48,72,72,46,46,48,48), beside=F, ylim=c(0,.5), las=1, col=c("grey"), axes=F, names.arg=c(expression("VL"["div"]),expression("L"["div"]),expression("M"["div"]),expression("H"["div"])), cex.names=1.4)
axis(2, lwd=2.5, las=1, cex.axis=1.4)
title(main=expression("T"["end"]),ylab="Proportion of lines", cex.lab=1.4, cex.main=2.5)
text(bp, summary(Telostspo$Cross_type)/c(48,72,46,48)+.04,  labels=c("48","72","46","48"), pos=3, cex=1.4)
#mtext("iii", side=3, line=1.2, adj=-.2, cex=2)
#mtext("iii", side=3, line=1.2, adj=-.3, cex=2)

dev.off()


####Figure 2 - Spore photos, FRS data and Autodiploidization####
###create the layout of the figure
pdf("Figure_2.pdf", height=7, width=24.5)
fig = matrix(c(rep(1,4),2,4,6,3,5,0, rep(7,4)), byrow=F, nrow=2, ncol=7)
layout(fig)

#####Insert images of the spores as Panel A
img = readJPEG("Spore_img2.jpg")

plot.new()
par(mar=c(3,6,5,0), xpd=NA)
plot.window(xlim=c(0,1), ylim=c(0,1))
axis(1, at=c(.3, .7), labels=c(expression("T"["ini"]), expression("T"["end"])), lwd=0, cex.axis=3)
rasterImage(img, .125,0,.875,1.14)
text(0,.18, "FRS\n > 0 \n(L1)", cex=2)
text(0,.58, "FRS\n = 0 \n(L2)", cex=2)
text(0,.95, "FRS\n < 0 \n(VL1)", cex=2)
mtext("A", side=3, adj=0, line=1.5, cex=3)

####Plot the histogram data for each cross type and the FRS densisty lines

par(mar=c(5,5,4,4), lty=1, xpd=0)
plot.new()
plot.window(xlim=c(-1,1), ylim=c(0,1))
axis(1, at=c(-1,-.5 , 0, .5, 1), lwd=2.5, cex.axis=1.5)
axis(2, lwd=2.5, las=1, cex.axis=1.5)			
rect(frs_VL$breaks[-length(frs_VL$breaks)],0,frs_VL$breaks[-1],frs_VL$density*.1, lwd=2)
text(0,.95, paste("P =", round(as.numeric(T0TeVL[3]),3)), cex=2)
par(xpd=NA)
polygon(c(seq(-1,1,1), seq(1,-1,-1)), c(rep(1.05,3), rep(1.27,3)), col=red, lty=0)
title(main=paste("VL mitotic (",length(which(VL$T0_Tend_scores != "NA")),")", sep=""), ylab="Proportion of lines (%)", xlab="FRS", cex.main=2.8, cex.lab=2, col.main="White")
mtext("B", side=3, line=0.5, adj=-.5, cex=3)
par(xpd=0)

plot.new()
plot.window(xlim=c(-1,1), ylim=c(0,1))
par(mar=c(5,5,4,4), lty=1)
axis(1, at=c(-1,-.5 , 0, .5, 1), lwd=2.5, cex.axis=1.5)
axis(2, lwd=2.5, las=1, cex.axis=1.5)			
rect(frs_H$breaks[-length(frs_M$breaks)],0,frs_H$breaks[-1],frs_H$density*.1, lwd=2)
text(0,.95, paste("P =", round(as.numeric(T0TeH[3]),3)), cex=2)
arrows(.65, .3, .65, .08, length=.1, lwd=2, col=red)
par(xpd=NA)
polygon(c(seq(-1,1,1), seq(1,-1,-1)), c(rep(1.05,3), rep(1.27,3)), col=black, lty=0)
title(main=paste("H mitotic (",length(which(H$T0_Tend_scores != "NA")),")", sep=""), ylab="Proportion of lines (%)", xlab="FRS", cex.main=2.8, cex.lab=2, col.main="White")
par(xpd=0)

plot.new()
plot.window(xlim=c(-1,1), ylim=c(0,1))
par(mar=c(5,5,4,4), lty=1)
axis(1, at=c(-1,-.5 , 0, .5, 1), lwd=2.5, cex.axis=1.5)
axis(2, lwd=2.5, las=1, cex.axis=1.5)			
rect(frs_L$breaks[-length(frs_L$breaks)],0,frs_L$breaks[-1],frs_L$density*.1, lwd=2)
text(0,.95, paste("P =", round(as.numeric(T0TeL[3]),3)), cex=2)
arrows(.75, .3, .75, .08, length=.1, lwd=2, col=red)
arrows(.55, .3, .55, .08, length=.1, lwd=2, col=red)
arrows(.45, .3, .45, .08, length=.1, lwd=2, col=red)
arrows(.35, .3, .35, .08, length=.1, lwd=2, col=red)
arrows(.25, .3, .25, .08, length=.1, lwd=2, col=red)

par(xpd=NA)
polygon(c(seq(-1,1,1), seq(1,-1,-1)), c(rep(1.05,3), rep(1.27,3)), col=blue, lty=0)
title(main=paste("L mitotic (",length(which(L$T0_Tend_scores != "NA")),")", sep=""), ylab="Proportion of lines (%)", xlab="FRS", cex.main=2.8, cex.lab=2, col.main="White")
par(xpd=0)

plot.new()
plot.window(xlim=c(-1,1), ylim=c(0,1))
par(mar=c(5,5,4,4), lty=1)
axis(1, at=c(-1,-.5 , 0, .5, 1), lwd=2.5, cex.axis=1.5)
axis(2, lwd=2.5, las=1, cex.axis=1.5)			
rect(frs_itc$breaks[-length(frs_itc$breaks)],0,frs_itc$breaks[-1],frs_itc$density*.1, lwd=2)
text(0,.95, paste("P < 0.01"), cex=2)
par(xpd=NA)
polygon(c(seq(-1,1,1), seq(1,-1,-1)), c(rep(1.05,3), rep(1.27,3)), col=black, lty=0)
title(main="ITC meiotic (12)", ylab="Proportion of lines (%)", xlab="FRS", cex.main=2.5, cex.lab=2, col.main="White")
par(xpd=0)

plot.new()
plot.window(xlim=c(-1,1), ylim=c(0,1))
par(mar=c(5,5,4,4), lty=1)
axis(1, at=c(-1,-.5 , 0, .5, 1), lwd=2.5, cex.axis=1.5)
axis(2, lwd=2.5, las=1, cex.axis=1.5)			
rect(frs_M$breaks[-length(frs_M$breaks)],0,frs_M$breaks[-1],frs_M$density*.1, lwd=2)
text(0,.95, paste("P =", round(as.numeric(T0TeM[3]),3)), cex=2)
arrows(.25, .3, .25, .08, length=.1, lwd=2, col=red)
par(xpd=NA)
polygon(c(seq(-1,1,1), seq(1,-1,-1)), c(rep(1.05,3), rep(1.27,3)), col=green, lty=0)
title(main=paste("M mitotic (",length(which(M$T0_Tend_scores != "NA")),")", sep=""), ylab="Proportion of lines (%)", xlab="FRS", cex.main=2.8, cex.lab=2, col.main="White")
par(xpd=0)


####FRS density plots
#plot.new()
#plot.window(xlim=c(-1,1), ylim=c(0,1))
#par(mar=c(5,5,4,4), lty=1)
#title(xlab="FRS",ylab="Proportion of lines (%)", cex.lab=2)
#axis(1, at=c(-1,-.5 , 0, .5, 1), lwd=2.5, cex.axis=1.5)
#axis(2, lwd=2.5, las=1, cex.axis=1.5)			
#legend(-1,1 , c("VL (48)", "L (72)", "M (46)", "H (48)","ITC (16)"), col=c(red, blue, green, black, black), lty=c(1,1,1,1,2), lwd=2.5, bty="n", cex=1.4)


#lines(frs_VL$breaks,c(frs_VL$density*.1,0), lwd=2, col=red)
#lines(frs_L$breaks,c(frs_L$density*.1,0), lwd=2, col=blue)
#lines(frs_M$breaks,c(frs_M$density*.1,0), lwd=2, col=green)
#lines(frs_H$breaks,c(frs_H$density*.1,0), lwd=2, col=black)
#lines(frs_itc$breaks,c(frs_itc$density*.1,0), lwd=2, lty=2, col=black)

#high recovery region#
#polygon(c(seq(0.30,.7,.1), seq(.7,0.30,-.1)), c(rep(.9,5), rep(0,5)), col=rgb(217/255,217/255,217/255, .3), lty=0)
#text(.5,.92, "High recovery", cex=1.4)


####Autodiploidization barplot
#par(lty=0, xpd=NA, mar=c(5,5,4,4))


#plot the barplot
#bp = barplot(avfert, ylim=c(0,1.2), names=c("Tini", "Tmid", "Tend","Tini", "Tmid", "Tend"), las=1, col="grey", cex.names=1.5, ylab="Fertility", cex.lab=1.5, cex.axis=1.5, lwd=2.5)
#axis(1, lty=1, lwd=2.5, labels=F, at=F)

#add segments to indicate standard deviations
#segments(bp, avfert+sdfert, bp,avfert-sdfert, lty=1,lwd=2)
#segments(bp-.1, avfert+sdfert, bp+.1,avfert+sdfert, lty=1,lwd=2)
#segments(bp-.1, avfert-sdfert, bp+.1,avfert-sdfert, lty=1,lwd=2)

#add the n over the std dev bars
#text(bp[1:3],max((avfert+sdfert)[1:3])+.1, labels=paste("(", autodip_n[1:3], ")", sep=""), cex=1.4)
#text(bp[4:6],max((avfert+sdfert)[4:6])+.1, labels=paste("(", autodip_n[4:6], ")", sep=""), cex=1.4)

#add segments for the ID of the categories
#segments(0,-.05, 3.7,-.05, lty=1, lwd=2.5, col="black")
#segments(3.8,-.05, 7.5,-.05, lty=1, lwd=2.5)

#add text in margin
#mtext("parental", side=1, line=3.4, adj=0.2, cex=1.2)
#mtext("autodiploidized", side=1, line=3.4, adj=.9, cex=1.2)

#dev.off()

#####boxplot version
par(xpd=NA, mar=c(6,6,4,4))
boxplot(Autodip_data[3:8], lwd=3, lty=1, frame.plot=0, ylim=c(0,1), axes=0, outline=0)
stripchart(Autodip_data[3:8], vertical=T, add=T, pch=16, method="jitter", cex=2)
axis(1, lwd=3, labels=NA, at=c(1,2,3,4,5,6), cex.axis=2)
axis(1, lwd=0, line=1, labels= c(expression("T"["ini"]), expression("T"["mid"]), expression("T"["end"]), expression("T"["ini"]), expression("T"["mid"]), expression("T"["end"])), at=c(1,2,3,4,5,6), cex.axis=2.5)
axis(2, lwd=3, las=1, cex.axis=2)
text(c(1,2,3,4,5,6), c(rep(.025,6)), labels=paste("(", autodip_n[1:6], ")", sep=""), cex=2.5)
title(ylab="Spore survival (%)", cex.lab=2, line=4)

segments(.65,-.13, 3.5,-.13, lty=1, lwd=3)
segments(3.7,-.13, 6.4,-.13, lty=1, lwd=3)

mtext("parental", side=1, line=4.5, adj=0.25, cex=1.5)
mtext("autodiploidized", side=1, line=4.5, adj=.85, cex=1.5)
mtext("C", side=3, line=0.5, adj=-.15, cex=3)
dev.off()


####Figure 4 - Slopegraphs integrating Fertility and Ploidy

#pdf("Figure_4_colorborder.pdf", height=9, width=16.2)
pdf("Figure_4_colorborder_pval.pdf", height=9, width=16.2)
par(mfrow=c(1,4), lty=1, mar=c(4,7,4,1.5))

####SpBxSpB
####Call new plot
plot.new()
plot.window(xlim=c(-.05,1.05), ylim=c(0,1))
axis(1, at=c(0, .5, 1), lwd=2.5, cex.axis=2.5, labels=NA)
axis(1, labels=c(expression("T"["ini"]), expression("T"["mid"]), expression("T"["end"])), at=c(0, .5, 1), cex.axis=3, line=1, lwd=0)
axis(2, lwd=2.5, las=1, cex.axis=2.5)			
legend(0.40,1.05 , c("VL1", "VL2"), col=c(rgb(0,0,0,0.4), rgb(0,0,0,0.4)), lty=c(1,2), bty="n", cex=2, lwd=2)
legend(0,1.05, c("2n", "3n", "4n"), pch=c(1, 16, 16), col=c(grey,orange,purple), cex=2, bty="n", pt.cex=2.5)
#Draw Polygon outside of graph area (remove "#" for display) 
par(xpd=NA)
#grey version
#polygon(c(seq(-0.05,1.05,.1), seq(1.05,-0.05,-.1)), c(rep(1.04,12), rep(1.1,12)), col=rgb(217/255,217/255,217/255), lty=0)
#title(main="VL", ylab="Fertility", cex.main=4, cex.lab=2)
#Color version
polygon(c(seq(-0.05,1.05,.1), seq(1.05,-0.05,-.1)), c(rep(1.04,12), rep(1.2,12)), col=red, lty=0)
title(main=expression("VL"["div"]),cex.main=4, cex.lab=2.5, col.main="White")
mtext(side=2,"Fertility", line=4.5, cex=2)
par(xpd=0)

####Loop for generating the different curves (one per line of the specified subset)
for(i in 1:nrow(VL)){
	#add random spread of the X values for better vizualization of all the data
	jit = jitter(rep(0, nrow(VL)))
	#Use line number "i" data to draw the segments and points
	lineage = VL[i,]
	#if the strain contain H it the line will be full
	if(grepl("H",lineage$Strain) == TRUE){
		ltype = 2
	#if not it will be dashed
	} else{
		ltype = 1
}
	#if the lineage gets to P35, the line fill be full
	#if(lineage$Last_glycerol == "P35"){
		#ltype = 1
	#if not it will be dashed
	#}else{
		#ltype = 2		
#}
	#if diploid at Tini, use empty grey circles
	if(lineage$Ploidy_tini == 2){
		pchini= 1
		colini = rgb(160/255,160/255,160/255)
	#if triploid use full orange circles
	} else if (lineage$Ploidy_tini == 3){
		pchini= 16
		colini = rgb(255/255,128/255,0/255)
	#if tetraploid use full purple circles
	} else if (lineage$Ploidy_tini == 4){
		pchini= 16
		colini = rgb(117/255,112/255,179/255)
}
	#if diploid at Tmid, use empty grey circles
	if(lineage$Ploidy_tmid == 2){
		pchmid= 1
		colmid = rgb(160/255,160/255,160/255)
	#if triploid use full orange circles
	} else if (lineage$Ploidy_tmid == 3){
		pchmid= 16
		colmid = rgb(255/255,128/255,0/255)
	#if tetraploid use full purple circles
	} else if (lineage$Ploidy_tmid == 4){
		pchmid= 16
		colmid = rgb(117/255,112/255,179/255)

}
	#if diploid at Tend, use empty grey circles
	if(lineage$Ploidy_tend == 2){
		pchend= 1
		colend = rgb(160/255,160/255,160/255)
	#if triploid use full orange circles
	} else if (lineage$Ploidy_tend == 3){
		pchend= 16
		colend = rgb(255/255,128/255,0/255)
	#if tetraploid full purple circles
	} else if (lineage$Ploidy_tend == 4){
		pchend= 16
		colend = rgb(117/255,112/255,179/255)

}
	#draw the segments between Tini and Tmid
	strain = which(row.names(pvalue_fisher_corr) == VL[i,1])
	if(pvalue_fisher_corr[strain, 6] <= .05 & is.na(pvalue_fisher_corr[strain, 6]) == "FALSE" ){
	segments(jit[i],lineage$t_ini, jit[i]+0.5, lineage$t_mid, lwd=2.5, col=rgb(0,0,0,0.7), lty=ltype)
	}else{
	segments(jit[i],lineage$t_ini, jit[i]+0.5, lineage$t_mid, lwd=1, col=rgb(0,0,0,.2), lty=ltype)
	}
	#draw the segments between Tmid and Tend
	if(pvalue_fisher_corr[strain, 6] <= .05 & is.na(pvalue_fisher_corr[strain, 6]) == "FALSE" ){
	segments(jit[i]+.5, lineage$t_mid, jit[i]+1, lineage$t_end, lwd=2.5, col=rgb(0,0,0,.7), lty=ltype)
	}else{
	segments(jit[i]+.5, lineage$t_mid, jit[i]+1, lineage$t_end, lwd=1, col=rgb(0,0,0,.2), lty=ltype)
	}
	#draw the point at Tini
	points (jit[i], lineage$t_ini, cex=2.5, pch=pchini, col=colini)
	#draw the point at Tmid
	points (jit[i]+.5, lineage$t_mid, cex=2.5, pch=pchmid, col=colmid)
	#draw the point at Tend
	points (jit[i]+1, lineage$t_end, cex=2.5, pch=pchend, col=colend)


}

####SpBxSpC
plot.new()
plot.window(xlim=c(-.05,1.05), ylim=c(0,1))
axis(1, at=c(0, .5, 1), lwd=2.5, cex.axis=2.5, labels=NA)
axis(1, labels=c(expression("T"["ini"]), expression("T"["mid"]), expression("T"["end"])), at=c(0, .5, 1), cex.axis=3, line=1, lwd=0)
#axis(2, lwd=2.5, las=1, cex.axis=3)			
legend(0.4,1.05 , c("L1", "L2"), col=c(rgb(0,0,0,0.4), rgb(0,0,0,0.4)), lty=c(1,2), lwd=2, bty="n", cex=2)
#legend(0,.985 , c("Lost"), col= rgb(0,0,0,.4), lty=2, lwd=2, bty="n", cex=2)
legend(0,1.05, c("2n", "3n", "4n"), pch=c(1, 16, 16), col=c(grey,orange,purple), cex=2, bty="n", pt.cex=2.5)
par(xpd=NA)
#polygon(c(seq(-0.05,1.05,.1), seq(1.05,-0.05,-.1)), c(rep(1.04,12), rep(1.1,12)), col=rgb(217/255,217/255,217/255), lty=0)
#title(main=expression("L"["div"]), ylab="Fertility", cex.main=4, cex.lab=2)
polygon(c(seq(-0.05,1.05,.1), seq(1.05,-0.05,-.1)), c(rep(1.04,12), rep(1.2,12)), col=blue, lty=0)
title(main="L", cex.main=4, cex.lab=2, col.main="White")
par(xpd=0)


for(i in 1:nrow(L)){

	jit = jitter(rep(0, nrow(L)))
	lineage = L[i,]
	if(grepl("A",lineage$Strain) == TRUE){
		ltype = 1
	} else{
		ltype = 2
}
	#if(lineage$Last_glycerol == "P35"){
		#ltype = 1
	#}else{
		#ltype = 2		
#}

	if(lineage$Ploidy_tini == 2){
		pchini= 1
		colini = grey
	} else if (lineage$Ploidy_tini == 3){
		pchini= 16
		colini = orange
	} else if (lineage$Ploidy_tini == 4){
		pchini= 16
		colini = purple
}
	if(lineage$Ploidy_tmid == 2){
		pchmid= 1
		colmid = grey
	} else if (lineage$Ploidy_tmid == 3){
		pchmid= 16
		colmid = orange
	} else if (lineage$Ploidy_tmid == 4){
		pchmid= 16
		colmid = purple

}
	if(lineage$Ploidy_tend == 2){
		pchend= 1
		colend = grey
	} else if (lineage$Ploidy_tend == 3){
		pchend= 16
		colend = orange
	} else if (lineage$Ploidy_tend == 4){
		pchend= 16
		colend = purple

}
	strain = which(row.names(pvalue_fisher_corr) == L[i,1])
	if(pvalue_fisher_corr[strain, 6] <= .05 & is.na(pvalue_fisher_corr[strain, 6]) == "FALSE" ){
	segments(jit[i],lineage$t_ini, jit[i]+0.5, lineage$t_mid, lwd=2.5, col=rgb(0,0,0,0.7), lty=ltype)
	}else{
	segments(jit[i],lineage$t_ini, jit[i]+0.5, lineage$t_mid, lwd=1, col=rgb(0,0,0,.2), lty=ltype)
	}
	if(pvalue_fisher_corr[strain, 6] <= .05 & is.na(pvalue_fisher_corr[strain, 6]) == "FALSE" ){
	segments(jit[i]+.5, lineage$t_mid, jit[i]+1, lineage$t_end, lwd=2.5, col=rgb(0,0,0,.7), lty=ltype)
	}else{
	segments(jit[i]+.5, lineage$t_mid, jit[i]+1, lineage$t_end, lwd=1, col=rgb(0,0,0,.2), lty=ltype)
	}
	points (jit[i], lineage$t_ini, cex=2.5, pch=pchini, col=colini)
	points (jit[i]+.5, lineage$t_mid, cex=2.5, pch=pchmid, col=colmid)
	points (jit[i]+1, lineage$t_end, cex=2.5, pch=pchend, col=colend)


}



####SpBxSpA
plot.new()
plot.window(xlim=c(-.05,1), ylim=c(0,1))
axis(1, at=c(0, .5, 1), lwd=2.5, cex.axis=2.5, labels=NA)
axis(1, labels=c(expression("T"["ini"]), expression("T"["mid"]), expression("T"["end"])), at=c(0, .5, 1), cex.axis=3, line=1, lwd=0)
#axis(2, lwd=2.5, las=1, cex.axis=3)			
legend(0.4,1.05 , c("M1", "M2"), col=c(rgb(0,0,0,0.4), rgb(0,0,0,0.4)), lty=c(1,2), lwd=2, bty="n", cex=2)
#legend(0,.985 , c("Lost"), col= rgb(0,0,0,.4), lty=2, lwd=2, bty="n", cex=2)
legend(0,1.05, c("2n", "3n", "4n"), pch=c(1, 16, 16), col=c(grey,orange,purple), cex=2, bty="n", pt.cex=2.5)
par(xpd=NA)
#polygon(c(seq(-0.05,1.05,.1), seq(1.05,-0.05,-.1)), c(rep(1.04,12), rep(1.1,12)), col=rgb(217/255,217/255,217/255), lty=0)
#title(main="M", ylab="Fertility", cex.main=4, cex.lab=2)
polygon(c(seq(-0.05,1.05,.1), seq(1.05,-0.05,-.1)), c(rep(1.04,12), rep(1.2,12)), col=green, lty=0)
title(main=expression("M"["div"]), cex.main=4, cex.lab=2, col.main="white")
par(xpd=0)


for(i in 1:nrow(M)){

	jit = jitter(rep(0, nrow(M)))
	lineage = M[i,]
	if(grepl("B",lineage$Strain) == TRUE){
		ltype = 1
	} else{
		ltype = 2
}
	#if(lineage$Last_glycerol == "P35"){
		#ltype = 1
	#}else{
		#ltype = 2		
#}

	
	if(lineage$Ploidy_tini == 2){
		pchini= 1
		colini = grey
	} else if (lineage$Ploidy_tini == 3){
		pchini= 16
		colini = orange
	} else if (lineage$Ploidy_tini == 4){
		pchini= 16
		colini = purple
}
	if(lineage$Ploidy_tmid == 2){
		pchmid= 1
		colmid = grey
	} else if (lineage$Ploidy_tmid == 3){
		pchmid= 16
		colmid = orange
	} else if (lineage$Ploidy_tmid == 4){
		pchmid= 16
		colmid = purple

}
	if(lineage$Ploidy_tend == 2){
		pchend= 1
		colend = grey
	} else if (lineage$Ploidy_tend == 3){
		pchend= 16
		colend = orange
	} else if (lineage$Ploidy_tend == 4){
		pchend= 16
		colend = purple

}
	strain = which(row.names(pvalue_fisher_corr) == M[i,1])
	if(pvalue_fisher_corr[strain, 6] <= .05 & is.na(pvalue_fisher_corr[strain, 6]) == "FALSE" ){
	segments(jit[i],lineage$t_ini, jit[i]+0.5, lineage$t_mid, lwd=2.5, col=rgb(0,0,0,0.7), lty=ltype)
	}else{
	segments(jit[i],lineage$t_ini, jit[i]+0.5, lineage$t_mid, lwd=1, col=rgb(0,0,0,.2), lty=ltype)
	}
	if(pvalue_fisher_corr[strain, 6] <= .05 & is.na(pvalue_fisher_corr[strain, 6]) == "FALSE" ){
	segments(jit[i]+.5, lineage$t_mid, jit[i]+1, lineage$t_end, lwd=2.5, col=rgb(0,0,0,.7), lty=ltype)
	}else{
	segments(jit[i]+.5, lineage$t_mid, jit[i]+1, lineage$t_end, lwd=1, col=rgb(0,0,0,.2), lty=ltype)
	}
	points (jit[i], lineage$t_ini, cex=2.5, pch=pchini, col=colini)
	points (jit[i]+.5, lineage$t_mid, cex=2.5, pch=pchmid, col=colmid)
	points (jit[i]+1, lineage$t_end, cex=2.5, pch=pchend, col=colend)
}

####SpBxScer
plot.new()
plot.window(xlim=c(-.05,1.05), ylim=c(0,1))
axis(1, at=c(0, .5, 1), lwd=2.5, cex.axis=2.5, labels=NA)
axis(1, labels=c(expression("T"["ini"]), expression("T"["mid"]), expression("T"["end"])), at=c(0, .5, 1), cex.axis=3, line=1, lwd=0)
#axis(2, lwd=2.5, las=1, cex.axis=3)			
legend(0.4,1.05 , c("H1", "H2"), col=c(rgb(0,0,0,0.4), rgb(0,0,0,0.4)), lty=c(1,2), lwd=2, bty="n", cex=2)
#legend(0,.985 , c("Lost"), col= rgb(0,0,0,.4), lty=2, lwd=2, bty="n", cex=2)
legend(0,1.05, c("2n", "3n", "4n"), pch=c(1, 16, 16), col=c(grey,orange,purple), cex=2, bty="n", pt.cex=2.5)
par(xpd=NA)
#polygon(c(seq(-0.05,1.05,.1), seq(1.05,-0.05,-.1)), c(rep(1.04,12), rep(1.1,12)), col=rgb(217/255,217/255,217/255), lty=0)
#title(main="H", ylab="Fertility", cex.main=4, cex.lab=2)
polygon(c(seq(-0.05,1.05,.1), seq(1.05,-0.05,-.1)), c(rep(1.04,12), rep(1.2,12)), col=black, lty=0)
title(main=expression("H"["div"]), cex.main=4, cex.lab=2, col.main="White")
par(xpd=0)


for(i in 1:nrow(H)){

	jit = jitter(rep(0, nrow(H)))
	lineage = H[i,]
	if(grepl("C",lineage$Strain) == TRUE){
		ltype = 1
	} else{
		ltype = 2
}
	#if(lineage$Last_glycerol == "P35"){
		#ltype = 1
	#}else{
		#ltype = 2		
#}

	if(lineage$Ploidy_tini == 2){
		pchini= 1
		colini = grey
	} else if (lineage$Ploidy_tini == 3){
		pchini= 16
		colini = orange
	} else if (lineage$Ploidy_tini == 4){
		pchini= 16
		colini = purple
}
	if(lineage$Ploidy_tmid == 2){
		pchmid= 1
		colmid = grey
	} else if (lineage$Ploidy_tmid == 3){
		pchmid= 16
		colmid = orange
	} else if (lineage$Ploidy_tmid == 4){
		pchmid= 16
		colmid = purple

}
	if(lineage$Ploidy_tend == 2){
		pchend= 1
		colend = grey
	} else if (lineage$Ploidy_tend == 3){
		pchend= 16
		colend = orange
	} else if (lineage$Ploidy_tend == 4){
		pchend= 16
		colend = purple

}
	strain = which(row.names(pvalue_fisher_corr) == H[i,1])
	if(pvalue_fisher_corr[strain, 6] <= .05 & is.na(pvalue_fisher_corr[strain, 6]) == "FALSE" ){
	segments(jit[i],lineage$t_ini, jit[i]+0.5, lineage$t_mid, lwd=2.5, col=rgb(0,0,0,0.7), lty=ltype)
	}else{
	segments(jit[i],lineage$t_ini, jit[i]+0.5, lineage$t_mid, lwd=1, col=rgb(0,0,0,.2), lty=ltype)
	}
	if(pvalue_fisher_corr[strain, 6] <= .05 & is.na(pvalue_fisher_corr[strain, 6]) == "FALSE" ){
	segments(jit[i]+.5, lineage$t_mid, jit[i]+1, lineage$t_end, lwd=2.5, col=rgb(0,0,0,.7), lty=ltype)
	}else{
	segments(jit[i]+.5, lineage$t_mid, jit[i]+1, lineage$t_end, lwd=1, col=rgb(0,0,0,.2), lty=ltype)
	}
	points (jit[i], lineage$t_ini, cex=2.5, pch=pchini, col=colini)
	points (jit[i]+.5, lineage$t_mid, cex=2.5, pch=pchmid, col=colmid)
	points (jit[i]+1, lineage$t_end, cex=2.5, pch=pchend, col=colend)
}

dev.off()


#################################
######Supplemetary figures#######
#################################

#### Fig S1 - Generations per passage -
#subset only considereing good YPD
gen_good = subset(generation, generation$YPD == "good")

#T-test for comparing the timepoints
Ti= subset(gen_good$generation_number, gen_good$Time=="ini")
Te= subset(gen_good$generation_number, gen_good$Time=="end")

t.test(Ti,Te)

#Welch Two Sample t-test

#data:  Ti and Te
#t = 0.12353, df = 60.886, p-value = 0.9021
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
# -0.3351109  0.3792380
#sample estimates:
#mean of x mean of y 
# 22.13010  22.10803 

#ANOVA for the effect of ploidy on generation number
Ploidy_effect=aov(gen_good$generation_number~gen_good$Ploidy)
#TukeyHSD on the fit
PloidyTHSD = TukeyHSD(Ploidy_effect)

#Tukey multiple comparisons of means
   #95% family-wise confidence level

#Fit: aov(formula = gen_good$generation_number ~ gen_good$Ploidy)

#$`gen_good$Ploidy`
#              diff       lwr        upr     p adj
#3n-2n -0.722307462 -1.186441 -0.2581744 0.0011873
#4n-2n -0.714160834 -1.328151 -0.1001705 0.0187932
#4n-3n  0.008146628 -0.688053  0.7043463 0.9995641


pdf("Fig S1.pdf", height=6, width=8)

par(mar=c(5,5,4,4), mfrow=c(1,2))
boxplot(gen_good$generation_number~gen_good$Time, axes=F, lwd=2.5, lty=1, ylim=c(20,24))
title(ylab= "Generations per passage", cex.lab=1.5)
axis(1, lwd=0, at=c(1,2), labels=c(expression("T"["ini"]), expression("T"["end"])), cex.axis=1.5)
axis(2, lwd=2.5, las=1, cex.axis=1.5)
mtext("A", side=3, line=1.2, adj=-.45, cex=2)
segments(1,23.8,2,23.8, lwd=2.5)
segments(c(1,2),c(23.85,23.85),c(1,2),c(23.75,23.75), lwd=2.5)
text(1.5, 23.9, "n.s.", cex=1.5)


boxplot(gen_good$generation_number~gen_good$Ploidy, axes=F, lwd=2.5, lty=1, ylim=c(20,24))
title(ylab= "Generations per passage", xlab="Ploidy", cex.lab=1.5)
axis(1, lwd=0, at=c(1,2,3), labels=c("2n", "3n", "4n"), cex.axis=1.5)
axis(2, lwd=2.5, las=1, cex.axis=1.5)
mtext("B", side=3, line=1.2, adj=-.45, cex=2)
segments(c(1,1,2),c(23.8,24,23),c(2,3,3),c(23.8,24,23), lwd=2.5)
segments(c(1,2,1,3,2,3),c(23.85,23.85,24.05,24.05,23.05,23.05),c(1,2,1,3,2,3),c(23.75,23.75, 23.95, 23.95, 22.95, 22.95), lwd=2.5)
text(c(1.5,2,2.5), c(23.9,24.1,23.1), c("**","*","n.s."), cex=1.5)

dev.off()


#### Fig. S2 - Survival curves with different groupings -


pdf("Fig S2.pdf", height=5, width=13.33)
par(mfrow=c(1,3), mar=c(5,6,4,4), lty=1)

##Per cross
plot.new()
plot.window(xlim=c(0,35), ylim=c(.5,1))
title(xlab="Passage", ylab="Proportion of surviving lines (%)", cex.lab=2)
axis(1, cex.axis=1.5, lwd=2.5)
axis(2, las=1, cex.axis=1.5, lwd=2.5)
lines(fit_cross, lwd=2, col=c(blue, green, black, blue, green, black, red, red), lty=c(1,1,1,2,2,2,1,2))
#lines(fit_cross, lwd=2, col=c(rgb(0,0,(139/255),.01), rgb(0,(139/255),0,.01), rgb(0,0,0,.01), rgb(0,0,(139/255),.01), rgb(0,(139/255),0,.01) ,rgb(0,0,0,.01), rgb((237/255),0,0,.01),rgb((237/255),0,0,.01)) , lty=c(1,1,1,1), conf.times=seq(0,36,.01), conf.cap=0)
abline(v=c(0,1,seq(4,31,3), 35), lwd=6, col=rgb(160/255,160/255,160/255,0.2), lty=1)
par(lty=0)
legend(2,.66, title="Crosses", ncol=2L, c("VL1 (48)","VL2 (48)" , "L1 (96)","L2 (96)", "M1 (96)","M2 (96)", "H1 (96)", "H2 (96)"), lty=c(2,1,1,2,1,2,1,2),lwd=2, col=c(red, red, blue, blue, green, green, black, black), bg="white", cex=1.5)
#legend(15,.65, c("Set 1", "Set 2"),lty=c(1,2) , bg="white", lwd=2, cex=1.5)
mtext("A", side=3, line=0.5, adj=-.15, cex=3)


##Per replication block

par(lty=1)
plot.new()
plot.window(xlim=c(0,35), ylim=c(.5,1))
title(xlab="Passage", ylab="Proportion of surviving lines (%)", cex.lab=2)
axis(1, cex.axis=1.5, lwd=2.5)
axis(2, las=1, cex.axis=1.5, lwd=2.5)
lines(fit_block, lwd=2, col=c(black, black, rgb(160/255,160/255,160/255,0.7), rgb(160/255,160/255,160/255,0.7), red,red), lty=c(1,2,1,2,1,2))
#lines(fit_block, lwd=2, col=c(rgb(0,0,0,.01), rgb(0,0,0,.01), rgb(160/255,160/255,160/255,.01), rgb(160/255,160/255,160/255,.01), rgb((237/255),0,0,.01),rgb((237/255),0,0,.01)) , lty=c(1,1,1,1), conf.times=seq(0,36,.01), conf.cap=0)
abline(v=c(0,1,seq(4,31,3), 35), lwd=6, col=rgb(160/255,160/255,160/255,0.2))
par(lty=0)
legend(2,.7, title="Replication block",c("L1-M1-H1 (1-64)", "L1-M1-H1 (65-96)", "L2-M2-H2 (1-64)", "L2-M2-H2 (65-96)", "VL1 (1-48)", "VL2 (1-48)"), lty=c(1,2,1,2,2,1), col=c("black", "black", rgb(160/255,160/255,160/255,0.7), rgb(160/255,160/255,160/255,0.7), red,red), lwd=2, bg="white", cex=1.3)
mtext("B", side=3, line=0.5, adj=-.15, cex=3)


##Per group (L-H, VL-M)

par(lty=1)
plot.new()
plot.window(xlim=c(0,35), ylim=c(.5,1))
title(xlab="Passage", ylab="Proportion of surviving lines (%)", cex.lab=2)
axis(1, cex.axis=1.5, lwd=2.5)
axis(2, las=1, cex.axis=1.5, lwd=2.5)
lines(fit_group, lwd=2, col=c(rgb(160/255,160/255,160/255), black), lty=c(1,1))
lines(fit_group, lwd=2, col=c(rgb(160/255,160/255,160/255,.01), rgb(0,0,0,.01)), conf.times=seq(0,36,.01), conf.cap=0)
abline(v=c(0,1,seq(4,31,3), 35), lwd=6, col=rgb(160/255,160/255,160/255,0.2), lty=1)
par(lty=0)
legend(2,.7, c("VL-M (288)", "L-H (384)"), fill=c("black", "grey"), bg="white", cex=1.5)
mtext("C", side=3, line=0.5, adj=-.15, cex=3)
par(lty=1)

dev.off()

####Fig. S3 - spore viability per group -

Div_fert=aov(Fertility~Divergence, data=fertility)
TukeyHSD(Div_fert)[1]

#H2-H1    0.037485256 -0.043678652  0.11864916 8.548102e-01
#L2-L1    0.050946168 -0.007759479  0.10965182 1.438471e-01
#M2-M1   -0.166275362 -0.241923651 -0.09062707 1.983513e-09
#VL2-VL1  0.265856968  0.195716493  0.33599744 4.255382e-10


pdf("Fig S3.pdf", height=7, width=6)
par(mar=c(5,5,4,4))
boxplot(VL1$t_ini, L$t_ini, M1$t_ini, H$t_ini, axes=F, lwd=2, lty=1, ylim=c(0,1))
boxplot(VL2$t_ini, L$t_ini, M2$t_ini, H$t_ini, axes=F, lwd=2, lty=1, border=c("grey","black","grey","black") , add=T)
title(xlab= "Crosses", ylab="Spore viability", cex.lab=1.5)
axis(1, lwd=2, labels= c("VL", "L", "M", "H"), at=c(1:4), cex.axis=1.5)
axis(2, lwd=2, las=2, cex.axis=1.5)
abline(h=c(.385, .438, .552), lwd=2, lty=2, col=c(blue, green, red))
text(rep(4.4,3),c(.574, .46, .407), c("VL", "M", "L"), col=c(red, green, blue), cex=1.2)
dev.off()

####Fig. S4



####Fig S5 - ITC slope graph and tetrad types -

pdf("Fig S5.pdf", height=4, width=12)

par(mfrow=c(1,3), lty=1, mar=c(5,5,4,4))

plot.new()
plot.window(xlim=c(0,12), ylim=c(0,1))
axis(1, labels=c("ITC 0", "ITC 1", "ITC 6", "ITC 12"),  at=c(0,1,6,12), lwd=2.5, cex.axis=1.4, las=2)
axis(2, lwd=2.5, las=1, cex.axis=1.4)
legend("topleft" , c("L1", "L2"), col=c(rgb(1,0,0,0.4), rgb(0,0,0,0.4)), lty=1, lwd=2, bty="n", cex=1.2)
par(xpd=NA)
polygon(c(seq(-.5,12.5,1), seq(12.5,-.5,-1)), c(rep(1.06,14), rep(1.23,14)), col=black, lty=0)
title(main="ITC", ylab="Fertility", col.main="white", cex.main=4, cex.lab=2)
mtext("A", side=3, line=1.2, adj=-.15, cex=2)

par(xpd=0)


for(i in 1:nrow(ITC_data)){

	lineage= ITC_data[i,]
	if(grepl("A",lineage$lineage) == TRUE){
		color= rgb(1,0,0,0.4)
	} else{
		color = rgb(0,0,0,0.4)
}

	segments(0,lineage$P0, 1,lineage$P1, col=color, lwd=2)
	segments(1,lineage$P1, 6,lineage$P6, col=color, lwd=2)
	segments(6,lineage$P6, 12,lineage$P12, col=color, lwd=2)
	points(c(0,1,6,12), c(lineage$P0,lineage$P1,lineage$P6,lineage$P12), pch=1)
}	



#sporecounts ITC

ITC_spores=aggregate(ITC_table[,5:9], by=list(ITC_table$ITC_number, ITC_table$Divergence), FUN=sum)

ITC_L1 = t(ITC_spores[c(1:4),c(3:7)])
colnames(ITC_L1) = c("ITC_0", "ITC_1", "ITC_6", "ITC_12")
ITC_L1= apply(ITC_L1, 2, FUN=prop.table)

ITC_L2 = t(ITC_spores[c(5:8),c(3:7)])
colnames(ITC_L2) = c("ITC_0", "ITC_1", "ITC_6", "ITC_12")
ITC_L2= apply(ITC_L2, 2, FUN=prop.table)

par(lty=0)

barplot(ITC_L1, col=c("red", "orange", "yellow", "green", "darkgreen"), cex.axis=1.5, cex.names=1.5, lwd=2.5, las=1, names.arg=c("ITC 0", "ITC 1","ITC 6","ITC 12"))
par(xpd=NA)
polygon(c(seq(0.1,4.9,.1), seq(4.9,0.1,-.1)), c(rep(1.02,49), rep(1.18,49)), col=black, lty=0)
title(main="ITC L1", ylab="Tetrad proportion (%)", col.main="white", cex.main=4, cex.lab=2)
legend(1,-.10, c("0 spores", "1 spore", "2 spores", "3 spores", "4 spores"), fill=c("red", "orange", "yellow", "green", "darkgreen"), bg="white", box.col="white", ncol=5, cex=1.5)
mtext("B", side=3, line=1.2, adj=-.15, cex=2)
par(xpd=0)


barplot(ITC_L2, col=c("red", "orange", "yellow", "green", "darkgreen"), cex.axis=1.5, cex.names=1.5, las=1, lwd=2.5, names.arg=c("ITC 0", "ITC 1","ITC 6","ITC 12"))
par(xpd=NA)
polygon(c(seq(0.1,4.9,.1), seq(4.9,0.1,-.1)), c(rep(1.02,49), rep(1.18,49)), col=black, lty=0)
title(main="ITC L2",ylab="Tetrad proportion (%)", col.main="white", cex.main=4, cex.lab=2)
mtext("C", side=3, line=1.2, adj=-.15, cex=2)
par(xpd=0)


dev.off()


#### Fig S6 - Spore counts evolution lines -
evo_spo = subset(fertility, fertility$Fertility != "NA" & fertility$Sporulation == "yes")

evo_spores=aggregate(evo_spo[,9:13], by=list(evo_spo$Time_point, evo_spo$Cross_type), FUN=sum)

tetrad_VL = t(evo_spores[c(4:6),c(3:7)])
colnames(tetrad_VL) = c("Tini", "Tmid", "Tend")
tetrad_VL= apply(tetrad_VL, 2, FUN=prop.table)

tetrad_L = t(evo_spores[c(7:9),c(3:7)])
colnames(tetrad_L) = c("Tini", "Tmid", "Tend")
tetrad_L= apply(tetrad_L, 2, FUN=prop.table)

tetrad_M = t(evo_spores[c(1:3),c(3:7)])
colnames(tetrad_M) = c("Tini", "Tmid", "Tend")
tetrad_M= apply(tetrad_M, 2, FUN=prop.table)

tetrad_H = t(evo_spores[c(10:12),c(3:7)])
colnames(tetrad_H) = c("Tini", "Tmid", "Tend")
tetrad_H= apply(tetrad_H, 2, FUN=prop.table)

#####plots

pdf("Fig S6.pdf", height= 10, width=10)
par(mar=c(5,5,4,4), lty=0, mfrow=c(2,2))


barplot(tetrad_VL, col=c("red", "orange", "yellow", "green", "darkgreen"), cex.axis=1.5, cex.names=2, las=1, lwd=2.5, names.arg=c(expression("T"["ini"]), expression("T"["mid"]), expression("T"["end"])))
par(xpd=NA)
polygon(c(seq(0.1,3.7,.1), seq(3.7,0.1,-.1)), c(rep(1.02,37), rep(1.18,37)), col=red, lty=0)
title(main="VL", ylab="Tetrad proportion (%)", col.main="white", cex.main=4, cex.lab=1.5)
legend(1,-.10, c("0 spores", "1 spore", "2 spores", "3 spores", "4 spores"), fill=c("red", "orange", "yellow", "green", "darkgreen"), bg="white", box.col="white", ncol=5, cex=1.4)
mtext("A", side=3, line=1.2, adj=-.15, cex=3)
par(xpd=0)

barplot(tetrad_L, col=c("red", "orange", "yellow", "green", "darkgreen"), cex.axis=1.5, cex.names=2, las=1, lwd=2.5, names.arg=c(expression("T"["ini"]), expression("T"["mid"]), expression("T"["end"])))
par(xpd=NA)
polygon(c(seq(0.1,3.7,.1), seq(3.7,0.1,-.1)), c(rep(1.02,37), rep(1.18,37)), col=blue, lty=0)
title(main="L", ylab="Tetrad proportion (%)", col.main="white", cex.main=4, cex.lab=1.5)
mtext("B", side=3, line=1.2, adj=-.15, cex=3)
par(xpd=0)

barplot(tetrad_M, col=c("red", "orange", "yellow", "green", "darkgreen"), cex.axis=1.5, cex.names=2, las=1, lwd=2.5, names.arg=c(expression("T"["ini"]), expression("T"["mid"]), expression("T"["end"])))
par(xpd=NA)
polygon(c(seq(0.1,3.7,.1), seq(3.7,0.1,-.1)), c(rep(1.02,37), rep(1.18,37)), col=green, lty=0)
title(main="M", ylab="Tetrad proportion (%)", col.main="white", cex.main=4, cex.lab=1.5)
mtext("C", side=3, line=1.2, adj=-.15, cex=3)
par(xpd=0)

barplot(tetrad_H, col=c("red", "orange", "yellow", "green", "darkgreen"), cex.axis=1.5, cex.names=2, las=1, lwd=2.5, names.arg=c(expression("T"["ini"]), expression("T"["mid"]), expression("T"["end"])))
par(xpd=NA)
polygon(c(seq(0.1,3.7,.1), seq(3.7,0.1,-.1)), c(rep(1.02,37), rep(1.18,37)), col=black, lty=0)
title(main="H", ylab="Tetrad proportion (%)", col.main="white", cex.main=4, cex.lab=1.5)
mtext("D", side=3, line=1.2, adj=-.15, cex=3)
par(xpd=0)


dev.off()

###Fig S7 - Fertility of diploid vs triploids -



#Data for the timepoints
dip_ini = L$t_ini[c(which(L$Ploidy_tini == 2))]
dip_mid =L$t_mid[c(which(L$Ploidy_tmid == 2))]
dip_end =L$t_end[c(which(L$Ploidy_tend == 2))]

trip_ini = L$t_ini[c(which(L$Ploidy_tini == 3))]
trip_mid =L$t_mid[c(which(L$Ploidy_tmid == 3))]
trip_end =L$t_end[c(which(L$Ploidy_tend == 3))]

#Beanplot Figure
pdf("Fig S7.pdf", height=5, width=5)
par(mar=c(5,5,4,4))
beanplot(dip_ini,trip_ini, dip_mid, trip_mid, dip_end, trip_end, side="both", overallline="median", col=list(c("black", "grey", "white", "grey"), c("grey","grey" ,"grey" ,"black" )), ylim=c(0,1), names=c("Tini", "Tmid", "Tend"), what=c(0,1,1,0), axes=F, border=NA)
axis(1, at=c(1,2,3), labels=c(expression("T"["ini"]), expression("T"["mid"]), expression("T"["end"])), cex.axis=1.4, lwd=2.5)
axis(2, las=1, lwd=2.5, cex.axis=1.4)
title(ylab="Fertility", cex.lab=1.4)
par(lty=0)
legend(1.5,1, c("Diploids", "Triploids"), fill=c(black,"grey"), bty="n")
par(lty=1)
segments(1,0,1,.7, lwd=2)
segments(2,0,2,.7, lwd=2)
segments(3,0,3,.7, lwd=2)

dev.off()

####Fig S8



####Fig S9



####Fig S10




##### Supplementary figure 22 - Figure effect of yeast extract -
genYPD = generation[c(47:80),]

####Plot image of YPD batches as panel A
pdf("Fig S12.pdf", height=4, width=12.489)

#plan the layout of the figure
fig = matrix(c(rep(1,6),rep(2,4)), byrow=F, nrow=2, ncol=5)
layout(fig)
img2 = readJPEG("Figure_Sup_ade2_coloration.jpg")

par(mar=c(4,0,4,4))
plot.new()
plot.window(xlim=c(0,2), ylim=c(0,1.0675))
rasterImage(img2, 0,0,2,1.0675)
mtext("A", side=3, line=1.2, adj=0.01, cex=2)

par(mar=c(4,5,4,4))
boxplot(genYPD$generation_number~genYPD$YPD, axes=F, lwd=2.5, lty=1, ylim=c(18,23))
title(ylab= "Generations per passage", cex.lab=1.5)
axis(1, lwd=0, at=c(1,2), labels=c("Lot# VM806753 728", "Lot# VM781753 704"), cex.axis=1.5)
axis(2, lwd=2.5, las=1, cex.axis=1.5)
mtext("B", side=3, line=1.2, adj=-.1, cex=2)

dev.off()


