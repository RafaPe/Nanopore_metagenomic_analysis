setwd("../") #Set working directory

knitr::opts_chunk$set(echo = TRUE)
library(readr)
BC01_16S <- read_delim("aligned2/BC01_filt.fasta_align.paf", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
BC02_16S <- read_delim("aligned2/BC02_filt.fasta_align.paf", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
BC03_16S <- read_delim("aligned2/BC03_filt.fasta_align.paf", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
BC04_16S <- read_delim("aligned2/BC04_filt.fasta_align.paf", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

dfz1<-BC01_16S[,1:17]
dfz2<-BC02_16S[,1:17]
dfz3<-BC03_16S[,1:17]
dfz4<-BC04_16S[,1:17]

colnames(dfz1)<-c("Query","Q_length","Q_start","Q_end","Strand","op","T_length","T_start","T_end","N_res_matches","Align_block","MapQ","NM","ms","AS","nn","P_S")
colnames(dfz2)<-c("Query","Q_length","Q_start","Q_end","Strand","op","T_length","T_start","T_end","N_res_matches","Align_block","MapQ","NM","ms","AS","nn","P_S")
colnames(dfz3)<-c("Query","Q_length","Q_start","Q_end","Strand","op","T_length","T_start","T_end","N_res_matches","Align_block","MapQ","NM","ms","AS","nn","P_S")
colnames(dfz4)<-c("Query","Q_length","Q_start","Q_end","Strand","op","T_length","T_start","T_end","N_res_matches","Align_block","MapQ","NM","ms","AS","nn","P_S")


dfz1$AS <- as.numeric(gsub("AS:i:", "", dfz1$AS))
dfz2$AS <- as.numeric(gsub("AS:i:", "", dfz2$AS))
dfz3$AS <- as.numeric(gsub("AS:i:", "", dfz3$AS))
dfz4$AS <- as.numeric(gsub("AS:i:", "", dfz4$AS))
#-----------------------------------------------------
df_z12<-cbind.data.frame(dfz1$Query,dfz1$op,dfz1$N_res_matches,dfz1$Align_block,dfz1$AS,dfz1$MapQ)
df_z22<-cbind.data.frame(dfz2$Query,dfz2$op,dfz2$N_res_matches,dfz2$Align_block,dfz2$AS,dfz2$MapQ)
df_z32<-cbind.data.frame(dfz3$Query,dfz3$op,dfz3$N_res_matches,dfz3$Align_block,dfz3$AS,dfz3$MapQ)
df_z42<-cbind.data.frame(dfz4$Query,dfz4$op,dfz4$N_res_matches,dfz4$Align_block,dfz4$AS,dfz4$MapQ)
#-----------------------------------------------------
colnames(df_z12)<-c("Query","op","N_res_matches","Align_block","AS","MapQ")
colnames(df_z22)<-c("Query","op","N_res_matches","Align_block","AS","MapQ")
colnames(df_z32)<-c("Query","op","N_res_matches","Align_block","AS","MapQ")
colnames(df_z42)<-c("Query","op","N_res_matches","Align_block","AS","MapQ")
#-----------------------------------------------------
unique(df_z12$Query)
unique(df_z22$Query)
unique(df_z32$Query)
unique(df_z42$Query)
#-----------------------------------------------------
per.match<-(df_z12$N_res_matches/df_z12$Align_block)*100
df_2_modz1<-cbind.data.frame(df_z12,per.match)
per.match<-(df_z22$N_res_matches/df_z22$Align_block)*100
df_2_modz2<-cbind.data.frame(df_z22,per.match)
per.match<-(df_z32$N_res_matches/df_z32$Align_block)*100
df_2_modz3<-cbind.data.frame(df_z32,per.match)
per.match<-(df_z42$N_res_matches/df_z42$Align_block)*100
df_2_modz4<-cbind.data.frame(df_z42,per.match)
#-----------------------------------------------------
tax<-read.csv("files/rrn_tax.txt",sep = "\t")
#buscarv_tax<-merge(df_2,Tax, by="op") ##revision
buscarv_taxz1<-merge(df_2_modz1,tax, by="op")
buscarv_taxz2<-merge(df_2_modz2,tax, by="op")
buscarv_taxz3<-merge(df_2_modz3,tax, by="op")
buscarv_taxz4<-merge(df_2_modz4,tax, by="op")
#-----------------------------------------------------
colnames(buscarv_taxz1)<-c("op","Query","N_res_matches","Align_block","AS","MapQ", "Matching","Tax")
colnames(buscarv_taxz2)<-c("op","Query","N_res_matches","Align_block","AS","MapQ", "Matching","Tax")
colnames(buscarv_taxz3)<-c("op","Query","N_res_matches","Align_block","AS","MapQ", "Matching","Tax")
colnames(buscarv_taxz4)<-c("op","Query","N_res_matches","Align_block","AS","MapQ", "Matching","Tax")
#-----------------------------------------------------
##For the 16S rRNA gene
df_3z1<-subset(buscarv_taxz1, Align_block > 999)
df_3z2<-subset(buscarv_taxz2, Align_block > 999)
df_3z3<-subset(buscarv_taxz3, Align_block > 999)
df_3z4<-subset(buscarv_taxz4, Align_block > 999)

library(doBy)
df_4z1<- summaryBy(AS ~ Query + Tax + Matching + MapQ, data=df_3z1, FUN=max)
df_4z2<- summaryBy(AS ~ Query + Tax + Matching + MapQ, data=df_3z2, FUN=max)
df_4z3<- summaryBy(AS ~ Query + Tax + Matching + MapQ, data=df_3z3, FUN=max)
df_4z4<- summaryBy(AS ~ Query + Tax + Matching + MapQ, data=df_3z4, FUN=max)
#-----------------------------------------------------
df_5z1<-df_4z1[!(duplicated(df_4z1$Query) | duplicated(df_4z1$Query, fromLast = TRUE)), ]  # unique tax records
df_5z2<-df_4z2[!(duplicated(df_4z2$Query) | duplicated(df_4z2$Query, fromLast = TRUE)), ]  # unique tax records
df_5z3<-df_4z3[!(duplicated(df_4z3$Query) | duplicated(df_4z3$Query, fromLast = TRUE)), ]  # unique tax records
df_5z4<-df_4z4[!(duplicated(df_4z4$Query) | duplicated(df_4z4$Query, fromLast = TRUE)), ]  # unique tax records
#-----------------------------------------------------
df_dupz1<-df_4z1[(duplicated(df_4z1$Query) | duplicated(df_4z1$Query, fromLast = TRUE)), ] # replicates
df_dupz2<-df_4z2[(duplicated(df_4z2$Query) | duplicated(df_4z2$Query, fromLast = TRUE)), ] # replicates
df_dupz3<-df_4z3[(duplicated(df_4z3$Query) | duplicated(df_4z3$Query, fromLast = TRUE)), ] # replicates
df_dupz4<-df_4z4[(duplicated(df_4z4$Query) | duplicated(df_4z4$Query, fromLast = TRUE)), ] # replicates
#-----------------------------------------------------
df_dup2z1<- df_dupz1[order(df_dupz1$Query, -abs(df_dupz1$AS) ), ] # sort by Query & AS value (highest to lowest)
df_dup2z2<- df_dupz2[order(df_dupz2$Query, -abs(df_dupz2$AS) ), ] # sort by Query & AS value (highest to lowest)
df_dup2z3<- df_dupz3[order(df_dupz3$Query, -abs(df_dupz3$AS) ), ] # sort by Query & AS value (highest to lowest)
df_dup2z4<- df_dupz4[order(df_dupz4$Query, -abs(df_dupz4$AS) ), ] # sort by Query & AS value (highest to lowest)
#-----------------------------------------------------
df_6z1<-df_dup2z1[ !duplicated(df_dup2z1$Query), ]    
df_6z2<-df_dup2z2[ !duplicated(df_dup2z2$Query), ] 
df_6z3<-df_dup2z3[ !duplicated(df_dup2z3$Query), ] 
df_6z4<-df_dup2z4[ !duplicated(df_dup2z4$Query), ] 
# take the first row within each Query
#-----------------------------------------------------
df_7z1 <-rbind.data.frame(df_5z1,df_6z1) 
write.csv(df_7z1, file = "files/r_generated/Z1_16S_Qscores.csv")     # Quality scores file
df_7z2 <-rbind.data.frame(df_5z2,df_6z2) 
write.csv(df_7z2, file = "files/r_generated/Z2_16S_Qscores.csv")     # Quality scores file
df_7z3 <-rbind.data.frame(df_5z3,df_6z3) 
write.csv(df_7z3, file = "files/r_generated/Z3_16S_Qscores.csv")     # Quality scores file
df_7z4 <-rbind.data.frame(df_5z4,df_6z4) 
write.csv(df_7z4, file = "files/r_generated/Z4_16S_Qscores.csv")     # Quality scores file

library(reshape2)
Tax_tablez1<-dcast(data=melt(df_7z1$Tax, id.vars="Tax"), value ~ .)
colnames(Tax_tablez1) <- c("Tax","Counts")
rel.abz1 <- Tax_tablez1$Counts/sum(Tax_tablez1$Counts)
readsz1<-sum(Tax_tablez1$Counts)
readsz1
Tax_tablez2<-dcast(data=melt(df_7z2$Tax, id.vars="Tax"), value ~ .)
colnames(Tax_tablez2) <- c("Tax","Counts")
rel.abz2 <- Tax_tablez2$Counts/sum(Tax_tablez2$Counts)
readsz2<-sum(Tax_tablez2$Counts)
readsz2
Tax_tablez3<-dcast(data=melt(df_7z3$Tax, id.vars="Tax"), value ~ .)
colnames(Tax_tablez3) <- c("Tax","Counts")
rel.abz3 <- Tax_tablez3$Counts/sum(Tax_tablez3$Counts)
readsz3<-sum(Tax_tablez3$Counts)
readsz3
Tax_tablez4<-dcast(data=melt(df_7z4$Tax, id.vars="Tax"), value ~ .)
colnames(Tax_tablez4) <- c("Tax","Counts")
rel.abz4 <- Tax_tablez4$Counts/sum(Tax_tablez4$Counts)
readsz4<-sum(Tax_tablez4$Counts)
readsz4
#-----------------------------------------------------
Tax_table_RAz1<-cbind(Tax_tablez1, rel.abz1,"Z1_16S")
colnames(Tax_table_RAz1) <- c("Tax","Counts","RA","SampleID")
write.csv(Tax_table_RAz1, file = "files/r_generated/Z1_16S_Tax_summ.csv")      # Taxa summary file
Tax_table_RAz2<-cbind(Tax_tablez2, rel.abz2,"Z2_16S")
colnames(Tax_table_RAz2) <- c("Tax","Counts","RA","SampleID")
write.csv(Tax_table_RAz2, file = "files/r_generated/Z2_16S_Tax_summ.csv")      # Taxa summary file
Tax_table_RAz3<-cbind(Tax_tablez3, rel.abz3,"Z3_16S")
colnames(Tax_table_RAz3) <- c("Tax","Counts","RA","SampleID")
write.csv(Tax_table_RAz3, file = "files/r_generated/Z3_16S_Tax_summ.csv")      # Taxa summary file
Tax_table_RAz4<-cbind(Tax_tablez4, rel.abz4,"Z4_16S")
colnames(Tax_table_RAz4) <- c("Tax","Counts","RA","SampleID")
write.csv(Tax_table_RAz4, file = "files/r_generated/Z4_16S_Tax_summ.csv")      # Taxa summary file

library(ggplot2)
plot_z1 <- read_delim("files/r_generated/Z1_16S_Tax_summ.csv", ",", escape_double = FALSE, col_names = TRUE, trim_ws = TRUE)
plot1 <- ggplot(plot_z1, aes(x = SampleID, y = Counts, fill = Tax,)) + geom_bar(aes(), stat="identity", position="stack") 
#-----------------------------------------------------
plot_z2 <- read_delim("files/r_generated/Z2_16S_Tax_summ.csv", ",", escape_double = FALSE, col_names = TRUE, trim_ws = TRUE)
plot2 <- ggplot(plot_z2, aes(x = SampleID, y = Counts, fill = Tax,)) + geom_bar(aes(), stat="identity", position="stack")
#-----------------------------------------------------
plot_z3 <- read_delim("files/r_generated/Z3_16S_Tax_summ.csv", ",", escape_double = FALSE, col_names = TRUE, trim_ws = TRUE)
plot3 <- ggplot(plot_z3, aes(x = SampleID, y = Counts, fill = Tax,)) + geom_bar(aes(), stat="identity", position="stack")
#-----------------------------------------------------
plot_z4<- read_delim("files/r_generated/Z4_16S_Tax_summ.csv", ",", escape_double = FALSE, col_names = TRUE, trim_ws = TRUE)
plot4 <- ggplot(plot_z4, aes(x = SampleID, y = Counts, fill = Tax,)) + geom_bar(aes(), stat="identity", position="stack")
#-----------------------------------------------------
#plot_zF<- read_delim("allF.csv", ",", escape_double = FALSE, col_names = TRUE, trim_ws = TRUE)
#plotF <- ggplot(plot_zF, aes(x = SampleID, y = Counts, fill = Tax,)) + geom_bar(aes(), stat="identity", position="stack")
#-------------Z1----------------------------------------
plot1+ scale_fill_distiller(palette = 'Spectral')
plot1+ scale_fill_gradient()
library(viridis)
plot1 + scale_color_viridis(discrete = TRUE)
plot1 + scale_color_gradientn(colours = rainbow(5))
#save plot to working directory
ggsave(plot1,filename="files/r_generated/BC01_rnn.png",height=9,width=12,units="in",dpi=200)
plot2 <- ggplot(plot_z2, aes(x = SampleID, y = Counts, fill = Tax,)) + geom_bar(aes(), stat="identity", position="stack")
#-------------Z2----------------------------------------
plot2+ scale_fill_distiller(palette = 'Spectral')
plot2+ scale_fill_gradient()
library(viridis)
plot2 + scale_color_viridis(discrete = TRUE)
plot2 + scale_color_gradientn(colours = rainbow(5))
#save plot to working directory
ggsave(plot2,filename="files/r_generated/BC02_rnn.png",height=9,width=12,units="in",dpi=200)
#-------------Z3----------------------------------------
plot3+ scale_fill_distiller(palette = 'Spectral')
plot3+ scale_fill_gradient()
library(viridis)
plot3 + scale_color_viridis(discrete = TRUE)
plot3 + scale_color_gradientn(colours = rainbow(5))
#save plot to working directory
ggsave(plot3,filename="files/r_generated/BC03_rnn.png",height=9,width=12,units="in",dpi=200)
#-------------Z4----------------------------------------
plot4+ scale_fill_distiller(palette = 'Spectral')
plot4+ scale_fill_gradient()
library(viridis)
plot4 + scale_color_viridis(discrete = TRUE)
plot4 + scale_color_gradientn(colours = rainbow(5))
#save plot to working directory
ggsave(plot4,filename="files/r_generated/BC04_rnn.png",height=9,width=12,units="in",dpi=200)
#-------------all----------------------------------------
#plotF+ scale_fill_distiller(palette = 'Spectral')
#plotF+ scale_fill_gradient()
#library(viridis)
#plotF + scale_color_viridis(discrete = TRUE)
#plotF + scale_color_gradientn(colours = rainbow(5))
#save plot to working directory
#ggsave(plotF,filename="ZF_SILVA.png",height=9,width=12,units="in",dpi=200)

