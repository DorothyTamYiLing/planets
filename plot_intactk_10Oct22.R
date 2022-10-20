install.packages("shape")
library("shape")

setwd("")
mytable<-read.table("")

plot(NA,NA, xlim=c(0,10), ylim=c(0,10))
Arrows(1,7,1.25,7, arr.type="simple",arr.length=1.25)
#the 4 numbers: arrow start at x-axis, arrow start at y-axis, arrow end at x-axis, arrow end at y-axis

#install.packages("ggplot2")
library("ggplot2")
#install.packages("ggforce")
library("ggforce")
#install.packages("ggpubr")
library("ggpubr")

arrows <- data.frame(
  stringsAsFactors = FALSE,
  x1 = rep(c(0), 4), 
  y1 = rep(c(0), 4), 
  x2 = c(0.03547301, 0.10231616, 0.10231616, 0.08409951
  ), 
  y2 = c(-0.122964349, 0.467737581, 0.009871989, 0.005768529),
  mycat=1:4
)

p <- ggplot(data = arrows, aes(x = x1, y = y1, xend = x2, yend = y2, color=mycat)) + 
  geom_link(arrow = grid::arrow(length = grid::unit(0.2, 'cm')))+
  scale_colour_gradientn(name = "mycat", 
                        colours=c("darkred","orange","yellow","green"))

p <- ggplot(data = arrows, aes(x = x1, y = y1, xend = x2, yend = y2, color=mycat)) + 
  geom_link(arrow = grid::arrow(length = grid::unit(0.2, 'cm')))+
  scale_colour_gradient(name = "mycat", 
                        low = "blue", high = "red")

#for clus1 clus2 GWAS only
#install.packages("ggplot2")
library("ggplot2")
#install.packages("ggforce")
library("ggforce")
#install.packages("ggpubr")
library("ggpubr")

setwd("/Users/dorothytam/Desktop/")
mytable<-read.delim("myintactk_out.txt",header=T,sep="\t")
#select kmers that are fwd in "0" genomes and rev in "1" genomes
myx<-mytable[which(mytable$fwdk_0gen_prop>0.6 & mytable$revk_1gen_prop>0.6),]
myx<-myx[order(myx$fwdk_0gen_med),]
#remove rows that are duplicated in myx$fwdk_0gen_med
mykeeprow<-c()
mycoorlist<-unique(myx$fwdk_0gen_med)
for (i in 1:length(mycoorlist)){
  mykeeprow<-c(mykeeprow,which(myx$fwdk_0gen_med==mycoorlist[i])[1])
}
myx_unqiue<-myx[mykeeprow,]
#indexing the kmers for colour gradients
myx_unqiue$mycol_index<-1:nrow(myx_unqiue)

myx_unqiue$x<-myx_unqiue$fwdk_0gen_med+15
myx_unqiue$y<-myx_unqiue$revk_1gen_med-15

p <- ggplot(data = myx_unqiue) + 
  geom_link(aes(x = fwdk_0gen_med, y = 0, xend = x, yend = 0, color=mycol_index),arrow = grid::arrow(length = grid::unit(0.6, 'cm')))+
  geom_link(aes(x = revk_1gen_med, y = 1, xend = y, yend = 1, color=mycol_index),arrow = grid::arrow(length = grid::unit(0.6, 'cm')))+
  scale_colour_gradientn(name = "mycol_index", 
                         colours=c("darkred","orange","red","blue","chartreuse3"))+
  scale_x_continuous(breaks = seq(1, 4300000, by = 100000))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#for PRN GWAS
setwd("/Users/dorothytam/Desktop/")
mytable<-read.table("PRN_468_tree_merge200GWAS_exbadchisq_myintactk_out.txt",header=T)
colnames(mytable)

#plotting all the intact kmers here

mytable<-mytable[order(mytable$fwdk_0gen_med),]

#indexing the kmers for colour gradients
mytable$mycol_index<-1:nrow(mytable)

#select kmers that are "fwd" in majority of the "0" genomes
my0<-mytable[which(mytable$fwdk_0gen_prop>0.6),]
colnames(my0)<-c("kmer","fwdk_med","revk_med")
  
#remove rows that are duplicated in myx$fwdk_0gen_med
mykeeprow<-c()
mycoorlist<-unique(my0$fwdk_med)
for (i in 1:length(mycoorlist)){
  mykeeprow<-c(mykeeprow,which(my0$fwdk_med==mycoorlist[i])[1])
}
my0_unqiue<-my0[mykeeprow,]

#select "1" genomes colummns
my1<-mytable[,c("kmer","fwdk_1gen_med","revk_1gen_med")]
colnames(my1)<-c("kmer","fwdk_med","revk_med")

#remove rows that are duplicated in myx$fwdk_0gen_med
mykeeprow<-c()
mycoorlist<-unique(my1$fwdk_med)
for (i in 1:length(mycoorlist)){
  mykeeprow<-c(mykeeprow,which(my1$fwdk_med==mycoorlist[i])[1])
}
my1_unqiue<-my1[mykeeprow,]

my0_unqiue$fwdk_med_end<-my0_unqiue$fwdk_med+15
my0_unqiue$revk_med_end<-my0_unqiue$revk_med+15
my1_unqiue$fwdk_med_end<-my1_unqiue$fwdk_med+15
my1_unqiue$revk_med_end<-my1_unqiue$revk_med+15

my0_unqiue$gen<-"0"
my1_unqiue$gen<-"1"

#indexing the kmers for colour gradients,one kmer one colour
myplotall<-rbind(my0_unqiue,my1_unqiue)
myplotall$mycol_index<-"NA"
mylistcol<-unique(myplotall$kmer)
for (i in 1:length(mylistcol)){
  myplotall[which(myplotall$kmer==mylistcol[i]),"mycol_index"]<-as.numeric(i)
}

#install.packages("ggplot2")
library("ggplot2")
#install.packages("ggforce")
library("ggforce")
#install.packages("ggpubr")
library("ggpubr")

myplotall$mycol_index<-as.numeric(myplotall$mycol_index)

p <- ggplot(data = myplotall) + 
  geom_link(aes(x = fwdk_med, y = gen, xend = fwdk_med_end, yend = gen, color=mycol_index),arrow = grid::arrow(length = grid::unit(0.9, 'cm'),size = 5))+
  geom_link(aes(x = revk_med_end, y = gen, xend = revk_med, yend = gen, color=mycol_index),arrow = grid::arrow(length = grid::unit(0.9, 'cm')))+
  scale_colour_gradientn(name = "mycol_index", 
                         colours=c("orange","red","blue","chartreuse3"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_x_continuous(breaks = seq(1, 4300000, by = 100000))


