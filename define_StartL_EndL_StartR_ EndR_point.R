#defining StartL, EndL, StartR, EndR from blast output

######New Scheme######
#odd rows refer to the first blast hit of the kmer/genome combination
colnames(data_row_odd)<-c("query_o","subject_o","identity_o","alig_len_o","mismatches_o","gap_o","qstart_o","qend_o","sStart_o","sEnd_o","evalue_o","bitscore_o","label_o")
#even rows refer to the second blast hit of the kmer/genome combination
colnames(data_row_even)<-c("query_e","subject_e","identity_e","alig_len_e","mismatches_e","gap_e","qstart_e","qend_e","sStart_e","sEnd_e","evalue_e","bitscore_e","label_e")

mymerge is the merge table with odd and even blast hit rows pasted side to side

#_o refers to columns in the odd rows
#_e refers to columns in the even rows

#qstart_o sStart_o
#qend_o sEnd_o
#qstart_e sStart_e
#qend_e sEnd_e

#StartL EndL StartR EndR

#if even rows > odd rows, qend_o > qstart_o, qend_e > qstart_e
#qstart_o < qend_o < qstart_e < qend_e
mymerge_1_row<-which(mymerge$qstart_o < mymerge$qend_o & mymerge$qend_o < mymerge$qstart_e & mymerge$qstart_e < mymerge$qend_e)
mymerge[mymerge_1_row,"StartL"]<-mymerge[mymerge_1_row,"sStart_o"] 
mymerge[mymerge_1_row,"EndL"]<-mymerge[mymerge_1_row,"sEnd_o"] 
mymerge[mymerge_1_row,"StartR"]<-mymerge[mymerge_1_row,"sStart_e"] 
mymerge[mymerge_1_row,"EndR"]<-mymerge[mymerge_1_row,"sEnd_e"] 

#if even rows > odd rows, qend_o < qstart_o, qend_e < qstart_e
#qend_o < qstart_o < qend_e < qstart_e 
mymerge_2_row<-which(mymerge$qend_o < mymerge$qstart_o & mymerge$qstart_o < mymerge$qend_e & mymerge$qend_e < mymerge$qstart_e)
mymerge[mymerge_2_row,"StartL"]<-mymerge[mymerge_2_row,"sEnd_o"] 
mymerge[mymerge_2_row,"EndL"]<-mymerge[mymerge_2_row,"sStart_o"] 
mymerge[mymerge_2_row,"StartR"]<-mymerge[mymerge_2_row,"sEnd_e"] 
mymerge[mymerge_2_row,"EndR"]<-mymerge[mymerge_2_row,"sStart_e"]

#if even rows > odd rows, qend_o > qstart_o, qend_e < qstart_e
#qstart_o < qend_o < qend_e < qstart_e
mymerge_3_row<-which(mymerge$qstart_o < mymerge$qend_o & mymerge$qend_o < mymerge$qend_e & mymerge$qend_e < mymerge$qstart_e)
mymerge[mymerge_3_row,"StartL"]<-mymerge[mymerge_3_row,"sStart_o"] 
mymerge[mymerge_3_row,"EndL"]<-mymerge[mymerge_3_row,"sEnd_o"] 
mymerge[mymerge_3_row,"StartR"]<-mymerge[mymerge_3_row,"sEnd_e"] 
mymerge[mymerge_3_row,"EndR"]<-mymerge[mymerge_3_row,"sStart_e"]

#if even rows > odd rows, qend_o < qstart_o, qend_e > qstart_e
#qend_o < qstart_o < qstart_e < qend_e 
mymerge_4_row<-which(mymerge$qend_o < mymerge$qstart_o & mymerge$qstart_o < mymerge$qstart_e & mymerge$start_e < mymerge$qend_e)
mymerge[mymerge_4_row,"StartL"]<-mymerge[mymerge_4_row,"sEnd_o"] 
mymerge[mymerge_4_row,"EndL"]<-mymerge[mymerge_4_row,"sStart_o"] 
mymerge[mymerge_4_row,"StartR"]<-mymerge[mymerge_4_row,"sStart_e"] 
mymerge[mymerge_4_row,"EndR"]<-mymerge[mymerge_4_row,"sEnd_e"]

#if even rows < odd rows, qend_o > qstart_o, qend_e > qstart_e
#qstart_e < qend_e < qstart_o < qend_o
mymerge_5_row<-which(mymerge$qstart_e < mymerge$qend_e & mymerge$qend_e < mymerge$qstart_o & mymerge$qstart_o < mymerge$qend_o)
mymerge[mymerge_5_row,"StartL"]<-mymerge[mymerge_5_row,"sStart_e"] 
mymerge[mymerge_5_row,"EndL"]<-mymerge[mymerge_5_row,"sEnd_e"] 
mymerge[mymerge_5_row,"StartR"]<-mymerge[mymerge_5_row,"sStart_o"] 
mymerge[mymerge_5_row,"EndR"]<-mymerge[mymerge_5_row,"sEnd_o"]

#if even rows < odd rows, qend_o < qstart_o, qend_e < qstart_e
#qend_e < qstart_e < qend_o < qstart_o
mymerge_6_row<-which(mymerge$qend_e < mymerge$qstart_e & mymerge$qstart_e < mymerge$qend_o & mymerge$qend_o < mymerge$qstart_o)
mymerge[mymerge_6_row,"StartL"]<-mymerge[mymerge_6_row,"sEnd_e"] 
mymerge[mymerge_6_row,"EndL"]<-mymerge[mymerge_6_row,"sStart_e"] 
mymerge[mymerge_6_row,"StartR"]<-mymerge[mymerge_6_row,"sEnd_o"] 
mymerge[mymerge_6_row,"EndR"]<-mymerge[mymerge_6_row,"sStart_o"]

#if even rows < odd rows, qend_o > qstart_o, qend_e < qstart_e
#qend_e < qstart_e < qstart_o < qend_o  
mymerge_7_row<-which(mymerge$qend_e < mymerge$qstart_e & mymerge$qstart_e < mymerge$qstart_o & mymerge$qstart_o < mymerge$qend_o)
mymerge[mymerge_7_row,"StartL"]<-mymerge[mymerge_7_row,"sEnd_e"] 
mymerge[mymerge_7_row,"EndL"]<-mymerge[mymerge_7_row,"sStart_e"] 
mymerge[mymerge_7_row,"StartR"]<-mymerge[mymerge_7_row,"sStart_o"] 
mymerge[mymerge_7_row,"EndR"]<-mymerge[mymerge_7_row,"sEnd_o"]

#if even rows < odd rows, qend_o < qstart_o, qend_e > qstart_e
#qstart_e < qend_e < qend_o < qstart_o
mymerge_8_row<-which(mymerge$qstart_e < mymerge$qend_e & mymerge$qend_e < mymerge$qend_o & mymerge$qend_o < mymerge$qstart_o)
mymerge[mymerge_8_row,"StartL"]<-mymerge[mymerge_8_row,"sStart_e"] 
mymerge[mymerge_8_row,"EndL"]<-mymerge[mymerge_8_row,"sEnd_e"] 
mymerge[mymerge_8_row,"StartR"]<-mymerge[mymerge_8_row,"sEnd_o"] 
mymerge[mymerge_8_row,"EndR"]<-mymerge[mymerge_8_row,"sStart_o"]

#############################


######old scheme###########

#defining StartL and EndL

k.len=200

#for rows where qstart_o==1
#myqstart_o_1<-which(mymerge$qstart_o==1)
#mymerge[myqstart_o_1,"StartL"]<-mymerge[myqstart_o_1,"sStart_o"] 
#mymerge[myqstart_o_1,"EndL"]<-mymerge[myqstart_o_1,"sEnd_o"] 

#for rows where qend_o==1
#myqend_o_1<-which(mymerge$qend_o==1)
#mymerge[myqend_o_1,"StartL"]<-mymerge[myqend_o_1,"sEnd_o"] 
#mymerge[myqend_o_1,"EndL"]<-mymerge[myqend_o_1,"sStart_o"] 

#for rows where qstart_e==1
#myqstart_e_1<-which(mymerge$qstart_e==1)
#mymerge[myqstart_e_1,"StartL"]<-mymerge[myqstart_e_1,"sStart_e"] 
#mymerge[myqstart_e_1,"EndL"]<-mymerge[myqstart_e_1,"sEnd_e"] 

#for rows where qend_e==1
#myqend_e_1<-which(mymerge$qend_e==1)
#mymerge[myqend_e_1,"StartL"]<-mymerge[myqend_e_1,"sEnd_e"] 
#mymerge[myqend_e_1,"EndL"]<-mymerge[myqend_e_1,"sStart_e"] 

#head(mymerge[which(mymerge$qend_e==k.len),])

#defining StartR and EndR
#for rows where qstart_o==k.len
#myqstart_o_klen<-which(mymerge$qstart_o==k.len)
#mymerge[myqstart_o_klen,"EndR"]<-mymerge[myqstart_o_klen,"sStart_o"] 
#mymerge[myqstart_o_klen,"StartR"]<-mymerge[myqstart_o_klen,"sEnd_o"] 

#for rows where qend_o==k.len
#myqend_o_klen<-which(mymerge$qend_o==k.len)
#mymerge[myqend_o_klen,"EndR"]<-mymerge[myqend_o_klen,"sEnd_o"] 
#mymerge[myqend_o_klen,"StartR"]<-mymerge[myqend_o_klen,"sStart_o"] 

#for rows where qstart_e==k.len
#myqstart_e_klen<-which(mymerge$qstart_e==k.len)
#mymerge[myqstart_e_klen,"EndR"]<-mymerge[myqstart_e_klen,"sStart_e"] 
#mymerge[myqstart_e_klen,"StartR"]<-mymerge[myqstart_e_klen,"sEnd_e"] 

#for rows where qend_e==k.len
#myqend_e_klen<-which(mymerge$qend_e==k.len)
#mymerge[myqend_e_klen,"StartR"]<-mymerge[myqend_e_klen,"sStart_e"] 
#mymerge[myqend_e_klen,"EndR"]<-mymerge[myqend_e_klen,"sEnd_e"] 

######################
