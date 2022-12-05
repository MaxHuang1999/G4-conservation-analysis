##This script was used for calculating the G4score for each g4seq.

###extract from the program G4Hunter
###################################################################

require(S4Vectors)

G4translate <- function(x)		
  # x is a Rle of a sequence
{
  xres=x
  runValue(xres)[runValue(x)=='C' & runLength(x)>3] <- -4
  runValue(xres)[runValue(x)=='C' & runLength(x)==3] <- -3
  runValue(xres)[runValue(x)=='C' & runLength(x)==2] <- -2
  runValue(xres)[runValue(x)=='C' & runLength(x)==1] <- -1
  runValue(xres)[runValue(x)=='G' & runLength(x)>3] <- 4
  runValue(xres)[runValue(x)=='G' & runLength(x)==3] <- 3
  runValue(xres)[runValue(x)=='G' & runLength(x)==2] <- 2
  runValue(xres)[runValue(x)=='G' & runLength(x)==1] <- 1
  runValue(xres)[runValue(x)!='C' & runValue(x)!='G'] <- 0
  Rle(as.numeric(xres))
}

G4Hscore <- function(y)		
  # y can be DNAString or a DNAStringSet or a simple char string.
{
  y2 <- Rle(strsplit(as.character(y),NULL)[[1]])
  y3 <- G4translate(y2)
  abs(round(mean(y3),3))
}

ScoreCalculate=function(inputfile,outputfile){
    g4df=read.csv2(file=inputfile,sep=',',header = F,stringsAsFactors= F)
    colnames(g4df)=c("AccessionNumber","pattern","Start","End","Seq")
    score_content=as.vector(sapply(g4df$Seq, G4Hscore))
    g4df$G4HunterScore=score_content
    write.table(g4df,file = outputfile,quote=F,sep=",",col.names = T,row.names = F)
}

args=commandArgs(T)
ScoreCalculate(args[1],args[2])

#setwd("/home/zhluo/Project/lizhen_lncRNA/Allvirus/shuffle_GNG/")
#ScoreCalculate("GNG_adjust_pattern_longseq.csv","GNG_adjust_pattern_longseq_g4score.csv")
