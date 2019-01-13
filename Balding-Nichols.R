setwd("/Users/miguelguardado/Desktop/SFSU/BIO_358")
AAFrequnecyTable=read.table('Allele-Freqs-PopStats.txt',sep=",")
AAGenotypeTable=read.csv('Genotypes-PopStats.csv')

#Get the observed Frequnecy of one STR to use to compare with the expected frequencies
AAGenotypeTable=as.matrix(AAGenotypeTable)
AAGenotypeTable=AAGenotypeTable[1:342,3:28]
Names<-data.matrix(AAFrequnecyTable$V1)
STR1=data.matrix(AAFrequnecyTable$V2)
STR2=data.matrix(AAFrequnecyTable$V3)
STR3=data.matrix(AAFrequnecyTable$V4)
STR4=data.matrix(AAFrequnecyTable$V5)
STR5=data.matrix(AAFrequnecyTable$V6)
STR6=data.matrix(AAFrequnecyTable$V7)
STR7=data.matrix(AAFrequnecyTable$V8)
STR8=data.matrix(AAFrequnecyTable$V9)
STR9=data.matrix(AAFrequnecyTable$V10)
STR10=data.matrix(AAFrequnecyTable$V11)
STR11=data.matrix(AAFrequnecyTable$V12)
STR12=data.matrix(AAFrequnecyTable$V13)
STR13=data.matrix(AAFrequnecyTable$V14)
STRArr<-cbind(Names,STR1,STR2,STR3,STR4,STR5,STR6,STR7,STR8,STR9,STR10,STR11,STR12,STR13)
STRArr<-STRArr[5:106,1:14]

Obsmean<-mean(as.numeric(STR[103,1:13]))

CSF1POfreq=data.matrix(AAFrequnecyTable$V2)
CSF1POfreq=cbind(refvar,CSF1POfreq)
obsfreq=as.numeric(AAFrequnecyTable[103,2])

#Algorithim for Calculating Genotype frequnecy under Hardy-Weinberg
HWEProb<- function(GenoType,FrequencyTable){
  LocusProb=1
  geno1<-as.numeric(FrequencyTable[which(as.numeric(FrequencyTable[,1])==as.numeric(GenoType[1])),2])
  geno2<-as.numeric(FrequencyTable[which(as.numeric(FrequencyTable[,1])==as.numeric(GenoType[2])),2])
  if(as.numeric(GenoType[1])==as.numeric(GenoType[2])){  
   LocusProb=geno1*geno2
    }else{
   LocusProb=2*geno1*geno2
 }
 return(LocusProb)
}


#Algorithim for calculating Genotype frequnecy under Balding Nichols
BNDProb<-function(GenoType,FrequencyTable,degree){
  Locusprob=1
  geno1<-as.numeric(FrequencyTable[which(as.numeric(FrequencyTable[,1])==as.numeric(GenoType[1])),2])
  geno2<-as.numeric(FrequencyTable[which(as.numeric(FrequencyTable[,1])==as.numeric(GenoType[2])),2])
  if(as.numeric(GenoType[1])==as.numeric(Genotype[2])){
    Locusprob=(degree*geno1)+(1-degree)*((geno1)^2)
  }else{
    Locusprob=(1-degree)*2*geno1*geno2
  } 
  return(Locusprob)
}

#Simulate though the first STR to calculate the expected genotype frequency using the two models.
GenoFact1<-data.matrix(AAGenotypeTable$CSF1PO)
GenoFact2<-data.matrix(AAGenotypeTable$X)
HWEprob=rep(0,342)
for(i in 2:344){
  Genotype<-c(GenoFact1[i],GenoFact2[i])
  HWEprob[i-1]=HWEProb(Genotype,CSF1POfreq)
}

STRHWE=rep(0,13)
STRBND=rep(0,13)
STRBND2=rep(0,13)
STRBND4=rep(0,13)
STRBND6=rep(0,13)
SREBND8=rep(0,13)
for(STR in 1:13){
  HWEprob=rep(0,342)
  BND<-rep(0,342)
  BND2<-rep(0,342)
  BND4<-rep(0,342)
  BND6<-rep(0,342)
  BND8<-rep(0,342)
  
  GenoFact1=AAGenotypeTable[1:342,STR*2-1]
  GenoFact2=AAGenotypeTable[1:342,STR*2]
  STRFreq<-cbind(STRArr[1:97,1],STRArr[1:97,STR+1])
for(i in 1:342){
  Genotype<-c(GenoFact1[i],GenoFact2[i])
  HWEprob[i]=HWEProb(Genotype,STRFreq)
  BND[i]<-BNDProb(Genotype,STRFreq,.10)
  BND2[i] <- BNDProb(Genotype,STRFreq,.02)
  BND4[i] <- BNDProb(Genotype,STRFreq,.04)
  BND6[i] <- BNDProb(Genotype,STRFreq,.06)
  BND8[i] <- BNDProb(Genotype,STRFreq,.08) 

}

HWEProbMean[STR]<-mean(HWEprob)
BNDMean[STR]<-mean(BND)
BND2Mean[STR]<-mean(BND2)
BND4Mean[STR]<-mean(BND4)
BND6Mean[STR]<-mean(BND6)
BND8Mean[STR]<-mean(BND8)
}
FinHweMean<-(Obsmean-mean(HWEProbMean))^2
FinBNDMean<-(Obsmean-mean(BNDMean))^2
FinBND2Mean<-(Obsmean-mean(BND2Mean))^2
FinBND4Mean<-(Obsmean-mean(BND4Mean))^2
FinBND6Mean<-(Obsmean-mean(BND6Mean))^2
FinBND8Mean<-(Obsmean-mean(BND8Mean))^2

color<-c("red",'blue',"green","black","orange","purple")

FinData<-c(FinHweMean,FinBND2Mean,FinBND4Mean,FinBND6Mean,FinBND8Mean,FinBNDMean)
plot((FinData),main="Means Squared Test of Observed and Exprected Genotype Frequnecies",xlab="Simulations Ran",ylab="(Observed vs Expected)^2",col=color)
legend("topleft", 
       legend=c("Hardy-Weignberg","Balding-Nichols(Theta=.02)","Balding-Nichols(Theta=.04)","Balding-Nichols(Theta=.06)","Balding-Nichols(Theta=.08)","Balding-Nichols(Theta=.10)"), 
       col=color, lty=1:6, cex=0.8)


