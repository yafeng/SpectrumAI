library(MSnbase)
library(mzR)
library(stringr)
options(stringsAsFactors = FALSE)

#Specfify the folder where synthetic and experimental spectra (mzML files) are stored
synpep_mzml_path=""
endopep_mzml_path=""

#set working directory which contains the PSM table
setwd("")
#Set the file name of the PSM table 
infile1=""
infile2=""

##double check if the dataframe has column names "SpecFile" , "Peptide", "Charge", "ScanNum", "SpecEValue", "Precursor"
DF1=read.table(infile1,header=TRUE,sep="\t",header=T,comment.char = "",quote = "")
DF2=read.table(infile2,header=TRUE,sep="\t",header=T,comment.char = "",quote = "")

# a list of peptides to draw mnirror image.
# modifications on the peptides must be noted in same way as in the above PSM table.
peptide_list=c()

pdf("Synthetic.Endogenous.peptide.MS2spectrum.mirrorImage.pdf",width=12,height=7,useDingbats = FALSE)

for (peptide in peptide_list){
  
  df.endopep=DF1[DF1$Peptide==peptide,]

  # choose the best scored PSM to draw mirror image
  scanNum.endo = as.integer(df.endopep[which.min (df.endopep$SpecEValue),]$ScanNum)
  precMass.endo = as.integer(df.endopep[which.min (df.endopep$SpecEValue),]$Precursor)
  precCharge.endo = as.integer(df.endopep[which.min (df.endopep$SpecEValue),]$Charge)
  mzml_file.endo=df.endopep[which.min (df.endopep$SpecEValue),]$SpecFile
  

  df.synpep=DF2[DF2$Peptide==Peptide & DF2$Charge == precCharge.endo, ]
  
  if (nrow(df.synpep)==0){sprintf ("%s no corresponding PSMs of synthetic peptide ",peptide);next}
  
  scanNum.syn=df.synpep[which.min(df.synpep$SpecEValue),]$ScanNum
  precMass.syn=df.synpep[which.min(df.synpep$SpecEValue),]$Precursor
  mzml_file.syn=df.synpep[which.min(df.synpep$SpecEValue),]$SpecFile
  
  
  spectraFile1 = paste0(endopep_mzml_path,mzml_file.endo) #endogenous peptide
  spectraFile2 = paste0(synpep_mzml_path,mzml_file.syn) # synthetic peptide
  
  rawdata1 <- openMSfile(spectraFile1,verbose=T)
  rawdata2 <- openMSfile(spectraFile2,verbose=T)
  
  exp_peaks1<-peaks(rawdata1,scan=scanNum.endo)
  exp_peaks2<-peaks(rawdata2,scan=scanNum.syn)
  
  #remove TMT reporter ions peaks 
  exp_peaks1 <-exp_peaks1[exp_peaks1[,1]>131.2,]
  exp_peaks2 <-exp_peaks2[exp_peaks2[,1]>131.2,]
  
  endopep_spectrum=new("Spectrum2",intensity=exp_peaks1[,2],mz=exp_peaks1[,1],centroided=TRUE,
                       precScanNum=scanNum.endo,precursorMz=precMass,precursorCharge=precCharge)
  
  synpep_spectrum=new("Spectrum2",intensity=exp_peaks2[,2],mz=exp_peaks2[,1],centroided=TRUE,
                      precScanNum=scanNum.syn,precursorMz=precMass.syn,precursorCharge=precCharge)
  
  #make mirror image plot, top is endogenous peptide spectrum, bottome is synthetic peptide spectrum
  
  pep.seq = gsub("[^A-Z]","",peptide)
  
  if (grepl("+15.995",peptide)){
    plot(endopep_spectrum,synpep_spectrum,tolerance=0.04,relative=FALSE,
         sequences=c(pep.seq,pep.seq),modifications=c(M=15.995,C=57.021,Nterm=229.163,K=229.163),
         z=seq(1,2,by=1),neutralLoss=NULL,peaks.cex=0,peaks.lwd=1)
  }else {
    plot(endopep_spectrum,synpep_spectrum,tolerance=0.04,relative=FALSE,
         sequences=c(pep.seq,pep.seq),modifications=c(C=57.021,Nterm=229.163,K=229.163),
         z=seq(1,2,by=1),neutralLoss=NULL,peaks.cex=0,peaks.lwd=1)
  }
  legend("topright", "endogenous")
  legend("bottomright", "synthetic")
}

dev.off()
