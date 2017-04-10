library(mzR)
library(MSnbase)
library(stringr)
options(stringsAsFactors = FALSE)

#set working directory which contains the PSM table
setwd("")

mzml_path=""

infile1=""

##double check if the dataframe has column names "SpecFile" , "Peptide", "Charge", "ScanNum", "SpecEValue", "Precursor"
DF = read.table(infile1,header=T,comment.char = "",quote="",sep = "\t")
DF$Sequence=gsub("[^A-Z]","",DF$Peptide)

Oxidation = "+15,995"
Frag.ions.tolerance= 0.02 # unit Dalton

pdf("PSM.annotated.spectra.pdf",width=12,height=7,useDingbats = FALSE)
for (i in 1:nrow(DF)){
  spectra_file=as.character(DF[i,]$SpectraFile)
  mzml_file=paste0(mzml_path,spectra_file)
  ScanNum=as.integer(DF[i,]$ScanNum)
  peptide=as.character(DF[i,]$Peptide)
  seq=DF[i,]$Sequence
  precMass=DF[i,]$Precursor
  precCharge=DF[i,]$Charge
  
  
  rawdata <- openMSfile(mzml_file,verbose=T)
  
  exp_peaks<-peaks(rawdata,scan=ScanNum)
  
  exp_spectrum=new("Spectrum2",intensity=exp_peaks[,2],mz=exp_peaks[,1],centroided=TRUE,
               precScanNum=ScanNum,precursorMz=precMass,precursorCharge=precCharge)

  
  if (grepl(Oxidation,peptide)){
    
    theoretical_peaks=calculateFragments(seq,modifications=c(M=15.995,C=57.021,Nterm=229.163,K=229.163),
                                        z=1,neutralLoss=NULL)
    theoretical_peaks$intensity=10000
    theoretical_spectrum = new("Spectrum2",intensity=theoretical_peaks$intensity,mz=theoretical_peaks$mz,centroided=TRUE)
    
    plot(exp_spectrum,theoretical_spectrum,tolerance=Frag.ions.tolerance,relative=FALSE,
         sequences=c(seq,seq),modifications=c(M=15.995,C=57.021,Nterm=229.163,K=229.163),
         z=1,neutralLoss=NULL,peaks.cex=0,peaks.lwd=1)
  }else {
    
    theoretical_peaks=calculateFragments(seq,modifications=c(C=57.021,Nterm=229.163,K=229.163),
                                         z=1,neutralLoss=NULL)
    theoretical_peaks$intensity=10000
    theoretical_spectrum = new("Spectrum2",intensity=theoretical_peaks$intensity,mz=theoretical_peaks$mz,centroided=TRUE)
    
    plot(exp_spectrum,theoretical_spectrum,tolerance=Frag.ions.tolerance,relative=FALSE,
         sequences=c(seq,seq),modifications=c(C=57.021,Nterm=229.163,K=229.163),
         z=1,neutralLoss=NULL,peaks.cex=0,peaks.lwd=1)
  }
  legend("topright", "experimental")
  legend("bottomright", "theoretical")
}

dev.off()
