library(mzR)
library(MSnbase)
library(stringr)
options(stringsAsFactors = FALSE)

InspectSpectrum <- function (DF){
    DF$Sequence=gsub("[^A-Z]","",DF$Peptide)
    DF$peptide_length=nchar(DF$Sequence)

    DF$status="skiped"

    DF$ions_support="NO"
    DF$support_ions=""
    DF$sum.supportions.intensity=0

    DF$flanking_ions_support="NO"
    DF$flanking_ions=""
    DF$sum.flanking.ions.intensity=0

    DF$matched_ions=""
    DF$sum.matchedions.intensity=0
    DF$sum.fragmentions.intensity=0
    DF$maxintensity=0
    DF$average_intensity=0
    DF$median_intensity=0

    Spectra_list= vector(mode="list", length=length(unique(as.character(DF[,spectra_file_column]))))
    names(Spectra_list)=unique(as.character(DF[,spectra_file_column]))

    for (i in 1:nrow(DF)){
        spectra_file=as.character(DF[i,]$SpectraFile)
        mzml_file=paste0(mzml_path,spectra_file)
        ScanNum=as.integer(DF[i,]$ScanNum)
        peptide=as.character(DF[i,]$Peptide)
        sub_pos=DF[i,]$sub_pos
        seq=DF[i,]$Sequence
        product_ions_charge=as.integer(DF[i,]$Charge)-1
        if (product_ions_charge>3){product_ions_charge=2} #set maximum charge of product ions
        if (product_ions_charge<1){product_ions_charge=1} #set minimum charge of product ions
        
        if (is.na(sub_pos)){next}
        if (sub_pos==0){next}
        if (is.null(Spectra_list[[spectra_file]])){
            Spectra_list[[spectra_file]]=openMSfile(mzml_file,verbose=T)
        }
        
        exp_peaks<-peaks(Spectra_list[[spectra_file]],scan=ScanNum)
        exp_ms2_spectrum=new("Spectrum2",intensity=exp_peaks[,2],mz=exp_peaks[,1])
        DF[i,]$status="checked"
        
        # use neutralLoss=NULL to disable all neutralLoss
        if (grepl(Oxidation_notation,peptide)){
            match_ions=calculateFragments(seq,exp_ms2_spectrum,z=seq(1,product_ions_charge,1),relative=TRUE,tolerance=Frag.ions.tolerance,
            modifications=c(M=15.995,C=57.021,Nterm=229.163,K=229.163),
            neutralLoss=defaultNeutralLoss(disableAmmoniaLoss=c("K", "R","N","Q")),verbose=F)
        }else {
            match_ions=calculateFragments(seq,exp_ms2_spectrum,z=seq(1,product_ions_charge,1),relative=TRUE,tolerance=Frag.ions.tolerance,
            modifications=c(C=57.021,Nterm=229.163,K=229.163),
            neutralLoss=defaultNeutralLoss(disableAmmoniaLoss=c("K", "R","N","Q")),verbose=F)
        }
        
        if (nrow(match_ions)==0){next}
        DF[i,]$matched_ions=paste(unique(match_ions$ion),collapse = ",")
        
        maxintensity=max(intensity(exp_ms2_spectrum))
        average_intensity=mean(intensity(exp_ms2_spectrum))
        median_intensity=median(intensity(exp_ms2_spectrum))
        
        supportions_intensity=0
        ions_support="NO"
        supportions=""
        
        for (j in 1:nrow(match_ions)){
            type=match_ions[j,]$type
            pos=match_ions[j,]$pos
            ion=match_ions[j,]$ion
            if (type=="b" & pos>=sub_pos){
                ions_support="YES"
                supportions_intensity=supportions_intensity+match_ions[j,]$intensity
                supportions=paste0(supportions,ion,",")
            }else if (type=="y" & pos>nchar(seq)-sub_pos){
                ions_support="YES"
                supportions_intensity=supportions_intensity+match_ions[j,]$intensity
                supportions=paste0(supportions,ion,",")
            }
        }
        
        DF[i,]$ions_support=ions_support
        DF[i,]$support_ions=supportions
        
        DF[i,]$sum.supportions.intensity=supportions_intensity
        DF[i,]$sum.matchedions.intensity=sum(match_ions$intensity)
        DF[i,]$sum.fragmentions.intensity=sum(intensity(exp_ms2_spectrum))
        
        DF[i,]$maxintensity=maxintensity
        DF[i,]$average_intensity=average_intensity
        DF[i,]$median_intensity=median_intensity
        
        #check if it is a noise peak or isotope peak supporting mutant ions
        if (DF[i,]$sum.supportions.intensity < DF[i,]$median_intensity){DF[i,]$ions_support <- "NO"}
        
        flanking_ions_left=c()
        flanking_ions_right=c()
        flanking_ions=c()
        flanking_ions_support="NO"
        n1=DF[i,]$peptide_length
        n2=sub_pos
        if (n2 ==1){
            flanking_ions=c("b1",paste0("y",as.character(n1-1)))
            if ("b1" %in% match_ions$ion){
                flanking_ions_support="YES"
                
            }else if (paste0("y",as.character(n1-1)) %in% match_ions$ion) {
                flanking_ions_support="YES"
            }
        }else if (n2 == n1){
            flanking_ions=c("y1",paste0("b",as.character(n1-1)))
            if ("y1" %in% match_ions$ion){
                flanking_ions_support="YES"
            }else if (paste0("b",as.character(n1-1)) %in% match_ions$ion) {
                flanking_ions_support="YES"
            }
        }else {
            flanking_ions_left=c(paste0("b",as.character(n2-1)))
            flanking_ions_left=c(flanking_ions_left,paste0("y",as.character(n1-n2+1)))
            
            flanking_ions_right=c(paste0("b",as.character(n2)))
            flanking_ions_right=c(flanking_ions_right,paste0("y",as.character(n1-n2)))
            
            flanking_ions=union(flanking_ions_left,flanking_ions_right)
            if ( length(intersect(flanking_ions_left,match_ions$ion))>=1 &
            length(intersect(flanking_ions_right,match_ions$ion))>=1){
                flanking_ions_support="YES"
            }
        }
        
        DF[i,]$flanking_ions_support=flanking_ions_support
        DF[i,]$flanking_ions=paste(flanking_ions,collapse = ",")
        DF[i,]$sum.flanking.ions.intensity=sum(match_ions[match_ions$ion %in% flanking_ions,]$intensity)
        
        #fragmentation is not preferable at Cterm side of proline, so just to loosen the rules
        if (grepl("P",substr(seq, sub_pos-1, sub_pos))){
            DF[i,]$flanking_ions_support=DF[i,]$ions_support
        }
        if (DF[i,]$sum.flanking.ions.intensity < DF[i,]$median_intensity){DF[i,]$flanking_ions_support <- "NO"}
        
    }
}

#set your working directory
setwd()

#set one or multiple paths in which raw files are located
mzml_path_vector=c()

#set one or multiple PSM tables which should be in same order as the path vector
#For example, path of raw spectra in the first PSM file should be found in first path location in the path vector, and so on
infile_vector=c()
#set corresponding output file name for each of input PSM file, in same order
outfile_vector=c()

Frag.ions.tolerance= 10e-6 # 10 ppm tolerance for MS2 fragment ions mass accuracy.
Oxidation_notation = "+15.995"  # define how oxidation on methionine noted in the Peptide Column

start.time <- Sys.time()
start.time
for (n in 1:length(mzml_path_vector)){
  mzml_path=mzml_path_vector[n]
  infile=infile_vector[n]
  outfile=outfile_vector[n]
  
  df.psm=read.table(infile,sep="\t",header=T,comment.char = "",quote = "")
  #Before running the next command, double check the header names in the input PSM table
  #The df.psm dataframe should have at least the following columns with exactly same names (the order can be different): 
  # "SpectraFile", "ScanNum", "Peptide", "Charge", "sub_pos" 
  InspectSpectrum(df.psm)
  write.table(DF,outfile,sep="\t",quote=F,row.names=F)
}

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
