library(mzR)
library(protViz)
library(stringr)


# Usage:
# RScript SpectrumAI.R /path/to/mzmls/ /path/to/psmtable.txt specfile_colnr /path/to/out.file
# E.g:
# Rscript SpectrumAI.R /home/user/mzmldir /home/user/psms.txt 1 /home/user/specAI_output.txt

# If using interactive environment such as RStudio, simply set the following variable to True,
# and uncomment and change the args below it
use.interactive = F
# args = c('/path/to/mzmls/', '/path/to/psmtable.txt', '/path/to/out.file')


if (use.interactive) {
       source('./Spectra_functions.R')
} else {
        args = commandArgs(trailingOnly = F)  # For scripted use
        # Get script file location when running RScript
        fileflag <- "--file="
        script.fullpath <- sub(fileflag, "", args[grep(fileflag, args)])
        script.dir <- dirname(script.fullpath)
        specfuncfile <- file.path(script.dir, "Spectra_functions.R")
        source(specfuncfile)
        cmargs = commandArgs(trailingOnly = T)
        mzml_path= cmargs[1]
        infile_name = cmargs[2]
        outfile_name = cmargs[3]
}


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

    Spectra_list= vector(mode="list", length=length(unique(as.character(DF[,]$SpectraFile))))
    names(Spectra_list)=unique(as.character(DF[,]$SpectraFile))

    for (i in 1:nrow(DF)){
        spectra_file=as.character(DF[i,]$SpectraFile)
        mzml_file=file.path(mzml_path,spectra_file)
        ScanNum=as.integer(DF[i,]$ScanNum)
        peptide=as.character(DF[i,]$Peptide)
        sub_pos=as.integer(DF[i,]$sub_pos)
        seq=DF[i,]$Sequence
        
        if (is.null(Spectra_list[[spectra_file]])){
            Spectra_list[[spectra_file]]=openMSfile(mzml_file,verbose=T)
        }
        
        exp_peaks<-peaks(Spectra_list[[spectra_file]],scan=ScanNum)
        predicted_peaks = predict_MS2_spectrum(Peptide =  as.character(DF[i,]$Peptide))
        match_ions = match_exp2predicted(exp_peaks, predicted_peaks, tolerance =Frag.ions.tolerance, relative = relative )
        
        maxintensity=max(exp_peaks[,2])
        average_intensity=mean(exp_peaks[,2])
        median_intensity=median(exp_peaks[,2])
        
        DF[i,]$sum.fragmentions.intensity=sum(exp_peaks[,2])
        DF[i,]$maxintensity=maxintensity
        DF[i,]$average_intensity=average_intensity
        DF[i,]$median_intensity=median_intensity
        
        if (nrow(match_ions)==0){next}
        DF[i,]$matched_ions=paste(unique(match_ions$ion),collapse = ",")
        DF[i,]$sum.matchedions.intensity=sum(match_ions$intensity)
        
        if (is.na(sub_pos)){next}
        if (sub_pos==0){next}
        if (sub_pos>DF[i,]$peptide_length){next}
        
        DF[i,]$status="checked"
        
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
            flanking_ions=intersect(flanking_ions,match_ions$ion)
            if (length(flanking_ions)>0){
                flanking_ions_support="YES"
            }
        }else if (n2 == n1){
            flanking_ions=c("y1",paste0("b",as.character(n1-1)))
            flanking_ions=intersect(flanking_ions,match_ions$ion)
            if (length(flanking_ions)>0){
                flanking_ions_support="YES"
            }
        }else {
            flanking_ions_left=c(paste0("b",as.character(n2-1)))
            flanking_ions_left=c(flanking_ions_left,paste0("y",as.character(n1-n2+1)))
            
            flanking_ions_right=c(paste0("b",as.character(n2)))
            flanking_ions_right=c(flanking_ions_right,paste0("y",as.character(n1-n2)))
            
            
            flanking_ions_left=intersect(flanking_ions_left,match_ions$ion)
            flanking_ions_right=intersect(flanking_ions_right,match_ions$ion)
            
            flanking_ions=union(flanking_ions_left,flanking_ions_right)
            if ( length(flanking_ions_left)>0 & length(flanking_ions_left)>0){
                flanking_ions_support="YES"
            }
        }
        
        DF[i,]$flanking_ions_support=flanking_ions_support
        DF[i,]$flanking_ions=paste(flanking_ions,collapse = ",")
        DF[i,]$sum.flanking.ions.intensity=sum(match_ions[match_ions$ion %in% flanking_ions,]$intensity)
        
        if (DF[i,]$sum.flanking.ions.intensity < DF[i,]$median_intensity){DF[i,]$flanking_ions_support <- "NO"}
        
        #fragmentation is not preferable at Cterm side of proline, so only require supporting ions
        if (grepl("P",substr(seq, sub_pos-1, sub_pos))){
            DF[i,]$flanking_ions_support=DF[i,]$ions_support
        }
    }
    return(DF)
}

Frag.ions.tolerance= 0.02 # 0.02 Da tolerance for MS2 fragment ions mass accuracy.
relative=FALSE

# or you can use ppm threshold
# Frag.ions.tolerance= 10 # 10 ppm tolerance for MS2 fragment ions mass accuracy.
# relative=TRUE

start.time <- Sys.time()
start.time

df.psm=read.table(infile_name,sep="\t",header=T,comment.char = "",quote = "")
  #Before running the next command, double check the header names in the input PSM table
  #The df.psm dataframe should have at least the following columns with exactly same names (the order can be different): 
  # "SpectraFile", "ScanNum", "Peptide",  "sub_pos" 

df.output = InspectSpectrum(df.psm)
write.table(df.output,outfile_name,sep="\t",quote=F,row.names=F)

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
