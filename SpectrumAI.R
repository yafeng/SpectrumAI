library(mzR)
library(protViz)
library(stringr)

options(stringsAsFactors = FALSE)

predict_MS2_spectrum <- function (Peptide){
  pepMSGF=gsub("[^A-Z]","",Peptide)
  size = nchar(pepMSGF)
  
  pepMSGFMods=Peptide
  pepMSGFMods=str_replace_all(pepMSGFMods,pattern = "([\\+,0-9,:,\\.]+)" ,replacement = "\\1:" )
  
  l='none'
  while (l != pepMSGFMods){
    l=pepMSGFMods
    pepMSGFMods=sub("([^\\+,0-9,:])([^\\+,0-9,:])", "\\1:\\2", pepMSGFMods, perl=TRUE)  
  }
  pepMSGFMods=t(str_split(pepMSGFMods,pattern = ":",simplify = T))
  pepMSGFMods=data.frame(str_split_fixed(pepMSGFMods,pattern = "\\+",2), stringsAsFactors = F)
  pepMSGFMods[pepMSGFMods[,2]=="",2]=0
  pepMSGFMods[,2]=  as.double(pepMSGFMods[,2])
  pepMSGFMods[2,2] = pepMSGFMods[1,2] + pepMSGFMods[2,2]
  
  pepMSGFMods=pepMSGFMods[nchar(pepMSGFMods[,1])>0,]
  pepMSGFWeights = protViz::aa2mass(pepMSGF)[[1]]
  pepMSGFWeights =  pepMSGFWeights + t(as.double(pepMSGFMods[,2]))
  
  
  ions=fragmentIon(pepMSGFWeights)
  ions <- data.frame(mz=c(ions[[1]]$b,ions[[1]]$y), 
                     ion=c(paste0(rep("b",size),1:size),paste0(rep("y",size),1:size)),
                     type=c(rep("b",size),rep("y",size)),
                     pos =rep(1:size,2),
                     z=rep(1,size*2))
  
  if (grepl("K+229.163",Peptide,fixed = TRUE)){ ions[nrow(ions),1] = ions[nrow(ions),1] - 229.163}
  
  return (ions)
}

match_exp2predicted <- function (exp_peak,pred_peak,tolerance = 0.02, relative=FALSE){
  pred_peak$error = apply(pred_peak,1,function(x) min(abs(as.numeric(x[1])-exp_peak[,1])))
  pred_peak$intensity = apply(pred_peak,1,function(x) exp_peak[which.min(abs(as.numeric(x[1]) - exp_peak[,1])),2])
  pred_peak$ppm = round(pred_peak$error/pred_peak$mz*1000000,2)
  
  if (relative){
    match_ions = pred_peak[pred_peak$ppm <tolerance,]
  }else{match_ions = pred_peak[pred_peak$error <tolerance,]
  }
  
  return (match_ions)
}


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
}

#set your working directory
setwd()

#set one or multiple paths in which raw files are located
mzml_path= ""

#set one or multiple PSM tables which should be in same order as the path vector
#For example, path of raw spectra in the first PSM file should be found in first path location in the path vector, and so on
infile_name =""
#set corresponding output file name for each of input PSM file, in same order
outfile_name =""

Frag.ions.tolerance= 0.02 # 0.02 Da tolerance for MS2 fragment ions mass accuracy.
relative=FALSE

# or you can use ppm threshold
# Frag.ions.tolerance= 10 # 10 ppm tolerance for MS2 fragment ions mass accuracy.
# relative=TRUE

start.time <- Sys.time()
start.time

df.psm=read.table(infile,sep="\t",header=T,comment.char = "",quote = "")
  #Before running the next command, double check the header names in the input PSM table
  #The df.psm dataframe should have at least the following columns with exactly same names (the order can be different): 
  # "SpectraFile", "ScanNum", "Peptide",  "sub_pos" 

InspectSpectrum(df.psm)
write.table(df.psm,outfile,sep="\t",quote=F,row.names=F)

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
