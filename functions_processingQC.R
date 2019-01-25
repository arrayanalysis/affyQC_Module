#=============================================================================#
# ArrayAnalysis - affyAnalysisQC                                              #
# a tool for quality control and pre-processing of Affymetrix array data      #
#                                                                             #
# Copyright 2010-2011 BiGCaT Bioinformatics                                   #
#                                                                             #
# Licensed under the Apache License, Version 2.0 (the "License");             #
# you may not use this file except in compliance with the License.            #
# You may obtain a copy of the License at                                     #
#                                                                             #
# http://www.apache.org/licenses/LICENSE-2.0                                  #
#                                                                             #
# Unless required by applicable law or agreed to in writing, software         #
# distributed under the License is distributed on an "AS IS" BASIS,           #
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.    #
# See the License for the specific language governing permissions and         #
# limitations under the License.                                              #
#=============================================================================#

##################
## getArrayType ##
##################

# determine array type used
#--------------------------

getArrayType <- function(Data) {

  aType <- "PMMM"  

  # Test whether the dataset is of aType "PM-only"
  mismatches <- mm(Data[,1])  

  if(is.null(mismatches)) {
    # mm does not exist
    aType <- "PMonly"
  } else {
    if(sum(is.na(mismatches))>=length(mismatches)/2){ 
      # mm is always NA or there are more NA values in the mm probes than 
	  # defined values (assuming these would just be controls)
      aType <- "PMonly"
    } else {
      matches <- pm(Data[,1])
      notNA <- !is.na(mismatches) & !is.na(matches)
      if(sum(mismatches[notNA]!=matches[notNA])==0){
        # MM contains a copy of PM, which indicates a PMonly array
        aType <- "PMonly"
      }
    }
  }

  cat("The arrays are determined to contain", ifelse(aType=="PMMM",
  "perfect match and mismatch probes\n", "perfect match probes only\n"))
  
  return(aType)
}


#######################
## addStandardCDFenv ##
#######################

addStandardCDFenv <- function(Data, overwrite=FALSE) {
  #if overwrite is FALSE a cdf environment will be kept if already loaded
  #if overwrite is TRUE it will always be overwritten (unless none is found)
  #the first option would be the regular one, used to add cdf environments where
  #   automatic loading failed
  #the second could be used to set back updated cdf files to the standard ones

  #start with r0 just in case this could exist
  rev <- 0

  #initial value
  CDFenv <- 0

  # recall which cdfName was added, in case no other one is found (set back even
  #   if it does not exist)
  presetCDF <- Data@cdfName

  #check whether environment is already correct
  suppressWarnings(try(CDFenv <- getCdfInfo(Data),TRUE))

  #try the annotation plus cdf as a cdf name
  if ((class(CDFenv)!="environment") | overwrite) {
    CDFenv <- 0 #needed for cases where overwrite is TRUE, but CDFenv already
                # was an environment
    Data@cdfName<-paste(Data@annotation,".cdf",sep="")
    suppressWarnings(try(CDFenv <- getCdfInfo(Data),TRUE))
    #if no succes also try without a dot
    if (class(CDFenv)!="environment") {
      Data@cdfName<-paste(Data@annotation,"cdf",sep="")
      suppressWarnings(try(CDFenv <- getCdfInfo(Data),TRUE))
    }
  }

  # don't run the loop if CDF is already known up front, or correct one has been
  # found
  while((class(CDFenv)!="environment") & (rev < 10)) {
    Data@cdfName<-paste(Data@annotation,".r",rev,"cdf",sep="")
    suppressWarnings(try(CDFenv <- getCdfInfo(Data),TRUE))
    rev <- rev + 1
  }

  if ((class(CDFenv)!="environment")) {
    Data@cdfName <- presetCDF
    warning("could not automatically retrieve CDF environment for this chip type - object kept as is")
  }

  cat("current cdf environment loaded:",Data@cdfName,"\n")
  return(Data)
}


######################
## addUpdatedCDFenv ##
######################

addUpdatedCDFenv <- function(Data, species=NULL, type="ENSG") {
  # note: this function will add an updated cdf environment to the data object
  # and will overwrite a possible already loaded environment, unless no updated
  # cdf environment is found

  # developer's note: it may be of interest to find out whether available
  # species and types can be retrieved automatically from the brainarray website

  if(is.null(species) || (species=="")) stop("The species must be provided")

  types <- c("ENTREZG","REFSEQ","ENSG","ENSE","ENST","VEGAG","VEGAE","VEGAT",
          "TAIRG","TAIRT","UG","MIRBASEF","MIRBASEG")
  if(!tolower(type) %in% tolower(types)) {
    stop("selected type not valid, select from ", paste(types, collapse=" "))
  } else {
    type <- types[match(tolower(type),tolower(types))]
  }

  spp <- c("Ag","At","Bt","Ce","Cf","Dr","Dm","Gg","Hs","MAmu","Mm","Os","Rn",
           "Sc","Sp","Ss")
  names(spp) <- c("Anopheles gambiae","Arabidopsis thaliana","Bos taurus",
           "Caenorhabditis elegans","Canis familiaris", "Danio rerio",
           "Drosophila melanogaster","Gallus gallus","Homo sapiens",
           "Macaca mulatta","Mus musculus", "Oryza sativa","Rattus norvegicus",
           "Saccharomyces cerevisiae","Schizosaccharomyces pombe","Sus scrofa")
  if(tolower(species) %in% tolower(names(spp)))
       species <- spp[tolower(names(spp))==tolower(species)]
  if(!tolower(species) %in% tolower(spp)) {
    stop("selected species not valid, select from:\n",
          paste(names(spp), collapse="\n"), "\nor abbreviated as ",
          paste(spp,collapse=" "))
  } else {
    species <- spp[match(tolower(species),tolower(spp))]
  }

  #initial value
  CDFenv <- 0

  # recall which cdfName was added, in case no updated one is found (set back
  # even if it does not exist)
  presetCDF <- Data@cdfName

  #try to find updated cdf file of choice***
  print(Data@cdfName<-paste(Data@annotation,species,type,sep="_"))

  suppressWarnings(try(CDFenv <- getCdfInfo(Data),TRUE))
  #try without a version number
  print(Data@cdfName<-paste(gsub("v[0-9]$","",Data@annotation),species,type,sep="_"))
  suppressWarnings(try(CDFenv <- getCdfInfo(Data),TRUE)) 
  

  #if it hasn't loaded, try to download ***
  if ((class(CDFenv)!="environment")) {
    print("***")
    print(Data@annotation)
    print(species)
    print(type)
    install.packages(tolower(paste(Data@annotation,species,type,"cdf",sep="")),
      repos="http://brainarray.mbni.med.umich.edu/bioc")
    suppressWarnings(try(CDFenv <- getCdfInfo(Data),TRUE))
  }

  #if it hasn't loaded, try to download without version number
  if ((class(CDFenv)!="environment")) {
    install.packages(tolower(paste(gsub("v[0-9]$","",Data@annotation),species,type,"cdf",sep="")),
      repos="http://brainarray.mbni.med.umich.edu/bioc")
    suppressWarnings(try(CDFenv <- getCdfInfo(Data),TRUE))
  }
  
  if ((class(CDFenv)!="environment")) {
    Data@cdfName <- presetCDF
    warning("Could not automatically retrieve CDF environment for this chip type - object kept as is")
  }

  cat("current cdf environment loaded:",Data@cdfName,"\n")
  return(Data)
}


####################
## colorsByFactor ##
####################

#create colors for the plots and the legends
#-------------------------------------------

colorsByFactor <- function(experimentFactor) {
  
  #check whether a factor has been provided
  if(class(experimentFactor)!="factor") stop("Parameter 'experimentFactor' must be of class 'factor'")

  if(length(levels(experimentFactor))==1) {
    #if there is only one group (or no groups are provided) take equally spread colors over the rainbow palette
    plotColors <- rainbow(length(experimentFactor),s=.8,v=.7)
	#set group legend color to white, as there is not a specific group color
	legendColors <- "white"
  } else {
    #compute the number of colors needed for each class
    tab.tmp <- table(experimentFactor)

    #set the two extreme colors for each class
    colors.light <- rainbow(length(levels(experimentFactor)),s=1-sapply(tab.tmp,min,5)*.1)
    colors.dark <- rainbow(length(levels(experimentFactor)),v=1-sapply(tab.tmp,min,5)*.14)

    #create the colors to plot, and colors for the legend (average one per experimental group)
    plotColors <- NULL
    legendColors <- NULL
    for(l in 1:length(levels(experimentFactor))) {
      colorFun <- colorRampPalette(c(colors.light[l],colors.dark[l]))
      tmpColors <- colorFun(tab.tmp[l])
      plotColors[experimentFactor==levels(experimentFactor)[l]] <- tmpColors
      legendColors[l] <- tmpColors[ceiling(length(tmpColors)/2)]
    }
  }
  return(list(plotColors=plotColors,legendColors=legendColors))

}


###################
## deduceSpecies ##
###################

#try to find the species of the chiptype
#---------------------------------------

deduceSpecies <- function(descr=NULL) {

  organism <- ""
  
  if(!is.null(descr) && (descr!="")) {
    try({lib <- paste(descr,".db",sep="")
      eval(parse("",-1,paste("require(",lib,")",sep="")))
      organism <- eval(parse("",-1,paste(descr,"ORGANISM",sep="")))},TRUE)
    if(organism=="") {
      descr <- tolower(descr)
      if(length(grep("^nugomm",descr)) > 0) organism <- "Mus musculus"
      if(length(grep("^nugohs",descr)) > 0) organism <- "Homo Sapiens"
      if(length(grep("^hgu133plus2",descr)) > 0) organism <- "Homo Sapiens"      
      if(length(grep("^hugene",descr)) > 0) organism <- "Homo Sapiens"
      if(length(grep("^mogene",descr)) > 0) organism <- "Mus musculus"
      if(length(grep("^ragene",descr)) > 0) organism <- "Rattus norvegicus"
    }
  }
  
  return(organism)
}


###################
## normalizeData ##
###################

#normalize the data set
#----------------------

normalizeData <- function(Data, normMeth="", perGroup=FALSE, experimentFactor=NULL, 
  customCDF=TRUE, species=NULL, CDFtype=NULL, aType=NULL, isOligo = FALSE, WIDTH=1000, HEIGHT=1414) {

  if((normMeth=="") || is.null(normMeth)) stop("normMeth, the requested normalization method, must be provided")
  normMeth <- toupper(normMeth)
  
	if(customCDF) {
		if(is.null(species) || species=="") stop("When customCDF is required, the species must be provided")
		if(is.null(CDFtype) || CDFtype=="") stop("When customCDF is required, the CDFtype must be provided")
	}
	if(perGroup) {
		if(is.null(experimentFactor)) stop("When normalization per group is requested, experimentFactor must be provided")
	}
  if((normMeth=="PLIER") && (is.null(aType))) stop("When selecting PLIER normalization, aType must be provided")
  if((normMeth=="GCRMA") && (is.null(aType))) stop("When selecting GCRMA normalization, aType must be provided")
  
	#if customCDF option is chosen, apply to copy of Data, in order not to change the original data object
	Data.copy <- Data
	if(customCDF){
		print ("Change CDF before pre-processing")
		if(!isOligo)
			Data.copy <- addUpdatedCDFenv(Data.copy, species, CDFtype)
		else
			warning("CUSTOM CDF NOT AVAILABLE FOR OLIGO PACKAGE. Standard CDF kept!")
	}
	print ("Pre-processing is running")

  nGroups <- 1
  if(perGroup) {
    nGroups <- max(1,length(levels(experimentFactor)))
    if(nGroups==1) warning("normalization per group requested, but no groups indicated in data set")
  }
	
  #if per group normalization required, or a method selected that does not return an ExpressionSet object,
  #make a model of class ExpressionSet to paste real values in, use the relatively fast RMA method
  #note that binding of ExpressionSet objects is NOT possible
  if((nGroups>1)) { # || (normMeth=="MAS5")) {
    normData <- rma(Data.copy)
    exprs(normData)[] <- NA
  }
  
  for(group in 1:nGroups) {
    if(nGroups==1) {
      Data.tmp <- Data.copy
    } else {
      Data.tmp <- Data.copy[,experimentFactor==(levels(experimentFactor)[group])]
    }
    switch(normMeth, 
      "MAS5" = {
      #doesn't work
      normData.tmp <- mas5(Data.tmp) 
      },
      "GCRMA" = {
      if(customCDF) {
        #probe library needed, first try whether this has been intalled, otherwise do so
	#***
        probeLibrary <- tolower(paste(Data@annotation,species,CDFtype,"probe",sep=""))
        loaded <- suppressWarnings(try(eval(parse("",-1,paste("library(",probeLibrary,")", sep=""))),TRUE))
        if(class(loaded)=="try-error") {
          install.packages(probeLibrary, repos="http://brainarray.mbni.med.umich.edu/bioc")
        }
      }
      if(aType == "PMMM") ntype = "fullmodel"
      if(aType == "PMonly") ntype = "affinities" # good results if most of the genes are not expressed
	  normData.tmp <- gcrma(Data.tmp, type=ntype, fast = FALSE)
      },
      "RMA" = {
      normData.tmp <- rma(Data.tmp)
      },
      "PLIER" = {
      if(aType == "PMMM") ntype = "together"
      if(aType == "PMonly") ntype = "pmonly"
      normData.tmp <- justPlier(Data.tmp, normalize=TRUE, norm.type = ntype)
      }
    )
    if(nGroups==1) {
        normData <- normData.tmp
		if(normMeth=="MAS5") exprs(normData)<-log2(exprs(normData))
    } else {
      try(
	  if(normMeth=="MAS5"){
		exprs(normData)[,match(sampleNames(normData.tmp), sampleNames(normData))] <- log2(exprs(normData.tmp))
	  }else{
		exprs(normData)[,match(sampleNames(normData.tmp), sampleNames(normData))] <- exprs(normData.tmp)
	  },TRUE)
    }
    rm(normData.tmp, Data.tmp)
  }

  #create an 'inter sheet'  
  png(file="Cover_2.png", width=WIDTH, height=HEIGHT)
  plot(c(0,2), type = 'n', ann = FALSE, axes = FALSE, 
     frame.plot = TRUE, xlim = c(0, 2), ylim = c(0,2))
  text(1,1,"Pre-processing of Raw Data\n\n\n",cex=3)
  
  #isOligo test
  result <- try(  text(1,1,paste("\n\nMethod: ",normMeth,"\nAnnotation: ",Data.copy@cdfName),cex=2.5), silent=TRUE)
  if (class(result) == "try-error") {
	  text(1,1,paste("\n\nMethod: ",normMeth,"\nAnnotation: ",Data.copy@annotation),cex=2.5)
  }
  #isOligo test end

  if(perGroup) text(1,1,paste("\n\n\n\n\nNormalization per experimental group"),cex=2.5)
  dev.off()
  
	rm(Data.copy)
	return(normData)
}

#####################
## computePMAtable ##
#####################

#prepare and return a table of PMA calls
#---------------------------------------

computePMAtable <- function(Data, customCDF=TRUE, species=NULL, CDFtype=NULL) {

	if(customCDF) {
		if(is.null(species)) stop("When customCDF is used, the species must be provided")
		if(is.null(CDFtype)) stop("When customCDF is used, the CDFtype must be provided")
  }

  Data.copy <- Data  
  if(customCDF) {
    Data.copy <- addUpdatedCDFenv(Data.copy, species, CDFtype)
	}
  
  PMAtable <- NULL
  
  #ocmpute the calls using the detection.p.val function from simpleaffy
  try(PMAtable <- detection.p.val(Data.copy)$call,TRUE)
  
  if(!is.null(PMAtable)) {
    #add column of IDs to table
    PMAtable <- cbind(rownames(PMAtable),PMAtable)

    #remove "_at" from custom probeset IDs to get to real ID
    if(customCDF) {
      control_rows <- grep("affx",tolower(PMAtable[,1]))
      if(length(control_rows) > 0) {
        PMAtable[-control_rows,1] <- substring(PMAtable[-control_rows,1],1,nchar(PMAtable[-control_rows,1])-3)
      } else {
        PMAtable[,1] <- substring(PMAtable[,1],1,nchar(PMAtable[,1])-3)
      }
    }
    
    #add column names to PMAtable
    colnames(PMAtable)[1] <- ifelse(customCDF,paste(CDFtype,"_ID",sep=""),"Probeset_ID")
  } else {
    warning("PMA table could not be computed for this arraytype")
  }
  
  return(PMAtable)
}


#########################
## createNormDataTable ##
#########################

#prepare and export the normalized data table
#--------------------------------------------

createNormDataTable <- function(normData, customCDF=NULL, species=NULL, CDFtype=NULL) {

  if(is.null(customCDF)) stop("The customCDF parameter must be provided")
	if(customCDF) {
		if(is.null(CDFtype)) stop("When customCDF is used, the CDFtype must be provided")
		if(species=="" || is.null(species)) {
		  warning("Species has not been set and custom cdf requested, attempting to deduce species for chip type")
		  species <- deduceSpecies(rawData@annotation)
		}
		if(species=="" || is.null(species)){
			warning("Could not define species; the CDF will not be changed")
			customCDF<-FALSE
		}
	}
  
  #add column of IDs and normalized data to normDataTable
  normDataTable<-cbind(rownames(exprs(normData)),exprs(normData))

  #remove "_at" from custom probeset IDs to get to real ID
  if(customCDF) {
    control_rows <- grep("affx",tolower(normDataTable[,1]))
    if(length(control_rows)>0) {
      normDataTable[-control_rows,1]<-substring(normDataTable[-control_rows,1],1,nchar(normDataTable[-control_rows,1])-3)
    } else {
      normDataTable[,1]<-substring(normDataTable[,1],1,nchar(normDataTable[,1])-3)
    }
  }
  
  #add column names to normDataTable
  colnames(normDataTable)<-c(ifelse(customCDF,paste(CDFtype,"_ID",sep=""),"Probeset_ID"),colnames(exprs(normData)))

  #add gene name and description in case ensembl IDs have been used (otherwise there is no 1 to 1 mapping)
  if(customCDF && CDFtype=="ENSG") {
    #load gene name and description annotations
    library(biomaRt)

    spName <- ""
    if(species=="Ag" || species=="Anopheles gambiae") spName <- "agambiae"
    if(species=="At" || species=="Arabidopsis thaliana") spName <- "athaliana"
    if(species=="Bt" || species=="Bos taurus") spName <- "btaurus"
    if(species=="Ce" || species=="Caenorhabditis elegans") spName <- "celegans"
    if(species=="Cf" || species=="Canis familiaris") spName <- "cfamiliaris"
    if(species=="Dr" || species=="Danio rerio") spName <- "drerio"
    if(species=="Dm" || species=="Drosophila melanogaster") spName <- "dmelanogaster"
    if(species=="Gg" || species=="Gallus gallus") spName <- "ggallus"
    if(species=="Hs" || species=="Homo sapiens") spName <- "hsapiens"
    if(species=="MAmu" || species=="Macaca mulatta") spName <- "mmulatta"
    if(species=="Mm" || species=="Mus musculus") spName <- "mmusculus"
    if(species=="Os" || species=="Oryza sativa") spName <- "osativa"
    if(species=="Rn" || species=="Rattus norvegicus") spName <- "rnorvegicus"
    if(species=="Sc" || species=="Saccharomyces cerevisiae") spName <- "scerevisiae"
    if(species=="Sp" || species=="Schizosaccharomyces pombe") spName <- "spombe"
    if(species=="Ss" || species=="Sus scrofa") spName <- "sscrofa"
    
    try(ensembl <- useMart("ensembl", dataset = paste(spName,"_gene_ensembl",sep="")))
    if(exists("ensembl")) {
      try(annotationTable<-getBM(attributes=c("ensembl_gene_id","external_gene_name","description"), mart=ensembl, uniqueRows=TRUE),TRUE)
    }
    if(exists("annotationTable")) {
      normDataTable <- as.data.frame(normDataTable,stringsAsFactors=FALSE)
      suppressWarnings(normDataTable<-cbind(normDataTable,annotationTable[match(normDataTable[,1],annotationTable[,1]),2:(dim(annotationTable)[2])]))
      normDataTable[,2:(dim(exprs(normData))[2]+1)] <- apply(normDataTable[,2:(dim(exprs(normData))[2]+1),drop=FALSE],2,as.numeric)
    } else {
      warning("No gene names and annotation could be retrieved from BioMart for this species or no connection could be established, gene information not added to normalized data table")
    }
  }
    
  return(normDataTable)
}
