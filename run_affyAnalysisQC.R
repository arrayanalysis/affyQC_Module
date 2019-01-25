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

version_nb <- "1.0.0"
cat("Script run using R version ",R.Version()$major,".",R.Version()$minor,
  " and affyAnalysisQC version_",version_nb,"\n",sep="")

#set memory to maximum on windows 32bit machines
if(length(grep("w32",R.Version()$os,fixed=TRUE))>0) memory.size(4095)

###############################################################################
# Load R libraries and affyAnalysisQC functions  							  #
###############################################################################
require("affy", quietly = TRUE)
require("affycomp", quietly = TRUE)
require("affyPLM", quietly = TRUE)
require("affypdnn", quietly = TRUE)
require("bioDist", quietly = TRUE)
require("simpleaffy", quietly = TRUE)
require("affyQCReport", quietly = TRUE)
require("plier", quietly = TRUE)
if(exists("samplePrep")) require("yaqcaffy", quietly = TRUE)
require("gdata", quietly = TRUE) #trim function
require("gplots", quietly = TRUE) #heatmap.2 functions
print("Libraries have been loaded")

reload <- function() {
source(paste(SCRIPT.DIR,"functions_processingQC.R",sep=""))
source(paste(SCRIPT.DIR,"functions_imagesQC.R",sep=""))
print ("Functions have been loaded")
}
reload();

if(!exists("rawData")) { # this is the case when the script is run locally
  setwd(DATA.DIR)

  ##-- check if is oligo, set global variable isOligo (turn this into function?!?!) (remove global variable!!!!!!!!!)
  ##copied, inserted, and modified (as indicated) from processFilesQC_web.R - temp fix

  #print(DATA.DIR)
  filenameindex <- grep("[.CEL]$", dir(), fixed=FALSE, ignore.case=TRUE)[1] #modified line
  
  print(dir()[filenameindex])

  if(!is.na(filenameindex)){
     #print("going to check oligo")     
     filename <- dir()[filenameindex];
     res <- tryCatch( ReadAffy(filenames=filename) , warning = function(e){e$message}, error = function(e){e$message})
     print(res)
     if(class(res)=="AffyBatch")
	isOligo <- FALSE
     else
        isOligo <- grepl("oligo", res, ignore.case=TRUE)
     #print(isOligo)
  } else {
     stop("No .cel file found in the specified directory");
  }

  if(isOligo){
	  #source("http://bioconductor.org/biocLite.R")
	  require("oligo") #asklars
          print("reading affy with oligo")
	  #print(regexpr("[\\.CEL^]", dir(), fixed=FALSE, ignore.case=TRUE))
          #rawData <- read.celfiles( regexpr("[\\.CEL^]", dir(), fixed=FALSE, ignore.case=TRUE)  )
	  celFiles <- list.celfiles()
	  print(celFiles)
          rawData <- read.celfiles( grep("[\\.CEL^]$", dir(), fixed=FALSE, ignore.case=TRUE,value=TRUE)  )#modified line
  } else {
    if(exists("prefOligo")) {
      if(prefOligo){
        require("oligo") #asklars
        celFiles <- list.celfiles()
        print(celFiles)
        try(rawData <- read.celfiles( grep("[\\.CEL^]$", dir(), fixed=FALSE, ignore.case=TRUE,value=TRUE)  ),TRUE) #modified line
		if(exists("rawData")) isOligo <- TRUE
      }
    }
    if(!exists("rawData")) {
	    print("reading affy with regular lib")
      rawData <- ReadAffy()
    } else {
      print("data read with oligo")
    }
  }
  ##--
  
  print("Raw data have been loaded in R")
}

if(!exists("libdir")) { # libdir exists only for GenePattern usage
  setwd(SCRIPT.DIR)
  setwd(WORK.DIR)
}

# Make sure that the CDF environment works
if(!isOligo){
	rawData <- addStandardCDFenv(rawData)   # if already works, won't be changed
	# Verify the array type (PMMM or PMonly)
	aType <- getArrayType(rawData)
} else {
	aType <- "PMonly"
}

# When refName does not exist, use the empty string
if(!exists("refName")) refName <- ""

###############################################################################
# Create array groups and array names                                         #
###############################################################################

if(arrayGroup!=""){
  # Information is available: groups will be created
  # 1- read the arrayGroup file and trim spaces
  # 2- define the array names and classes (experimentFactor)
  if(!exists("DESC.DIR")) DESC.DIR <- ""
  
  descfile <- paste(DESC.DIR, arrayGroup, sep="")
  extension<-strsplit(descfile,"\\.")
  extension<-paste(".",extension[[1]][length(extension[[1]])],sep="")  
  description = NULL;
  switch(extension,
         ".txt" = description<-as.data.frame(apply(read.delim(descfile, fill = FALSE, as.is=TRUE),2,trimws)),
         ".csv" = description<-as.data.frame(apply(read.csv(descfile, fill = FALSE, as.is=TRUE),2,trimws)),
         ".xls" = {library(gdata); description<-as.data.frame(apply(read.xls(descfile, as.is=TRUE),2,trimws))},
         ".xlsx" = {library(gdata); description<-as.data.frame(apply(read.xls(descfile, as.is=TRUE),2,trimws))}
	)
  if(is.null(description)) stop(paste("extension",extension,"not recognised"))
  
 # description <- trim(read.table(paste(DESC.DIR, arrayGroup, sep=""), 
	 # 	  header = TRUE, as.is = TRUE, sep="\t"))
	
  if(length(grep(".CEL",toupper(colnames(description)[1]),
    ignore.case = TRUE))>0) {
    stop(paste("The description file may not contain a header, as the first",
     	"column header seems to be a CEL file name"))
  }
  file_order <- match(description[,1],sampleNames(rawData))
  if(sum(is.na(file_order)) > 0) stop("file names in data directory and file names in description file do not match")
  if(length(unique(file_order)) < length(file_order)) stop("file names in description file are not unique")
  rawData <- rawData[,file_order]

  sampleNames(rawData)<- as.character(description[,2])
  experimentFactor <- factor(description[,3])
  
  # if required reorder the arrays according to group levels in order to keep 
  # groups together in all plots
  if(reOrder) {
    rawData <- rawData[,order(experimentFactor)]
    experimentFactor <- experimentFactor[order(experimentFactor)]
  }
} else {
  # No information: arrays will be computed/colored independently  
  sampleNames(rawData) <- as.character(sampleNames(rawData))
  experimentFactor <- factor(rep(1, length(sampleNames(rawData))))
  description <- cbind(sampleNames(rawData),sampleNames(rawData),
    experimentFactor)
  colnames(description) <- c("ArrayDataFile","SourceName","FactorValue")
}

# Create colorset for the array groups
#-------------------------------------
colList <- colorsByFactor(experimentFactor)
plotColors <- colList$plotColors
legendColors <- colList$legendColors
rm(colList)

# Create symbolset for the array groups
#--------------------------------------
plotSymbols <- 18-as.numeric(experimentFactor)
legendSymbols <- sort(plotSymbols, decreasing=TRUE)

###############################################################################
# Define display parameters for the images			                          #
###############################################################################

WIDTH <- 1000
HEIGHT <- 1414
POINTSIZE <- 24
if(!exists("maxArray")) maxArray <- 41

###############################################################################
# Calculate the indicator values and begin the report                         #
###############################################################################

#create a cover sheet for the report to be created later
#and create a page indicating the naming and grouping used
coverAndKeyPlot(description, refName,WIDTH=WIDTH,HEIGHT=HEIGHT)

#create a table with several QC indicators
if(samplePrep || ratio || hybrid || percPres || bgPlot || scaleFact) {

  # The indicators are calculated only for PM-MM arrays as the calculation
  # based on MAS5 does not work for PM-only arrays

  quality <- NULL
  try(quality <- qc(rawData),TRUE) # calculate Affymetrix quality data for PMMM
  if(is.null(quality)) {
    warning("Plots based on the simpleaffy qc function cannot be created for this chip type")
  }
  
  if(samplePrep) {    
    # find the data 
    try(yack <- yaqc(rawData),TRUE)
    if(exists("yack")) {
      spnames<-rownames(yack@morespikes[grep("(lys|phe|thr|dap).*3", # only 3' 
      rownames(yack@morespikes), ignore.case = TRUE),])
      sprep<-t(yack@morespikes[spnames,])
    } else {
      sprep <- NULL
      warning("Plots based on the yaqc function cannot be created for this chip type")
    }    
    
    try({calls<-detection.p.val(rawData)$call
    lys<-calls[rownames(calls)[grep("lys.*3",rownames(calls),ignore.case=TRUE)],]
    rm(calls)},TRUE)
    if(!exists("lys")) {
      lys <- NULL
      warning("Plots based on the detection.p.val function cannot be created for this chip type")
    }else{
		if(length(lys) > length(sampleNames(rawData))) { lys<-lys[1,] }
    }
  }
  
  QCtablePlot(rawData,quality,sprep,lys,samplePrep=samplePrep,ratio=ratio,
      hybrid=hybrid,percPres=percPres,bgPlot=bgPlot,scaleFact=scaleFact,
	  WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE)
}


###############################################################################
# Raw data Quality Control graphs                                             #
###############################################################################

print("Graphs ready to be computed")

# 1.1 Sample prep controls
#-------------------------

if(samplePrep && !is.null(sprep) && !is.null(lys)) {
  print ("   plot sample prep controls"  )
  samplePrepPlot(rawData,sprep,lys,plotColors,
     WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=maxArray)
}
 
# 1.2 3'/5' ratio - only for PM-MM arrays
#----------------------------------------

if(ratio && !is.null(quality)) {
  print ("   plot beta-actin & GAPDH 3'/5' ratio")
  ratioPlot(rawData,quality=quality,experimentFactor,plotColors,legendColors,
     WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=maxArray)
}

# 1.3 RNA degradation plot
#-------------------------

if(degPlot) {
  print ("   plot degradation plot (skipped for oligo based analysis)"  )
  try(RNAdegPlot(rawData,plotColors=plotColors,
     WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=maxArray),TRUE)
}

###############################################################################
# 2.1 Spike-in controls - only for PM-MM arrays
#----------------------------------------------

if(hybrid && !is.null(quality)) {
  print ("   plot spike-in hybridization controls"  )
  hybridPlot(rawData,quality=quality,plotColors,
     WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=maxArray)
}

# 2.2 Background intensities - only for PM-MM arrays
#---------------------------------------------------

if(bgPlot && !is.null(quality)) {
  print ("   plot background intensities"  )
  backgroundPlot(rawData,quality=quality,experimentFactor,plotColors,legendColors,
     WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=maxArray)
}

# 2.3 Percent present - only for PM-MM arrays
#---------------------------------------------

if(percPres && !is.null(quality)) {
  print ("   plot percent present"  )
  percPresPlot(rawData,quality=quality,experimentFactor,plotColors,legendColors,
     WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=maxArray)
}

# 2.4 Table of PMA-calls based on the MAS5 algorithm - only for PM-MM arrays
#---------------------------------------------------------------------------

if(PMAcalls) {
  if(customCDF) {
    if(species=="") {
      warning("Species has not been set and custom cdf requested, attempting to deduce species for chip type")
      species <- deduceSpecies(rawData@annotation)
    }
	if(species!=""){
		PMAtable <- computePMAtable(rawData,customCDF,species,CDFtype)
	}else{
		warning("Could not define species; the CDF will not be changed")
		PMAtable <- computePMAtable(rawData,customCDF)
	}
  } else {
    PMAtable <- computePMAtable(rawData,customCDF)
  }
  if(!is.null(PMAtable)) {
    write.table(PMAtable, "PMAtable.txt", sep="\t", row.names=FALSE, 
	  col.names=TRUE, quote=FALSE)
  }
}

# 2.5 Pos and Neg control distribution
#-------------------------------------

if(posnegDistrib) {
  print ("   plot pos & neg control distribution  (skipped for oligo based analysis)"  )
  try(PNdistrPlot(rawData,
     WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE),TRUE)
}

# 2.6 affx control profiles and boxplot
#--------------------------------------

if(controlPlot) {
  print ("   plot control profiles and/or boxplots (skipped for oligo based analysis)")
  try(controlPlots(rawData,plotColors,experimentFactor,legendColors,
     WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=maxArray),TRUE)
}

###############################################################################
# 3.1.1 Scale factor - only for PM-MM arrays
#-------------------------------------------

if(scaleFact && !is.null(quality)) {
  print ("   plot scale factors")
  scaleFactPlot(rawData,quality=quality,experimentFactor,plotColors,
     legendColors,WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,
	 MAXARRAY=maxArray)
}

# 3.1.2 Boxplot of raw log-intensities
#-------------------------------------

if(boxplotRaw){
  print ("   plot boxplot for raw intensities")
  boxplotFun(Data=rawData, experimentFactor, plotColors, legendColors,
     WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=maxArray)
}

# 3.1.3 Density histogram of raw log-intensities
#-----------------------------------------------

if(densityRaw){
  print ("   plot density histogram for raw intensities")
  densityFun(Data=rawData, plotColors,
     WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=maxArray)
  #densityFunUnsmoothed(Data=rawData, plotColors,
  #  WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=maxArray)
}

# 3.2.1 MA-plot or raw data
#--------------------------


if(MARaw){
  print ("   MA-plots for raw intensities")
  maFun(Data=rawData, experimentFactor, perGroup=(MAOption1=="group"), 
     aType=aType,WIDTH=WIDTH,HEIGHT=HEIGHT,MAXARRAY=maxArray)
}

# 3.3.1 Plot of the array layout
#-------------------------------

if(layoutPlot) {
  print ("   plot array reference layout")
  plotArrayLayout(rawData,aType,
     WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE)
}

# 3.3.2 Pos and Neg control Position
#-----------------------------------

if(posnegCOI && !isOligo){  
  print ("   Pos/Neg COI")
  PNposPlot(rawData,WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE)
}

if(!isOligo){ #start of !isOligo

# 3.3.3.1 Create PLM object
#--------------------------

# fit a probe level model on the raw data, used by nuse and rle plot as well
  rawData.pset <- NULL
  if(spatialImage || PLMimage || Nuse || Rle) {
  print ("   Fit a probe level model (PLM) on the raw data")  
    rawData.pset <- fitPLM(rawData)                     
  }


# 3.3.3.2 Spatial images
#---------------------

if(spatialImage) {   
  print ("   2D virtual images")
  valtry<-try(spatialImages(rawData, Data.pset=rawData.pset, TRUE,FALSE,FALSE,FALSE, 
	          WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE),
		      silent=TRUE)
  if(class(valtry)=="try-error") {
	print("      Use array.image instead of spatialImages function")
	if(length(sampleNames(rawData))>6){
		# Usage of a median array is interesting when there are enough arrays
		array.image(rawData,WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE)
	}else{
		# Usage when few arrays in dataset (one page for 3 arrays -> max: 2 pages)
		array.image(rawData,relative=FALSE,col.mod=4,symm=TRUE,WIDTH=WIDTH,
		  HEIGHT=HEIGHT,POINTSIZE=POINTSIZE)
	}
  }
}

# 3.3.3.3 PLM images
#---------------------

if(PLMimage) {  
  print ("   Complete set of 2D PLM images")
  valtry<-try(spatialImages(rawData, Data.pset=rawData.pset, TRUE, TRUE, TRUE, TRUE,
	            WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=maxArray), 
				silent=TRUE)
  if(class(valtry)=="try-error") {
	print("      Could not create the PLM images.")
  }
}


# 3.4.1 NUSE
#-----------

if(Nuse){
  print ("   NUSE boxplot")
  nuseFun(rawData, Data.pset=rawData.pset, experimentFactor, plotColors, 
     legendColors,WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,
	 MAXARRAY=maxArray)
}

# 3.4.2 RLE
#----------

if(Rle){          
  print ("   RLE boxplot")
  rleFun(rawData, Data.pset=rawData.pset, experimentFactor, plotColors, 
     legendColors,WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,
	 MAXARRAY=maxArray)
}

} #end of isOligo

###############################################################################
# 4.1 Correlation Plot  of raw data
#----------------------------------

if(correlRaw){
  print ("   Correlation plot of raw data")
  correlFun(Data=rawData, experimentFactor=experimentFactor, legendColors=legendColors,
     WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=maxArray)  
}

# 4.2 PCA analysis of raw data
#-----------------------------

if(PCARaw){
  print("   PCA analysis of raw data")
  pcaFun(Data=rawData, experimentFactor=experimentFactor, 
	plotColors=plotColors, legendColors=legendColors, plotSymbols=plotSymbols,
	legendSymbols=legendSymbols, namesInPlot=((max(nchar(sampleNames(rawData)))<=10)&&
	(length(sampleNames(rawData))<=(maxArray/2))),WIDTH=WIDTH,HEIGHT=HEIGHT,
	POINTSIZE=POINTSIZE)
}

# 4.3 Hierarchical Clustering of raw data
#-----------------------------------------

if(clusterRaw){
  print ("   Hierarchical clustering of raw data") 
  clusterFun(Data=rawData, experimentFactor=experimentFactor,
   clusterOption1=clusterOption1, clusterOption2=clusterOption2,
   plotColors=plotColors, legendColors=legendColors,
   plotSymbols=plotSymbols, legendSymbols=legendSymbols,
   WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=maxArray) 
}

###############################################################################
# Pre-processing                                                              #
###############################################################################

if (aType == "PMonly") {
  if (normMeth == "MAS5") {
    warning("MAS5 cannot be applied to PMonly arrays. Changed MAS5 to PLIER")
    normMeth <- "PLIER"
  }
#  if (normMeth == "GCRMA") {
#    warning("GCRMA cannot be applied to PMonly arrays. Changed GCRMA to RMA")
#    normMeth <- "RMA"
#  }  
}

if(normMeth!="" && normMeth!="none") {
  if(customCDF) { 
    if(species=="") {
      warning("Species has not been set and custom cdf requested, attempting to deduce species for chip type")
      species <- deduceSpecies(rawData@annotation)
    }	
	if(species!=""){
		normData <- normalizeData(rawData,normMeth,perGroup=(normOption1=="group"), 
		  experimentFactor, aType=aType, customCDF, species, CDFtype, isOligo, WIDTH=WIDTH,
		  HEIGHT=HEIGHT)
	}else{
		warning("Could not define species; the CDF will not be changed")
		normData <- normalizeData(rawData,normMeth,perGroup=(normOption1=="group"), 
		  experimentFactor, aType=aType, customCDF, isOligo, WIDTH=WIDTH,HEIGHT=HEIGHT)
	}

  } else {
    normData <- normalizeData(rawData,normMeth,perGroup=(normOption1=="group"), 
	  experimentFactor, aType=aType, customCDF,WIDTH=WIDTH,HEIGHT=HEIGHT)
  }
}

if((boxplotNorm || densityNorm || MANorm || correlNorm || clusterNorm || 
    PCANorm) && ((normMeth=="") || (normMeth=="none"))) {
  warning("One or more QC plots of normalized data requested, but no normalization selected, plots will be omitted")
} else {

###############################################################################
# Evaluation of the pre-processing                                            #
###############################################################################

# 5.1 Make a Box-plot of the normalized data
#-------------------------------------------

  if(boxplotNorm){
    print ("   plot boxplot for normalized intensities") 
    boxplotFun(Data=normData, experimentFactor, plotColors, legendColors, 
	  normMeth=normMeth,WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,
	  MAXARRAY=maxArray)
  }

# 5.2 Make a Density histogram of the normalized data
#----------------------------------------------------

  if(densityNorm){
    print ("   plot density histogram for normalized intensities")
    densityFun(Data=normData, plotColors, normMeth=normMeth,
      WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=maxArray)
    #densityFunUnsmoothed(Data=normData, plotColors, normMeth=normMeth,
    #  WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=maxArray)
  }

# 5.3 Make separate MA-plots for each group on normalized data
#-------------------------------------------------------------	

  if(MANorm){
    print ("   MA-plots for normalized intensities") 
    maFun(Data=normData, experimentFactor, perGroup=(MAOption1=="group"), 
	 normMeth=normMeth,WIDTH=WIDTH,HEIGHT=HEIGHT,MAXARRAY=maxArray)
  }

# 5.4 Make correlation plots on normalized data
#----------------------------------------------

  if(correlNorm){
    print ("   Correlation plot of normalized data") 
    correlFun(Data=normData, normMeth=normMeth, experimentFactor=experimentFactor, legendColors=legendColors,
     WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=maxArray)
  }

# 5.5 PCA analysis of normalized data
# -----------------------------------

  if(PCANorm){
    print("   PCA graph for normalized data")
    pcaFun(Data=normData, experimentFactor=experimentFactor,normMeth=normMeth, 
	  plotColors=plotColors, legendColors=legendColors, plotSymbols=plotSymbols,
	  legendSymbols=legendSymbols, namesInPlot=((max(nchar(sampleNames(rawData)))<=10)&&
	  (length(sampleNames(rawData))<=(maxArray/2))),WIDTH=WIDTH,HEIGHT=HEIGHT,
	  POINTSIZE=POINTSIZE)
  }

# 5.6 Make hierarchical clustering on normalized data
#----------------------------------------------------

  if(clusterNorm){
    print ("   Hierarchical clustering of normalized data") 
    clusterFun(Data=normData, experimentFactor=experimentFactor,
    clusterOption1=clusterOption1, clusterOption2=clusterOption2,
    normMeth=normMeth, plotColors = plotColors, legendColors = legendColors,
    plotSymbols=plotSymbols, legendSymbols=legendSymbols,
    WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=maxArray)
  }
}

###############################################################################
# Prepare the output data                          				              #
###############################################################################

# Export the normalized data 

if((normMeth=="") || (normMeth=="none")) {
  warning("No normalization selected, normalized data table not saved")
} else {  
  print("Saving normalized data table")

  if(isOligo)
   normDataTable <- createNormDataTable(normData, customCDF=FALSE, species, CDFtype)
  else
   normDataTable <- createNormDataTable(normData, customCDF=(sum(featureNames(normData)!=featureNames(rawData)[1:length(featureNames(normData))])>0), species, CDFtype)
    
  #output normalised expression data to file
  refName <- sub("(_\\d{4}-\\d{2}-\\d{2}_\\d{2}-\\d{2}_\\d{2})", "", refName)  
  normFileName <- paste(normMeth,"NormData_",refName,".txt",sep="")
  print(paste("Normalized data table:", normFileName))
  write.table(normDataTable, normFileName, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
}

# clean R: (or quit without saving the environment...)
# rm(list = ls())

print("I am on test mode")
warning("test mode");
write("prints to stdout", stdout())


if(!is.null(warnings())) warnings()
