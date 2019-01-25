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

###############################################################################
# Read and check array data and define repositories					          #
###############################################################################

# This script is called only from the WEBPORTAL 

###############################################################################
# Set directories 		                          				              #
###############################################################################

correctDIR <- function(d) { 
  lastChar <- substr(d,nchar(d),nchar(d))
  if((lastChar != "/") && (lastChar != "\\")) d <- paste(d,"/",sep="")
  return(d)
}

SCRIPT.DIR <- getwd()
if(exists("SCRIPT.DIR")) SCRIPT.DIR <- correctDIR(SCRIPT.DIR)

dirbase <- strsplit(SCRIPT.DIR,"/", fixed=TRUE)[[1]]
dirbase <- paste(dirbase[1:(length(dirbase)-1)],collapse="/")

WORK.DIR <- paste(dirbase,"/temp/",refName,"_QC", sep="") # directory where the results are computed (this assumes it relates to script.dir! To be adapted
if(exists("WORK.DIR")) WORK.DIR <- correctDIR(WORK.DIR)

if(!file.exists(WORK.DIR)) dir.create(WORK.DIR)

DATA.DIR <- paste(WORK.DIR,"DATA.DIR/",sep="") # CEL files are loaded using 'ReadAffy' then the repository is deleted
dir.create(DATA.DIR) # DATA.DIR is a relative path name
if(exists("DATA.DIR")) DATA.DIR <- correctDIR(DATA.DIR)

###############################################################################
# Unzip CEL files		                          				              #
###############################################################################

setwd(WORK.DIR)
print(paste("Zip file: ",rawdataZip))

#print(dir(WORK.DIR))

if(!file.exists(rawdataZip)) {
  stop("Execution halted! CANNOT FIND THE INPUT ZIP FILE CONTAINING THE DATA!") 
}

unzip(rawdataZip, exdir = DATA.DIR) 

###############################################################################
# Copy CEL files in DATA.DIR and remove files from WORK.DIR					  # 
###############################################################################
# The content of the Zip file is verified:
setwd(DATA.DIR)
# Case 1: the .CEL files are in WORK.DIR : we move them in the DATA.DIR file
celfile <- list.files(pattern = ".CEL",ignore.case = TRUE)
if(length(celfile) > 0 ){
#	listfile <- list.files()
#	file.copy(listfile, DATA.DIR)
#	unlink(listfile)    
#	setwd(DATA.DIR)
}else{  
# Case 2: no .CEL files were found
	stop("Execution halted! CANNOT FIND CEL FILES IN THE INPUT ZIP FILE!") 
}

print("Raw data ready to be loaded in R")

###############################################################################
# Load array data in R                          						      #
###############################################################################
#load the data
require("affy", quietly = TRUE)
#print(paste("R version:",R.Version()$version.string))
#print(paste("affy version:",sessionInfo()$otherPkgs$affy$Version))

  ##-- check if is oligo, set global variable isOligo (turn this into function?!?!) (remove global variable!!!!!!!!!)

  #print(DATA.DIR)
  filenameindex <- grep("[.CEL$]", dir(), fixed=FALSE, ignore.case=TRUE)[1]
  
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
          rawData <- read.celfiles( celFiles ) # this should be improved to get only the cell files
  }else{
	  print("reading affy with regular lib")
          rawData <- ReadAffy()
  }
  ##--

###############################################################################
# Clean space if the usage is the webportal        						      #
###############################################################################
#setwd(SCRIPT.DIR)
setwd(WORK.DIR)
unlink(DATA.DIR, recursive = TRUE)
rm(DATA.DIR,rawdataZip) #,listfile) 

###############################################################################
print("Raw data have been loaded in R")
setwd(SCRIPT.DIR)
