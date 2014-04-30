#=============================================================================#
# ArrayAnalysis                                                               #
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

# Called by QC module 
check_chiptype <- function(filename=NULL) {

  if(is.null(filename)) stop("data parameter not provided");
  
  library(affy)
  
  d <- ReadAffy(filenames=filename)
  
  source("functions_processingQC.R")
  
  aType <- getArrayType(d)
  
  if(!(is.null(d@annotation)) && (d@annotation!="")) {
    descr <- d@annotation
  } else {
    d <- addStandardCDFenv(d)
    if(!(is.null(d@annotation)) && (d@annotation!="")) {
      descr <- d@annotation
    } else {
      descr <- d@cdfName
    }
  }
  
  organism <- deduceSpecies(d@annotation)
  
  print(aType)
  
  print(descr)
  
  print(organism)
}

