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
# Get input parameters                          				              #
###############################################################################

affyAnalysisQC<-function(...) {


options(warn = 1, keep.source = TRUE, error = 
  quote({ 
    cat("Environment:\n", file=stderr()); 

    # TODO: setup option for dumping to a file (?)
    # Set `to.file` argument to write this to a file for post-mortem debugging    
    dump.frames();  # writes to last.dump

    #
    # Debugging in R
    #   http://www.stats.uwo.ca/faculty/murdoch/software/debuggingR/index.shtml
    #
    # Post-mortem debugging
    #   http://www.stats.uwo.ca/faculty/murdoch/software/debuggingR/pmd.shtml
    #
    # Relation functions:
    #   dump.frames
    #   recover
    # >>limitedLabels  (formatting of the dump with source/line numbers)
    #   sys.frame (and associated)
    #   traceback
    #   geterrmessage
    #
    # Output based on the debugger function definition.

    n <- length(last.dump)
    calls <- names(last.dump)
    cat(paste("  ", 1L:n, ": ", calls, sep = ""), sep = "\n", file=stderr())
    cat("\n", file=stderr())

    if (!interactive()) {
      q()
    }
  }))



#tryCatch({

# GenePattern command line of this type: 
#          <R2.5> <libdir>myscript.R myfunction -i<input.file> 
# Input parameters formatted for GenePattern and arrayanalysis.org
# Default values managed by GP and web input forms
  
	args <- list(...)
	for(i in 1:length(args)) {
		flag <- substring(args[[i]], 0, 2)
		value <- substring(args[[i]], 3, nchar(args[[i]]))
		if(flag=='-d'){
			rawdataZip = value
		}	
		if(flag=='-A'){
			refName = value
		}				
		if(flag=='-g'){
			arrayGroup = value
		}				
		if(flag=='-m'){
			maxArray = value
		}
		if(flag=='-G'){
			reOrder = value
		}
		if(flag=='-s'){
			samplePrep = value
		}
		if(flag=='-r'){
			ratio = value
		}
		if(flag=='-e'){
			degPlot = value
		}
		if(flag=='-h'){
			hybrid = value
		}
		if(flag=='-b'){
			bgPlot = value
		}
		if(flag=='-p'){
			percPres = value
		}		  
		if(flag=='-P'){
			PMAcalls = value
		}
		if(flag=='-n'){
			posnegDistrib = value
		}
		if(flag=='-H'){
			controlPlot = value
		}
		if(flag=='-f'){
			scaleFact = value
		}
		if(flag=='-x'){
			boxplotRaw = value
		}
		if(flag=='-X'){
			boxplotNorm = value
		}
		if(flag=='-y'){
			densityRaw = value
		}
		if(flag=='-Y'){
			densityNorm = value
		}
		if(flag=='-k'){
			MARaw = value
		}
		if(flag=='-K'){
			MANorm = value
		}
		if(flag=='-j'){
			MAOption1 = value
		}
		if(flag=='-F'){
			layoutPlot = value
		}		
		if(flag=='-N'){
			posnegCOI = value
		}
		if(flag=='-R'){
			spatialImage = value
		}
		if(flag=='-W'){
			PLMimage = value
		}
		if(flag=='-u'){
			Nuse = value
		}
		if(flag=='-a'){
			Rle = value
		}
		if(flag=='-c'){
			correlRaw = value
		}
		if(flag=='-C'){
			correlNorm = value
		}
		if(flag=='-t'){
			PCARaw = value  
		}
		if(flag=='-T'){
			PCANorm = value  
		}
		if(flag=='-o'){
			clusterRaw = value
		}
		if(flag=='-O'){
			clusterNorm = value
		}
		if(flag=='-v'){
			clusterOption1 = value
		}
		if(flag=='-w'){
			clusterOption2 = value
		}
		if(flag=='-z'){
			normMeth = value
		}     
		if(flag=='-J'){
			normOption1 = value
		} 
		if(flag=='-l'){
			customCDF = value
		}                                            
		if(flag=='-L'){
			CDFtype = value
		}                                            
		if(flag=='-S'){
			species = value
		} 
	}


  source("setParametersQC_web.R",local=TRUE)
  source("processFilesQC_web.R",local=TRUE)
  source("run_affyAnalysisQC.R",local=TRUE)


#}, error = function(e){ traceback()  });
   
# DESCRIPTION : see affyAnalysisQC.R
}
