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

#for compatibility with R local script, set every variable to a boolean depending on whether it exists
# only to be run when the code is called from the webportal or GenePattern

if(!exists("rawdataZip")) rawdataZip <- ""
if(!exists("refName")) refName <- ""
if(!exists("arrayGroup")) arrayGroup <- ""
reOrder <- exists("reOrder")
if(!exists("maxArray")) maxArray <- 41
layoutPlot <- exists("layoutPlot")
controlPlot <- exists("controlPlot")
samplePrep <- exists("samplePrep")
ratio <- exists("ratio")
degPlot <- exists("degPlot")
hybrid <- exists("hybrid")
percPres <- exists("percPres")
posnegDistrib <- exists("posnegDistrib")
bgPlot <- exists("bgPlot")
scaleFact <- exists("scaleFact")
boxplotRaw <- exists("boxplotRaw")
boxplotNorm <- exists("boxplotNorm")
densityRaw <- exists("densityRaw")
densityNorm <- exists("densityNorm")
MARaw <- exists("MARaw")
MANorm <- exists("MANorm")
if(!exists("MAOption1")) MAOption1 <- ""
spatialImage <- exists("spatialImage")
PLMimage <- exists("PLMimage")
posnegCOI <- exists("posnegCOI")
Nuse <- exists("Nuse")
Rle <- exists("Rle")
correlRaw <- exists("correlRaw")
correlNorm <- exists("correlNorm")
clusterRaw <- exists("clusterRaw")
clusterNorm <- exists("clusterNorm")
if(!exists("clusterOption1")) clusterOption1 <- ""
if(!exists("clusterOption2")) clusterOption2 <- ""
PCARaw <- exists("PCARaw")
PCANorm <- exists("PCANorm")        
PMAcalls <- exists("PMAcalls")
if(!exists("normMeth")) normMeth <- ""
if(!exists("normOption1")) normOption1 <- ""
customCDF <- exists("customCDF")
if(!exists("CDFtype")) CDFtype <- ""
if(!exists("species")) species <- ""

print ("Parameters have been registered")
