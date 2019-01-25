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
# Set directories 		                          				              #
###############################################################################

### Change these Paths if needed ###
DATA.DIR <- "D:/Microarrays/cBITE data/Total data set collection/StudyID_14 Piel nuclear confinement/raw-data"
SCRIPT.DIR <- "D:\\R-stuff\\ArrayAnalysis\\affyQC_Module oligo version"
# If you want to use the last version, uncomment the following line after 
# updating the milestone number:	
# SCRIPT.DIR <- "http://svn.bigcat.unimaas.nl/arrayanalysis/tags/version_1.0.0/src/"
WORK.DIR <- "D:/Microarrays/cBITE data/Total data set collection/StudyID_14 Piel nuclear confinement/raw-data/normalization"

####################################

correctDIR <- function(d) { 
  lastChar <- substr(d,nchar(d),nchar(d))
  if((lastChar != "/") && (lastChar != "/")) d <- paste(d,"/",sep="")
  return(d)
}
if(exists("DATA.DIR")) DATA.DIR <- correctDIR(DATA.DIR)
if(exists("SCRIPT.DIR")) SCRIPT.DIR <- correctDIR(SCRIPT.DIR)
if(exists("WORK.DIR")) WORK.DIR <- correctDIR(WORK.DIR)

###############################################################################
# Set input parameters                          				              #
###############################################################################

# Usage from R local install is to run the function in the repository 
# containing all CEL files (already unziped)

prefOligo <- TRUE
arrayGroup <- "D:\\Microarrays\\cBITE data\\Total data set collection\\StudyID_14 Piel nuclear confinement\\raw-data\\description.txt" #see comments below
reOrder <- TRUE
layoutPlot <- TRUE
controlPlot <- TRUE
samplePrep <- TRUE
ratio <- TRUE
degPlot <- TRUE
hybrid <- TRUE
percPres <- TRUE
posnegDistrib <- TRUE
bgPlot <- TRUE
scaleFact <- TRUE
boxplotRaw <- TRUE
boxplotNorm <- TRUE
densityRaw <- TRUE
densityNorm <- TRUE
MARaw <- TRUE
MANorm <- TRUE
MAOption1 <- "dataset" #see comments below
spatialImage <- TRUE
PLMimage <- FALSE
posnegCOI <- TRUE
Nuse <- TRUE
Rle <- TRUE
correlRaw <- TRUE
correlNorm <- TRUE
clusterRaw <- TRUE
clusterNorm <- TRUE
clusterOption1 <- "Spearman" #see comments below
clusterOption2 <- "ward.D2" #see comments below
PCARaw <- TRUE  
PCANorm <- TRUE     
PMAcalls <- FALSE
normMeth <- "RMA" #see comments below
normOption1 <- "dataset" #see comments below
customCDF <- FALSE
CDFtype <- "ENSG" #see comments below
species <- "" #see comments below
 
source(paste(SCRIPT.DIR,"run_affyAnalysisQC.R",sep=""),local=TRUE)

###############################################################################
# PARAMETER DESCRIPTION                         				              #
###############################################################################

# "automated calls" means usage from GenePattern and arrayanalysis.org
# Note that flag letters are only used in these automated calls
  
#     prefOligo = prefer Oligo for chiptypes for which both Affy and Oligo can be used
# d = rawdataZip = zip file with the .CEL data  (usage in automated calls only)
# D = libdir = global variable sent by GenePattern (usage in GenePattern call 
#        only)
# g = arrayGroup = description file describing the array names and experimental 
#        groups
# G = reOrder = boolean for whether the arrays have to be ordered per group in 
#        the plots FALSE keeps the order of the description file, TRUE reorders 
#		 per group
# m = maxArray = parameter adapting the image diplay if more than 'maxArray' 
#        arrays in the dataset 
# s = samplePrep = boolean for Sample prep controls 
# r = ratio = boolean for 3’/5’ for b-actin and GAPDH 
# e = degPlot = boolean for DNA degration plot 
# h = hybrid = boolean for Spike-in controls 
# b = bgPlot = boolean for Background intensity 
# p = percPres = boolean for Percent present 
# P = PMAcalls = boolean for Present/Marginal/Absent calls using MAS5
# n = posnegDistrib = boolean for +and - controls distribution
# H = controlPlot = boolean for plots of the AFFX controls on the arrays
# f = scaleFact = boolean for Scale factor 
# x = boxplotRaw = boolean for Raw boxplot of log-intensity 
# X = boxplotNorm = boolean for Norm boxplot of log-intensity 
# y = densityRaw = boolean for Raw density histrogram 
# Y = densityNorm = boolean for Norm density histrogram 
# k = MARaw = boolean for Raw MA-plot 
# K = MANorm = boolean for Norm MA-plot 
# j = MAOption1 = two possible values: "group" or "dataset"
# F = layoutPlot = boolean for plot of the array layout
# N = posnegCOI = boolean for + and – controls COI plot 
# R = spatialImage = boolean for 2D images
# W = PLMimage = boolean for 2D PLM plots 
# u = NUSE = boolean for NUSE
# a = RLE = boolean for RLE
# c = correlRaw = boolean for Raw correlation plot 
# C = correlNorm = boolean for Norm correlation plot 
# t = PCARaw = boolean for PCA analysis of raw data
# T = PCANorm = boolean for PCA analysis of normalized data    
# o = clusterRaw = boolean for Raw hierarchical clustering 
# O = clusterNorm = boolean for Norm hierarchical clustering 
# v = clusterOption1 = possibles values for Distance: "Spearman", "Pearson" or 
#        "Euclidian"
# w = clusterOption2 = possible values for Tree: "ward", "single",
#        "complete", "average", "mcquitty", "median" or centroid".    
# z = normMeth = poosible values for Data pre-processing: "RMA", "GCRMA", 
#        "PLIER", or "none"
# J = normOption1 = two possible values: "group" or "dataset"
# l = customCDF = boolean for a custom CDF for the pre-processed data 
# L = CDFtype = annotation format (default: ENSG), possibilities: "ENTREZG",
#        "REFSEQ","ENSG","ENSE","ENST","VEGAG","VEGAE","VEGAT",
#        "TAIRG","TAIRT","UG","MIRBASEF","MIRBASEG"
# S = species. It is required when customCDF is called. Possibilities: 
#        abbreviations: "Ag","At","Bt","Ce","Cf","Dr","Dm","Gg","Hs","MAmu",
#          "Mm","Os","Rn","Sc","Sp","Ss"
#        or full names: "Anopheles gambiae","Arabidopsis thaliana","Bos taurus",
#          "Caenorhabditis elegans","Canis familiaris", "Danio rerio",
#          "Drosophila melanogaster","Gallus gallus","Homo sapiens",
#          "Macaca mulatta","Mus musculus", "Oryza sativa","Rattus norvegicus",
#          "Saccharomyces cerevisiae","Schizosaccharomyces pombe","Sus scrofa"
# U = return image containing R data objects (not implemented yet)
