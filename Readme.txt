######***********************######
Part 0  BasicSpatialAnalysis
#The basic analysis method  for SPT data in Seurat.

##First time use,install packages (Seuart)
##Source functions and packages
##Adjust the working directory to the directory where the function script is placed
##Data read and save
##Graphing, a preliminary look at the data distribution
##Citing single-cell sequencing dimensionality reduction methods
##dimensionality reduction
##Look at the gene expression of subpopulations.
##Look at the gene expression of subgroups after grouping according to pathological information

######***********************######
Part 1 DataIntegration&Merge
#The two methods to integrate data.

######***********************######
Part 1 SPATAanalysis
#transform from Seurat To Spata, and analysis in SPATA.

######***********************######
Part 1 SpatialGetPercentPlotSplit
#The method to get split percent plot in SPT data.

######***********************######
Part 2 SpatialGenesetScore
#The method to get Geneset Score in SPT data.

######***********************######
Part 3 ContourAnalysis
#Extract gene exp and spot coordinates from each seurat object, to get Contour plot.

######***********************######
Part 3 CalculateGenesetExpression
#A function use global variable genes.with.exp to filter genesets from GenesetContourMaps.

######***********************######
Part 4 SpatialGetNearlyCell
#The distance between selected cells and all tumor cells were calculated.
#The shortest distance were used to classify cells surrounding the tumor foci.

######***********************######
Part 4 SpatialInfercnv
#SpatialInfercnv is used to explore tumor SPT data to identify evidence for large-scale chromosomal copy number variations.

######***********************######
Part 4 SpatialMonocle
#SpatialMonocle performs differential expression and time-series analysis for SPT experiments.

######***********************######
Part 4 SpatialNichenet
#SpatialNicheNet: modeling intercellular communication by linking ligands to target genes.

######***********************######
Part 5 CellLabelTransfer
#CellLabelTransfer is a method to find anchors between a reference and query object, and use the anchors to transfer data from the reference to query object

######***********************######
Part 5 MakePredictionDotPlot
#The method to get dotplot in cell label transfer data.

######***********************######
Part 6 Cibersort
#CIBERSORT analysis tool is a gene expression-based arithmetic, which uses a series of bar code genetic expression results (a ¡°signature matrix¡± of 547 genes) to assess data of immunocyte constituents from bulk cancer specimens.

######***********************######
Part 6 MCPcounter
#Estimating the population abundance of tissue-infiltrating immune and stromal cell populations using gene expression.