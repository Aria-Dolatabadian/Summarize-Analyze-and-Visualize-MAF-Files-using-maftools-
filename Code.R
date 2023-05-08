#https://www.bioconductor.org/packages/devel/bioc/vignettes/maftools/inst/doc/maftools.html

library(maftools)

#path to TCGA LAML MAF file

laml.maf = ("tcga_laml.maf.gz")
laml.clin = ("tcga_laml_annot.tsv")
laml = read.maf(maf = laml.maf, clinicalData = laml.clin)

#Typing laml shows basic summary of MAF file.
laml

#Shows sample summry.
getSampleSummary(laml)
#Shows gene summary.
getGeneSummary(laml)
#shows clinical data associated with samples
getClinicalData(laml)
#Shows all fields in MAF
getFields(laml)
#Writes maf summary to an output file with basename laml.
write.mafSummary(maf = laml, basename = 'laml')

plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

#oncoplot for top ten mutated genes.
oncoplot(maf = laml, top = 10)

laml.titv = titv(maf = laml, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = laml.titv)

#lollipop plot for DNMT3A, which is one of the most frequent mutated gene in Leukemia.
lollipopPlot(
  maf = laml,
  gene = 'DNMT3A',
  AACol = 'Protein_Change',
  showMutationRate = TRUE,
  labelPos = 882
)

plotProtein(gene = "TP53", refSeqID = "NM_000546")

#Rainfall plots

brca <- system.file("extdata", "brca.maf.gz", package = "maftools")
brca = read.maf(maf = brca, verbose = FALSE)

rainfallPlot(maf = brca, detectChangePoints = TRUE, pointSize = 0.4)

laml.mutload = tcgaCompare(maf = laml, cohortName = 'Example-LAML', logscale = TRUE, capture_size = 50)

#Plotting VAF
plotVaf(maf = laml, vafCol = 'i_TumorVAF_WU')

