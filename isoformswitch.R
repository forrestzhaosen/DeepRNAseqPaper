if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
BiocManager::install(version='devel')
BiocManager::install("IsoformSwitchAnalyzeR")


library(IsoformSwitchAnalyzeR)
remotes::install_github("TobiTekath/DTUrtle")


salmonQuant <- importIsoformExpression(
  parentDir = system.file("extdata/",package="IsoformSwitchAnalyzeR"))

myDesign <- data.frame(
  sampleID = colnames(salmonQuant$abundance)[-1],
  condition = gsub('_.*', '', colnames(salmonQuant$abundance)[-1])
)

aSwitchList <- importRdata(
  isoformCountMatrix   = salmonQuant$counts,
  isoformRepExpression = salmonQuant$abundance,
  designMatrix         = myDesign,
  isoformExonAnnoation = system.file("extdata/example.gtf.gz"             , package="IsoformSwitchAnalyzeR"),
  isoformNtFasta       = system.file("extdata/example_isoform_nt.fasta.gz", package="IsoformSwitchAnalyzeR"),
  fixStringTieAnnotationProblem = TRUE,
  showProgress = FALSE
)
