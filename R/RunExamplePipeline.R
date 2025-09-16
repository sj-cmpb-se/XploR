#' Run XploR Example Pipeline
#'
#' Runs the XploR pipeline on bundled example data and writes results to a temporary directory (or user-specified directory).
#' @param out_dir Output directory for results. Defaults to a temporary directory.
#' @return Invisibly returns the output directory path.
#' @export
RunExamplePipeline <- function(out_dir) {
  # Locate example data
  seg <- system.file("extdata", "example_female.called.igv.seg", package = "XploR")
  cov <- system.file("extdata", "example_female.counts", package = "XploR")
  ai <- system.file("extdata", "example_female_dragen.tumor.ballele.counts.gz", package = "XploR")
  cr <- system.file("extdata", "example_female.denoisedCR.tsv", package = "XploR")
  pon <- system.file("extdata", "PON_AI.Rdata", package = "XploR")
  cytoband <- system.file("extdata", "hg19_cytoBand.dat", package = "XploR")
  whitelist_bed <- system.file("extdata", "hg19_gatk_female_pon_whitelist.bed", package = "XploR")
  blacklist_bed <- system.file("extdata", "hg19_gatk_female_pon_blacklist.bed", package = "XploR")
  whitelist_edge <- system.file("extdata","hg19_detectable_edges.txt",package = "XploR")
  gene <- system.file("extdata", "RefSeqCurated.genePred.gene_region.txt", package = "XploR")
  ai_pon <- system.file("extdata", "PON_AI.Rdata", package = "XploR")
  gender <- "female"
  prefix <- "example"


# Segmentation
print("Run segmentation.......")
segmentation <- tryCatch({
  RunAIsegmentation(
    seg = seg,
    cov = cov,
    ai = ai,
    gender = gender,
    out_dir = out_dir,
    prefix = prefix,
    ai_pon = ai_pon,
    aitype = "dragen"
  )
  TRUE
}
, error = function(e) {
  message("ERROR during segmentation: ", conditionMessage(e))
  quit(save = "no", status = 1, runLast = FALSE) # Exits the workflow with error status
})

if(segmentation){
  print("Segmentation is DONE.")
}

## Modeling
print("Run model likelihood calculation and selection...............")
modeling <- tryCatch({
  RunModelLikelihood(
    seg = paste0(out_dir,"/",prefix,"_GATK_AI_segment.tsv"),
    out_dir = out_dir,
    prefix = prefix,
    gender = gender,
    modelminprobes = 20,
    modelminAIsize = 5000000,
    minsf = 0.4,
    callcov = 0.3,
    thread = 6)
  TRUE

}, error = function(e) {
  message("ERROR during likelihood calculation: ", conditionMessage(e))
  quit(save = "no", status = 1, runLast = FALSE) # Exits the workflow with error status
})

if(modeling){
  print("Estimate model is DONE.")
}

print("Run annotation segments................")
anno <- tryCatch({
  AnnotateSegments(
    input = paste0(out_dir,"/",prefix,"_final_calls.tsv"),
    out_dir = out_dir,
    prefix = prefix,
    cytoband = cytoband,
    whitelist_edge = whitelist_edge,
    gene = gene)
  TRUE
}, error = function(e) {
  message("ERROR during segment annotation: ", conditionMessage(e))
  quit(save = "no", status = 1, runLast = FALSE) # Exits the workflow with error status
})

if (anno) {
  print("CNV segment annotation DONE.")
}

print("Generating CNV plot.............")

plot_result <- tryCatch({
  RunPlotCNV(
    seg = paste0(out_dir,"/",prefix,"_CNV_annotation.tsv"),
    cr = cr,
    ballele = ai,
    ai_binsize = 100000,
    cov_binsize = 100000,
    whitelist = whitelist_bed,
    gender = gender,
    out_dir = out_dir,
    prefix = prefix,
    aitype = "dragen"
  )
  TRUE  # Return TRUE if successful
}, error = function(e) {
  message("ERROR during CNV plot generation: ", conditionMessage(e))
  quit(save = "no", status = 1, runLast = FALSE) # Exits the workflow with error status
})

if (plot_result) {
  print("CNV plot generation DONE.")
}

print("Generating AI segment quality file............")
bafqc <- tryCatch({
  bafqc_result <- BafQC(
    annofile = paste0(out_dir,"/",prefix,"_CNV_annotation.tsv"),
    out_dir = out_dir,
    prefix = prefix)
  TRUE
},error = function(e) {
  message("ERROR during generate quality control table: ", conditionMessage(e))
  quit(save = "no", status = 1, runLast = FALSE) # Exits the workflow with error status
})

if (bafqc) {
  print("Segment QC is DONE.")
}

#RerunCNV(
#  seg = paste0(out_dir,"/",prefix,"_GATK_AI_segment.tsv"),
#  input = paste0(out_dir,"/",prefix,"_top_likelihood_calls.tsv"),
#  models = "paste0(out_dir,"/",prefix,"_Models_likelihood.tsv"),
#  call = "paste0(out_dir,"/",prefix,"_final_calls.tsv"),
#  gender = gender,
#  mode = "region",
#  chromosome = "5",
#  start = "63986201",
#  out_file = out_file,
#  callcov = 0.3
#)
if( segmentation && modeling && anno && plot_result && bafqc ){
  print(paste("Example pipeline run complete. Results are in:", out_dir))
}
invisible(out_dir)
}
