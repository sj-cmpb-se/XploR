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
  cytoband <- system.file("extdata", "hg19_cytoBand.dat", package = "XploR")
  whitelist_bed <- system.file("extdata", "hg19_gatk_female_pon_whitelist.bed", package = "XploR")
  blacklist_bed <- system.file("extdata", "hg19_gatk_female_pon_blacklist.bed", package = "XploR")
  whitelist_edge <- system.file("extdata","hg19_detectable_edges.txt",package = "XploR")
  gene <- system.file("extdata", "RefSeqCurated.genePred.gene_region.txt", package = "XploR")
  gender <- "female"
  prefix <- "example"
# generate reference files, the edge table is same in male and female only Y chrom is different, so in the reference file i only kept male version
#PonProcess(pon_file = pon_file,
#           blacklist_bed = blacklist_bed,
#           whitelist_bed = whitelist_bed,
#           cytoband = cytoband,
#           detectable_edge = detectable_edge,
#           gender = gender
#           )


# Segmentation
print("Run segmentation.......")
tryCatch({
  RunAIsegmentation(
    seg = seg,
    cov = cov,
    ai = ai,
    gender = gender,
    out_dir = out_dir,
    prefix = prefix,
    aibinsize = 500000,
    mergeai = 0.15,
    mergecov = 0.2,
    minaisize = 1000000,
    snpmin = 7,
    minsnpcallaicutoff = 10,
    mergecovminsize = 500000,
    aitype = "dragen"
  )
}
, error = function(e) {
  message("ERROR during segmentation: ", conditionMessage(e))
  quit(save = "no", status = 1, runLast = FALSE) # Exits the workflow with error status
})

## Modeling
print("Run model likelihood calculation and selection...............")
tryCatch({
  RunModelLikelihood(
    seg = paste0(out_dir,"/",prefix,"_GATK_AI_segment.tsv"),
    out_dir = out_dir,
    prefix = prefix,
    gender = gender,
    lambda = 1,
    gamma = 1,
    epsilon = 0.01,
    modelminprobes = 20,
    modelminAIsize = 5000000,
    minsf = 0.4,
    callcov = 0.3,
    thread = 6)

}, error = function(e) {
  message("ERROR during likelihood calculation: ", conditionMessage(e))
  quit(save = "no", status = 1, runLast = FALSE) # Exits the workflow with error status
})

print("Run annotation segments................")
tryCatch({
  AnnotateSegments(
    input = paste0(out_dir,"/",prefix,"_final_calls.tsv"),
    out_dir = out_dir,
    prefix = prefix,
    cytoband = cytoband,
    whitelist_edge = whitelist_edge,
    gene = gene)
}, error = function(e) {
  message("ERROR during segment annotation: ", conditionMessage(e))
  quit(save = "no", status = 1, runLast = FALSE) # Exits the workflow with error status
})

print("Generating CNV plot.............")
tryCatch({
  RunPlotCNV(
    seg = paste0(out_dir,"/",prefix,"_CNV_annotation.tsv"),
    cr =cr,
    ballele = ai,
    ai_binsize = 100000,
    cov_binsize = 100000,
    whitelist = whitelist_bed,
    gender = gender,
    out_dir = out_dir,
    prefix = prefix,
    aitype = "dragen"
  )
  }, error = function(e) {
    message("ERROR during generate CNV plot: ", conditionMessage(e))
    quit(save = "no", status = 1, runLast = FALSE) # Exits the workflow with error status
})

print("Generating AI segment quality file............")
tryCatch({
  BafQC(
    annofile = paste0(out_dir,"/",prefix,"_CNV_annotation.tsv"),
    out_dir = out_dir,
    prefix = prefix)
},error = function(e) {
  message("ERROR during generate quality control table: ", conditionMessage(e))
  quit(save = "no", status = 1, runLast = FALSE) # Exits the workflow with error status
})


#SummarizeMapQC(
#  input_dir = daragen_output,
#  prefix = prefix,
#  out_dir = out_dir
#)

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


message("Example pipeline run complete. Results are in: ", out_dir)
invisible(out_dir)
}
