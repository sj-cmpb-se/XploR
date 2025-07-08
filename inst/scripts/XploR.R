#!/usr/bin/env Rscript
library(optparse)
library(XploR)

# ... (define option_list and parse args as before) ...

SummarizeMapQC(input_dir, prefix, output_dir)

runAIsegmentation(
  seg = opt$seg,
  cov = opt$cov,
  ai = opt$ai,
  gender = opt$gender,
  out_dir = opt$out_dir,
  prefix = opt$prefix,
  aibinsize = opt$aibinsize,
  mergeai = opt$mergeai,
  mergecov = opt$mergecov,
  minaisize = opt$minaisize,
  snpmin = opt$snpmin,
  minsnpcallaicutoff = opt$minsnpcallaicutoff,
  mergecovminsize = opt$mergecovminsize
)


runModelLikelihood(
    seg = opt$seg,
    out_dir = opt$out_dir,
    prefix = opt$prefix,
    gender = opt$gender,
    lambda = 1,
    gamma = 1,
    epsilon = 0.01,
    modelminprobes = 20,
    modelminAIsize = 5000000,
    minsf = 0.4,
    callcov = 0.3,
    thread = 4,
    diploidweight = 0.5
)


AnnotateSegments(
  input,
  out_dir,
  prefix,
  cytoband,
  whitelist_edge,
  gene)


runPlotCNV(
    seg, coverage, ballele,
    ai_binsize = 100000, cov_binsize = 100000,
    whitelist, gender = "female",
    out_dir, prefix
)

CoverageQC(annofile)


rerunCalling(
    seg, input, models, call, gender,
    dicovsf = NULL, purity = NULL,
    chromosome = NULL, start = NULL, end = NULL,
    mode = "model",
    out_file
)
