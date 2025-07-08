#' Run Model Likelihood and Copy Number Calling Workflow
#'
#' This function runs the complete workflow for likelihood calculation and model selection.
#'
#' @param seg Path to the combined segment file.
#' @param out_dir Output directory.
#' @param prefix Output file prefix.
#' @param gender Sample gender ("male" or "female").
#' @param lambda, gamma, epsilon, modelminprobes, modelminAIsize, minsf, callcov, thread, diploidweight, ratio Numeric parameters for the workflow (see script for details).
#' @param ... Additional arguments (reserved for future use).
#' @return Invisibly returns the output file paths.
#' @importFrom data.table fread
#' @importFrom dplyr filter mutate arrange rowwise ungroup select group_by summarise desc left_join
#' @importFrom tidyr replace_na
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster clusterEvalQ clusterExport
#' @export
runModelLikelihood <- function(
    seg, out_dir, prefix, gender,
    lambda = 1, gamma = 1, epsilon = 0.01, modelminprobes = 20,
    modelminAIsize = 5000000, minsf = 0.4, callcov = 0.3, thread = 4,
    diploidweight = 0.5, ...
) {
  # 1. Read and preprocess input files
  seg <- data.table::fread(seg)
  diploid_cov <- 100
  chr_levels <- c(as.character(1:22), "X", "Y")
  seg <- seg %>%
    dplyr::rowwise() %>%
    dplyr::mutate(Segcov = SegmentMeanToOriCov(gender = gender,
                                               chromosome = Chromosome,
                                               diploid_cov = diploid_cov,
                                               SM = Segment_Mean)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Chromosome = factor(Chromosome, levels = chr_levels)) %>%
    dplyr::arrange(Chromosome, Start) %>%
    dplyr::mutate(index = as.character(row_number()))

  # 2. Filter and prepare data
  autosome <- seg %>%
    dplyr::filter(!Chromosome %in% c("X", "Y")) %>%
    dplyr::mutate(Purity_estimate = ifelse(
      size >= modelminAIsize & FILTER != "FAILED" & MAF_Probes >= modelminprobes,
      "Include", "Exclude"
    ))
  observeddata <- data.frame(
    index = autosome$index,
    Segcov = autosome$Segcov,
    BAF = autosome$MAF,
    Tag = autosome$Purity_estimate,
    k = categorize_K(K = autosome$MAF_Probes)
  )

  # 3. Model source
  opt <- list(modelminprobes = modelminprobes, modelminAIsize = modelminAIsize, minsf = minsf, callcov = callcov, gender = gender)
  model_source <- ModelSource(seg, opt)

  # 4. Estimate variance
  var_sf <- EstimateSFVariance(x = observeddata[which(abs(observeddata$BAF - 0.5) <= 0.05), "Segcov"])

  # 5. Generate combinations
  sf_range <- round(seq(0.3, 3, length.out = 55), 3)
  purity_range <- round(seq(0.1, 1, length.out = 30), 3)
  purity_sf <- expand.grid(mu = sf_range, rho = purity_range)

  # 6. Calculate likelihoods in parallel
  num_cores <- thread
  cl <- doParallel::makeCluster(num_cores)
  doParallel::registerDoParallel(cl)
  chunks <- split(purity_sf, rep(1:num_cores, each = ceiling(nrow(purity_sf) / num_cores)))
  doParallel::clusterExport(cl, c(
    "RunCallikelihood", "AssignPriors", "Calccf",
    "Callikelihood", "CalSegmentLikelihood", "GenerateCombinations",
    "observeddata", "var_sf", "opt", "purity_sf", "chunks"
  ))
  doParallel::clusterEvalQ(cl, {
    library(dplyr); library(reshape2); library(tidyr)
  })
  results <- foreach::foreach(
    chunk_df = chunks,
    .combine = function(...) {
      res_list <- list(...)
      valid_res <- Filter(function(x) !is.null(x) && nrow(x) > 0, res_list)
      if (length(valid_res) > 0) do.call(rbind, valid_res) else NULL
    },
    .maxcombine = length(chunks),
    .multicombine = TRUE,
    .errorhandling = "pass",
    .packages = c("dplyr", "reshape2", "tidyr")
  ) %dopar% {
    tryCatch({
      result <- RunCallikelihood(
        purity_sf = chunk_df,
        data = observeddata,
        sigma_C = var_sf,
        opt = opt
      )
      return(result)
    }, error = function(e) {
      message("Chunk Error [", Sys.getpid(), "]: ", conditionMessage(e))
      NULL
    })
  }
  doParallel::stopCluster(cl)

  # 7. Rank and select calls
  results <- results %>%
    dplyr::mutate(ccf_BAF = CcfLOH(BAF = BAF, rho = rho, minor = minor, major = major))
  write.table(results, file = file.path(out_dir, paste0(prefix, "_likelihood_raw.tsv")), row.names = FALSE, quote = FALSE, sep = "\t")
  top_likelihood_rows <- SelectCallpersegment(results = results)
  write.table(top_likelihood_rows, file = file.path(out_dir, paste0(prefix, "_top_likelihood_calls.tsv")), row.names = FALSE, quote = FALSE, sep = "\t")

  # 8. Summarize and plot
  if (model_source == "Coverage + BAF") {
    return_models <- SelectFinalModel(
      top_likelihood_rows = top_likelihood_rows,
      groupinfo = observeddata,
      opt = opt,
      model_source = model_source
    )
    PlotModel(data = return_models$models,
              opt = opt,
              max_L_mu = return_models$Final_model$Final_mu,
              max_L_rho = return_models$Final_model$Final_rho)
  } else if (model_source == "Coverage") {
    return_models <- EstimatePurityCov(seg = seg)
    PlotCovDisCN(dis = return_models$dis_df,
                 opt = opt,
                 purity = return_models$Final_model$Final_rho,
                 min_dis = return_models$min_dis)
  } else {
    return_models <- list(Final_model = list(Final_mu = 1, Final_rho = 1))
  }

  raw_call <- ExtractCall(
    df = top_likelihood_rows,
    max_L_mu = return_models$Final_model$Final_mu,
    max_L_rho = return_models$Final_model$Final_rho,
    seg = seg
  )
  refined_call <- RefineCalls(
    df = raw_call,
    max_L_mu = return_models$Final_model$Final_mu,
    max_L_rho = return_models$Final_model$Final_rho,
    opt = opt
  )
  final_call <- RefineCallsSecond(
    df = refined_call,
    results = results,
    final_mu = return_models$Final_model$Final_mu,
    final_rho = return_models$Final_model$Final_rho,
    opt = opt
  )
  final_call$Model_source <- model_source
  final_call$rho <- return_models$Final_model$Final_rho
  final_call$mu <- return_models$Final_model$Final_mu
  write.table(final_call, file = file.path(out_dir, paste0(prefix, "_final_calls.tsv")), row.names = FALSE, quote = FALSE, sep = "\t")
  if (!is.null(return_models$dis_df)) {
    write.table(return_models$dis_df, file = file.path(out_dir, paste0(prefix, "_Models_dis.tsv")), row.names = FALSE, quote = FALSE, sep = "\t")
  }
  invisible(list(
    likelihood_raw = file.path(out_dir, paste0(prefix, "_likelihood_raw.tsv")),
    top_likelihood_calls = file.path(out_dir, paste0(prefix, "_top_likelihood_calls.tsv")),
    final_calls = file.path(out_dir, paste0(prefix, "_final_calls.tsv"))
  ))
}
