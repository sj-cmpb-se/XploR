#' Run Model Likelihood and Copy Number Calling Workflow
#'
#' This function runs the complete workflow for likelihood calculation and model selection.
#'
#' @param seg Path to the combined segment file.
#' @param out_dir Output directory.
#' @param prefix Output file prefix.
#' @param gender Sample gender ("male" or "female").
#' @param lambda, gamma, epsilon, modelminprobes, modelminAIsize, minsf, callcov, thread, callcovcutoff, callaicutoff, minsnpcallaicutoff Numeric parameters for the workflow (see script RunCallikelihood for details).
#' @return Invisibly returns the output file paths.
#' @importFrom data.table fread
#' @importFrom dplyr filter mutate arrange rowwise ungroup select group_by summarise desc left_join
#' @importFrom tidyr replace_na
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster clusterEvalQ clusterExport
#' @export
RunModelLikelihood <- function(
    seg, out_dir, prefix, gender,
    lambda = 1, gamma = 1, epsilon = 0.01, modelminprobes = 20,
    modelminAIsize = 5000000, minsf = 0.4, callcov = 0.3, thread = 4,
    ...
) {

  print("Run likelihood Parameters are:")
  print(paste0("Segment file is: ", seg))
  print(paste0("Output dir is: ", out_dir))
  print(paste0("Prefix is: ", prefix))
  print(paste0("Gender is: ", gender))
  print(paste0("Lambda is: ", lambda))
  print(paste0("gamma is: ", gamma))
  print(paste0("epsilon is: ", epsilon))
  print(paste0("modelminprobes is: ", modelminprobes))
  print(paste0("modelminAIsize is: ", modelminAIsize))
  print(paste0("minsf is: ", minsf ))
  print(paste0("callcov is: ", callcov ))
  print(paste0("thread is: ", thread ))


  # 1. Read and preprocess input files
  print(paste0("Reading segment file: ", seg))
  seg <- data.table::fread(seg)
  diploid_cov <- 100
  chr_levels <- c(as.character(1:22), "X", "Y")
  print(paste0("Processing segment file: nrow=", nrow(seg)))
  seg <- seg %>%
    dplyr::rowwise() %>%
    dplyr::mutate(Segcov = SegmentMeanToOriCov(gender = gender,
                                               chromosome = Chromosome,
                                               diploid_cov = diploid_cov,
                                               SM = Segment_Mean)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Chromosome = factor(Chromosome, levels = chr_levels)) %>%
    dplyr::arrange(Chromosome, Start) %>%
    dplyr::mutate(index = as.character(dplyr::row_number()))

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
    MAF = autosome$MAF,
    Tag = autosome$Purity_estimate,
    k = CategorizeK(K = autosome$MAF_Probes)
  )

  # 3. Model source

  model_source <- ModelSource(seg_df = seg, modelminprobes = modelminprobes )
  print(paste0("Estimated model source is: ", model_source))
  # 4. Estimate variance of diploid region
  var_sf <- EstimateVariance(x = observeddata[which(abs(observeddata$MAF - 0.5) <= 0.05 ), "Segcov"])

  # 5. Generate combinations
  sf_range <- round(seq(0.3, 3, length.out = 55), 3)
  purity_range <- round(seq(0.1, 1, length.out = 30), 3)
  purity_sf <- expand.grid(mu = sf_range, rho = purity_range)
  print(paste0("Calculating likelihood for possible models: n=", nrow(purity_sf)))
  # 6. Calculate likelihoods in parallel
  num_cores <- thread
  cl <- parallel::makeCluster(num_cores)
  doParallel::registerDoParallel(cl)
  chunks <- split(purity_sf, rep(1:num_cores, each = ceiling(nrow(purity_sf) / num_cores)))
 # parallel::clusterExport(cl, c(
#    "observeddata", "var_sf", "purity_sf", "chunks","lambda",
#    "gamma", "epsilon", "modelminprobes", "modelminAIsize",
#    "minsf", "callcov", "thread"
#  ))
  parallel::clusterEvalQ(cl, { library(XploR)})
  results <- foreach(
    chunk_df = chunks,
    .combine = function(...) {
      res_list <- list(...)
      valid_res <- Filter(function(x) !is.null(x) && nrow(x) > 0, res_list)
      if (length(valid_res) > 0) do.call(rbind, valid_res) else NULL
    },
    .maxcombine = length(chunks),
    .multicombine = TRUE,
    .errorhandling = "pass",
    .packages = c("XploR"),
    .export = c(
      "observeddata", "var_sf", "purity_sf", "chunks","lambda",
      "gamma", "epsilon", "modelminprobes", "modelminAIsize",
      "minsf", "callcov", "thread"
    )
  ) %dopar% {
    tryCatch({
      result <- RunCallikelihood(
        purity_sf = chunk_df,
        data = observeddata,
        sigma_C = var_sf,
        lambda = lambda,
        gamma = gamma,
        epsilon = epsilon
      )
      message("Worker finished on PID: ", Sys.getpid())
      return(result)
    }, error = function(e) {
      list(
        error = TRUE,
        pid = Sys.getpid(),
        error_message = conditionMessage(e),
        error_class = class(e),
        error_call = deparse(conditionCall(e))
      )
    })
  }
  parallel::stopCluster(cl)
  if (is.list(results) && any(sapply(results, function(x) is.list(x) && !is.null(x$error) && x$error))) {
    print("Some workers failed:")
    print(Filter(function(x) is.list(x) && !is.null(x$error) && x$error, results))
  }
  if( ! is.null(results) ){
    print(paste0("Likelihood calculation DONE."))
  }else{
    stop("Likelihood calculation failed: 'results' is NULL. Check previous error messages for details.")

  }
  # 7. Rank and select calls
  results <- results %>%
    dplyr::mutate(ccf_MAF = CcfLOH(MAF = MAF, rho = rho, minor = minor, major = major))
  write.table(results, file = file.path(out_dir, paste0(prefix, "_likelihood_raw.tsv")), row.names = FALSE, quote = FALSE, sep = "\t")
  print(paste0("Raw likelihood calculation results is saved at: ", paste0(out_dir,"/", paste0(prefix, "_likelihood_raw.tsv")) ))
  top_likelihood_rows <- SelectCallpersegment(results = results)
  write.table(top_likelihood_rows, file = file.path(out_dir, paste0(prefix, "_top_likelihood_calls.tsv")), row.names = FALSE, quote = FALSE, sep = "\t")
  print(paste0("Top likelihood allelic combinations of each segment are saved at: ", paste0(out_dir,"/", paste0(prefix, "_top_likelihood_calls.tsv")) ))

   # 8. Summarize and plot
  print("Selecting final models.")
  if (model_source == "Coverage + MAF") {
    return_models <- SelectFinalModel(
      results = results,
      callcov = callcov,
      modelminprobes= modelminprobes,
      modelminAIsize = modelminAIsize,
      minsf = minsf,
      top_likelihood_rows = top_likelihood_rows,
      groupinfo = observeddata,
      model_source = model_source,
      prefix = prefix,
      gender = gender,
      out_dir = out_dir,
      seg = seg
    )

    PlotModel(data = return_models$models,
              prefix = prefix,
              out_dir = out_dir,
              max_L_mu = return_models$Final_model$Final_mu,
              max_L_rho = return_models$Final_model$Final_rho)

    PlotModelcluster(
      models = models,
      refined_calls = refined_calls,
      out_dir = out_dir,
      prefix = prefix

    )


  } else if (model_source == "Coverage") {
    return_models <- EstimatePurityCov(seg = seg,gender = gender)
    PlotCovDisCN(dis = return_models$dis_df,
                 prefix = prefix,
                 out_dir = out_dir,
                 purity = return_models$Final_model$Final_rho,
                 min_dis = return_models$min_dis)
  } else {
    return_models <- list(Final_model = list(Final_mu = 1, Final_rho = 1))
  }
  print("Extracting calls under final model...........")
  raw_call <- ExtractCall(
    df = top_likelihood_rows,
    max_L_mu = return_models$Final_model$Final_mu,
    max_L_rho = return_models$Final_model$Final_rho,
    seg = seg
  )
  print("Refining calls...............")
  refined_call <- RefineCalls(
    df = raw_call,
    max_L_mu = return_models$Final_model$Final_mu,
    max_L_rho = return_models$Final_model$Final_rho,
    gender = gender
  )
  final_call <- RefineCallsSecond(
    df = refined_call,
    results = results,
    final_mu = return_models$Final_model$Final_mu,
    final_rho = return_models$Final_model$Final_rho,
    callcov = callcov,
    gender = gender
  )
  print(paste0("Reporting final calls at: ", paste0(out_dir,"/",paste0(prefix, "_final_calls.tsv"))))
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
