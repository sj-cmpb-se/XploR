#' Plot Coverage Distance to Integer Copy Number vs. Purity
#'
#' Plots the distance to integer copy number (\code{dis_integer_CN}) against purity (\code{rho}) and highlights the minimum distance.
#'
#' @param dis Data frame or tibble. Must include columns \code{rho} and \code{dis_integer_CN}.
#' @param opt List of options. Must include \code{out_dir} and \code{prefix}.
#' @param purity Numeric. The selected purity value to highlight.
#' @param min_dis Numeric. The minimum distance to integer copy number to highlight.
#'
#' @return Invisibly returns the ggplot object (or \code{FALSE} if plotting failed).
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point theme_minimal labs ggsave
#' @examples
#' # PlotCovDisCN(dis, opt, purity = 0.7, min_dis = 2)
#'
#' @export
PlotCovDisCN <- function(dis, opt, purity, min_dis ){
  p <- tryCatch({
    displot <- ggplot2::ggplot(dis, aes(x = rho, y = dis_integer_CN)) +
      ggplot2::geom_line(color = "grey", linewidth = 1) +
      ggplot2::geom_point(color = "lightblue", size = 3) +
      ggplot2::geom_point( aes(x = purity, y = min_dis ), color = "blue", size = 3 ) +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "",
           x = "Purity",
           y = "Distance to CN")
    TRUE
  }, error = function(e) {
    message <- conditionMessage(e)
    print(paste("An error occurred while trying to plot: ", message))
    FALSE
  })

  if(p) {
    outFile <- paste0( opt$prefix, "_CNV_Coverage_purity_estimate.png")

    ggplot2::ggsave(filename = outFile,
           plot = displot,
           device = "png",
           path = opt$out_dir,
           width = 8,height = 8,dpi = 600, bg = 'white' )
    print(paste0("Plotting is DONE. The plot is saved at: ", opt$out_dir,"/", outFile))

  }}



#' Plot Model Total Log-Likelihoods as a Dot Plot
#'
#' Plots the total log-likelihood of each model as a dot plot and marks the tier 1 index as a vertical line.
#'
#' @param models Data frame or tibble. Must include a \code{total_log_likelihood} column.
#' @param opt List of options. Must include \code{out_dir} and \code{prefix}.
#' @param tier1_index Integer. The index (rank) of the tier 1 model to highlight.
#'
#' @return Invisibly returns the ggplot object (or \code{FALSE} if plotting failed).
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_vline ylab xlab theme_bw ggsave
#' @examples
#' # PlotModelDot(models, opt, tier1_index = 10)
#'
#' @export
PlotModelDot<- function(models, opt, tier1_index) {
  a <- tryCatch({
    modeldot <- ggplot2::ggplot(models, aes(y = total_log_likelihood, x = 1:nrow(models))) +
      ggplot2::geom_point(size = 2, color = "blue", alpha = 0.5 ) +
      ggplot2::geom_vline(xintercept = tier1_index, linetype = "dashed", color = "grey") +
      ggplot2::ylab(" Total likelihood ") +
      ggplot2::xlab(" Rank ") +
      ggplot2::theme_bw()
    TRUE
  }, error = function(e) {
    message <- conditionMessage(e)
    print(paste("An error occurred while trying to plot: ", message))
    FALSE
  })

  if(a) {
    outFile <- paste0( opt$prefix, "_CNV_likelihood_dot_plot.png")

    ggplot2::ggsave(filename = outFile,
           plot = modeldot,
           device = "png",
           path = opt$out_dir,
           width = 8,height = 8,dpi = 600, bg = 'white' )
    print(paste0("Plotting is DONE. The plot is saved at: ", opt$out_dir,"/", outFile))
  }
}



#' Plot Model Likelihood and Distance Clusters as an Interactive Heatmap
#'
#' Plots an interactive heatmap of total likelihood cluster assignments, with the best model highlighted, and saves as an HTML file.
#'
#' @param df Data frame or tibble. Must include columns \code{mu}, \code{rho}, \code{total_likelihood_cluster}, \code{diploid_distance_cluster}, \code{nondiploid_distance_cluster}.
#' @param Final_model List. Output from \code{\link{ClusterModels}}, must include \code{Final_mu} and \code{Final_rho}.
#' @param opt List of options. Must include \code{out_dir} and \code{prefix}.
#'
#' @return Invisibly returns the interactive plotly object (or \code{FALSE} if plotting failed).
#'
#' @importFrom ggplot2 ggplot aes geom_tile geom_point scale_fill_gradient ylab xlab theme_bw element_blank
#' @importFrom plotly ggplotly
#' @importFrom htmlwidgets saveWidget
#' @examples
#' # PlotLikelihoodAndDistance(df, Final_model, opt)
#'
#' @export

PlotLikelihoodAndDistance <- function(df , Final_model, opt){
  a <- tryCatch({
    model_plot <- ggplot2::ggplot(df ) +
      ggplot2::geom_tile(aes(y = mu,
                    x = rho,
                    fill = total_likelihood_cluster,
                    text = paste("SF: ", mu, "<br>",
                                 "Purity: ", rho, "<br>",
                                 "Total log likelihood cluster:", total_likelihood_cluster, '<br>',
                                 "diploid_distance_cluster:",diploid_distance_cluster,'<br>',
                                 "nondiploid_distance_cluster:", nondiploid_distance_cluster,'<br>' )),
                size = 2, alpha = 0.5 ) +
      ggplot2::scale_fill_gradient(low = 'yellow',high = 'navy') +
      ggplot2::geom_point(aes( y = Final_model$Final_mu,
                      x= Final_model$Final_rho,
                      text = paste("SF: ", mu, "<br>",
                                   "Purity: ", rho, "<br>",
                                   "Total log likelihood cluster:", total_likelihood_cluster, '<br>',
                                   "diploid_distance_cluster:",diploid_distance_cluster,'<br>',
                                   "nondiploid_distance_cluster:", nondiploid_distance_cluster,'<br>' )),
                 color = "red",size = 3) +
      ggplot2::ylab(" Size factor ") +
      ggplot2::xlab(" Tumor purity ") +
      ggplot2::theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
      ) + theme_bw()

    TRUE
  }, error = function(e) {
    message <- conditionMessage(e)
    print(paste("An error occurred while trying to plot: ", message))
    FALSE
  })

  if(a) {
    outFile <- paste0( opt$out_dir,"/", opt$prefix, "_CNV_likelihood_distance.html")
    model_plot <- plotly::ggplotly(model_plot, tooltip = "text")
    htmlwidgets::saveWidget(model_plot, file = outFile )

    print(paste0("Plotting is DONE. The plot is saved at: ", outFile))
  }
}


#' Plot Model Total Log-Likelihoods as a Heatmap
#'
#' Plots a heatmap of total log-likelihoods for (mu, rho) combinations and highlights the selected model.
#'
#' @param data Data frame or tibble. Must include columns \code{rho}, \code{mu}, and \code{total_log_likelihood_before_refine}.
#' @param opt List of options. Must include \code{out_dir} and \code{prefix}.
#' @param max_L_mu Numeric. The selected value of mutation multiplicity (\code{mu}) to highlight.
#' @param max_L_rho Numeric. The selected value of tumor purity (\code{rho}) to highlight.
#'
#' @return Invisibly returns the ggplot object (or \code{FALSE} if plotting failed).
#'
#' @importFrom ggplot2 ggplot aes geom_tile geom_point scale_fill_gradientn xlab ylab theme_bw ggsave
#' @importFrom scales rescale
#' @examples
#' # PlotModel(data, opt, max_L_mu = 1, max_L_rho = 0.7)
#'
#' @export
PlotModel <- function( data, opt, max_L_mu, max_L_rho ){
  data$total_log_likelihood <- data$total_log_likelihood_before_refine
  p <- tryCatch({
    range <- quantile(data$total_log_likelihood, probs = c(0, 0.33, 0.67, 1), na.rm = T)
    color_palette <- c("navy", "navy", "navy",  "yellow")

    # Create the plot
    model_plot <- ggplot2::ggplot(data, aes(x = rho, y = mu, fill = total_log_likelihood)) +
      ggplot2::geom_tile() +
      ggplot2::geom_point(aes(x = max_L_rho, y = max_L_mu), color = "lightblue", size = 5) +
      ggplot2::scale_fill_gradientn(colors = color_palette, values = scales::rescale(range) ) +
      ggplot2::xlab("Tumor Purity") +
      ggplot2::ylab("Adjusted Factor for Diploid Coverage") +
      ggplot2::theme_bw()
    TRUE
  }, error = function(e) {
    message <- conditionMessage(e)
    print(paste("An error occurred while trying to plot: ", message))
    FALSE
  })

  if(p) {
    outFile <- paste0( opt$prefix, "_CNV_Modelplot.pdf")

    ggsave(filename = outFile,
           plot = model_plot,
           device = "pdf",
           path = opt$out_dir,
           width = 8,height = 8,dpi = 600, bg = 'white' )
    print(paste0("Plotting is DONE. The plot is saved at: ", opt$out_dir,"/", outFile))
  }

}
