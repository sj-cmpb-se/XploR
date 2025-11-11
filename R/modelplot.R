#' Plot Coverage Distance to Integer Copy Number vs. Purity
#'
#' Plots the distance to integer copy number (\code{dis_integer_CN}) as a function of purity (\code{rho}), with the minimum distance point highlighted.
#'
#' @param dis Data frame or tibble. Must include columns \code{rho} and \code{dis_integer_CN}.
#' @param purity Numeric. The selected purity value to highlight.
#' @param min_dis Numeric. The minimum distance to integer copy number to highlight.
#' @param prefix Character. Output file prefix.
#' @param out_dir Character. Output directory for saving the plot.
#'
#' @return Invisibly returns the ggplot object (or \code{FALSE} if plotting failed). The plot is saved as a PNG file in the output directory.
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point theme_minimal labs ggsave
#' @examples
#' # PlotCovDisCN(dis, purity = 0.7, min_dis = 2, prefix = "Sample1", out_dir = "results/")
#'
#' @export
PlotCovDisCN <- function(dis, purity, min_dis, prefix, out_dir ){
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
    outFile <- paste0( prefix, "_CNV_Coverage_purity_estimate.png")

    ggplot2::ggsave(filename = outFile,
           plot = displot,
           device = "png",
           path = out_dir,
           width = 8,height = 8,dpi = 600, bg = 'white' )
    print(paste0("Plotting is DONE. The plot is saved at: ", out_dir,"/", outFile))

  }}



#' Plot Model Total Log-Likelihoods as a Dot Plot
#'
#' Plots the total log-likelihood of each model as a dot plot, with the tier 1 index highlighted as a vertical dashed line.
#'
#' @param models Data frame or tibble. Must include a \code{total_log_likelihood} column.
#' @param prefix Character. Output file prefix.
#' @param out_dir Character. Output directory for saving the plot.
#' @param tier1_index Integer. The index (rank) of the tier 1 model to highlight.
#'
#' @return Invisibly returns the ggplot object (or \code{FALSE} if plotting failed). The plot is saved as a PNG file in the output directory.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_vline ylab xlab theme_bw ggsave
#' @examples
#' # PlotModelDot(models, prefix = "Sample1", out_dir = "results/", tier1_index = 10)
#'
#' @export
PlotModelDot<- function(models, prefix, out_dir, tier1_index) {
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
    outFile <- paste0( prefix, "_CNV_likelihood_dot_plot.png")

    ggplot2::ggsave(filename = outFile,
           plot = modeldot,
           device = "png",
           path = out_dir,
           width = 8,height = 8,dpi = 600, bg = 'white' )
    print(paste0("Plotting is DONE. The plot is saved at: ", out_dir,"/", outFile))
  }
}



#' Plot Model Total Log-Likelihoods as a Heatmap
#'
#' Plots a heatmap of total log-likelihoods for (mu, rho) combinations and highlights the selected model.
#'
#' @param data Data frame or tibble. Must include columns \code{rho}, \code{mu}, and \code{total_log_likelihood_before_refine}.
#' @param prefix Character. Output file prefix.
#' @param out_dir Character. Output directory for saving the plot.
#' @param max_L_mu Numeric. The selected value of mutation multiplicity (\code{mu}) to highlight.
#' @param max_L_rho Numeric. The selected value of tumor purity (\code{rho}) to highlight.
#'
#' @return Invisibly returns the ggplot object (or \code{FALSE} if plotting failed).
#'
#' @importFrom ggplot2 ggplot aes geom_tile geom_point scale_fill_gradientn xlab ylab theme_bw ggsave
#' @importFrom scales rescale
#' @importFrom grDevices colorRampPalette
#'
#' @examples
#' # PlotModel(data, opt, max_L_mu = 1, max_L_rho = 0.7)
#'
#' @export
PlotModel <- function( data, prefix, out_dir, max_L_mu, max_L_rho ){
  data$total_log_likelihood <- data$total_log_likelihood_before_refine

  p <- tryCatch({


    tier1_min <- data %>% dplyr::filter( !is.na(Tier1))
    tier1_min <- min(tier1_min$total_log_likelihood, na.rm = T)
    median_ll <- median(data$total_log_likelihood, na.rm = T)
    if(tier1_min > median_ll){
      range <- c(min(data$total_log_likelihood, na.rm = T), median_ll, tier1_min, max(data$total_log_likelihood, na.rm = T))
      my_ramp <- grDevices::colorRampPalette(c("navy", "yellow", "orange"))
      color_palette <- c(
        "navy",
        my_ramp(4)[2],
        my_ramp(4)[3],
        "orange"
      )
      }else if( tier1_min < median_ll){
      range <- c(min(data$total_log_likelihood, na.rm = T), tier1_min, median_ll, max(data$total_log_likelihood, na.rm = T))
      my_ramp <- grDevices::colorRampPalette(c("navy", "yellow", "orange"))
      color_palette <- c(
        "navy",
        my_ramp(4)[2],
        my_ramp(4)[3],
        "orange"
      )
      } else{
         range <- c(min(data$total_log_likelihood, na.rm = T), median_ll, max(data$total_log_likelihood, na.rm = T))
         my_ramp <- grDevices::colorRampPalette(c("navy", "yellow", "orange"))
         color_palette <- c(
           "navy",
           my_ramp(3)[2],
           "orange"
         )
       }


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
    outFile <- paste0( prefix, "_CNV_Modelplot.pdf")

    ggsave(filename = outFile,
           plot = model_plot,
           device = "pdf",
           path = out_dir,
           width = 8,height = 8,dpi = 600, bg = 'white' )
    print(paste0("Plotting is DONE. The plot is saved at: ", out_dir,"/", outFile))
  }

}


#' Plot all CNV calls under each Tier1 models
#'
#' Generates an interactive plot (HTML) of clustered CNV models, showing copy number profiles for Tier 1 models. The y-axis represents model IDs (clustered by CNV profile similarity), the x-axis shows genomic coordinates, and color indicates copy number (CN). Tooltip text provides model summary statistics.
#'
#' @param models Data frame of model likelihood results. Must include columns \code{mu}, \code{rho}, \code{Tier1}, and summary statistics (\code{total_likelihood_cluster}, \code{diploid_distance_cluster}, \code{nondiploid_distance_cluster}).
#' @param refined_calls Data frame of segment calls. Must include columns \code{Chromosome}, \code{Start}, \code{End}, \code{CN}, \code{mu}, and \code{rho}.
#' @param out_dir Character. Output directory for the HTML plot.
#' @param prefix Character. Prefix for output file name.
#'
#' @return Invisibly returns \code{TRUE} if the plot was successfully created and saved, otherwise returns \code{FALSE}.
#'
#' @details
#' The function filters for Tier 1 models, clusters them by CNV profile similarity, and generates an interactive heatmap using \code{ggplot2} and \code{plotly}. The plot is saved as an HTML file in \code{out_dir}.
#'
#' @importFrom dplyr filter mutate select left_join
#' @importFrom tidyr pivot_wider
#' @importFrom ggplot2 ggplot geom_segment scale_color_gradientn facet_grid theme_minimal labs theme element_blank element_text element_line
#' @importFrom plotly ggplotly
#' @importFrom htmlwidgets saveWidget
#' @importFrom stats hclust dist
#'
#' @examples
#' \dontrun{
#' PlotModelcluster(
#'   models = models_df,
#'   refined_calls = calls_df,
#'   out_dir = "results/",
#'   prefix = "sample1"
#' )
#' }
#'
#' @export
PlotModelcluster<- function( models, refined_calls, out_dir, prefix ){
  # only include tier1 models
  models <- models %>%
    dplyr::filter( !is.na(Tier1)) %>%
    dplyr::mutate( id = paste( mu,rho, sep=":"))

  calls <- refined_calls %>%
    dplyr::select(Chromosome, Start, End, CN, mu, rho) %>%
    dplyr::mutate( id = paste(mu, rho, sep = ":")) %>%
    dplyr::filter( id %in% models$id) %>%
    dplyr::select(-mu,-rho) %>%
    dplyr::mutate( pos_id = paste( Chromosome, Start, End , sep=":")) %>%
    dplyr::ungroup()


  # cluster calls
  calls_matrix <- calls %>%
    tidyr::pivot_wider(
      id_cols = id,
      names_from = pos_id,
      values_from = CN,
      values_fill = NA )

  row_ids <- calls_matrix$id
  calls_matrix <- as.matrix(calls_matrix[,-1])

  row_order <- hclust(dist(calls_matrix))$order
  ordered_ids <- row_ids[row_order]
  calls$id <- factor(calls$id, levels = ordered_ids)

  max_cn <- max(calls$CN,na.rm = T)
  max_cn <- ifelse(max_cn >6,6,max_cn)
  calls <- calls %>% dplyr::mutate(CN = ifelse(CN > max_cn, max_cn,CN)) %>% dplyr::left_join(models,by = "id")

  p <- tryCatch({
    ggplot2::ggplot() +
      ggplot2::geom_segment(
        data = calls,
        ggplot2::aes(
          x = Start, xend = End, y = id, yend = id, color = CN,
          text = paste(
            "SF: ", mu, "<br>",
            "Purity: ", rho, "<br>",
            "CN: ", CN, "<br>",
            "Total log likelihood cluster:", total_likelihood_cluster, "<br>",
            "diploid_distance_cluster:", diploid_distance_cluster, "<br>",
            "nondiploid_distance_cluster:", nondiploid_distance_cluster, "<br>"
          )
        ),
        linewidth = 1
      ) +
      ggplot2::scale_color_gradientn(
        colors = c("blue", "white", "red"),
        values = scales::rescale(c(0, 2, max_cn)),
        limits = c(0, max_cn),
        breaks = 0:max_cn,
        name = "Copy Number"
      ) +
      ggplot2::facet_grid(cols = vars(Chromosome), scales = "free_x", space = "free_x") +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = prefix, y = "SF:Purity") +
      ggplot2::theme(
        panel.spacing = grid::unit(0, "lines"),
        strip.background = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(face = "bold", size = 10),
        plot.title = ggplot2::element_text(hjust = 0.5),
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.line.y.left = ggplot2::element_line(color = "black"),
        axis.ticks.y = ggplot2::element_line(color = "black"),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        legend.position = "bottom"
      )
  }, error = function(e) {
    message("An error occurred while trying to plot: ", conditionMessage(e))
    return(NULL)
  })
  if (!is.null(p)) {
    outFile <- paste0(out_dir, "/", prefix, "_Tier1_models.html")
    model_plot <- plotly::ggplotly(p, tooltip = "text")
    htmlwidgets::saveWidget(model_plot, file = outFile, selfcontained = TRUE)
  }

}
