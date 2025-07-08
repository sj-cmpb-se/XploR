#' Cytoband reference table
#'
#' This dataset contains cytoband annotations for the human genome.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{chrom}{Chromosome name (character).}
#'   \item{start}{Start coordinate (integer).}
#'   \item{end}{End coordinate (integer).}
#'   \item{band}{Cytoband name (character).}
#'   \item{gieStain}{Giemsa stain (character).}
#' }
#' @source UCSC Table Browser
"cytoband_table"

#' Detectable Edge table for cytoband annotation follows the 5MB rule.
#'
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{Chrom}{Chromosome name (character).}
#'   \item{p_chromStart}{P arm chromosome Start coordinates (integer).}
#'   \item{p_chromEnd}{P arm chromosome End coordinates (integer).}
#'   \item{p_first_name}{P arm first cytoband name (character).}
#'   \item{p_last_name}{P arm last cytoband name (character).}
#'   \item{q_chromStart}{Q arm chromosome Start coordinates (integer).}
#'   \item{q_chromEnd}{Q arm chromosome End coordinates (integer).}
#'   \item{q_first_name}{Q arm first cytoband name (character).}
#'   \item{q_last_name}{Q arm last cytoband name (character).}
#' }
#' @source Internal analysis
"edge_table"

#' Non mapable regions of the genome
#'
#' A table listing genomic intervals known to be problematic for mapping.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{chrom}{Chromosome name (character).}
#'   \item{start}{Start coordinate (integer).}
#'   \item{end}{End coordinate (integer).}
#'   \item{tag}{ Could be ignored (character).}
#' }
#' @source DRAGEN4.2 CNV output *.cnv.excluded_intervals.bed.gz file
"blacklist_table"

#' RefSeq gene annotations
#'
#' RefSeq annotations including gene names and start and end.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{chrom}{Chromosome name (character).}
#'   \item{start}{Start coordinate (integer).}
#'   \item{end}{End coordinate (integer).}
#'   \item{genesymbl}{Gene name (character).}
#' }
#' @source RefSeq
"gene_annotation"


#' Non blacklist region for each cytobands binned by 300nt for female.
#'
#' File used in plot CNV profile to keep only non-blacklist region.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{cytoband}{cytoband name (character).}
#'   \item{coordinates}{ coordinates binned 300nt (integer).}
#' }
#' @source DRAGEN
"band_plottable_female"


#' Non blacklist region for each cytobands binned by 300nt for male.
#'
#' File used in plot CNV profile to keep only non-blacklist region.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{cytoband}{cytoband name (character).}
#'   \item{coordinates}{ coordinates binned 300nt (integer).}
#' }
#' @source DRAGEN
"band_plottable_male"
