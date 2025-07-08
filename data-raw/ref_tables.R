## code to prepare `ref_tables` dataset goes here

cytoband_table <- read.table("data-raw/hg19_cytoBand.dat", stringsAsFactors = FALSE,header = T)
edge_table <- read.table("data-raw/hg19_detectable_edges.txt", stringsAsFactors = FALSE, header = T)
blacklist_table <- read.table("data-raw/extended_centromere_blacklist.bed", stringsAsFactors = FALSE, header = F )
colnames(blacklist_table) <- c("chrom","start","end","tag")
gene_annotation <- read.table("data-raw/RefSeqCurated.genePred.gene_region.txt", sep = "\t",stringsAsFactors = FALSE)
colnames(gene_annotation) <- c("chrom","start","end","genesymbl")
band_plottable_female <- read.csv("data-raw/band_positions_female.csv.gz", stringsAsFactors = FALSE,header = F)
colnames(band_plottable_female) <- c("cytoband","pos")
band_plottable_male <- read.csv("data-raw/band_positions_male.csv.gz", stringsAsFactors = FALSE,header = F)
colnames(band_plottable_male) <- c("cytoband","pos")

usethis::use_data(
  cytoband_table,
  edge_table,
  blacklist_table,
  gene_annotation,
  band_plottable_female,
  band_plottable_male,
  overwrite = TRUE
)
