#' Calculate Median Coverage for a Genomic Region
#'
#' @param chrom Chromosome name (character)
#' @param start Start position (numeric)
#' @param end End position (numeric)
#' @param cov Data frame recording coverage with columns CONTIG, START, END, COUNT
#'
#' @return Median coverage (numeric)
#' @importFrom dplyr filter
#' @export
CalMedian <- function( chrom, start, end, cov){

  tmp <- cov %>%
    dplyr::filter( CONTIG == chrom & START >= start & END <= end )
  return(median(tmp$COUNT,na.rm = T))
}


#' Calculate Baseline Coverage
#'
#' @param chrom Chromosome
#' @param cr Segment mean
#' @param count Coverage count
#' @return Baseline coverage (numeric)
#' @export
CalbaselineCov <- function( chrom, cr, count ){
  baseline <- count/(2^cr)
  return(baseline)

}

#' Fix Segment Mean Based on Gender Discrepancy
#'
#' @param sm Segment mean
#' @param gatkgender Gender from GATK
#' @param pipeline_gender Gender from clinical info
#' @return Adjusted segment mean
#' @export
FixsegmentMean<- function( sm, gatkgender, pipeline_gender ){
  if(gatkgender == "male" && pipeline_gender == "female"){
    sm <- sm-1
  }else if(gatkgender == "female" && pipeline_gender == "male" ) { sm <- sm + 1 }
  return(sm)
}


#' Run automatic gender detection
#'
#' Detects sample gender based on coverage and segmentation data, and returns the predicted gender.
#' If gender information is missing or set to "unknown", the function will infer gender automatically.
#'
#' @param cov Data frame. Coverage data.
#' @param seg Character. Path to segmentation data file.
#' @param gender Character. Provided gender information; use "unknown" to trigger automatic detection.
#'
#' @return Character. Predicted gender ("male" or "female").
#'
#' @importFrom dplyr filter mutate
#' @importFrom data.table fread
#' @export
Runcheckgender <- function( cov, seg, gender, out_dir, prefix){

  seg_df <- data.table::fread(seg) %>%
    dplyr::mutate(size = End - Start)

  # Step 1: Gender check
  gender = tolower(gender)
  # Read coverage file
  lines <- readLines(cov)
  filtered_lines <- grep("^\\s*@", lines, invert = TRUE, value = TRUE)
  cov_df <- data.table::fread(text = filtered_lines)
  seg_df <- CheckGender(cov = cov_df, seg = seg_df, gender = gender)
  gender <- as.character(seg_df$pipeline_gender[1])
  outFile <- file.path(out_dir, paste0(prefix, "_gender.txt"))
  write.table(gender,file=outFile,quote = F,row.names = F)
  return(gender)
}



#' Check and Adjust Gender in GATK Segmentation Data
#'
#' @param cov Coverage data frame
#' @param seg Segmentation data frame
#' @param gender Gender from clinical information (character, "male", "female" or "unknown")
#' @return Modified segmentation seg data frame with gender information
#' @importFrom dplyr filter rowwise mutate
#' @importFrom stats median
#' @export
CheckGender <- function( cov, seg, gender ){
  ## check what is the baseline cov is used in autosome and x and Y separately.

  seg <- seg %>%
    dplyr::rowwise() %>%
    dplyr::mutate( Count = CalMedian( chrom = Chromosome, start = Start, end = End, cov = cov )) %>%
    #filter( Num_Probes >= 500 | Chromosome %in% c("X","Y")) %>%
    dplyr::rowwise() %>%
    dplyr::mutate( Baseline_cov = CalbaselineCov( chrom = Chromosome, cr = Segment_Mean, count = Count ))

  autobasecov <- seg %>%
    dplyr::filter( ! Chromosome %in% c("X","Y") ) %>%
    dplyr::filter( Num_Probes >= 500 )
 # autobasecov <- median(autobasecov$Baseline_cov)
  autobasecov <- sum(autobasecov$Baseline_cov * autobasecov$Num_Probes) /sum( autobasecov$Num_Probes)
  xbasecov <- seg %>% dplyr::filter( Chromosome == "X" )
  xbasecov <- sum(xbasecov$Baseline_cov * xbasecov$Num_Probes) /sum( xbasecov$Num_Probes)
   ybasecov <- seg %>% dplyr::filter( Chromosome == "Y" )

  if( nrow(ybasecov) > 0){
    ybasecov <- (ybasecov$Baseline_cov * ybasecov$Num_Probes) /sum( ybasecov$Num_Probes)
    print( paste0("The Y baseline cov/autosome is: ", ybasecov/autobasecov))

    print( paste0("Y chromosome baseline cov is: ", ybasecov))
  }

  print( paste0("The autosome baseline cov is: ", autobasecov))
  print( paste0("X chromosome baseline cov is: ", xbasecov))

  print( paste0("The X baseline cov/autosome is: ", xbasecov/autobasecov))
  if( xbasecov <= autobasecov/1.5 ){
    gatkgender <- "male"
  }else{gatkgender <- "female"}

  seg$gatk_gender <- gatkgender
  seg$Segment_Mean_raw <- seg$Segment_Mean

  if( gender == "unknown"){
    message("Automatic gender assginded as : ", gatkgender)
    seg$pipeline_gender <- gatkgender
  }else if (gender != "unknown" && gatkgender == gender){
    message("GATK gender is correct!")
    seg$pipeline_gender <- gender
  }else if (gender != "unknown" && gatkgender != gender ){
    message("GATK gender is wrong, fixing the X chromosome from the seg file.")
    seg$pipeline_gender <- gender
    seg <- seg %>%
      dplyr::rowwise() %>%
      dplyr::mutate(Segment_Mean = ifelse( gatk_gender != gender & Chromosome == "X",
                                    FixsegmentMean( sm = Segment_Mean,
                                                    gatkgender = gatk_gender,
                                                    pipeline_gender = gender ), Segment_Mean ))

  }
  return(seg)
}




#' Check Whether Two Segments Should Be Merged (Copy Number)
#'
#' Determines whether two adjacent copy number segments should be merged based on the difference in their segment means
#' and whether their calls are identical.
#'
#' @param cur_row A data frame row (list or tibble row) representing the current segment. Must have \code{Segment_Mean} and \code{Call}.
#' @param next_row A data frame row (list or tibble row) representing the next segment. Must have \code{Segment_Mean} and \code{Call}.
#' @param mergecov Numeric threshold for segment mean difference
#'
#' @return Logical value: \code{TRUE} if the two segments should be merged, \code{FALSE} otherwise.
#'
#' @export
MergeSegCheck <- function(cur_row,next_row, mergecov){
  # check merge conditions
  # ai segment diff <= mergeai

  if(abs( as.numeric(cur_row$Segment_Mean) - as.numeric(next_row$Segment_Mean)) <= mergecov &
     cur_row$Call == next_row$Call ){
    re <- TRUE
  }else{ re <- FALSE}
  return(re)
}


#' Merge Adjacent Copy Number Segments Based on Similarity Criteria
#'
#' Iteratively merges adjacent rows in a data frame of copy number segments if they meet criteria defined by \code{MergeSegCheck}.
#' For merged segments, the segment with the greater number of probes provides the call, segment mean, count, and baseline coverage.
#'
#' @param df A data frame or tibble of copy number segments, with columns such as Sample, Chromosome, Start, End, Num_Probes, Call, Segment_Mean, etc.
#' @param mergecov Numeric threshold for merging of segment mean value.
#'
#' @return A data frame or tibble with merged segments and updated segment information.
#'
#' @details
#' This function uses \code{\link{MergeSegCheck}} to determine if two adjacent segments should be merged.
#'
#' @importFrom dplyr arrange
#' @export
MergeSegRow <- function(df, mergecov) {

  if(nrow(df) > 1 ){
    i <- 1
    while ( i < ( nrow(df) ) ) {

      cur_row <- df[i,]
      next_row <- df[ i+1,]

      if (  MergeSegCheck(cur_row = cur_row, next_row = next_row, mergecov) ) {
        new_df <- data.frame(
          Sample = cur_row$Sample,
          Chromosome = cur_row$Chromosome,
          Start = cur_row$Start,
          End = next_row$End,
          Num_Probes = cur_row$Num_Probes + next_row$Num_Probes,
          Call = ifelse( cur_row$Num_Probes > next_row$Num_Probes, cur_row$Call, next_row$Call),
          Segment_Mean = ifelse( cur_row$Num_Probes > next_row$Num_Probes,
                                 cur_row$Segment_Mean,  next_row$Segment_Mean),
          size = next_row$End - cur_row$Start,
          Count = ifelse( cur_row$Num_Probes > next_row$Num_Probes, cur_row$Count, next_row$Count),
          Baseline_cov = ifelse( cur_row$Num_Probes > next_row$Num_Probes, cur_row$Baseline_cov, next_row$Baseline_cov),
          gatk_gender = cur_row$gatk_gender,
          pipeline_gender = cur_row$pipeline_gender,
          Segment_Mean_raw = ifelse( cur_row$Num_Probes > next_row$Num_Probes, cur_row$Segment_Mean_raw, next_row$Segment_Mean_raw)
        )
        df <- rbind(new_df, df[-c(i,i+1),])
        df <- df %>% dplyr::arrange(by=Chromosome,Start,na.last = T)
        i <- 1
      }else{i <- i+1}
    }
  }
  return(df)
}



