#' Find bad tbpB sequences
#'
#' @param Org_id Organism to query: GAS, PNEUMO or GONO
#'
#' @details Scans through MG-MAST downloaded file tbpB.fasta to identify erroneous alleles
#' that do not start CGTCTGAA.
#' Creates a bad_tbpB.csv file will bad allele ids that is manually combined with bad_porB_alleles.csv to create
#' the main bad_alleles.csv list which is the combined bad tbpB and porB allele ids.
#' This combined list is used to filter 100% BLAST matches during NG-MAST sequence typing.
#'
#'
#'
#' @return A table frame containing the results of the query
#' @export


tbpB_scan <- function(Org_id) {

#Code-----------------------------------------####
library(plyr)
library(dplyr)
library(stringr)
#---------------------------------------------
  setwd("L:\\GC_PAHO\\Whole_Genome_Sequencing\\WGS_Typing\\GONO\\NGMAST_R\\temp")

  FileName <- "L:\\GC_PAHO\\Whole_Genome_Sequencing\\WGS_Typing\\GONO\\NGMAST_R\\allele_lkup_dna\\tbpB.fasta"
  con <- file(FileName, open="r")
  linn <- readLines(con)
  close(con)

  # scan for each new line that does not start with CGTCTGAA,
  # copy allele number to bad_alleles table,
  # or change to "XXXXX"
  bad_tpbB_alleles.df <-data.frame(Allele = character(), stringsAsFactors = FALSE)
  #bad_porB_alleles.df <- read.csv("L:\\GC_PAHO\\Whole_Genome_Sequencing\\WGS_Typing\\GONO\\NGMAST_R\\temp\\bad_porB_alleles.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

  for (i in seq(1, length(linn), 2))
  {
    if (str_detect(linn[i], ">"))
    {
      AlleleLine <- unlist(linn[i])
      AlleleParts <- strsplit(AlleleLine, ">")
      AlleleLine2 <- unlist(AlleleParts)
      Allele_id <- AlleleLine2[2]
      Sequence_Line <- unlist(linn[i+1])
      Start_Line <- substr(Sequence_Line, 1, 8)
      if (Start_Line != "CGTCTGAA")
      {
        cat(Allele_id, " Bad Sequence\n")
        bad_allele.df <- data_frame(Allele = Allele_id)
        bad_tpbB_alleles.df <- rbind(bad_tpbB_alleles.df, bad_allele.df)
      }
    }
  }

  #bad_alleles.df <- rbind(bad_tpbB_alleles.df, bad_porB_alleles.df)
  write.csv(bad_tbpB_alleles.df, "L:\\GC_PAHO\\Whole_Genome_Sequencing\\WGS_Typing\\GONO\\NGMAST_R\\temp\\bad_tbpB_alleles.csv", row.names = F)
  #write.csv(bad_alleles.df, "L:\\GC_PAHO\\Whole_Genome_Sequencing\\WGS_Typing\\GONO\\NGMAST_R\\temp\\bad_alleles.csv", row.names = F)



cat("\n\nDONE! ... Output - bad_alleles.csv created.")


return(bad_alleles.df)

}
