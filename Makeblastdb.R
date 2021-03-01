# Makeblastdb.R indexes Blast fasta databases

#' Index BLAST fasta datases of lookup data
#'
#' Takes Organism, Test, and Locus and create BLAST databases
#' @param Org_id Organism to query: GAS, PNEUMO or GONO
#' @param TestID Sample number associated with contig.fasta file
#' @param locus Locus_name to build blast database
#' @param curr_work_dir Start up directory from pipeline project to locate system file structure
#' @return A table frame containing the results of the query
#' @export

Index_pipeline <- function(Org_id, Test_id, LocusID, curr_work_dir){

  #--------------------------------------------------------------------------------------------------------

  #Org_id <- "PNEUMO"                     #GAS, PNEUMO or GONO
  #Test_id <- "AMR"              #AMR, TOXINS, VIRULENCE, NGSTAR use MASTER for 16S id
  #LocusID <- "pbp1a"                   #L:\\GC_PAHO\\Whole Genome Sequencing\\WGS Typing\\GAS\\Wamr_R\\temp\\loci.csv"

  #------------------------------------------------------------------------------------------------------------
  #------------------------------------------------------------------------------------------------------------
  # get directory structure
  curr_dir <- curr_work_dir
  dir_file <- paste(curr_work_dir, "DirectoryLocations.csv", sep="")
  #cat("\n directory string passed to function: ", curr_dir, "\n")
  Directories.df <- as_tibble(read.csv(dir_file, header = TRUE, sep = ",", stringsAsFactors = FALSE))
  Directories_org.df <- filter(Directories.df, OrgID == Org_id)
  local_dir <- Directories_org.df$LocalDir
  #setwd(local_dir)
  SampList <- paste(local_dir, "list.csv", sep = "")
  local_output_dir <- paste(local_dir, "Output\\", sep = "")
  local_temp_dir <- paste(local_dir, "temp\\", sep = "")
  system_dir <- Directories_org.df$SystemDir
  ContigsDir <- Directories_org.df$ContigsDir

  switch(Test_id,
         AMR={test_dir <- "Wamr_R\\"},
         AMR_2={test_dir <- "Wamr_R\\"},
         AMR_LW={test_dir <- "Wamr_R\\"},
         TOXINS={test_dir <- "Toxins_R\\"},
         EMM={test_dir <- "emm_R\\"},
         VIRULENCE={test_dir <- "Virulence_R\\"},
         SERO={test_dir <- "Serotype_R\\"},
         MASTER={test_dir <- "Master_Blaster_R\\"},
         MLST={test_dir <- "MLST_R\\"},
         rRNA23S={test_dir <- "NGSTAR_R\\"},
         NGMAST={test_dir <- "NGMAST_R\\"},
         NGSTAR={test_dir <- "NGSTAR_R\\"}
  )
  LkupDir <- paste(system_dir, Org_id, "\\", test_dir, sep = "")
  #------------------------------------------------------------------------------------------------------------

if(LocusID == "list")
{
  LocusList <- paste(LkupDir, "temp\\loci.csv", sep = "")
  LocusList.df <- as_tibble(read.csv(LocusList, header = TRUE, sep = ",", stringsAsFactors = FALSE))
}else
{
  LocusList.df <- tibble(LocusID)
}

SizeList <- dim(LocusList.df)
NumLoci <- SizeList[1]
q<-1
for(q in 1L:NumLoci)
{
  locus <- as.character(LocusList.df[q,1])
  LocusLkupDNA <- paste(LkupDir, "allele_lkup_dna\\", locus, ".fasta", sep = "")
  if(file.exists(LocusLkupDNA))
  {
    #BlastFormatCommand <- paste("formatdb -i ", LocusLkupDNA, " -p F", sep = "")
    BlastFormatCommand <- paste("makeblastdb -in ", LocusLkupDNA, " -dbtype nucl", sep = "")

    try(system(BlastFormatCommand))
  } #else {cat("\nBlast Error\n")}
}

cat("\n\nDONE.  BLAST databases indexed! \n")

done_signal.df <- tibble(Output = "Blast databases indexed.")
return(done_signal.df)


}
