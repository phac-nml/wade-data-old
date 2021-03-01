#' run MasterBlastR first, then this script outputs csv to LabWareUpload_GAS_TOXINS.csv in Output folder
#'
#' Then run this analysis to combine data the full amr profile to upload to LabWare.
#'
#'
#' @param Org_id Organism to query: GAS, PNEUMO or GONO
#' @param curr_work_dir Start up directory from pipeline project to locate system file structure
#' @return A table frame containing the results of the query
#' @export
#'
#'

# Org_id <- "PNEUMO"
# curr_work_dir <- "L:\\GC_PAHO\\Whole_Genome_Sequencing\\WGS_Typing\\"


labware_pneumo_virulence <- function(Org_id, curr_work_dir) {

  #------------------------------------------------------------------------------------------------------------
  # get directory structure
  dir_file <- paste(curr_work_dir, "DirectoryLocations.csv", sep="")
  Directories.df <- as_tibble(read.csv(dir_file, header = TRUE, sep = ",", stringsAsFactors = FALSE))
  Directories_org.df <- filter(Directories.df, OrgID == Org_id)
  local_dir <- Directories_org.df$LocalDir
  SampList <- paste(local_dir, "list.csv", sep = "")
  local_output_dir <- paste(local_dir, "Output\\", sep = "")
  local_temp_dir <- paste(local_dir, "temp\\", sep = "")
  system_dir <- Directories_org.df$SystemDir
  ContigsDir <- Directories_org.df$ContigsDir


Output.df <- as_tibble(read.csv(paste(local_output_dir, "output_profile_PNEUMO_VIRULENCE.csv", sep = ""),
                                header = TRUE, sep = ",", stringsAsFactors = FALSE))

Output.df$SampleProfile[Output.df$SampleProfile == ""] <- "Error"

Size.df <- dim(Output.df)
NumSamples <- Size.df[1]
NumLoci <- ((Size.df[2]-2) / 8)

if (NumLoci > 1)
{
LabWare.df <- tibble(
                                 Output.df$SampleNo,
                                 Output.df$PI1_result,
                                 Output.df$PI2_result,
                                 Output.df$lytA_result,
                                 Output.df$ply_result,
                                 Output.df$rpoB_allele,
                                 Output.df$rpoB_mutations,
                                 Output.df$rpoB_PctWT,
                                 Output.df$SampleProfile
                          )
headers <- list("SampleNo", "PI1", "PI2", "lytA", "ply", "rpoB", "rpoB_ID", "rpoB_match", "Profile")

names(LabWare.df) <- headers
#View(LabWare.df)
write.csv(LabWare.df, paste(local_output_dir, "LabWareUpload_PNEUMO_VIRULENCE.csv", sep = ""), quote = FALSE, row.names = FALSE)
}else  # else if a single loci was run return the original output.csv
{
  LabWare.df <- Output.df
}


return(LabWare.df)

}

