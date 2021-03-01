#' Run AMR first
#' Then run this analysis to combine data the full amr profile to upload to LabWare.
#'
#'
#' @param Org_id Organism to query: GAS, PNEUMO or GONO
#' @param curr_work_dir Start up directory from pipeline project to locate system file structure
#' @return A table frame containing the results of the query
#' @export
#'
#'

#Org_id <- "GBS"

labware_gbs_amr <- function(Org_id, curr_work_dir) {

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

Output.df <- as_tibble(read.csv(paste(local_output_dir, "output_profile_GBS_AMR.csv", sep = ""),
                                header = TRUE, sep = ",", stringsAsFactors = FALSE))

Size.df <- dim(Output.df)
NumSamples <- Size.df[1]
NumLoci <- ((Size.df[2]-2) / 8)

if (NumLoci > 1)
{

sepr <- " - "

m <- 1

for (m in 1L:NumSamples)  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
{
  molec_profile <- ""
  lw_comments <- ""

  lw_CurrSampleNo <- as.character(Output.df[m, "SampleNo"])

  lw_ermA <- as.character(Output.df[m, "ermA_result"])
  if (lw_ermA == "POS")
  {
    if (molec_profile == "")
    {molec_profile <- "ermA"}else
      {molec_profile <- paste(molec_profile, sepr, "ermA", sep = "")}
  }

  lw_ermB <- as.character(Output.df[m, "ermB_result"])
  if (lw_ermB == "POS")
  {
    if (molec_profile == "")
    {molec_profile <- "ermB"}else
    {molec_profile <- paste(molec_profile, sepr, "ermB", sep = "")}
  }

  lw_ermT <- as.character(Output.df[m, "ermT_result"])
  if (lw_ermT == "POS")
  {
    if (molec_profile == "")
    {molec_profile <- "ermT"}else
    {molec_profile <- paste(molec_profile, sepr, "ermT", sep = "")}
  }

  lw_mefAE <- as.character(Output.df[m, "mefAE_result"])
  if (lw_mefAE == "POS")
  {
    if (molec_profile == "")
    {molec_profile <- "mefAE"}else
    {molec_profile <- paste(molec_profile, sepr, "mefAE", sep = "")}
  }


  lw_gyrA <-   as.character(Output.df[m, "gyrA_mutations"])


  if (lw_gyrA == "???" | lw_gyrA == "x" | lw_gyrA == "" | is.na(lw_gyrA))
  {
    lw_gyrA <- "Err"
    lw_comments <- paste(lw_comments, "gyrA curation error.", sep = "")
  }


  if (lw_gyrA != "WT")
  {
    if (molec_profile == "")
    {molec_profile <- paste("gyrA ", lw_gyrA, sep = "")}else
    {molec_profile <- paste(molec_profile, sepr, "gyrA ", lw_gyrA, sep = "")}
  }

  lw_gyrA_result <-   as.character(Output.df[m, "gyrA_result"])
  if (lw_gyrA_result == "NEG")
  {
    lw_gyrA <- "no gene"
    lw_comments <- paste(lw_comments, "gyrA missing.", sep = "")
  }

  lw_parC <-   as.character(Output.df[m, "parC_mutations"])


  if (lw_parC == "???" | lw_parC == "x" | lw_parC == "" | is.na(lw_parC))
  {
    lw_parC <- "Err/Err/Err"
    lw_comments <- paste(lw_comments, "parC curation error.", sep = "")
  }


  #if (lw_parC != "WT")
  #{

  lw_parC_parts <- unlist(strsplit(lw_parC, "/"))
  lw_parC_S79 <-  lw_parC_parts[1]
  lw_parC_S80 <- lw_parC_parts[2]
  lw_parC_D83 <- lw_parC_parts[3]


  lw_parC_prof <- NA

  #if (lw_parC_D86 != "WT")
  #{
  #  if (is.na(lw_parC_prof))
  #  {lw_parC_prof <- paste("parC ", lw_parC_D86)} else
  #  {lw_parC_prof <- paste(lw_parC_prof, "/", lw_parC_D86, sep = "")}
  #}

  #lw_parC_prof <- ""
  if (lw_parC_S79 != "WT")
  {
    if (is.na(lw_parC_prof))
    {lw_parC_prof <- paste("parC ", lw_parC_S79, sep = "")} else
    {lw_parC_prof <- paste(lw_parC_prof, "/", lw_parC_S79, sep = "")}
  }

  if (lw_parC_S80 != "WT")
  {
    if (is.na(lw_parC_prof))
    {lw_parC_prof <- paste("parC ", lw_parC_S80, sep = "")} else
    {lw_parC_prof <- paste(lw_parC_prof, "/", lw_parC_S80, sep = "")}
  }

  if (lw_parC_D83 != "WT")
  {
    if (is.na(lw_parC_prof))
    {lw_parC_prof <- paste("parC ", lw_parC_D83, sep = "")} else
    {lw_parC_prof <- paste(lw_parC_prof, "/", lw_parC_D83, sep = "")}
  }

  if (!is.na(lw_parC_prof))
  {
    if (molec_profile == "")
    {molec_profile <- lw_parC_prof} else
    {molec_profile <- paste(molec_profile, sepr, lw_parC_prof, sep = "")}
  }

  #}

  lw_parC_result <-   as.character(Output.df[m, "parC_result"])
  if (lw_parC_result == "NEG")
  {
    lw_parC <- "no gene"
    lw_comments <- paste(lw_comments, "parC missing.", sep = "")
  }

  lw_tetM <- as.character(Output.df[m, "tetM_result"])
  if (lw_tetM == "POS")
  {
    if (molec_profile == "")
    {molec_profile <- "tetM"}else
    {molec_profile <- paste(molec_profile, sepr, "tetM", sep = "")}
  }

  lw_tetO <- as.character(Output.df[m, "tetO_result"])
  if (lw_tetO == "POS")
  {
    if (molec_profile == "")
    {molec_profile <- "tetO"}else
    {molec_profile <- paste(molec_profile, sepr, "tetO", sep = "")}
  }

  lw_cat <- as.character(Output.df[m, "cat_result"])
  if (lw_cat == "POS")
  {
    if (molec_profile == "")
    {molec_profile <- "cat"}else
    {molec_profile <- paste(molec_profile, sepr, "cat", sep = "")}
  }


  # lw_dfrF <- as.character(Output.df[m, "dfrF_result"])
  #
  # if (lw_dfrF == "POS")
  # {
  #   if (molec_profile == "")
  #   {molec_profile <- "dfrF"}else
  #   {molec_profile <- paste(molec_profile, sepr, "dfrF", sep = "")}
  # }
  # lw_dfrG <- as.character(Output.df[m, "dfrG_result"])
  # if (lw_dfrG == "POS")
  # {
  #   if (molec_profile == "")
  #   {molec_profile <- "dfrG"}else
  #   {molec_profile <- paste(molec_profile, sepr, "dfrG", sep = "")}
  # }

  # lw_folA <-   as.character(Output.df[m, "folA_mutations"])
  # if (lw_folA == "???" | lw_folA == "x" | lw_folA == "" | is.na(lw_folA)) {lw_folA <- "Err"}
  #
  #
  # if (lw_folA != "WT")
  # {
  #   if (molec_profile == "")
  #   {molec_profile <- paste("folA ", lw_folA, sep = "")}else
  #   {molec_profile <- paste(molec_profile, sepr, "folA ", lw_folA, sep = "")}
  # }
  #
  # lw_folP <-   as.character(Output.df[m, "folP_mutations"])
  # if (lw_folP == "???" | lw_folP == "x" | lw_folP == "" | is.na(lw_folP)) {lw_folP <- "Err"}
  #
  #
  #
  # if (lw_folP != "WT")
  # {
  #   if (molec_profile == "")
  #   {molec_profile <- paste("folP ", lw_folP, sep = "")}else
  #   {molec_profile <- paste(molec_profile, sepr, "folP ", lw_folP, sep = "")}
  # }


  lw_pbp2x <-   as.character(Output.df[m, "pbp2x_mutations"])


  if (lw_pbp2x == "???" | lw_pbp2x == "x" | lw_pbp2x == "" | is.na(lw_pbp2x))
  {
    lw_pbp2x <- "Err/Err"
    lw_comments <- paste(lw_comments, "pbp2x curation error.", sep = "")
  }

  if (lw_pbp2x != "WT/WT")
  {

    lw_pbp2x_parts <- unlist(strsplit(lw_pbp2x, "/"))
    lw_pbp2x_400 <-  lw_parC_parts[1]
    lw_pbp2x_552 <- lw_parC_parts[2]

    lw_pbp2x_prof <- NA

    if (lw_pbp2x_400 != "WT")
    {
      if (is.na(lw_pbp2x_prof))
      {lw_pbp2x_prof <- paste("pbp2x ", lw_pbp2x_400, sep = "")} else
      {lw_pbp2x_prof <- paste(lw_pbp2x_prof, "/", lw_pbp2x_400, sep = "")}
    }

    if (lw_pbp2x_552 != "WT")
    {
      if (is.na(lw_pbp2x_prof))
      {lw_pbp2x_prof <- paste("pbp2x ", lw_pbp2x_552, sep = "")} else
      {lw_pbp2x_prof <- paste(lw_pbp2x_prof, "/", lw_pbp2x_552, sep = "")}
    }


    if (!is.na(lw_pbp2x_prof))
    {
      if (molec_profile == "")
      {molec_profile <- lw_pbp2x_prof} else
      {molec_profile <- paste(molec_profile, sepr, lw_pbp2x_prof, sep = "")}
    }

  }


  lw_pbp2x_result <-   as.character(Output.df[m, "pbp2x_result"])
  if (lw_pbp2x_result == "NEG")
  {
    lw_pbp2x <- "no gene"
    lw_comments <- paste(lw_comments, "pbp2x missing.", sep = "")
  }


  if (molec_profile == "") {molec_profile <- "Wild Type"}


  #----------------------------------------------------------------------  INTERPRETATIONS

  amr_profile <- "Susceptible"

  if (str_detect(molec_profile, paste(c("ermA", "ermB", "ermT", "mefAE"),collapse = '|')))
  {
    ery <- "Resistant"
  }else{ery <- "Susceptible"}

  if (str_detect(molec_profile, "cat"))
  {
    chl <- "Resistant"
  }else{chl <- "Susceptible"}

  if (str_detect(molec_profile, "ermB"))
  {
    cli <- "Resistant"
  }else if (str_detect(molec_profile, paste(c("ermA", "ermT"),collapse = '|')))
  {
    cli <- "Inducible"
  } else {cli <-"Susceptible"}


  cip <- "Undetermined"
  lev <- "Undetermined"
  if ( (str_detect(molec_profile, "gyrA Err")) | (str_detect(molec_profile, "parC Err")) )
  {cip <- "Error"
  amr_profile <- "Error"
  }else
  {
    if (str_detect(molec_profile, paste(c("gyrA", "parC"),collapse = '|')))
    {
      cip <- "Resistant"
    }else{cip <- "Susceptible"}

    if (str_detect(molec_profile, paste(c("gyrA S81L", "S79F", "S79Y", "D83Y", "D83G"),collapse = '|')))
    {
      lev <- "Resistant"
    }else{lev <- "Susceptible"}
  }


  if (str_detect(molec_profile, paste(c("tetM", "tetO"),collapse = '|')))
  {
    tet <- "Resistant"
  }else{tet <- "Susceptible"}

  # if ( (str_detect(molec_profile, "folA Err")) | (str_detect(molec_profile, "folP Err")) )
  # {sxt <- "Error"
  # amr_profile <- "Error"
  # }else
  #
  # {
  #   if (str_detect(molec_profile, paste(c("dfrG", "dfrF", "folA", "folP"),collapse = '|')))
  #   {
  #
  #     if (str_detect(molec_profile, paste(c("dfrG", "dfrF"),collapse = '|')))
  #     {
  #     sxt <- "Resistant"
  #     } else if (str_detect(molec_profile, paste(c("folA", "folP"),collapse = '|')))
  #     {
  #       sxt <- "Intermediate"
  #       if ((str_detect(molec_profile, "folA")) & (str_detect(molec_profile, "folP")))
  #       {
  #       sxt <- "Resistant"
  #       }
  #     }
  #   }else{sxt <- "Susceptible"}
  #
  # }

  if (str_detect(molec_profile, "pbp2x"))
  {
    pen <- "Decreased Susceptiblity"
  }else{pen <- "Susceptible"}

  #--------------------------------------------------------------  MAKE AMR PROFILE
  if (amr_profile != "Error")
  {

  amr_profile <- "Susceptible"
  sepr2 <- "/"

  if (ery == "Resistant")
  {
    if (amr_profile == "Susceptible")
    {amr_profile <- "ERY-R"} else
    {amr_profile <- paste(amr_profile, sepr2, "ERY-R", sep = "")}
  }
  if (cli == "Resistant")
  {
    if (amr_profile == "Susceptible")
    {amr_profile <- "CLI-R"} else
    {amr_profile <- paste(amr_profile, sepr2, "CLI-R", sep = "")}
  }
  if (cli == "Inducible")
  {
    if (amr_profile == "Susceptible")
    {amr_profile <- "CLI-Ind"} else
    {amr_profile <- paste(amr_profile, sepr2, "CLI-Ind", sep = "")}
  }

  if (chl == "Resistant")
  {
    if (amr_profile == "Susceptible")
    {amr_profile <- "CHL-R"} else
    {amr_profile <- paste(amr_profile, sepr2, "CHL-R", sep = "")}
  }
  if (cip == "Resistant")
  {
    if (amr_profile == "Susceptible")
    {amr_profile <- "CIP-R"} else
    {amr_profile <- paste(amr_profile, sepr2, "CIP-R", sep = "")}
  }
  if (lev == "Resistant")
  {
    if (amr_profile == "Susceptible")
    {amr_profile <- "LEV-R"} else
    {amr_profile <- paste(amr_profile, sepr2, "LEV-R", sep = "")}
  }
  if (tet == "Resistant")
  {
    if (amr_profile == "Susceptible")
    {amr_profile <- "TET-R"} else
    {amr_profile <- paste(amr_profile, sepr2, "TET-R", sep = "")}
  }

  # if (sxt == "Resistant")
  # {
  #   if (amr_profile == "Susceptible")
  #   {amr_profile <- "SXT-R"} else
  #   {amr_profile <- paste(amr_profile, sepr2, "SXT-R", sep = "")}
  # }
  # if (sxt == "Intermediate")
  # {
  #   if (amr_profile == "Susceptible")
  #   {amr_profile <- "SXT-I"} else
  #   {amr_profile <- paste(amr_profile, sepr2, "SXT-I", sep = "")}
  # }

  if (pen == "Decreased Susceptiblity")
  {
    if (amr_profile == "Susceptible")
    {amr_profile <- "PEN-DS"} else
    {amr_profile <- paste(amr_profile, sepr2, "PEN-DS", sep = "")}
  }

  }

  if (lw_ermA == "Sample_Err")
  {
    lw_comments <- "Contig file not found"
       ery <- "Error"
       chl <- "Error"
       cip <- "Error"
       cli <- "Error"
       tet <- "Error"
       pen <- "Error"
  }



  #-------------------------------------------------------------- New Uploader structure
  LabWare_Sample.df <- tibble(
                                 lw_CurrSampleNo,
                                 lw_ermA,
                                 lw_ermB,
                                 lw_ermT,
                                 lw_mefAE,
                                 lw_gyrA,
                                 lw_parC,
                                 lw_tetM,
                                 lw_tetO,
                                 lw_cat,
                                 lw_pbp2x,
                                 molec_profile,
                                 ery,
                                 chl,
                                 cip,
                                 cli,
                                 tet,
                                 pen,
                                 amr_profile,
                                 lw_comments
                                 )


#--------------------------------------------------------------


  if(m==1)  #if first sample make one row profile table, otherwise add new row to table
  {
  LabWare.df <- tibble(LabWare_Sample.df)
  }else
  {
  LabWare.df <- rbind(LabWare.df, LabWare_Sample.df)
  }

}

lw_output_bad.df <- filter(LabWare.df, amr_profile == "Error")
#lw_output_good.df <- filter(LabWare.df, amr_profile != "Error")

write.csv(LabWare.df, paste(local_output_dir, "LabWareUpload_GBS_AMR.csv", sep = ""), quote = FALSE,  row.names = FALSE)
write.csv(lw_output_bad.df, paste(local_output_dir, "LabWareUpload_GBS_AMR_bad.csv", sep = ""), quote = FALSE,  row.names = FALSE)
#write.csv(lw_output_good.df, "LabWareUpload_GBS_AMR_good.csv", quote = FALSE,  row.names = FALSE)


}else  # else if a single loci was run return the original output.csv
{
  LabWare.df <- Output.df
}

return(LabWare.df)


}
