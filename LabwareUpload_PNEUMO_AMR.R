#' Run AMR first, then run the 23S allele counts,
#' Then run this analysis to combine data from AMR, 23S rRNA
#' to prepare full amr profile to upload to LabWare.
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


labware_pneumo_amr <- function(Org_id, curr_work_dir) {


  #------------------------------------------------------------------------------------------------------------
  # get directory structure
  dir_file <- paste(curr_work_dir, "DirectoryLocations.csv", sep="")
  Directories.df <- read.csv(dir_file, header = TRUE, sep = ",")
  Directories.df <- as_tibble(Directories.df)
  Directories_org.df <- filter(Directories.df, OrgID == Org_id)
  local_dir <- Directories_org.df$LocalDir
  SampList <- paste(local_dir, "list.csv", sep = "")
  local_output_dir <- paste(local_dir, "Output\\", sep = "")
  local_temp_dir <- paste(local_dir, "temp\\", sep = "")
  system_dir <- Directories_org.df$SystemDir
  ContigsDir <- Directories_org.df$ContigsDir


cat(paste("\n\n", "Formatting output table...", "\n\n", sep = ""))

AMR_Output.df <- as_tibble(read.csv(paste(local_output_dir, "output_profile_PNEUMO_AMR.csv", sep = ""),
                          header = TRUE, sep = ",", stringsAsFactors = FALSE))
#lw_Output.df <- as_tibble(AMR_Output.df)
Size.df <- dim(AMR_Output.df)
NumSamples <- Size.df[1]
NumLoci <- ((Size.df[2]-3) / 6)

if (NumLoci >1)
{

rRNA23S_Output.df <- as_tibble(read.csv(paste(local_output_dir, "output_profile_23S.csv", sep = ""),
                              header = TRUE, sep = ",", stringsAsFactors = FALSE))
#rRNA23S_Output.df <- as_tibble(rRNA23S_Output.df)
Combined_Output.df <- full_join(AMR_Output.df, rRNA23S_Output.df, by = "SampleNo")

Size.df <- dim(Combined_Output.df)
NumSamples <- Size.df[1]

pen_mic.df <- as_tibble(read.csv(paste(system_dir, "PNEUMO\\Wamr_R\\temp\\inc_mic_penicillin.csv", sep = ""),
                       header = TRUE, sep = ",", stringsAsFactors = FALSE))
cro_mic.df <- as_tibble(read.csv(paste(system_dir, "PNEUMO\\Wamr_R\\temp\\inc_mic_ceftriaxone.csv", sep = ""),
                       header = TRUE, sep = ",", stringsAsFactors = FALSE))
cfm_mic.df <- as_tibble(read.csv(paste(system_dir, "PNEUMO\\Wamr_R\\temp\\inc_mic_cefuroxime.csv", sep = ""),
                       header = TRUE, sep = ",", stringsAsFactors = FALSE))
ery_mic.df <- as_tibble(read.csv(paste(system_dir, "PNEUMO\\Wamr_R\\temp\\inc_mic_erythromycin.csv", sep = ""),
                       header = TRUE, sep = ",", stringsAsFactors = FALSE))
cla_mic.df <- as_tibble(read.csv(paste(system_dir, "PNEUMO\\Wamr_R\\temp\\inc_mic_clarithromycin.csv", sep = ""),
                       header = TRUE, sep = ",", stringsAsFactors = FALSE))
cli_mic.df <- as_tibble(read.csv(paste(system_dir, "PNEUMO\\Wamr_R\\temp\\inc_mic_clindamycin.csv", sep = ""),
                       header = TRUE, sep = ",", stringsAsFactors = FALSE))
lev_mic.df <- as_tibble(read.csv(paste(system_dir, "PNEUMO\\Wamr_R\\temp\\inc_mic_levofloxacin.csv", sep = ""),
                       header = TRUE, sep = ",", stringsAsFactors = FALSE))
mox_mic.df <- as_tibble(read.csv(paste(system_dir, "PNEUMO\\Wamr_R\\temp\\inc_mic_moxifloxacin.csv", sep = ""),
                       header = TRUE, sep = ",", stringsAsFactors = FALSE))
tet_mic.df <- as_tibble(read.csv(paste(system_dir, "PNEUMO\\Wamr_R\\temp\\inc_mic_tetracycline.csv", sep = ""),
                       header = TRUE, sep = ",", stringsAsFactors = FALSE))
dox_mic.df <- as_tibble(read.csv(paste(system_dir, "PNEUMO\\Wamr_R\\temp\\inc_mic_doxycycline.csv", sep = ""),
                       header = TRUE, sep = ",", stringsAsFactors = FALSE))
sxt_mic.df <- as_tibble(read.csv(paste(system_dir, "PNEUMO\\Wamr_R\\temp\\inc_mic_trimeth_sulfa.csv", sep = ""),
                       header = TRUE, sep = ",", stringsAsFactors = FALSE))
chl_mic.df <- as_tibble(read.csv(paste(system_dir, "PNEUMO\\Wamr_R\\temp\\inc_mic_chloramphenicol.csv", sep = ""),
                       header = TRUE, sep = ",", stringsAsFactors = FALSE))
sepr <- " - "

m <- 1
ErrorFound <- FALSE

for (m in 1L:NumSamples)  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Build Molecular Profile
{
  molec_profile <- NA
  lw_comments <- ""

  lw_CurrSampleNo <- as.character(Combined_Output.df[m, "SampleNo"])

  lw_ermB <- as.character(Combined_Output.df[m, "ermB_result"])
  lw_ermB_allele <- as.character(Combined_Output.df[m, "ermB_allele"])  #check for NF(Not Found found in allele lkup)

  if (lw_ermB != "Sample_Err")
  {
  lw_ermB_mut <- as.character(Combined_Output.df[m, "ermB_mutations"])  #check for ermB(S) -- Susceptible

    if (lw_ermB == "POS")
      {
      molec_profile <- "ermB"

      if (lw_ermB_mut == "S")
        {
        lw_ermB <- "NEG"
        molec_profile <- "ermB(S)"
        }
      if (lw_ermB_allele == "NF")
      {
        lw_ermB <- "NF"
        lw_comments <- paste(lw_comments, "ermB curation error.", sep="")
      }

      }

    #this next test needs to be outside of the above loop to work for some reason???
    if (lw_ermB_mut == "S" | is.na(lw_ermB_mut))
    {
    lw_ermB <- "NEG"
    }
  }


  lw_ermTR <- as.character(Combined_Output.df[m, "ermTR_result"])
      if (lw_ermTR == "POS")
      {
        if (is.na(molec_profile))
        {molec_profile <- "ermTR"} else
        {molec_profile <- paste(molec_profile, sepr, "ermTR", sep = "")}

      }

  lw_mefAE <- as.character(Combined_Output.df[m, "mefAE_result"])
  if (lw_mefAE == "POS")
  {
    if (is.na(molec_profile))
    {molec_profile <- "mefAE"} else
    {molec_profile <- paste(molec_profile, sepr, "mefAE", sep = "")}

  }


  lw_23S_prof <- NA

  if (is.na(Combined_Output.df[m, "A2059G"]))
  {
    lw_23S_A2059G <- 0L
    lw_23S_prof <- "23S Err"
    lw_comments <- paste(lw_comments, "23S VCF file not found.", sep="")
  }else
    {lw_23S_A2059G <- as.integer(Combined_Output.df[m, "A2059G"])}

  if (is.na(Combined_Output.df[m, "C2611T"]))
  {
    lw_23S_C2611T <- 0L
    lw_23S_prof <- "23S Err"
  }else
  {lw_23S_C2611T <- as.integer(Combined_Output.df[m, "C2611T"])}



  #if (is.na(lw_23S_A2059G)) {lw_23S_A2059G <- "Err"}
  #if (is.na(lw_23S_C2611T)) {lw_23S_C2611T <- "Err"}

  if (lw_23S_A2059G != 0L)
  {
    if (is.na(lw_23S_prof))
    {lw_23S_prof <- paste("23S rRNA A2059G ", lw_23S_A2059G, "/4", sep = "")} else
    {lw_23S_prof <- paste(ls_23S_prof, "/A2059G ", lw_23S_A2059G, "/4",  sep = "")}
  }

  if (lw_23S_C2611T != 0L)
  {
    if (is.na(lw_23S_prof))
    {lw_23S_prof <- paste("23S rRNA C2611T ", lw_23S_C2611T, "/4", sep = "")} else
    {lw_23S_prof <- paste(lw_23S_prof, "/C2611T ", lw_23S_C2611T, "/4",  sep = "")}
  }

  if (!is.na(lw_23S_prof))
  {
    if (is.na(molec_profile))
    {molec_profile <- lw_23S_prof} else
    {molec_profile <- paste(molec_profile, sepr, lw_23S_prof, sep = "")}
  }else(lw_23S_prof <- "")



  lw_tetM <- as.character(Combined_Output.df[m, "tetM_result"])
  if (lw_tetM == "POS")
  {
    if (is.na(molec_profile))
    {molec_profile <- "tetM"} else
    {molec_profile <- paste(molec_profile, sepr, "tetM", sep = "")}
  }

  lw_tetO <- as.character(Combined_Output.df[m, "tetO_result"])
  if (lw_tetO == "POS")
  {
    if (is.na(molec_profile))
    {molec_profile <- "tetO"} else
    {molec_profile <- paste(molec_profile, sepr, "tetO", sep = "")}
  }

  lw_cat <- as.character(Combined_Output.df[m, "cat_result"])
  if (lw_cat == "POS")
  {
    if (is.na(molec_profile))
    {molec_profile <- "cat"} else
    {molec_profile <- paste(molec_profile, sepr, "cat", sep = "")}
  }

  lw_folA <-   as.character(Combined_Output.df[m, "folA_mutations"])
  lw_folA_result <- as.character(Combined_Output.df[m, "folA_result"])
  if (lw_folA == "???" | lw_folA == "x" | lw_folA == "" | is.na(lw_folA))
  {
    lw_folA <- "Err"
    lw_comments <- paste(lw_comments, "folA curation error.", sep="")
  }

  if (lw_folA != "WT")
  {
    if (is.na(molec_profile))
    {molec_profile <- paste("folA ", lw_folA, sep = "")}else
    {molec_profile <- paste(molec_profile, sepr, "folA ", lw_folA, sep = "")}
  }

  lw_folP <-   as.character(Combined_Output.df[m, "folP_mutations"])
  lw_folP_result <- as.character(Combined_Output.df[m, "folP_result"])
  if (lw_folP == "???" | lw_folP == "x" | lw_folP == "" | is.na(lw_folP))
  {
    lw_folP <- "Err"
    lw_comments <- paste(lw_comments, "folP curation error.", sep="")
  }

  if (lw_folP != "WT")
  {
    if (is.na(molec_profile))
    {molec_profile <- paste("folP ", lw_folP, sep = "")}else
    {molec_profile <- paste(molec_profile, sepr, "folP ", lw_folP, sep = "")}
  }

  lw_pbp1a <-   as.character(Combined_Output.df[m, "pbp1a_mutations"])
  lw_pbp1a_result <-   as.character(Combined_Output.df[m, "pbp1a_result"])

  if (lw_pbp1a == "???" | lw_pbp1a == "x" | lw_pbp1a == "" | is.na(lw_pbp1a))
  {
    lw_pbp1a <- "Err/Err/Err/Err"
    lw_comments <- paste(lw_comments, "pbp1a curation error.", sep="")
    }

  if (lw_pbp1a != "WT/WT/WT/WT")
  {
    if (is.na(molec_profile))
    {molec_profile <- paste("pbp1a ", lw_pbp1a, sep = "")}else
    {molec_profile <- paste(molec_profile, sepr, "pbp1a ", lw_pbp1a, sep = "")}
  }

  lw_pbp2b <-   as.character(Combined_Output.df[m, "pbp2b_mutations"])
  lw_pbp2b_result <-   as.character(Combined_Output.df[m, "pbp2b_result"])
  if (lw_pbp2b == "???" | lw_pbp2b == "x" | lw_pbp2b == "" | is.na(lw_pbp2b))
  {
    lw_pbp2b <- "Err/Err/Err/Err"
    lw_comments <- paste(lw_comments, "pbp2b curation error.", sep="")
    }

  if (lw_pbp2b != "WT/WT/WT/WT")
  {
    if (is.na(molec_profile))
    {molec_profile <- paste("pbp2b ", lw_pbp2b, sep = "")}else
    {molec_profile <- paste(molec_profile, sepr, "pbp2b ", lw_pbp2b, sep = "")}
  }

  lw_pbp2x <-   as.character(Combined_Output.df[m, "pbp2x_mutations"])
  lw_pbp2x_result <-   as.character(Combined_Output.df[m, "pbp2x_result"])
  if (lw_pbp2x == "???" | lw_pbp2x == "x" | lw_pbp2x == "" | is.na(lw_pbp2x))
  {
    lw_pbp2x <- "Err/Err/Err/Err"
    lw_comments <- paste(lw_comments, "pbp2x curation error.", sep="")
    }

  if (lw_pbp2x != "WT/WT/WT/WT")
  {
    if (is.na(molec_profile))
    {molec_profile <- paste("pbp2x ", lw_pbp2x, sep = "")}else
    {molec_profile <- paste(molec_profile, sepr, "pbp2x ", lw_pbp2x, sep = "")}
  }


  lw_pbp_profile <- paste("PBP" ,as.character(Combined_Output.df[m, "pbp1a_allele"]),
                        as.character(Combined_Output.df[m, "pbp2b_allele"]),
                        as.character(Combined_Output.df[m, "pbp2x_allele"]), sep = ":")


  lw_gyrA <- as.character(Combined_Output.df[m, "gyrA_mutations"])
  lw_gyrA_result <- as.character(Combined_Output.df[m, "gyrA_result"])
  if (lw_gyrA == "???" | lw_gyrA == "x" | lw_gyrA == "" | is.na(lw_gyrA))
  {
    lw_gyrA <- "Err"
    lw_comments <- paste(lw_comments, "gyrA curation error.", sep="")
    }

  if (lw_gyrA != "WT")
  {
    if (is.na(molec_profile))
    {molec_profile <- paste("gyrA ", lw_gyrA, sep = "")}else
    {molec_profile <- paste(molec_profile, sepr, "gyrA ", lw_gyrA, sep = "")}
  }


  lw_parC <-   as.character(Combined_Output.df[m, "parC_mutations"])
  lw_parC_result <- as.character(Combined_Output.df[m, "parC_result"])


  if (lw_parC == "???" | lw_parC == "x" | lw_parC == "" | is.na(lw_parC))
  {
    lw_parC <- "Err/Err/Err"
    lw_comments <- paste(lw_comments, "parC curation error.", sep="")
    }

  #if (lw_parC != "WT/WT/WT")
  #{

    lw_parC_parts <- unlist(strsplit(lw_parC, "/"))
    lw_parC_S79 <-  lw_parC_parts[1]
    lw_parC_D83 <- lw_parC_parts[2]
    lw_parC_N91 <- lw_parC_parts[3]


    lw_parC_prof <- NA


    if (lw_parC_S79 != "WT")
    {
      if (is.na(lw_parC_prof))
      {lw_parC_prof <- paste("parC ", lw_parC_S79, sep = "")} else
      {lw_parC_prof <- paste(lw_parC_prof, "/", lw_parC_S79, sep = "")}
    }
    if (lw_parC_D83 != "WT")
    {
      if (is.na(lw_parC_prof))
      {lw_parC_prof <- paste("parC ", lw_parC_D83, sep = "")} else
      {lw_parC_prof <- paste(lw_parC_prof, "/", lw_parC_D83, sep = "")}
    }
    if (lw_parC_N91 != "WT")
    {
      if (is.na(lw_parC_prof))
      {lw_parC_prof <- paste("parC ", lw_parC_N91, sep = "")} else
      {lw_parC_prof <- paste(lw_parC_prof, "/", lw_parC_N91, sep = "")}
    }

    if (!is.na(lw_parC_prof))
    {
      if (is.na(molec_profile))
      {molec_profile <- lw_parC_prof} else
      {molec_profile <- paste(molec_profile, sepr, lw_parC_prof, sep = "")}
    }

  #}

  if ((molec_profile == "") | is.na(molec_profile)) {molec_profile <- "Wild Type"}




  #-----------------------------------------------------Fluorquinolone mics and interpretatons

  lw_gyrA_value_S81F <- 0
  lw_gyrA_value_S81Y <- 0
  lw_gyrA_value_S81L <- 0
  lw_parC_value_S79_any <- 0
  lw_parC_value_D83_any <- 0

  if (lw_gyrA == "S81F")
  {lw_gyrA_value_S81F <- 1}
  if (lw_gyrA == "S81Y")
  {lw_gyrA_value_S81Y <- 1}
  if (lw_gyrA == "S81L")
  {lw_gyrA_value_S81L <- 1}

  if (lw_parC_S79 != "WT") #any mutation
  {lw_parC_value_S79_any <- 1}
  if (lw_parC_D83 != "WT") #any mutation
  {lw_parC_value_D83_any <- 1}


  lev_mic_inc <-  round(-0.218 +
                          (lw_gyrA_value_S81F * 2.028)+
                          (lw_gyrA_value_S81Y * 1.564)+
                          (lw_gyrA_value_S81L * 3.564)+
                          (lw_parC_value_S79_any * 1.654)+
                          (lw_parC_value_D83_any * 0.834)
                          )

  lev_mic <- round( (2^lev_mic_inc), digits = 3)
  lev_mic2 <- paste(lev_mic, " ug/ml", sep = "")

  if (lev_mic < 2)
  {
    lev_mic2 <- paste("<= 1 ug/ml", sep = "")
  }else if (lev_mic > 16)
  {
    lev_mic2 <- paste(">= 32 ug/ml", sep = "")
  }

  lev_interp <- lev_mic.df$Interp[lev_mic.df$MIC == lev_mic]


  mox_mic_inc <-  round(-2.819 +
                          (lw_gyrA_value_S81F * 3.130)+
                          (lw_gyrA_value_S81Y * 3.907)+
                          (lw_gyrA_value_S81L * 4.907)+
                          (lw_parC_value_S79_any * 0.911)
                          )

  mox_mic <- round( (2^mox_mic_inc), digits = 3)
  mox_mic2 <- paste(mox_mic, " ug/ml", sep = "")

  if (mox_mic < 0.25)
  {
    mox_mic2 <- paste("<= 0.125 ug/ml", sep = "")
  }else if (mox_mic > 16)
  {
    mox_mic2 <- paste(">= 32 ug/ml", sep = "")
  }

  mox_interp <- mox_mic.df$Interp[mox_mic.df$MIC == mox_mic]

  if (lw_gyrA == "Err" | lw_parC == "Err/Err/Err" | lw_gyrA == "short" | lw_parC == "short")
  {
    lev_mic2 <- "Error"
    lev_interp <- "Error"
    mox_mic2 <- "Error"
    mox_interp <- "Error"
  }

  if (lw_gyrA_result=="NEG")
  {
    lw_gyrA <- "no gene"
    lw_comments <- paste(lw_comments, "gyrA missing.", sep="")
    }
  if (lw_parC_result=="NEG")
    {
    lw_parC <- "no gene"
    lw_parC_S79 <- "no gene"
    lw_parC_D83 <- "no gene"
    lw_parC_N91 <- "no gene"
    lw_comments <- paste(lw_comments, "parC missing.", sep="")
    }


  #---------------------------------------------------------------------  Penicillin  MICs and interpretaions

  lw_pbp1a_parts <- unlist(strsplit(lw_pbp1a, "/")) #get motif regions
  lw_pbp1a_1 <-  lw_pbp1a_parts[1]
  lw_pbp1a_2 <-  lw_pbp1a_parts[2]
  lw_pbp1a_3 <-  lw_pbp1a_parts[3]
  lw_pbp1a_4 <-  lw_pbp1a_parts[4]

  lw_pbp2b_parts <- unlist(strsplit(lw_pbp2b, "/")) #get motif regions
  lw_pbp2b_1 <-  lw_pbp2b_parts[1]
  lw_pbp2b_2 <-  lw_pbp2b_parts[2]
  lw_pbp2b_3 <-  lw_pbp2b_parts[3]
  lw_pbp2b_4 <-  lw_pbp2b_parts[4]


  lw_pbp2x_parts <- unlist(strsplit(lw_pbp2x, "/"))
  lw_pbp2x_1 <-  lw_pbp2x_parts[1]
  lw_pbp2x_2 <-  lw_pbp2x_parts[2]
  lw_pbp2x_3 <-  lw_pbp2x_parts[3]
  lw_pbp2x_4 <-  lw_pbp2x_parts[4]

  lw_pbp1a_motif1_value_any <- 0
  lw_pbp1a_motif1_value_SAMK <- 0
  lw_pbp1a_motif1_value_SSMK <- 0
  lw_pbp1a_motif2_value_any <- 0
  lw_pbp1a_motif3_value_any <- 0
  lw_pbp1a_motif4_value_any <- 0

  lw_pbp2b_motif1_value_any <- 0
  lw_pbp2b_motif2_value_any <- 0
  lw_pbp2b_motif3_value_any <- 0
  lw_pbp2b_motif4_value_any <- 0

  lw_pbp2x_motif1_value_any <- 0
  lw_pbp2x_motif1_value_SAFK <- 0
  lw_pbp2x_motif1_value_other <- 0
  lw_pbp2x_motif2_value_any <- 0
  lw_pbp2x_motif3_value_any <- 0
  lw_pbp2x_motif3_value_EDT <- 0
  lw_pbp2x_motif3_value_KEA <- 0
  lw_pbp2x_motif4_value_VKSG <- 0
  lw_pbp2x_motif4_value_any <- 0


  if (lw_pbp1a_1 != "WT") #any mutation
  {lw_pbp1a_motif1_value_any <- 1}

  if (lw_pbp1a_1 == "SAMK")
  {lw_pbp1a_motif1_value_SAMK <- 1}

  if (lw_pbp1a_1 == "SSMK")
  {lw_pbp1a_motif1_value_SSMK <- 1}

  if (lw_pbp1a_4 != "WT") #any mutation
  {lw_pbp1a_motif4_value_any <- 1}

  if (lw_pbp2b_2 != "WT") #any mutation
  {lw_pbp2b_motif2_value_any <- 1}

  if (lw_pbp2b_3 != "WT") #any mutation
  {lw_pbp2b_motif3_value_any <- 1}

  if (lw_pbp2x_1 == "SAFK")
  {
    lw_pbp2x_motif1_value_SAFK <- 1
  }else
    {
      if (lw_pbp2x_1 != "WT")
      {lw_pbp2x_motif1_value_other <- 1}
    }

  if (lw_pbp2x_3 == "EDT")
  {lw_pbp2x_motif3_value_EDT <- 1}

  if (lw_pbp2x_3 == "KEA")
  {lw_pbp2x_motif3_value_EDT <- 1}

  if (lw_pbp2x_4 == "VKSG")
  {lw_pbp2x_motif4_value_VKSG <- 1}

  # Intercept	-4.609850206
  # pbp1a_1_Any	1.546643583
  # pbp1a_4_Any	0.949074251
  # pbp2b_2_Any	1.202527621
  # pbp2b_3_Any	0.356298924
  # pbp2x_1_SAFK	1.625871441
  # pbp2x_3_EDT	1.547626017
  # pbp2x_3_KEA	0.675981078
  # pbp2x_4_VKSG	0.753453309


  pen_mic_inc <-  round(-4.6099 +
    (lw_pbp1a_motif1_value_any * 1.5466)+
    (lw_pbp1a_motif4_value_any * 0.9491)+
    (lw_pbp2b_motif2_value_any * 1.2025)+
    (lw_pbp2b_motif3_value_any * 0.3563)+
    (lw_pbp2x_motif1_value_SAFK * 1.6259)+
    (lw_pbp2x_motif3_value_EDT * 1.5476)+
    (lw_pbp2x_motif3_value_KEA * 0.6760)+
    (lw_pbp2x_motif4_value_VKSG * 0.7536))

  pen_mic <- round( (2^pen_mic_inc), digits = 3)
  pen_mic2 <- paste(pen_mic, " ug/ml", sep = "")

  if (pen_mic < 0.032)
  {
    pen_mic2 <- paste("<= 0.03 ug/ml", sep = "")
  }else if (pen_mic > 3)
  {
    pen_mic2 <- paste(">= 4 ug/ml", sep = "")
  }


  pen_interp <- pen_mic.df$Interp[pen_mic.df$MIC == pen_mic]



  cro_mic_inc <-  round(-2.938 +
                          (lw_pbp1a_motif1_value_any * 1.425)+
                          (lw_pbp2x_motif1_value_SAFK * 2.677)+
                          (lw_pbp2x_motif1_value_other * 0.161)+
                          (lw_pbp2x_motif4_value_VKSG * 1.085))

  #cro_mic <- as.character(round( (2^cro_mic_inc), digits = 3))
  cro_mic <- round( (2^cro_mic_inc), digits = 3)
  cro_mic2 <- paste(cro_mic, " ug/ml", sep = "")

  if (cro_mic < 0.13)
  {
    cro_mic2 <- paste("<= 0.125 ug/ml", sep = "")
  }else if (cro_mic > 3)
  {
    cro_mic2 <- paste(">= 4 ug/ml", sep = "")
  }


  cro_interp <- cro_mic.df$Interp[cro_mic.df$MIC == cro_mic]

  cfm_mic_inc <-  round(-1.018 +
                          (lw_pbp1a_motif1_value_SAMK * 1.509)+
                          (lw_pbp1a_motif1_value_SSMK * 2.170)+
                          (lw_pbp2x_motif1_value_SAFK * 2.322)+
                          (lw_pbp2x_motif1_value_other * 0.256)+
                          (lw_pbp2x_motif4_value_VKSG * 1.026))

  #cfm_mic <- as.character(round( (2^cfm_mic_inc), digits = 3))
  cfm_mic <- round( (2^cfm_mic_inc), digits = 3)

  cfm_mic2 <- paste(cfm_mic, " ug/ml", sep = "")

  if (cfm_mic < 0.51)
  {
    cfm_mic2 <- paste("<= 0.5 ug/ml", sep = "")
  }else if (cfm_mic > 15)
  {
    cfm_mic2 <- paste(">= 16 ug/ml", sep = "")
  }


  cfm_interp <- cfm_mic.df$Interp[cfm_mic.df$MIC == cfm_mic]

  if ((lw_pbp1a == "Err/Err/Err") | (lw_pbp2b == "Err/Err/Err") | (lw_pbp2x == "Err/Err/Err"))
  {
    pen_mic2 <- "Error"
    cro_mic2 <- "Error"
    cfm_mic2 <- "Error"
    pen_interp <- "Error"
    cro_interp <- "Error"
    cfm_interp <- "Error"
  }

  if(lw_pbp1a_result=="NEG")
  {
    lw_pbp1a <- "no gene"
    lw_comments <- paste(lw_comments, "pbp1a missing.", sep="")
    }
  if(lw_pbp2b_result=="NEG")
  {
    lw_pbp2b <- "no gene"
    lw_comments <- paste(lw_comments, "pbp2b missing.", sep="")
    }
  if(lw_pbp2x_result=="NEG")
  {
    lw_pbp2x <- "no gene"
    lw_comments <- paste(lw_comments, "pbp2x missing.", sep="")
    }


  #-----------------------------------------------------------------Trimethoprim/Sulfamethoxazole
  lw_folA_value <- 0
  lw_folP_value <- 0

  if (lw_folA == "I100L")
  {lw_folA_value <- 1}
  if (lw_folP == "INS")
  {lw_folP_value <- 1}


  sxt_mic_inc <- round(-2.265 +
                         (lw_folA_value * 2.113)+
                         (lw_folP_value * 2.668)
                       )
  sxt_mic <-  round( (2^sxt_mic_inc), digits = 3)
  sxt_mic2 <- paste(sxt_mic, " ug/ml", sep = "")

  if (sxt_mic < 0.51)
  {
    sxt_mic2 <- paste("<= 0.5/9.5 ug/ml", sep = "")
  }else if (sxt_mic > 3)
  {
    sxt_mic2 <- paste(">= 4/76 ug/ml", sep = "")
  }

  sxt_interp <- sxt_mic.df$Interp[sxt_mic.df$MIC == sxt_mic]

  if ((lw_folA == "Err") | (lw_folP == "Err"))
  {
    sxt_mic2 <- "Error"
    sxt_interp <- "Error"
  }

  if(lw_folA_result=="NEG")
  {
    lw_folA <- "no gene"
    lw_comments <- paste(lw_comments, "folA missing.", sep="")
    }
  if(lw_folP_result=="NEG")
  {
    lw_folP <- "no gene"
    lw_comments <- paste(lw_comments, "folP missing.", sep="")
    }


  #------------------------------------------------------------------------- MACROLIDEs     MICs and Interpretations

  # lw_23S_A2059G set above
  # lw_23s_C2611T set above
  lw_ermB_value <- 0
  lw_ermTR_value <- 0
  lw_mefAE_value <- 0

  if (lw_ermB == "POS")
  {lw_ermB_value <- 1} # if ermB(S), then lw_ermB is set to NEG above
  if (lw_ermTR == "POS")
  {lw_ermTR_value <- 1}
  if (lw_mefAE == "POS")
  {lw_mefAE_value <- 1}

  # if both ermB and mef are positive, only use ermB value to calculate
  if ((lw_ermB == "POS") && (lw_mefAE == "POS"))
  {
    lw_mefAE_value <- 0
  }

  ery_mic_inc <-  round(-2.975 +
                          (lw_23S_A2059G * 1.993)+
                          (lw_ermB_value * 7.680)+
                          (lw_mefAE_value * 4.808)
                        )

  ery_mic <- round( (2^ery_mic_inc), digits = 5)
  ery_mic2 <- paste(ery_mic, " ug/ml", sep = "")

  if (ery_mic < 0.25)
  {
    ery_mic2 <- paste("<= 0.125 ug/ml", sep = "")
  }else if (ery_mic > 16)
  {
    ery_mic2 <- paste(">= 32 ug/ml", sep = "")
  }


  ery_interp <- ery_mic.df$Interp[ery_mic.df$MIC == ery_mic]



  cla_mic_inc <-  round(-4.984 +
                          (lw_23S_A2059G * 1.819)+
                          (lw_23S_C2611T * 1.246)+
                          (lw_ermB_value * 10.820)+
                          (lw_mefAE_value * 6.064)
  )

  cla_mic <- round( (2^cla_mic_inc), digits = 5)
  cla_mic2 <- paste(cla_mic, " ug/ml", sep = "")

  if (cla_mic < 0.6)
  {
    cla_mic2 <- paste("<= 0.3 ug/ml", sep = "")
  }else if (cla_mic > 32)
  {
    cla_mic2 <- paste(">= 64 ug/ml", sep = "")
  }

  cla_interp <- cla_mic.df$Interp[cla_mic.df$MIC == cla_mic]


  cli_mic_inc <-  round(-2.8145 +
                          (lw_23S_A2059G * 0.456)+
                          (lw_ermB_value * 9.048)
  )

  cli_mic <- round( (2^cli_mic_inc), digits = 5)

  cli_mic2 <- paste(cli_mic, " ug/ml", sep = "")

  if (cli_mic < 0.25)
  {
    cli_mic2 <- paste("<= 0.125 ug/ml", sep = "")
  }else if (cli_mic > 32)
  {
    cli_mic2 <- paste(">= 64 ug/ml", sep = "")
  }


  cli_interp <- cli_mic.df$Interp[cli_mic.df$MIC == cli_mic]
  if (lw_ermTR == "POS")
    {
    cli_interp <- "Inducible"
  }


  if (lw_ermB == "Sample_Err")
  {
    ery_mic2 <- "Sample_Err"
    ery_interp <- "Sample_Err"
    cli_mic2 <- "Sample_Err"
    cli_interp <- "Sample_Err"
    cla_mic2 <- "Sample_Err"
    cla_interp <- "Sample_Err"
    azi_mic2 <- "Sample_Err"
    azi_interp <- "Sample_Err"
  }

  if (lw_ermB == "NF")
  {
    ery_mic2 <- "Error"
    ery_interp <- "Curation Error"
    cli_mic2 <- "Error"
    cli_interp <- "Curation Error"
    cla_mic2 <- "Error"
    cla_interp <- "Curation Error"
    azi_mic2 <- "Error"
    azi_interp <- "Curation Error"
  }


  if (lw_23S_prof == "23S Err")
  {
    lw_23S_A2059G <- 99L
    lw_23S_C2611T <- 99L
    ery_mic2 <- "Error"
    ery_interp <- "No VCF"
    cli_mic2 <- "Error"
    cli_interp <- "No VCF"
    cla_mic2 <- "Error"
    cla_interp <- "No VCF"
    azi_mic2 <- "Error"
    azi_interp <- "No VCF"
  }
  #-----------------------------------------------------------------------------------------


  chl_mic <- 2
  if (lw_cat == "POS"){chl_mic <- 16} #list in increasing mic order

  if (chl_mic == 2)
  {
    chl_mic2 <- paste("<= 2 ug/ml", sep = "")
  }else if (chl_mic == 16)
  {
    chl_mic2 <- paste(">= 16 ug/ml", sep = "")
  }

  chl_interp <- chl_mic.df$Interp[chl_mic.df$MIC == chl_mic]

  if (lw_cat == "Sample_Err")
  {
    chl_mic2 <- "Sample_Err"
    chl_interp <- "Sample_Err"
  }

  #-----------------------------------------------------------------------------------------


  tet_mic <- 1
  if (lw_tetM == "POS"){tet_mic <- 16} #list in increasing mic order

  if (tet_mic == 1)
  {
    tet_mic2 <- paste("<= 1 ug/ml", sep = "")
  }else if (tet_mic == 16)
  {
    tet_mic2 <- paste(">= 16 ug/ml", sep = "")
  }

  tet_interp <- tet_mic.df$Interp[tet_mic.df$MIC == tet_mic]



  dox_mic <- 0.25
  if (lw_tetM == "POS"){dox_mic <- 4} #list in increasing mic order

  if (dox_mic == 0.25)
  {
    dox_mic2 <- paste("<= 0.25 ug/ml", sep = "")
  }else if (dox_mic == 4)
  {
    dox_mic2 <- paste(">= 4 ug/ml", sep = "")
  }


  dox_interp <- dox_mic.df$Interp[dox_mic.df$MIC == dox_mic]

  if (lw_tetM == "Sample_Err")
  {
    tet_mic2 <- "Sample_Err"
    tet_interp <- "Sample_Err"
    dox_mic2 <- "Sample_Err"
    dox_interp <- "Sample_Err"

  }


  #--------------------------------------------------------------  MAKE AMR PROFILE

  amr_profile <- "Susceptible"


  if ( (str_detect(molec_profile, "folA Err")) | (str_detect(molec_profile, "folP Err")) )
  {
  amr_profile <- "Error"
  sxt_interp <- "Error"
  }


  if ( (str_detect(molec_profile, "gyrA Err")) | (str_detect(molec_profile, "parC Err")) )
  {lev_interp <- "Error"
  mox_interp <- "Error"
  amr_profile <- "Error"
  }

  if (str_detect(molec_profile, "Err"))
  {
    amr_profile <- "Error"
  }


  if (lw_ermB == "Sample_Err")
    {
    amr_profile <- "Sample Error"
    }


  if (amr_profile == "Susceptible")
  {

    #amr_profile <- "Susceptible"
    sepr2 <- "/"

    if (pen_interp == "Intermediate")
    {
      if (amr_profile == "Susceptible")
      {amr_profile <- "PEN-I"} else
      {amr_profile <- paste(amr_profile, sepr2, "PEN-I", sep = "")}
    }
    if (pen_interp == "Resistant")
    {
      if (amr_profile == "Susceptible")
      {amr_profile <- "PEN-R"} else
      {amr_profile <- paste(amr_profile, sepr2, "PEN-R", sep = "")}
    }


    if (cro_interp == "Intermediate")
    {
      if (amr_profile == "Susceptible")
      {amr_profile <- "CRO-I"} else
      {amr_profile <- paste(amr_profile, sepr2, "CRO-I", sep = "")}
    }
    if (cro_interp == "Resistant")
    {
      if (amr_profile == "Susceptible")
      {amr_profile <- "CRO-R"} else
      {amr_profile <- paste(amr_profile, sepr2, "CRO-R", sep = "")}
    }

    if (cfm_interp == "Intermediate")
    {
      if (amr_profile == "Susceptible")
      {amr_profile <- "CFM-I"} else
      {amr_profile <- paste(amr_profile, sepr2, "CFM-I", sep = "")}
    }
    if (cfm_interp == "Resistant")
    {
      if (amr_profile == "Susceptible")
      {amr_profile <- "CFM-R"} else
      {amr_profile <- paste(amr_profile, sepr2, "CFM-R", sep = "")}
    }


    if (ery_interp == "Resistant")
    {
      if (amr_profile == "Susceptible")
      {amr_profile <- "ERY-R"} else
      {amr_profile <- paste(amr_profile, sepr2, "ERY-R", sep = "")}
    }
    if (cla_interp == "Resistant")
    {
      if (amr_profile == "Susceptible")
      {amr_profile <- "CLA-R"} else
      {amr_profile <- paste(amr_profile, sepr2, "CLA-R", sep = "")}
    }

    if (cli_interp == "Resistant")
    {
      if (amr_profile == "Susceptible")
      {amr_profile <- "CLI-R"} else
      {amr_profile <- paste(amr_profile, sepr2, "CLI-R", sep = "")}
    }
    if (cli_interp == "Inducible")
    {
      if (amr_profile == "Susceptible")
      {amr_profile <- "CLI-Ind"} else
      {amr_profile <- paste(amr_profile, sepr2, "CLI-Ind", sep = "")}
    }

    if (chl_interp == "Resistant")
    {
      if (amr_profile == "Susceptible")
      {amr_profile <- "CHL-R"} else
      {amr_profile <- paste(amr_profile, sepr2, "CHL-R", sep = "")}
    }



    if (lev_interp == "Resistant")
    {
      if (amr_profile == "Susceptible")
      {amr_profile <- "LEV-R"} else
      {amr_profile <- paste(amr_profile, sepr2, "LEV-R", sep = "")}
    }

    if (mox_interp == "Resistant")
    {
      if (amr_profile == "Susceptible")
      {amr_profile <- "MOX-R"} else
      {amr_profile <- paste(amr_profile, sepr2, "MOX-R", sep = "")}
    }

    if (tet_interp == "Resistant")
    {
      if (amr_profile == "Susceptible")
      {amr_profile <- "TET-R"} else
      {amr_profile <- paste(amr_profile, sepr2, "TET-R", sep = "")}
    }
    if (dox_interp == "Resistant")
    {
      if (amr_profile == "Susceptible")
      {amr_profile <- "DOX-R"} else
      {amr_profile <- paste(amr_profile, sepr2, "DOX-R", sep = "")}
    }

    if (sxt_interp == "Resistant")
    {
      if (amr_profile == "Susceptible")
      {amr_profile <- "SXT-R"} else
      {amr_profile <- paste(amr_profile, sepr2, "SXT-R", sep = "")}
    }
    if (sxt_interp == "Intermediate")
    {
      if (amr_profile == "Susceptible")
      {amr_profile <- "SXT-I"} else
      {amr_profile <- paste(amr_profile, sepr2, "SXT-I", sep = "")}
    }

  }

  if (lw_ermB == "Sample_Err") {lw_comments <- "Contig file not found"}

  #---------------------------------------------------------------------New LabWare Uploader structure
  sample_data.df <- tibble(lw_CurrSampleNo, lw_pbp1a, lw_pbp2b, lw_pbp2x,
                               lw_23S_A2059G, lw_23S_C2611T,  lw_ermB, lw_ermTR, lw_mefAE,
                               lw_folA, lw_folP,
                               lw_gyrA, lw_parC,
                               lw_tetM, lw_tetO, lw_cat, molec_profile,
                               pen_mic2, pen_interp,
                               cro_mic2, cro_interp,
                               cfm_mic2, cfm_interp,
                               azi_mic2 = "Not Done", azi_interp = "Not Done",
                               ery_mic2, ery_interp,
                               cla_mic2, cla_interp,
                               cli_mic2, cli_interp,
                               chl_mic2, chl_interp,
                               lev_mic2, lev_interp,
                               mox_mic2, mox_interp,
                               tet_mic2, tet_interp,
                               dox_mic2, dox_interp,
                               sxt_mic2, sxt_interp,
                               amr_profile,
                               lw_comments)
  #----------------------------------------------------------------------



  if(m==1)
  {
    lw_Output.df <- tibble(sample_data.df)
  }else
  {
    lw_Output.df <- bind_rows(lw_Output.df, sample_data.df)
  }



}#end samples loop

#View(lw_Output.df)

lw_Output_bad.df <- filter(lw_Output.df, molec_profile == "Err")
#lw_Output_good.df <- filter(lw_Output.df, molec_profile != "Err")

write.csv(lw_Output.df, paste(local_output_dir, "LabWareUpload_PNEUMO_AMR.csv", sep = ""), quote = FALSE, row.names = FALSE)
#write.csv(lw_Output_good.df, paste(local_output_dir, "LabWareUpload_PNEUMO_AMR_good.csv", sep = ""), quote = FALSE, row.names = FALSE)
write.csv(lw_Output_bad.df, paste(local_output_dir, "LabWareUpload_PNEUMO_AMR_bad.csv", sep = ""), quote = FALSE, row.names = FALSE)



}#end if more than 1 locus

 #else {
  #       lw_Output.df <- AMR_Output.df
  #       }

return(lw_Output.df)

}  #end function



