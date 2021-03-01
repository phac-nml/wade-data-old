#' NG-MAST MLST pipeline for WGS assemblies
#'
#' Takes Organism, Sample Number, Locus, and queries a contig.fasta file vs. the NG-MAST Blast database
#' @param Org_id Organism to query: GONO
#' @param SampleNo Sample number (or "list" to use a list of sample numbers) associated with contig.fasta file
#' @param locus The locus to query, (or "list" to use a list of alleles)
#' @param curr_work_dir Start up directory from pipeline project to locate system file structure
#' @return A table frame containing the results of the query
#' @export


# Org_id <- "GONO"
# SampleNo <- "40561"
# locus <- "list"
# curr_work_dir <- "L:\\GC_PAHO\\Whole_Genome_Sequencing\\WGS_Typing\\"


NGMAST_pipeline <- function(Org_id, SampleNo, locus, curr_work_dir) {
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

  #------------------------------------------------------------------------------------------------------------

  Lkup_Dir <- paste(system_dir, "GONO\\NGMAST_R\\allele_lkup_dna\\", sep = "")
  Tmp_Dir <- paste(system_dir, "GONO\\NGMAST_R\\temp\\", sep = "")
  Loci_List <- paste(Tmp_Dir, "loci.csv", sep = "")
  Profiles <- paste(Tmp_Dir, "profiles.csv", sep = "")

  Variable <- NA
  Allele <- ""
  AlleleNo <- ""

  #--------------------------------------------------------------------------------------Setup locus list table ####
  LocusList.df <- as_tibble(read.csv(Loci_List, header = TRUE, sep = ",", stringsAsFactors = FALSE))
  if(locus != "list")
  {
    LocusList.df <- filter(LocusList.df, LocusName == locus)
  }
  Size.df <- dim(LocusList.df)
  NumLoci <- Size.df[1]

  #wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww  Index BLAST lookup files
  # for(q in 1L:NumLoci)
  # {
  #   CurrLocus <- as.character(LocusList.df[q,1])
  #   LocusLkupDNA <- paste(Lkup_Dir, CurrLocus, ".fasta", sep = "")
  #   if(file.exists(LocusLkupDNA))
  #   {
  #     #BlastFormatCommand <- paste("formatdb -i ", LocusLkupDNA, " -p F", sep = "")
  #     BlastFormatCommand <- paste("makeblastdb -in ", LocusLkupDNA, " -dbtype nucl", sep = "")
  #     #try(system(BlastFormatCommand))
  #     shell(BlastFormatCommand)
  #   }
  #
  # }
  #wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww

  cat("\n\n", "NG-MAST \n", "Sample No.", "\n", sep = "")

  #------------------------------------------------------------------------------------Set up Sample list table
  if(SampleNo == "list")
  {
    SampleList.df <- as_tibble(read.csv(SampList, header = TRUE, sep = ",", stringsAsFactors = FALSE))
  }else
  {
    SampleList.df <- tibble(SampleNo, Variable)
  }

  Size.df <- dim(SampleList.df)
  NumSamples <- Size.df[1]

  df.bad_alleles <- as_tibble(read.csv(paste(Tmp_Dir, "bad_alleles.csv",sep = ""),
                             header = TRUE, sep = ",", stringsAsFactors = FALSE))

  #----------------------------------------------------------------------------------------- Load sequence type profiles
  profiles.df <- as_tibble(read.csv(Profiles, header = TRUE, sep = ",", stringsAsFactors = FALSE))
  #------------------------------------------------------------------------------------------

  for (m in 1L:NumSamples)  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Loop through sample table
  {
    Allele <- ""
    AlleleNo <- ""

    CurrSampleNo <- as.character(SampleList.df[m, "SampleNo"])
    QueryFile <- paste(ContigsDir, CurrSampleNo, ".fasta", sep = "")
    DestFile <- paste(local_temp_dir, "queryfile.fasta", sep = "")
    if (!file.copy(QueryFile, DestFile, overwrite = T))
    {
      #stop("Sample number not found in contigs directory!")
      #return(SampleList.df)
      Allele <- "Sample Number Error"
      AlleleNo <- "Sample_Err"

    }


      for (n in 1L:NumLoci)  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< loop through loci table
      {
        CurrLocus <- as.character(LocusList.df[n, "LocusName"])
        CurrLocusLen <- as.integer(LocusList.df[n, "size"])
      if (AlleleNo != "Sample_Err")
      {

        LocusLkupDNA <- paste(Lkup_Dir, CurrLocus, ".fasta", sep = "")
        if (!(file.exists(LocusLkupDNA)))
        {
          cat("Locus error for: ", locus, ".... file not found!")
          sample_error.df <- tibble(Loucs_ID = locus, Output = "File Not Found")
          return(sample_error.df)
        }

        Blast_Out_File <- paste(local_temp_dir, "blastout.txt", sep = "")
        #BlastCommand <- paste("blastall -i ./temp/queryfile.fasta -d ", LocusLkupDNA, " -p blastn -o ./temp/blastout.txt -v 0 -b 10 -e 10e-175 -m 8 -F F")
        BlastCommand <- paste("blastn -query ", DestFile, " -db ", LocusLkupDNA, " -out ", Blast_Out_File, " -num_alignments 10 -evalue 10e-50 -outfmt 6")

        shell(BlastCommand)

        #FileName <- paste(Tmp_Dir, "blastout.txt", sep = "")

        info = file.info(Blast_Out_File)
        if(info$size == 0)
        {
          Allele <- "No gene present"
          AlleleNo <- "0"
        }else
        {
          df.blastout <- as_tibble(read.csv(Blast_Out_File, header = FALSE, sep = "\t", stringsAsFactors = FALSE))
          names(df.blastout) <- c("SampleNo", "Allele", "Ident", "Align", "Mismatches", "Gaps", "SampleStart", "SampleEnd", "AlleleStart", "AlleleEnd", "eValue", "bit")
          #arrange(df.blastout, desc(Ident)) #this was causing errors, don't think we need this sorted anyway?
          df.blastout100 <- filter(df.blastout, Ident == 100 & Align == CurrLocusLen)
          if (CurrLocus == "porB"){str_start <- 5}else{str_start <- 6}
          df.blastout100 <- arrange(df.blastout100, SampleStart)
          df.blastout100 <- anti_join(df.blastout100, df.bad_alleles, by = "Allele")


          df.blastout100$Allele2 <- as.numeric(substr(df.blastout100$Allele, str_start, 10))
          #df.blastout100 <- arrange(df.blastout100, Allele2)

          dfSize <- nrow(df.blastout100)
          # if (dfSize > 0)
          # {
          #   Allele <- df.blastout100$Allele[1]
          #   AlleleParts <- unlist(strsplit(Allele, "_"))
          #   AlleleNo <- AlleleParts[2]
          # }else
          # {
          #   Allele <- "Not Found"
          #   AlleleNo <- "?"
          # }
          Allele <- ""
          Allele1 <- ""
          Allele2 <- ""
          AlleleNo <- ""

          if (dfSize == 1)  #only 1 unique 100% match
          {
            Allele <- df.blastout100$Allele[1]
            AlleleParts <- unlist(strsplit(Allele, "_"))
            AlleleNo <- AlleleParts[2]
          }

          if (dfSize > 1) #more than 1 100% match
          {
            Allele1 <- df.blastout100$Allele[1]
            AlleleParts1 <- unlist(strsplit(Allele1, "_"))
            AlleleNo1 <- AlleleParts1[2]

            Allele2 <- df.blastout100$Allele[2]
            AlleleParts2 <- unlist(strsplit(Allele2, "_"))
            AlleleNo2 <- AlleleParts2[2]

            Allele <- paste(Allele1, "/", Allele2, sep = "")
            AlleleNo <- paste(AlleleNo1, "/", AlleleNo2, sep = "")
          }

          if (dfSize == 0)
          {
            Allele <- "Not Found"
            AlleleNo <- "?"
          }


        }

      }#end if not sample error

        #-----------------------if this is the first locus, add headers to output, else cbind next locus to output.
        if(n==1)
        {
          Output.df <- tibble(CurrSampleNo, AlleleNo)
          names(Output.df) <- c("Sample", CurrLocus)
        }else
        {
          LocusOutput.df <- tibble(AlleleNo)
          names(LocusOutput.df) <- c(CurrLocus)
          Output.df <- bind_cols(Output.df, LocusOutput.df)
        }



      }#end locus loop


    if(locus == "list")
    {
      #--------------------------------------------------- lookup profiles
      profile2.df <- tibble(profiles.df)
      p<-1L
      for (p in 1L:NumLoci)
      {
        if (!empty(profile2.df))
        {
          profile2.df <- filter(profile2.df, profile2.df[,p] == as.character(Output.df[1,p+1]))
        }
      }

      if (empty(profile2.df))
      {
        ST <- NA
        MLSTtype.df <- tibble(ST)
      }else
      {
        MLSTtype.df <- select(profile2.df, ST)
      }

      Output.df <- bind_cols(Output.df, MLSTtype.df)

    }else
    {
      MLSTtype.df <- tibble(Output.df[2])
    }



    #-------------------- if this is the first sample, copy output of first, else rowbind to add next sample profiles to output
    if(m==1)
    {
      SampleOutput.df <- tibble(Output.df)
    }else
    {
      SampleOutput.df <- bind_rows(SampleOutput.df, Output.df)
    }

    if (locus == "list")
    {
      cat(CurrSampleNo, "\t por-", Output.df$porB[1], "\t tbpB-", Output.df$tbpB[1], "\t ST-", MLSTtype.df$ST[1], "\n", sep = "")

    }else
    {
      cat(CurrSampleNo, "\t", locus, "-", as.character(MLSTtype.df[1]), "\n", sep = "")

    }

  } #end sample loop



  if (NumLoci == 1L)
  {
    #View(df.blastout)
    write.csv(df.blastout, paste(local_output_dir, "output_profile_GONO_NGMAST.csv", sep = ""), row.names = F)

    cat("\n\nDone! output_profile_GONO_NGMAST.csv created. \n\n\n")


    return(df.blastout)

  }else
  {
    #View(SampleOutput.df)
    SampleOutput2.df <- tibble(SampleOutput.df)
    SampleOutput2.df$porB <- paste("porB-", SampleOutput2.df$porB, sep = "")
    SampleOutput2.df$tbpB <- paste("tbpB-", SampleOutput2.df$tbpB, sep = "")
    SampleOutput2.df$ST <- paste("ST-", SampleOutput2.df$ST, sep = "")

    SampleOutput2_good.df <- filter(SampleOutput2.df, !is.na(ST))
    SampleOutput2_bad.df <- filter(SampleOutput2.df, is.na(ST))

    write.csv(SampleOutput2.df, paste(local_output_dir, "output_profile.csv", sep = ""), row.names = F)

    write.csv(SampleOutput2.df, paste(local_output_dir, "LabWareUpload_GONO_NGMAST.csv", sep = ""), quote = FALSE,  row.names = FALSE)
    write.csv(SampleOutput2_good.df, paste(local_output_dir, "LabWareUpload_GONO_NGMAST_good.csv", sep = ""), quote = FALSE,  row.names = FALSE)
    write.csv(SampleOutput2_bad.df, paste(local_output_dir, "LabWareUpload_GONO_NGMAST_bad.csv", sep = ""), quote = FALSE,  row.names = FALSE)


    write.csv(SampleOutput2.df, paste(local_output_dir, "output_profile_GONO_NGMAST.csv", sep = ""), quote = FALSE, row.names = FALSE)

    cat("\n\nDone! LabWareUpload_GONO_NGMAST.csv created. \n\n\n")

    return(SampleOutput2.df)

  }




} #end of function
