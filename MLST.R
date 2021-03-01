#' MLST pipeline for WGS assemblies
#'
#' Takes Organism, Sample Number, Locus, and a Variable at queries a contig.fasta file
#' @param Org_id Organism to query: GAS, PNEUMO or GONO
#' @param SampleNo Sample number associated with contig.fasta file
#' @param locus The locus to query, or enter list to use a list of alleles
#' @param curr_work_dir Start up directory from pipeline project to locate system file structure
#' @return A table frame containing the results of the query
#' @export

# Org_id <- "GONO"
# SampleNo <- "40911"
# locus <- "list"
# curr_work_dir <- "L:\\GC_PAHO\\Whole_Genome_Sequencing\\WGS_Typing\\"

MLST_pipeline <- function(Org_id, SampleNo, locus, curr_work_dir) {
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

  Lkup_Dir <- paste(system_dir, Org_id, "\\MLST_R\\allele_lkup_dna\\", sep = "")
  Tmp_Dir <- paste(system_dir, Org_id, "\\MLST_R\\temp\\", sep = "")
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

  cat("\n\n", "Sample No.", "\n", sep = "")

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

  #----------------------------------------------------------------------------------------- Load sequence type profiles
  profiles.df <-as_tibble(read.csv(Profiles, header = TRUE, sep = ",", stringsAsFactors = FALSE))
  #------------------------------------------------------------------------------------------
  m<-1L
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

      n<-1L
      for (n in 1L:NumLoci)  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< loop through loci table
      {
        CurrLocus <- as.character(LocusList.df[n, "LocusName"])
        CurrLocusLen <- as.integer(LocusList.df[n, "size"])
      if (AlleleNo != "Sample_Err")
      {

        LocusLkupDNA <- paste(Lkup_Dir, CurrLocus, ".fasta", sep = "")
        if (!(file.exists(LocusLkupDNA)))
        {
          cat("Locus error for: ", CurrLocus, ".fasta file not found!", sep = "")
          sample_error.df <- tibble(Loucs_ID = CurrLocus, Output = "File Not Found")
          return(sample_error.df)
        }

        Blast_Out_File <- paste(local_temp_dir, "blastout.txt", sep = "")
        #BlastCommand <- paste("blastall -i ./temp/queryfile.fasta -d ", LocusLkupDNA, " -p blastn -o ./temp/blastout.txt -v 0 -b 10 -e 10e-175 -m 8 -F F")
        BlastCommand <- paste("blastn -query ", DestFile, " -db ", LocusLkupDNA, " -out ", Blast_Out_File, " -num_alignments 10 -evalue 10e-175 -outfmt 6")

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
          if ((Org_id == "PNEUMO" & CurrLocus == "xpt") | (Org_id == "GBS")) #xpt has variable lengths in lookup table seep xpt_113
          {
            df.blastout100 <- filter(df.blastout, Ident == 100)
          }else
          {
            df.blastout100 <- filter(df.blastout, Ident == 100 & Align == CurrLocusLen)
          }
          #df.blastout100 <- filter(df.blastout, Ident == 100 & Align == CurrLocusLen)
          dfSize <- nrow(df.blastout100)
          if (dfSize > 0)
          {
            Allele <- df.blastout100$Allele[1]
            AlleleParts <- unlist(strsplit(Allele, "_"))
            AlleleNo <- AlleleParts[2]
          }else
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
      #snvTable1.df <- unite(snvTable.df, SNVprofile, 1:27, sep = "", remove = TRUE)
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
      cat(CurrSampleNo, "\t", "ST-", MLSTtype.df$ST[1], "\n", sep = "")

    }else
    {
      cat(CurrSampleNo, "\t", locus, "-", MLSTtype.df$ST[1], "\n", sep = "")

    }

  } #end sample loop


  cat("Done!", sep ="")

  if (NumLoci == 1L)
  {
    #View(df.blastout)
    write.csv(df.blastout, paste(local_output_dir, "output_profile_MLST.csv", sep = ""), row.names = F)

    cat("\n\nDone!\n\n\n")


    return(df.blastout)

  }else
  {
    #View(SampleOutput.df)
    #write.csv(SampleOutput.df, "./temp/output_profile.csv", row.names = F)
    #SampleOutput.df <- read.csv("C:\\WGS_Typing\\Output\\LabWareUpload_GONO_MLST.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE) # this line is only for debugging
    SampleOutput_good.df <- filter(SampleOutput.df, !is.na(ST))
    SampleOutput_bad.df <- filter(SampleOutput.df, is.na(ST))


    write.csv(SampleOutput.df, paste(local_output_dir, "LabWareUpload_", Org_id, "_MLST.csv", sep = ""), quote = FALSE,  row.names = FALSE)
    write.csv(SampleOutput_good.df, paste(local_output_dir, "LabWareUpload_", Org_id, "_MLST_good.csv", sep = ""), quote = FALSE,  row.names = FALSE)
    write.csv(SampleOutput_bad.df, paste(local_output_dir, "LabWareUpload_", Org_id, "_MLST_bad.csv", sep = ""), quote = FALSE,  row.names = FALSE)

    write.csv(SampleOutput.df, paste(local_output_dir, "output_profile_MLST.csv", sep = ""), quote = FALSE, row.names = FALSE)

    cat("\n\nDone!\n\n\n")

    return(SampleOutput.df)

  }




} #end of function
