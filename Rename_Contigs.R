#Rename Assemblies and move to main contigs folder
#R Studio scripting on local version of R Studio
#2018-05-25
#Walter Demczuk

#' Rename files in contigs_new to remove _ _ numbers from file names from Galaxy:
#'
#'
#' @param Org_id Organism to query: GAS, PNEUMO or GONO
#' @return A table frame containing the results of the query
#' @export
#'
#'

rename_contigs <- function(Org_id) {

library(dplyr)
library(stringr)

switch(Org_id,
       GAS={
            ContigsDir <- "W:\\Projects\\Project_GC_WalterD\\MiSeq\\Streptococcus\\GAS\\contigs_new"
            WinCommand <- "ren *__* *A.*"

            },


        PNEUMO={
                ContigsDir <- "W:\\Projects\\Project_GC_WalterD\\MiSeq\\Streptococcus\\PNEUMO\\contigs_new"
                WinCommand <- "ren *__* *P.*"

                },

       GROUP_CG={
                ContigsDir <- "W:\\Projects\\Project_GC_WalterD\\MiSeq\\Streptococcus\\GROUP_CG\\contigs_new"
                WinCommand <- "ren *__* *CG.*"

                },
       GONO={
         ContigsDir <- "W:\\Projects\\Project_GC_WalterD\\MiSeq\\Gonorrhoea\\contigs_new"
         WinCommand <- "ren *__*.fasta ?????.fasta"

       },

       {
         Switch_entry <- "NOTHING"
       }
)


setwd(ContigsDir)
shell(WinCommand, intern = TRUE)


shell("move /Y *.fasta ..\\contigs\\")

}
