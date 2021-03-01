#Rename Assemblies and move to main contigs folder
#R Studio scripting on local version of R Studio
#2018-05-25
#Walter Demczuk


#' Rename files in VCF_new to remove _ _ numbers from file names from Galaxy:
#'
#'
#' @param Org_id Organism to query: GAS, PNEUMO or GONO
#' @return A table frame containing the results of the query
#' @export
#'
#'


rename_vcfs <- function(Org_id) {

library(dplyr)
library(stringr)

switch(Org_id,
       GONO={
         ContigsDir <- "W:\\Projects\\Project_GC_WalterD\\MiSeq\\Gonorrhoea\\vcf_new"
         setwd(ContigsDir)
         WinCommand <- "ren *__*.vcf ?????.vcf"
         shell(WinCommand, intern = TRUE)

         shell("move /Y W:\\Projects\\Project_GC_WalterD\\MiSeq\\Gonorrhoea\\vcf_new\\*.vcf W:\\Projects\\Project_GC_WalterD\\MiSeq\\Gonorrhoea\\vcf\\")
       },
       PNEUMO={
         ContigsDir <- "W:\\Projects\\Project_GC_WalterD\\MiSeq\\Streptococcus\\Pneumo\\vcf_new"
         setwd(ContigsDir)
         WinCommand <- "ren *__*.vcf *P.vcf"

         shell(WinCommand, intern = TRUE)

         shell("move /Y W:\\Projects\\Project_GC_WalterD\\MiSeq\\Streptococcus\\Pneumo\\vcf_new\\*.vcf W:\\Projects\\Project_GC_WalterD\\MiSeq\\Streptococcus\\Pneumo\\vcf\\")
       },
       {
         Switch_entry <- "OTHERWISE DO NOTHING"
       }
)




}
