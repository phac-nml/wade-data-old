#WADE Shiny app

# Set up a local directory to ensure outputs and temporary files are stored locally and not on the network
# so that file sharing issues can be avoided.
# C:/WGS_Typing
# with two subdirectories:
# C:/WGS_Typing/Output
# C:/WGS_Typing/temp
# The multiple sample list file is located at: C:/WGS_Typing/list.csv
# list.csv must have the following structure:

# SampleNo	      Variable
# SC20-3242-B	    4 ug/ml
# SC20-3243-B	    8 ug/ml


#To install Biostrings copy and paste next 2 lines into R console
# Might need to upgrade R to version 4?

# https://bioconductor.org/install/

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.11")

# BiocManager::install("ggtree", "Biostrings")


#Download the latest "BLAST+" software from NCBI:
# https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download

# if you get errors with the makeblastdb program:
# got to Windows Settings and search for "Environmental Variables"
# opens the "System Properties" dialogue box, click the "Environmental Variables" button.
# in the 'User Variables for..." box, click "New..." button.
#   Variable Name: BLASTDB_LMDB_MAP_SIZE
#   Variable Value: 1000000
# note that the value is one million, 1,000,000

library(plyr)
library(dplyr)
library(tidyverse)
library(tidyselect)
library(stringr)
library(Biostrings)
library(shiny)
library(DT)
library(WGSpipeline)


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< set location of directory structure mapping table
#curr_work_dir <- "C:\\WGS_Typing\\"
#curr_work_dir <- "C:\\WADE\\"
curr_work_dir <- "L:\\GC_PAHO\\Whole_Genome_Sequencing\\WGS_Typing\\"


# Define UI ---------------------------------------------------------------------------------------
ui <- fluidPage(

  img(src = "Logo.png"),

  #titlePanel(strong("Strep/STI WGS Molecular Typing Pipeline")),

  sidebarLayout(position = "left",
    sidebarPanel(
      #helpText("Molecular typing profiles from WGS contig assemblies."),

      selectInput("Org",
                  label = h3("Choose an Organism"),
                  choices = list("GAS",
                                 "GBS",
                                 "PNEUMO",
                                 "GONO",
                                 "Curator"
                                 ),
                  selected = "GAS"),




      conditionalPanel(
        condition = "input.Org == 'GONO'",
        radioButtons("test", h3("Choose an analysis:"),
                     choices = list("AMR profile" = "AMR",
                                    "NG-STAR Type" = "NGSTAR",
                                    "23S rRNA Alleles" = "rRNA23S",
                                    "LabWare AMR profile" = "AMR_LW",
                                    "LabWare AMR profile ALL" = "AMR_ALL",
                                    "MLST Type" = "MLST",
                                    "NG-MAST Type" = "NGMAST",
                                    "MasterBlastR" = "MASTER",
                                    "LabWare Metrics" = "LW_METRICS"

                                    #"ARG-ANNOT/Resfinder/CARD" = "AMR_DB",
                                    #"VFDB (Virulence Factor Database" = "VFDB",
                                    #"Rename Contigs" = "CONTIGS",
                                    #"Rename 23S rRNA VCFs" = "VCF"


                                    ),
                     selected = "AMR")
      ),

      conditionalPanel(
        condition = "input.Org == 'GAS'",
        radioButtons("test2", h3("Choose an analysis:"),
                    choices = list("AMR profile" = "AMR",
                                   "Toxin profile" = "TOXINS",
                                   "MLST Type" = "MLST",
                                   "Virulence Factors" = "VIRULENCE",
                                   "emm Typing" = "EMM",
                                   "MasterBlastR" = "MASTER",
                                   "ARG-ANNOT/Resfinder/CARD" = "AMR_DB",
                                   "VFDB (Virulence Factor Database" = "VFDB",
                                   "LabWare Metrics" = "LW_METRICS"

                                   #"Rename Contigs" = "CONTIGS"
                    ),
                    selected = "TOXINS")
      ),


      conditionalPanel(
        condition = "input.Org == 'GBS'",
        radioButtons("test4", h3("Choose an analysis:"),
                     choices = list(
                                    "AMR profile" = "AMR",
                                    "MLST Type" = "MLST",
                                    "Serotyping CPS" = "SERO",
                                    "Serotyping SNPs" = "SERO2",
                                    "MasterBlastR" = "MASTER",
                                    "ARG-ANNOT/Resfinder/CARD" = "AMR_DB",
                                    "VFDB (Virulence Factor Database" = "VFDB",
                                    "LabWare Metrics" = "LW_METRICS"

                                    #"Rename Contigs" = "CONTIGS"
                     ),
                     selected = "SERO")
      ),



      conditionalPanel(
        condition = "input.Org == 'PNEUMO'",
        radioButtons("test3", h3("Choose an analysis:"),
                     choices = list("AMR + 23S rRNA profile" = "AMR_ALL",
                                    "AMR factors only" = "AMR",
                                    "23S rRNA Alleles" = "rRNA23S",
                                    "LabWare AMR profile" = "AMR_LW",
                                    "Serotype" = "SERO",
                                    "MLST Type" = "MLST",
                                    "Virulence Factors" = "VIRULENCE",
                                    "MasterBlastR" = "MASTER",
                                    "ARG-ANNOT/Resfinder/CARD" = "AMR_DB",
                                    "VFDB (Virulence Factor Database" = "VFDB",
                                    "LabWare Metrics" = "LW_METRICS"

                                    #"Rename Contigs" = "CONTIGS",
                                    #"Rename 23S rRNA VCFs" = "VCF"
                     ),
                     selected = "SERO")
      ),


      conditionalPanel(
        condition = "input.Org == 'Curator'",
        radioButtons("test5", h3("Choose an analysis:"),
                     choices = list(
                                    "Tree Sort Order" = "TREE_ORDER",
                                    "Tree Rename Leaves" = "TREE_RENAME",
                                    "Remove duplicate sequences" = "REMOVE_DUPLICATES"
                     ),
                     selected = "TREE_ORDER")
      ),

      conditionalPanel(
        condition = "input.test5 == 'REMOVE_DUPLICATES' && input.Org == 'Curator'",
        textInput("allele", h3("Enter starting allele no."), value = "1")
      ),

      conditionalPanel(
        condition = "input.Org != 'Curator'",
        textInput("locus", h3("Enter a locus to query or \"list\" for default loci list"),
                  value = "list")
      ),

#      textInput("locus", h3("Enter a locus to query or \"list\" for default loci list"),
#                value = "list"),

      conditionalPanel(
        condition = "input.Org != 'Curator'",
        textInput("sample", h3("Enter sample number or \"list\" for mulitple samples"),
                  value = "list")
      ),

#      textInput("sample", h3("Enter sample number or \"list\" for mulitple samples"),
#                value = "list"),


      actionButton("action", "Go"),
      actionButton("openxls", "Output"),
      conditionalPanel(
        condition = "input.Org != 'Curator'",
        actionButton("index", "MakeBlastdb")
      )

      #actionButton("index", "MakeBlastdb")

      #br(),
      #br(),
      #submitButton("Submit")


      ),



    mainPanel(
      #img(src = "Logo.png"),
      titlePanel(strong("Strep/STI WGS Analysis and Detection of Molecular Markers (WADE)")),
      textOutput("selected_Org"),
      textOutput("selected_test"),
      textOutput("selected_test2"),
      textOutput('selected_test3'),
      textOutput("selected_test4"),
      textOutput("entered_locus"),
      textOutput("entered_sample"),
      textOutput("button_value"),
      br(),
      textOutput("error_text"),
      textOutput("Status"),

      DT::dataTableOutput("profile_table")
      )

    )

)

# Define server logic -------------------------------------------------------------------------------
server <- function(input, output) {

  #if (input$sample == "") {input$sample <- "list"}
  #if (input$locus == "") {input$locus <- "list"}

  # output$selected_Org <- renderText({
  #   paste("The Organism you have selected is: ", input$Org)
  # })
  #
  # output$selected_test <- renderText({
  #   paste("The GONO test you have selected is: ", input$test)
  # })
  #
  # output$selected_test2 <- renderText({
  #   paste("The GAS test you have selected is: ", input$test2)
  # })
  #
  # output$selected_test3 <- renderText({
  #   paste("The PNEUMO test you have selected is: ", input$test3)
  # })
  #
  # output$entered_locus <- renderText({
  #   paste("The locus you have entered is: ", input$locus)
  # })
  #
  # output$entered_sample <- renderText({
  #   paste("The sample you have entered is: ", input$sample)
  # })



  observeEvent(input$action,
  {
    #if ((input$Org == "GONO") & (input$test == "NGSTAR") & (input$locus != "list"))
    #{
    #  input$test <- "NGSTAR_locus"
    #}

    if (input$Org == "Curator")
    {
      switch(input$test5,
             TREE_ORDER={output.df <- nwk_sort_order(input$Org, curr_work_dir)},
             TREE_RENAME={output.df <- nwk_rename_samples(input$Org, curr_work_dir)},
             REMOVE_DUPLICATES = {output.df <- remove_duplicate_fasta(input$Org, input$allele, curr_work_dir)},

             {output.df <- MASTER_pipeline(input$Org, input$test, input$sample, input$locus, curr_work_dir)}
      )
    }


    if (input$Org == "GONO")
      {


    switch(input$test,
           MLST={output.df <- MLST_pipeline(input$Org, input$sample, input$locus, curr_work_dir)},
           NGSTAR={output.df <- NGSTAR_MLST_pipeline(input$Org, input$sample, input$locus, curr_work_dir)},
           NGMAST={output.df <- NGMAST_pipeline(input$Org, input$sample, input$locus, curr_work_dir)},
           rRNA23S={output.df <- rRNA23S_pipeline(input$Org, input$sample, curr_work_dir)},
           CONTIGS={output.df <- rename_contigs(input$Org)},
           VCF={output.df <- rename_vcfs(input$Org)},
           AMR_DB=(output.df <- AMR_DATABASE_pipeline(input$Org, input$sample, curr_work_dir)),
           VFDB=(output.df <- VFDB_pipeline(input$Org, input$sample, curr_work_dir)),
           AMR_LW={output.df <- labware_gono_amr(input$Org, curr_work_dir)},
           LW_METRICS={output.df <- metrics(input$Org, curr_work_dir)},


           AMR_ALL={output.df <- MASTER_pipeline(input$Org, input$test, input$sample, "list", curr_work_dir)
                    output.df <- NGSTAR_MLST_pipeline(input$Org, input$sample, input$locus, curr_work_dir)
                   output.df <- rRNA23S_pipeline(input$Org, input$sample, curr_work_dir)
                    output.df <- labware_gono_amr(input$Org, curr_work_dir)
           },


           {output.df <- MASTER_pipeline(input$Org, input$test, input$sample, input$locus, curr_work_dir)}
    )
    }

    if (input$Org == "GAS")
    {


      switch(input$test2,
             MLST={output.df <- MLST_pipeline(input$Org, input$sample, input$locus, curr_work_dir)},
             AMR_DB=(output.df <- AMR_DATABASE_pipeline(input$Org, input$sample, curr_work_dir)),
             EMM=(output.df <- EMM_pipeline(input$Org, input$sample, input$locus, curr_work_dir)),
             VFDB=(output.df <- VFDB_pipeline(input$Org, input$sample, curr_work_dir)),
             CONTIGS={output.df <- rename_contigs(input$Org)},
             LW_METRICS={output.df <- metrics(input$Org, curr_work_dir)},
             AMR={
                 output.df <- MASTER_pipeline(input$Org, input$test2, input$sample, input$locus, curr_work_dir)
                 output.df <- labware_gas_amr(input$Org, curr_work_dir)
                 },
             TOXINS={
                    output.df <- MASTER_pipeline(input$Org, input$test2, input$sample, input$locus, curr_work_dir)
                    output.df <- labware_gas_toxins(input$Org, curr_work_dir)
                    },

             {output.df <- MASTER_pipeline(input$Org, input$test2, input$sample, input$locus, curr_work_dir)}
      )
    }


    if (input$Org == "GBS")
    {


      switch(input$test4,
             MLST={output.df <- MLST_pipeline(input$Org, input$sample, input$locus, curr_work_dir)},
             AMR={
               output.df <- MASTER_pipeline(input$Org, input$test4, input$sample, input$locus, curr_work_dir)
               output.df <- labware_gbs_amr(input$Org, curr_work_dir)
             },
             AMR_DB=(output.df <- AMR_DATABASE_pipeline(input$Org, input$sample, curr_work_dir)),
             VFDB=(output.df <- VFDB_pipeline(input$Org, input$sample, curr_work_dir)),
             CONTIGS={output.df <- rename_contigs(input$Org)},
             LW_METRICS={output.df <- metrics(input$Org, curr_work_dir)},
             SERO=(output.df <- SEROTYPE_pipeline(input$Org, input$sample, input$locus, input$test4, curr_work_dir)),
             SERO2=(output.df <- SEROTYPE_pipeline(input$Org, input$sample, input$locus, input$test4, curr_work_dir)),
             {output.df <- MASTER_pipeline(input$Org, input$test4, input$sample, input$locus, curr_work_dir)}
      )
    }




    if (input$Org == "PNEUMO")
    {
      switch(input$test3,
             AMR_ALL={
               output.df <- MASTER_pipeline(input$Org, input$test3, input$sample, "list", curr_work_dir)
               output.df <- rRNA23S_pipeline(input$Org, input$sample, curr_work_dir)
               output.df <- labware_pneumo_amr(input$Org, curr_work_dir)
             },
             AMR={output.df <- MASTER_pipeline(input$Org, input$test3, input$sample, input$locus, curr_work_dir)},
             MLST={output.df <- MLST_pipeline(input$Org, input$sample, input$locus, curr_work_dir)},
             rRNA23S={output.df <- rRNA23S_pipeline(input$Org, input$sample, curr_work_dir)},
             AMR_LW={output.df <- labware_pneumo_amr(input$Org, curr_work_dir)},
             AMR_DB=(output.df <- AMR_DATABASE_pipeline(input$Org, input$sample, curr_work_dir)),
             VFDB=(output.df <- VFDB_pipeline(input$Org, input$sample, curr_work_dir)),
             CONTIGS={output.df <- rename_contigs(input$Org)},
             VCF={output.df <- rename_vcfs(input$Org)},
             #SERO=(output.df <- PneumoCaT_pipeline(input$Org, input$sample)),
             SERO=(output.df <- SEROTYPE_pipeline(input$Org, input$sample, input$locus, input$test3, curr_work_dir)),
             LW_METRICS={output.df <- metrics(input$Org, curr_work_dir)},
             VIRULENCE={
                        output.df <- MASTER_pipeline(input$Org, input$test3, input$sample, input$locus, curr_work_dir)
                        output.df <- labware_pneumo_virulence(input$Org, curr_work_dir)
                        },

             {output.df <- MASTER_pipeline(input$Org, input$test3, input$sample, input$locus, curr_work_dir)}
      )
    }

     dir_file <- paste(curr_work_dir, "DirectoryLocations.csv", sep="")
     #cat("\n directory string passed to function: ", curr_dir, "\n")
     Directories.df <- read.csv(dir_file, header = TRUE, sep = ",", stringsAsFactors = FALSE)
     Directories_org.df <- filter(Directories.df, OrgID == input$Org)
     local_dir <- Directories_org.df$LocalDir


     unlink(paste(local_dir, "Output\\output_profile.csv", sep = ""))
     write.csv(output.df, paste(local_dir, "Output\\output_profile.csv", sep = ""), row.names = F)
     #"C:\WGS_Typing\Output\output_profile.csv"



    output$profile_table <- renderDataTable({output.df})




  })

  observeEvent(input$openxls,
               {
                 dir_file <- paste(curr_work_dir, "DirectoryLocations.csv", sep="")
                 #cat("\n directory string passed to function: ", curr_dir, "\n")
                 Directories.df <- read.csv(dir_file, header = TRUE, sep = ",", stringsAsFactors = FALSE)
                 Directories_org.df <- filter(Directories.df, OrgID == input$Org)
                 local_dir <- Directories_org.df$LocalDir

                 shell.exec(paste(local_dir, "Output\\output_profile.csv", sep = ""))
               })

  observeEvent(input$index,
               {
                 switch(input$Org,
                        GONO = {runblast <- Index_pipeline(input$Org, input$test, input$locus, curr_work_dir)},
                        GAS = {runblast <- Index_pipeline(input$Org, input$test2, input$locus, curr_work_dir)},
                        GBS = {runblast <- Index_pipeline(input$Org, input$test4, input$locus, curr_work_dir)},
                        PNEUMO = {runblast <- Index_pipeline(input$Org, input$test3, input$locus, curr_work_dir)}

                 )


               })


}

# Run the app -------------------------------------------------------------------------------------
shinyApp(ui = ui, server = server)

#--------------------------------------------------------------------------------------------------
# Code to get the console output displayed -- it works but only displayed after results are complete
# not displayed as they scroll through the console.
#
# withConsoleRedirect <- function(containerId, expr) {
#   # Change type="output" to type="message" to catch stderr
#   # (messages, warnings, and errors) instead of stdout.
#   txt <- capture.output(results <- expr, type = "output")
#   if (length(txt) > 0) {
#     insertUI(paste0("#", containerId), where = "beforeEnd",
#              ui = paste0(txt, "\n", collapse = "")
#     )
#   }
#   results
# }
#
# # Example usage
#
# ui <- fluidPage(
#   pre(id = "console")
# )
#
# server <- function(input, output, session) {
#   observe({
#     invalidateLater(1000)
#
#     withConsoleRedirect("console", {
#       str(cars)
#     })
#   })
# }
#--------------------------------------------------------------------------------------------------
# how to set up tabs with a different sidebar in each tab
# -- would be good for the curator functions or the phlogenetic tree annotation scripts.
# library(shiny)
# library(plotly)
#
# shinyApp(
#   ui = fluidPage(
#     tabsetPanel(
#       tabPanel("Map", fluid = TRUE,
#                sidebarLayout(
#                  sidebarPanel(selectInput("Country", "Select Country", choices = "", selected = "")),
#                  mainPanel(
#                    htmlOutput("Attacks")
#                  )
#                )
#       ),
#       tabPanel("plot", fluid = TRUE,
#                sidebarLayout(
#                  sidebarPanel(sliderInput("year", "Year:", min = 1968, max = 2009, value = 2009, sep='')),
#                  mainPanel(fluidRow(
#                    column(7,  plotlyOutput("")),
#                    column(5, plotlyOutput(""))
#                  )
#                  )
#                )
#       )
#     )
#   ),
#   server = function(input, output) {
#
#   }
# )
