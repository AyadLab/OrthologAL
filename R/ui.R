library(shiny)
library(Seurat)
library(biomaRt)
library(data.table)
library(dplyr)
library(bslib)
library(DT)
library(ggplot2)
library(viridis)

ui <- fluidPage(

  tags$script(src = "https://kit.fontawesome.com/070e476711.js"),
  titlePanel("OrthologAL", windowTitle = "OrthologAL"),
  hr(),
  # br(),
  tags$head(tags$style(HTML(
    ".nav.nav-pills.nav-stacked > .active > a, .nav.nav-pills.nav-stacked > .active > a:hover {
    background-color: #000000;
  }

  .well {
      min-height: 20px;
      max-width: 5000px;
      padding: 19px;
      margin-bottom: 20px;
      background-color: #ffffff;
      border: 1px solid #ffffff;
      border-radius: 4px;
      -webkit-box-shadow: inset 0 1px 1px rgba(0,0,0,.05);
      box-shadow: inset 0 1px 1px rgba(0,0,0,.05);
      font-family: 'sans-serif', Arial Rounded MT Bold;

  }

                            "))),
  navlistPanel(well = T, fluid = T, widths = c(2, 8),
               # tags$head(tags$style(HTML(".tab-content { height: 83vh; overflow-y: auto !important; }" ))),
               tabPanel(tags$div(
                 tags$i(class = "fa-sharp fa-solid fa-desktop"),
                 tags$span("- Overview"),
                 tags$style(type = "text/css", "li a{color:#000000;}")
               ),
               h3("About OrthologAL"),
               br(),
               p("OrthologAL is a no-code solution for converting non-human gene expression data to that of corresponding human orthologs in Seurat objects"),
               br(),
               p("The OrthologAL package leverages the data mining tool biomaRt to access different gene sets from ENSEMBL, and facilitates the interaction with these servers."),
               br(),
               p("Researchers can effortlessly input their Seurat object, a standard datatype format for single-cell RNA sequencing (scRNAseq), single-nuclei RNA sequencing (snRNAseq), or spatial transcriptomics (stRNAseq) data of any species, and
           OrthologAL will output a human-gene converted Seurat object for download."),
               br(),
               p("Uniquely, OrthologAL can process dual-species model data, such as PDXs which may contain host cells/tissue as well as human tumor."),
               br(),
               p("In the OrthologAL manuscript (Chowdary et al.), we benchmark OrthologAL conversion data effiency in single-cell, single-nuclei, and spatial transcriptomics datasets from cancer and spinal cord injury models."),
               br(),
               tags$img(src = "Figure_1.png", width = 1000, height = 900),
               br(),
               br(),
               div(
               a(href="https://github.com/AyadLab", "Keep up-to-date with the most recent release of OrthologAL at The Ayad lab GitHub", style = "font-size: 18px; font-weight: bold;"),
               style = "text-align: center;"),
               hr(),
               ),

               tabPanel(tags$div(
                 tags$i(class="fa-sharp fa-solid fa-magnifying-glass-chart"),
                 tags$span("- Run OrthologAL")
               ), #put contents for actual application here
               ######################################################################################################################

               conditionalPanel(condition = "!input.acceptAgreement",
                                h1("OrthologAL User Agreement"),
                                hr(),
                                br(),
                                h1("License Agreement"),
                                # hr(),
                                br(), # Note this is a draft... I dont know if this is right. Not a lawyer!
                                p(strong("1. The Board of Trustees of the Georgetown University (“Georgetown”) provides OrthologAL software and code (“Service”) free of charge for non-commercial use only. Use of the Service by any commercial entity for any purpose, including research, is prohibited.")),
                                p(strong("2. By using the Service, you agree to be bound by the terms of this Agreement. Please read it carefully.")),
                                p(strong("3. You agree not to use the Service for commercial advantage, or in the course of for-profit activities. You agree not to use the Service on behalf of any organization that is not a non-profit organization. Commercial entities wishing to use this Service should contact Georgetown University’s Office of Technology Licensing.")),
                                p(strong("4. THE SERVICE IS OFFERED “AS IS”, AND, TO THE EXTENT PERMITTED BY LAW, GEORGETOWN MAKES NO REPRESENTATIONS AND EXTENDS NO WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED. GEORGETOWN SHALL NOT BE LIABLE FOR ANY CLAIMS OR DAMAGES WITH RESPECT TO ANY LOSS OR OTHER CLAIM BY YOU OR ANY THIRD PARTY ON ACCOUNT OF, OR ARISING FROM THE USE OF THE SERVICE. YOU HEREBY AGREE TO DEFEND AND INDEMNIFY GEORGETOWN, ITS TRUSTEES, EMPLOYEES, OFFICERS, STUDENTS, AGENTS, FACULTY, REPRESENTATIVES, AND VOLUNTEERS (“GEORGETOWN INDEMNITEES”) FROM ANY LOSS OR CLAIM ASSERTED AGAINST GEORGETOWN INDEMNITEES ARISING FROM YOUR USE OF THE SERVICE.")),
                                p(strong("5. All rights not expressly granted to you in this Agreement are reserved and retained by GEORGETOWN or its licensors or content providers. This Agreement provides no license under any patent.")),
                                p(strong("6. You agree that this Agreement and any dispute arising under it is governed by the laws of the District of Columbia, United States of America, applicable to agreements negotiated, executed, and performed within the DISTRICT OF COLUMBIA")),
                                p(strong("7. Subject to your compliance with the terms and conditions set forth in this Agreement, GEORGETOWN grants you a revocable, non-exclusive, non-transferable right to access and make use of the Service.")),
                                # hr(),
                                br(),
                                p(em("Do you accept the terms and conditions in this agreement?")),
                                fluidRow(column(width = 4, actionButton(inputId = "acceptAgreement", label = "Accept")), column(width = 4, actionButton(inputId = "rejectAgreement", label = "Reject and close"))),
                                conditionalPanel(condition = "input.acceptAgreement",
                                                 "Please scroll down for instructions."),
                                # hr(),
                                br()),
                                # agreement accepted, run app
                                conditionalPanel(condition = "input.acceptAgreement",
                                                 titlePanel("OrthologAL", windowTitle = "OrthologAL"),
                                                 sidebarLayout(
                                                   sidebarPanel(
                                                     fileInput("file", "Choose RDS File"),
                                                     selectInput("species", "Select Species", choices = c("Mouse", "Zebrafish","Rat","Human","Custom"), selected = "Mouse"),
                                                     conditionalPanel(
                                                       condition = "input.species == 'Custom'",
                                                       tags$div(textInput("customEnsemblId", "Enter Ensembl ID", placeholder = "e.g., mmusculus_gene_ensembl"), class = "text-input"),
                                                       tags$div(textInput("customAttributes", "Enter Attributes", placeholder = "e.g., mgi_symbol"), class = "text-input"),
                                                       tags$div(textInput("customFilters", "Enter Filters (Optional)", placeholder = "e.g., mgi_symbol"), class = "text-input")
                                                     ),
                                                     selectInput("Selected_assay", "Select Assay", choices = c("RNA", "SCT", "Spatial", "Integrated","alra"), selected = "RNA"),
                                                     #selectInput("Select_model", "Select Data type", choices = c("scRNAseq","snRNAseq","spatial transcriptomics","patient derived xenograft or PDX"),selected = "No Selection required"),
                                                     selectInput("Select_model", "Run in PDX mode?", choices = c("Patient Derived Xenograft (PDX)","No Selection required"),selected = "No Selection required"),
                                                     #using bootstrap to make it more app like %structure%
                                                     div(class = "form-group",
                                                         conditionalPanel(condition = 'output.seuratLoaded',
                                                                          actionButton("convertButton", "Convert", class = "btn btn-primary btn-block")
                                                         )
                                                     ),
                                                     uiOutput("download_visibile_in_main_page"),
                                                   ),
                                                   mainPanel(
                                                     code("status"),
                                                     uiOutput("status"),
                                                     br(),
                                                     navset_pill(
                                                       #title = "OUTPUT",
                                                       # Panel with plot ----
                                                       nav_panel("Genetype Plot", plotOutput("gene_type")),
                                                       nav_panel("% Match Plot", plotOutput("pieChart")),
                                                       #nav_panel("Genes Detected",plotOutput("genesdetected")),
                                                       # Panel with table ----
                                                       nav_panel("Table", DTOutput("geneTable"))
                                                     ),
                                                     conditionalPanel(condition = 'output.genes_list_ready',
                                                                      downloadButton("download_geneslist", "Download Genes Mapping list")
                                                     )


                                                   )
                                                 )

                                ) # end condition acccepted agreement
               ), # end app tab

               tabPanel(tags$div(
                 tags$i(class = "fa-brands fa-github"),
                 tags$span("- Github")
               ), # put contents for R package installation here
               h1("The OrthologAL R Package"),
               br(),
               div(
                 a(href="https://github.com/AyadLab", "Keep up-to-date with the most recent release of OrthologAL at The Ayad lab GitHub", style = "font-size: 18px; font-weight: bold;"),
                 style = "text-align: center;"),

               ),

               tabPanel(tags$div(
                 tags$i(class="fa-sharp fa-solid fa-envelope"),
                 tags$span("- Contact")
               ), #put contents for Contact here
               splitLayout(
                 wellPanel(
                   # hr(),
                   br(),
                   div(
                     p(strong("Rishika Chowdary, MS")),
                     p(em("Bioinformatician")),
                     style = "text-align: center;"
                   ),
                   hr(),
                   div(
                     img(src = "rishika.png", width = 185, height = 230),
                     style = "text-align: center;"
                     ),
                   br(), br(),
                   hr()
                 ),
                 wellPanel(
                   # hr(),
                   br(),
                   div(
                     p(strong("Robert K. Suter, PhD")),
                     p(em("Assistant Professor, LCCC")),
                     style = "text-align: center;"
                   ),
                   hr(),
                   # img(src = "rks_headshot_2022.png", width = 185*1.25, height = 160*1.25),
                   div(
                     tags$img(
                       src = "rks_headshot_2022.png",
                       width = 170*1.4, height = 160*1.4,
                       alt = "Suter lab"
                     ), style = "text-align: center;"
                   ),
                   br(), br(),
                   div(
                     a(href="mailto:rks82@georgetown.edu", "Contact"),
                     br(),
                     hr(),
                     a(href="https://suterlab.com", "The Suter Lab"),
                     style = "text-align: center;"
                   ),
                   hr()
                 ),
                 wellPanel(
                   # hr(),
                   br(),
                   div(
                     p(strong("Nagi G. Ayad, PhD")),
                     p(em("Professor, LCCC")),
                     style = "text-align: center;"
                     ),
                   hr(),
                   div(
                     img(src = "nagi_headshot.png", width = 100*1.4, height = 160*1.4),
                     style = "text-align: center;"
                     ),
                   br(), br(),
                   div(
                     a(href="mailto:na853@georgetown.edu", "Contact"),
                     br(),
                     hr(),
                     a(href="https://sites.google.com/georgetown.edu/ayadlab/home", "The Ayad Lab"),
                     style = "text-align: center;"
                   ),
                   hr()
                 )
               )
               )
  )
)
