library(shiny)
library(Seurat)
library(biomaRt)
library(data.table)
library(dplyr)
library(bslib)
library(DT)
library(ggplot2)
library(viridis)

options(shiny.maxRequestSize = 80000 * 1024^2) #to increase the upload size of seurat object

server <- function(input, output) {
  convertedData <- reactiveVal(NULL)

  obj <- reactive({
    req(input$file)
    readRDS(input$file$datapath)
  })

  # rks
  sl <- reactive({
    if (class(obj()) == "Seurat"){
      return(NULL)
    }
  })
  output$seuratLoaded <- reactive({
    return(is.null(sl()))
  })
  outputOptions(output, 'seuratLoaded', suspendWhenHidden = FALSE)

  # load object as reactive for input into observe event

  observeEvent(input$convertButton, {
    obj <- reactiveVal()
    req(input$file)
    obj <- readRDS(input$file$datapath) # obj() #
    assay <- input$Selected_assay
    print(assay)
    get_counts_matrix <- function(obj, assay) {
      version <- as.character(Version(object = obj))
      if (startsWith(version, "5")) {
        # Check if layers$counts exists
        if (!is.null(obj[[assay]]$counts)) {
          return(obj[[assay]]$counts)
        } else {
          stop("The counts matrix could not be found in the layers slot.")
        }
      } else if (startsWith(version, "3") || startsWith(version, "4")) {
        # Check if @counts exists
        if (!is.null(obj[[assay]]@counts)) {
          return(obj[[assay]]@counts)
        } else {
          stop("The counts matrix could not be found in the counts slot.")
        }
      } else {
        stop("Unsupported version.")
      }
    }

    counts_matrix <- get_counts_matrix(obj, assay)
    genes <- rownames(counts_matrix)
    gene_all <- sub("^hg38-|^mm10-", "", genes)
    if (input$species != "Custom") {
      # Use attributes and filters based on selected species
      species_lookup <- data.frame(
        Mouse = c(ensembl_id = "mmusculus_gene_ensembl", attributes = 'mgi_symbol', filters = 'mgi_symbol'),
        Human = c(ensembl_id = "hsapiens_gene_ensembl", attributes = 'ensembl_gene_id', filters = 'ensembl_gene_id'),
        Zebrafish = c(ensembl_id = "drerio_gene_ensembl", attributes = 'zfin_id_symbol', filters = 'zfin_id_symbol'),
        Rat = c(ensembl_id = "rnorvegicus_gene_ensembl", attributes = 'rgd_symbol', filters = 'rgd_symbol'),
        stringsAsFactors = FALSE
      )
      species_info <- species_lookup[[input$species]]
      print(species_info)
    } else {
      # Use custom attributes and filters
      ensembl_id <- req(input$customEnsemblId)
      attributes <- req(input$customAttributes)
      filters <- input$customFilters
      # If filters is empty, set it to the same value as attributes
      filters <- if (nzchar(filters)) filters else attributes
      # Create a dataframe for custom species info
      species_info <- data.frame(
        ensembl_id = ensembl_id,
        attributes = attributes,
        filters = filters,
        stringsAsFactors = FALSE
      )
    }
    #We get the species gene list from biomart server using useEnsembl function and converted variable helps us to create a dataframe of genes names equivalent of species used and human gene
    ######## tested on different biomart servers ########
    #mart.species <- useEnsembl("ensembl", species_info[[1]], mirror = 'useast', host = "https://dec2021.archive.ensembl.org")
    # mart.human <- useEnsembl("ensembl", "hsapiens_gene_ensembl", mirror = 'useast', host = "https://dec2021.archive.ensembl.org")
    mart.species <- useEnsembl("ensembl", species_info[[1]], mirror = 'useast',host = "https://nov2020.archive.ensembl.org")
    mart.human <- useEnsembl("ensembl", "hsapiens_gene_ensembl", mirror = 'useast', host = "https://nov2020.archive.ensembl.org")
    converted <- biomaRt::getLDS(
      attributes =  c(species_info[[2]],"gene_biotype"),
      filters = species_info[[3]],
      values = as.character(gene_all),
      mart = mart.species,
      attributesL = c('hgnc_symbol'),
      martL = mart.human,
      uniqueRows = T
    )
    print(class(converted))
    print(str(converted))
    species_symbol <- function(attr) {
      parts <- strsplit(attr, "_")[[1]]
      formatted <- paste0(toupper(parts[1]), ".symbol")
      return(formatted)
    }
    #species_sym gives us the gene symbol of species using the function species_symbol which is strip split function#####################################################################
    species_sym <- species_symbol(species_info[[2]])
    #Selecting the PDOX model on the app which has two species information (human and mouse/rat/zebrafish) in the seurat object,but we just need to convert the species data into human##########
    if (input$Select_model == "Yes") {
      print("PDOX model to convert species to human gene set successful......")
      #Necessary to paste the gene symbols here, as the current PDOX model objects have these symbols attached to them to recognize the MM10/HG38 IDENTIFIER##################################################################
      converted$MGI.symbol <- paste0("mm10-",converted$MGI.symbol)
      converted$HGNC.symbol <- paste0("hg38-", converted$HGNC.symbol)
      hasspecies <- which(rownames(counts_matrix) %in% converted[[species_sym]])
      tmp.counts <- get_counts_matrix(obj,assay)[hasspecies,]
    }
    else {
      print("Running in 'normal' mode, if input data is dual-species, please select to run in PDX mode!")
      #tmp.counts <- obj[[assay]]@counts
      length(rownames(obj))
      length(converted[[species_sym]])
      genes_present_converted <- which(rownames(tryCatch(obj[[assay]]$counts, error = function(e) NULL) %||% obj[[assay]]@counts) %in% converted[[species_sym]])
      tmp.counts <- get_counts_matrix(obj,assay)[genes_present_converted,]
    }
    species_genes <- getBM(
      attributes = c(species_info[[2]],"ensembl_gene_id", "gene_biotype"),
      uniqueRows = TRUE,
      mart = mart.species
    )
    symbol_id <- species_info[[2]]
    species_converted_hg <-  biomaRt::getLDS(
      attributes =  c(species_info[[2]],"gene_biotype","ensembl_gene_id"),
      filters = species_info[[3]],
      values = species_genes[[symbol_id]],
      mart = mart.species,
      attributesL = c('hgnc_symbol'),
      martL = mart.human,
      uniqueRows = T
    )
    converted_unique_h <- species_converted_hg[!duplicated(species_converted_hg$HGNC.symbol), ]
    gene_classification <- as.data.frame(table(converted$Gene.type))
    colnames(gene_classification) <- c("Gene_Type", "Freq")
    gene_classification_DB <- as.data.frame(table(converted_unique_h$Gene.type))
    colnames(gene_classification_DB) <- c("Gene_Type", "Freq")

    output$gene_type <- renderPlot({
      graph1 <- ggplot(gene_classification_DB, aes(x = "", y = Freq, fill = Gene_Type)) +
        geom_bar(width = 1, stat = "identity") +
        coord_polar(theta = "y") +
        scale_fill_viridis(discrete = TRUE, option = "turbo") +
        theme_void() + labs(title = "Species DB Distribution ") +  theme(
          plot.title = element_text(hjust = 0.5, size = 15, face = "bold", margin = margin(b = 10)),
          legend.title = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 12))

      graph2 <- ggplot(gene_classification, aes(x = "", y = Freq, fill = Gene_Type)) +
        geom_bar(width = 1, stat = "identity") +
        coord_polar(theta = "y") +
        scale_fill_viridis(discrete = TRUE, option = "turbo") +
        theme_void() +
        labs(title = "Dataset converted Distribution")  +  theme(
          plot.title = element_text(hjust = 0.5, size = 15, face = "bold", margin = margin(b = 10)),
          legend.title = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 12))
      cowplot::plot_grid(plotlist = list(graph1,graph2),ncol = 2)
    })
    species_hg_class <- as.data.frame(table(converted_unique_h$Gene.type))
    colnames(species_hg_class) <- c("Gene_Type", "Freq")
    converted_unique <- converted[!duplicated(converted$HGNC.symbol), ]
    ortho_class_Data <- as.data.frame(table(converted_unique$Gene.type))
    colnames(ortho_class_Data) <- c("Gene_Type", "Freq")
    dataset_pco <- ortho_class_Data[ortho_class_Data$Gene_Type == "protein_coding", "Freq"]
    biomart_ortho_pco_mouse <- species_hg_class[species_hg_class$Gene_Type == "protein_coding", "Freq"] #17,620
    matched <- dataset_pco/biomart_ortho_pco_mouse * 100
    unmatched <- 100 - matched
    pie_data <- data.frame(
      category = c("Matched", "Unmatched"),
      count = c(matched,unmatched)
    )
    # # Calculate intersected genes
    matched_genes <- intersect(rownames(obj), converted[[species_sym]])
    unmatched_genes <- setdiff(rownames(obj), matched_genes)
    all_genes <- c(matched_genes, unmatched_genes)
    type <- c(rep('matched', length(matched_genes)), rep('unmatched', length(unmatched_genes)))
    data_df <- data.frame(gene = all_genes, type = type)
    rv <- reactiveValues(data = data_df)
    output$pieChart <- renderPlot({
      ggplot(pie_data, aes(x = "", y = count, fill = category)) +
        geom_bar(stat = "identity", width = 0.9) +
        coord_polar(theta = "y") +
        theme_void() +
        geom_text(aes(label = paste0(round(count, 2), "%")),
                  position = position_stack(vjust = 0.5),
                  color = "white", size = 5, fontface = "bold") +
        labs(title = "Protein Coding genes overlap \n between Biomart DB and data uploaded") +
        theme(
          plot.title = element_text(hjust = 0.5, size = 15, face = "bold", margin = margin(b = 10)),
          legend.title = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 12),
          legend.position = "right",
          legend.box.background = element_rect(color = "grey", size = 0.5),
          legend.box.margin = margin(6, 6, 6, 6)
        )
    })
    output$geneTable <- renderDT(
      data_df ,
      options = list(
        paging = TRUE,
        pageLength = 10,
        autoWidth = TRUE,
        server = TRUE,           # Use server-side processing
        dom = 'Bfrtip'
      ),
      selection = 'single',
      filter = 'bottom',
      rownames = FALSE
    )
    output$download_geneslist <- downloadHandler(
      filename = function() {
        paste("geneslist", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        write.csv(data.frame(data_df), file, row.names = FALSE)
      }
    )

    # condition for genes list download button...
    output$genes_list_ready <- reactive({
      return(!is.null(data_df))
    })
    outputOptions(output, 'genes_list_ready', suspendWhenHidden = FALSE)

    converted[[species_sym]] <- as.character(converted[[species_sym]])
    converted$HGNC.symbol <- as.character(converted$HGNC.symbol)
    #mapping required to convert species rownames from the object to human and use these tmp.counts to create a new seurat object which only has human gene information
    rownames(tmp.counts) <- make.unique(plyr::mapvalues(
      x = as.character(rownames(tmp.counts)),
      from = as.character(converted[[species_sym]]),
      to = converted$HGNC.symbol,
      warn_missing = FALSE
    ))
    updated_assay_name <- paste0(assay, "_ortho")
    tmp.counts <- as(tmp.counts, "dgCMatrix")
    tmp <- CreateSeuratObject(counts = tmp.counts,assay = updated_assay_name)
    tmp@meta.data <- obj@meta.data
    tmp@reductions <- obj@reductions
    tmp@assays[[assay]] <- obj@assays[[assay]]
    # Add images if assay is "Spatial"
    if (assay == "Spatial") {
      tmp@images <- obj@images
    }
    #tmp <- NormalizeData(tmp)
    #tmp <- ScaleData(tmp)
    convertedData(tmp)
  })
  output$status <- renderUI({
    if (!is.null(convertedData())) {
      tags$span("Conversion completed. You can now download the OrthologAL converted Seurat object.", style = "color: green;")
    } else {
      tags$span("Upload an RDS file. Once loaded, click 'Convert' to start.", style = "color: blue;")
    }
  })
  output$download_visibile_in_main_page <- renderUI({
    if (!is.null(convertedData())) {
      downloadButton("downloadButton", "Download Converted Data", class = "btn btn-success btn-block mt-3")
    }
  })
  output$downloadButton <- downloadHandler(
    filename = function() {
      if (!is.null(convertedData()))
        paste("OrthologAL_", input$file$name)
    },
    content = function(file) {
      if (!is.null(convertedData()))
        saveRDS(convertedData(), file)
    }
  )
}
