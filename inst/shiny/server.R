library(shiny)
library(Seurat)
library(biomaRt)
library(data.table)
library(dplyr)
library(bslib)
library(DT)
library(ggplot2)
library(viridis)
library(shinyFiles)

options(shiny.maxRequestSize = 80000 * 1024^2) # to increase the upload size of seurat object

server <- function(input, output, session) {
  convertedData <- reactiveVal(NULL)
  uploaded_filename <- reactiveVal("converted_data.rds")
  
  # --- shinyFiles Setup for UPLOADING ---
  volumes <- getVolumes()() # Ensure you call the function to get the actual volumes list
  shinyFileChoose(input, "file_server", roots = volumes, session = session, filetypes = c('', 'rds', 'RDS'))
  
  output$selected_server_file <- renderText({
    if (is.integer(input$file_server)) {
      "No file selected"
    } else {
      as.character(parseFilePaths(volumes, input$file_server)$datapath[1])
    }
  })
  
  # --- shinyFiles Setup for SAVING ---
  shinyFileSave(input, "save_server_btn", roots = volumes, session = session)
  # ------------------------
  
  # Dynamically load the Seurat object based on which upload method the user selected
  obj <- reactive({
    if (input$upload_source == "local") {
      req(input$file_local)
      uploaded_filename(input$file_local$name)
      readRDS(input$file_local$datapath)
    } else {
      req(!is.integer(input$file_server))
      file_info <- parseFilePaths(volumes, input$file_server)
      req(nrow(file_info) > 0)
      
      file_path <- as.character(file_info$datapath[1])
      req(file.exists(file_path)) 
      
      uploaded_filename(as.character(file_info$name[1]))
      
      tryCatch({
        readRDS(file_path)
      }, error = function(e) {
        showNotification(paste("Could not read RDS file:", e$message), type = "error", duration = 15)
        return(NULL)
      })
    }
  })
  
  sl <- reactive({
    if (inherits(obj(), "Seurat")){
      return(NULL)
    }
  })
  
  output$seuratLoaded <- reactive({
    return(is.null(sl()))
  })
  outputOptions(output, 'seuratLoaded', suspendWhenHidden = FALSE)
  
  observeEvent(input$convertButton, {
    
    # 1. Safely extract the Seurat object from the reactive function
    seurat_obj <- obj() 
    req(seurat_obj)
    
    assay <- input$Selected_assay
    print(assay)
    
    get_counts_matrix <- function(seurat_data, target_assay) {
      version <- as.character(Version(object = seurat_data))
      if (startsWith(version, "5")) {
        if (!is.null(seurat_data[[target_assay]]$counts)) {
          return(seurat_data[[target_assay]]$counts)
        } else {
          stop("The counts matrix could not be found in the layers slot.")
        }
      } else if (startsWith(version, "3") || startsWith(version, "4")) {
        if (!is.null(seurat_data[[target_assay]]@counts)) {
          return(seurat_data[[target_assay]]@counts)
        } else {
          stop("The counts matrix could not be found in the counts slot.")
        }
      } else {
        stop("Unsupported version.")
      }
    }
    
    counts_matrix <- get_counts_matrix(seurat_obj, assay)
    genes <- rownames(counts_matrix)
    gene_all <- sub("^hg38-|^mm10-", "", genes)
    
    if (input$species != "Custom") {
      species_lookup <- data.frame(
        Mouse = c(
          ensembl_id = "mmusculus_gene_ensembl",
          attributes = 'mgi_symbol',
          filters = 'mgi_symbol',
          filename = "ortho_df_Mouse_Human" 
        ),
        Zebrafish = c(
          ensembl_id = "drerio_gene_ensembl",
          attributes = 'zfin_id_symbol',
          filters = 'zfin_id_symbol',
          filename = "ortho_df_Zebrafish_Human"     
        ),
        Rat = c(
          ensembl_id = "rnorvegicus_gene_ensembl",
          attributes = 'rgd_symbol',
          filters = 'rgd_symbol',
          filename = "ortho_df_Rat_Human"            
        ),
        stringsAsFactors = FALSE
      )
      species_info <- species_lookup[[input$species]]
      print(species_info)
    } else {
      ensembl_id <- req(input$customEnsemblId)
      attributes <- req(input$customAttributes)
      filters <- input$customFilters
      filters <- if (nzchar(filters)) filters else attributes
      species_info <- data.frame(
        ensembl_id = ensembl_id,
        attributes = attributes,
        filters = filters,
        stringsAsFactors = FALSE
      )
    }
    
    target_object_name <- as.character(species_info[4])
    
    tryCatch({
      local_rda_path <- file.path("data", paste0(target_object_name, ".rda"))
      if (file.exists(local_rda_path)) {
        load(local_rda_path) 
      }
      
      master_ref <- get(target_object_name)
      print(paste0("Successfully loaded .rda object: ", target_object_name))
      
    }, error = function(e) {
      stop(paste0("Critical Error: Could not find dataset '", target_object_name, "'. Make sure it is saved as an .rda file in the data/ folder."))
    })
    
    species_symbol <- function(attr) {
      parts <- strsplit(attr, "_")[[1]]
      formatted <- paste0(toupper(parts[1]), ".symbol")
      return(formatted)
    }
    
    species_sym <- species_symbol(species_info[[2]])
    converted <- master_ref[master_ref[[species_sym]] %in% as.character(gene_all), ]
    converted <- converted[!duplicated(converted[[species_sym]]), ]
    
    if (input$Select_model == "Yes") {
      print("PDOX model to convert species to human gene set successful......")
      converted$MGI.symbol <- paste0("mm10-",converted$MGI.symbol)
      converted$HGNC.symbol <- paste0("hg38-", converted$HGNC.symbol)
      
      hasspecies <- which(rownames(counts_matrix) %in% converted[[species_sym]])
      tmp.counts <- counts_matrix[hasspecies, ]
      
    } else {
      print("Running in 'normal' mode, if input data is dual-species, please select to run in PDX mode!")
      
      genes_present_converted <- which(rownames(counts_matrix) %in% converted[[species_sym]])
      tmp.counts <- counts_matrix[genes_present_converted, ]
    }
    
    species_genes <- master_ref
    species_converted_hg <- master_ref
    
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
    #dataset_pco <- ortho_class_Data[ortho_class_Data$Gene_Type == "protein-coding", "Freq"]
    #biomart_ortho_pco_mouse <- species_hg_class[species_hg_class$Gene_Type == "protein-coding", "Freq"] 
    dataset_pco <- sum(ortho_class_Data[ortho_class_Data$Gene_Type %in% c("protein-coding", "protein_coding"), "Freq"])
    biomart_ortho_pco_mouse <- sum(species_hg_class[species_hg_class$Gene_Type %in% c("protein-coding", "protein_coding"), "Freq"])
    matched <- dataset_pco/biomart_ortho_pco_mouse * 100
    unmatched <- 100 - matched
    pie_data <- data.frame(
      category = c("Matched", "Unmatched"),
      count = c(matched,unmatched)
    )
    
    matched_genes <- intersect(rownames(seurat_obj), converted[[species_sym]])
    unmatched_genes <- setdiff(rownames(seurat_obj), matched_genes)
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
        server = TRUE,
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
    
    output$genes_list_ready <- reactive({
      return(!is.null(data_df))
    })
    outputOptions(output, 'genes_list_ready', suspendWhenHidden = FALSE)
    
    converted[[species_sym]] <- as.character(converted[[species_sym]])
    converted$HGNC.symbol <- as.character(converted$HGNC.symbol)
    
    rownames(tmp.counts) <- make.unique(plyr::mapvalues(
      x = as.character(rownames(tmp.counts)),
      from = as.character(converted[[species_sym]]),
      to = converted$HGNC.symbol,
      warn_missing = FALSE
    ))
    updated_assay_name <- paste0(assay, "_ortho")
    tmp.counts <- as(tmp.counts, "dgCMatrix")
    tmp <- CreateSeuratObject(counts = tmp.counts, assay = updated_assay_name)
    tmp@meta.data <- seurat_obj@meta.data
    tmp@reductions <- seurat_obj@reductions
    tmp@assays[[assay]] <- seurat_obj@assays[[assay]]
    
    if (assay == "Spatial") {
      tmp@images <- seurat_obj@images
    }
    
    convertedData(tmp)
  })
  
  output$status <- renderUI({
    if (!is.null(convertedData())) {
      tags$span("Conversion completed. You can now download the OrthologAL converted Seurat object.", style = "color: green;")
    } else {
      tags$span("Upload an RDS file. Once loaded, click 'Convert' to start.", style = "color: blue;")
    }
  })
  
  # --- DYNAMIC DOWNLOAD UI ---
  output$download_visibile_in_main_page <- renderUI({
    if (!is.null(convertedData())) {
      wellPanel(
        h5("Choose Download Destination:"),
        radioButtons("download_destination", label = NULL, 
                     choices = c("Local Computer" = "local", "Save to Server" = "server"), 
                     inline = TRUE),
        
        # Show standard download button if Local is selected
        conditionalPanel(
          condition = "input.download_destination == 'local'",
          downloadButton("downloadButton", "Download Converted Data", class = "btn btn-success btn-block mt-3")
        ),
        
        # Show shinyFiles save button if Server is selected
        conditionalPanel(
          condition = "input.download_destination == 'server'",
          shinySaveButton("save_server_btn", "Save to Server", "Save file as...", 
                          filetype = list(RDS = "rds"), class = "btn btn-primary btn-block mt-3")
        )
      )
    }
  })
  
  # --- LOCAL DOWNLOAD ---
  output$downloadButton <- downloadHandler(
    filename = function() {
      if (!is.null(convertedData()))
        paste0("OrthologAL_", uploaded_filename())
    },
    content = function(file) {
      if (!is.null(convertedData()))
        saveRDS(convertedData(), file)
    }
  )
  
  # --- SERVER SAVE ---
  observeEvent(input$save_server_btn, {
    req(!is.integer(input$save_server_btn))
    file_info <- parseSavePath(volumes, input$save_server_btn)
    req(nrow(file_info) > 0)
    
    save_path <- as.character(file_info$datapath[1])
    
    tryCatch({
      # Save the file to the chosen server path
      showNotification("Saving to server, please wait...", type = "default", duration = 3)
      saveRDS(convertedData(), save_path)
      showNotification(paste("Success! File saved to:", save_path), type = "message", duration = 10)
    }, error = function(e) {
      showNotification(paste("Error saving file:", e$message), type = "error", duration = 10)
    })
  })
}
