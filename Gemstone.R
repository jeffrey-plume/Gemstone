options(shiny.maxRequestSize=500*1024^2) 
options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx12288m"))
library(shiny)
library(tidyverse)

ui <- fluidPage(
  
  titlePanel("Gemstone"),
  shinybusy::add_busy_spinner(spin = "fading-circle"),
  shinyjs::useShinyjs(),
  
  sidebarLayout(
    sidebarPanel(
      
      numericInput(inputId = 'aid', label = 'Target Assay ID', value = 1190),
      
      actionButton('load_data', 'Load Descriptors & Fingerprints', style = 'background-color:green; font-weight:bold; color:white'),
      
      # Preprocessing controls
      checkboxInput('remove_na', 'Impute means of missing data', value = TRUE),
      checkboxInput('remove_zeroVar', 'Remove near zero variance data', value = TRUE),
      
      selectInput('scale', 'Scaling Method', choices = c('None', 'Standardization', 'Min-Max Scaling')),
      checkboxInput('outlier_removal', 'Remove Outliers (Z-score ±3)', value = FALSE),

      numericInput('folds', 'Number of CV Folds', value = 5, min = 2, max = 10),
      numericInput('tuneLength', 'Tune Length', value = 3, min = 1, max = 10),
      
      # Sampling controls for handling imbalanced data
      radioButtons('sampling_method', 'Sampling Method for Imbalanced Data',
                   choices = c('None', 'Oversampling', 'Undersampling', 'Both'),
                   selected = 'None'),
      
      
      # Action buttons
      actionButton('train', 'Train Model', style = 'background-color:green; font-weight:bold; color:white'),
      HTML('<br>'),
      HTML('<br>'),
      
      fileInput(inputId = 'unknowns_in', label = 'Upload Unknowns in SDF Format', placeholder = "repurposing.sdf"),
      actionButton('unknowns_pred', 'Test Unknowns', style = 'background-color:blue; font-weight:bold; color:white'),

      downloadButton('downloadData', 'Download Predictions')
      
    ),

    
    mainPanel(
      

        shinyjs::hidden(
          div(id = "advanced",
              
              HTML('<br>'),
              textOutput('assay_description')
          )
        ),
        a(id = "toggleAdvanced", "Show Description"),
        
      HTML('<br>'),
      HTML('<br>'),
      
      div(DT::dataTableOutput(outputId = 'overview'), style = 'overflow-x: auto; white-space: nowrap; text-overflow: ellipsis;'),
      uiOutput('error'),
      
      HTML('<br>'),
      
      div(DT::dataTableOutput(
        outputId = 'prop'), style = 'overflow-x: auto; white-space: nowrap; text-overflow: ellipsis;'),
      HTML('<br>'),
      
      div(DT::dataTableOutput(outputId = 'corr'), style = 'overflow-x: auto; white-space: nowrap; text-overflow: ellipsis;'),
      
      HTML('<br>'),
      
      div(DT::dataTableOutput(outputId = 'prop_corr'), style = 'overflow-x: auto; white-space: nowrap; text-overflow: ellipsis;'),
      HTML('<br>'),
      
      div(DT::dataTableOutput(outputId = 'summary'), style = 'overflow-x: auto; white-space: nowrap; text-overflow: ellipsis;'),
      HTML('<br>'),
      
      plotOutput(outputId = 'roc'),
      HTML('<br>'),
      
      plotOutput(outputId = 'imp_out'),
      HTML('<br>'),
      
      div(DT::dataTableOutput(outputId = 'imp_table'), style = 'overflow-x: auto; white-space: nowrap; text-overflow: ellipsis;'),
      HTML('<br>'),
      
      div(DT::dataTableOutput(outputId = 'prediction'), style = 'overflow-x: auto; white-space: nowrap; text-overflow: ellipsis;'),
      HTML('<br>'),
      
      div(DT::dataTableOutput(outputId = 'unknowns_out'), style = 'overflow-x: auto; white-space: nowrap; text-overflow: ellipsis;'),
      HTML('<br>')
    )
  )
)

server <- function(input, output, session) {
  
  rv <- reactiveValues(index = NULL,
                       data_loaded = FALSE,
                       assay_data = NULL,
                       training_data = NULL,
                       fit = NULL)
  
  
  shinyjs::onclick("toggleAdvanced",
                   shinyjs::toggle(id = "advanced", anim = TRUE))   
  
  
  retrieve_data <- reactive({
    req(input$aid)

    
    assay_data <- tryCatch({
      rpubchem::get.assay(as.numeric(input$aid)) %>%
        dplyr::filter(!is.na(PUBCHEM.CID))
    }, error = function(e) {
      showNotification("Failed to retrieve assay data. Please check the target AID.", type = "error")
      return(NULL)
    })
    
    return(assay_data)
  }) %>% debounce(1000)
  
  observeEvent(input$load_data, {
    req(retrieve_data())
    
    showNotification("Loading data and generating fingerprints...", type = "message", duration = NULL)

    
    tryCatch({
      assay_data <- isolate(retrieve_data())
      
      cl <- parallel::makeCluster(parallel::detectCores() - 1)
      doParallel::registerDoParallel(cl)
      
      # Check if 'PUBCHEM.EXT.DATASOURCE.SMILES' column exists
      if ('PUBCHEM.EXT.DATASOURCE.SMILES' %in% colnames(assay_data)) {
        # Use SMILES strings from the assay data
        mol <- rcdk::parse.smiles(assay_data$PUBCHEM.EXT.DATASOURCE.SMILES)
      } else {
        # Retrieve SMILES using CIDs
        cids <- map(seq(ceiling(length(assay_data$PUBCHEM.CID) / 500)), 
                    ~str_c(
                      assay_data$PUBCHEM.CID[((.x - 1) * 500 + 1):ifelse(length(assay_data$PUBCHEM.CID) > (.x * 500), (.x * 500), length(assay_data$PUBCHEM.CID))], 
                      collapse = ","))
        
        url <- map(cids, ~paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/", .x, "/property/canonicalsmiles/csv"))
        smiles <- unlist(map(url, ~read_csv(httr::content(httr::GET(.x), "text", encoding = "UTF-8"))$CanonicalSMILES))
        
        mol <- rcdk::parse.smiles(as.vector(smiles))
      }

      # Generate molecular descriptors (e.g., Rule of 5)
      descriptors <- do.call(cbind, map(
        c('extractDrugExtendedComplete', 'extractDrugApol', 'extractDrugTPSA', 'extractDrugHBondAcceptorCount', 
          'extractDrugHBondDonorCount', 'extractDrugWeight', 'extractDrugXLogP'), 
        function(f) eval(parse(text = paste0('Rcpi::', f, '(mol)'))))) %>%
        data.frame() %>%
        set_names(str_replace(colnames(.), "(X)(\\d)", "ECFP\\2"))
      
      parallel::stopCluster(cl)
      
      # Combine fingerprints and descriptors into a single dataset
      rv$training_data <- data.frame(cbind(PUBCHEM.ACTIVITY.OUTCOME = as.factor(assay_data$PUBCHEM.ACTIVITY.OUTCOME), descriptors)) 
      rv$data_loaded <- TRUE

      showNotification("Data and fingerprints loaded successfully!", type = "message")
      
    }, error = function(e) {
      showNotification("Error loading data or generating fingerprints. Please try again.", type = "error")
    })
  })
  # Function to generate training data from molecular structures using rcdk

  
  # Preprocess the training data
  preprocess_data <- reactive({
    req(rv$training_data)
    
    data <- rv$training_data # Generate the molecular feature dataset
    scales <- list()
    
    # Apply the sampling method chosen by the user
    if (input$sampling_method == 'Oversampling') {
      data <- ROSE::ovun.sample(PUBCHEM.ACTIVITY.OUTCOME ~ ., data = data, method = "over")$data
    } else if (input$sampling_method == 'Undersampling') {
      data <- ROSE::ovun.sample(PUBCHEM.ACTIVITY.OUTCOME ~ ., data = data, method = "under")$data
    } else if (input$sampling_method == 'Both') {
      data <- ROSE::ovun.sample(PUBCHEM.ACTIVITY.OUTCOME ~ ., data = data, method = "both")$data
    } else {
      data <- data  # No sampling applied
    }
    
    
    if (input$remove_zeroVar) {
      # impute means for numeric columns
      data <- data[,-caret::nearZeroVar(data)]
    }
    
    index <- caret::createDataPartition(data$PUBCHEM.ACTIVITY.OUTCOME, p = 0.75, list = FALSE)
    rv$index <- index
    
    if (input$remove_na) {
      # impute means for numeric columns
      data[, sapply(data, is.numeric)] <- map(data[, sapply(data, is.numeric)], ~ifelse(is.na(.x), median(.x[index], na.rm = T), .x))
      
    }
    
    
    if (input$scale == 'Standardization') {
      scales["center"] <- 'center'
      scales["scale"] <- 'scale'
      
    } else if (input$scale == 'Min-Max Scaling') {
      scales["range"] <- 'range'
    } else {
      scales <- NULL
    }
      
    if (!is.null(scales) & length(scales) > 0) {
      scales <- str_c(scales, collapse = ',')
      
      preProc  <- caret::preProcess(data[,sapply(data, is.numeric)])
      rv$preProc <- preProc
      data[,sapply(data, is.numeric)] <- predict(preProc, data[,sapply(data, is.numeric)])
    }

    if (input$outlier_removal) {
      data <- data %>%
        filter(apply(select_if(., is.numeric), 1, function(x) all(abs(scale(x)) < 3)))
    }

    return(data)
  })
  
  
  # Function to retrieve assay description
  retrieve_assay_description <- reactive({
    req(input$aid)
    
    Sys.sleep(1)
    
    description <- tryCatch({
      rpubchem::get.assay.desc(input$aid)
    }, error = function(e) {
      showNotification("Failed to retrieve assay description. Please check the target AID.", type = "error")
      return(NULL)
    })
    
    return(description)
  })
  
  # Output the assay description
  output$assay_description <- renderText({
    req(retrieve_assay_description())
    paste("Bioassay Description: ", retrieve_assay_description())
  })
  

  
  observeEvent(input$train, {
    req(preprocess_data())
    processed_data <- preprocess_data()

    withProgress(message = 'Training model...', value = 0, {

      incProgress(0.3)  # Update progress bar
      
      index = rv$index
      # Data partitioning
      train <- processed_data[index, ]
      test <- processed_data[-index, ]
      assay_data <- retrieve_data()
      
      cl <- parallel::makeCluster(parallel::detectCores() - 1)
      doParallel::registerDoParallel(cl)
      
      tryCatch({
        
        ctrl <- caret::trainControl(
          method = 'repeatedcv', 
          number = input$folds, 
          repeats = 3, 
          classProbs = TRUE,
          savePredictions = 'all')
        
        fit <- caret::train(PUBCHEM.ACTIVITY.OUTCOME ~ .,
                            data = train,
                            method = 'rf', 
                            trControl = ctrl, 
                            na.action = na.omit,
                            tuneLength = input$tuneLength, 
                            metric = 'Accuracy')
        
        incProgress(0.6)  # Update progress bar
        
        output$summary <- DT::renderDT({
          DT::datatable(fit$results, extensions = 'Buttons', filter = 'top', options = list(dom = 'Bfrtlp', buttons = c('copy', 'csv', 'excel')))
        })
        
        output$roc <- renderPlot({
          plot(fit, type = 'b')
        })
        
        importance <- caret::varImp(fit)
        output$imp_out <- renderPlot({
          plot(importance, top = 20)
        })
        
        probs <- predict(fit, test, type = 'prob') %>%
          set_names(paste("Probs", colnames(.)))
        
        pred <- cbind(assay_data[-index, sapply(assay_data, is.character)], probs)
        
        output$prediction <- DT::renderDT({
          DT::datatable(pred, extensions = 'Buttons', options = list(dom = 'Bfrtlp', buttons = c('copy', 'csv', 'excel')))
        })
        
        rv$fit <- fit
        
        incProgress(1)  # Complete progress
      }, finally = {
        parallel::stopCluster(cl)
      })
    })
  })
  


  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("predictions-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(pred, file)
    }
  )

  output$overview <- DT::renderDT({
    req(retrieve_data())
    DT::datatable(retrieve_data())
  })
  
  output$prop <- DT::renderDT({
    req(preprocess_data())
    
    
    DT::datatable(preprocess_data())
  },
  server = FALSE)


    
  # Preprocess and predict the uploaded unknown data
  
  observeEvent(input$unknowns_pred, {
    req(input$unknowns_in)
    req(rv$fit)
    
    fit <- rv$fit
    preprocessed <- preprocess_data()
    
    showNotification("Processing unknowns...", type = "message", duration = NULL)
    
    tryCatch({
      mol <- rcdk::load.molecules(input$unknowns_in$datapath)
      
      unknown_data <- do.call(cbind, map(
        c('extractDrugExtendedComplete', 'extractDrugApol', 'extractDrugTPSA', 'extractDrugHBondAcceptorCount', 
          'extractDrugHBondDonorCount', 'extractDrugWeight', 'extractDrugXLogP'), 
        function(f) eval(parse(text = paste0('Rcpi::', f, '(mol)'))))) %>%
        data.frame() %>%
        set_names(str_replace(colnames(.), "(X)(\\d)", "ECFP\\2"))
      
      unknown_data <- dplyr::select(unknown_data, any_of(colnames(preprocessed)))
      
      if (exists("rv$preProc")) {
        preProc <- rv$preProc
        unknown_data <- predict(preProc, unknown_data)
      }
      
      predictions <- predict(fit, unknown_data, type = "prob") %>%
        print()

      output$unknowns_out <- DT::renderDT({
        DT::datatable(predictions, extensions = 'Buttons', options = list(dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel')))
      })
      
      showNotification("Predictions generated successfully!", type = "message")
      
    }, error = function(e) {
      showNotification("Error processing unknowns. Please check the file format.", type = "error")
    })
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)

