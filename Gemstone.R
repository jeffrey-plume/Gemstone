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
      
      textInput(inputId = 'aid', label = 'Target Assay ID', value = "588834"),
      
      # Preprocessing controls
      checkboxInput('remove_na', 'Remove Rows with Missing Data', value = TRUE),
      selectInput('scale', 'Scaling Method', choices = c('None', 'Standardization', 'Min-Max Scaling')),
      checkboxInput('outlier_removal', 'Remove Outliers (Z-score ±3)', value = FALSE),

      numericInput('folds', 'Number of CV Folds', value = 5, min = 2, max = 10),
      numericInput('tuneLength', 'Tune Length', value = 3, min = 1, max = 10),

      # Action buttons
      actionButton('train', 'Train Model', style = 'background-color:blue; font-weight:bold; color:white'),
      fileInput(inputId = 'unknowns_in', label = 'Upload Unknowns in SDF Format', placeholder = "repurposing.sdf"),
      actionButton('unknowns_pred', 'Test Unknowns', style = 'background-color:blue; font-weight:bold; color:white'),

      downloadButton('downloadData', 'Download Predictions')
      
    ),

    
    mainPanel(
      

      
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
  
  
  retrieve_data <- reactive({
    req(input$aid)
    assay_data <- tryCatch({
      rpubchem::get.assay(input$aid) %>%
        dplyr::filter(!is.na(PUBCHEM.EXT.DATASOURCE.SMILES))
    }, error = function(e) {
      showNotification("Failed to retrieve assay data. Please check the target AID.", type = "error")
      return(NULL)
    })
    
    return(assay_data)
  })
  
  # Function to generate training data from molecular structures using rcdk
  generate_training_data <- reactive({
    req(retrieve_data())
    assay_data <- retrieve_data()
    
    # Parse SMILES strings to molecular structures
    mol <- rcdk::parse.smiles(assay_data$PUBCHEM.EXT.DATASOURCE.SMILES)
    
    # Generate molecular fingerprints
    fingerprints <- Rcpi::extractDrugExtendedComplete(mol) %>%
      data.frame() %>%
      set_names(paste("ECFP", 1:ncol(.)))
    
    # Generate molecular descriptors (e.g., Rule of 5)
    descriptors <- do.call(cbind, lapply(
      c('extractDrugApol', 'extractDrugTPSA', 'extractDrugHBondAcceptorCount', 
        'extractDrugHBondDonorCount', 'extractDrugWeight', 'extractDrugXLogP'), 
      function(f) eval(parse(text = paste0('Rcpi::', f, '(mol)'))))) %>%
      data.frame()
    
    # Combine fingerprints and descriptors into a single dataset
    training_data <- cbind(fingerprints, descriptors)
    
    return(training_data)
  })
  
  # Preprocess the training data
  preprocess_data <- reactive({
    data <- generate_training_data()  # Generate the molecular feature dataset
    req(data)
    
    # Apply user-specified preprocessing steps
    if (input$remove_na) {
      data <- na.omit(data)
    }
    
    if (input$scale != 'None') {
      num_columns <- sapply(data, is.numeric)
      if (input$scale == 'Standardization') {
        data[num_columns] <- scale(data[num_columns])
      } else if (input$scale == 'Min-Max Scaling') {
        data[num_columns] <- apply(data[num_columns], 2, function(x) (x - min(x)) / (max(x) - min(x)))
      }
    }
    
    if (input$outlier_removal) {
      data <- data %>%
        filter(apply(select_if(., is.numeric), 1, function(x) all(abs(scale(x)) < 3)))
    }
    
    return(data)
  })
  
  observeEvent(input$train, {
    req(preprocess_data())
    
    withProgress(message = 'Training model...', value = 0, {
      
      assay_data = retrieve_data()
      processed_data <- preprocess_data()
        
      incProgress(0.3)  # Update progress bar
      
      # Data partitioning
      index <- caret::createDataPartition(assay_data$PUBCHEM.ACTIVITY.OUTCOME, p = 0.75, list = FALSE)
      train <- processed_data[index, ]
      test <- processed_data[-index, ]
      
  
      cl <- parallel::makeCluster(parallel::detectCores() - 1)
      doParallel::registerDoParallel(cl)
      
      tryCatch({
        
        ctrl <- caret::trainControl(
          method = 'repeatedcv', 
          number = input$folds, 
          repeats = 3, 
          classProbs = TRUE,
          savePredictions = 'all')
        
        fit <- caret::train(
          x = train, 
          y = assay_data$PUBCHEM.ACTIVITY.OUTCOME[index], 
          method = 'rf', 
          trControl = ctrl, 
          tuneLength = input$tuneLength, 
          metric = 'Accuracy')
        
        saveRDS(fit, file = "trained_model.rds")
        
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
        
        pred <- predict(fit, test, type = 'prob')
        output$prediction <- DT::renderDT({
          DT::datatable(pred, extensions = 'Buttons', options = list(dom = 'Bfrtlp', buttons = c('copy', 'csv', 'excel')))
        })
        
        incProgress(1)  # Complete progress
      }, finally = {
        stopCluster(cl)
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

  
}

# Run the application 
shinyApp(ui = ui, server = server)

