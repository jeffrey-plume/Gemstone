options(shiny.maxRequestSize=500*1024^2) 
options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx12288m"))
library(shiny)
library(tidyverse)
library(future)
library(DBI)
library(RSQLite)
library(future)
library(promises)

# Create a SQLite database connection
con <- dbConnect(RSQLite::SQLite(), "models.db")

## Create a table to store models if it doesn't exist
dbExecute(con, "CREATE TABLE IF NOT EXISTS models (aid INTEGER PRIMARY KEY, model BLOB)")

# Function to save the model to the database
save_model_to_db <- function(model, aid) {
  # Serialize the model into a raw vector
  serialized_model <- serialize(model, connection = NULL)
  
  # Convert the serialized model to a raw vector before inserting it
  dbExecute(con, "INSERT OR REPLACE INTO models (aid, model) VALUES (?, ?)",
            params = list(aid, list(serialized_model)))
}


# Function to load the model from the database
load_model_from_db <- function(aid) {
  model_blob <- dbGetQuery(con, "SELECT model FROM models WHERE aid = ?", params = list(aid))
  
  if (nrow(model_blob) == 0) {
    return(NULL)  # No model found for this AID
  }
  
  unserialized_model <- unserialize(model_blob$model[[1]])
  return(unserialized_model)
}

# Close the database connection when the app is stopped
onStop(function() {
  dbDisconnect(con)
})

plan(multisession)

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
      actionButton('unknowns_pred', 'Test Unknowns', style = 'background-color:green; font-weight:bold; color:white'),

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
  
  rv <- reactiveValues(index = NULL, data_loaded = FALSE, assay_data = NULL, training_data = NULL, fit = NULL)
  
  shinyjs::onclick("toggleAdvanced", shinyjs::toggle(id = "advanced", anim = TRUE))
  
  retrieve_data <- reactive({
    req(input$aid)
    assay_data <- tryCatch({
      rpubchem::get.assay(as.numeric(input$aid)) %>%
        dplyr::filter(!is.na(PUBCHEM.CID))
    }, error = function(e) {
      shinyalert::shinyalert("Error", "Failed to retrieve assay data. Please check the target AID.", type = "error")
      return(NULL)
    })
    return(assay_data)
  }) %>% debounce(1000)
  
  
  observeEvent(input$aid, {
    aid <- input$aid
    
    # Load the model from the database if it exists
    model <- load_model_from_db(aid)
    
    if (!is.null(model)) {
      rv$fit <- model
      showNotification("Loaded pre-trained model for AID: ", type = "message")
    } else {
      rv$fit <- NULL
      showNotification("No pre-trained model found for AID. Please train the model.", type = "error")
    }
  })
  
  observeEvent(input$load_data, {
    req(retrieve_data())
    
    showNotification("Loading data and generating fingerprints...", type = "message", duration = NULL)
    
    future({
      assay_data <- isolate(retrieve_data())
      if ('PUBCHEM.EXT.DATASOURCE.SMILES' %in% colnames(assay_data)) {
        mol <- rcdk::parse.smiles(assay_data$PUBCHEM.EXT.DATASOURCE.SMILES)
      } else {
        cids <- map(seq(ceiling(length(assay_data$PUBCHEM.CID) / 500)), 
                    ~str_c(assay_data$PUBCHEM.CID[((.x - 1) * 500 + 1):min(.x * 500, length(assay_data$PUBCHEM.CID))], 
                           collapse = ","))
        url <- map(cids, ~paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/", .x, "/property/canonicalsmiles/csv"))
        smiles <- unlist(map(url, ~read_csv(httr::content(httr::GET(.x), "text", encoding = "UTF-8"))$CanonicalSMILES))
        mol <- rcdk::parse.smiles(as.vector(smiles))
      }
      descriptors <- do.call(cbind, map(
        c('extractDrugExtendedComplete', 'extractDrugApol', 'extractDrugTPSA', 'extractDrugHBondAcceptorCount', 
          'extractDrugHBondDonorCount', 'extractDrugWeight', 'extractDrugXLogP'), 
        function(f) eval(parse(text = paste0('Rcpi::', f, '(mol)'))))) %>%
        data.frame() %>%
        set_names(str_replace(colnames(.), "(X)(\\d)", "ECFP\\2")) 
      return(data.frame(cbind(PUBCHEM.ACTIVITY.OUTCOME = as.factor(assay_data$PUBCHEM.ACTIVITY.OUTCOME), descriptors)))
    }) %...>% {
      rv$training_data <- .
      rv$data_loaded <- TRUE
      showNotification("Data and fingerprints loaded successfully!", type = "message")
    } %...!% {
      shinyalert::shinyalert("Error", "Failed to load data or generate fingerprints.", type = "error")
    }
  })
        

  
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
    
    # Capture reactive inputs before starting future
    folds <- input$folds
    tuneLength <- input$tuneLength
    aid <- input$aid
    

    index <- caret::createDataPartition(processed_data$PUBCHEM.ACTIVITY.OUTCOME, p = 0.75, list = FALSE)
    rv$index <- index
    train <- processed_data[index, ]
    test <- processed_data[-index, ]
    
    # Asynchronously train the model using future
    future({
      tryCatch({
        ctrl <- caret::trainControl(method = 'repeatedcv', 
                                    number = folds,  # Use the captured value
                                    repeats = 3, 
                                    classProbs = TRUE, 
                                    savePredictions = 'all')
        
        fit <- caret::train(PUBCHEM.ACTIVITY.OUTCOME ~ ., 
                            data = train, 
                            method = 'rf',
                            trControl = ctrl, 
                            tuneLength = tuneLength,  # Use the captured value
                            metric = 'Accuracy')
        
        
        fit  # Return the model
      }, error = function(e) {
        message("Error during model training: ", e$message)
        stop(e)  # Re-throw the error to the outer handler
      })
    }, seed = TRUE) %...>% {
      # Handle the successful result (when the future is resolved)
      rv$fit <- .
      save_model_to_db(., aid)
      
      output$summary <- DT::renderDT({
        DT::datatable(rv$fit$results, 
                  extensions = 'Buttons', 
                  options = list(
                    dom = 'Bfrtip', 
                    buttons = c('copy', 'csv', 'excel'),
                    serverSide = TRUE  # Enable server-side processing
                  ))
      }, server = TRUE)
      output$roc <- renderPlot({
        plot(rv$fit, type = 'b')
      })
      importance <- caret::varImp(rv$fit)
      output$imp_out <- renderPlot({
        plot(importance, top = 20)
      })
      showNotification("Model training completed successfully!", type = "message")
    } %...!% {
      # Handle the error case (when the future encounters an error)
      showNotification("Failed to train the model.", type = "error")
    } %>% finally(function() {
      # Ensure progress is completed
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
      
      row.names(unknown_data) <- unlist(map(mol, ~rcdk::get.property(.x, 'smiles')))

      unknown_data <- dplyr::select(unknown_data, any_of(colnames(preprocessed)))
      
      if (exists("rv$preProc")) {
        preProc <- rv$preProc
        unknown_data <- predict(preProc, unknown_data)
      }
      
      predictions <- predict(fit, unknown_data, type = "prob") 
      
      row.names(predictions) <- row.names(unknown_data)

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

