# **Gemstone** - Predictive Modeling with Shiny

**Gemstone** (Good Enough Model, SomeTimes, Often, NEver) is a Shiny-based web application designed for predictive modeling on biological assay data using machine learning algorithms. This app allows users to train models, save and retrieve models from a database, and test new compound predictions.

## **Table of Contents**

1. [Features](#features)
2. [Installation](#installation)
3. [Usage](#usage)
4. [Model Storage](#model-storage)
5. [Technologies](#technologies)
6. [Contributing](#contributing)
7. [License](#license)

---

## **Features**

- **Load Data:** Load chemical compound descriptors and fingerprints from PubChem Assay IDs (AID).
- **Preprocessing:** Options for preprocessing data, including handling missing values, removing near-zero variance features, scaling, and outlier removal.
- **Model Training:** Train machine learning models using Random Forest and other algorithms with user-defined cross-validation settings.
- **Imbalanced Data Handling:** Option to oversample, undersample, or use both techniques to handle imbalanced datasets.
- **Model Storage:** Automatically save and retrieve trained models to/from a SQLite database, so that models can be reused across sessions.
- **Unknown Compound Prediction:** Test unknown compounds by uploading files in SDF format and making predictions based on a trained model.
- **Download Predictions:** Export the prediction results as CSV files for further analysis.

---

## **Installation**

### **Prerequisites**

Ensure you have the following installed on your machine:
- R (version 4.0 or higher)
- RStudio (optional but recommended)
- The following R packages:
  - `shiny`
  - `tidyverse`
  - `future`
  - `DBI`
  - `RSQLite`
  - `rcdk`
  - `Rcpi`
  - `DT`
  - `shinybusy`
  - `shinyjs`

To install the required packages, run the following R command:

```r
install.packages(c("shiny", "tidyverse", "future", "DBI", "RSQLite", "rcdk", "Rcpi", "DT", "shinybusy", "shinyjs"))
```

### **Clone the Repository**

```bash
git clone https://github.com/yourusername/gemstone.git
cd gemstone
```

### **Run the Application**

To run the Shiny application, use the following R command:

```r
library(shiny)
runApp("app.R")
```

This will launch the application in your default web browser.

---

## **Usage**

1. **Input Target Assay ID (AID):** Enter the PubChem Assay ID in the provided input box.
2. **Load Descriptors and Fingerprints:** Click the 'Load Descriptors & Fingerprints' button to retrieve and process data from PubChem.
3. **Preprocessing Options:** Choose options for imputing missing data, scaling, and handling outliers.
4. **Model Training:** Click 'Train Model' to train a machine learning model. The model will be saved to the SQLite database.
5. **Testing Unknowns:** Upload a compound in SDF format, and click 'Test Unknowns' to generate predictions based on the trained model.
6. **Download Predictions:** After testing unknowns, you can download the predictions as a CSV file.

---

## **Model Storage**

### **Save Model to Database**

Trained models are automatically serialized and saved to a local SQLite database (`models.db`). When training a model for a specific AID, it will be stored, allowing you to reuse the model in future sessions.

```r
save_model_to_db(model, aid)
```

### **Retrieve Model from Database**

If a model for a specific AID already exists, it will be automatically loaded when the AID is entered.

```r
retrieved_model <- load_model_from_db(aid)
```

---

## **Technologies**

- **Shiny**: A web application framework for R.
- **SQLite**: A lightweight database to store models for retrieval across sessions.
- **Caret**: A machine learning package in R used for training models.
- **rcdk & Rcpi**: R packages used for generating chemical descriptors and fingerprints.

---

## **Contributing**

Contributions are welcome! To contribute:

1. Fork this repository.
2. Create a new branch: `git checkout -b my-feature-branch`.
3. Make your changes and commit them: `git commit -m 'Add new feature'`.
4. Push to the branch: `git push origin my-feature-branch`.
5. Submit a pull request.

Please make sure to update tests as appropriate.

---

## **License**

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## **Credits**

This project relies on several R packages and datasets from PubChem. Special thanks to the R community for providing the tools and resources that made this application possible.

---

### **Contact**

For any inquiries or issues, please open an issue in this repository, or reach out to the project maintainers via GitHub.

---

