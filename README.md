# License agreement
1. The Board of Trustees of the Georgetown University provides OrthologAL software and code (“Service”) free of charge for non-commercial use only.Use of the Service by any commercial entity for any purpose, including research, is prohibited.
2. By using the Service, you agree to be bound by the terms of this Agreement. Please read it carefully.
3. You agree not to use the Service for commercial advantage, or in the course of for-profit activities.You agree not to use the Service on behalf of any organization that is not a non-profit organization. Commercial entities wishing to use this Service should
contact Georgetown University Office of Technology Commercialization.
4. THE SERVICE IS OFFERED “AS IS”, AND, TO THE EXTENT PERMITTED BY LAW, Georgetown MAKES NO REPRESENTATIONS AND EXTENDS NO WARRANTIES OF ANY KIND,EITHER EXPRESS OR IMPLIED. GEORGETOWN SHALL NOT BE LIABLE FOR ANY CLAIMS OR DAMAGES WITH RESPECT TO ANY LOSS OR OTHER CLAIM BY YOU OR ANY THIRD PARTY ON ACCOUNT OF, OR ARISING FROM THE USE OF THE SERVICE.
YOU HEREBY AGREE TO DEFEND AND INDEMNIFY GEORGETOWN, ITS TRUSTEES, EMPLOYEES, OFFICERS, STUDENTS, AGENTS, FACULTY, REPRESENTATIVES, AND VOLUNTEERS (“GEORGETOWN INDEMNITEES”) FROM ANY LOSS OR CLAIM ASSERTED AGAINST GEORGETOWN INDEMNITEES ARISING FROM YOUR USE OF THE SERVICE.
5. All rights not expressly granted to you in this Agreement are reserved and retained by Georgetown or its licensors or content providers. This Agreement provides no license under any patent.
6. You agree that this Agreement and any dispute arising under it is governed by the laws of the State of Washington DC, United States of America, applicable to agreements negotiated, executed, and performed within DC.
7. Subject to your compliance with the terms and conditions set forth in this Agreement, Georgetown grants you a revocable, non-exclusive, non-transferable right to access and make use of the Service.

# OrthologAL Package


## Table of Contents

1. [Introduction](#introduction)
2. [Features](#features)
3. [prerequisites](#prerequisites)
4. [Installation](#installation)
5. [Orthogonal App](#Orthogonal_App)
6. [Contributing](#contributing)
7. [Contact](#contact)
8. [License](#license)

## 1. Introduction
Researchers can input their single-cell or other high-dimensional gene expression data from any species, and OrthologAL will output a human ortholog-converted dataset for download and use.
To demonstrate the utility of this application, we characterized orthologous conversion in single-cell, single-nuclei, and spatial transcriptomic data derived from common pre-clinical models, including genetically engineered mouse models of medulloblastoma, patient-derived orthotopic xenografts of medulloblastoma, and mouse and rat models of spinal cord injury.
We show that OrthologAL can convert these data efficiently while preserving the dimensional architecture of the transcriptomic expression data. 
OrthologAL will be broadly useful for applying pre-clinical, high-dimensional transcriptomics data to functional small molecule predictions using existing human-annotated databases.

## 2. Features

- User friendly
- Standard Datatype format (Seurat)
- Can use single-cell,single-nuclei and spatial data
- Multiple assays can be used with this package
- PDX harmonization
- Quality control metrics

### 3. Prerequisites
Dependencies required 

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("biomaRt")
install.packages("shiny")
install.packages("dplyr")
install.packages("data.table")
install.packages('DT')
install.packages("ggplot2")
install.packages("bslib")
remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)
```
## 4.Installation

### Step 1: Clone the GitHub Repository
You can clone the repository using the git clone in your terminal, or download the ZIP folder from GitHub.
```git
git clone "https://github.com/AyadLab/OrthologAL.git" 
cd OrthologAL-main/ 
```
### Step 2: Install the package 
```r
#package installation
setwd("~/OrthologAL-main")
devtools::install()
```
### Step 3: Run the Shiny App
```r
library(shiny)
library(OrthologAL)
#Run the app
OrthologAL::RunOrthologAL()
```

### 5. OrthologAl App 
  
### Step 1: Upload the standard datatype format Seurat in the shiny app
![Screenshot 2024-10-25 at 2 56 30 PM](https://github.com/user-attachments/assets/1ad8b909-b898-4c5c-a3ce-48581fd035d3)


### Step 2: Select the species supported by the seurat object uploaded in step 1
- Select the species,by default we have mouse,rat,human and zebrafish species data
- There is a custom option to access other species through these identifiers Ensembl ID ( which extracts the information of genes based on the species code),
  attribute(gene symbol) acts a query  to look through the Ensembl data,to get all the genes associated with the species and filter is optional.

![Screenshot 2024-10-25 at 2 59 19 PM](https://github.com/user-attachments/assets/3a7b4f87-6a2c-47d7-8b42-ee566b65f380)


### Step 3: Select the assay supported by the seurat object uploaded in step 1

- Spatial,RNA,SCT and integrated assays are available here.

<img width="249" alt="Screenshot 2024-11-01 at 9 29 11 AM" src="https://github.com/user-attachments/assets/441d1914-f520-40a3-b2f2-4b51958086b2">

### Step 4: Select the ANimal model

- By default,no selection required here unless we are trying to use  Patient Derived Xenograft (PDX) model.
![Screenshot 2024-10-25 at 3 01 28 PM](https://github.com/user-attachments/assets/8f48486f-000f-4c14-8a81-2a27e1e39675)

### Step 5: Click on Convert
- Status changes to green(conversion done)

### Step 6 : Outputs
- Example data (Assay : RNA and species : RAT )
[![Download Data File](https://img.shields.io/badge/download-data--file-green)](https://drive.google.com/drive/folders/1icVieksEhdIUTEqkVSHZEfQxhfmKxU3m?usp=sharing)

- We get the pie chart which shows Biomart species genes with human orthologous gene Distribution and converted data gene distribution
![Artboard 1](https://github.com/user-attachments/assets/84704a2a-498c-4041-8a92-650797567329)
- Also shows protein coding genes overlap between BioMart Database and data uploaded
  ![Screenshot 2024-10-25 at 3 05 40 PM](https://github.com/user-attachments/assets/544ace45-db5c-4a60-9863-185148d20cc5)
- The system generates a summary table displaying both successfully matched genes and those without positive matches. This table is available to download as a .CSV file from within the application.

![Screenshot 2024-10-25 at 3 06 06 PM](https://github.com/user-attachments/assets/f8a77769-f008-4c84-b3a6-f43eaedc1641)


### 6. Contributing 
Rishika Chowdary, MS, Robert K. Suter, Ph.D., Matthew D’Antuono, BS, Cynthia Gomes, Ph.D., Joshua Stein, BA, 
Ki-Bum Lee, Ph.D., Jae K Lee, Ph.D., Nagi G. Ayad, Ph.D.

### 7. Contact
[na853@georgetown.edu ](mailto:na853@georgetown.edu) 



