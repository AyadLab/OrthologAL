# OrthologAL Package

A package is needed for converting genes from different species to human while preserving the dimensional architecture of the transcriptomic expression data. OrthologAL will be broadly useful for applying pre-clinical, high-dimensional transcriptomics data to functional small molecule predictions using existing human-annotated databases.
allow direct comparison using pharmacogenomic databases.

## Table of Contents

1. [Introduction](#introduction)
2. [Features](#features)
3. [prerequisites](#prerequisites)
4. [Installation](#installation)
5. [OrthologAL App](#Orthogonal_App)
6. [Contributing](#contributing)
7. [Contact](#contact)
8. [License](#license)

## 1. Introduction

Orthologous gene sets are common in drug discovery but there are no direct means of comparing disease signatures from different species. A method is needed for converting genes from different species to the human annotated genome to allow direct comparison using pharmaco- genomic databases. BioMart server is a free, scalable, large database supporting  large data resources including Ensembl, UniProt, HapMap, Wormbase, and Gramene. From the BioMart server, we can access different species gene sets. Our app OrthogonAL leverages the BioMart server to access gene sets from multiple species and are extracted through common identifiers like ENSEMBL ID, gene symbol, and attribute.These common identifiers help us to find the species genes from their database, and subsequently map them with human genes. OrthogonAL can deal with different species, animal models and allow the easy use of Seurat objects with spatial, single-cell, or single-nuclei datasets.

## 2. Features

- User friendly
- Standard Datatype format (Seurat)
- Can use single-cell,single-nuclei and spatial data
- Multiple assays can be used with this package 

### 3.Prerequisites
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
library(shiny)
library(Seurat)
library(biomaRt)
library(data.table)
library(dplyr)
library(ggplot2)
library(DT)
library(bslib)

```
## 4.Installation

### Step 1: Clone the GitHub Repository
You can clone the repository using the git clone in your terminal, or download the ZIP folder from GitHub.
```git
git clone "https://github.com/AyadLab/OrthogonAL.git" 
cd OrthogonAL-main/ 
```
### Step 2: Install the package 
```r
#package installation
setwd("~/OrthogconAL-main")
devtools::install()
```
### Step 3: Run the Shiny App
```r
library(shiny)
library(OrthologAL)
#Run the app
OrthologAL::RunOrthologAL()
```

### 5.OrthologAL App Use case 
  
### Step 1: Upload the standard datatype format Seurat in the shiny app
![Screenshot 2024-10-25 at 2 56 30 PM](https://github.com/user-attachments/assets/d3b743b2-11cf-4f24-bcd0-9c8ab0e1405e)


### Step 2: Select the species supported by the seurat object uploaded in step 1

![Screenshot 2024-10-25 at 2 59 19 PM](https://github.com/user-attachments/assets/94debb81-9177-41dc-9531-55b0bf6bde64)

- Select the species,by default we have mouse,human and zebrafish species data
- There is a custom option to access other species through these identifiers Ensembl ID ( which extracts the information of genes based on the species code),attribute(gene symbol) acts a query  to look through the Ensembl data,to get all the genes associated with the species and filter is optional.
<img width="474" alt="Screenshot 2024-05-14 at 2 01 43 PM" src="https://github.com/AyadLab/OrthogonAL/assets/43522813/47dc3ad1-8afc-4378-b1de-b6bfb3cc8970">


### Step 3: Select the assay supported by the seurat object uploaded in step 1

- Spatial,RNA,SCT and integrated assays are available here.
<img width="474" alt="Screenshot 2024-05-14 at 2 02 04 PM" src="https://github.com/AyadLab/OrthogonAL/assets/43522813/ca12e73c-42bf-4723-b218-0a018b56f633">

### Step 4: Select dual species model
![Screenshot 2024-10-25 at 3 01 28 PM](https://github.com/user-attachments/assets/78f976f9-fc98-4612-b647-8c012a229d85)

- By default,no selection required here unless we are trying to use Patient Derived Xenograft(PDX) model.

### Step 5: Click on Convert
- Status changes to green(conversion done)
<img width="507" alt="Screenshot 2024-05-14 at 2 03 51 PM" src="https://github.com/AyadLab/OrthogonAL/assets/43522813/759cf4f8-eb33-4ff2-838d-593fd91f4802">

### Step 6 : Outputs
- Example data (Assay : RNA and species : RAT )
[![Download Data File](https://img.shields.io/badge/download-data--file-green)](https://drive.google.com/drive/folders/1icVieksEhdIUTEqkVSHZEfQxhfmKxU3m?usp=sharing)

- Distribution of transcripts present in the Rat SCI scRNA-seq dataset (left) and distribution of genes present after orthologous conversion using OrthologAL (right).
![Screenshot 2024-10-25 at 3 05 03 PM](https://github.com/user-attachments/assets/bf947c38-2200-471a-b363-d52f9b488e50)
- Pie chart showing the proportion of human orthologous genes profiled by the rat SCI dataset 
![Screenshot 2024-10-25 at 3 05 40 PM](https://github.com/user-attachments/assets/49639cae-9aae-42ef-8c05-4e3a529cd4d2)
- The system generates a summary table displaying both successfully matched genes and those without positive matches. This table is available to download as a .CSV file from within the application.
![Screenshot 2024-10-25 at 3 06 06 PM](https://github.com/user-attachments/assets/1bce8fab-3ac6-4a94-9d22-bff8a139fcfb)

### 6.OrthologAL App Screenshot
![Screenshot 2024-10-25 at 3 12 22 PM](https://github.com/user-attachments/assets/4329e9ff-887d-42d8-8f90-c1855f5d031c)

### 7.Contributing 
Rishika Chowdary, MS,Robert K. Suter, Ph.D., Matthew D’Antuono, BS, 
Cynthia Gomes, Ph.D., Joshua Stein, BA,
Ki-Bum Lee, Ph.D., Jae K Lee, Ph.D., Nagi G. Ayad, Ph.D.

### 8.Contact
vc514@georgetown.edu


