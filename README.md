# OrthologAL Package

A package is needed for converting genes from different species to human to allow direct comparison using pharmacogenomic databases.

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
library(OrthogonAL)
#Run the app
OrthogonAL::OrthoConvo()
```

### 5.OrthologAL App
<img width="521" alt="Screenshot 2024-05-14 at 2 00 23 PM" src="https://github.com/AyadLab/OrthogonAL/assets/43522813/6ba1f970-7652-48af-9601-a94a469ea133">
  
### Step 1: Upload the standard datatype format Seurat in the shiny app
<img width="483" alt="Screenshot 2024-05-14 at 2 01 08 PM" src="https://github.com/AyadLab/OrthogonAL/assets/43522813/505149ee-3503-40a6-a741-99782570786b">


### Step 2: Select the species supported by the seurat object uploaded in step 1
- Select the species,by default we have mouse,human and zebrafish species data

<img width="474" alt="Screenshot 2024-05-14 at 2 01 30 PM" src="https://github.com/AyadLab/OrthogonAL/assets/43522813/5e90e438-7651-4610-a40b-f210ffae3ccf">

- There is a custom option to access other species through these identifiers Ensembl ID ( which extracts the information of genes based on the species code),attribute(gene symbol) acts a query  to look through the Ensembl data,to get all the genes associated with the species and filter is optional.
<img width="474" alt="Screenshot 2024-05-14 at 2 01 43 PM" src="https://github.com/AyadLab/OrthogonAL/assets/43522813/47dc3ad1-8afc-4378-b1de-b6bfb3cc8970">


### Step 3: Select the assay supported by the seurat object uploaded in step 1

- Spatial,RNA,SCT and integrated assays are available here.
<img width="474" alt="Screenshot 2024-05-14 at 2 02 04 PM" src="https://github.com/AyadLab/OrthogonAL/assets/43522813/ca12e73c-42bf-4723-b218-0a018b56f633">

### Step 4: Select the ANimal model

- By default,no selection required here unless we are trying to use  Patient derived orthotopic xenograft (PDOX) model.

<img width="470" alt="Screenshot 2024-05-14 at 2 02 24 PM" src="https://github.com/AyadLab/OrthogonAL/assets/43522813/bb42ba4b-b530-48e3-b2a2-a9d4c025e698">

### Step 5: Click on Convert
- Status changes to green(conversion done)
<img width="507" alt="Screenshot 2024-05-14 at 2 03 51 PM" src="https://github.com/AyadLab/OrthogonAL/assets/43522813/759cf4f8-eb33-4ff2-838d-593fd91f4802">

### Step 6 : Outputs
- Example data (Assay : Spatial and species : Mouse)
[![Download Data File](https://img.shields.io/badge/download-data--file-green)](https://drive.google.com/drive/folders/1icVieksEhdIUTEqkVSHZEfQxhfmKxU3m?usp=sharing)

- We get the pie chart (which shows which genes are matched and unmatched with biomaRT genes list) and download the csv list
<img width="659" alt="Screenshot 2024-06-21 at 1 19 54 PM" src="https://github.com/AyadLab/OrthogonAL/assets/43522813/f55363eb-bf0e-4f0c-9ee7-578db22fbfb8">
<img width="772" alt="Screenshot 2024-06-21 at 1 20 21 PM" src="https://github.com/AyadLab/OrthogonAL/assets/43522813/7dd5830e-4307-4c96-8f35-1642c792ef1e">



### 6.Contributing 
Rishika Chowdary, MS,Robert K. Suter, Ph.D., Matthew D’Antuono, BS, 
Cynthia Gomes, Ph.D., Joshua Stein, BA,
Ki-Bum Lee, Ph.D., Jae K Lee, Ph.D., Nagi G. Ayad, Ph.D.

### 7.Contact
vc514@georgetown.edu


