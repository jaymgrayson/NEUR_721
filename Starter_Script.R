# This script is to get you ready for Wed's class
# Place cursor at the lefthand edge of line 5
# Then click the Run button
# This will run the script one line at a time
rm(list=ls())
if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install()
BiocManager::install("flowCore", version = "3.8")
BiocManager::install("flowTrans", version = "3.8")
install.packages("tidyverse")
install.packages("Rtsne")

# Set working directory to where the data file you downloaded is