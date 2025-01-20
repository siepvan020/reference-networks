#!/usr/bin/env Rscript

library(reshape2)
library(ggplot2)
library(tidyverse)

setwd("/div/pythagoras/u1/siepv/siep/Analysis")

motif <- read.delim("data/priors/motif_prior_names_2024.tsv", header = FALSE, sep = "\t")

blood_cd4 <- readRDS("Rdata/temp_cd4_positive_t_cell_bloodScorpionOutput.Rdata")
blood_platelet <- readRDS("Rdata/temp_platelet_bloodScorpionOutput.Rdata")

# Outdegree stuff
blood_cd4_outdegree <- rowSums(blood_cd4$regNet)
hist(blood_cd4_outdegree, breaks = 25, main = "Histogram of blood_cd4 - Outdegree", xlab = "Outdegree")
top_regulators_blood_cd4 <- sort(blood_cd4_outdegree, decreasing = TRUE)[1:10]
as.data.frame(top_regulators_blood_cd4)

blood_platelet_outdegree <- rowSums(blood_platelet$regNet)
hist(blood_platelet_outdegree, breaks = 25, main = "Histogram of blood_platelet - Outdegree", xlab = "Outdegree")
top_regulators_blood_platelet <- sort(blood_platelet_outdegree, decreasing = TRUE)[1:10]
as.data.frame(top_regulators_blood_platelet)




# Indegree stuff
blood_cd4_indegree <- colSums(blood_cd4$regNet)
hist(blood_cd4_indegree, breaks = 25, main = "Histogram of blood_cd4 - Indegree", xlab = "Indegree")
top_genes_blood_cd4 <- sort(blood_cd4_indegree, decreasing = TRUE)[1:10]
as.data.frame(top_genes_blood_cd4)

blood_platelet_indegree <- colSums(blood_platelet$regNet)
hist(blood_platelet_indegree, breaks = 25, main = "Histogram of blood_platelet - Indegree", xlab = "Indegree")
top_genes_blood_platelet <- sort(blood_platelet_indegree, decreasing = TRUE)[1:10]
as.data.frame(top_genes_blood_platelet)




rs_blood_cd4 <- reshape2::melt(as.matrix(blood_cd4$regNet))
temp <- dplyr::left_join(motif, rs_blood_cd4, by = c("V1" = "Var1", "V2" = "Var2"))

ggplot(rs_blood_cd4, aes(x = value)) +
    geom_histogram(binwidth = 0.1, fill = "blue", color = "#000000") +
    labs(title = "Histogram of blood_cd4", x = "Value", y = "Frequency")



lung_cd4 <- readRDS("Rdata/temp_cd4_positive_t_cell_lungScorpionOutput.Rdata")
lung_pneumo <- readRDS("Rdata/temp_type_ii_pneumocyte_lungScorpionOutput.Rdata")

# Outdegree stuff
lung_cd4_outdegree <- rowSums(lung_cd4$regNet)
hist(lung_cd4_outdegree, breaks = 25, main = "Histogram of lung_cd4 - Outdegree", xlab = "Outdegree")
top_regulators_lung_cd4 <- sort(lung_cd4_outdegree, decreasing = TRUE)[1:10]
as.data.frame(top_regulators_lung_cd4)

lung_pneumo_outdegree <- rowSums(lung_pneumo$regNet)
hist(lung_pneumo_outdegree, breaks = 25, main = "Histogram of lung_pneumo - Outdegree", xlab = "Outdegree")
top_regulators_lung_pneumo <- sort(lung_pneumo_outdegree, decreasing = TRUE)[1:10]
as.data.frame(top_regulators_lung_pneumo)



# Indegree stuff
lung_cd4_indegree <- colSums(lung_cd4$regNet)
hist(lung_cd4_indegree, breaks = 25, main = "Histogram of lung_cd4 - Indegree", xlab = "Indegree")
top_genes_lung_cd4 <- sort(lung_cd4_indegree, decreasing = TRUE)[1:10]
as.data.frame(top_genes_lung_cd4)

lung_pneumo_indegree <- colSums(lung_pneumo$regNet)
hist(lung_pneumo_indegree, breaks = 25, main = "Histogram of lung_pneumo - Indegree", xlab = "Indegree")
top_genes_lung_pneumo <- sort(lung_pneumo_indegree, decreasing = TRUE)[1:10]
as.data.frame(top_genes_lung_pneumo)


# blood_cd8 <- readRDS("Rdata/temp_cd8_positive_t_cell_bloodScorpionOutput.Rdata")
# lung_cd8 <- readRDS("Rdata/temp_cd8_positive_t_cell_lungScorpionOutput.Rdata")

# blood_mono <- readRDS("Rdata/temp_monocyte_bloodScorpionOutput.Rdata")
# lung_mono <- readRDS("Rdata/temp_monocyte_lungScorpionOutput.Rdata")

# blood_platelet <- readRDS("Rdata/temp_platelet_bloodScorpionOutput.Rdata")
# lung_pneumo <- readRDS("Rdata/temp_type_ii_pneumocyte_lungScorpionOutput.Rdata")


rs_lung_cd4 <- reshape2::melt(as.matrix(lung_cd4$regNet))



ggplot(rs_lung_cd4, aes(x = value)) +
    geom_histogram(binwidth = 0.1, fill = "red", color = "#00000000") +
    labs(title = "Histogram of lung_cd4", x = "Value", y = "Frequency")
