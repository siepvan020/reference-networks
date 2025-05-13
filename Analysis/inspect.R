



setwd("/div/pythagoras/u1/siepv/siep/Analysis/Rdata")

blood_platelet <- readRDS("temp_platelet_bloodScorpionOutput.Rdata")
lung_pneu <- readRDS("temp_type_ii_pneumocyte_lungScorpionOutput.Rdata")

blood_cd4 <- readRDS("temp_cd4_positive_t_cell_bloodScorpionOutput.Rdata")
lung_cd4 <- readRDS("temp_cd4_positive_t_cell_lungScorpionOutput.Rdata")

blood_cd8 <- readRDS("temp_cd8_positive_t_cell_bloodScorpionOutput.Rdata")
lung_cd8 <- readRDS("temp_cd8_positive_t_cell_lungScorpionOutput.Rdata")

blood_mono <- readRDS("temp_monocyte_bloodScorpionOutput.Rdata")
lung_mono <- readRDS("temp_monocyte_lungScorpionOutput.Rdata")

v1_indegrees <- data.frame(gene = blood_platelet$regNet@Dimnames[[2]])

v1_indegrees$blood_platelet <- colSums(as.matrix(blood_platelet$regNet))
v1_indegrees$lung_pneu <- colSums(as.matrix(lung_pneu$regNet))
v1_indegrees$blood_cd4 <- colSums(as.matrix(blood_cd4$regNet))
v1_indegrees$lung_cd4 <- colSums(as.matrix(lung_cd4$regNet))
v1_indegrees$blood_cd8 <- colSums(as.matrix(blood_cd8$regNet))
v1_indegrees$lung_cd8 <- colSums(as.matrix(lung_cd8$regNet))
v1_indegrees$blood_mono <- colSums(as.matrix(blood_mono$regNet))
v1_indegrees$lung_mono <- colSums(as.matrix(lung_mono$regNet))


v1_cor_matrix <- cor(v1_indegrees[, -1])
