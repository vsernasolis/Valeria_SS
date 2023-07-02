# Honors Thesis Ratfish Code
#Revised May 03, 2023 VSS

#Load relevant packages
library(geomorph) #v.4.0.3
library(Morpho) #v.2.9
library(nlme) #v.3.1-157 

#Set working directory
setwd("~/Desktop/Ratfish.R") #Valeria
setwd("~/Desktop/Current Projects/Ratfish/RCode-Ratfish/Final") #Brandon

#Do L1 vs TL and L1 vs TW once for each structure 
####Adult Clasper Tenacula Data####
maleClasperShape <- read.morphologika("claspermorphologikano11_unscaled.txt") #Load male shape data
      maleClasperGPA <- gpagen(maleClasperShape, 
                               ProcD = FALSE) #Run Generalized Procrustes Analysis using bending energy criterion
      
      #Change to clasperMaleData
      claspermaleData <- read.csv("clasper2mratfish.csv") #Load in data with specimenID, L1 (total length), TL, TW, AG
      
      #male PCA to examine shape trends in morphospace
      maleClasperPCA <- gm.prcomp(maleClasperGPA$coords)
      maleClasperPCASum <- summary(maleClasperPCA) #Gives you breakdown of variance and eigenvalues
      
      #maleClasperPCASum$PC.summary and maleClasperPCA$x are exported as supplemental file XX
      
      #Set a color palette for the PCA based on L1 where larger males are red and smaller males are blue
      pal = colorRampPalette(c("blue", "red"))
      maleClasperGPA$col <- pal(10)[as.numeric(cut(log10(claspermaleData$L1), breaks = 10))] #This breaks into 10 colors
      
      #FIGURE XX
      maleClasperPCAPlot <- plot(maleClasperPCA,
                                 xlab = paste("Principal Component 1 ", "(", sep = "", 
                                              paste(round(maleClasperPCASum$PC.summary$Comp1[2]*100, digits = 2), "%", ")", sep = "")),
                                 ylab = paste("Principal Component 2 ", "(", sep = "", 
                                              paste(round(maleClasperPCASum$PC.summary$Comp2[2]*100, digits = 2), "%", ")", sep = "")),
                                 pch = 21,
                                 col = "black",
                                 bg = maleClasperGPA$col,
                                 cex = 4,
                                 cex.lab = 1.2)
      text(maleClasperPCAPlot$PC.points[ , 2] ~ maleClasperPCAPlot$PC.points[ , 1], 
           labels = claspermaleData$Order, cex= 1, col = c("white"))
      
      #generate a color gradient for an inset in the plot
      plot(rep(1,100),col=(pal(100)), pch=19,cex=10)
      
      #Make Shape Extremes with Point Clouds using Morpho package
      procClasperMale <- procSym(maleClasperShape)
      plot(procClasperMale$PCscores[,2] ~ procClasperMale$PCscores[,1]) #replot in morpho just to double check axes did not flip
      
      #Plot at ends of PC1 axes, negative is purple and positive is green, Magnified by 2 to accentuate trends
      palette(c("white", "green", "purple"))
      posPC1claspermale <- restoreShapes(2*sd(procClasperMale$PCscores[,1]), 
                                         procClasperMale$PCs[,1], 
                                         procClasperMale$mshape)
      negPC1claspermale <- restoreShapes(-2*sd(procClasperMale$PCscores[,1]), 
                                         procClasperMale$PCs[,1], 
                                         procClasperMale$mshape)
      deformGrid3d(posPC1claspermale, 
                   negPC1claspermale, 
                   ngrid = 0, 
                   lines = FALSE) #Type no
      
      #and again for PC2
      posPC2claspermale <- restoreShapes(2*sd(procClasperMale$PCscores[,2]), 
                                         procClasperMale$PCs[,2], 
                                         procClasperMale$mshape)
      negPC2claspermale <- restoreShapes(-2*sd(procClasperMale$PCscores[,2]), 
                                         procClasperMale$PCs[,2], 
                                         procClasperMale$mshape)
      deformGrid3d(posPC2claspermale, 
                   negPC2claspermale, 
                   ngrid = 0, 
                   lines = FALSE) #Type no
      
      #Make a geomorph dataframe for using functions in geomorph package 
      maleClasperGDF <-geomorph.data.frame(shape = maleClasperGPA$coords,
                                           cs = log10(maleClasperGPA$Csize),
                                           L1 = log10(claspermaleData$L1),
                                           TL = log10(claspermaleData$TL),
                                           TW = log10(claspermaleData$TW)) #Note that Csize is log10 transformed for downstream analyses
      
      #Examine relationship between L1 and clasper centroid size (log10 transformed)
      clasperDFMale <- data.frame(log10(claspermaleData$L1),
                                  as.numeric(log10(maleClasperGPA$Csize)),
                                  log10(claspermaleData$TL),
                                  log10(claspermaleData$TW))
      
      colnames(clasperDFMale) <- c("MaleL1", "MaleCS", "TL", "TW") #Name columns
      
      
#Look for relationships between L1 and CS and testes dimensions using gls      
      #L1CSFit a gls model and plot the relationship
      clasperL1CSFitMale <- gls(MaleCS ~ MaleL1, 
                         data = clasperDFMale)
      
      plot(MaleCS ~ MaleL1, 
           data = clasperDFMale,
           pch = 21,
           col = "black",
           bg = maleClasperGPA$col,
           cex = 4,
           cex.lab = 1.2,
           xlim = c(1.55, 1.64), #NEED TO CHANGE ACROSS VARIABLES
           ylim = c(2.78, 2.97), #NEED TO CHANGE ACROSS VARIABLES
           xlab = "log10 (L1)",
           ylab = "log10 (Centroid Size)")
      abline(a = coef(clasperL1CSFitMale)[1],
             b = coef(clasperL1CSFitMale)[2],
             lw = 3)
      
      summary(clasperL1CSFitMale) #Check if the confidence intervals include 1 to assess allometric trajectory
      confint(clasperL1CSFitMale)
      
      #CSTLFit a gls model and plot the relationship
      clasperCSTLFitMale <- gls(MaleCS ~ TL, 
                         data = clasperDFMale)
      
      plot(MaleCS ~ TL, 
           data = clasperDFMale,
           pch = 21,
           col = "black",
           bg = maleClasperGPA$col,
           cex = 4,
           cex.lab = 1.2,
           #xlim = c(1,3), #NEED TO CHANGE ACROSS VARIABLES
           ylim = c(2.78, 2.97), #NEED TO CHANGE ACROSS VARIABLES
           xlab = "log10 (TL)",
           ylab = "log10 (Centroid Size)")
      abline(a = coef(clasperCSTLFitMale)[1],
             b = coef(clasperCSTLFitMale)[2],
             lw = 3)
      
      summary(clasperCSTLFitMale) #Check if the confidence intervals include 1 to assess allometric trajectory
      confint(clasperCSTLFitMale)
      
      #CSTWFit a gls model and plot the relationship
      clasperCSTWFitMale <- gls(MaleCS ~ TW, 
                         data = clasperDFMale)
      
      plot(MaleCS ~ TW, 
           data = clasperDFMale,
           pch = 21,
           col = "black",
           bg = maleClasperGPA$col,
           cex = 4,
           cex.lab = 1.2,
           xlim = c(0.05, 0.35), #NEED TO CHANGE ACROSS VARIABLES
           ylim = c(2.75, 3), #NEED TO CHANGE ACROSS VARIABLES
           xlab = "log10 (TW)",
           ylab = "log10 (Centroid Size)")
      abline(a = coef(clasperCSTWFitMale)[1],
             b = coef(clasperCSTWFitMale)[2],
             lw = 3)
      
      summary(clasperCSTWFitMale) #Check if the confidence intervals include 1 to assess allometric trajectory
      confint(clasperCSTWFitMale)
      
      #L1TLFit a gls model and plot the relationship
      clasperL1TLFitMale <- gls(TL ~ MaleL1, 
                                data = clasperDFMale)
      
      plot(TL ~ MaleL1, 
           data = clasperDFMale,
           pch = 21,
           col = "black",
           bg = maleClasperGPA$col,
           cex = 4,
           cex.lab = 1.2,
           xlim = c(1.55,1.64), #NEED TO CHANGE ACROSS VARIABLES
           ylim = c(0.2,0.6), #NEED TO CHANGE ACROSS VARIABLES
           xlab = "log10 (L1)",
           ylab = "log10 (TL)")
      abline(a = coef(clasperL1TLFitMale)[1],
             b = coef(clasperL1TLFitMale)[2],
             lw = 3)
      
      summary(clasperL1TLFitMale) #Check if the confidence intervals include 1 to assess allometric trajectory
      confint(clasperL1TLFitMale)
      
      #L1TWFit a gls model and plot the relationship
      clasperL1TWFitMale <- gls(TW ~ MaleL1, 
                                data = clasperDFMale)
      
      plot(TW ~ MaleL1, 
           data = clasperDFMale,
           pch = 21,
           col = "black",
           bg = maleClasperGPA$col,
           cex = 4,
           cex.lab = 1.2,
           xlim = c(1.55,1.64), #NEED TO CHANGE ACROSS VARIABLES
           ylim = c(0,0.4), #NEED TO CHANGE ACROSS VARIABLES
           xlab = "log10 (L1)",
           ylab = "log10 (TW)")
      abline(a = coef(clasperL1TWFitMale)[1],
             b = coef(clasperL1TWFitMale)[2],
             lw = 3)
      
      summary(clasperL1TWFitMale) #Check if the confidence intervals include 1 to assess allometric trajectory
      confint(clasperL1TWFitMale)
#Put it all together in a table
claspCS <- matrix(, nrow = 3, ncol = 3)
      colnames(claspCS) <- c("p", "lowerCI", "upperCI")
      rownames(claspCS) <- c("L1", "TL", "TW")
      
      claspCS[[1, 1]] <- round(summary(clasperL1CSFitMale)$tTable[2, 4], 3)
      claspCS[[1, 2]] <- round(confint(clasperL1CSFitMale)[2, 1], 3)
      claspCS[[1, 3]] <- round(confint(clasperL1CSFitMale)[2, 2], 3)
      
      claspCS[[2, 1]] <- round(summary(clasperCSTLFitMale)$tTable[2, 4], 3)
      claspCS[[2, 2]] <- round(confint(clasperCSTLFitMale)[2, 1], 3)
      claspCS[[2, 3]] <- round(confint(clasperCSTLFitMale)[2, 2], 3)
      
      claspCS[[3, 1]] <- round(summary(clasperCSTWFitMale)$tTable[2, 4], 3)
      claspCS[[3, 2]] <- round(confint(clasperCSTWFitMale)[2, 1], 3)
      claspCS[[3, 3]] <- round(confint(clasperCSTWFitMale)[2, 2], 3)
    
claspCS
      
      
claspL1 <- matrix(, nrow = 3, ncol = 3)
colnames(claspL1) <- c("p", "lowerCI", "upperCI")
rownames(claspL1) <- c("L1", "TL", "TW")

claspL1[[1, 1]] <- round(summary(clasperL1CSFitMale)$tTable[2, 4], 3)
claspL1[[1, 2]] <- round(confint(clasperL1CSFitMale)[2, 1], 3)
claspL1[[1, 3]] <- round(confint(clasperL1CSFitMale)[2, 2], 3)

claspL1[[2, 1]] <- round(summary(clasperL1TLFitMale)$tTable[2, 4], 3)
claspL1[[2, 2]] <- round(confint(clasperL1TLFitMale)[2, 1], 3)
claspL1[[2, 3]] <- round(confint(clasperL1TLFitMale)[2, 2], 3)

claspL1[[3, 1]] <- round(summary(clasperL1TWFitMale)$tTable[2, 4], 3)
claspL1[[3, 2]] <- round(confint(clasperL1TWFitMale)[2, 1], 3)
claspL1[[3, 3]] <- round(confint(clasperL1TWFitMale)[2, 2], 3)
      
claspL1
#Look for relationships between shape and L1, CS, and testes dimensions using Procrustes ANOVAs         
      #Do a Procrustes ANOVA of shape and L1
      clasperMaleFit <- procD.lm(shape ~ L1,
                                 data = maleClasperGDF,
                                 iter = 1000)
      clasperMaleFitSum <- summary(clasperMaleFit) 
      
      #And with CS of the claspers as in the males
      clasperMaleFit2 <- procD.lm(shape ~ cs,
                                  data = maleClasperGDF,
                                  iter = 1000)
      clasperMaleFitSum2 <- summary(clasperMaleFit2) 
      
      #Do a Procrustes ANOVA of shape and TL
      clasperMaleFit3 <- procD.lm(shape ~ TL,
                                  data = maleClasperGDF,
                                  iter = 1000)
      clasperMaleFitSum3 <- summary(clasperMaleFit3) 
      
      #Do a Procrustes ANOVA of shape and TW
      clasperMaleFit4 <- procD.lm(shape ~ TW,
                                  data = maleClasperGDF,
                                  iter = 1000)
      clasperMaleFitSum4 <- summary(clasperMaleFit4) 
      
#Move relevant info into a table 
claspShapeMat <- matrix(, nrow = 4, ncol = 2)
      colnames(claspShapeMat) <- c("R^2", "p")
      rownames(claspShapeMat) <- c("L1", "CS", "TL", "TW")
      
      claspShapeMat[[1, 1]] <- round(clasperMaleFitSum$table$Rsq[[1]], 3)
      claspShapeMat[[1, 2]] <- round(clasperMaleFitSum$table$`Pr(>F)`[[1]], 3)
      claspShapeMat[[2, 1]] <- round(clasperMaleFitSum2$table$Rsq[[1]], 3)
      claspShapeMat[[2, 2]] <- round(clasperMaleFitSum2$table$`Pr(>F)`[[1]], 3)
      claspShapeMat[[3, 1]] <- round(clasperMaleFitSum3$table$Rsq[[1]], 3)
      claspShapeMat[[3, 2]] <- round(clasperMaleFitSum3$table$`Pr(>F)`[[1]], 3)
      claspShapeMat[[4, 1]] <- round(clasperMaleFitSum4$table$Rsq[[1]], 3)
      claspShapeMat[[4, 2]] <- round(clasperMaleFitSum4$table$`Pr(>F)`[[1]], 3)
      
write.csv(claspShapeMat, "claspShapeMat.csv") #Export file for thesis
      
####Juvenile Claspers#### 
#These data were aligned separately from the adult claspers so they are not directly comparable with the adults.
      jclasperShape <- read.morphologika("juvclasper2morphologika_unscaled.txt") #Load male shape data
      jclasperGPA <- gpagen(jclasperShape, 
                            ProcD = FALSE) #Run Generalized Procrustes Analysis using bending energy criterion
      
      #Change to clasperMaleData
      jmaleData <- read.csv("jmratfish.csv") #Load in data with specimenID, L1 (total length), TL, TW, AG without HYCO010M
      
      #male PCA to examine shape trends in morphospace
      jclasperPCA <- gm.prcomp(jclasperGPA$coords)
      jclasperPCASum <- summary(jclasperPCA) #Gives you breakdown of variance and eigenvalues
      
      #maleClasperPCASum$PC.summary and maleClasperPCA$x are exported as supplemental file XX
      
      #Set a color palette for the PCA based on L1 where larger males are red and smaller males are blue
      pal = colorRampPalette(c("blue", "red"))
      jclasperGPA$col <- pal(10)[as.numeric(cut(log10(jmaleData$L1), breaks = 10))] #This breaks into 10 colors
      
      #FIGURE XX
      jclasperPCAPlot <- plot(jclasperPCA,
                              xlab = paste("Principal Component 1 ", "(", sep = "", 
                                           paste(round(jclasperPCASum$PC.summary$Comp1[2]*100, digits = 2), "%", ")", sep = "")),
                              ylab = paste("Principal Component 2 ", "(", sep = "", 
                                           paste(round(jclasperPCASum$PC.summary$Comp2[2]*100, digits = 2), "%", ")", sep = "")),
                              pch = 21,
                              col = "black",
                              bg = jclasperGPA$col,
                              cex = 4,
                              cex.lab = 1.2)
      text(jclasperPCAPlot$PC.points[ , 2] ~ jclasperPCAPlot$PC.points[ , 1], 
           labels = jmaleData$Order, cex= 1, col = c("white"))
      
      #generate a color gradient for an inset in the plot
      plot(rep(1,100),col=(pal(100)), pch=19,cex=10)
      
      #Make Shape Extremes with Point Clouds using Morpho package
      procjclasperMale <- procSym(jclasperShape)
      plot(procjclasperMale$PCscores[,2] ~ procjclasperMale$PCscores[,1])
      
      #Plot at ends of PC1 axes, negative is purple and positive is green, Magnified by 2 to accentuate trends
      palette(c("white", "green", "purple"))
      posPC1jclaspermale <- restoreShapes(2*sd(procjclasperMale$PCscores[,1]), 
                                          procjclasperMale$PCs[,1], 
                                          procjclasperMale$mshape)
      negPC1jclaspermale <- restoreShapes(-2*sd(procjclasperMale$PCscores[,1]), 
                                          procjclasperMale$PCs[,1], 
                                          procjclasperMale$mshape)
      deformGrid3d(posPC1jclaspermale, 
                   negPC1jclaspermale, 
                   ngrid = 0, 
                   lines = FALSE) #Type no
      
      #and again for PC2
      posPC2jclaspermale <- restoreShapes(2*sd(procjclasperMale$PCscores[,2]), 
                                          procjclasperMale$PCs[,2], 
                                          procjclasperMale$mshape)
      negPC2jclaspermale <- restoreShapes(-2*sd(procjclasperMale$PCscores[,2]), 
                                          procjclasperMale$PCs[,2], 
                                          procjclasperMale$mshape)
      deformGrid3d(posPC2jclaspermale, 
                   negPC2jclaspermale, 
                   ngrid = 0, 
                   lines = FALSE) #Type no
      
      #Make a geomorph dataframe for using functions in geomorph package 
      jclasperGDF <-geomorph.data.frame(shape = jclasperGPA$coords,
                                        cs = log10(jclasperGPA$Csize),
                                        L1 = log10(jmaleData$L1),
                                        TL = log10(jmaleData$TL),
                                        TW = log10(jmaleData$TW)) #Note that Csize is log10 transformed for downstream analyses
      
      #Examine relationship between L1 and clasper centroid size (log10 transformed)
      jclasperDFMale <- data.frame(log10(jmaleData$L1),
                                   as.numeric(log10(jclasperGPA$Csize)),
                                   log10(jmaleData$TL),
                                   log10(jmaleData$TW))
      
      colnames(jclasperDFMale) <- c("MaleL1", "MaleCS", "TL", "TW") #Name columns
      
      #L1CSFit a gls model and plot the relationship
      jclasperL1CSFitMale <- gls(MaleCS ~ MaleL1, 
                                 data = jclasperDFMale)
      
      plot(MaleCS ~ MaleL1, 
           data = jclasperDFMale,
           pch = 21,
           col = "black",
           bg = jclasperGPA$col,
           cex = 4,
           cex.lab = 1.2,
           xlim = c(1.46,1.58), #NEED TO CHANGE ACROSS VARIABLES
           ylim = c(2,2.5), #NEED TO CHANGE ACROSS VARIABLES
           xlab = "log10 (L1)",
           ylab = "log10 (Centroid Size)")
      abline(a = coef(jclasperL1CSFitMale)[1],
             b = coef(jclasperL1CSFitMale)[2],
             lw = 3)
      
      summary(jclasperL1CSFitMale) #Check if the confidence intervals include 1 to assess allometric trajectory
      confint(jclasperL1CSFitMale)
      
      #CSTLFit a gls model and plot the relationship
      jclasperCSTLFitMale <- gls(MaleCS ~ TL, 
                                 data = jclasperDFMale)
      
      plot(MaleCS ~ TL, 
           data = jclasperDFMale,
           pch = 21,
           col = "black",
           bg = jclasperGPA$col,
           cex = 4,
           cex.lab = 1.2,
           xlim = c(-0.4,0.08), #NEED TO CHANGE ACROSS VARIABLES
           ylim = c(2,2.5), #NEED TO CHANGE ACROSS VARIABLES
           xlab = "log10(TL)",
           ylab = "log10 (Centroid Size)")
      abline(a = coef(jclasperCSTLFitMale)[1],
             b = coef(jclasperCSTLFitMale)[2],
             lw = 3)
      
      summary(jclasperCSTLFitMale) #Check if the confidence intervals include 1 to assess allometric trajectory
      confint(jclasperCSTLFitMale)
      
      #CSTWFit a gls model and plot the relationship
      jclasperCSTWFitMale <- gls(MaleCS ~ TW, 
                                 data = jclasperDFMale)
      
      plot(MaleCS ~ TW, 
           data = jclasperDFMale,
           pch = 21,
           col = "black",
           bg = jclasperGPA$col,
           cex = 4,
           cex.lab = 1.2,
           xlim = c(-0.8,0), #NEED TO CHANGE ACROSS VARIABLES
           ylim = c(2.1,2.45), #NEED TO CHANGE ACROSS VARIABLES
           xlab = "log10(TW)",
           ylab = "log10 (Centroid Size)")
      abline(a = coef(jclasperCSTWFitMale)[1],
             b = coef(jclasperCSTWFitMale)[2],
             lw = 3)
      
      summary(jclasperCSTWFitMale) #Check if the confidence intervals include 1 to assess allometric trajectory
      confint(jclasperCSTWFitMale)
      
      #L1TLFit a gls model and plot the relationship
      jclasperL1TLFitMale <- gls(TL ~ MaleL1, 
                                data = jclasperDFMale)
      
      plot(TL ~ MaleL1, 
           data = jclasperDFMale,
           pch = 21,
           col = "black",
           bg = jclasperGPA$col,
           cex = 4,
           cex.lab = 1.2,
           xlim = c(1.47,1.58), #NEED TO CHANGE ACROSS VARIABLES
           ylim = c(-0.5,0.1), #NEED TO CHANGE ACROSS VARIABLES
           xlab = "log10 (L1)",
           ylab = "log10 (TL)")
      abline(a = coef(jclasperL1TLFitMale)[1],
             b = coef(jclasperL1TLFitMale)[2],
             lw = 3)
      
      summary(jclasperL1TLFitMale) #Check if the confidence intervals include 1 to assess allometric trajectory
      confint(jclasperL1TLFitMale)
      
      #L1TWFit a gls model and plot the relationship
      jclasperL1TWFitMale <- gls(TW ~ MaleL1, 
                                data = jclasperDFMale)
      
      plot(TW ~ MaleL1, 
           data = jclasperDFMale,
           pch = 21,
           col = "black",
           bg = jclasperGPA$col,
           cex = 4,
           cex.lab = 1.2,
           xlim = c(1.47,1.58), #NEED TO CHANGE ACROSS VARIABLES
           ylim = c(-1,0.2), #NEED TO CHANGE ACROSS VARIABLES
           xlab = "log10 (L1)",
           ylab = "log10 (TW)")
      abline(a = coef(jclasperL1TWFitMale)[1],
             b = coef(jclasperL1TWFitMale)[2],
             lw = 3)
      
      summary(jclasperL1TWFitMale) #Check if the confidence intervals include 1 to assess allometric trajectory
      confint(jclasperL1TWFitMale)
      
#Put it all together in a table
      jClaspCS <- matrix(, nrow = 3, ncol = 3)
      colnames(jClaspCS) <- c("p", "lowerCI", "upperCI")
      rownames(jClaspCS) <- c("L1", "TL", "TW")
      
      jClaspCS[[1, 1]] <- round(summary(jclasperL1CSFitMale)$tTable[2, 4], 3)
      jClaspCS[[1, 2]] <- round(confint(jclasperL1CSFitMale)[2, 1], 3)
      jClaspCS[[1, 3]] <- round(confint(jclasperL1CSFitMale)[2, 2], 3)
      
      jClaspCS[[2, 1]] <- round(summary(jclasperCSTLFitMale)$tTable[2, 4], 3)
      jClaspCS[[2, 2]] <- round(confint(jclasperCSTLFitMale)[2, 1], 3)
      jClaspCS[[2, 3]] <- round(confint(jclasperCSTLFitMale)[2, 2], 3)
      
      jClaspCS[[3, 1]] <- round(summary(jclasperCSTWFitMale)$tTable[2, 4], 3)
      jClaspCS[[3, 2]] <- round(confint(jclasperCSTWFitMale)[2, 1], 3)
      jClaspCS[[3, 3]] <- round(confint(jclasperCSTWFitMale)[2, 2], 3)
      
      #Do a Procrustes ANOVA of shape and L1
      jclasperFit <- procD.lm(shape ~ L1,
                              data = jclasperGDF,
                              iter = 1000)
      jClaspFitSum <- summary(jclasperFit) 
      
      #And with CS of the claspers as in the males
      jclasperFit2 <- procD.lm(shape ~ cs,
                               data = jclasperGDF,
                               iter = 1000)
      jClaspFitSum2 <- summary(jclasperFit2) 
      
      #Do a Procrustes ANOVA of shape and TL
      jclasperFit3 <- procD.lm(shape ~ TL,
                               data = jclasperGDF,
                               iter = 1000)
      jClaspFitSum3 <- summary(jclasperFit3) 
      
      #Do a Procrustes ANOVA of shape and TW
      jclasperFit4 <- procD.lm(shape ~ TW,
                               data = jclasperGDF,
                               iter = 1000)
      jClaspFitSum4 <- summary(jclasperFit4)  
      
#Move relevant info into a table 
      jClaspShapeMat <- matrix(, nrow = 4, ncol = 2)
      colnames(jClaspShapeMat) <- c("R^2", "p")
      rownames(jClaspShapeMat) <- c("L1", "CS", "TL", "TW")
      
      jClaspShapeMat[[1, 1]] <- round(jClaspFitSum$table$Rsq[[1]], 3)
      jClaspShapeMat[[1, 2]] <- round(jClaspFitSum$table$`Pr(>F)`[[1]], 3)
      jClaspShapeMat[[2, 1]] <- round(jClaspFitSum2$table$Rsq[[1]], 3)
      jClaspShapeMat[[2, 2]] <- round(jClaspFitSum2$table$`Pr(>F)`[[1]], 3)
      jClaspShapeMat[[3, 1]] <- round(jClaspFitSum3$table$Rsq[[1]], 3)
      jClaspShapeMat[[3, 2]] <- round(jClaspFitSum3$table$`Pr(>F)`[[1]], 3)
      jClaspShapeMat[[4, 1]] <- round(jClaspFitSum4$table$Rsq[[1]], 3)
      jClaspShapeMat[[4, 2]] <- round(jClaspFitSum4$table$`Pr(>F)`[[1]], 3)
      

    
####Adult Prepelvic Tenacula Data####    
    malePPTShape <- read.morphologika("pptmorphologikaNo11_unscaled.txt") #Load male shape data
    malePPTGPA <- gpagen(malePPTShape, 
                         ProcD = FALSE) #Run GPA using bending energy criterion
    
    PPTmaleData <- read.csv("ppt2mratfish.csv") #Load in data with specimenID, L2cm (length 2 cm)
    
    #male PCA to examine shape trends in morphospace
    malePPTPCA <- gm.prcomp(malePPTGPA$coords)
    malePPTPCASum <- summary(malePPTPCA) #Gives you breakdown of variance and eigenvalues
    
    #malePPTPCASum$PC.summary and malePPTPCA$x are exported as supplemental file XX
    
    #Set a color palette for the PCA based on L1 where larger males are red and smaller males are blue
    pal = colorRampPalette(c("blue", "red"))
    malePPTGPA$col <- pal(10)[as.numeric(cut(log10(PPTmaleData$L1), breaks = 10))]
    
    #FIGURE XX
    malePPTPCAPlot <- plot(malePPTPCA,
                           xlab = paste("Principal Component 1 ", "(", sep = "",
                                        paste(round(malePPTPCASum$PC.summary$Comp1[2]*100, digits = 2), "%", ")", sep = "")),
                           ylab = paste("Principal Component 2 ", "(", sep = "", 
                                        paste(round(malePPTPCASum$PC.summary$Comp2[2]*100, digits = 2), "%", ")", sep = "")),
                           pch = 21,
                           col = "black",
                           bg = malePPTGPA$col,
                           cex = 4,
                           cex.lab = 1.2)
    text(malePPTPCAPlot$PC.points[ , 2] ~ malePPTPCAPlot$PC.points[ , 1], 
         labels = PPTmaleData$Order, cex= 1, col = c("white"))
    
    #Make Shape Extremes with Point Clouds using Morpho package
    procPPTMale <- procSym(malePPTShape)
    plot(procPPTMale$PCscores[,2] ~ procPPTMale$PCscores[,1])
    
    #Plot at ends of PC1 axes, negative is purple and positive is green, Magnified by 2 to accentuate trends
    palette(c("white", "green", "purple"))
    posPC1PPTmale <- restoreShapes(2*sd(procPPTMale$PCscores[,1]), 
                                   procPPTMale$PCs[,1], 
                                   procPPTMale$mshape)
    negPC1PPTmale <- restoreShapes(-2*sd(procPPTMale$PCscores[,1]), 
                                   procPPTMale$PCs[,1], 
                                   procPPTMale$mshape)
    deformGrid3d(posPC1PPTmale, 
                 negPC1PPTmale, 
                 ngrid = 0, 
                 lines = FALSE) #Type no
    
    #and again for PC2
    posPC2PPTmale <- restoreShapes(2*sd(procPPTMale$PCscores[,2]), 
                                   procPPTMale$PCs[,2], 
                                   procPPTMale$mshape)
    negPC2PPTmale <- restoreShapes(-2*sd(procPPTMale$PCscores[,2]), 
                                   procPPTMale$PCs[,2], 
                                   procPPTMale$mshape)
    deformGrid3d(posPC2PPTmale, 
                 negPC2PPTmale, 
                 ngrid = 0, 
                 lines = FALSE) #Type no
    
    #Make a geomorph dataframe for using functions in geomorph package 
    malePPTGDF <-geomorph.data.frame(shape = malePPTGPA$coords,
                                     cs = log10(malePPTGPA$Csize),
                                     L1 = log10(PPTmaleData$L1),
                                     TL = log10(PPTmaleData$TL),
                                     TW = log10(PPTmaleData$TW))
    
    #Make dataframe to compare log10 L1 and log10CS
    PPTDFMale <- data.frame(log10(PPTmaleData$L1),
                            as.numeric(log10(malePPTGPA$Csize)),
                            log10(PPTmaleData$TL),
                            log10(PPTmaleData$TW))
    colnames(PPTDFMale) <- c("MaleL1", "MaleCS", "TL", "TW")
    
    #Fit a GLS model comparing CS and L1
    PPTL1CSFitMale <- gls(MaleCS ~ MaleL1, 
                       data = PPTDFMale)
    
    plot(MaleCS ~ MaleL1, 
         data = PPTDFMale,
         pch = 21,
         col = "black",
         bg = malePPTGPA$col,
         cex = 4,
         cex.lab = 1.2,
         xlim = c(1.575, 1.64), #NEED TO CHANGE ACROSS VARIABLES
         ylim = c(2.45, 2.7), #NEED TO CHANGE ACROSS VARIABLES
         xlab = "log10 (L1)",
         ylab = "log10 (Centroid Size)")
    abline(a = coef(PPTL1CSFitMale)[1],
           b = coef(PPTL1CSFitMale)[2],
           lw = 3)
    
    summary(PPTL1CSFitMale) #Check if the confidence intervals include 1 to assess allometric trajectory
    confint(PPTL1CSFitMale)
    
    #CSTLFit a gls model and plot the relationship
    PPTCSTLFitMale <- gls(MaleCS ~ TL, 
                       data = PPTDFMale)
    
    plot(MaleCS ~ TL, 
         data = PPTDFMale,
         pch = 21,
         col = "black",
         bg = malePPTGPA$col,
         cex = 4,
         cex.lab = 1.2,
         xlim = c(0.3,0.52), #NEED TO CHANGE ACROSS VARIABLES
         ylim = c(2.45, 2.7), #NEED TO CHANGE ACROSS VARIABLES
         xlab = "log10 (TL)",
         ylab = "log10 (Centroid Size)") 
    abline(a = coef(PPTCSTLFitMale)[1],
           b = coef(PPTCSTLFitMale)[2],
           lw = 3)
    
    summary(PPTCSTLFitMale) #Check if the confidence intervals include 1 to assess allometric trajectory
    confint(PPTCSTLFitMale)
    
    #CSTWFit a gls model and plot the relationship
    PPTCSTWFitMale <- gls(MaleCS ~ TW, 
                       data = PPTDFMale)
    
    plot(MaleCS ~ TW, 
         data = PPTDFMale,
         pch = 21,
         col = "black",
         bg = malePPTGPA$col,
         cex = 4,
         cex.lab = 1.2,
         xlim = c(0.07, 0.35), #NEED TO CHANGE ACROSS VARIABLES
         ylim = c(2.47, 2.65), #NEED TO CHANGE ACROSS VARIABLES
         xlab = "log10 (TW)",
         ylab = "log10 (Centroid Size)")
    abline(a = coef(PPTCSTWFitMale)[1],
           b = coef(PPTCSTWFitMale)[2],
           lw = 3)
    
    summary(PPTCSTWFitMale) #Check if the confidence intervals include 1 to assess allometric trajectory
    confint(PPTCSTWFitMale)
    
    
#Put it all together in a table
PPTCS <- matrix(, nrow = 3, ncol = 3)
    colnames(PPTCS) <- c("p", "lowerCI", "upperCI")
    rownames(PPTCS) <- c("L1", "TL", "TW")
    
    PPTCS[[1, 1]] <- round(summary(PPTL1CSFitMale)$tTable[2, 4], 3)
    PPTCS[[1, 2]] <- round(confint(PPTL1CSFitMale)[2, 1], 3)
    PPTCS[[1, 3]] <- round(confint(PPTL1CSFitMale)[2, 2], 3)
    
    PPTCS[[2, 1]] <- round(summary(PPTCSTLFitMale)$tTable[2, 4], 3)
    PPTCS[[2, 2]] <- round(confint(PPTCSTLFitMale)[2, 1], 3)
    PPTCS[[2, 3]] <- round(confint(PPTCSTLFitMale)[2, 2], 3)
    
    PPTCS[[3, 1]] <- round(summary(PPTCSTWFitMale)$tTable[2, 4], 3)
    PPTCS[[3, 2]] <- round(confint(PPTCSTWFitMale)[2, 1], 3)
    PPTCS[[3, 3]] <- round(confint(PPTCSTWFitMale)[2, 2], 3)
    
    #Do a Procrustes ANOVA of shape and L1
    PPTmaleFit <- procD.lm(shape ~ L1,
                           data = malePPTGDF,
                           iter = 1000)
    PPTMaleFitSum <- summary(PPTmaleFit) 
    
    #And with CS of the PPT 
    PPTmaleFit2 <- procD.lm(shape ~ cs,
                            data = malePPTGDF,
                            iter = 1000)
    PPTMaleFitSum2 <- summary(PPTmaleFit2) 
    #And with TL of the PPT    
    PPTmaleFit3 <- procD.lm(shape ~ TL,
                            data = malePPTGDF,
                            iter = 1000)
    PPTMaleFitSum3 <- summary(PPTmaleFit3) 
    #And with TW of the PPT    
    PPTmaleFit4 <- procD.lm(shape ~ TW,
                            data = malePPTGDF,
                            iter = 1000)
    PPTMaleFitSum4 <- summary(PPTmaleFit4) 
    
#Move relevant info into a table 
PPTShapeMat <- matrix(, nrow = 4, ncol = 2)
    colnames(PPTShapeMat) <- c("R^2", "p")
    rownames(PPTShapeMat) <- c("L1", "CS", "TL", "TW")
    
    PPTShapeMat[[1, 1]] <- round(PPTMaleFitSum$table$Rsq[[1]], 3)
    PPTShapeMat[[1, 2]] <- round(PPTMaleFitSum$table$`Pr(>F)`[[1]], 3)
    PPTShapeMat[[2, 1]] <- round(PPTMaleFitSum2$table$Rsq[[1]], 3)
    PPTShapeMat[[2, 2]] <- round(PPTMaleFitSum2$table$`Pr(>F)`[[1]], 3)
    PPTShapeMat[[3, 1]] <- round(PPTMaleFitSum3$table$Rsq[[1]], 3)
    PPTShapeMat[[3, 2]] <- round(PPTMaleFitSum3$table$`Pr(>F)`[[1]], 3)
    PPTShapeMat[[4, 1]] <- round(PPTMaleFitSum4$table$Rsq[[1]], 3)
    PPTShapeMat[[4, 2]] <- round(PPTMaleFitSum4$table$`Pr(>F)`[[1]], 3)
    
    
####Juvenile Pre-pelvic Tenaculum####
    jmalePPTShape <- read.morphologika("jmpptmorphologika_unscaled.txt") #Load male shape data
    jmalePPTGPA <- gpagen(jmalePPTShape, 
                          ProcD = FALSE) #Run Generalized Procrustes Analysis using bending energy criterion
    
    #Change to clasperMaleData
    jmaleData <- read.csv("jmratfish.csv") #Load in data with specimenID, L1 (total length), TL, TW, AG without HYCO010M
    
    #male PCA to examine shape trends in morphospace
    jmalePPTPCA <- gm.prcomp(jmalePPTGPA$coords)
    jmalePPTPCASum <- summary(jmalePPTPCA) #Gives you breakdown of variance and eigenvalues
    
    #maleClasperPCASum$PC.summary and maleClasperPCA$x are exported as supplemental file XX
    
    #Set a color palette for the PCA based on L1 where larger males are red and smaller males are blue
    pal = colorRampPalette(c("blue", "red"))
    jmalePPTGPA$col <- pal(10)[as.numeric(cut(log10(jmaleData$L1), breaks = 10))] #This breaks into 10 colors
    
    #FIGURE XX
    jmalePPTPCAPlot <- plot(jmalePPTPCA,
                            xlab = paste("Principal Component 1 ", "(", sep = "", 
                                         paste(round(jmalePPTPCASum$PC.summary$Comp1[2]*100, digits = 2), "%", ")", sep = "")),
                            ylab = paste("Principal Component 2 ", "(", sep = "", 
                                         paste(round(jmalePPTPCASum$PC.summary$Comp2[2]*100, digits = 2), "%", ")", sep = "")),
                            pch = 21,
                            col = "black",
                            bg = jmalePPTGPA$col,
                            cex = 4,
                            cex.lab = 1.2)
    text(jmalePPTPCAPlot$PC.points[ , 2] ~ jmalePPTPCAPlot$PC.points[ , 1], 
         labels = jmaleData$Order, cex= 1, col = c("white"))
    
    #generate a color gradient for an inset in the plot
    plot(rep(1,100),col=(pal(100)), pch=19,cex=10)
    
    #Make Shape Extremes with Point Clouds using Morpho package
    procjPPTMale <- procSym(jmalePPTShape)
    plot(procjPPTMale$PCscores[,2] ~ procjPPTMale$PCscores[,1])
    
    #Plot at ends of PC1 axes, negative is purple and positive is green, Magnified by 2 to accentuate trends
    palette(c("white", "green", "purple"))
    posPC1jPPTmale <- restoreShapes(2*sd(procjPPTMale$PCscores[,1]), 
                                    procjPPTMale$PCs[,1], 
                                    procjPPTMale$mshape)
    negPC1jPPTmale <- restoreShapes(-2*sd(procjPPTMale$PCscores[,1]), 
                                    procjPPTMale$PCs[,1], 
                                    procjPPTMale$mshape)
    deformGrid3d(posPC1jPPTmale, 
                 negPC1jPPTmale, 
                 ngrid = 0, 
                 lines = FALSE) #Type no
    
    #and again for PC2
    posPC2jPPTmale <- restoreShapes(2*sd(procjPPTMale$PCscores[,2]), 
                                    procjPPTMale$PCs[,2], 
                                    procjPPTMale$mshape)
    negPC2jPPTmale <- restoreShapes(-2*sd(procjPPTMale$PCscores[,2]), 
                                    procjPPTMale$PCs[,2], 
                                    procjPPTMale$mshape)
    deformGrid3d(posPC2jPPTmale, 
                 negPC2jPPTmale, 
                 ngrid = 0, 
                 lines = FALSE) #Type no
    
    #Make a geomorph dataframe for using functions in geomorph package 
    jmalePPTGDF <-geomorph.data.frame(shape = jmalePPTGPA$coords,
                                      cs = log10(jmalePPTGPA$Csize),
                                      L1 = log10(jmaleData$L1),
                                      TL = log10(jmaleData$TL),
                                      TW = log10(jmaleData$TW)) #Note that Csize is log10 transformed for downstream analyses
    
    #Examine relationship between L1 and clasper centroid size (log10 transformed)
    jmalePPTDFMale <- data.frame(log10(jmaleData$L1),
                                 as.numeric(log10(jmalePPTGPA$Csize)),
                                 log10(jmaleData$TL),
                                 log10(jmaleData$TW))
    
    colnames(jmalePPTDFMale) <- c("MaleL1", "MaleCS", "TL", "TW") #Name columns
    
    #L1CSFit a gls model and plot the relationship
    jmPPTL1CSFitMale <- gls(MaleCS ~ MaleL1, 
                            data = jmalePPTDFMale)
    
    plot(MaleCS ~ MaleL1, 
         data = jmalePPTDFMale,
         pch = 21,
         col = "black",
         bg = jmalePPTGPA$col,
         cex = 4,
         cex.lab = 1.2,
         xlim = c(1.47, 1.58), #NEED TO CHANGE ACROSS VARIABLES
         ylim = c(2.2, 2.6), #NEED TO CHANGE ACROSS VARIABLES
         xlab = "log10 (L1)",
         ylab = "log10 (Centroid Size)")
    abline(a = coef(jmPPTL1CSFitMale)[1],
           b = coef(jmPPTL1CSFitMale)[2],
           lw = 3)
    
    summary(jmPPTL1CSFitMale) #Check if the confidence intervals include 1 to assess allometric trajectory
    confint(jmPPTL1CSFitMale)
    
    #CSTLFit a gls model and plot the relationship
    jmPPTCSTLFitMale <- gls(MaleCS ~ TL, 
                            data = jmalePPTDFMale)
    
    plot(MaleCS ~ TL, 
         data = jmalePPTDFMale,
         pch = 21,
         col = "black",
         bg = jmalePPTGPA$col,
         cex = 4,
         cex.lab = 1.2,
         xlim = c(-0.4,0.1), #NEED TO CHANGE ACROSS VARIABLES
         ylim = c(2.2,2.6), #NEED TO CHANGE ACROSS VARIABLES
         xlab = "log10(TL)",
         ylab = "log10 (Centroid Size)")
    abline(a = coef(jmPPTCSTLFitMale)[1],
           b = coef(jmPPTCSTLFitMale)[2],
           lw = 3)
    
    summary(jmPPTCSTLFitMale) #Check if the confidence intervals include 1 to assess allometric trajectory
    confint(jmPPTCSTLFitMale)
    
    #CSTWFit a gls model and plot the relationship
    jmPPTCSTWFitMale <- gls(MaleCS ~ TW, 
                            data = jmalePPTDFMale)
    
    plot(MaleCS ~ TW, 
         data = jmalePPTDFMale,
         pch = 21,
         col = "black",
         bg = jmalePPTGPA$col,
         cex = 4,
         cex.lab = 1.2,
         xlim = c(-0.8,0), #NEED TO CHANGE ACROSS VARIABLES
         ylim = c(2.2,2.6), #NEED TO CHANGE ACROSS VARIABLES
         xlab = "log10(TW)",
         ylab = "log10 (Centroid Size)")
    abline(a = coef(jmPPTCSTWFitMale)[1],
           b = coef(jmPPTCSTWFitMale)[2],
           lw = 3)
    
    summary(jmPPTCSTWFitMale) #Check if the confidence intervals include 1 to assess allometric trajectory
    confint(jmPPTCSTWFitMale)
   
    
#Put it all together in a table
    jPPTCS <- matrix(, nrow = 3, ncol = 3)
    colnames(jPPTCS) <- c("p", "lowerCI", "upperCI")
    rownames(jPPTCS) <- c("L1", "TL", "TW")
    
    jPPTCS[[1, 1]] <- round(summary(jmPPTL1CSFitMale)$tTable[2, 4], 3)
    jPPTCS[[1, 2]] <- round(confint(jmPPTL1CSFitMale)[2, 1], 3)
    jPPTCS[[1, 3]] <- round(confint(jmPPTL1CSFitMale)[2, 2], 3)
    
    jPPTCS[[2, 1]] <- round(summary(jmPPTCSTLFitMale)$tTable[2, 4], 3)
    jPPTCS[[2, 2]] <- round(confint(jmPPTCSTLFitMale)[2, 1], 3)
    jPPTCS[[2, 3]] <- round(confint(jmPPTCSTLFitMale)[2, 2], 3)
    
    jPPTCS[[3, 1]] <- round(summary(jmPPTCSTWFitMale)$tTable[2, 4], 3)
    jPPTCS[[3, 2]] <- round(confint(jmPPTCSTWFitMale)[2, 1], 3)
    jPPTCS[[3, 3]] <- round(confint(jmPPTCSTWFitMale)[2, 2], 3)
    
    
    #Do a Procrustes ANOVA of shape and L1
    jmalePPTFit <- procD.lm(shape ~ L1,
                            data = jmalePPTGDF,
                            iter = 1000)
    jPPTMaleFitSum <- summary(jmalePPTFit) 
    
    #And with CS of the claspers as in the males
    jmalePPTFit2 <- procD.lm(shape ~ cs,
                             data = jmalePPTGDF,
                             iter = 1000)
    jPPTMaleFitSum2 <- summary(jmalePPTFit2) 
    
    #Do a Procrustes ANOVA of shape and TL
    jmalePPTFit3 <- procD.lm(shape ~ TL,
                             data = jmalePPTGDF,
                             iter = 1000)
    jPPTMaleFitSum3 <- summary(jmalePPTFit3) 
    
    #Do a Procrustes ANOVA of shape and TW
    jmalePPTFit4 <- procD.lm(shape ~ TW,
                             data = jmalePPTGDF,
                             iter = 1000)
    jPPTMaleFitSum4 <- summary(jmalePPTFit4) 

#Move relevant info into a table 
    jPPTShapeMat <- matrix(, nrow = 4, ncol = 2)
    colnames(jPPTShapeMat) <- c("R^2", "p")
    rownames(jPPTShapeMat) <- c("L1", "CS", "TL", "TW")
    
    jPPTShapeMat[[1, 1]] <- round(jPPTMaleFitSum$table$Rsq[[1]], 3)
    jPPTShapeMat[[1, 2]] <- round(jPPTMaleFitSum$table$`Pr(>F)`[[1]], 3)
    jPPTShapeMat[[2, 1]] <- round(jPPTMaleFitSum2$table$Rsq[[1]], 3)
    jPPTShapeMat[[2, 2]] <- round(jPPTMaleFitSum2$table$`Pr(>F)`[[1]], 3)
    jPPTShapeMat[[3, 1]] <- round(jPPTMaleFitSum3$table$Rsq[[1]], 3)
    jPPTShapeMat[[3, 2]] <- round(jPPTMaleFitSum3$table$`Pr(>F)`[[1]], 3)
    jPPTShapeMat[[4, 1]] <- round(jPPTMaleFitSum4$table$Rsq[[1]], 3)
    jPPTShapeMat[[4, 2]] <- round(jPPTMaleFitSum4$table$`Pr(>F)`[[1]], 3)
    
    ####For all PPT####
    
    allmalePPTShape <- read.morphologika("allpptmorphologikaINT_unscaled.txt") #Load male shape data
    allmalePPTGPA <- gpagen(allmalePPTShape, 
                           ProcD = FALSE) #Run Generalized Procrustes Analysis using bending energy criterion
    
    #Change to clasperMaleData
    allmaleData <- read.csv("allINT.csv") #Load in data with specimenID, L1 (total length), TL, TW, AG without HYCO010M
    
    #male PCA to examine shape trends in morphospace
    allmalePPTPCA <- gm.prcomp(allmalePPTGPA$coords)
    allmalePPTPCASum <- summary(allmalePPTPCA) #Gives you breakdown of variance and eigenvalues
    
    #maleClasperPCASum$PC.summary and maleClasperPCA$x are exported as supplemental file XX
    
    #Set a color palette for the PCA based on L1 where larger males are red and smaller males are blue
    pal = colorRampPalette(c("blue", "red"))
    allmalePPTGPA$col <- pal(10)[as.numeric(cut(log10(allmaleData$L1), breaks = 10))] #This breaks into 10 colors
    
    #FIGURE XX
    allmalePPTPCAPlot <- plot(allmalePPTPCA,
                             xlab = paste("Principal Component 1 ", "(", sep = "", 
                                          paste(round(allmalePPTPCASum$PC.summary$Comp1[2]*100, digits = 2), "%", ")", sep = "")),
                             ylab = paste("Principal Component 2 ", "(", sep = "", 
                                          paste(round(allmalePPTPCASum$PC.summary$Comp2[2]*100, digits = 2), "%", ")", sep = "")),
                             pch = 21,
                             col = "black",
                             bg = allmalePPTGPA$col,
                             cex = 4,
                             cex.lab = 1.2)
    text(allmalePPTPCAPlot$PC.points[ , 2] ~ allmalePPTPCAPlot$PC.points[ , 1], 
         labels = allmaleData$Order, cex= 1, col = c("white"))
    
    #generate a color gradient for an inset in the plot
    plot(rep(1,100),col=(pal(100)), pch=19,cex=10)
    
    #Make Shape Extremes with Point Clouds using Morpho package
    procallPPTMale <- procSym(allmalePPTShape)
    plot(procallPPTMale$PCscores[,2] ~ procallPPTMale$PCscores[,1])
    
    #Plot at ends of PC1 axes, negative is purple and positive is green, Magnified by 2 to accentuate trends
    palette(c("white", "green", "purple"))
    posPC1allPPTmale <- restoreShapes(2*sd(procallPPTMale$PCscores[,1]), 
                                     procallPPTMale$PCs[,1], 
                                     procallPPTMale$mshape)
    negPC1allPPTmale <- restoreShapes(-2*sd(procallPPTMale$PCscores[,1]), 
                                     procallPPTMale$PCs[,1], 
                                     procallPPTMale$mshape)
    deformGrid3d(posPC1allPPTmale, 
                 negPC1allPPTmale, 
                 ngrid = 0, 
                 lines = FALSE) #Type no
    
    #and again for PC2
    posPC2allPPTmale <- restoreShapes(2*sd(procallPPTMale$PCscores[,2]), 
                                     procallPPTMale$PCs[,2], 
                                     procallPPTMale$mshape)
    negPC2allPPTmale <- restoreShapes(-2*sd(procallPPTMale$PCscores[,2]), 
                                     procallPPTMale$PCs[,2], 
                                     procallPPTMale$mshape)
    deformGrid3d(posPC2allPPTmale, 
                 negPC2allPPTmale, 
                 ngrid = 0, 
                 lines = FALSE) #Type no
    
    #Make a geomorph dataframe for using functions in geomorph package 
    allmalePPTGDF <-geomorph.data.frame(shape = allmalePPTGPA$coords,
                                       cs = log10(allmalePPTGPA$Csize),
                                       L1 = log10(allmaleData$L1),
                                       TL = log10(allmaleData$TL),
                                       TW = log10(allmaleData$TW)) #Note that Csize is log10 transformed for downstream analyses
    
    #Examine relationship between L1 and clasper centroid size (log10 transformed)
    allmalePPTDFMale <- data.frame(log10(allmaleData$L1),
                                  as.numeric(log10(allmalePPTGPA$Csize)),
                                  log10(allmaleData$TL),
                                  log10(allmaleData$TW))
    
    colnames(allmalePPTDFMale) <- c("MaleL1", "MaleCS", "TL", "TW") #Name columns
    
    #L1CSFit a gls model and plot the relationship
    allPPTL1CSFitMale <- gls(MaleCS ~ MaleL1, 
                            data = allmalePPTDFMale)
    
    plot(MaleCS ~ MaleL1, 
         data = allmalePPTDFMale,
         pch = 21,
         col = "black",
         bg = allmalePPTGPA$col,
         cex = 4,
         cex.lab = 1.2,
         xlim = c(1.46, 1.65), #NEED TO CHANGE ACROSS VARIABLES
         ylim = c(1.3,3), #NEED TO CHANGE ACROSS VARIABLES
         xlab = "log10 (L1)",
         ylab = "log10 (Centroid Size)")
    abline(a = coef(allPPTL1CSFitMale)[1],
           b = coef(allPPTL1CSFitMale)[2],
           lw = 3)
    
    summary(allPPTL1CSFitMale) #Check if the confidence intervals include 1 to assess allometric trajectory
    confint(allPPTL1CSFitMale)
    
    #CSTLFit a gls model and plot the relationship
    allPPTCSTLFitMale <- gls(MaleCS ~ TL, 
                            data = allmalePPTDFMale)
    
    plot(MaleCS ~ TL, 
         data = allmalePPTDFMale,
         pch = 21,
         col = "black",
         bg = allmalePPTGPA$col,
         cex = 4,
         cex.lab = 1.2,
         xlim = c(-0.4,0.6), #NEED TO CHANGE ACROSS VARIABLES
         ylim = c(2.2,2.7), #NEED TO CHANGE ACROSS VARIABLES
         xlab = "log10(TL)",
         ylab = "log10 (Centroid Size)")
    abline(a = coef(allPPTCSTLFitMale)[1],
           b = coef(allPPTCSTLFitMale)[2],
           lw = 3)
    
    summary(allPPTCSTLFitMale) #Check if the confidence intervals include 1 to assess allometric trajectory
    confint(allPPTCSTLFitMale)
    
    #CSTWFit a gls model and plot the relationship
    allPPTCSTWFitMale <- gls(MaleCS ~ TW, 
                            data = allmalePPTDFMale)
    
    plot(MaleCS ~ TW, 
         data = allmalePPTDFMale,
         pch = 21,
         col = "black",
         bg = allmalePPTGPA$col,
         cex = 4,
         cex.lab = 1.2,
         xlim = c(-0.8,0.4), #NEED TO CHANGE ACROSS VARIABLES
         ylim = c(2.2,2.8), #NEED TO CHANGE ACROSS VARIABLES
         xlab = "log10(TW)",
         ylab = "log10 (Centroid Size)")
    abline(a = coef(allPPTCSTWFitMale)[1],
           b = coef(allPPTCSTWFitMale)[2],
           lw = 3)
    
    summary(allPPTCSTWFitMale) #Check if the confidence intervals include 1 to assess allometric trajectory
    confint(allPPTCSTWFitMale)
    
    #L1TLFit a gls model and plot the relationship
    allPPTL1TLFitMale <- gls(TL ~ MaleL1, 
                              data = allmalePPTDFMale)
    
    plot(TL ~ MaleL1, 
         data = allmalePPTDFMale,
         pch = 21,
         col = "black",
         bg = allmalePPTGPA$col,
         cex = 4,
         cex.lab = 1.2,
         xlim = c(1.45,1.65), #NEED TO CHANGE ACROSS VARIABLES
         ylim = c(-0.5,0.6), #NEED TO CHANGE ACROSS VARIABLES
         xlab = "log10 (L1)",
         ylab = "log10 (TL)")
    abline(a = coef(allPPTL1TLFitMale)[1],
           b = coef(allPPTL1TLFitMale)[2],
           lw = 3)
    
    summary(allPPTL1TLFitMale) #Check if the confidence intervals include 1 to assess allometric trajectory
    confint(allPPTL1TLFitMale)
    
    #L1TWFit a gls model and plot the relationship
    allPPTL1TWFitMale <- gls(TW ~ MaleL1, 
                              data = allmalePPTDFMale)
    
    plot(TW ~ MaleL1, 
         data = allmalePPTDFMale,
         pch = 21,
         col = "black",
         bg = allmalePPTGPA$col,
         cex = 4,
         cex.lab = 1.2,
         xlim = c(1.45,1.65), #NEED TO CHANGE ACROSS VARIABLES
         ylim = c(-1,0.5), #NEED TO CHANGE ACROSS VARIABLES
         xlab = "log10 (L1)",
         ylab = "log10 (TW)")
    abline(a = coef(allPPTL1TWFitMale)[1],
           b = coef(allPPTL1TWFitMale)[2],
           lw = 3)
    
    summary(allPPTL1TWFitMale) #Check if the confidence intervals include 1 to assess allometric trajectory
    confint(allPPTL1TWFitMale)
    
    #Put it all together in a table
    allPPTCS <- matrix(, nrow = 3, ncol = 3)
    colnames(allPPTCS) <- c("p", "lowerCI", "upperCI")
    rownames(allPPTCS) <- c("L1", "TL", "TW")
    
    allPPTCS[[1, 1]] <- round(summary(allPPTL1CSFitMale)$tTable[2, 4], 3)
    allPPTCS[[1, 2]] <- round(confint(allPPTL1CSFitMale)[2, 1], 3)
    allPPTCS[[1, 3]] <- round(confint(allPPTL1CSFitMale)[2, 2], 3)
    
    alLPPTCS[[2, 1]] <- round(summary(allPPTCSTLFitMale)$tTable[2, 4], 3)
    allPPTCS[[2, 2]] <- round(confint(allPPTCSTLFitMale)[2, 1], 3)
    allPPTCS[[2, 3]] <- round(confint(allPPTCSTLFitMale)[2, 2], 3)
    
    allPPTCS[[3, 1]] <- round(summary(allPPTCSTWFitMale)$tTable[2, 4], 3)
    allPPTCS[[3, 2]] <- round(confint(allPPTCSTWFitMale)[2, 1], 3)
    allPPTCS[[3, 3]] <- round(confint(allPPTCSTWFitMale)[2, 2], 3)
    
    
    #Do a Procrustes ANOVA of shape and L1
    allmalePPTFit <- procD.lm(shape ~ L1,
                             data = allmalePPTGDF,
                             iter = 1000)
    allPPTMaleFitSum <- summary(allmalePPTFit) 
    
    #And with CS of the claspers as in the males
    allmalePPTFit2 <- procD.lm(shape ~ cs,
                              data = allmalePPTGDF,
                              iter = 1000)
    allPPTMaleFitSum2 <- summary(allmalePPTFit2) 
    
    #Do a Procrustes ANOVA of shape and TL
    allmalePPTFit3 <- procD.lm(shape ~ TL,
                              data = allmalePPTGDF,
                              iter = 1000)
    allPPTMaleFitSum3 <- summary(allmalePPTFit3) 
    
    #Do a Procrustes ANOVA of shape and TW
    allmalePPTFit4 <- procD.lm(shape ~ TW,
                              data = allmalePPTGDF,
                              iter = 1000)
    allPPTMaleFitSum4 <- summary(allmalePPTFit4) 
    
    #Move relevant info into a table 
    allPPTShapeMat <- matrix(, nrow = 4, ncol = 2)
    colnames(allPPTShapeMat) <- c("R^2", "p")
    rownames(allPPTShapeMat) <- c("L1", "CS", "TL", "TW")
    
    allPPTShapeMat[[1, 1]] <- round(allPPTMaleFitSum$table$Rsq[[1]], 3)
    allPPTShapeMat[[1, 2]] <- round(allPPTMaleFitSum$table$`Pr(>F)`[[1]], 3)
    allPPTShapeMat[[2, 1]] <- round(allPPTMaleFitSum2$table$Rsq[[1]], 3)
    allPPTShapeMat[[2, 2]] <- round(allPPTMaleFitSum2$table$`Pr(>F)`[[1]], 3)
    allPPTShapeMat[[3, 1]] <- round(allPPTMaleFitSum3$table$Rsq[[1]], 3)
    allPPTShapeMat[[3, 2]] <- round(allPPTMaleFitSum3$table$`Pr(>F)`[[1]], 3)
    allPPTShapeMat[[4, 1]] <- round(allPPTMaleFitSum4$table$Rsq[[1]], 3)
    allPPTShapeMat[[4, 2]] <- round(allPPTMaleFitSum4$table$`Pr(>F)`[[1]], 3)
    
    
    allPPTShapeMat

    
    
#Include age as a factor
    agePPTFit <- procD.lm(shape ~ L1 + Age,
                          data = allPPTGDF,
                          iter = 1000)
    summary(agePPTFit)
    
    agePPTFit2 <- procD.lm(shape ~ Age,
                           data = allPPTGDF,
                           iter = 1000)
    summary(agePPTFit2)
    
    agePPTFit3 <- procD.lm(shape ~ cs + Age,
                           data = allPPTGDF,
                           iter = 1000)
    summary(agePPTFit3)
    
    agePPTFit4 <- procD.lm(shape ~ cs,
                           data = allPPTGDF,
                           iter = 1000)
    summary(agePPTFit4)
    
    
    
    
####Adult Frontal Tenaculum####    
    maleFTShape <- read.morphologika("ft2morphologikaNo11_unscaled.txt") #Load male shape data
    maleFTGPA <- gpagen(maleFTShape, 
                        ProcD = FALSE) #Run GPA using bending energy criterion
    
    FTmaleData <- read.csv("ftmratfishNo.11.csv") #Load in data with specimenID, L2cm (length 2 cm)
    
    #male PCA to examine shape trends in morphospace
    maleFTPCA <- gm.prcomp(maleFTGPA$coords)
    maleFTPCASum <- summary(maleFTPCA) #Gives you breakdown of variance and eigenvalues
    
    #maleFTPCASum$PC.summary and maleFTPCA$x are exported as supplemental file XX
    
    #Set a color palette for the PCA based on L1 where larger males are red and smaller males are blue
    pal = colorRampPalette(c("blue", "red"))
    maleFTGPA$col <- pal(10)[as.numeric(cut(log10(FTmaleData$L1), breaks = 10))]
    
    #FIGURE XX
    maleFTPCAPlot <- plot(maleFTPCA,
                          xlab = paste("Principal Component 1 ", "(", sep = "", 
                                       paste(round(maleFTPCASum$PC.summary$Comp1[2]*100, digits = 2), "%", ")", sep = "")),
                          ylab = paste("Principal Component 2 ", "(", sep = "", 
                                       paste(round(maleFTPCASum$PC.summary$Comp2[2]*100, digits = 2), "%",")", sep = "")),
                          pch = 21,
                          col = "black",
                          bg = maleFTGPA$col,
                          cex = 4,
                          cex.lab = 1.2)
    text(maleFTPCAPlot$PC.points[ , 2] ~ maleFTPCAPlot$PC.points[ , 1], 
         labels = FTmaleData$Order, cex= 1, col = c("white"))
    
    #Make Shape Extremes with Point Clouds using Morpho package
    procFTMale <- procSym(maleFTShape)
    plot(procFTMale$PCscores[,2] ~ procFTMale$PCscores[,1])
    
    #Plot at ends of PC1 axes, negative is purple and positive is green, Magnified by 2 to accentuate trends
    palette(c("white", "green", "purple"))
    posPC1FTmale <- restoreShapes(2*sd(procFTMale$PCscores[,1]), 
                                  procFTMale$PCs[,1], 
                                  procFTMale$mshape)
    negPC1FTmale <- restoreShapes(-2*sd(procFTMale$PCscores[,1]), 
                                  procFTMale$PCs[,1], 
                                  procFTMale$mshape)
    deformGrid3d(posPC1FTmale, 
                 negPC1FTmale, 
                 ngrid = 0, 
                 lines = FALSE) #type no
    
    #and again for PC2
    posPC2FTmale <- restoreShapes(2*sd(procFTMale$PCscores[,2]), 
                                  procFTMale$PCs[,2], 
                                  procFTMale$mshape)
    negPC2FTmale <- restoreShapes(-2*sd(procFTMale$PCscores[,2]), 
                                  procFTMale$PCs[,2], 
                                  procFTMale$mshape)
    deformGrid3d(posPC2FTmale,
                 negPC2FTmale, 
                 ngrid = 0, 
                 lines = FALSE) #type no
    
    #Make a geomorph dataframe for using functions in geomorph package 
    maleFTGDF <-geomorph.data.frame(shape = maleFTGPA$coords,
                                    cs = log10(maleFTGPA$Csize),
                                    L1 = log10(FTmaleData$L1),
                                    TL = log10(FTmaleData$TL),
                                    TW = log10(FTmaleData$TW))
    
    #Make dataframe to compare L1 and FT Csize
    FTTDFMale <- data.frame(log10(FTmaleData$L1),
                            as.numeric(log10(maleFTGPA$Csize)), 
                            log10(FTmaleData$TL),
                            log10(FTmaleData$TW))
    
    colnames(FTTDFMale) <- c("MaleL1", "MaleCS", "TL", "TW")
    
    
    #Fit GLS of CS and L1
    FTL1CSFitMale <- gls(MaleCS ~ MaleL1, 
                       data = FTTDFMale)
    
    plot(MaleCS ~ MaleL1, 
         data = FTTDFMale,
         pch = 21,
         col = "black",
         bg = maleFTGPA$col,
         cex = 4,
         cex.lab = 1.2,
         xlim = c(1.55,1.65), #NEED TO CHANGE ACROSS VARIABLES
         ylim = c(1.63,1.8), #NEED TO CHANGE ACROSS VARIABLES
         xlab = "log10 (L1)",
         ylab = "log10 (Centroid Size)")
    abline(a = coef(FTL1CSFitMale)[1],
           b = coef(FTL1CSFitMale)[2],
           lw = 3)
    
    summary(FTL1CSFitMale) #Check if the confidence intervals include 1 to assess allometric trajectory
    confint(FTL1CSFitMale)
    
    #CSTLFit a gls model and plot the relationship
    FTCSTLFitMale <- gls(MaleCS ~ TL, 
                       data = FTTDFMale)
    
    plot(MaleCS ~ TL, 
         data = FTTDFMale,
         pch = 21,
         col = "black",
         bg = maleFTGPA$col,
         cex = 4,
         cex.lab = 1.2,
         xlim = c(0.25, 0.53), #NEED TO CHANGE ACROSS VARIABLES
         ylim = c(1.64, 1.8), #NEED TO CHANGE ACROSS VARIABLES
         xlab = "log10 (TL)",
         ylab = "log10 (Centroid Size)") 
    abline(a = coef(FTCSTLFitMale)[1],
           b = coef(FTCSTLFitMale)[2],
           lw = 3)
    
    summary(FTCSTLFitMale) #Check if the confidence intervals include 1 to assess allometric trajectory
    confint(FTCSTLFitMale)
    
    #CSTWFit a gls model and plot the relationship
    FTCSTWFitMale <- gls(MaleCS ~ TW, 
                       data = FTTDFMale)
    
    plot(MaleCS ~ TW, 
         data = FTTDFMale,
         pch = 21,
         col = "black",
         bg = maleFTGPA$col,
         cex = 4,
         cex.lab = 1.2,
         xlim = c(0.05, 0.35), #NEED TO CHANGE ACROSS VARIABLES
         ylim = c(1.65, 1.8), #NEED TO CHANGE ACROSS VARIABLES
         xlab = "log10 (TW)",
         ylab = "log10 (Centroid Size)")
    abline(a = coef(FTCSTWFitMale)[1],
           b = coef(FTCSTWFitMale)[2],
           lw = 3)
    
    summary(FTCSTWFitMale) #Check if the confidence intervals include 1 to assess allometric trajectory
    confint(FTCSTWFitMale)
    
#Put it all together in a table
    FTCSMat <- matrix(, nrow = 3, ncol = 3)
    colnames(FTCSMat) <- c("p", "lowerCI", "upperCI")
    rownames(FTCSMat) <- c("L1", "TL", "TW")
    
    FTCSMat[[1, 1]] <- round(summary(FTL1CSFitMale)$tTable[2, 4], 3)
    FTCSMat[[1, 2]] <- round(confint(FTL1CSFitMale)[2, 1], 3)
    FTCSMat[[1, 3]] <- round(confint(FTL1CSFitMale)[2, 2], 3)
    
    FTCSMat[[2, 1]] <- round(summary(FTCSTLFitMale)$tTable[2, 4], 3)
    FTCSMat[[2, 2]] <- round(confint(FTCSTLFitMale)[2, 1], 3)
    FTCSMat[[2, 3]] <- round(confint(FTCSTLFitMale)[2, 2], 3)
    
    FTCSMat[[3, 1]] <- round(summary(FTCSTWFitMale)$tTable[2, 4], 3)
    FTCSMat[[3, 2]] <- round(confint(FTCSTWFitMale)[2, 1], 3)
    FTCSMat[[3, 3]] <- round(confint(FTCSTWFitMale)[2, 2], 3)
    
#Procrustes ANOVA of shape ~ L1, CS, TL, TW
    #Do a Procrustes ANOVA of FT shape and L2
    FTmaleFit <- procD.lm(shape ~ L1,
                          data = maleFTGDF,
                          iter = 1000)
    FTmaleFitSum <- summary(FTmaleFit) 
    
    #And with CS of the FT as in the males
    FTmaleFit2 <- procD.lm(shape ~ cs,
                           data = maleFTGDF,
                           iter = 1000)
    FTmaleFitSum2 <- summary(FTmaleFit2) 
    
    #And with TL of the male
    FTmaleFit3 <- procD.lm(shape ~ TL,
                           data = maleFTGDF,
                           iter=1000)
    FTmaleFitSum3 <- summary(FTmaleFit3)
  
    #And with TW of the male 
    FTmaleFit4 <- procD.lm(shape ~ TW,
                           data = maleFTGDF,
                           iter=1000)
    FTmaleFitSum4 <- summary(FTmaleFit4)

#Move relevant info into a table 
    FTShapeMat <- matrix(, nrow = 4, ncol = 2)
    colnames(FTShapeMat) <- c("R^2", "p")
    rownames(FTShapeMat) <- c("L1", "CS", "TL", "TW")
    
    FTShapeMat[[1, 1]] <- round(FTmaleFitSum$table$Rsq[[1]], 3)
    FTShapeMat[[1, 2]] <- round(FTmaleFitSum$table$`Pr(>F)`[[1]], 3)
    FTShapeMat[[2, 1]] <- round(FTmaleFitSum2$table$Rsq[[1]], 3)
    FTShapeMat[[2, 2]] <- round(FTmaleFitSum2$table$`Pr(>F)`[[1]], 3)
    FTShapeMat[[3, 1]] <- round(FTmaleFitSum3$table$Rsq[[1]], 3)
    FTShapeMat[[3, 2]] <- round(FTmaleFitSum3$table$`Pr(>F)`[[1]], 3)
    FTShapeMat[[4, 1]] <- round(FTmaleFitSum4$table$Rsq[[1]], 3)
    FTShapeMat[[4, 2]] <- round(FTmaleFitSum4$table$`Pr(>F)`[[1]], 3)
    
####Juvenile Frontal Tenaculum####
    jmaleFTShape <- read.morphologika("juvft_unscaled.txt") #Load male shape data
    jmaleFTGPA <- gpagen(jmaleFTShape, 
                          ProcD = FALSE) #Run Generalized Procrustes Analysis using bending energy criterion
    
    #Change to clasperMaleData
    jmaleData <- read.csv("jmratfish.csv") #Load in data with specimenID, L1 (total length), TL, TW, AG without HYCO010M
    
    #male PCA to examine shape trends in morphospace
    jmaleFTPCA <- gm.prcomp(jmaleFTGPA$coords)
    jmaleFTPCASum <- summary(jmaleFTPCA) #Gives you breakdown of variance and eigenvalues
    
    #maleClasperPCASum$PC.summary and maleClasperPCA$x are exported as supplemental file XX
    
    #Set a color palette for the PCA based on L1 where larger males are red and smaller males are blue
    pal = colorRampPalette(c("blue", "red"))
    jmaleFTGPA$col <- pal(10)[as.numeric(cut(log10(jmaleData$L1), breaks = 10))] #This breaks into 10 colors
    
    #FIGURE XX
    jmaleFTPCAPlot <- plot(jmaleFTPCA,
                            xlab = paste("Principal Component 1 ", "(", sep = "", 
                                         paste(round(jmaleFTPCASum$PC.summary$Comp1[2]*100, digits = 2), "%", ")", sep = "")),
                            ylab = paste("Principal Component 2 ", "(", sep = "", 
                                         paste(round(jmaleFTPCASum$PC.summary$Comp2[2]*100, digits = 2), "%", ")", sep = "")),
                            pch = 21,
                            col = "black",
                            bg = jmaleFTGPA$col,
                            cex = 4,
                            cex.lab = 1.2)
    text(jmaleFTPCAPlot$PC.points[ , 2] ~ jmaleFTPCAPlot$PC.points[ , 1], 
         labels = jmaleData$Order, cex= 1, col = c("white"))
    
    #generate a color gradient for an inset in the plot
    plot(rep(1,100),col=(pal(100)), pch=19,cex=10)
    
    #Make Shape Extremes with Point Clouds using Morpho package
    procjFTMale <- procSym(jmaleFTShape)
    plot(procjFTMale$PCscores[,2] ~ procjFTMale$PCscores[,1])
    
    #Plot at ends of PC1 axes, negative is purple and positive is green, Magnified by 2 to accentuate trends
    palette(c("white", "green", "purple"))
    posPC1jFTmale <- restoreShapes(2*sd(procjFTMale$PCscores[,1]), 
                                    procjFTMale$PCs[,1], 
                                    procjFTMale$mshape)
    negPC1jFTmale <- restoreShapes(-2*sd(procjFTMale$PCscores[,1]), 
                                    procjFTMale$PCs[,1], 
                                    procjFTMale$mshape)
    deformGrid3d(posPC1jFTmale, 
                 negPC1jFTmale, 
                 ngrid = 0, 
                 lines = FALSE) #Type no
    
    #and again for PC2
    posPC2jFTmale <- restoreShapes(2*sd(procjFTMale$PCscores[,2]), 
                                    procjFTMale$PCs[,2], 
                                    procjFTMale$mshape)
    negPC2jFTmale <- restoreShapes(-2*sd(procjFTMale$PCscores[,2]), 
                                    procjFTMale$PCs[,2], 
                                    procjFTMale$mshape)
    deformGrid3d(posPC2jFTmale, 
                 negPC2jFTmale, 
                 ngrid = 0, 
                 lines = FALSE) #Type no
    
    #Make a geomorph dataframe for using functions in geomorph package 
    jmaleFTGDF <-geomorph.data.frame(shape = jmaleFTGPA$coords,
                                      cs = log10(jmaleFTGPA$Csize),
                                      L1 = log10(jmaleData$L1),
                                      TL = log10(jmaleData$TL),
                                      TW = log10(jmaleData$TW)) #Note that Csize is log10 transformed for downstream analyses
    
    #Examine relationship between L1 and clasper centroid size (log10 transformed)
    jmaleFTDFMale <- data.frame(log10(jmaleData$L1),
                                 as.numeric(log10(jmaleFTGPA$Csize)),
                                 log10(jmaleData$TL),
                                 log10(jmaleData$TW))
    
    colnames(jmaleFTDFMale) <- c("MaleL1", "MaleCS", "TL", "TW") #Name columns
    
    #L1CSFit a gls model and plot the relationship
    jmFTL1CSFitMale <- gls(MaleCS ~ MaleL1, 
                            data = jmaleFTDFMale)
    
    plot(MaleCS ~ MaleL1, 
         data = jmaleFTDFMale,
         pch = 21,
         col = "black",
         bg = jmaleFTGPA$col,
         cex = 4,
         cex.lab = 1.2,
         xlim = c(1.45, 1.58), #NEED TO CHANGE ACROSS VARIABLES
         ylim = c(-1.1,-0.6),#NEED TO CHANGE ACROSS VARIABLES
         xlab = "log10 (L1)",
         ylab = "log10 (Centroid Size)")
    abline(a = coef(jmFTL1CSFitMale)[1],
           b = coef(jmFTL1CSFitMale)[2],
           lw = 3)
    
    summary(jmFTL1CSFitMale) #Check if the confidence intervals include 1 to assess allometric trajectory
    confint(jmFTL1CSFitMale)
    
    #CSTLFit a gls model and plot the relationship
    jmFTCSTLFitMale <- gls(MaleCS ~ TL, 
                            data = jmaleFTDFMale)
    
    plot(MaleCS ~ TL, 
         data = jmaleFTDFMale,
         pch = 21,
         col = "black",
         bg = jmaleFTGPA$col,
         cex = 4,
         cex.lab = 1.2,
         xlim = c(-0.4,0.1), #NEED TO CHANGE ACROSS VARIABLES
         ylim = c(-1.1,-0.6), #NEED TO CHANGE ACROSS VARIABLES
         xlab = "log10(TL)",
         ylab = "log10 (Centroid Size)")
    abline(a = coef(jmFTCSTLFitMale)[1],
           b = coef(jmFTCSTLFitMale)[2],
           lw = 3)
    
    summary(jmFTCSTLFitMale) #Check if the confidence intervals include 1 to assess allometric trajectory
    confint(jmFTCSTLFitMale)
    
    #CSTWFit a gls model and plot the relationship
    jmFTCSTWFitMale <- gls(MaleCS ~ TW, 
                            data = jmaleFTDFMale)
    
    plot(MaleCS ~ TW, 
         data = jmaleFTDFMale,
         pch = 21,
         col = "black",
         bg = jmaleFTGPA$col,
         cex = 4,
         cex.lab = 1.2,
         xlim = c(-0.8,0), #NEED TO CHANGE ACROSS VARIABLES
         ylim = c(-1.2,-0.6), #NEED TO CHANGE ACROSS VARIABLES
         xlab = "log10(TW)",
         ylab = "log10 (Centroid Size)")
    abline(a = coef(jmFTCSTWFitMale)[1],
           b = coef(jmFTCSTWFitMale)[2],
           lw = 3)
    
    summary(jmFTCSTWFitMale) #Check if the confidence intervals include 1 to assess allometric trajectory
    confint(jmFTCSTWFitMale)
    
    
    #Put it all together in a table
    jFTCS <- matrix(, nrow = 3, ncol = 3)
    colnames(jFTCS) <- c("p", "lowerCI", "upperCI")
    rownames(jFTCS) <- c("L1", "TL", "TW")
    
    jPPTCS[[1, 1]] <- round(summary(jmFTL1CSFitMale)$tTable[2, 4], 3)
    jPPTCS[[1, 2]] <- round(confint(jmFTL1CSFitMale)[2, 1], 3)
    jPPTCS[[1, 3]] <- round(confint(jmFTL1CSFitMale)[2, 2], 3)
    
    jPPTCS[[2, 1]] <- round(summary(jmFTCSTLFitMale)$tTable[2, 4], 3)
    jPPTCS[[2, 2]] <- round(confint(jmFTCSTLFitMale)[2, 1], 3)
    jPPTCS[[2, 3]] <- round(confint(jmFTCSTLFitMale)[2, 2], 3)
    
    jPPTCS[[3, 1]] <- round(summary(jmFTCSTWFitMale)$tTable[2, 4], 3)
    jPPTCS[[3, 2]] <- round(confint(jmFTCSTWFitMale)[2, 1], 3)
    jPPTCS[[3, 3]] <- round(confint(jmFTCSTWFitMale)[2, 2], 3)
    
    
    #Do a Procrustes ANOVA of shape and L1
    jmaleFTFit1 <- procD.lm(shape ~ L1,
                            data = jmaleFTGDF,
                            iter = 1000)
    jFTMaleFitSum <- summary(jmaleFTFit1) 
    
    #And with CS of the claspers as in the males
    jmaleFTFit2 <- procD.lm(shape ~ cs,
                             data = jmaleFTGDF,
                             iter = 1000)
    jFTMaleFitSum2 <- summary(jmaleFTFit2) 
    
    #Do a Procrustes ANOVA of shape and TL
    jmaleFTFit3 <- procD.lm(shape ~ TL,
                             data = jmaleFTGDF,
                             iter = 1000)
    jFTMaleFitSum3 <- summary(jmaleFTFit3) 
    
    #Do a Procrustes ANOVA of shape and TW
    jmaleFTFit4 <- procD.lm(shape ~ TW,
                             data = jmaleFTGDF,
                             iter = 1000)
    jFTMaleFitSum4 <- summary(jmaleFTFit4) 
    
    #Move relevant info into a table 
    jFTShapeMat <- matrix(, nrow = 4, ncol = 2)
    colnames(jFTShapeMat) <- c("R^2", "p")
    rownames(jFTShapeMat) <- c("L1", "CS", "TL", "TW")
    
    jFTShapeMat[[1, 1]] <- round(jFTMaleFitSum$table$Rsq[[1]], 3)
    jFTShapeMat[[1, 2]] <- round(jFTMaleFitSum$table$`Pr(>F)`[[1]], 3)
    jFTShapeMat[[2, 1]] <- round(jFTMaleFitSum2$table$Rsq[[1]], 3)
    jFTShapeMat[[2, 2]] <- round(jFTMaleFitSum2$table$`Pr(>F)`[[1]], 3)
    jFTShapeMat[[3, 1]] <- round(jFTMaleFitSum3$table$Rsq[[1]], 3)
    jFTShapeMat[[3, 2]] <- round(jFTMaleFitSum3$table$`Pr(>F)`[[1]], 3)
    jFTShapeMat[[4, 1]] <- round(jFTMaleFitSum4$table$Rsq[[1]], 3)
    jFTShapeMat[[4, 2]] <- round(jFTMaleFitSum4$table$`Pr(>F)`[[1]], 3)
    
    

####General allometric trajectories####
#(using clasper datasheet, but does not relate to claspers)
    
    #For Juveniles look at Testes length versus L1
    #Fit a gls model and plot the relationship
    juvTLL1FitMale <- gls(TL ~ MaleL1, 
                          data = jclasperDFMale)
    
    plot(TL ~ MaleL1, 
         data = jclasperDFMale,
         pch = 21,
         col = "black",
         bg = jclasperGPA$col,
         cex = 4,
         cex.lab = 1.2,
         xlim = c(1.46, 1.58), #NEED TO CHANGE ACROSS VARIABLES
         ylim = c(0.3, 1.2), #NEED TO CHANGE ACROSS VARIABLES
         xlab = "log10 (L1)",
         ylab = "Testes Length")
    abline(a = coef(juvTLL1FitMale)[1],
           b = coef(juvTLL1FitMale)[2],
           lw = 3)
    
    summary(juvTLL1FitMale) #Check if the confidence intervals include 1 to assess allometric trajectory
    confint(juvTLL1FitMale)
    
    #For Juveniles look at Testes width versus L1
    #Fit a gls model and plot the relationship
    juvTWL1FitMale <- gls(TW ~ MaleL1, 
                          data = jclasperDFMale)
    
    plot(TW ~ MaleL1, 
         data = jclasperDFMale,
         pch = 21,
         col = "black",
         bg = jclasperGPA$col,
         cex = 4,
         cex.lab = 1.2,
         xlim = c(1.46, 1.58), #NEED TO CHANGE ACROSS VARIABLES
         ylim = c(0, 0.9), #NEED TO CHANGE ACROSS VARIABLES
         xlab = "log10 (L1)",
         ylab = "Testes Width")
    abline(a = coef(juvTWL1FitMale)[1],
           b = coef(juvTWL1FitMale)[2],
           lw = 3)
    
    summary(juvTWL1FitMale) #Check if the confidence intervals include 1 to assess allometric trajectory
    confint(juvTWL1FitMale)
    
    #For Adults look at Testes length versus L1
    #Fit a gls model and plot the relationship
    adtTLL1FitMale <- gls(TL ~ MaleL1, 
                          data = clasperDFMale)
    
    plot(TL ~ MaleL1, 
         data = clasperDFMale,
         pch = 21,
         col = "black",
         bg = maleClasperGPA$col,
         cex = 4,
         cex.lab = 1.2,
         xlim = c(1.55, 1.64), #NEED TO CHANGE ACROSS VARIABLES
         ylim = c(0.25, 0.54), #NEED TO CHANGE ACROSS VARIABLES
         xlab = "log10 (L1)",
         ylab = "log10 (TL)")
    abline(a = coef(adtTLL1FitMale)[1],
           b = coef(adtTLL1FitMale)[2],
           lw = 3)
    
    summary(adtTLL1FitMale) #Check if the confidence intervals include 1 to assess allometric trajectory
    confint(adtTLL1FitMale)
    
    #Fit a gls model and plot the relationship
    adtTWL1FitMale <- gls(TW ~ MaleL1, 
                          data = clasperDFMale)
    
    plot(TW ~ MaleL1, 
         data = clasperDFMale,
         pch = 21,
         col = "black",
         bg = maleClasperGPA$col,
         cex = 4,
         cex.lab = 1.2,
         xlim = c(1.55, 1.64), #NEED TO CHANGE ACROSS VARIABLES
         ylim = c(0.08, 0.33), #NEED TO CHANGE ACROSS VARIABLES
         xlab = "log10 (L1)",
         ylab = "log10 (TW)")
    abline(a = coef(adtTWL1FitMale)[1],
           b = coef(adtTWL1FitMale)[2],
           lw = 3)
    
    summary(adtTWL1FitMale) #Check if the confidence intervals include 1 to assess allometric trajectory
    confint(adtTWL1FitMale)
    
