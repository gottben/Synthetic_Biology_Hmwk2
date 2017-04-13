# Uploads library data on nif gene expression levels and pathway yield - don't forget to replace filePath with the path to library_yields.csv on your system
setwd("C:\\Users\\gottb\\OneDrive\\Documents\\Graduate_School\\Synthetic Biology\\Hmwk\\Hmwk2")
myData <- read.csv("library_yields.csv")

# Log 10 transformation of data to address nonlinearity and heteroscedasticity

transformData <- log10(myData)

# Multiple linear regression on transformed data

# The gene expression terms on the right hand side of the equation (and all other equations in this script) should match the genes tested by your library

# You may also add interaction terms of the form nifB:nifH, nifN:nifE:nifM, etc. to the equation

fit <- lm(yield ~ nifB + nifN + nifE + nifM + nifU + nifH + nifD + nifS + nifK, data=transformData)

# Displays summary of regression model, including coefficients, p-values, R^2, etc.

summary(fit)

# Creates residual plots for regression model

layout(matrix(c(1,2,3,4),2,2))

plot(fit)

# Fits many regression models to transformed data and ranks them by adjusted R^2 - see suggested reading 3

library(leaps) 

leaps <- regsubsets(yield ~ (nifB + nifN + nifE + nifM + nifU + nifH + nifD + nifS + nifK)^2, data=transformData, nbest=6, nvmax=11, force.in=c(1,2,3,4,5,6))

plot(leaps, scale="adjr2")

# Uses regression model to predict yield for a new set of gene expression levels

newData <- data.frame(nifB = log10(5000), nifN = log10(6000), nifE = log10(7000), nifM = log10(3000), nifU = log10(5000), nifH = log10(4000), nifD = log10(5000), nifS = log10(6000), nifK = log10(7000))

predictedYield <- 10^predict(fit, newData)