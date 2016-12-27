# Clear the workspace 
rm(list=ls())

# Set Working Directory 
path = "F:/OPIM-5503-Data Analytics using R-SECB13-1168/MultiVariateAnalysis/"
setwd(path)


# install.packages("car")  
# install.packages("RColorBrewer")
# install.packages("MASS")
# install.packages("splitstackshape")

library(splitstackshape)
library(car)
library(RColorBrewer)
library(MASS)

# Read Multivariate Data into R
wineData <- read.csv("wineData.csv",header = F)

# Summary Statistics 
summary(wineData)

# Stratified Sampling (80/20 Ratio)
StrataDataList <- stratified(indt = wineData,group = "V1",size = 0.2,bothSets = T)
TestData <- data.frame(StrataDataList[[1]])
TrainData <- data.frame(StrataDataList[[2]])
remove(StrataDataList)
wineData <- TrainData


# Scatter Plot Matrix 
scatterplotMatrix(wineData)


# Scatter Plot along with Wine Group Labels 
plot(wineData$V4,wineData$V5)
text(wineData$V4,wineData$V5,wineData$V1,cex = 0.7,pos = 4,col = "red")

# Profile Plot 
makeProfilePlot <- function(mylist,names)
{
  require(RColorBrewer)
  # find out how many variables we want to include
  TotalVariables <- length(mylist)
  # choose 'TotalVariables' random colors
  colors <- brewer.pal(TotalVariables,"Set1")
  # find out the minimum and maximum values of the variables:
  mymin <- 1e+20
  mymax <- 1e-20
  for (i in 1:TotalVariables)
  {
    vectori <- mylist[[i]]
    mini <- min(vectori)
    maxi <- max(vectori)
    if (mini < mymin) { mymin <- mini }
    if (maxi > mymax) { mymax <- maxi }
  }
  # plot the variables
  for (i in 1:TotalVariables)
  {
    vectori <- mylist[[i]]
    namei <- names[i]
    colori <- colors[i]
    if (i == 1) { plot(vectori,col=colori,type="l",ylim=c(mymin,mymax)) }
    else { points(vectori, col=colori,type="l") }
    lastxval <- length(vectori)
    lastyval <- vectori[length(vectori)]
    text((lastxval-10),(lastyval),namei,col="black",cex=0.9)
  }
}

names <- c("V2","V3","V4","V5","V6")
mylist <- list(wineData$V2,wineData$V3,wineData$V4,wineData$V5,wineData$V6)
makeProfilePlot(mylist,names)


names <- c("V7","V8","V9","V10")
mylist <- list(wineData$V7,wineData$V8,wineData$V9,wineData$V10)
makeProfilePlot(mylist,names)

names <- c("V11","V12","V13","V14")
mylist <- list(wineData$V11,wineData$V12,wineData$V13,wineData$V14)
makeProfilePlot(mylist,names)


# Calculate Summary Statistics 
sapply(wineData[2:14], sd)
sapply(wineData[2:14], mean)

# Calculate Summary Statistics for a particular wine group
WineType2 <- wineData[wineData$V1 == "2",]
sapply(WineType2[2:14], sd)
sapply(WineType2[2:14], mean)


####################################################################################################################
# Print Mean and SD By Wine Group
####################################################################################################################

# install.packages("knitr")
library(knitr)

printMeanAndSdByGroup <- function(variables,groupvariable)
{
  require(knitr)  

  # find the names of the variables
  variablenames <- c(names(groupvariable),names(as.data.frame(variables)))
  
  # within each group, find the mean of each variable
  groupvariable <- groupvariable[,1] # ensures groupvariable is not a list
  means <- aggregate(as.matrix(variables) ~ groupvariable, FUN = mean)
  names(means) <- variablenames
  print(kable(means))
  cat("Table 1: Means")
  
  # within each group, find the standard deviation of each variable:
  sds <- aggregate(as.matrix(variables) ~ groupvariable, FUN = sd)
  names(sds) <- variablenames
  print(kable(sds))
  cat("Table 2: Standard Deviations")
  
  # within each group, find the number of samples:
  samplesizes <- aggregate(as.matrix(variables) ~ groupvariable, FUN = length)
  names(samplesizes) <- variablenames
  print(kable(samplesizes))
  cat("Table 3: Sample sizes")
}


printMeanAndSdByGroup(wineData[2:14],wineData[1])


####################################################################################################################
# Calculate Winthin Group Variance Function
####################################################################################################################

calcWithinGroupsVariance <- function(variable,groupvariable)
{
  
  # find out how many values the group variable can take
  groupvariable2 <- as.factor(groupvariable[[1]])
  levels <- levels(groupvariable2)
  numlevels <- length(levels)
  # get the mean and standard deviation for each group:
  numtotal <- 0
  denomtotal <- 0
  for (i in 1:numlevels)
  {
    leveli <- levels[i]
    levelidata <- variable[groupvariable==leveli,]
    levelilength <- length(levelidata)
    # get the standard deviation for group i:
    sdi <- sd(levelidata)
    numi <- (levelilength - 1)*(sdi * sdi)
    denomi <- levelilength
    numtotal <- numtotal + numi
    denomtotal <- denomtotal + denomi
  }
  # calculate the within-groups variance
  Vw <- numtotal / (denomtotal - numlevels)
  return(Vw)
}



####################################################################################################################
# Calculate Between Group Variance Function
####################################################################################################################

calcBetweenGroupsVariance <- function(variable,groupvariable)
{
  # find out how many values the group variable can take
  groupvariable2 <- as.factor(groupvariable[[1]])
  levels <- levels(groupvariable2)
  numlevels <- length(levels)
  # calculate the overall grand mean:
  grandmean <- mean(as.vector(variable[,1]))
  
  # get the mean and standard deviation for each group:
  numtotal <- 0
  denomtotal <- 0
  for (i in 1:numlevels)
  {
    leveli <- levels[i]
    levelidata <- variable[groupvariable==leveli,]
    levelilength <- length(levelidata)
    # get the mean and standard deviation for group i:
    meani <- mean(levelidata)
    sdi <- sd(levelidata)
    numi <- levelilength * ((meani - grandmean)^2)
    denomi <- levelilength
    numtotal <- numtotal + numi
    denomtotal <- denomtotal + denomi
  }
  # calculate the between-groups variance
  Vb <- numtotal / (numlevels - 1)
  Vb <- Vb[[1]]
  return(Vb)
}

####################################################################################################################
# Calculate Separations Funciton
####################################################################################################################

calcSeparations <- function(variables,groupvariable)
{
  # find out how many variables we have
  variables <- as.data.frame(variables)
  numvariables <- length(variables)
  # find the variable names
  variablenames <- colnames(variables)
  # calculate the separation for each variable
  for (i in 1:numvariables)
  {
    variablei <- variables[i]
    variablename <- variablenames[i]
    Vw <- calcWithinGroupsVariance(variablei, groupvariable)
    Vb <- calcBetweenGroupsVariance(variablei, groupvariable)
    sep <- Vb/Vw
    print(paste("variable",variablename,"Vw=",Vw,"Vb=",Vb,"separation=",sep))
  }
}

####################################################################################################################
# Calculate Separations
####################################################################################################################

calcSeparations(wineData[2:14],wineData[1])

####################################################################################################################
# Calculate Highly Correlated Variables Function
####################################################################################################################

mosthighlycorrelated <- function(mydataframe,numtoreport){

  require(Hmisc)
  
  # Find Correlations 
  corMatrix <- rcorr(as.matrix(mydataframe))$r
  # set the correlations on the diagonal or lower triangle to zero,
  diag(corMatrix) <- 0
  corMatrix[lower.tri(corMatrix)] <- 0
  # flatten the matrix into a dataframe
  corrFm <- as.data.frame(as.table(corMatrix))
  # assign column names
  names(corrFm) <- c("First.Variable", "Second.Variable","Correlation")
  
  # Find p-values related to Correlations 
  pValMatrix <- rcorr(as.matrix(mydataframe))$P
  # set the p-values on the diagonal or lower triangle to zero,
  diag(pValMatrix) <- 0
  pValMatrix[lower.tri(pValMatrix)] <- 0
  # flatten the matrix into a dataframe 
  pValFm <- as.data.frame(as.table(pValMatrix))
  # assign column names
  names(pValFm) <- c("First.Variable", "Second.Variable","p-value")
  
  
  # Create a single frame 
  fm <- data.frame(fm$First.Variable,fm$Second.Variable,fm$Correlation,fm1$`p-value`)
  # assign column names
  names(fm) <- c("First.Variable", "Second.Variable","Correlation","p-value")
  # Format p-value column
  fm$`p-value` <- format(round(fm$`p-value`, 4), nsmall = 4)
  
  
  # sort and print the top n correlations
  head(fm[order(abs(fm$Correlation),decreasing=T),],n=numtoreport)
  
}


####################################################################################################################
# Calculate Highly Correlated Variables Function
####################################################################################################################

library(Hmisc)
mosthighlycorrelated(wineData[2:14], 10)


####################################################################################################################
# Standardising Variables 
####################################################################################################################

standardisedconcentrations <- as.data.frame(scale(wineData[2:14]))

# Standardised Variables will have standard deviation 0 and mean 1 
sapply(standardisedconcentrations,sd)
sapply(standardisedconcentrations,mean)

####################################################################################################################
# Principal Componant Analysis 
####################################################################################################################

wineData.pca <- prcomp(standardisedconcentrations)
summary(wineData.pca)

screeplot(wineData.pca, type="lines")

# The most obvious change in slope in the scree plot occurs at component 4, which is the "elbow" of the scree
# plot. Therefore, it cound be argued based on the basis of the scree plot that the first three components should be
# retained.

(wineData.pca$sdev)^2

# We see that the variance is above 1 for principal components 1, 2, and 3 (which have variances 4.71, 2.50, and
# 1.45, respectively). Therefore, using Kaiser's criterion, we would retain the first three principal components.


####################################################################################################################
# Loadings for the Principal Components
####################################################################################################################

wineData.pca$rotation[,1]

# This means that the first principal component is a linear combination of the variables: -0.144*Z2 + 0.245*Z3 +
# 0.002*Z4 + 0.239*Z5 - 0.142*Z6 - 0.395*Z7 - 0.423*Z8 + 0.299*Z9 -0.313*Z10 + 0.089*Z11 - 0.297*Z12 -
# 0.376*Z13 - 0.287*Z14, where Z2, Z3, Z4...Z14 are the standardised versions of the variables V2, V3, V4...V14
# (that each have mean of 0 and variance of 1).

# The first principal component has highest (in absolute value) loadings for V8 (-0.423), V7 (-0.395), V13 (-0.376),
# V10 (-0.313), V12 (-0.297), V14 (-0.287), V9 (0.299), V3 (0.245), and V5 (0.239). The loadings for V8, V7,
# V13, V10, V12 and V14 are negative, while those for V9, V3, and V5 are positive. Therefore, an interpretation
# of the first principal component is that it represents a contrast between the concentrations of V8, V7, V13, V10,
# V12, and V14, and the concentrations of V9, V3 and V5.

# Note that the square of the loadings sum to 1, as this is a constraint used in calculating the loadings:
# sum((wine.pca$rotation[,1])^2)

# First Principal Componant 
wineData.pca$x[,1]


wine.pca$rotation[,2]

# This means that the second principal component is a linear combination of the variables: 0.484*Z2 + 0.225*Z3 +
# 0.316*Z4 - 0.011*Z5 + 0.300*Z6 + 0.065*Z7 - 0.003*Z8 + 0.029*Z9 + 0.039*Z10 + 0.530*Z11 - 0.279*Z12 -
# 0.164*Z13 + 0.365*Z14, where Z1, Z2, Z3...Z14 are the standardised versions of variables V2, V3, ... V14 that
# each have mean 0 and variance 1.

# The second principal component has highest loadings for V11 (0.530), V2 (0.484), V14 (0.365), V4 (0.316), V6
# (0.300), V12 (-0.279), and V3 (0.225). The loadings for V11, V2, V14, V4, V6 and V3 are positive, while the
# loading for V12 is negative. Therefore, an interpretation of the second principal component is that it represents a
# contrast between the concentrations of V11, V2, V14, V4, V6 and V3, and the concentration of V12. Note that
# the loadings for V11 (0.530) and V2 (0.484) are the largest, so the contrast is mainly between the concentrations
# of V11 and V2, and the concentration of V12.


####################################################################################################################
# Scatterplots of the Principal Components
####################################################################################################################

plot(wineData.pca$x[,1],wineData.pca$x[,2]) # make a scatterplot
text(wineData.pca$x[,1],wineData.pca$x[,2], wineData$V1, cex=0.7, pos=4, col="red") # add labels


# Above, we interpreted the first principal component as a contrast between the concentrations of V8, V7, V13,
# V10, V12, and V14, and the concentrations of V9, V3 and V5. We can check whether this makes sense in terms
# of the concentrations of these chemicals in the different cultivars, by printing out the means of the standardised
# concentration variables in each cultivar, using the "printMeanAndSdByGroup()" function (see above)

printMeanAndSdByGroup(standardisedconcentrations,wineData[1])

# In cultivar 1, the mean values of V8 (0.954), V7 (0.871), V13 (0.769), V10 (0.539), V12 (0.458) and V14 (1.171) are very high
# compared to the mean values of V9 (-0.577), V3 (-0.292) and V5 (-0.736). In cultivar 3, the mean values of V8
# (-1.249), V7 (-0.985), V13 (-1.307), V10 (-0.764), V12 (-1.202) and V14 (-0.372) are very low compared to the
# mean values of V9 (0.688), V3 (0.893) and V5 (0.575). Therefore, it does make sense that principal component 1
# is a contrast between the concentrations of V8, V7, V13, V10, V12, and V14, and the concentrations of V9, V3
# and V5; and that principal component 1 can separate cultivar 1 from cultivar 3.

# Similarly PC2 seperates Groyp 2 from Groups 1 & 3 

# The purpose of principal component analysis is to find the best low-dimensional representation of the variation in
# a multivariate data set. For example, in the wine data set, we have 13 chemical concentrations describing wine
# samples from three cultivars. By carrying out a principal component analysis, we found that most of the variation
# in the chemical concentrations between the samples can be captured using the first two principal components,
# where each of the principal components is a particular linear combination of the 13 chemical concentrations.


####################################################################################################################
# Linear Discriminant Analysis
####################################################################################################################

# The purpose of linear discriminant analysis (LDA) is to find the linear combinations of the original variables
# (the 13 chemical concentrations here) that gives the best possible separation between the groups (wine cultivars
# here) in our data set. Linear discriminant analysis is also known as "canonical discriminant analysis", or simply"discriminant analysis".


# If we want to separate the wines by cultivar, the wines come from three different cultivars, so the number of groups
# (G) is 3, and the number of variables is 13 (13 chemicals' concentrations; p = 13). The maximum number of useful
# discriminant functions that can separate the wines by cultivar is the minimum of G-1 and p, and so in this case it
# is the minimum of 2 and 13, which is 2. Thus, we can find at most 2 useful discriminant functions to separate the
# wines by cultivar, using the 13 chemical concentration variables.


wineData.lda <- lda(wineData$V1 ~ wineData$V2 + wineData$V3 + wineData$V4 + wineData$V5 + wineData$V6 + wineData$V7 + wineData$V8
                     +wineData$V9 + wineData$V10 + wineData$V11 + wineData$V12 + wineData$V13 + wineData$V14)

wineData.lda

wineData.lda$scaling[,1]

# This means that the first discriminant function is a linear combination of the variables: -0.403*V2 + 0.165*V3 -
# 0.369*V4 + 0.155*V5 - 0.002*V6 + 0.618*V7 - 1.661*V8 - 1.496*V9 + 0.134*V10 + 0.355*V11 - 0.818*V12
# - 1.158*V13 - 0.003*V14, where V2, V3, ... V14 are the concentrations of the 14 chemicals found in the wine
# samples.


# Calculate First LD Values 
wineData.lda.values <- predict(wineData.lda, wineData[2:14])
wineData.lda.values$x[,1]

# Calculate LDA &  First LD Values for Test Data

TestData.lda <- lda(TestData$V1 ~ TestData$V2 + TestData$V3 + TestData$V4 + TestData$V5 + TestData$V6 + TestData$V7 + TestData$V8
                    +TestData$V9 + TestData$V10 + TestData$V11 + TestData$V12 + TestData$V13 + TestData$V14)


TestData.lda.values <- predict(TestData.lda,TestData[2:14])
TestData.lda.values$x[,1]

#########################################################################################################################
# LDA on Group-Standardised Variables
#########################################################################################################################

# In linear discriminant analysis, the standardised version of an input variable is defined so that it has mean zero
# and within-groups variance of 1. Thus, we can calculate the "group-standardised" variable by subtracting the
# mean from each value of the variable, and dividing by the within-groups standard deviation


groupStandardise <- function(variables, groupvariable)
{
  # find out how many variables we have
  variables <- as.data.frame(variables)
  numvariables <- length(variables)
  # find the variable names
  variablenames <- colnames(variables)
  # calculate the group-standardised version of each variable
  for (i in 1:numvariables)
  {
    variablei <- variables[i]
    variablei_name <- variablenames[i]
    variablei_Vw <- calcWithinGroupsVariance(variablei, groupvariable)
    variablei_mean <- mean(variablei[,1])
    variablei_new <- (variablei - variablei_mean)/(sqrt(variablei_Vw))
    data_length <- nrow(variablei)
    if (i == 1) { variables_new <- data.frame(row.names=seq(1,data_length)) }
    variables_new[`variablei_name`] <- variablei_new
  }
  return(variables_new)
}


groupstandardisedconcentrations <- groupStandardise(wineData[2:14], wineData[1])


wineData.lda2 <- lda(wineData$V1 ~ groupstandardisedconcentrations$V2 + groupstandardisedconcentrations$V3 + groupstandardisedconcentrations$V4 + groupstandardisedconcentrations$V5 + groupstandardisedconcentrations$V6 + groupstandardisedconcentrations$V7 + groupstandardisedconcentrations$V8 + groupstandardisedconcentrations$V9 + groupstandardisedconcentrations$V10 + groupstandardisedconcentrations$V11 + groupstandardisedconcentrations$V12 + groupstandardisedconcentrations$V13 + groupstandardisedconcentrations$V14)


wineData.lda2



# In the first discriminant function calculated for the group-standardised variables, the largest loadings (in absolute)
# value are given to V8 (-0.871), V11 (0.537), V13 (-0.464), V14 (-0.464), and V5 (0.438). The loadings for V8,
# V13 and V14 are negative, while those for V11 and V5 are positive. Therefore, the discriminant function seems
# to represent a contrast between the concentrations of V8, V13 and V14, and the concentrations of V11 and V5.
# We saw above that the individual variables which gave the greatest separations between the groups were V8
# (separation 233.93), V14 (207.92), V13 (189.97), V2 (135.08) and V11 (120.66). These were mostly the same
# variables that had the largest loadings in the linear discriminant function (loading for V8: -0.871, for V14: -0.464, for V13: -0.464, for V11: 0.537).
# We found above that variables V8 and V11 have a negative between-groups covariance (-60.41) and a positive
# within-groups covariance (0.29). When the between-groups covariance and within-groups covariance for two
# variables have opposite signs, it indicates that a better separation between groups can be obtained by using a linear
# combination of those two variables than by using either variable on its own.
# Thus, given that the two variables V8 and V11 have between-groups and within-groups covariances of opposite
# signs, and that these are two of the variables that gave the greatest separations between groups when used individually,
# it is not surprising that these are the two variables that have the largest loadings in the first discriminant
# function.


wineData.lda.values2 <- predict(wineData.lda2, groupstandardisedconcentrations)
wineData.lda.values2$x[,1]

# We can see that although the loadings are different for the first discriminant functions calculated using unstandardised
# and group-standardised data, the actual values of the first discriminant function are the same.

calcSeparations(wineData.lda.values$x,wineData[1])
# the total separation is the sum of these, which is (794.652200566216+361.241041493455=1155.893)
# 1155.89, rounded to two decimal places. Therefore, the "percentage separation" achieved by the first discriminant
# function is (794.652200566216*100/1155.893=) 68.75%, and the percentage separation achieved by the second
# discriminant function is (361.241041493455*100/1155.893=) 31.25%.



# The "proportion of trace" that is printed when you type "wineData.lda2" OR "wineData.lda" (the variable returned by the lda() function) is
# the percentage separation achieved by each discriminant function. For example, for the wine data we get the same
# values as just calculated (68.75% and 31.25%):
wineData.lda2

#########################################################################################################
# A Stacked Histogram of the LDA Values
#########################################################################################################

ldahist(data = wineData.lda.values$x[,1], g=wineData$V1)

# We can see from the histogram that cultivars 1 and 3 are well separated by the first discriminant function, since
# the values for the first cultivar are between -6 and -1, while the values for cultivar 3 are between 2 and 6, and so
# there is no overlap in values.

#However, the separation achieved by the linear discriminant function on the training set may be an overestimate.
#To get a more accurate idea of how well the first discriminant function separates the groups, we would need to see
#a stacked histogram of the values for the three cultivars using some unseen "test set", that is, using a set of data
#that was not used to calculate the linear discriminant function.

#We see that the first discriminant function separates cultivars 1 and 3 very well, but does not separate cultivars 1
#and 2, or cultivars 2 and 3, so well.

#We therefore investigate whether the second discriminant function separates those cultivars, by making a stacked
#histogram of the second discriminant function's values:

ldahist(data = wineData.lda.values$x[,2], g=wineData$V1)

# We see that the second discriminant function separates cultivars 1 and 2 quite well, although there is a little overlap
# in their values. Furthermore, the second discriminant function also separates cultivars 2 and 3 quite well, although
# again there is a little overlap in their values so it is not perfect.


##########################################################################################################################
# Scatterplots of the Discriminant Functions
##########################################################################################################################
plot(wineData.lda.values$x[,1],wineData.lda.values$x[,2])
text(wineData.lda.values$x[,1],wineData.lda.values$x[,2],wineData$V1,cex=0.7,pos=4,col="red")


##########################################################################################################################
# Allocation Rules and Misclassification Rate
##########################################################################################################################

printMeanAndSdByGroup(wineData.lda.values$x,wineData[1])

# We find that the mean value of the first discriminant function is -3.42248851 for cultivar 1, -0.07972623 for
#cultivar 2, and 4.32473717 for cultivar 3. The mid-way point between the mean values for cultivars 1 and 2 is
#(-3.42248851-0.07972623)/2=-1.751107, and the mid-way point between the mean values for cultivars 2 and 3 is
#(-0.07972623+4.32473717)/2 = 2.122505.
#Therefore, we can use the following allocation rule:
#  . if the first discriminant function is <= -1.751107, predict the sample to be from cultivar 1
#  . if the first discriminant function is > -1.751107 and <= 2.122505, predict the sample to be from cultivar 2
#  . if the first discriminant function is > 2.122505, predict the sample to be from cultivar 3




calcAllocationRuleAccuracy <- function(ldavalue, groupvariable, cutoffpoints)
{
  # find out how many values the group variable can take
  groupvariable2 <- as.factor(groupvariable[[1]])
  levels <- levels(groupvariable2)
  numlevels <- length(levels)
  # calculate the number of true positives and false negatives for each group
  numlevels <- length(levels)
  for (i in 1:numlevels)
  {
    leveli <- levels[i]
    levelidata <- ldavalue[groupvariable==leveli]
    # see how many of the samples from this group are classified in each group
    for (j in 1:numlevels)
    {
      levelj <- levels[j]
      if (j == 1)
      {
        cutoff1 <- cutoffpoints[1]
        cutoff2 <- "NA"
        results <- summary(levelidata <= cutoff1)
      }
      else if (j == numlevels)
      {
        cutoff1 <- cutoffpoints[(numlevels-1)]
        cutoff2 <- "NA"
        results <- summary(levelidata > cutoff1)
      }
      else
      {
        cutoff1 <- cutoffpoints[(j-1)]
        cutoff2 <- cutoffpoints[(j)]
        results <- summary(levelidata > cutoff1 & levelidata <= cutoff2)
      }
      trues <- results["TRUE"]
      trues <- trues[[1]]
      print(paste("Number of samples of group",leveli,"classified as group",levelj," : ",
                  trues,"(cutoffs:",cutoff1,",",cutoff2,")"))
    }
  }
}

calcAllocationRuleAccuracy(wineData.lda.values$x[,1], wineData[1], c(-1.751107, 2.122505))

calcAllocationRuleAccuracy(TestData.lda.values$x[,1], TestData[1], c(-1.751107, 2.122505))


# There are 3+5+1=9 wine samples that are misclassified, out of (56+3+5+65+1+48=) 178 wine samples: 3 samples
# from cultivar 1 are predicted to be from cultivar 2, 5 samples from cultivar 2 are predicted to be from cultivar
# 1, and 1 sample from cultivar 2 is predicted to be from cultivar 3. Therefore, the misclassification rate is 9/178,
# or 5.1%. The misclassification rate is quite low, and therefore the accuracy of the allocation rule appears to be
# relatively high.

# However, this is probably an underestimate of the misclassification rate, as the allocation rule was based on this
# data (this is the "training set"). If we calculated the misclassification rate for a separate "test set" consisting of data
# other than that used to make the allocation rule, we would probably get a higher estimate of the misclassification
# rate.

