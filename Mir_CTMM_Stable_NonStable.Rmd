---
title: 'Predictive Power of Monocyte Expressed microRNA on CAD: Stable vs Non Stable'
author: "Huayu"
output:
  html_document:
    df_print: paged
number_sections: yes
---


```{r setup}
knitr::opts_chunk$set(echo = TRUE)
library(knitr)
library(ggplot2)
library(reshape2)
library(caret)
library(doParallel)
library(pROC)


library(caTools)
#library(pls) # pls
library(MASS) #lda
library(naivebayes) # naivebayes
library(ipred); library(plyr); library(e1071) # bagged tree
library(randomForest) #rf
library(gbm) # gbm
#library(xgboost) #xgboost
```

#Data Input, preparation and descriptives
##Import the dataframe of patient parameter and NGC expression

```{r}
ctmmMir <- read.delim(file = "CTMM_miR_all.txt", stringsAsFactors = FALSE)
str(ctmmMir[,1:6])
kable(head(ctmmMir[,1:6]))
```

```{r}
names(ctmmMir)
```


##Split into responce and NGC expression
Split data to clinical parameters and NGC expression
```{r}
res <- ctmmMir[,1:121]
row.names(res) <- ctmmMir$CC_number
ftrExprs <- ctmmMir[,122:213]
row.names(ftrExprs) <- ctmmMir$CC_number
head(names(res))
head(names(ftrExprs))
```



##Descriptives of the data

Missing Values (Only one missing for confirmed diagnosis)
```{r}
numNA <- sapply(ctmmMir, function(x) sum(is.na(x)))
numNA[!(numNA == 0)]
```

Different types of variables (continuous, categrical or descriptive)
```{r}
levelsVar <- sapply(X = res, FUN = function(x) { nlevels(factor(x))})


catVar <- levelsVar <10
numVar <- !catVar &  sapply(X = res, FUN = function(x) { is.numeric(x)})
charVar <- !(catVar | numVar)
```

Categrical variables are:
```{r}
names(res)[catVar]
```

Continuous variables are:
```{r}
names(res)[numVar]
```

Descriptive variables are:
```{r}
names(res)[charVar]
```

Tables of categorical variables
```{r}
catTables <- sapply(res[catVar], function(x) table(as.factor(x)))
lapply(catTables, function(x) x)
```


Distributions of numeric variables
```{r, fig.height = 10, fig.width = 7}
numVarGG <- melt(res[,numVar])
numDist <- ggplot(data = numVarGG) +
            geom_histogram(mapping = aes(x = value), bins = 15, na.rm = TRUE) +
            facet_wrap(facets = ~variable, ncol = 4, scales = "free")
            
numDist

```

##Distribution of Expression of NGC

```{r, fig.height = 12, fig.width = 7}
ftrExprsGG <- melt(ftrExprs)
ftrDist <- ggplot(data = ftrExprsGG) +
            geom_histogram(mapping = aes(x = value), bins = 15) +
            facet_wrap(facets = ~variable, ncol = 7, scales = "free") +
            theme(
              axis.text.x = element_text(size = 8, angle = 45, hjust = 1, vjust =  1),
              axis.text.y = element_text(size = 8))
ftrDist

```

Box plot of expression of all NGCs

```{r, fig.height = 5, fig.width = 10}
ftrBox <- ggplot(data = ftrExprsGG) +
            geom_boxplot(mapping = aes(x = variable, y = value)) +
            theme(
              axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1)
            )
ftrBox
```

#Data Imputation for NGC expression

```{r}
ftrExprsImp <- preProcess(x = ftrExprs, method = "bagImpute")
  
ftrExprs <- predict(object = ftrExprsImp, newdata = ftrExprs)
```


#Data paraperation for model fitting (Confirmed diagnosis as the responce)
##Remove NA values from the response

All diseased will be put in the same group

```{r}
ftrExprsNarm <- ftrExprs[ctmmMir$Conf_Diagn %in% 1:3, ]
resNarm <- res[ctmmMir$Conf_Diagn %in% 1:3, ]
#Put the named response into a vector
resNamed <- resNarm$Conf_Diagn
resNamed[resNamed == 1] <- 0
resNamed[resNamed == 2] <- 1
resNamed[resNamed == 3] <- 1
resNamed <- factor(resNamed, levels = c(0,1), labels = c("Stable_CVD", "Unstable_CVD"))
nrow(ftrExprsNarm)
table(resNamed)
```

##Violin plot of features vs response

Violin plot is more handy since the response is categorical. Mean and sd are indicated in red color.

```{r}
ftrExprsNarmGG <- melt(ftrExprsNarm)
ftrExprsNarmGG$response <- resNamed
```

```{r, fig.height = 14 , fig.width = 13}
ftrResScat <- ggplot(data = ftrExprsNarmGG, mapping = aes(x = response, y = value)) +
                geom_violin() +
                stat_summary(fun.data = mean_sdl, geom = "pointrange", size = 0.3, color = "red", fun.args = list(mult = 1)) +
                facet_wrap(facets = ~variable, ncol = 11, scales = "free") +
                theme(
                  axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45)
                )
ftrResScat
```


##Vacanoplot



```{r, fig.height = 5, fig.width = 4}
ftrParam <- data.frame(Pvalue = -log10(sapply(ftrExprsNarm, function(x) t.test(x ~ resNamed)$p.value)))
ftrParam$FC <- log2(sapply(ftrExprsNarm, function(x) t.test(x ~ resNamed)$estimate[2]/t.test(x ~ resNamed)$estimate[1]))
ftrParam$label <- NA
ftrParam$label[ftrParam$Pvalue > -log10(0.05) & abs(ftrParam$FC) > 0.0125] <- row.names(ftrParam)[ftrParam$Pvalue > -log10(0.05) & abs(ftrParam$FC) > 0.0125]
ftrParam$color <- as.numeric(!is.na(ftrParam$label)) + 1

ggplot(ftrParam) +
  geom_point(mapping = aes(x = FC, y = Pvalue), color =  ftrParam$color) +
  geom_text(mapping = aes(x = FC, y = Pvalue, label = label), color = ftrParam$color, hjust = -0.1, vjust = 0.1, angle = 45, na.rm = TRUE, size = 3) +
  geom_hline(mapping = aes(yintercept = -log10(0.05)), linetype = 2) +
  geom_vline(xintercept = 0.0125, linetype = 2) +
  geom_vline(xintercept = -0.0125, linetype = 2) +
  coord_cartesian(xlim = c(-0.06, 0.06), ylim = c(0, 4)) +
  xlab("log(Fold Change)") +
  ylab("-log(p-value)")
```

##P-value and Expression Level
Does the p-value rely on expression of guidance cues?

```{r, fig.height = 5, fig.width = 4}
ftrParam$meanExprs <- sapply(ftrExprsNarm, mean)

pCut <- 0.2

ftrParam$labelExprs <- NA
ftrParam$labelExprs[ftrParam$Pvalue > -log10(pCut)] <- row.names(ftrParam)[ftrParam$Pvalue > -log10(pCut)]
ftrParam$colorExprs <- as.numeric(!is.na(ftrParam$labelExprs)) + 1

ggplot(ftrParam) +
  geom_point(mapping = aes(x = meanExprs, y = Pvalue), color = ftrParam$colorExprs) +
  geom_text(mapping = aes(x = meanExprs, y = Pvalue, label = labelExprs), color = ftrParam$colorExprs, hjust = -0.1, vjust = 0.1, angle = 45, na.rm = TRUE, size = 3) +
  geom_hline(yintercept = -log10(pCut), linetype = 2) +
  coord_cartesian(ylim = c(0,3.5))
```

Distribution of -log10 p
```{r}
hist(ftrParam$Pvalue)
```



```{r}
ftrIn <- ftrParam$Pvalue > -log10(pCut)
sum(ftrIn)
```

```{r}
ftrExprsIn <- ftrExprsNarm[,ftrIn]
names(ftrExprsIn)
```

##Correlation between features

```{r, fig.height = 6, fig.width = 6}
ftrCorr <- cor(ftrExprsIn)
heatmap(ftrCorr, revC = TRUE, scale = "none")

findCorrelation(ftrCorr, cutoff = 0.8)
```

##Filtering by correlations of all features



#Confounding factors

Age and Sex are very common confounding factors in biomarker studies. We examine them also here. 

##Examining age
```{r}
ggplot(resNarm, aes(x = resNamed, y = age)) +
  geom_violin(color = "red") +
  geom_dotplot(binaxis = "y", dotsize = 0.5, stackdir =  "center") +
  geom_text(mapping = aes(x = 1, y = 90, label = paste("p =", round(t.test(resNarm$age ~ resNamed)$p.value, 3))))
```

Age is not very different amoung stable and unstable group. Therefore, age is unlikely to be a confounding factor in this study.

```{r, fig.height = 14, fig.width = 7}
ftrExprsNarmGG$age <- resNarm$age
ggplot(ftrExprsNarmGG, mapping = aes(x = age, y = value)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm") +
  facet_wrap(facets = ~variable, scales = "free", ncol = 7) +
                theme(
                  axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45, size = 6)
                )
```

Age doesn't correlate much with ftr expression, so that age is again not a confounding factor in this modeling. But age will be included anyway to cancel any confounding effect. 

##Examining sex
```{r, fig.height = 14, fig.width = 7}
sexRes <- data.frame(sex = factor(resNarm$Sex, 0:1, c("female", "male")))
sexRes$res <- resNamed
table(sexRes)
```


```{r}
chisq.test(table(sexRes))
```

Sex is not significantly biased to any group. Sex will also be included anyway to cancell confounding effect.

```{r, fig.height = 16, fig.width = 7}
ftrExprsNarmGG$Sex <- sexRes$sex
ggplot(ftrExprsNarmGG, aes(x = Sex, y = value)) +
  geom_violin() +
  stat_summary(fun.data = mean_sdl, geom = "pointrange", size = 0.3, color = "red", fun.args = list(mult = 1)) +
  facet_wrap(facets = ~variable, ncol = 7, scales = "free") +
  theme(
        axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45)
      )
```



#Data partition using the samples with a valid confirmed diagnostics, settings of model parameters

##Data partition

```{r}
set.seed(4000)
trainingSamples <- createDataPartition(y = resNamed, p = 0.8, list = FALSE)
length(trainingSamples)

resTr <- resNamed[trainingSamples]
resTest <- resNamed[-trainingSamples]

table(resTr)
table(resTest)

#Age and Sex will be included
fea <- cbind(age = as.numeric(resNarm$age), Sex = factor(resNarm$Sex), ftrExprsIn)
fea$Sex <- as.numeric(fea$Sex) -1
feaTr <- fea[trainingSamples, ] 
feaTest <- fea[-trainingSamples, ]



feaTr[1:6, 1:6]
feaTest[1:6, 1:6]


```

##Model Training parameters (cross validation ...)

10 fold cross-validation used

```{r}
  ctrl = trainControl(method = "repeatedcv", repeats = 10, classProbs = TRUE)
```



#Model Fitting


##Function for model report


```{r}
modelReport <- function(model = NULL, obsTr = resTr, obsTest = resTest, newTr = feaTr, newTest = feaTest)
{
  predTr <- predict(object = model, newdata = newTr)
  predTest <- predict(object = model, newdata = newTest)
  predTrProb <- predict(object = model, newdata = newTr, type = "prob")
  predTestProb <- predict(object = model, newdata = newTest, type = "prob")
  
  feaImp <- varImp(model)$importance
  feaImp$fea <- row.names(feaImp)
  colnames(feaImp)[1] <- "Imp"
  feaImp <- feaImp[order(feaImp$Imp, decreasing = T), ]
  feaImp$fea <- factor(feaImp$fea, levels = row.names(feaImp))
  
  feaImpVis <- ggplot(data = feaImp[1:10,], mapping = aes(x = fea, y = Imp)) +
      geom_col() +
      xlab(label = NULL) +
      ylab(label = paste("Variable Importance", model$modelInfo$label)) +
      theme(
        axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45)
      )
    
  cmGG <- expand.grid(y = -1:-12, x = 1:7)
cmGG$label <- c(NA, NA, "Training Set", "Prediction", "Stable_CVD", "Unstable_CVD",
NA, "Training Set", "Accuracy", "Sensitivity", "Specificity", "Kappa",
"Confusion Matrix and Statistics", NA, "Reference", "Stable_CVD", "Training TN", "Training FP",
NA, NA, "Training Accuracy", "Training Sensitivity", "Training Specificity", "Training Kappa",
NA, NA, NA, "Unstable_CVD", "Training FN", "Training TP",
NA, NA, NA, NA, NA, NA,
NA, NA, NA, NA, NA, NA,
NA, NA, NA, NA, NA, NA,
NA, NA, "Test Set", "Prediction", "Stable_CVD", "Unstable_CVD",
NA, "Test Set", "Accuracy", "Sensitivity", "Specificity", "Kappa",
NA, NA, "Reference", "Stable_CVD", "Test TN", "Test FP",
NA, NA, "Test Accuracy", "Test Sensitivity", "Test Specificity", "Test Kappa",
NA, NA, NA, "Unstable_CVD", "Test FN", "Test TP",
NA, NA, NA, NA, NA, NA)

cmGG$label[c(5, 16, 53, 64)] <- as.character(model$levels[1])
cmGG$label[c(6, 28, 54, 76)] <- as.character(model$levels[2])

cmTr <- confusionMatrix(data = predTr, reference = obsTr, positive = levels(resTest)[2])
cmTest <- confusionMatrix(data = predTest, reference = obsTest, positive = levels(resTest)[2])

cmGG$label[17] <- cmTr$table[1,1]
cmGG$label[18] <- cmTr$table[2,1]
cmGG$label[29] <- cmTr$table[1,2]
cmGG$label[30] <- cmTr$table[2,2]
cmGG$label[65] <- cmTest$table[1,1]
cmGG$label[66] <- cmTest$table[2,1]
cmGG$label[77] <- cmTest$table[1,2]
cmGG$label[78] <- cmTest$table[2,2]
cmGG$label[21] <- round(cmTr$overall[1], 4)
cmGG$label[22] <- round(cmTr$byClass[1], 4)
cmGG$label[23] <- round(cmTr$byClass[2], 4)
cmGG$label[24] <- round(cmTr$overall[2], 4)
cmGG$label[69] <- round(cmTest$overall[1], 4)
cmGG$label[70] <- round(cmTest$byClass[1], 4)
cmGG$label[71] <- round(cmTest$byClass[2], 4)
cmGG$label[72] <- round(cmTest$overall[2], 4)
  
cmPlot <- ggplot(cmGG) +
  geom_text(aes(x, y, label = label), hjust = 1, na.rm = TRUE) +
  coord_cartesian(xlim = c(0,8), ylim = c(0,-13)) +
  geom_segment(mapping = aes(x = -0.05, y = -2.5, xend = 7.1, yend = -2.5), size = 1) + #Top border
  geom_segment(mapping = aes(x = -0.05, y = -4.5, xend = 7.1, yend = -4.5)) + #Middel border
  geom_segment(mapping = aes(x = -0.05, y = -6.5, xend = 7.1, yend = -6.5), size = 1) + #Bottom border
  geom_segment(mapping = aes(x = 1.1, y = -3.5, xend = 3.1, yend = -3.5)) + #Left reference
  geom_segment(mapping = aes(x = 5.1, y = -3.5, xend = 7.1, yend = -3.5)) + #Right reference
  geom_segment(mapping = aes(x = -0.05, y = -7.5, xend = 2.1, yend = -7.5), size = 1) +
  geom_segment(mapping = aes(x = -0.05, y = -8.5, xend = 2.1, yend = -8.5)) +
  geom_segment(mapping = aes(x = -0.05, y = -12.5, xend = 2.1, yend = -12.5), size = 1) +
  geom_segment(mapping = aes(x = 4.1, y = -7.5, xend = 6.1, yend = -7.5), size = 1) +
  geom_segment(mapping = aes(x = 4.1, y = -8.5, xend = 6.1, yend = -8.5)) +
  geom_segment(mapping = aes(x = 4.1, y = -12.5, xend = 6.1, yend = -12.5), size = 1) +
  
 theme(line = element_blank(),
       text = element_blank(),
       title = element_blank(),
       plot.background = element_blank(),
       panel.background = element_blank()
       )
  
  rocCurve <- roc(response = obsTest, predictor = predTestProb[,1])
  
  return(list(
              model,
              feaImpVis,
              cmTr,
              cmTest,
              plot(rocCurve, legacy.axes = TRUE),
              cmPlot
              )
         )
}

```


##Boosted Logistic regression model
```{r, fig.height = 3, fig.height = 6}
registerDoParallel(detectCores())
ModelLogit <- train(
    x = feaTr, 
    y = resTr, 
    method = "LogitBoost", 
    preProcess = c("center","scale"), 
    trControl = ctrl,
    metric = "Kappa",
    tuneLength = 20
    )
stopImplicitCluster()
modelReport(model = ModelLogit)
```


##Linear Discriminant Analysis
```{r, fig.height = 3, fig.height = 6}
ModelLda <- train(
    x = feaTr, 
    y = resTr, 
    method = "lda", 
    preProcess = c("center","scale"), 
    metric = "Kappa",
    trControl = ctrl
    )
modelReport(model = ModelLda)
```


##K nearest Neighbours

```{r, fig.height = 3, fig.height = 6}
ModelKnn <- train(
    x = feaTr, 
    y = resTr, 
    method = "knn", 
    preProcess = c("center","scale"), 
    trControl = ctrl,
    metric = "Kappa",
    tuneLength = 10
    )
modelReport(model = ModelKnn)
```


##Naive Bayes

```{r, fig.height = 3, fig.height = 6}

nbGrid <- expand.grid(laplace = 1:3, usekernel = c(T,F), adjust = 1:3)
registerDoParallel(detectCores())
ModelNb <- train(
    x = feaTr, 
    y = resTr, 
    method = "naive_bayes", 
    preProcess = c("center","scale"), 
    metric = "Kappa",
    trControl = ctrl, 
    tuneGrid = nbGrid
    )
stopImplicitCluster()
modelReport(model = ModelNb)
```


##Bagged CART

```{r, fig.height = 3, fig.height = 6}
registerDoParallel(detectCores())
ModelBagCart <- train(
    x = feaTr, 
    y = resTr, 
    method = "treebag", 
    preProcess = c("center","scale"), 
    metric = "Kappa",
    trControl = ctrl
    )
stopImplicitCluster()
modelReport(model = ModelBagCart)
```


##Random Forest

```{r, fig.height = 3, fig.height = 6}
registerDoParallel(detectCores())
ModelRf <- train(
    x = feaTr, 
    y = resTr, 
    method = "rf", 
    preProcess = c("center","scale"), 
    trControl = ctrl,
    metric = "Kappa",
    tuneLength = 10
    )
stopImplicitCluster()
modelReport(model = ModelRf)
```


##Stochastic Gradient Boosting

```{r, fig.height = 3, fig.height = 6}
gbmGrid <- expand.grid(interaction.depth = c(3, 5),
                       n.trees = c(100, 300, 500),
                       shrinkage = c(0.01, 0.1),
                       n.minobsinnode = c(5, 10))
registerDoParallel(detectCores())
ModelGbm <- train(
    x = feaTr, 
    y = resTr, 
    method = "gbm", 
    preProcess = c("center","scale"), 
    trControl = ctrl,
    metric = "Kappa",
    tuneGrid = gbmGrid
    )
stopImplicitCluster()
modelReport(model = ModelGbm)
```


#Data summarizations

Comparison of All tested models

Information retrieval

```{r, fig.keep = "none"}
modelBestTune <- function(model = NULL)
{
  bestTuneParam <- as.data.frame(model$bestTune)
  resultTable <- model$results
  n <- ncol(bestTuneParam)
  bestTuneParam <- as.data.frame(bestTuneParam[,names(resultTable)[1:n]])
  
  bestTune <- rep(TRUE, nrow(resultTable))
  for (i in 1:n)
    bestTune <- bestTune & (resultTable[,i] == bestTuneParam[1,i])
  bestModel <- as.matrix(resultTable[bestTune, (n+1):ncol(resultTable)])
  row.names(bestModel) <- model$modelInfo$label
  
  reportBest <- modelReport(model)
  
  bestModel <- cbind(bestModel, 
                    t(as.matrix(reportBest[[3]]$overall[1:2])), 
                    t(as.matrix(reportBest[[3]]$byClass[1:2])),
                    t(as.matrix(reportBest[[4]]$overall[1:2])), 
                    t(as.matrix(reportBest[[4]]$byClass[1:2])),
                    auc = as.numeric(reportBest[[5]]$auc))
  colnames(bestModel)[5:8] <- paste("Tr", colnames(bestModel)[5:8], sep = "_")
  colnames(bestModel)[9:12] <- paste("Test", colnames(bestModel)[9:12], sep = "_")
  return(bestModel)
}


modelList <- mget(ls(pattern = "Model"))
modelSum <- do.call(rbind, lapply(X = modelList, FUN = modelBestTune))
modelSum <- data.frame(modelSum)
modelSum$modelName <- row.names(modelSum)

```
##Names and abbreviations of Models

```{r}
kable(cbind(Name = row.names(modelSum), Abbreviation = as.character(modelSum$modelName)))
```

##Accuracy in Cross Validation/Training/Test
```{r, fig.width= 10, fig.height= 5}
ggplot(data = modelSum) +
  geom_point(mapping = aes(x = modelName, y = Accuracy, color = "cv"), size = 2) +
  geom_errorbar(mapping = aes(x = modelName, 
              ymax = Accuracy + AccuracySD,
              ymin = Accuracy - AccuracySD,
              color = "cv"
              ), width = 0.3, show.legend = FALSE) +
  geom_point(mapping = aes(x = modelName, y = Tr_Accuracy, color = "Training"), size = 2) +
  geom_point(mapping = aes(x = modelName, y = Test_Accuracy, color = "Test"), size = 2) +
  geom_hline(yintercept = 282/345, linetype = 3 ) +
  xlab("Models") +
  ylab("Accuracy") + 
  coord_flip() +
  theme_bw()
```

##Kappa

```{r, fig.width= 10, fig.height= 5}
ggplot(data = modelSum) +
  geom_point(mapping = aes(x = modelName, y = Kappa, color = "cv"), size = 2) +
  geom_errorbar(mapping = aes(x = modelName, 
               ymax = Kappa + KappaSD,
              ymin = Kappa - KappaSD,
              color = "cv"
              ), width = 0.3) +
  geom_point(mapping = aes(x = modelName, y = Tr_Kappa, color = "Training"), size = 2) +
  geom_point(mapping = aes(x = modelName, y = Test_Kappa, color = "Test"), size = 2) +
  xlab("Models") +
  ylab("Kappa") + 
  coord_flip() +
  theme_bw()
```

##Sensitivity
```{r, fig.width= 10, fig.height= 5}
ggplot(data = modelSum) +
  geom_point(mapping = aes(x = modelName, y = Tr_Sensitivity, color = "Sensitivity Training"), size = 2) +
  geom_point(mapping = aes(x = modelName, y = Test_Sensitivity, color = "Sensitivity Test"), size = 2) +
  xlab("Models") +
  ylab("Sensitivity") + 
  coord_flip() +
  theme_bw()
```


##Specificity


```{r, fig.width= 10, fig.height= 5}
ggplot(data = modelSum) +
  geom_point(mapping = aes(x = modelName, y = Tr_Specificity, color = "Specificity Training"), size = 2) +
  geom_point(mapping = aes(x = modelName, y = Test_Specificity, color = "Specificity Test"), size = 2) +
  xlab("Models") +
  ylab("Specificity") + 
  coord_flip() +
  theme_bw()
```

##The Specificity-Sensitivity Plot

```{r, , fig.width= 7, fig.height= 5}
ggplot(data = modelSum) +
  geom_point(mapping = aes(x = 1 - Tr_Specificity, y = Tr_Sensitivity, color = "Training"), size = 2) +
  geom_point(mapping = aes(x = 1 - Test_Specificity, y = Test_Sensitivity, color = "Test"), size = 2) +
  geom_text(mapping = aes(x = 1 - Tr_Specificity, y = Tr_Sensitivity, color = "Training", label = modelName), size = 3, angle = 45, hjust = 0, show.legend = FALSE) +
  geom_text(mapping = aes(x = 1 - Test_Specificity, y = Test_Sensitivity, color = "Test", label = modelName), size = 3, angle = 45, hjust = 0, show.legend = FALSE) + 
  xlab("1 - Specificity") +
  ylab("Sensitivity") +
  coord_cartesian(xlim = c(0,1)) +
  geom_abline(slope = 1, linetype = 3) +
  theme_bw(
  )
```

##Area under curve


```{r, fig.width= 10, fig.height= 5}
ggplot(data = modelSum) +
  geom_point(mapping = aes(x = modelName, y = auc), size = 2) +
  xlab("Models") +
  ylab("Area Under ROC curve") + 
  coord_flip() +
  theme_bw()
```

#Conclusion
None of the performance of the above models is satisfactory, meaning that ftrs are weak predictors of stable and unstable CVD. Gbm performs better than other models, with best Area Under ROC curve and best sensitivity in the test dataset. 

#Session info
```{r}
save.image(file = "./Rdata/Mir_Stable_NonStable_new.RData")
sessionInfo()
```
