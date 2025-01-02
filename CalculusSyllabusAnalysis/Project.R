rm(list = ls())
dev.off()
gc()
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
stopCluster(clust) # Free cores

library(ggfortify)
library(ggplot2)
library(caret)
library(quanteda)
library(doSNOW)
library(randomForest)
library(lsa)
library(readxl)
library(glmnet)
library(e1071)
library(MASS)
library(class)
library(naivebayes)
library(rpart.plot)
library(gridExtra)

clust = makeCluster(8, type = "SOCK")
registerDoSNOW(clust)  # Train in parallel

set.seed(200)

data = read_excel("dataset.xlsx")
data$Label = as.factor(data$Label) # Convert to factors

data$SUNY = as.factor(data$SUNY)# Convert to factors
data$State = as.factor(data$State)# Convert to factors
data$Admission = as.numeric(data$Admission)# Convert to number
data$CUNY = as.factor(data$CUNY)# Convert to factors
data$Length = nchar(data$Syllabus) # Text length as new feature
prop.table(table(data$Label)) # Show proportion
data$Credits = as.numeric(data$Credits)# Convert to number
data$Public = as.factor(data$Public)# Convert to factors
data$Community = as.factor(data$Community)

# Visualize distribution
# Good

x.1 = ggplot(data, aes(x = Length, fill = Label)) + theme_bw() +
  geom_density(position = "fill") +
  labs(y = "Density of Schools", x = "Text Length", 
       title = "Distribution of Text Length")

x.2 = ggplot(data, aes(x = Credits, fill = Label)) + theme_bw() +
  geom_bar(position = "fill") +
  labs(y = "Number of Schools", x = "Credit Hours", 
       title = "Distribution of Credit Hours")

x.3 = ggplot(data, aes(x = Admission, fill = Label)) + theme_bw() +
  geom_density(position = "fill") +
  labs(y = "Density of Schools", x = "Admission Rate", 
       title = "Distribution of Admission Rate")

x.4 = ggplot(data, aes(x = State, fill = Label)) + theme_bw() +
  geom_bar(position = "fill") +
  labs(y = "Number of Schools", x = "In State School", 
       title = "Distribution of In/Out State School")

# Bad

y.1 = ggplot(data, aes(x = CUNY, fill = Label)) + theme_bw() +
  geom_bar(position = "fill") +
  labs(y = "Number of Schools", x = "CUNY School", 
       title = "Distribution of CUNY/None CUNY")

y.2 = ggplot(data, aes(x = Public, fill = Label)) + theme_bw() +
  geom_bar(position = "fill") +
  labs(y = "Number of Schools", x = "Public School", 
       title = "Distribution of Public/Private")

y.3 = ggplot(data, aes(x = Community, fill = Label)) + theme_bw() +
  geom_bar(position = "fill") +
  labs(y = "Number of Schools", x = "Community School", 
       title = "Distribution of Community/None Community")

y.4 = ggplot(data, aes(x = SUNY, fill = Label)) + theme_bw() +
  geom_bar(position = "fill") +
  labs(y = "Number of Schools", x = "SUNY School", 
       title = "Distribution of SUNY/None SUNY")

grid.arrange(x.1, x.2, y.1, y.2, ncol=2, nrow = 2)

# Stratified split
split = createDataPartition(data$Label, times = 1, p = 0.7, list = FALSE)
train = data[split, ]
test = data[-split, ]
prop.table(table(train$Label))
prop.table(table(test$Label))

total.tokens = tokens(data$Syllabus, what = "word", remove_numbers = TRUE,
                      remove_punct = TRUE, remove_symbols = TRUE, remove_url = TRUE)
# Tokenization
total.tokens = tokens_tolower(total.tokens) # Lowercase
total.tokens = tokens_select(total.tokens, stopwords(), selection = "remove")
# Remove stop-word
total.tokens = tokens_wordstem(total.tokens, language = "english") # Stemming
total.tokens = tokens_ngrams(total.tokens, n = 1:3) # Add bigram and trigram
# For word-ordering 
total.frequency = dfm(total.tokens, tolower = FALSE)
colnames(total.frequency) = make.names(colnames(total.frequency))
# Make valid column names
total.frequency
total.column = colnames(total.frequency) # Total column

train.tokens = tokens(train$Syllabus, what = "word", remove_numbers = TRUE,
                     remove_punct = TRUE, remove_symbols = TRUE, remove_url = TRUE)
train.tokens = tokens_tolower(train.tokens)
train.tokens = tokens_select(train.tokens, stopwords(), selection = "remove")
train.tokens = tokens_wordstem(train.tokens, language = "english")
train.tokens = tokens_ngrams(train.tokens, n = 1:3)
train.frequency = dfm(train.tokens, tolower = FALSE)
colnames(train.frequency) = make.names(colnames(train.frequency))
# Make valid column names
train.frequency
train.frequency = dfm_match(train.frequency, total.column) # Match
train.frequency
train.matrix = as.matrix(train.frequency) # Document-frequency matrix

term.frequency = function(row){
  row / sum(row) # Normalization
} # Term frequency proportion
inverse.doc.frequency = function(col){
  size = length(col)
  count = length(which(col > 0))
  log(size / count) # Penalization
} # Inverse document-frequency
tf.idf = function(tf, idf){
  tf * idf
}

train.tf = apply(train.matrix, 1, term.frequency)
# Loop through rows
# Return transpose
# Term-frequency matrix
train.idf = apply(train.matrix, 2, inverse.doc.frequency)
# Loop through columns
# Still transpose
# Term-frequency matrix
train.tf.idf = apply(train.tf, 2, tf.idf, idf = train.idf)
# Loop though columns because transpose
# Still transpose
# Term-frequency matrix
train.tf.idf = t(train.tf.idf)
# Transpose back
# Document-frequency matrix

train.incomplete.column = which(!complete.cases(t(train.tf.idf)))
train.tf.idf[, train.incomplete.column] = rep(0.0, nrow(train.tf.idf))
# Fix incomplete cases

cv.folds = createMultiFolds(train$Label, k = 10, times = 3)
# 10-fold cross-validation repeated 3 times
cv.control = trainControl(method = "repeatedcv", number = 10, repeats = 3, 
                          index = cv.folds)
lambda.range = 10^seq(-1, -3, length = 100)

train.lasso = cv.glmnet(train.tf.idf, train$Label, family = "binomial",
                        alpha = 1, lambda = lambda.range, type.measure = "class")
train.lasso.column = train.tf.idf[, (which(coef(train.lasso, s = "lambda.min")
                                                 != 0) - 1)[-1]]
train.lasso.column = as.data.frame(train.lasso.column)
length((which(coef(train.lasso, s = "lambda.min")
              != 0) - 1)[-1])
train.tf.idf.lasso = data.frame(train.lasso.column,
                                  State = train$State, Admission = train$Admission,
                                  Length = train$Length, Credits = train$Credits)

test.tokens = tokens(test$Syllabus, what = "word", remove_numbers = TRUE,
                      remove_punct = TRUE, remove_symbols = TRUE,
                      split_hyphens = TRUE, remove_url = TRUE) # Tokenization
test.tokens = tokens_tolower(test.tokens) # Lowercase
test.tokens = tokens_select(test.tokens, stopwords(), selection = "remove")
# Remove stop-word
test.tokens = tokens_wordstem(test.tokens, language = "english") # Stemming

test.tokens = tokens_ngrams(test.tokens, n = 1:3) # Add bigrams
# For word-ordering 

test.frequency = dfm(test.tokens, tolower = FALSE)
# Document-feature frequency
test.frequency = dfm_match(test.frequency, total.column) # Match

test.frequency = dfm_match(test.frequency, total.column) # Match
test.frequency
test.matrix = as.matrix(test.frequency) # Document-frequency matrix

test.tf = apply(test.matrix, 1, term.frequency)
# Loop through rows
# Return transpose
# Term-frequency matrix
test.tf.idf = apply(test.tf, 2, tf.idf, idf = train.idf)
# Loop though columns because transpose
# Still transpose
# Term-frequency matrix
test.tf.idf = t(test.tf.idf)
# Transpose back
# Document-frequency matrix

test.incomplete.column = which(!complete.cases(t(test.tf.idf)))
test.tf.idf[, test.incomplete.column] = rep(0.0, nrow(test.tf.idf))
# Fix incomplete cases

test.lasso.column = test.tf.idf[, (which(coef(train.lasso, s = "lambda.min")
                                         != 0) - 1)[-1]]
test.lasso.column = as.data.frame(test.lasso.column)
test.tf.idf.lasso = data.frame(test.lasso.column,
                                State = test$State, Admission = test$Admission,
                                Length = test$Length, Credits = test$Credits)

ensemble.learning = function(train, test, train.label, test.label){
  
  train.data = data.frame(Label = train.label, train)
  k.accuracy = vector()
  
  # K-nearest neighbors
  for (i in 1:as.integer(sqrt(length(train.label)))){
    k.result = knn(train, test, train.label, k = i)
    k.accuracy = c(k.accuracy, mean(k.result == test.label))
  }
  print(as.integer(sqrt(length(train.label)))) # Range of K
  k.best = min(which(k.accuracy == max(k.accuracy)))
  print(k.best)
  k.result = knn(train, test, train.label, k = k.best)
  # Retrain
  print(mean(k.result == test.label))
  print(confusionMatrix(k.result, test.label))
  
  # Logistic regression
  logistic.fit = glm(Label~., data = train.data, family = "binomial")
  logistic.result = predict(logistic.fit, newdata = test, type = "response")
  logistic.result <- ifelse(logistic.result < 0.5, "223/224", "224/225")
  print(mean(logistic.result == test.label))
  logistic.result = as.factor(logistic.result)
  print(confusionMatrix(logistic.result, test.label))
  
  # Naive bayes
  naive.bayes = naive_bayes(Label~., data = train.data)
  bayes.result = predict(naive.bayes, newdata = test, type = "class")
  print(mean(bayes.result == test.label))
  print(confusionMatrix(bayes.result, test.label))
  
  # Neural network
  neural.network <- train(Label~., data = train.data, method = "nnet",
                          trControl = cv.control, tuneLength = 10, maxit = 100)
  print(neural.network)
  network.result = predict(neural.network, test)
  print(mean(network.result == test.label))
  print(confusionMatrix(network.result, test.label))
  
  # Single decision tree
  r.part <- train(Label~., data = train.data, method = "rpart",
                  trControl = cv.control, tuneLength = 20)
  print(r.part)
  decision.result = predict(r.part, test)
  print(mean(decision.result == test.label))
  print(r.part$finalModel$variable.importance)
  rpart.plot(r.part$finalModel)
  print(confusionMatrix(decision.result, test.label))
  
  # Random forests
  random.forests = train(Label~., data = train.data, method = "rf",
                         trControl = cv.control, tuneLength = 20,
                         importance = TRUE)
  print(random.forests)
  forests.result = predict(random.forests, test)
  print(mean(forests.result == test.label))
  print(confusionMatrix(forests.result, test.label))
  varImpPlot(random.forests$finalModel, main = "Importance of Word(s)")
  
  ensemble = data.frame(Nearest.Neighbors = k.result, Logistic = logistic.result,
                        Naive.Bayes = bayes.result,
                        Neural.Network = network.result,
                        Decision.Tree = decision.result,
                        Random.Forests = forests.result)
}

find.most = function(data){
  result = vector()
  row = dim(data)[1] # Find row number
  for (i in 1:row){
    temp = data[i, ]
    if (length(which(temp == "223/224")) > length(which(temp == "224/225")))
      result = c(result, "223/224")
    else
      result = c(result, "224/225")
  }
  return(result)
}


# Training
train.ensemble = ensemble.learning(train.lasso.column, train.lasso.column, train$Label,
                                   train$Label)
train.final = find.most(train.ensemble)
train.final = as.factor(train.final)
mean(train.final == train$Label)
confusionMatrix(train.final, train$Label)
# Confusion matrix


# Testing
test.ensemble = ensemble.learning(train.lasso.column, test.lasso.column, train$Label,
                                  test$Label)

test.final = find.most(test.ensemble)
test.final = as.factor(test.final)
mean(test.final == test$Label)
confusionMatrix(test.final, test$Label)
# Confusion matrix


# Training feature engineer
train.ensemble = ensemble.learning(train.tf.idf.lasso, train.tf.idf.lasso, train$Label,
                             train$Label)
train.final = find.most(train.ensemble)
train.final = as.factor(train.final)
mean(train.final == train$Label)
confusionMatrix(train.final, train$Label)
# Confusion matrix



# Testing feature engineer
test.ensemble = ensemble.learning(train.tf.idf.lasso, test.tf.idf.lasso, train$Label,
                             test$Label)

test.final = find.most(test.ensemble)
test.final = as.factor(test.final)
mean(test.final == test$Label)
confusionMatrix(test.final, test$Label)
# Confusion matrix

index = which(test.final != test$Label)

data.frame(School = test$School[index],
           Label = test$Label[index],
           Prediction = test.final[index])
# Label confilt prediction






total.frequency = as.matrix(total.frequency)


x = prcomp(total.frequency)
summary(x)
autoplot(x, data = data, colour = 'Label')
