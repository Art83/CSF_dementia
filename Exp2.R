library(caret)
library(doParallel)

n_cores <- detectCores()
cluster <- makeCluster(n_cores - 1)
registerDoParallel(cluster)

df <- readRDS("D:/AZ/CSF_dementia/proc/dataset_proc")

df <- df[df$status2 == "AZ+" | df$status2 == "control+",]

prots1 <- read.csv("D:/AZ/CSF_dementia/proc/top_50_list.csv", row.names = 1)[,1]
prots2 <- c("TBCA", "OTULIN", "CCL25", "CHCHD7", "S100A13", "LRRN1", "SPC25" )


#Split to training and testing
set.seed(42)
#rm(.Random.seed, envir=globalenv())
ind=sample(2, nrow(df), prob = c(.6, .4), replace = T)
ind_remove <- which(colnames(df) %in% c("RID","status",  "status2", "apoe", "APGEN2", "APGEN1" ))

training <- as.data.frame(scale(df[ind==1,-ind_remove]))
test <- as.data.frame(scale(df[ind==2,-ind_remove]))

training$target <- df[ind==1,"status2", drop=T]
training$target <- ifelse(training$target == "AZ+", "AZ", "Ctr")
test$target <- as.factor(df[ind==2,"status2", drop=T])
test$target <- ifelse(test$target == "AZ+", "AZ", "Ctr")


control <- trainControl(method="repeatedcv", number=3, repeats=3, 
                        summaryFunction = twoClassSummary,
                        classProbs = TRUE,
                        sampling = "up")


training <- training[,colnames(training) %in% c(prots1, "target")]
test <- test[,colnames(test) %in% c(prots1, "target")]


my_glm <- train(target ~ .,
                data = training,
                method = "rpart",
                metric = "ROC",
                tuneLength = 1000,
                #tuneGrid = grid,
                trControl = control)


predict_unseen <-predict(my_glm, test)
table_mat <- table(test$target, predict_unseen)
cm <- confusionMatrix(table_mat)
cm


predict_unseen <-predict(my_glm, test, type = 'prob')
vec <- predict_unseen[,1]
pROC::roc(test$target, vec, plot = T,print.auc=F,lwd =1,main="", percent = F, asp=NA)
