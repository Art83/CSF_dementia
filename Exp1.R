library(caret)
library(doParallel)

n_cores <- detectCores()
cluster <- makeCluster(n_cores - 1)
registerDoParallel(cluster)

df <- readRDS("D:/AZ/CSF_dementia/proc/dataset_proc")

df <- df[df$status2 == "AZ+" | df$status2 == "AZ-",]

prots1 <- read.csv("D:/AZ/CSF_dementia/proc/top_50_list.csv", row.names = 1)[,1]
prots2 <- c("TBCA", "OTULIN", "CCL25", "CHCHD7", "S100A13", "LRRN1", "SPC25" )

#Split to training and testing
set.seed(42)
#rm(.Random.seed, envir=globalenv())
ind=sample(2, nrow(df), prob = c(.7, .3), replace = T)
ind_remove <- which(colnames(df) %in% c("RID","status",  "status2", "apoe", "APGEN2", "APGEN1" ))

training <- as.data.frame(scale(df[ind==1,-ind_remove]))
test <- as.data.frame(scale(df[ind==2,-ind_remove]))

training$target <- df[ind==1,"status2", drop=T]
training$target <- ifelse(training$target == "AZ+", "AZ", "Ctr")
test$target <- as.factor(df[ind==2,"status2", drop=T])
test$target <- ifelse(test$target == "AZ+", "AZ", "Ctr")


# Information Value only on training
fs_training <- training
fs_training$target <- ifelse(fs_training$target == "AZ",1,0 )

rr <- Information::create_infotables(fs_training, y = "target")
IV <- rr$Summary
IV_proc <- IV[IV$IV > 0.3,] # only 1534 should be here at this point



# Trimming the datasets
training <- training[,colnames(training) %in% c(IV_proc$Variable, "target")]
test <- test[,colnames(test) %in% c(IV_proc$Variable, "target")]

control <- trainControl(method="repeatedcv", number=3, repeats=3, 
                        summaryFunction = twoClassSummary,
                        classProbs = TRUE,
                        sampling = "up")

colnames(training) <- make.names(colnames(training))
colnames(test) <- make.names(colnames(test))
IV_proc$Variable <- make.names(IV_proc$Variable)
for(feat in IV_proc$Variable){
  training_s <- training[,c(feat, "target")]
  test_s <- test[,c(feat, "target")]
  my_glm <- train(target ~ .,
                  data = training_s,
                  method = "rpart",
                  metric = "ROC",
                  tuneLength = 100,
                  trControl = control)
  
  
  predict_unseen <-predict(my_glm, test_s)
  table_mat <- table(test_s$target, predict_unseen)
  cm <- confusionMatrix(table_mat)
  if(cm$byClass[1] > 0.75 & cm$byClass[2] > 0.75){
    cat(feat,"\n")
  }
  
}


#TBCA
#OTULIN
#CCL25 
#CHCHD7 
#S100A13 
#LRRN1 
#SPC25 


iv_set <- IV_proc$Variable
iv_set <- iv_set[1:50]
write.csv(iv_set, "D:/AZ/proteomics/proc/top_50_list.csv")



# part B

training <- training[,colnames(training) %in% c(prots2, "target")]
test <- test[,colnames(test) %in% c(prots2, "target")]
control <- trainControl(method="repeatedcv", number=3, repeats=3, 
                        summaryFunction = twoClassSummary,
                        classProbs = TRUE,
                        sampling = "up")

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
