suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(SomaScan.db))


raw_df <- readRDS("D:/AZ/proteomics/proc/proteomics")






df_no_control <- raw_df[raw_df$status == "AZ",]
df_no_control <- raw_df[raw_df$status == "control",]
df_no_control <- raw_df[raw_df$status2 == "AZ+" | raw_df$status2 == "control+",]


colnames(df_no_control)[-c(1,7010:7014)] <- gsub("[.]", "-", gsub("X","", colnames(df_no_control)[-c(1,7010:7014)]))

ss <- select(SomaScan.db, keys = colnames(df_no_control)[-c(1,7010:7014)], columns = c("UNIPROT","SYMBOL"))

ss <- ss %>%
  filter(!duplicated(PROBEID, fromLast=F))%>%
  filter(!duplicated(UNIPROT, fromLast=F))

colnames(df_no_control)[-c(1,7010:7014)] <- ss[match(colnames(df_no_control)[-c(1,7010:7014)], ss$PROBEID), "SYMBOL"]


df_no_control <- df_no_control[,!is.na(colnames(df_no_control))]


cols_to_remove <- c("RID", "status", "status2", "apoe", "APGEN2", "APGEN1" )

ind_remove <- which(colnames(df_no_control) %in% cols_to_remove)


#PCA
pca <- prcomp(na.omit(scale(df_no_control[,rem])))


dataPCA <- data.frame(PC1 = pca$x[,1],PC2 = pca$x[,2], Condition = as.factor(df_no_control$status2) )

qplot(PC1, PC2, data = dataPCA, color = Condition) + geom_point(size=3) + theme_classic() + 
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size=14),
        legend.text = element_text(size=12))



pca.res <- FactoMineR::PCA(na.omit(scale(df_no_control[,-ind_remove])), graph = F)

contr <- data.frame(pca.res$var$cor)
#contr <- data.frame(pca.res$var$contrib)
contr_sort <- contr[order(abs(contr$Dim.1), decreasing = T),]
contr_sort_dataset1 <- contr_sort[abs(contr$Dim.1) >= 0.80, ]
contr_sort_dataset1 <- contr_sort[1:1000, ]
write.csv(row.names(contr_sort_dataset1), "D:/AZ/proteomics/proc/Dim1_proteins_1000.csv")


contr_sort_dataset2 <- contr_sort[c(1000:1050),]
contr_sort_dataset2$gene <- fdata_gse97760[match(rownames(contr_sort_dataset2), rownames(fdata_gse97760)), "GENE_SYMBOL"]


contr_sort <- contr[order(abs(contr$Dim.2), decreasing = T),]
contr_sort_dataset1 <- contr_sort[abs(contr$Dim.2) >= 0.80, ]
contr_sort_dataset1 <- contr_sort[1:1000, ]
write.csv(row.names(contr_sort_dataset1), "D:/AZ/proteomics/proc/Dim2_proteins_1000.csv")










li <- IV$Variable[1:30]
li <- gsub(" ", "", readClipboard())
li <- c("OTULIN", "CCL25", "CHCHD7", "S100A13", "LRRN1", "SPC25")
li <- c("HYKK")
li <- c("CXCL9","CXCL13","CXCL11","CXCL2","CXCL3","CXCL12","CXCL16","CXCL5","CXCL8","CXCL6","CXCL14","CXCL17")
li <- c("ATE1", "HYKK", "RDH16", "OTULIN", "CHCHD7", "RND1", "HAX1",
        "CTF1", "SNRPF", "LMO4", "MTMR7", "FOXO1", "CHFR", "ARL2", "AGR2", "PLEKHO2",
        "CABP7", "VPS29", "WNT10B", "KCTD2", "NFIA", "NUDT7", "TP53I11", "RNF215", "IRF6", "MAPK8")

df_ML <- df_no_control[,colnames(df_no_control) %in% c(li, cols_to_remove)]
df_ML <- df_no_control[,colnames(df_no_control) %in% c(IV_proc$Variable,cols_to_remove)]


colnames(df_ML) <-make.names(colnames(df_ML))


# ML
set.seed(42)
#rm(.Random.seed, envir=globalenv())
ind=sample(2, nrow(df_no_control), prob = c(.7, .3), replace = T)
ind_remove <- which(colnames(df_ML) %in% c("RID","status",  "status2", "apoe", "APGEN2", "APGEN1" ))

training <- as.data.frame(scale(df_ML[ind==1,-ind_remove]))
test <- as.data.frame(scale(df_ML[ind==2,-ind_remove]))

training$target <- df_ML[ind==1,"status2", drop=T]
training$target <- ifelse(training$target == "AZ+", "AZ", "Ctr")
test$target <- as.factor(df_ML[ind==2,"status2", drop=T])
test$target <- ifelse(test$target == "AZ+", "AZ", "Ctr")


training$target <- as.factor(df_ML[ind==1,"status2", drop=T])
training$target <- ifelse(training$target == "control+", "control1", "control2")
test$target <- as.factor(df_ML[ind==2,"status2", drop=T])
test$target <- ifelse(test$target == "control+", "control1", "control2")


training$target <- as.factor(df_ML[ind==1,"status2", drop=T])
training$target <- ifelse(training$target == "AZ+", "AZ", "control")
test$target <- as.factor(df_ML[ind==2,"status2", drop=T])
test$target <- ifelse(test$target == "AZ+", "AZ", "control")










control <- trainControl(method="repeatedcv", number=5, repeats=5, 
                        summaryFunction = twoClassSummary,
                        classProbs = TRUE,
                        sampling = "up")

cl <- makePSOCKcluster(15)
registerDoParallel(cl)


my_glm <- train(target ~ .,
                data = training,
                method = "rpart",
                metric = "ROC",
                #tuneLength = 100,
                #tuneGrid = grid,
                trControl = control)


predict_unseen <-predict(my_glm, test)
table_mat <- table(test$target, predict_unseen)
cm <- confusionMatrix(table_mat)
cm


plot(my_glm)
my_glm$bestTune
varImp(my_glm)

my_glm$finalModel

cm$byClass


first_rem <- c("OTULIN", "CCL25", "CHCHD7", "S100A13", "LRRN1", "SPC25")

rem <- c("TBCA", "CTF1")
rem <- c("TBCA", "CTF1", "ATE1", "HAX1")
rem <- c("TBCA", "CTF1", "ATE1", "HAX1", "SNRPF", "PLEKHO2", "AGR2")
rem <- c("TBCA", "CTF1", "ATE1", "HAX1", "SNRPF", "PLEKHO2", "AGR2", "HYKK", "RNF215")
rem <- c("TBCA", "CTF1", "ATE1", "HAX1", "SNRPF", "PLEKHO2", "AGR2", "HYKK", "RNF215", "VSNL1", "RND1")
rem <- c("TBCA", "CTF1", "ATE1", "HAX1", "SNRPF", "PLEKHO2", "AGR2", "HYKK", "RNF215", "VSNL1", "RND1", "RDH16", "PGK1", "MAPK8")

# Specificty in CART is < 0.7
rem <- c("TBCA", "CTF1", "ATE1", "HAX1", "SNRPF", "PLEKHO2", "AGR2", "HYKK", "RNF215", "VSNL1", "RND1", 
         "RDH16", "PGK1", "MAPK8", "MTMR7", "FOXO1", "WNT10B", "LMO4", "KCTD2")
rem <- c("TBCA", "CTF1", "ATE1", "HAX1", "SNRPF", "PLEKHO2", "AGR2", "HYKK", "RNF215", "VSNL1", "RND1", 
         "RDH16", "PGK1", "MAPK8", "MTMR7", "FOXO1", "WNT10B", "LMO4", "KCTD2", "CHFR", "CABP7", "TP53I11", 'ARL2', "VPS29", "IRF6", "NFIA")


training <- training[,!colnames(training) %in% rem]
test <- test[,!colnames(test) %in% rem]

iv_set <- IV_proc[IV_proc$IV > 1.1, "Variable"]



training <- training[,!colnames(training) %in% c(first_rem, iv_set)]
test <- test[,!colnames(test) %in% c(first_rem, iv_set)]

for(i in 2:7010){
  ts <- cor.test(df_no_control[df_no_control$status2 == "AZ+", i, drop=T],
              df_no_control[df_no_control$status2 == "AZ-", i, drop=T])
  if(ts$p.value <= 1e-20){
    cat(colnames(df_no_control)[i], ts$estimate, ts$p.value, "\n")
  }
  
}


# Information Value
fs_training <- training
fs_training$target <- ifelse(fs_training$target == "AZ",1,0 )

rr <- Information::create_infotables(fs_training, y = "target")
IV <- rr$Summary
IV_proc <- IV[IV$IV > 0.3,]



training <- training[,colnames(training) %in% c(IV_proc$Variable, "target")]
test <- test[,colnames(test) %in% c(IV_proc$Variable, "target")]

for(feat in IV_proc$Variable){

  my_glm <- train(target ~ .,
                  data = training,
                  method = "rpart",
                  metric = "ROC",
                  #tuneLength = 100,
                  #tuneGrid = grid,
                  trControl = control)
  
  
  predict_unseen <-predict(my_glm, test)
  table_mat <- table(test$target, predict_unseen)
  cm <- confusionMatrix(table_mat)
  if(cm$byClass[1] > 0.75 & cm$byClass[2] > 0.75){
    cat(feat,"\n")
  }
  
}


ggplot(df_no_control, aes(x = status2, y = LRRN1))+
  geom_boxplot()


for(feat in IV_proc$Variable){
  training_s <- training[,c(feat, "target")]
  test_s <- test[,c(feat, "target")]
  my_glm <- train(target ~ .,
                  data = training_s,
                  method = "rpart",
                  metric = "ROC",
                  #tuneLength = 100,
                  #tuneGrid = grid,
                  trControl = control)
  
  
  predict_unseen <-predict(my_glm, test_s)
  table_mat <- table(test_s$target, predict_unseen)
  cm <- confusionMatrix(table_mat)
  if(cm$byClass[1] > 0.75 & cm$byClass[2] > 0.75){
    cat(feat,"\n")
  }
  
}

iv_set <- IV_proc$Variable
iv_set <- iv_set[!iv_set %in% c(li, "TBCA")]
iv_set <- iv_set[1:50]

write.csv(iv_set, "D:/AZ/proteomics/proc/top_50_list.csv")
