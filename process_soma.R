suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(SomaScan.db))


raw_df <- readRDS("D:/AZ/CSF_dementia/proc/proteomics")


# 1-RID, 7010-7014 are "APGEN1"  "APGEN2"  "apoe"    "status"  "status2" respectively
cols_to_remove <- c("RID", "status", "status2", "apoe", "APGEN2", "APGEN1" )
ind_remove <- which(colnames(raw_df) %in% cols_to_remove)
colnames(raw_df)[-ind_remove] <- gsub("[.]", "-", gsub("X","", colnames(raw_df)[-ind_remove]))

# Annotation PROBEID -> Symbol
ss <- select(SomaScan.db, keys = colnames(raw_df)[-ind_remove], columns = c("UNIPROT","SYMBOL"))

ss <- ss %>%
  filter(!duplicated(PROBEID, fromLast=F))%>%
  filter(!duplicated(UNIPROT, fromLast=F))

colnames(raw_df)[-ind_remove] <- ss[match(colnames(raw_df)[-ind_remove], ss$PROBEID), "SYMBOL"]

#926 didn't map
sum(is.na(colnames(raw_df)))

raw_df<- raw_df[,!is.na(colnames(raw_df))]




#PCA
df_PCA <- raw_df[raw_df$status2 == "AZ+" | raw_df$status2 == "AZ-",]
cols_to_remove <- c("RID", "status", "status2", "apoe", "APGEN2", "APGEN1" )
ind_remove <- which(colnames(raw_df) %in% cols_to_remove)
pca <- prcomp(na.omit(scale(raw_df[,-ind_remove])))


dataPCA <- data.frame(PC1 = pca$x[,1],PC2 = pca$x[,2], Condition = as.factor(raw_df$status2) )

qplot(PC1, PC2, data = dataPCA, color = Condition) + geom_point(size=3) + theme_classic() + 
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size=14),
        legend.text = element_text(size=12))

saveRDS(raw_df, "D:/AZ/CSF_dementia/proc/dataset_proc")
