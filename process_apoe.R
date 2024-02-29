



apoeres <- read.csv("D:/AZ/proteomics/raw/apoe.csv", stringsAsFactors = F)


apoeres <- apoeres %>% 
  mutate(apoe = if_else(APGEN1 == 4 | APGEN2 == 4, 1,0) ) %>% 
  select(RID, APGEN1, APGEN2, apoe)

saveRDS(apoeres, "D:/AZ/proteomics/proc/apoe")


outcome <- read.csv("D:/AZ/aibl/adni/proc_outcome.csv", stringsAsFactors = F)


outcome <- outcome %>% 
  rename(status = DXCURRENT) %>% 
  filter(status != 2) %>%
  mutate(status = case_when(status == 1 ~ "control",
                            status == 3 ~ "AZ")) %>% 
  mutate(status = factor(status, levels = c("AZ", "control")))

saveRDS(outcome, "D:/AZ/proteomics/proc/outcome")


proteomics <- read.csv("D:/AZ/proteomics/raw/prot.csv", stringsAsFactors = F)

subst.mean <- function(x) replace(x, is.na(x), median(x, na.rm = TRUE))

prot_df <- proteomics %>%
  select(-c(EXAMDATE:PlateId)) %>%
  inner_join(apoeres, by="RID") %>%
  inner_join(outcome, by="RID") %>%
  mutate(status2 = case_when(apoe == 1 & status == "AZ" ~ "AZ+",
                              apoe == 0 & status == "AZ" ~ "AZ-",
                              apoe == 1 & status == "control" ~ "control+",
                              apoe == 0 & status == "control" ~ "control-",
                              .default = NA) ) %>%
  group_by(status2) %>%
  mutate(across(X10000.28:X9999.1, subst.mean))







saveRDS(prot_df, "D:/AZ/proteomics/proc/proteomics")

