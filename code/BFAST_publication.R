# 1. Box plot for cTF vs categorical Baseline Characteristics
library("rice")
library("ggplot2")
library("ggpubr")
library("dplyr")
library("tidyr")
library("survival")
library("survminer")
library('readxl')

# 1.0 Import sample (OSPL) data 
OSPL <- read.csv("~/ctdna/data_cleaning_20221219/cohortA/bfast_oSpl_cohortA.withDBIannotations.csv"
                 ,stringsAsFactors = TRUE)
OAPL <- read.csv("~/ctdna/data_cleaning_20221219/cohortA/bfast_oapl_cohortA.withDBIannotations.clean.csv"
                 ,stringsAsFactors = TRUE)
cTF <- OSPL %>%
  dplyr::select(PT, VISIT, cTF, assay) %>%
  dplyr::filter(VISIT %in% 'SCRN') %>%
  dplyr::arrange(PT)

# 1.1 ECOG
adsub <- rice_read("root/clinical_studies/RO5424802/CDT30119/BO29554/data_analysis/Cohort_A_BioPub2023/prod/outdata_vad/adsub.sas7bdat")
df <- adsub %>%
  select(SUBJID, PARAMCD, AVALC) %>%
  filter(PARAMCD %in% 'BECOG') %>%
  mutate(PT = as.numeric(sub(".*-","",SUBJID)),
         ECOG = as.factor(AVALC)) %>%
  arrange(PT) %>%
  dplyr::left_join(cTF, by = 'PT')
kruskal.test(log(df$cTF) ~ df$ECOG)
ggplot(df, aes(x = ECOG, y = log(cTF))) +
  geom_boxplot() +
  labs(x = "Baseline ECOG", y = "log(cTF)", title = "log(cTF) vs Baseline ECOG", subtitle = "Kruskal-Wallis, p = 0.007374") +
  theme_bw() #600x500

# 1.2 SLD by IRF
adtr <- rice_read("root/clinical_studies/RO5424802/CDT30119/BO29554/data_analysis/Cohort_A_BioPub2023/prod/outdata_vad/adtr.sas7bdat")
SLD <- adtr %>%
  filter(ABLFL %in% 'Y') %>%
  select(SUBJID, PARAMCD, AVALC) %>%
  filter(PARAMCD %in% 'SUMLDIAM') %>%
  mutate(PT = as.numeric(sub(".*-","",SUBJID)),
         AVALC = as.numeric(AVALC)) %>%
  arrange(PT) %>%
  dplyr::left_join(cTF, by = 'PT') %>%
  mutate(logcTF = log(cTF))
ggscatter(SLD, x = "AVALC", y = "logcTF", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Baseline Sum of Longest Diameter by IRF", ylab = "log(cTF)") +
          theme_bw() #600x500
cor.test(SLD$AVALC, 
         SLD$logcTF,
         method = "pearson")

# 1.3 Sum of Diameters by Investigator
SLD <- adtr %>%
  filter(ABLFL %in% 'Y') %>%
  select(SUBJID, PARAMCD, AVALC) %>%
  filter(PARAMCD %in% 'SD') %>%
  mutate(PT = as.numeric(sub(".*-","",SUBJID)),
         AVALC = as.numeric(AVALC)) %>%
  arrange(PT) %>%
  dplyr::left_join(cTF, by = 'PT') %>%
  mutate(logcTF = log(cTF))
ggscatter(SLD, y = "logcTF", x = "AVALC", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          ylab = "log(cTF)", xlab = "Baseline Sum of Diameters by Investigator") +
          theme_bw() #600x500
cor.test(SLD$AVALC, 
         log(SLD$cTF),
         method = "pearson")

# 1.4 AGE
AGE <- adsub %>%
  select(SUBJID, AAGE, PARAMCD) %>%
  filter(PARAMCD %in% 'BECOG') %>%
  mutate(PT = as.numeric(sub(".*-","",SUBJID)),
         AAGE = as.numeric(AAGE)) %>%
  arrange(PT) %>%
  dplyr::left_join(cTF, by = 'PT') %>%
  mutate(logcTF = log(cTF))
ggscatter(AGE, y = "logcTF", x = "AAGE", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          ylab = "log(cTF)", xlab = "Baseline Age")
cor.test(AGE$AAGE, 
         AGE$logcTF,
         method = "pearson")

# 1.5 BMI
df <- adsub %>%
  select(SUBJID, PARAMCD, AVAL) %>%
  filter(PARAMCD %in% 'BMISI') %>%
  mutate(PT = as.numeric(sub(".*-","",SUBJID))) %>%
  arrange(PT) %>%
  dplyr::left_join(cTF, by = 'PT')  %>%
  mutate(logcTF = log(cTF))
ggscatter(df, y = "logcTF", x = "AVAL", 
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "pearson",
          ylab = "log(cTF)", xlab = "BMI")

# 1.6 Tobacco Use History
df <- adsub %>%
  select(SUBJID, PARAMCD, AVALC) %>%
  filter(PARAMCD %in% 'BTOBAC') %>%
  mutate(PT = as.numeric(sub(".*-","",SUBJID))) %>%
  arrange(PT) %>%
  dplyr::left_join(cTF, by = 'PT')
kruskal.test(log(df$cTF) ~ df$AVALC)
ggplot(df, aes(x = AVALC, y = log(cTF))) +
  geom_boxplot() +
  labs(x = "Tobacco Use History", y = "log(cTF)", 
       title = "log(cTF) vs Tobacco Use History", subtitle = "Kruskal-Wallis, p = 0.04541") +
  theme_bw() #600x500

# 1.7 CNS
df <- adsub %>%
  select(SUBJID, PARAMCD, AVALC) %>%
  filter(PARAMCD %in% 'CNSFL') %>%
  mutate(PT = as.numeric(sub(".*-","",SUBJID))) %>%
  arrange(PT) %>%
  dplyr::left_join(cTF, by = 'PT')
wilcox.test(log(df$cTF) ~ df$AVALC)
ggplot(df, aes(x = AVALC, y = log(cTF))) +
  geom_boxplot() +
  labs(x = "CNS Disease Flag", y = "log(cTF)", 
       title = "log(cTF) vs CNS Disease Flag", subtitle = "Wilcoxon, p = 0.1271") +
  theme_bw() #600x500

# 1.8 EML4
df <- adsub %>%
  select(SUBJID, PARAMCD, AVALC) %>%
  filter(PARAMCD %in% 'EML4FL') %>%
  mutate(PT = as.numeric(sub(".*-","",SUBJID))) %>%
  arrange(PT) %>%
  dplyr::left_join(cTF, by = 'PT')
wilcox.test(log(df$cTF) ~ df$AVALC)
ggplot(df, aes(x = AVALC, y = log(cTF))) +
  geom_boxplot() +
  labs(x = "EML4", y = "log(cTF)", 
       title = "log(cTF) vs EML4", subtitle = "Wilcoxon, p = 0.6194") +
  theme_bw() #600x500

# 1.9 Histology
df <- adsub %>%
  select(SUBJID, PARAMCD, AVALC) %>%
  filter(PARAMCD %in% 'HIST') %>%
  mutate(PT = as.numeric(sub(".*-","",SUBJID))) %>%
  arrange(PT) %>%
  dplyr::left_join(cTF, by = 'PT')
kruskal.test(log(df$cTF) ~ df$AVALC)
ggplot(df, aes(x = AVALC, y = log(cTF))) +
  geom_boxplot() +
  labs(x = "Histology", y = "log(cTF)", 
       title = "log(cTF) vs Histology", subtitle = "Kruskal-Wallis, p = 0.3932") +
  theme_bw() #600x500

# 1.10 Investigator baseline number of mets
tu <- rice_read("root/clinical_studies/RO5424802/CDT30119/BO29554/data_analysis/Cohort_A_BioPub2023/data/sdtmv/tu.sas7bdat")
mets <- tu %>%
  filter(VISIT=="SCREENING" & TUEVAL=="INVESTIGATOR") %>%
  group_by(USUBJID,TUEVAL) %>%
  summarize(n=n()) %>%
  mutate(PT = as.numeric(sub(".*-","",USUBJID))) %>%
  arrange(PT) %>%
  dplyr::left_join(cTF, by = 'PT') %>%
  mutate(logcTF = log(cTF))
ggscatter(mets, y = "logcTF", x = "n", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          ylab = "log(cTF)", xlab = "Investigator baseline number of mets")

# 1.11 IRF baseline number of mets
mets <- tu %>%
  filter(VISIT=="SCREENING" & TUEVAL=="INDEPENDENT ASSESSOR" & TUACPTFL=="Y") %>%
  group_by(USUBJID,TUEVAL) %>%
  summarize(n=n()) %>%
  mutate(PT = as.numeric(sub(".*-","",USUBJID))) %>%
  arrange(PT) %>%
  dplyr::left_join(cTF, by = 'PT') %>%
  mutate(logcTF = log(cTF))
ggscatter(mets, y = "logcTF", x = "n", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          ylab = "log(cTF)", xlab = "IRF baseline number of mets")

# 1.12 Sex
df <- adsub %>%
  select(SUBJID, SEX) %>%
  mutate(PT = as.numeric(sub(".*-","",SUBJID)),
         SEX = as.factor(SEX)) %>%
  arrange(PT) %>%
  dplyr::left_join(cTF, by = 'PT') %>% 
  distinct()
median_by_group <- df %>%
  group_by(SEX) %>%
  summarize(median_cTF = median(cTF))
kruskal.test(log(df$cTF) ~ df$SEX)
ggplot(df, aes(x = SEX, y = log(cTF))) +
  geom_boxplot() +
  labs(x = "Sex", y = "log(cTF)", title = "log(cTF) vs Sex", subtitle = "Kruskal-Wallis, p = 0.3366") +
  theme_bw() #600x500

# 2.0 PFS/OS km plot
adtte <- rice_read("root/clinical_studies/RO5424802/CDT30119/BO29554/data_analysis/Cohort_A_BioPub2023/prod/outdata_vad/adtte.sas7bdat")

# 2.1.1 PFS Baseline cTF <1% vs >1%
OSPL_cTF <- OSPL %>% 
  select(PT, VISIT, cTF) %>% 
  filter((VISIT != "TRTDISC") & (!is.na(cTF))) %>% 
  pivot_wider(names_from = VISIT, values_from = cTF) %>% 
  mutate(cTF_above_1perc = (SCRN > 0.01)) %>% 
  mutate(PT = as.character(PT))
adtte_PFS <- adtte %>% 
  filter(PARAMCD == "PFS") %>% 
  select(SUBJID, AVAL, CNSR) %>% 
  separate(col=SUBJID, into=c("Study", "Site", "PT"), sep="-")
OSPL_cTF_adtte_PFS <- left_join(x = OSPL_cTF, y = adtte_PFS, by = "PT") %>%
  mutate(EVENT = CNSR == 0) %>%
  mutate(cTF_above_1perc = factor(cTF_above_1perc, levels=c(TRUE, FALSE)))
OSPL_cTF_adtte_PFS$cTF_above_1perc = relevel(OSPL_cTF_adtte_PFS$cTF_above_1perc, ref = "TRUE")
fit_survdiff <- survdiff(Surv(AVAL, event = EVENT) ~ cTF_above_1perc, data = OSPL_cTF_adtte_PFS)
pval <- pchisq(fit_survdiff$chisq, length(fit_survdiff$n)-1, lower.tail = FALSE)
fit_coxmodel <- coxph(Surv(AVAL, event = EVENT) ~ cTF_above_1perc, data = OSPL_cTF_adtte_PFS, ties = "exact")
summary(fit_coxmodel)$conf.int
fit_survfit <- survfit(Surv(AVAL, event = EVENT) ~ cTF_above_1perc, data = OSPL_cTF_adtte_PFS)
fit_survfit
km_plot <- ggsurvplot(fit_survfit, 
                      legend.labs = c("cTF ≥ 1%", "cTF < 1%"),
                      legend.position = "top",
                      pval = FALSE,
                      palette = c("red", "blue"),
                      risk.table = TRUE, risk.table.y.text.col = TRUE, risk.table.col = "strata",
                      xlab = "Months", ylab = "Progression-free Survival",
                      break.time.by = 5,
                      font.main = 16,
                      risk.table.title = NULL,      # Remove the risk table title
                      legend.title = NULL)
# HR=0.61 [0.36-1.03] p = 0.0609
km_plot$plot <- km_plot$plot + 
  geom_segment(aes(x = 0, xend = 41.5, y = 0.5, yend = 0.5), 
               linetype = "dashed", color = "darkgray", linewidth = 0.7) +
  geom_segment(aes(x = 18.4, xend = 18.4, y = 0, yend = 0.5), 
               linetype = "dashed", color = "red", linewidth = 0.7) +
  geom_segment(aes(x = 41.5, xend = 41.5, y = 0, yend = 0.5), 
               linetype = "dashed", color = "blue", linewidth = 0.7) +
  theme(legend.position = c(0.85, 0.85),
        legend.title = element_blank()) 
km_plot #1100*700
                  
# 2.1.2 OS Baseline cTF <1% vs >1%
adtte_OS <- adtte %>% 
  filter(PARAMCD == "OS") %>% 
  select(SUBJID, AVAL, CNSR) %>% 
  separate(col=SUBJID, into=c("Study", "Site", "PT"), sep="-")
OSPL_cTF_adtte_OS <- left_join(x = OSPL_cTF, y = adtte_OS, by = "PT") %>%
  mutate(EVENT = CNSR == 0) %>%
  mutate(cTF_above_1perc = factor(cTF_above_1perc, levels=c(TRUE, FALSE)))
OSPL_cTF_adtte_OS$cTF_above_1perc = relevel(OSPL_cTF_adtte_OS$cTF_above_1perc, ref = "TRUE")
fit_survdiff <- survdiff(Surv(AVAL, event = EVENT) ~ cTF_above_1perc, data = OSPL_cTF_adtte_OS)
pval <- pchisq(fit_survdiff$chisq, length(fit_survdiff$n)-1, lower.tail = FALSE)
fit_coxmodel <- coxph(Surv(AVAL, event = EVENT) ~ cTF_above_1perc, data = OSPL_cTF_adtte_OS, ties = "exact")
summary(fit_coxmodel)$conf.int
fit_survfit <- survfit(Surv(AVAL, event = EVENT) ~ cTF_above_1perc, data = OSPL_cTF_adtte_OS)
fit_survfit
km_plot <- ggsurvplot(fit_survfit, 
                      legend.labs = c("cTF ≥ 1%", "cTF < 1%"),
                      legend.position = "top",
                      pval = FALSE,
                      palette = c("red", "blue"),
                      risk.table = TRUE, risk.table.y.text.col = TRUE, risk.table.col = "strata",
                      xlab = "Months", ylab = "Overall Survival",
                      break.time.by = 5,
                      font.main = 16,
                      risk.table.title = NULL,      # Remove the risk table title
                      legend.title = NULL)
# HR=0.44 [0.22-0.85] p = 0.0118
km_plot$plot <- km_plot$plot + 
  geom_segment(aes(x = 0, xend = 65, y = 0.5, yend = 0.5), 
               linetype = "dashed", color = "darkgray", linewidth = 0.7) +
  geom_segment(aes(x = 39.5, xend = 39.5, y = 0, yend = 0.5), 
               linetype = "dashed", color = "red", linewidth = 0.7) +
  theme(legend.position = c(0.85, 0.85),
        legend.title = element_blank()) 
km_plot #1100*700

# 2.2.1 PFS Baseline cTF <10% vs >10%
OSPL_cTF <- OSPL %>% 
  select(PT, VISIT, cTF) %>% 
  filter((VISIT != "TRTDISC") & (!is.na(cTF))) %>% 
  pivot_wider(names_from = VISIT, values_from = cTF) %>% 
  mutate(cTF_above_10perc = (SCRN > 0.1)) %>% 
  mutate(PT = as.character(PT))
adtte_PFS <- adtte %>% 
  filter(PARAMCD == "PFS") %>% 
  select(SUBJID, AVAL, CNSR) %>% 
  separate(col=SUBJID, into=c("Study", "Site", "PT"), sep="-")
OSPL_cTF_adtte_PFS <- left_join(x = OSPL_cTF, y = adtte_PFS, by = "PT") %>%
  mutate(EVENT = CNSR == 0) %>%
  mutate(cTF_above_10perc = factor(cTF_above_10perc, levels=c(TRUE, FALSE)))
OSPL_cTF_adtte_PFS$cTF_above_10perc = relevel(OSPL_cTF_adtte_PFS$cTF_above_10perc, ref = "TRUE")
fit_survdiff <- survdiff(Surv(AVAL, event = EVENT) ~ cTF_above_10perc, data = OSPL_cTF_adtte_PFS)
pval <- pchisq(fit_survdiff$chisq, length(fit_survdiff$n)-1, lower.tail = FALSE)
fit_coxmodel <- coxph(Surv(AVAL, event = EVENT) ~ cTF_above_10perc, data = OSPL_cTF_adtte_PFS, ties = "exact")
summary(fit_coxmodel)$conf.int
fit_survfit <- survfit(Surv(AVAL, event = EVENT) ~ cTF_above_10perc, data = OSPL_cTF_adtte_PFS)
fit_survfit
km_plot <- ggsurvplot(fit_survfit, 
                      legend.labs = c("cTF ≥ 10%", "cTF < 10%"),
                      legend.position = "top",
                      pval = FALSE,
                      palette = c("red", "blue"),
                      risk.table = TRUE, risk.table.y.text.col = TRUE, risk.table.col = "strata",
                      xlab = "Months", ylab = "Progression-free Survival",
                      break.time.by = 5,
                      font.main = 16,
                      risk.table.title = NULL,      # Remove the risk table title
                      legend.title = NULL)
#surv_at_median_1 <- summary(fit_survfit, times = 12.7)$surv
#surv_at_median_2 <- summary(fit_survfit, times = 36.9)$surv
# HR=0.39 [0.20-0.76] p = 0.0041
km_plot$plot <- km_plot$plot + 
  geom_segment(aes(x = 0, xend = 36.9, y = 0.5, yend = 0.5), 
               linetype = "dashed", color = "darkgray", linewidth = 0.7) +
  geom_segment(aes(x = 12.7, xend = 12.7, y = 0, yend = 0.5), 
               linetype = "dashed", color = "red", linewidth = 0.7) +
  geom_segment(aes(x = 36.9, xend = 36.9, y = 0, yend = 0.5), 
               linetype = "dashed", color = "blue", linewidth = 0.7) +
  theme(legend.position = c(0.85, 0.85),
        legend.title = element_blank()) 
km_plot #1100*700

# 2.2.2 OS Baseline cTF <10% vs >10%
adtte_OS <- adtte %>% 
  filter(PARAMCD == "OS") %>% 
  select(SUBJID, AVAL, CNSR) %>% 
  separate(col=SUBJID, into=c("Study", "Site", "PT"), sep="-")
OSPL_cTF_adtte_OS <- left_join(x = OSPL_cTF, y = adtte_OS, by = "PT") %>%
  mutate(EVENT = CNSR == 0) %>%
  mutate(cTF_above_10perc = factor(cTF_above_10perc, levels=c(TRUE, FALSE)))
OSPL_cTF_adtte_OS$cTF_above_10perc = relevel(OSPL_cTF_adtte_OS$cTF_above_10perc, ref = "TRUE")
fit_survdiff <- survdiff(Surv(AVAL, event = EVENT) ~ cTF_above_10perc, data = OSPL_cTF_adtte_OS)
pval <- pchisq(fit_survdiff$chisq, length(fit_survdiff$n)-1, lower.tail = FALSE)
fit_coxmodel <- coxph(Surv(AVAL, event = EVENT) ~ cTF_above_10perc, data = OSPL_cTF_adtte_OS, ties = "exact")
summary(fit_coxmodel)$conf.int
fit_survfit <- survfit(Surv(AVAL, event = EVENT) ~ cTF_above_10perc, data = OSPL_cTF_adtte_OS)
fit_survfit
km_plot <- ggsurvplot(fit_survfit, 
                      legend.labs = c("cTF ≥ 10%", "cTF < 10%"),
                      legend.position = "top",
                      font.main = 16, 
                      pval = FALSE,
                      palette = c("red", "blue"),
                      risk.table = TRUE, risk.table.y.text.col = TRUE, risk.table.col = "strata",
                      xlab = "Months", ylab = "Overall Survival",
                      break.time.by = 5)
#HR=0.44 [0.22-0.85] p < 0.0001
km_plot$plot <- km_plot$plot + 
  geom_segment(aes(x = 0, xend = 62, y = 0.5, yend = 0.5), 
               linetype = "dashed", color = "darkgray", linewidth = 0.7) +
  geom_segment(aes(x = 18.2, xend = 18.2, y = 0, yend = 0.5), 
               linetype = "dashed", color = "red", linewidth = 0.7) +
  theme(legend.position = c(0.85, 0.85),
        legend.title = element_blank()) 
km_plot #1100*700


#########################################################################
# 2.3.1 PFS KM plot for EML4 fusion partner vs non-EML4 fusion partners #
#########################################################################
EMLVariant <- read.csv("~/ctdna/data_cleaning_20221219/cohortA/bo29554_ngs_dna_targdna_foundationact_17-567_n186_ALK_20200109.csv")
EML4_all <- EMLVariant %>% 
  filter(ptid != "#N/A") %>% 
  mutate(PT = as.character(ptid), EML4Variant = mutation_label) %>% 
  select(c(EML4Variant, PT)) %>%
  mutate(Type = "EML4")
adtte_PFS <- adtte %>% 
  filter(PARAMCD == "PFS") %>% 
  select(SUBJID, AVAL, CNSR) %>% 
  separate(col=SUBJID, into=c("Study", "Site", "PT"), sep="-")
EML4_adtte_km_pfs <- left_join(x = adtte_PFS, y = EML4_all, by = "PT") %>%
  mutate(EVENT = CNSR == 0) %>%
  replace(is.na(.), "others") %>%
  mutate(Type = factor(Type, levels=c("others", "EML4")))
fit_survdiff <- survdiff(Surv(AVAL, event = EVENT) ~ Type, data = EML4_adtte_km_pfs)
pval <- pchisq(fit_survdiff$chisq, length(fit_survdiff$n)-1, lower.tail = FALSE)
fit_coxmodel <- coxph(Surv(AVAL, event = EVENT) ~ Type, data = EML4_adtte_km_pfs, ties = "exact")
summary(fit_coxmodel)$conf.int
fit_survfit <- survfit(Surv(AVAL, event = EVENT) ~ Type, data = EML4_adtte_km_pfs)
fit_survfit
km_plot <- ggsurvplot(fit_survfit, 
                      legend.labs = c("EML4 fusion", "non-EML4 fusion"),
                      legend.position = "top",
                      pval = FALSE,
                      palette = c("red", "blue"),
                      risk.table = TRUE, risk.table.y.text.col = TRUE, risk.table.col = "strata",
                      xlab = "Months", ylab = "Progression-free Survival",
                      break.time.by = 5,
                      font.main = 16,
                      risk.table.title = NULL,      # Remove the risk table title
                      legend.title = NULL)
# HR=0.99 [0.49-2.02] p = 0.9790
km_plot$plot <- km_plot$plot + 
  geom_segment(aes(x = 0, xend = 33.3, y = 0.5, yend = 0.5), 
               linetype = "dashed", color = "darkgray", linewidth = 0.7) +
  geom_segment(aes(x = 33.3, xend = 33.3, y = 0, yend = 0.5), 
               linetype = "dashed", color = "red", linewidth = 0.7) +
  geom_segment(aes(x = 33.0, xend = 33.0, y = 0, yend = 0.5), 
               linetype = "dashed", color = "blue", linewidth = 0.7) +
  theme(legend.position = c(0.85, 0.85),
        legend.title = element_blank()) 
km_plot #1100*700

########################################################################
# 2.3.2 OS KM plot for EML4 fusion partner vs non-EML4 fusion partners #
########################################################################
adtte_OS <- adtte %>% 
  filter(PARAMCD == "OS") %>% 
  select(SUBJID, AVAL, CNSR) %>% 
  separate(col=SUBJID, into=c("Study", "Site", "PT"), sep="-")
EML4_adtte_km_os <- left_join(x = adtte_OS, y = EML4_all, by = "PT") %>%
  mutate(EVENT = CNSR == 0) %>%
  replace(is.na(.), "others") %>%
  mutate(Type = factor(Type, levels=c("others", "EML4")))
fit_survdiff <- survdiff(Surv(AVAL, event = EVENT) ~ Type, data = EML4_adtte_km_os)
pval <- pchisq(fit_survdiff$chisq, length(fit_survdiff$n)-1, lower.tail = FALSE)
fit_coxmodel <- coxph(Surv(AVAL, event = EVENT) ~ Type, data = EML4_adtte_km_os, ties = "exact")
summary(fit_coxmodel)$conf.int
fit_survfit <- survfit(Surv(AVAL, event = EVENT) ~ Type, data = EML4_adtte_km_os)
fit_survfit
km_plot <- ggsurvplot(fit_survfit, 
                      legend.labs = c("EML4 fusion", "non-EML4 fusion"),
                      legend.position = "top",
                      font.main = 16, 
                      pval = FALSE,
                      palette = c("red", "blue"),
                      risk.table = TRUE, risk.table.y.text.col = TRUE, risk.table.col = "strata",
                      xlab = "Months", ylab = "Overall Survival",
                      break.time.by = 5)
# HR=0.80 [0.35-1.82] p = 0.5964
km_plot$plot <- km_plot$plot + 
  geom_segment(aes(x = 0, xend = 65, y = 0.5, yend = 0.5), 
               linetype = "dashed", color = "darkgray", linewidth = 0.7) +
  geom_segment(aes(x = 56.8, xend = 56.8, y = 0, yend = 0.5), 
               linetype = "dashed", color = "red", linewidth = 0.7) +
  theme(legend.position = c(0.85, 0.85),
        legend.title = element_blank()) 
km_plot #1100*700

##########################################################################
# 2.4.1 PFS KM plot for EML4 fusion partner v1 vs EML4 fusion partner v3 #
##########################################################################
EML4 <- EMLVariant %>% 
  filter(ptid != "#N/A" & mutation_label %in% c("ALK_V1", "ALK_V3")) %>% 
  mutate(PT = as.character(ptid), EML4Variant = mutation_label) %>% 
  select(c(EML4Variant, PT)) 
EML4_adtte_km_pfs <- left_join(x = EML4, y = adtte_PFS, by = "PT") %>%
  mutate(EVENT = CNSR == 0) %>%
  mutate(EML4Variant = factor(EML4Variant, levels=c("ALK_V3", "ALK_V1")))
fit_survdiff <- survdiff(Surv(AVAL, event = EVENT) ~ EML4Variant, data = EML4_adtte_km_pfs)
pval <- pchisq(fit_survdiff$chisq, length(fit_survdiff$n)-1, lower.tail = FALSE)
fit_coxmodel <- coxph(Surv(AVAL, event = EVENT) ~ EML4Variant, data = EML4_adtte_km_pfs, ties = "exact")
summary(fit_coxmodel)$conf.int
fit_survfit <- survfit(Surv(AVAL, event = EVENT) ~ EML4Variant, data = EML4_adtte_km_pfs)
fit_survfit
km_plot <- ggsurvplot(fit_survfit, 
                      legend.labs = c("EML4 fusion v1", "EML4 fusion v3"),
                      legend.position = "top",
                      pval = FALSE,
                      palette = c("red", "blue"),
                      risk.table = TRUE, risk.table.y.text.col = TRUE, risk.table.col = "strata",
                      xlab = "Months", ylab = "Progression-free Survival",
                      break.time.by = 5,
                      font.main = 16,
                      risk.table.title = NULL,      # Remove the risk table title
                      legend.title = NULL)
# HR=0.79 [0.40-1.55] p = 0.4852
km_plot$plot <- km_plot$plot + 
  geom_segment(aes(x = 0, xend = 33.3, y = 0.5, yend = 0.5), 
               linetype = "dashed", color = "darkgray", linewidth = 0.7) +
  geom_segment(aes(x = 29.4, xend = 29.4, y = 0, yend = 0.5), 
               linetype = "dashed", color = "red", linewidth = 0.7) +
  geom_segment(aes(x = 33.3, xend = 33.3, y = 0, yend = 0.5), 
               linetype = "dashed", color = "blue", linewidth = 0.7) +
  theme(legend.position = c(0.85, 0.85),
        legend.title = element_blank()) 
km_plot #1100*700

#########################################################################
# 2.4.2 OS KM plot for EML4 fusion partner v1 vs EML4 fusion partner v3 #
#########################################################################
EML4_adtte_km_os <- left_join(x = EML4, y = adtte_OS, by = "PT") %>%
  mutate(EVENT = CNSR == 0) %>%
  mutate(EML4Variant = factor(EML4Variant, levels=c("ALK_V3", "ALK_V1")))
fit_survdiff <- survdiff(Surv(AVAL, event = EVENT) ~ EML4Variant, data = EML4_adtte_km_os)
pval <- pchisq(fit_survdiff$chisq, length(fit_survdiff$n)-1, lower.tail = FALSE)
fit_coxmodel <- coxph(Surv(AVAL, event = EVENT) ~ EML4Variant, data = EML4_adtte_km_os, ties = "exact")
summary(fit_coxmodel)$conf.int
fit_survfit <- survfit(Surv(AVAL, event = EVENT) ~ EML4Variant, data = EML4_adtte_km_os)
fit_survfit
km_plot <- ggsurvplot(fit_survfit, 
                      legend.labs = c("EML4 fusion v1", "EML4 fusion v3"),
                      legend.position = "top",
                      font.main = 16, 
                      pval = FALSE,
                      palette = c("red", "blue"),
                      risk.table = TRUE, risk.table.y.text.col = TRUE, risk.table.col = "strata",
                      xlab = "Months", ylab = "Overall Survival",
                      break.time.by = 5)
# HR=0.73 [0.31-1.72] p = 0.4699
km_plot$plot <- km_plot$plot + 
  geom_segment(aes(x = 0, xend = 65, y = 0.5, yend = 0.5), 
               linetype = "dashed", color = "darkgray", linewidth = 0.7) +
  theme(legend.position = c(0.85, 0.85),
        legend.title = element_blank()) 
km_plot #1100*700

###########################################
# 2.5.1 PFS KM plot for wtTP53 vs mutTP53 #
###########################################
OAPL_TP53 <- OAPL %>%
  filter(gene == "TP53" & somatic_status_functional_impact %in% c("known", "likely") & VISIT == "SCRN") %>%
  select(PT) %>%
  distinct() %>%
  mutate(PT = as.character(PT), Type = "TP53 commutation")
EML4_adtte_km_pfs <- left_join(x = adtte_PFS, y = OAPL_TP53, by = "PT") %>%
  mutate(EVENT = CNSR == 0) %>%
  replace(is.na(.), "TP53 wild type") %>%
  mutate(Type = factor(Type, levels=c("TP53 commutation", "TP53 wild type")))
fit_survdiff <- survdiff(Surv(AVAL, event = EVENT) ~ Type, data = EML4_adtte_km_pfs)
pval <- pchisq(fit_survdiff$chisq, length(fit_survdiff$n)-1, lower.tail = FALSE)
fit_coxmodel <- coxph(Surv(AVAL, event = EVENT) ~ Type, data = EML4_adtte_km_pfs, ties = "exact")
summary(fit_coxmodel)$conf.int
fit_survfit <- survfit(Surv(AVAL, event = EVENT) ~ Type, data = EML4_adtte_km_pfs)
fit_survfit
km_plot <- ggsurvplot(fit_survfit, 
                      legend.labs = c("Wild-type TP53", "Mutant TP53"),
                      legend.position = "top",
                      pval = FALSE,
                      palette = c("red", "blue"),
                      risk.table = TRUE, risk.table.y.text.col = TRUE, risk.table.col = "strata",
                      xlab = "Months", ylab = "Progression-free Survival",
                      break.time.by = 5,
                      font.main = 16,
                      risk.table.title = NULL,      # Remove the risk table title
                      legend.title = NULL)
# HR=0.59 [0.35-1.00] p = 0.0460
km_plot$plot <- km_plot$plot + 
  geom_segment(aes(x = 0, xend = 39.5, y = 0.5, yend = 0.5), 
               linetype = "dashed", color = "darkgray", linewidth = 0.7) +
  geom_segment(aes(x = 18.7, xend = 18.7, y = 0, yend = 0.5), 
               linetype = "dashed", color = "red", linewidth = 0.7) +
  geom_segment(aes(x = 39.5, xend = 39.5, y = 0, yend = 0.5), 
               linetype = "dashed", color = "blue", linewidth = 0.7) +
  theme(legend.position = c(0.85, 0.85),
        legend.title = element_blank()) 
km_plot #1100*700

##########################################
# 2.5.2 OS KM plot for wtTP53 vs mutTP53 #
##########################################
EML4_adtte_km_os <- left_join(x = adtte_OS, y = OAPL_TP53, by = "PT") %>%
  mutate(EVENT = CNSR == 0) %>%
  replace(is.na(.), "TP53 wild type") %>%
  mutate(Type = factor(Type, levels=c("TP53 commutation", "TP53 wild type")))
fit_survdiff <- survdiff(Surv(AVAL, event = EVENT) ~ Type, data = EML4_adtte_km_os)
pval <- pchisq(fit_survdiff$chisq, length(fit_survdiff$n)-1, lower.tail = FALSE)
fit_coxmodel <- coxph(Surv(AVAL, event = EVENT) ~ Type, data = EML4_adtte_km_os, ties = "exact")
summary(fit_coxmodel)$conf.int
fit_survfit <- survfit(Surv(AVAL, event = EVENT) ~ Type, data = EML4_adtte_km_os)
fit_survfit
km_plot <- ggsurvplot(fit_survfit, 
                      legend.labs = c("Wild-type TP53", "Mutant TP53"),
                      legend.position = "top",
                      font.main = 16, 
                      pval = FALSE,
                      palette = c("red", "blue"),
                      risk.table = TRUE, risk.table.y.text.col = TRUE, risk.table.col = "strata",
                      xlab = "Months", ylab = "Overall Survival",
                      break.time.by = 5)
# HR=0.41 [0.22-0.78] p = 0.0045
km_plot$plot <- km_plot$plot + 
  geom_segment(aes(x = 0, xend = 65, y = 0.5, yend = 0.5), 
               linetype = "dashed", color = "darkgray", linewidth = 0.7) +
  geom_segment(aes(x = 34, xend = 34, y = 0, yend = 0.5), 
               linetype = "dashed", color = "red", linewidth = 0.7) +
  theme(legend.position = c(0.85, 0.85),
        legend.title = element_blank()) 
km_plot #1100*700

#################################################################
# 2.6.1 PFS KM plot for ctDNA clearance at C3D1 vs no clearance #
#################################################################
OSPL_C3D1 <- OSPL %>%
  filter(VISIT=='C3D1' & assay == "F1LCdX") %>%
  select(PT, fmi_sample_id) # N=71 patients with Cycle 3 Day 1 (C3D1) data
OAPL_RA_ALK <- OAPL %>%
  filter(variant_type == "rearrangement" & gene=="ALK") %>%
  #filter(assay == "F1LCdX") %>%
  filter((PT %in% OSPL_C3D1$PT & VISIT == "SCRN") | fmi_sample_id %in% OSPL_C3D1$fmi_sample_id) %>%
  mutate(PT = as.character(PT))
BEP <- unique(OAPL_RA_ALK$PT)
no_clearance <- (OAPL_RA_ALK %>% filter(VISIT %in% 'C3D1'))$PT
OAPL_RA_ALK_adtte_km_pfs <- adtte_PFS %>%
  filter(PT %in% BEP) %>%
  mutate(clearance_C3D1 = !(PT %in% no_clearance)) %>%
  mutate(EVENT = CNSR == 0) %>%
  mutate(clearance_C3D1 = factor(clearance_C3D1, levels=c(FALSE, TRUE)))
fit_survdiff <- survdiff(Surv(AVAL, event = EVENT) ~ clearance_C3D1, data = OAPL_RA_ALK_adtte_km_pfs)
pval <- pchisq(fit_survdiff$chisq, length(fit_survdiff$n)-1, lower.tail = FALSE)
fit_coxmodel <- coxph(Surv(AVAL, event = EVENT) ~ clearance_C3D1, data = OAPL_RA_ALK_adtte_km_pfs, ties = "exact")
summary(fit_coxmodel)$conf.int
fit_survfit <- survfit(Surv(AVAL, event = EVENT) ~ clearance_C3D1, data = OAPL_RA_ALK_adtte_km_pfs)
fit_survfit
km_plot <- ggsurvplot(fit_survfit, 
                      legend.labs = c("ALK Detected at C3D1", "ALK Cleared by C3D1"),
                      legend.position = "top",
                      pval = FALSE,
                      palette = c("red", "blue"),
                      risk.table = TRUE, risk.table.y.text.col = TRUE, risk.table.col = "strata",
                      xlab = "Months", ylab = "Progression-free Survival",
                      break.time.by = 5,
                      font.main = 16,
                      risk.table.title = NULL,      # Remove the risk table title
                      legend.title = NULL)
# HR = 0.167 [0.062-0.451] p < 0.0001
km_plot$plot <- km_plot$plot + 
  geom_segment(aes(x = 0, xend = 36.8, y = 0.5, yend = 0.5), 
               linetype = "dashed", color = "darkgray", linewidth = 0.7) +
  geom_segment(aes(x = 7.1, xend = 7.1, y = 0, yend = 0.5), 
               linetype = "dashed", color = "red", linewidth = 0.7) +
  geom_segment(aes(x = 36.8, xend = 36.8, y = 0, yend = 0.5), 
               linetype = "dashed", color = "blue", linewidth = 0.7) +
  theme(legend.position = c(0.85, 0.85),
        legend.title = element_blank()) 
km_plot #1100*700


################################################################
# 2.6.2 OS KM plot for ctDNA clearance at C3D1 vs no clearance #
################################################################
OAPL_RA_ALK_adtte_km_os <- adtte_OS %>%
  filter(PT %in% BEP) %>%
  mutate(clearance_C3D1 = !(PT %in% no_clearance)) %>%
  mutate(EVENT = CNSR == 0) %>%
  mutate(clearance_C3D1 = factor(clearance_C3D1, levels=c(FALSE, TRUE)))
fit_survdiff <- survdiff(Surv(AVAL, event = EVENT) ~ clearance_C3D1, data = OAPL_RA_ALK_adtte_km_os)
pval <- pchisq(fit_survdiff$chisq, length(fit_survdiff$n)-1, lower.tail = FALSE)
fit_coxmodel <- coxph(Surv(AVAL, event = EVENT) ~ clearance_C3D1, data = OAPL_RA_ALK_adtte_km_os, ties = "exact")
summary(fit_coxmodel)$conf.int
fit_survfit <- survfit(Surv(AVAL, event = EVENT) ~ clearance_C3D1, data = OAPL_RA_ALK_adtte_km_os)
fit_survfit
km_plot <- ggsurvplot(fit_survfit, 
                      legend.labs = c("ALK Detected at C3D1", "ALK Cleared by C3D1"),
                      legend.position = "top",
                      pval = FALSE,
                      palette = c("red", "blue"),
                      risk.table = TRUE, risk.table.y.text.col = TRUE, risk.table.col = "strata",
                      xlab = "Months", ylab = "Overall Survival",
                      break.time.by = 5,
                      font.main = 16,
                      risk.table.title = NULL,      # Remove the risk table title
                      legend.title = NULL)
# HR = 0.155 [0.058-0.416] p < 0.0001
km_plot$plot <- km_plot$plot + 
  geom_segment(aes(x = 0, xend = 65, y = 0.5, yend = 0.5), 
               linetype = "dashed", color = "darkgray", linewidth = 0.7) +
  geom_segment(aes(x = 8.34, xend = 8.34, y = 0, yend = 0.5), 
               linetype = "dashed", color = "red", linewidth = 0.7) +
  theme(legend.position = c(0.85, 0.85),
        legend.title = element_blank()) 
km_plot #1100*700


#################################################################
# 2.7.1 PFS KM plot for Baseline ALK <median VAF vs >median VAF #
#################################################################
BFAST_cohortA_samples_SCRN_plasma_Jul_vaf <- read_excel("ctdna/DBI_20230508/BFAST_cohortA_samples_SCRN_plasma_Jul_vaf.xlsx")
alk <- BFAST_cohortA_samples_SCRN_plasma_Jul_vaf$percentReads_RE
median(alk) #0.6

wilcox.test(alk,
            alternative="two.sided",
            correct=TRUE,
            conf.int=TRUE,
            conf.level=0.95) # https://rcompanion.org/handbook/E_04.html
library(boot)
set.seed(1234)
Mboot = boot(alk,
             function(x,i) median(x[i]),
             R=100000)
boot.ci(Mboot,
        conf = 0.95,
        type = c("norm", "basic" ,"perc", "bca"))

SCRN_plasma <- BFAST_cohortA_samples_SCRN_plasma_Jul_vaf %>%
  rename(PT = `Subject ID`) %>%
  mutate(PT = as.character(PT),
         ALK_below_median = ifelse(percentReads_RE<0.6, TRUE, FALSE))
OAPL_RA_ALK_adtte_km_pfs <- left_join(x = SCRN_plasma, y = adtte_PFS, by = "PT") %>%
  mutate(EVENT = CNSR == 0) %>%
  mutate(ALK_below_median = factor(ALK_below_median, levels=c(FALSE, TRUE)))
fit_survdiff <- survdiff(Surv(AVAL, event = EVENT) ~ ALK_below_median, data = OAPL_RA_ALK_adtte_km_pfs)
pval <- pchisq(fit_survdiff$chisq, length(fit_survdiff$n)-1, lower.tail = FALSE)
fit_coxmodel <- coxph(Surv(AVAL, event = EVENT) ~ ALK_below_median, data = OAPL_RA_ALK_adtte_km_pfs, ties = "exact")
summary(fit_coxmodel)$conf.int
fit_survfit <- survfit(Surv(AVAL, event = EVENT) ~ ALK_below_median, data = OAPL_RA_ALK_adtte_km_pfs)
fit_survfit
km_plot <- ggsurvplot(fit_survfit, 
                      legend.labs = c("ALK VAF ≥ Median", "ALK VAF < Median"),
                      legend = "top",
                      main = "BEP n=77", font.main = 16, 
                      pval = FALSE,
                      palette = c("red", "blue"),
                      risk.table = TRUE, risk.table.y.text.col = TRUE, risk.table.col = "strata",
                      xlab = "Months", ylab = "Progression-free Survival",
                      break.time.by = 5)
# HR=0.67 [0.38-1.17] p=0.1578
km_plot$plot <- km_plot$plot + 
  geom_segment(aes(x = 0, xend = 36.2, y = 0.5, yend = 0.5), 
               linetype = "dashed", color = "darkgray", linewidth = 0.7) +
  geom_segment(aes(x = 22.2, xend = 22.2, y = 0, yend = 0.5), 
               linetype = "dashed", color = "red", linewidth = 0.7) +
  geom_segment(aes(x = 36.2, xend = 36.2, y = 0, yend = 0.5), 
               linetype = "dashed", color = "blue", linewidth = 0.7) +
  theme(legend.position = c(0.85, 0.85),
        legend.title = element_blank()) 
km_plot #1100*700


################################################################
# 2.7.2 OS KM plot for Baseline ALK <median VAF vs >median VAF #
################################################################
OAPL_RA_ALK_adtte_km_os <- left_join(x = SCRN_plasma, y = adtte_OS, by = "PT") %>%
  mutate(EVENT = CNSR == 0) %>%
  mutate(ALK_below_median = factor(ALK_below_median, levels=c(FALSE, TRUE)))
fit_survdiff <- survdiff(Surv(AVAL, event = EVENT) ~ ALK_below_median, data = OAPL_RA_ALK_adtte_km_os)
pval <- pchisq(fit_survdiff$chisq, length(fit_survdiff$n)-1, lower.tail = FALSE)
fit_coxmodel <- coxph(Surv(AVAL, event = EVENT) ~ ALK_below_median, data = OAPL_RA_ALK_adtte_km_os, ties = "exact")
summary(fit_coxmodel)$conf.int
fit_survfit <- survfit(Surv(AVAL, event = EVENT) ~ ALK_below_median, data = OAPL_RA_ALK_adtte_km_os)
fit_survfit
km_plot <- ggsurvplot(fit_survfit, 
                      legend.labs = c("ALK VAF ≥ Median", "ALK VAF < Median"),
                      legend = "top",
                      main = "BEP n=77", font.main = 16, 
                      pval = FALSE,
                      palette = c("red", "blue"),
                      risk.table = TRUE, risk.table.y.text.col = TRUE, risk.table.col = "strata",
                      xlab = "Months", ylab = "Overall Survival",
                      break.time.by = 5)
# HR=0.70 [0.36-1.38] p=0.3009
km_plot$plot <- km_plot$plot + 
  geom_segment(aes(x = 0, xend = 65, y = 0.5, yend = 0.5), 
               linetype = "dashed", color = "darkgray", linewidth = 0.7) +
  geom_segment(aes(x = 50, xend = 50, y = 0, yend = 0.5), 
               linetype = "dashed", color = "red", linewidth = 0.7) +
  theme(legend.position = c(0.85, 0.85),
        legend.title = element_blank()) 
km_plot #1100*700

####################################################################
# 3.0 The ALK fusion was not detected at treatment discontinuation #
####################################################################
TRTDISC <- OSPL %>%
  filter(VISIT %in% 'TRTDISC') %>%
  arrange(PT)
trtdisc <- unique(TRTDISC$PT)
no_ALK_TRTDISC <- OAPL %>%
  filter(VISIT %in% 'TRTDISC', variant_type %in% 'rearrangement', rearr_gene1 %in% 'ALK' | rearr_gene2 %in% 'ALK') %>%
  arrange(PT)
ALK_trtdisc <- unique(no_ALK_TRTDISC$PT)
34-21
trtdisc[!(trtdisc %in% ALK_trtdisc)]

####################
# 4.0 Longitudinal #
####################
pts <- c()
for (pt in unique(OSPL$PT)) {
  tmp <- OSPL %>% filter(PT %in% pt)
  if (all(c('SCRN', 'C3D1', 'TRTDISC') %in% tmp$VISIT)) {
    if (nrow(tmp) >= 5) {
      pts <- c(pts, pt)
    }
  }
}

adtr <- rice_read("root/clinical_studies/RO5424802/CDT30119/BO29554/data_analysis/Cohort_A_BioPub2023/prod/outdata_vad/adtr.sas7bdat"
                  ,prolong = TRUE)
adrs <- rice_read("root/clinical_studies/RO5424802/CDT30119/BO29554/data_analysis/Cohort_A_BioPub2023/prod/outdata_vad/adrs.sas7bdat")
rice_session_close()

# 4.1 SLD
SLD <- adtr %>%
  filter(PARAMCD %in% c('SD', 'SUMLDIAM')) %>%
  mutate(PT = as.numeric(sub(".*-","",SUBJID))) %>%
  filter(PT %in% 100087, AVALC != "NE") %>%
  arrange(ADY) %>%
  mutate(AVALC = as.numeric(AVALC),
         AVISIT = factor(AVISIT, levels = unique(AVISIT))) %>%
  select(PT, PARAM, AVALC, AVALU, AVISIT, ADT, ADY)
ggplot(SLD, aes(x = ADY, y = AVALC, group = PARAM, color = PARAM)) +
  geom_line() +
  geom_point(aes(shape = ifelse(AVALC == 0, "circle", "other")), show.legend = FALSE) +
  scale_shape_manual(values = c("circle" = 21, "other" = 16)) +
  labs(
    title = "100087",
    x = "Days from randomization",
    y = "SLD (mm)"
  ) +
  geom_hline(yintercept = 0, color = "black", linetype = "solid") +  # Highlight y=0
  geom_vline(xintercept = 0, color = "black", linetype = "solid") +  # Highlight x=0
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        legend.position = "top")

# 4.2 cTF
samples_date <- read.csv("~/bin/BFAST_cohortA_publication/samples_date.csv") %>%
  select(ACCSNM, BMSD) %>%
  distinct()
# Find the row with the minimum ADT in SLD
min_row <- SLD %>% filter(ADT == min(ADT))
# Extract the reference date and ADY
reference_ADT <- min_row$ADT[1]
reference_ADY <- min_row$ADY[1]
cTF <- OSPL %>%
  dplyr::select(PT, VISIT, cTF, assay, fmi_sample_id) %>%
  filter(PT %in% 100087) %>%
  left_join(samples_date, by = c("fmi_sample_id" = "ACCSNM")) %>%
  mutate(ADT = as.Date(as.character(BMSD), format = "%Y%m%d"),
         ADY = as.numeric(ADT - reference_ADT) + reference_ADY) %>%
  arrange(ADY)
ggplot(cTF, aes(x = ADY, y = cTF, group = PT)) +
  geom_line(color = "blue") +
  geom_point(aes(shape = factor(ifelse(cTF == 0, 1, 16))), color = "blue") +
  scale_shape_manual(
    values = c(1, 16),
    labels = c("Clearance", "Not Clearance")) +  # Define shapes for the factor levels
  labs(
    title = "100087",
    x = "Days from randomization",
    y = "cTF",
    shape = "cTF"
  ) +
  geom_hline(yintercept = 0, color = "black", linetype = "solid") +  # Highlight y=0
  geom_vline(xintercept = 0, color = "black", linetype = "solid") +  # Highlight x=0
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        legend.position = "top")

# 4.3 RE
RE <- OAPL %>%
  dplyr::select(PT, VISIT, mutID1, rearr_number_of_reads, somatic_status_functional_impact, fmi_sample_id) %>%
  filter(PT %in% 100087,
         mutID1 %in% "ALK_RE") %>%
  mutate(rearr_number_of_reads = as.numeric(as.character(rearr_number_of_reads)))
RE <- cTF %>%
  left_join(RE, by = c("fmi_sample_id")) %>%
  mutate(rearr_number_of_reads = ifelse(is.na(rearr_number_of_reads), 0, rearr_number_of_reads)) 
ggplot(RE, aes(x = ADY, y = rearr_number_of_reads)) +
  geom_line(color = "blue") +
  geom_point(aes(shape = factor(ifelse(rearr_number_of_reads == 0, 1, 16))), color = "blue") +
  scale_shape_manual(
    values = c(1, 16),
    labels = c("Clearance", "Not Clearance")) +  # Define shapes for the factor levels
  labs(
    title = "100087",
    x = "Days from randomization",
    y = "ALK Total Number of Rearrangement Reads",
    shape = "ALK Rearrangement"
  ) +
  geom_hline(yintercept = 0, color = "black", linetype = "solid") +  # Highlight y=0
  geom_vline(xintercept = 0, color = "black", linetype = "solid") +  # Highlight x=0
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        legend.position = "top")

# 4.4 ORR
df <- adrs %>%
  mutate(PT = as.numeric(sub(".*-","",SUBJID))) %>%
  filter(PT %in% 100087, PARAMCD %in% 'OVR') %>%
  arrange(ADY) %>%
  select(SUBJID, PARAM, AVALC, AVISIT)


df <- OSPL %>%
  filter(PT %in% 100076)
write.csv(df, file = "/gstore/home/jiangt17/bin/OSPL_100076.csv")
df <- OAPL %>%
  filter(PT %in% 100076)
write.csv(df, file = "/gstore/home/jiangt17/bin/OAPL_100076.csv")

###########
# 5.0 SLD #
###########
# 5.1 Percentage of patients larger than or equal to 100mm SLD vs less than 100mm SLD
adtr <- rice_read("root/clinical_studies/RO5424802/CDT30119/BO29554/data_analysis/Cohort_A_BioPub2023/prod/outdata_vad/adtr.sas7bdat")
SLD <- adtr %>%
  filter(ABLFL %in% 'Y') %>%
  select(SUBJID, PARAMCD, AVALC) %>%
  filter(PARAMCD %in% 'SD') %>%
  mutate(PT = as.numeric(sub(".*-","",SUBJID)),
         AVALC = as.numeric(AVALC)) %>%
  arrange(PT) %>%
  dplyr::left_join(cTF, by = 'PT') %>%
  mutate(logcTF = log(cTF))
SLD %>%
  summarise(
    greater_equal_100mm = sum(AVALC >= 100) / n() * 100,
    less_than_100mm = sum(AVALC < 100) / n() * 100
  )

# 5.2 PFS KM plot of patients ≥100mm SLD vs less than 100mm SLD 
adtte <- rice_read("root/clinical_studies/RO5424802/CDT30119/BO29554/data_analysis/Cohort_A_BioPub2023/prod/outdata_vad/adtte.sas7bdat")
SLD_100 <- SLD %>% 
  mutate(sld_above_100 = (AVALC >= 100)) %>% 
  mutate(PT = as.character(PT))
adtte_PFS <- adtte %>% 
  filter(PARAMCD == "PFS") %>% 
  select(SUBJID, AVAL, CNSR) %>% 
  separate(col=SUBJID, into=c("Study", "Site", "PT"), sep="-")
SLD_100_adtte_PFS <- left_join(x = SLD_100, y = adtte_PFS, by = "PT") %>%
  mutate(EVENT = CNSR == 0) %>%
  mutate(sld_above_100 = factor(sld_above_100, levels=c(TRUE, FALSE)))
SLD_100_adtte_PFS$sld_above_100 = relevel(SLD_100_adtte_PFS$sld_above_100, ref = "TRUE")
fit_survdiff <- survdiff(Surv(AVAL, event = EVENT) ~ sld_above_100, data = SLD_100_adtte_PFS)
pval <- pchisq(fit_survdiff$chisq, length(fit_survdiff$n)-1, lower.tail = FALSE)
fit_coxmodel <- coxph(Surv(AVAL, event = EVENT) ~ sld_above_100, data = SLD_100_adtte_PFS, ties = "exact")
summary(fit_coxmodel)$conf.int
fit_survfit <- survfit(Surv(AVAL, event = EVENT) ~ sld_above_100, data = SLD_100_adtte_PFS)
fit_survfit
km_plot <- ggsurvplot(fit_survfit, 
                      title = "BEP n=87 \nPFS HR=0.76 [0.44-1.32]",
                      legend = "none",
                      main = "BEP n=87", font.main = 16, 
                      pval = TRUE,
                      palette = c("red", "blue"),
                      risk.table = TRUE, risk.table.y.text.col = TRUE, risk.table.col = "strata",
                      xlab = "Months", ylab = "PFS")
km_plot #600x500

# 5.3 OS KM plot of patients ≥100mm SLD vs less than 100mm SLD 
adtte_OS <- adtte %>% 
  filter(PARAMCD == "OS") %>% 
  select(SUBJID, AVAL, CNSR) %>% 
  separate(col=SUBJID, into=c("Study", "Site", "PT"), sep="-")
SLD_100_adtte_OS <- left_join(x = SLD_100, y = adtte_OS, by = "PT") %>%
  mutate(EVENT = CNSR == 0) %>%
  mutate(sld_above_100 = factor(sld_above_100, levels=c(TRUE, FALSE)))
SLD_100_adtte_OS$sld_above_100 = relevel(SLD_100_adtte_OS$sld_above_100, ref = "TRUE")
fit_survdiff <- survdiff(Surv(AVAL, event = EVENT) ~ sld_above_100, data = SLD_100_adtte_OS)
pval <- pchisq(fit_survdiff$chisq, length(fit_survdiff$n)-1, lower.tail = FALSE)
fit_coxmodel <- coxph(Surv(AVAL, event = EVENT) ~ sld_above_100, data = SLD_100_adtte_OS, ties = "exact")
summary(fit_coxmodel)$conf.int
fit_survfit <- survfit(Surv(AVAL, event = EVENT) ~ sld_above_100, data = SLD_100_adtte_OS)
fit_survfit
km_plot <- ggsurvplot(fit_survfit, 
                      title = "BEP n=87 \nPFS HR=0.57 [0.30-1.08]",
                      legend = "none",
                      main = "BEP n=87", font.main = 16, 
                      pval = TRUE,
                      palette = c("red", "blue"),
                      risk.table = TRUE, risk.table.y.text.col = TRUE, risk.table.col = "strata",
                      xlab = "Months", ylab = "OS")
km_plot #600x500

# 5.4 Medians and ranges of cTF at baseline for patients ≥100mm SLD vs less than 100mm SLD
SLD_100 %>%
  group_by(sld_above_100) %>%
  summarise(
    median_value = median(cTF),
    max_value = max(cTF),
    min_value = min(cTF)
  )

# 5.5 Histogram and density plot of cTF at baseline in ITT
# Histogram
ggplot(SLD, aes(x = cTF)) +
  geom_histogram(binwidth = 0.01, fill = "blue", color = "black", alpha = 0.3) +
  labs(title = "Histogram of Baseline cTF (n=87)", x = "Baseline cTF", y = "Number of Patients") +
  theme_classic()

# Density Plot
ggplot(SLD, aes(x = cTF)) +
  geom_density(fill = "blue", alpha = 0.3) +
  labs(title = "Density Plot of Baseline cTF (n=87)", x = "Baseline cTF", y = "Density") +
  theme_classic()

# 6.1 At treatment discontinuation, median cTF with CI
median_cTF <- OSPL %>% filter(VISIT %in% 'TRTDISC')
length(unique(median_cTF$PT)) == nrow(median_cTF)
nrow(median_cTF)
median(median_cTF$cTF)

library(boot)
set.seed(1234)
Mboot = boot(median_cTF$cTF,
             function(x,i) median(x[i]),
             R=100000)
boot.ci(Mboot,
        conf = 0.95,
        type = c("norm", "basic" ,"perc", "bca"))

# 7.1
#################################################################
# 7.1.1 PFS #
#################################################################
adtte_PFS <- adtte %>% 
  filter(PARAMCD == "PFS") %>% 
  select(SUBJID, AVAL, CNSR) %>% 
  separate(col=SUBJID, into=c("Study", "Site", "PT"), sep="-") %>%
  mutate(EVENT = CNSR == 0)
fit_survfit <- survfit(Surv(AVAL, event = EVENT) ~ 1, data = adtte_PFS)
fit_survfit # Need to use Brookmeyer-Crowley method to replace this
km_plot <- ggsurvplot(fit_survfit, 
                      #title = "Kaplan-Meier Plot for Progression-Free Survival per RECIST v1.1 - Investigator (N=87)",
                      legend.title = "Treatment Group",
                      legend.labs = c("All"),
                      legend = "top",
                      main = "BEP n=87", font.main = 16, 
                      conf.int = FALSE,
                      palette = c("blue"),
                      risk.table = TRUE, risk.table.y.text.col = TRUE, risk.table.col = NULL,
                      risk.table.title = NULL, 
                      xlab = "Months", ylab = "Progression-free Survival",
                      break.time.by = 5)
surv_at_median <- summary(fit_survfit, times = 33)$surv
km_plot$plot <- km_plot$plot + 
  geom_segment(aes(x = 0, xend = 33, y = 0.5, yend = 0.5), 
               linetype = "dashed", color = "darkgray", linewidth = 0.7) +
  geom_segment(aes(x = 33, xend = 33, y = 0, yend = surv_at_median), 
               linetype = "dashed", color = "darkgray", linewidth = 0.7) +
  theme(legend.position = c(0.85, 0.85),
        legend.title = element_blank()) 
km_plot #1100*700

# 7.1.2 OS
adtte_OS <- adtte %>% 
  filter(PARAMCD == "OS") %>% 
  select(SUBJID, AVAL, CNSR) %>% 
  separate(col=SUBJID, into=c("Study", "Site", "PT"), sep="-") %>%
  mutate(EVENT = CNSR == 0)
fit_survfit <- survfit(Surv(AVAL, event = EVENT) ~ 1, data = adtte_OS)
fit_survfit
km_plot <- ggsurvplot(fit_survfit, 
                      #title = "Kaplan-Meier Plot for Overall Survival (N=87)",
                      legend.title = "Treatment Group",
                      legend.labs = c("All"),
                      legend = "top",
                      main = "BEP n=87", font.main = 16, 
                      conf.int = FALSE,
                      palette = c("blue"),
                      risk.table = TRUE, risk.table.y.text.col = TRUE, risk.table.col = NULL,
                      risk.table.title = NULL,
                      xlab = "Months", ylab = "Overall Survival",
                      break.time.by = 5)
km_plot$plot <- km_plot$plot + 
  theme(legend.position = c(0.85, 0.85),
        legend.title = element_blank()) 
km_plot #1100*700


# 8.1 Forest Plot
# Load necessary package
library(forestplot)
library(grid)

# Data for the rows, including the top 2 title rows and 8 bottom rows of information
table_data <- rbind(
  c("cTF cut off", " ", "cTF ≥ cut off", " ", " ", "cTF < cut off", " ", " ", " "),
  c(" ", "n", "Events", "Median", "n", "Events", "Median", "Hazard Ratio", 
    "95% CI"),
  c("0.5%", "65", "45", "22.2", "22", "13", "41.5", "0.62", "(0.34, 1.16)"),
  c("1%", "47", "34", "18.4", "40", "24", "41.5", "0.61", "(0.36, 1.03)"),
  c("1.5%", "37", "29", "14.8", "50", "29", "42.0", "0.48", "(0.29, 0.81)"),
  c("2%", "30", "24", "16.5", "57", "34", "42.0", "0.48", "(0.29, 0.82)"),
  c("2.5%", "27", "22", "14.6", "60", "36", "40.5", "0.42", "(0.25, 0.73)"),
  c("5%", "19", "15", "13.7", "68", "43", "38.9", "0.43", "(0.23, 0.78)"),
  c("10%", "13", "11", "12.7", "74", "47", "36.9", "0.39", "(0.20, 0.76)"),
  c("15%", "7", "7", "13.7", "80", "51", "36.8", "0.22", "(0.10, 0.52)")
)

# Hazard ratios and CIs for the plot
hr_values <- c(NA, NA, 0.62, 0.61, 0.48, 0.48, 0.42, 0.43, 0.39, 0.22)
lower_ci <- c(NA, NA, 0.34, 0.36, 0.29, 0.29, 0.25, 0.23, 0.20, 0.10)
upper_ci <- c(NA, NA, 1.16, 1.03, 0.81, 0.82, 0.73, 0.78, 0.76, 0.52)

# Create the forest plot
forestplot(
  labeltext = table_data,
  mean = hr_values, 
  lower = lower_ci, 
  upper = upper_ci,
  is.summary = c(TRUE, TRUE, rep(FALSE, 8)),
  xlab = "Hazard Ratio",
  col = fpColors(box = "royalblue", lines = "darkblue", zero = "black"),
  zero = 1,  # Reference line at HR = 1
  boxsize = 0.2,  # Size of the boxes
  lineheight = unit(1, "cm"),  # Adjust line height for rows
  xticks = c(0.1, 0.5, 1, 1.5, 2),  # Define the ticks on the x-axis
  hrzl_lines = list(
    "3" = gpar(lty = 1, columns = 1:9)  # Solid line between row 2 and 3 across all columns
  ),
  clip = c(0.1, 1.2),  # Updated range for the x-axis
  new_page = TRUE,  # Ensure a new page for the plot
  vert_lines = list(
    "3" = gpar(lty = 1, columns = 3:4)  # Solid line between column 3 and 4
  ),
  align = "c",  # Center align the text in each column
  cex.axis = 6  # Increase the font size of the x-axis labels
)

# 9.1 OncoPrint at End of Treatment
library(dplyr)
library(tidyr)
library(ComplexHeatmap)

# 1. Input alteration (OAPL) and response (adrs) data 
OAPL <- read.csv("~/ctdna/data_cleaning_20221219/cohortA/bfast_oapl_cohortA.withDBIannotations.clean.csv"
                 ,stringsAsFactors = TRUE)

# 2. Data manipulation
total_samples <- OSPL %>% filter(VISIT %in% 'TRTDISC')
OAPL_SCRN <- OAPL %>%
  filter(VISIT %in% 'SCRN')
OAPL_TRTDISC <- OAPL %>%
  filter(VISIT %in% 'TRTDISC')
tmp <- anti_join(OAPL_TRTDISC, OAPL_SCRN, by = "pt_mutID")

# Restrci to 'known' and 'likely'; top 20 genes 
OAPL_oncoprint <- tmp %>% 
  filter(variant_QC_pass, 
         somatic_status_functional_impact %in% c('known', 'likely'),
         #bep_62gene,
         VISIT %in% 'TRTDISC') %>%
  mutate(Alterations = paste0(alteration_type, "_", coding_type)) %>%
  select(PT, gene, Alterations) %>%
  distinct() %>%
  pivot_wider(names_from = gene, values_from = Alterations, values_fn = list)
OAPL_oncoprint <- data.frame(OAPL_oncoprint)
rownames(OAPL_oncoprint) <- OAPL_oncoprint$PT
OAPL_oncoprint <- OAPL_oncoprint[, -1]

OAPL_oncoprint_mat <- matrix("", nrow = nrow(OAPL_oncoprint), ncol = ncol(OAPL_oncoprint))
rownames(OAPL_oncoprint_mat) <- rownames(OAPL_oncoprint)
colnames(OAPL_oncoprint_mat) <- colnames(OAPL_oncoprint)
for (i in 1:nrow(OAPL_oncoprint)) {
  for (j in 1:ncol(OAPL_oncoprint)) {
    if (!is.null(OAPL_oncoprint[[j]][[i]])) {
      if (length(OAPL_oncoprint[[j]][[i]])==1) {
        OAPL_oncoprint_mat[i,j] <- OAPL_oncoprint[[j]][[i]]
      } else {
        OAPL_oncoprint_mat[i,j] <- paste(OAPL_oncoprint[[j]][[i]],collapse=";")
      }
    }
  }
}
OAPL_oncoprint_mat <- t(OAPL_oncoprint_mat)

# 3. OncoPrint Plot
table(paste0(OAPL$alteration_type, "_", OAPL$coding_type))
table(OAPL_oncoprint$Alterations)

col = c("CN_amplification" = "#EDC948", "RE_rearrangement" = "red", "SV_frameshift" = "#008000",
        "SV_missense" = "blue", "SV_nonframeshift" = "#79706E", "SV_nonsense" = "#9D7660",
        "SV_nonstop" = "#91DCEA", "SV_promoter" = "#e39802", "SV_splice" = "#bfbb60")
alter_fun = list(
  background = alter_graphic("rect", fill = "#CCCCCC"),   
  CN_amplification  = alter_graphic("rect", fill = "#EDC948"), #
  RE_rearrangement = alter_graphic("rect", fill = "red"), #
  SV_frameshift = alter_graphic("rect", height = 0.33, fill = "#008000"),
  SV_missense = alter_graphic("rect", height = 0.33, fill = "blue"), #
  SV_nonframeshift = alter_graphic("rect", fill = "#79706E"),
  SV_nonsense = alter_graphic("rect", fill = "#9D7660"),
  SV_nonstop = alter_graphic("rect", fill = "#91DCEA"),
  SV_promoter = alter_graphic("rect", fill = "#e39802"),
  SV_splice = alter_graphic("rect", fill = "#bfbb60")
)

heatmap_legend_param = list(title = "Alternations", 
                            at = c("CN_amplification", "RE_rearrangement", "SV_frameshift", "SV_missense",
                                   "SV_nonframeshift", "SV_nonsense", "SV_nonstop", "SV_promoter", "SV_splice"), 
                            labels = c("CN_amplification", "RE_rearrangement", "SV_frameshift", "SV_missense",
                                       "SV_nonframeshift", "SV_nonsense", "SV_nonstop", "SV_promoter", "SV_splice"))
oncoPrint(OAPL_oncoprint_mat,
          alter_fun = alter_fun, col = col,
          remove_empty_columns = TRUE, remove_empty_rows = TRUE,
          heatmap_legend_param = heatmap_legend_param)

# 9.2 bep 62 gene
OAPL_oncoprint <- tmp %>% 
  filter(variant_QC_pass, 
         somatic_status_functional_impact %in% c('known', 'likely'),
         bep_62gene,
         VISIT %in% 'TRTDISC') %>%
  mutate(Alterations = paste0(alteration_type, "_", coding_type)) %>%
  select(PT, gene, Alterations) %>%
  distinct() %>%
  pivot_wider(names_from = gene, values_from = Alterations, values_fn = list)
OAPL_oncoprint <- data.frame(OAPL_oncoprint)
rownames(OAPL_oncoprint) <- OAPL_oncoprint$PT
OAPL_oncoprint <- OAPL_oncoprint[, -1]

OAPL_oncoprint_mat <- matrix("", nrow = nrow(OAPL_oncoprint), ncol = ncol(OAPL_oncoprint))
rownames(OAPL_oncoprint_mat) <- rownames(OAPL_oncoprint)
colnames(OAPL_oncoprint_mat) <- colnames(OAPL_oncoprint)
for (i in 1:nrow(OAPL_oncoprint)) {
  for (j in 1:ncol(OAPL_oncoprint)) {
    if (!is.null(OAPL_oncoprint[[j]][[i]])) {
      if (length(OAPL_oncoprint[[j]][[i]])==1) {
        OAPL_oncoprint_mat[i,j] <- OAPL_oncoprint[[j]][[i]]
      } else {
        OAPL_oncoprint_mat[i,j] <- paste(OAPL_oncoprint[[j]][[i]],collapse=";")
      }
    }
  }
}
OAPL_oncoprint_mat <- t(OAPL_oncoprint_mat)

# 3. OncoPrint Plot
table(paste0(OAPL$alteration_type, "_", OAPL$coding_type))
table(OAPL_oncoprint$Alterations)

col = c("CN_amplification" = "#EDC948", "RE_rearrangement" = "red", "SV_frameshift" = "#008000",
        "SV_missense" = "blue", "SV_nonframeshift" = "#79706E", "SV_nonsense" = "#9D7660",
        "SV_nonstop" = "#91DCEA", "SV_promoter" = "#e39802", "SV_splice" = "#bfbb60")
alter_fun = list(
  background = alter_graphic("rect", fill = "#CCCCCC"),   
  CN_amplification  = alter_graphic("rect", fill = "#EDC948"), #
  RE_rearrangement = alter_graphic("rect", fill = "red"), #
  SV_frameshift = alter_graphic("rect", height = 0.33, fill = "#008000"),
  SV_missense = alter_graphic("rect", height = 0.33, fill = "blue"), #
  SV_nonframeshift = alter_graphic("rect", fill = "#79706E"),
  SV_nonsense = alter_graphic("rect", fill = "#9D7660"),
  SV_nonstop = alter_graphic("rect", fill = "#91DCEA"),
  SV_promoter = alter_graphic("rect", fill = "#e39802"),
  SV_splice = alter_graphic("rect", fill = "#bfbb60")
)

heatmap_legend_param = list(title = "Alternations", 
                            at = c("CN_amplification", "RE_rearrangement", "SV_frameshift", "SV_missense",
                                   "SV_nonframeshift", "SV_nonsense", "SV_nonstop", "SV_promoter", "SV_splice"), 
                            labels = c("CN_amplification", "RE_rearrangement", "SV_frameshift", "SV_missense",
                                       "SV_nonframeshift", "SV_nonsense", "SV_nonstop", "SV_promoter", "SV_splice"))
oncoPrint(OAPL_oncoprint_mat,
          alter_fun = alter_fun, col = col,
          remove_empty_columns = TRUE, remove_empty_rows = TRUE,
          heatmap_legend_param = heatmap_legend_param)

# 10.1 Density plot of baseline cTF
library(scales)
library(ggplot2)
library(dplyr)

base_cTF <- OSPL %>%
  filter(VISIT %in% 'SCRN') %>%
  select(cTF)

ggplot(base_cTF, aes(x = cTF)) + 
  geom_density(fill = "lightblue", alpha = 0.6) +
  labs(x = "Baseline cTF", y = "Number of Patients") +
  scale_x_continuous(labels = label_percent(accuracy = 1), 
                     breaks = pretty_breaks(n = 15)) +  # More tick marks with 'pretty_breaks'
  theme_classic()

