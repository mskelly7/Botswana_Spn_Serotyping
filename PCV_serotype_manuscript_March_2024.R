##### PCV Serotype and Antibiotic Susceptibility Code 03202024 ####
## Jillian Hurst, PhD ##
## jillian.hurst@duke.edu ##

set.seed(1234)
library(tidyverse)
library(rstatix)
library(ggplot2)
library(ggpubr)
library(furniture)
library(dplyr)
library(aod)
library(coin)
library(RColorBrewer)
library(ggthemes)
library(data.table)
library(Kendall)
library(scales)
library(data.table)
library(cowplot)
library(gridGraphics)
library(ggbreak)
library(ggrepel)


#Load Dataset
data <- read_csv("BMC_serotyping_dataset_112223.csv")
spec(data)

#Create a cohort that only includes the infant samples
cohort <- subset(data, Source == "INF") 
nrow(cohort) #there were 264 samples included in the cohort

#Number of unique sampleIDs
length(unique(cohort$SampleID)) #There is one sample with two isolates; 263 unique infants

#Group strains into vaccine serotypes
typeof(cohort$Serotype)
pcv13 <- list("1", "3", "4", "5", "6A", "6B", "6C", "7F", "9V", "14", "18C", "19A", "19F", "23F")
pcv15 <- list("22F", "33F")
pcv20 <- list("8", "10A", "11A", "12F", "15B")
cohort <- cohort %>% mutate(spneumo_type = case_when(
  (cohort$Serotype %in% pcv13) ~ "PCV13",
  (cohort$Serotype %in% pcv15) ~ "PCV15",
  (cohort$Serotype %in% pcv20) ~ "PCV20",
  .default = "NVT"))

cohort <- cohort %>% mutate(colonized_PCV13 = case_when(
  spneumo_type == "PCV13" ~ 1,
  .default = 0))
cohort <- cohort %>% mutate(colonized_PCV15 = case_when(
  spneumo_type == "PCV15" ~ 1,
  .default = 0))
cohort <- cohort %>% mutate(colonized_PCV20 = case_when(
  spneumo_type == "PCV20" ~ 1,
  .default = 0))
cohort <- cohort %>% mutate(colonized_NVT = case_when(
  spneumo_type == "NVT" ~ 1,
  .default = 0))

#Evaluate lytA density in PCV13 and non-PCV13 colonizing isolates
group_by(cohort, colonized_PCV13) %>%
  summarise(
    count = n(),
    median = median(inf_sp, na.rm = TRUE),
    IQR = IQR(inf_sp, na.rm = TRUE))
PCV13 <- cohort[cohort$colonized_PCV13 == 1,]
summary(PCV13$inf_sp)
Not_PCV13 <- cohort[cohort$colonized_PCV13 == 0,]
summary(Not_PCV13$inf_sp)

ggboxplot(cohort, x = "colonized_PCV13", y = "inf_sp", 
          color = "colonized_PCV13", palette = c("#00AFBB", "#E7B800"),
          ylab = "lytA", xlab = "Groups")

lytA_diff <- wilcox.test(inf_sp ~ colonized_PCV13, data = cohort,
                   exact = FALSE)
lytA_diff

##Summarize immutable characteristics of the cohort at the participant level (n=150)
##Include: Sex, Maternal HIV status, Birthweight, Low birthweight status
#Start by creating a table where infants have one entry only
colnames(cohort)
patient_level <- cohort %>%
  group_by(StudyID) %>%
  arrange(age_days) %>%
  filter(row_number()==1)
nrow(patient_level) #Confirm 150 unique infants in the cohort

#Table 1 Sex
table(patient_level$sex)
patient_level %>% 
  group_by(sex) %>% 
  summarise(Percentage= 100 * n()/nrow(.))

#Table 1 maternal HIV status
table(patient_level$mat_hiv)
patient_level %>% 
  group_by(mat_hiv) %>% 
  summarise( percent = 100 * n() / nrow(.) )

#Birthweight continuous table 1
summary(patient_level$bw)

#Table 1 low birth weight
table(patient_level$lbw)
patient_level %>% 
  group_by(lbw) %>% 
  summarise( percent = 100 * n() / nrow( . ) )


##Summarize data at the collection level for Table 1
summary(cohort$age_days)

#Table 1 Breastfeeding status
table(cohort$breastmilk)
cohort %>% 
  group_by(breastmilk) %>% 
  summarise( percent = 100 * n() / nrow( cohort ) )

#Table 1 Residence
table(cohort$residence)
cohort %>% 
  group_by(residence) %>% 
  summarise( percent = 100 * n() / nrow( cohort ) )

#Table 1 Sample collection year
table(cohort$year)
cohort %>% 
  group_by(year) %>% 
  summarise( percent = 100 * n() / nrow( cohort ) )

#Table 1 Sample collection season
table(cohort$season)
cohort %>% 
  group_by(season) %>% 
  summarise( percent = 100 * n() / nrow( cohort ) )

#Create new variable, pcv_doses, to show number of doses received at time of sample collection.
#NB: A dose only counts if it was received at least 14 days before sample collection
cohort <- cohort %>% mutate(pcv_doses = case_when(
  (pcv3 + 14) < age_days ~ 3,
  (pcv2 + 14) < age_days ~ 2,
  (pcv1 + 14) < age_days ~ 1,
  .default = 0))

table(cohort$pcv_doses)
cohort %>% 
  group_by(pcv_doses) %>% 
  summarise( percent = 100 * n() / nrow( cohort ) )
summary(cohort$pcv_doses)

no_doses_year <- cohort[cohort$pcv_doses == 0,]
table(no_doses_year$year)

#Create a new variable for days since last PCV13 dose
dose_3 <- cohort[cohort$pcv_doses == 3,]
dose_3$days_since_last_pcv <- dose_3$age_days - dose_3$pcv3
dose_3 <- dose_3[c("SampleID", "days_since_last_pcv")]
dose_2 <- cohort[cohort$pcv_doses == 2,]
dose_2$days_since_last_pcv <- dose_2$age_days - dose_2$pcv2
dose_2 <- dose_2[c("SampleID", "days_since_last_pcv")]
dose_1 <- cohort[cohort$pcv_doses == 1,]
dose_1$days_since_last_pcv <- dose_1$age_days - dose_1$pcv1
dose_1 <- dose_1[c("SampleID", "days_since_last_pcv")]
days_since_dose <- rbind(dose_3, dose_2, dose_1)
view(days_since_dose)
days_since_dose <- days_since_dose %>%
  group_by(SampleID) %>%
  slice_max(order_by = days_since_last_pcv, n = 1)
view(days_since_dose)
cohort <- merge(x=cohort, y=days_since_dose, by="SampleID", all.x = TRUE)
#Find and remove the duplicated rows
cohort <- cohort %>% distinct()

#Table 1 antibiotics 
table(cohort$inf_abx_any)
cohort %>% 
  group_by(inf_abx_any) %>% 
  summarise( percent = 100 * n() / nrow( cohort ) )
#Table 1 amoxicillin exposure
table(cohort$inf_abx_amox)
cohort %>% 
  group_by(inf_abx_amox) %>% 
  summarise( percent = 100 * n() / nrow( cohort ) )
#Table 1 TS exposure
table(cohort$inf_abx_cotrim)
cohort %>% 
  group_by(inf_abx_cotrim) %>% 
  summarise( percent = 100 * n() / nrow( cohort ) )
#Table 1 metronidazole exposure
table(cohort$inf_abx_metro)
cohort %>% 
  group_by(inf_abx_metro) %>% 
  summarise( percent = 100 * n() / nrow( cohort ) )

#Create a new variable to show if the infant has any respiratory virus diagnosis at the time of samples collection
cohort <- cohort %>% mutate(virus_any = case_when(
  cohort$inf_rv == "NEG" ~ "N",
  .default = "Y"))

#Table 1 Respiratory Virus infection
table(cohort$virus_any)
cohort %>% 
  group_by(virus_any) %>% 
  summarise(Percentage= 100 * n()/nrow(.))

#Get a list of values for the respiratory virus column
unique(cohort$inf_rv)
table(cohort$inf_rv)
#Top values are adenovirus (n=11) and enterovirus/rhinovirus (n=76). All other viruses have 5 or fewer occurrences and will be lumped together as other
cohort %>% 
  group_by(inf_rv) %>% 
  summarise( percent = 100 * n() / nrow( cohort ) )

#Create a new variable to lump the virus types
ER <- list("ER")
ADENO <- list("ADENO")
NEG <- list("NEG", NA)
cohort <- cohort %>% mutate(virus_type = case_when(
  (cohort$inf_rv %in% ER) ~ "ER",
  (cohort$inf_rv %in% ADENO) ~ "ADENO",
  (cohort$inf_rv %in% NEG) ~ "NEG",
  .default = "Other virus"))
table(cohort$virus_type)
cohort %>% 
  group_by(virus_type) %>% 
  summarise( percent = 100 * n() / nrow( cohort ) )

#Create a categorical variable that indicates infant month by age at time of sample collection
#Not included in table
levels( factor( cohort$month ) )
table(cohort$month)
cohort %>% 
  group_by(month) %>% 
  summarise(Percentage= 100 * n()/nrow(.))

#Identify infants from whom more than one pneumococcal isolate was obtained
repeat_spn <- cohort %>% group_by(StudyID) %>% filter(n()>1) 
unique(repeat_spn$StudyID)
length(unique(repeat_spn$StudyID))
table(repeat_spn$StudyID)

#Determine if any children in the cohort had more than one isolate identified in a single sample
repeat_sample <- cohort %>% group_by(SampleID) %>% filter(n()>1) 
nrow(repeat_sample) #Only one child had a sample that yielded 2 isolates.
repeat_sample$Serotype #The two serotypes were 15B and 6A

#Check the lytA density by number of PCV3 doses received (split between no doses and any doses, then 0-1 doses and 2+ doses)
lytA_vax <- cohort[,c("inf_sp", "pcv_doses")]
lytA_vax <- lytA_vax %>% mutate(doses_any = case_when(
  cohort$pcv_doses >0 ~ "Y",
  .default = "N"))
lytA_vax <- lytA_vax %>% mutate(doses_2 = case_when(
  cohort$pcv_doses >1 ~ "Y",
  .default = "N"))
Vax_1 <- lytA_vax[lytA_vax$doses_any == "Y",]
summary(Vax_1$inf_sp)
No_vax <- lytA_vax[lytA_vax$doses_any == "N",]
summary(No_vax$inf_sp)
group_by(lytA_vax, doses_any) %>%
  summarise(
    count = n(),
    median = median(inf_sp, na.rm = TRUE),
    IQR = IQR(inf_sp, na.rm = TRUE))
ggboxplot(lytA_vax, x = "doses_any", y = "inf_sp", 
          color = "doses_any", palette = c("#00AFBB", "#E7B800"),
          ylab = "lytA", xlab = "Groups")
lytA_vax_dose <- wilcox.test(inf_sp ~ doses_any, data = lytA_vax,
                    exact = FALSE)
lytA_vax_dose

group_by(lytA_vax, doses_2) %>%
  summarise(
    count = n(),
    median = median(inf_sp, na.rm = TRUE),
    IQR = IQR(inf_sp, na.rm = TRUE))
Vax_2 <- lytA_vax[lytA_vax$doses_2 == "Y",]
summary(Vax_2$inf_sp)
Vax0_1 <- lytA_vax[lytA_vax$doses_2 == "N",]
summary(Vax0_1$inf_sp)
ggboxplot(lytA_vax, x = "doses_2", y = "inf_sp", 
          color = "doses_2", palette = c("#00AFBB", "#E7B800"),
          ylab = "lytA", xlab = "Groups")
lytA_2doses <- wilcox.test(inf_sp ~ doses_2, data = lytA_vax,
                    exact = FALSE)
lytA_2doses

#Run a Cochrane-Armitage test to see if there is a significant change in the proportion of the different
#vaccine serotype groups over time
year <- cohort$year
spneumo_type <- cohort$spneumo_type
aim1_spneumo_type_year <- table(spneumo_type, year)
aim1_spneumo_type_year
sum(aim1_spneumo_type_year)
prop.table(aim1_spneumo_type_year,
           margin = 2)
Test_aim1_spneumotype = chisq_test(aim1_spneumo_type_year,
                          scores = list("year" = c(2016, 2017, 2018, 2019)))
Test_aim1_spneumotype
#Test for changes in PCV13 colonization
year <- cohort$year
colonized_PCV13 <- cohort$colonized_PCV13
aim1_PCV13_year <- table(colonized_PCV13, year)
aim1_PCV13_year
sum(aim1_PCV13_year)
prop.table(aim1_PCV13_year,
           margin = 2)
Test_aim1_PCV13_year = chisq_test(aim1_PCV13_year,
                                   scores = list("year" = c(2016, 2017, 2018, 2019)))
Test_aim1_PCV13_year
#Test for changes in PCV15 colonization
year <- cohort$year
colonized_PCV15 <- cohort$colonized_PCV15
aim1_PCV15_year <- table(colonized_PCV15, year)
aim1_PCV15_year
sum(aim1_PCV15_year)
prop.table(aim1_PCV15_year,
           margin = 2)
Test_aim1_PCV15_year = chisq_test(aim1_PCV15_year,
                                  scores = list("year" = c(2016, 2017, 2018, 2019)))
Test_aim1_PCV15_year
#Test for changes in PCV20 colonization
year <- cohort$year
colonized_PCV20 <- cohort$colonized_PCV20
aim1_PCV20_year <- table(colonized_PCV20, year)
aim1_PCV20_year
sum(aim1_PCV20_year)
prop.table(aim1_PCV20_year,
           margin = 2)
Test_aim1_PCV20_year = chisq_test(aim1_PCV20_year,
                                  scores = list("year" = c(2016, 2017, 2018, 2019)))
Test_aim1_PCV20_year
#Test for changes in NVT colonization
year <- cohort$year
colonized_NVT <- cohort$colonized_NVT
aim1_NVT_year <- table(colonized_NVT, year)
aim1_NVT_year
sum(aim1_NVT_year)
prop.table(aim1_NVT_year,
           margin = 2)
Test_aim1_NVT_year = chisq_test(aim1_NVT_year,
                                  scores = list("year" = c(2016, 2017, 2018, 2019)))
Test_aim1_NVT_year


#Figure 1 showing proportions of more of the individual serotypes by year
Serotype_summary <- table(cohort$Serotype, cohort$year)
rowSums(Serotype_summary)
spneumo_type_summary <- table(cohort$spneumo_type, cohort$year)
rowSums(spneumo_type_summary)

#Show all PCV13 serotypes and all non-vaccine serotypes with 5 or more detections: Non-vaccine (7C, 9N, 15A, 15C, 16F, 17F, 21, 23A, 23B, 34, 35B, 35F)#
#PCV20 show 11A, 15B and additional PCV20 (8, 10A, 12F)
#PCV15 show 22F, 33F 
#Everything else will be classified as non-vaccine serotypes (11D, 13, 15F, 18F, 20, 24, 28F, 29, 31, 33A, 35, 35A, 40, 44, 7A, 7F, 9L, NT)
SuppTab1 <- cohort[c("year", "Serotype")]
SuppTab1$sero_class[SuppTab1$Serotype == "1"] <- "1"
SuppTab1$sero_class[SuppTab1$Serotype == "3"] <- "3"
SuppTab1$sero_class[SuppTab1$Serotype == "4"] <- "4"
SuppTab1$sero_class[SuppTab1$Serotype == "5"] <- "5"
SuppTab1$sero_class[SuppTab1$Serotype == "6A"] <- "6A"
SuppTab1$sero_class[SuppTab1$Serotype == "6B"] <- "6B"
SuppTab1$sero_class[SuppTab1$Serotype == "6C"] <- "6C"
SuppTab1$sero_class[SuppTab1$Serotype == "7F"] <- "7F"
SuppTab1$sero_class[SuppTab1$Serotype == "9V"] <- "9V"
SuppTab1$sero_class[SuppTab1$Serotype == "14"] <- "14"
SuppTab1$sero_class[SuppTab1$Serotype == "18C"] <- "18C"
SuppTab1$sero_class[SuppTab1$Serotype == "19A"] <- "19A"
SuppTab1$sero_class[SuppTab1$Serotype == "19F"] <- "19F"
SuppTab1$sero_class[SuppTab1$Serotype == "23F"] <- "23F"
#Additional PCV15 serotypes
SuppTab1$sero_class[SuppTab1$Serotype == "22F"] <- "22F"
SuppTab1$sero_class[SuppTab1$Serotype == "33F"] <- "33F"
#Additional PCV20 serotypes
SuppTab1$sero_class[SuppTab1$Serotype == "11A"] <- "11A"
SuppTab1$sero_class[SuppTab1$Serotype == "15B"] <- "15B"
SuppTab1$sero_class[SuppTab1$Serotype == "8" | SuppTab1$Serotype == "10A" | SuppTab1$Serotype == "12F"] <- "Additional PCV-20"
#Highlighted non-vaccine serotypes: 7C, 9N, 15A, 15C, 16F, 17F, 21, 23A, 23B, 34, 35B, 35F
SuppTab1$sero_class[SuppTab1$Serotype == "7C"] <- "7C"
SuppTab1$sero_class[SuppTab1$Serotype == "9N"] <- "9N"
SuppTab1$sero_class[SuppTab1$Serotype == "15A"] <- "15A"
SuppTab1$sero_class[SuppTab1$Serotype == "15C"] <- "15C"
SuppTab1$sero_class[SuppTab1$Serotype == "16F"] <- "16F"
SuppTab1$sero_class[SuppTab1$Serotype == "17F"] <- "17F"
SuppTab1$sero_class[SuppTab1$Serotype == "21"] <- "21"
SuppTab1$sero_class[SuppTab1$Serotype == "23A"] <- "23A"
SuppTab1$sero_class[SuppTab1$Serotype == "23B"] <- "23B"
SuppTab1$sero_class[SuppTab1$Serotype == "34"] <- "34"
SuppTab1$sero_class[SuppTab1$Serotype == "35B"] <- "35B"
SuppTab1$sero_class[SuppTab1$Serotype == "35F"] <- "35F"
#Other non-vaccine serotypes
SuppTab1$sero_class[SuppTab1$Serotype == "11D" | SuppTab1$Serotype == "13" | SuppTab1$Serotype == "15F" | SuppTab1$Serotype == "18F" | SuppTab1$Serotype == "20" | SuppTab1$Serotype == "24" | SuppTab1$Serotype == "28F" | SuppTab1$Serotype == "29" | SuppTab1$Serotype == "31" | SuppTab1$Serotype == "33A" | SuppTab1$Serotype == "35" | SuppTab1$Serotype == "35A" | SuppTab1$Serotype == "40" | SuppTab1$Serotype == "44" | SuppTab1$Serotype == "7A" | SuppTab1$Serotype == "9L" | SuppTab1$Serotype == "NT"] <- "Other non-vaccine"
view(SuppTab1)
nrow(SuppTab1) #there are 264 data points for sero_class
SuppTab1_total <- table(SuppTab1$sero_class)
SuppTab1_total <- as.data.frame(SuppTab1_total)
colnames(SuppTab1_total) <- c("serotype", "frequency")
SuppTab1_total$prop <- SuppTab1_total$frequency/264
view(SuppTab1_total)

PCV13_serotypes <- dplyr::filter(SuppTab1_total, serotype %in% c("1", "3","4", "5", "6A", "6B", "6C", "7F", "9V", "14", "18C", "19A", "19F", "23F"))
sum(PCV13_serotypes$frequency)
sum(PCV13_serotypes$prop)
Nonvax_serotypes <- dplyr::filter(SuppTab1_total, serotype %in% c("7C", "9N", "15A", "15C", "16F", "17F", "21","23A", "23B", "34", "35B", "35F", "Other non-vaccine"))
sum(Nonvax_serotypes$frequency)
sum(Nonvax_serotypes$prop)
PCV15_serotypes<- dplyr::filter(SuppTab1_total, serotype %in% c("22F", "33F"))
sum(PCV15_serotypes$frequency)
sum(PCV15_serotypes$prop)
PCV20_serotypes<- dplyr::filter(SuppTab1_total, serotype %in% c("11A", "15B", "Additional PCV-20"))
sum(PCV20_serotypes$frequency)
sum(PCV20_serotypes$prop)

SuppTab1_freq <- table(SuppTab1$sero_class, SuppTab1$year)
SuppTab1_freq <- as.data.frame(SuppTab1_freq)
view(SuppTab1_freq)
colnames(SuppTab1_freq) <- c("serotype", "year", "frequency")
view(SuppTab1_freq)

y2016 <- SuppTab1_freq[which(SuppTab1_freq$year == 2016),]
sum(y2016$frequency) #sums to 51
y2016$prop <- y2016$frequency/51
sum(y2016$prop)  #sums to 1
view(y2016)

y2017 <- SuppTab1_freq[which(SuppTab1_freq$year == 2017),]
sum(y2017$frequency) #sums to 92
y2017$prop <- y2017$frequency/92
sum(y2017$prop)  #sums to 1
view(y2017)

y2018 <- SuppTab1_freq[which(SuppTab1_freq$year == 2018),]
sum(y2018$frequency) #sums to 93
y2018$prop <- y2018$frequency/93
sum(y2018$prop)  #sums to 1
view(y2018)

y2019 <- SuppTab1_freq[which(SuppTab1_freq$year == 2019),]
sum(y2019$frequency) #sums to 28
y2019$prop <- y2019$frequency/28
sum(y2019$prop)  #sums to 1
view(y2019)

#Combine proportion info from all years
Fig1_prop <- rbind(y2016, y2017, y2018, y2019)
view(Fig1_prop)
unique(Fig1_prop$serotype)
Fig1_prop_table <- Fig1_prop %>%
  group_by(serotype)
view(Fig1_prop_table)

#Figure 1
#Show the top serotypes with 10 or more detections: PCV13 (6A, 19A, 19F) vs. Non-vaccine (7C, 16F, 21, 23B, 35B)#
#Everything else will be classified as other PCV-13 serotypes, unique PCV-15 serotypes, unique PCV-20 serotypes, or other non-vaccine serotypes
Fig1 <- cohort[c("year", "Serotype")]
typeof(Fig1$Serotype)
Fig1$sero_class <- Fig1$Serotype

#Highlighted PCV13 = 6A, 19A, 19F
Fig1$sero_class[Fig1$Serotype == "6A"] <- "6A"
Fig1$sero_class[Fig1$Serotype == "19A"] <- "19A"
Fig1$sero_class[Fig1$Serotype == "19F"] <- "19F"
#Other PCV 13
Fig1$sero_class[Fig1$Serotype == "1" | Fig1$Serotype == "3" | Fig1$Serotype == "4" | Fig1$Serotype == "5" | Fig1$Serotype == "6B" | Fig1$Serotype == "6C" | Fig1$Serotype == "7F" | Fig1$Serotype == "9V" | Fig1$Serotype == "14" | Fig1$Serotype == "18C"|
                  Fig1$Serotype == "23F"] <- "Other PCV-13"
#Unique PCV15 serotypes
Fig1$sero_class[Fig1$Serotype == "22F" | Fig1$Serotype == "33F"] <- "Additional PCV-15"
#Unique PCV20 serotypes
Fig1$sero_class[Fig1$Serotype == "8" | Fig1$Serotype == "10A" | Fig1$Serotype == "11A" |Fig1$Serotype == "12F" | Fig1$Serotype == "15B"] <- "Additional PCV-20"
#Highlighted non-vaccine serotypes: 7C, 16F, 21, 23B, 35B
Fig1$sero_class[Fig1$Serotype == "7C"] <- "7C"
Fig1$sero_class[Fig1$Serotype == "16F"] <- "16F"
Fig1$sero_class[Fig1$Serotype == "21"] <- "21"
Fig1$sero_class[Fig1$Serotype == "23B"] <- "23B"
Fig1$sero_class[Fig1$Serotype == "35B"] <- "35B"
#Other non-vaccine serotypes
Fig1$sero_class[Fig1$Serotype == "11D" | Fig1$Serotype == "13" | Fig1$Serotype == "15A" | Fig1$Serotype == "15C" | Fig1$Serotype == "15F" | Fig1$Serotype == "17F" | Fig1$Serotype == "18F" | Fig1$Serotype == "20" | Fig1$Serotype == "23A" | Fig1$Serotype == "24" | Fig1$Serotype == "28F" | Fig1$Serotype == "29" | Fig1$Serotype == "31" | Fig1$Serotype == "33A" | Fig1$Serotype == "34" | Fig1$Serotype == "35" | Fig1$Serotype == "35A" | Fig1$Serotype == "35F" | Fig1$Serotype == "40" | Fig1$Serotype == "44" | Fig1$Serotype == "7A" | Fig1$Serotype == "9L" | Fig1$Serotype == "9N" | Fig1$Serotype == "NT"] <- "Other non-vaccine"
view(Fig1)
nrow(Fig1) #there are 264 data points for sero_class
Fig1_total <- table(Fig1$sero_class)
Fig1_total <- as.data.frame(Fig1_total)
colnames(Fig1_total) <- c("serotype", "frequency")
Fig1_total$prop <- Fig1_total$frequency/264
view(Fig1_total)
PCV13_serotypes <- dplyr::filter(Fig1_total, serotype %in% c("6A", "19A", "19F", "Other PCV-13"))
sum(PCV13_serotypes$frequency)
sum(PCV13_serotypes$prop)
Nonvax_serotypes <- dplyr::filter(Fig1_total, serotype %in% c("7C", "16F", "21", "23B", "35B", "Other non-vaccine"))
sum(Nonvax_serotypes$frequency)
sum(Nonvax_serotypes$prop)

Fig1_freq <- table(Fig1$sero_class, Fig1$year)
Fig1_freq <- as.data.frame(Fig1_freq)
view(Fig1_freq)
colnames(Fig1_freq) <- c("serotype", "year", "frequency")
Fig1_freq_year <- table(Fig1_freq$serotype, Fig1_freq$year, Fig1_freq$frequency)    
view(Fig1_freq_year)                          

y2016 <- Fig1_freq[which(Fig1_freq$year == 2016),]
sum(y2016$frequency) #sums to 51
y2016$prop <- y2016$frequency/51
sum(y2016$prop)  #sums to 1
view(y2016)

y2017 <- Fig1_freq[which(Fig1_freq$year == 2017),]
sum(y2017$frequency) #sums to 92
y2017$prop <- y2017$frequency/92
sum(y2017$prop)  #sums to 1
view(y2017)

y2018 <- Fig1_freq[which(Fig1_freq$year == 2018),]
sum(y2018$frequency) #sums to 93
y2018$prop <- y2018$frequency/93
sum(y2018$prop)  #sums to 1
view(y2018)

y2019 <- Fig1_freq[which(Fig1_freq$year == 2019),]
sum(y2019$frequency) #sums to 28
y2019$prop <- y2019$frequency/28
sum(y2019$prop)  #sums to 1
view(y2019)

#Combine proportion info from all years
Fig1_prop <- rbind(y2016, y2017, y2018, y2019)
view(Fig1_prop)
unique(Fig1_prop$serotype)

#Check unique serotype values to make sure we have the correct number of colors
unique(Fig1_prop$serotype)

#Create the figure 1 barchart
#color options: https://www.datanovia.com/en/blog/ggplot-colors-best-tricks-you-will-love/#change-colors-by-groups-ggplot-default-colors)
#annotate: https://ggplot2.tidyverse.org/reference/annotate.html)
#Creating our color palette (we need 12 colors). PCV13 (Other PCV13 and the individual strains), PCV15, PCV20, and non-vaccine strains
#Non-vaccine serotypes will be #3A488A, 7C will be #616CA1; 16F will be #8891B8, 21 will be #B0B5D0, 23B will be #D7DAE7, 35B will be #ebecf3
#PCV20 will be #DABD61; 
#PCV15 will be #D95F30
#Other PCV13 will be #BE3428FF; 6A will be #CB5C52; 19A will be #D8857E; 19F will be #E5ADA9
Fig1_colors <- c("#E5ADA9", "#D8857E", "#cb5c52", "#BE3428FF", "#D95F30", "#DABD61", "#ebecf3", "#D7DAE7", "#B0B5D0", "#8891B8", "#616CA1", "#3A488A")
 
#Creating custom x axis labels with 2 rows (year with n underneath)
xlab <- c("2016\n n=51", "2017\n n=92", "2018\n n=93", "2019\n n=28")

#Revised figure margins version
Fig1_bar <- Fig1_prop %>%
  mutate(serotype = fct_relevel(serotype, "6A", "19A", "19F", "Other PCV-13", "Additional PCV-15", "Additional PCV-20", "7C", "16F", "21", "23B", "35B", "Other non-vaccine")) %>%
  ggplot( aes(fill=serotype, x=year, y=frequency)) + geom_bar(stat= "identity", position = "fill", width = 0.9) +
  scale_fill_manual(values = Fig1_colors) +
  scale_y_continuous(labels = waiver()) +  
  scale_x_discrete(labels=xlab) +
  labs(x="", y="Proportion of isolates", fill = "")+
  theme_classic()+
  theme(plot.title = element_text(size = 26, hjust = 0.5),
    legend.title = NULL,
    legend.text = element_text(size = 26),
    legend.position = "right",
    legend.box.spacing = unit(0, "pt"),
    axis.text = element_text(size = 26),
    axis.title = element_text(size = 26),
    plot.margin = margin(0, 0, 0, 0, "cm")) 
Fig1_bar

###Create pie charts to show relationship between number of vaccine doses received and serotype colonization###
view(cohort)
library(ggrepel)
serotype_pies <- table(cohort$pcv_doses, cohort$spneumo_type)
serotype_pies <- as.data.frame(serotype_pies)
colnames(serotype_pies) <- c("doses", "serotype", "frequency")
serotype_pies$sero_class[serotype_pies$serotype == "NVT"] <- "Non-vaccine"
serotype_pies$sero_class[serotype_pies$serotype == "PCV13"] <- "PCV-13"
serotype_pies$sero_class[serotype_pies$serotype == "PCV15"] <- "Additional PCV-15"
serotype_pies$sero_class[serotype_pies$serotype == "PCV20"] <- "Additional PCV-20"
serotype_pies = subset(serotype_pies, select = -c(serotype))
names(serotype_pies)[names(serotype_pies) == "sero_class"] <- "Serotype"
view(serotype_pies)
pie_colors <- c("#f7f1dd", "#f1c3b2", "#8e99d0", "#e99d96")

#Code for zero dose pie
zero_doses <- serotype_pies[which(serotype_pies$doses == 0),]
sum(zero_doses$frequency) #sums to 51
zero_doses$prop <- zero_doses$frequency/51
sum(zero_doses$prop)  #sums to 1
zero_doses$percent = zero_doses$prop*100
zero_doses <- as.data.frame(zero_doses)
colnames(zero_doses)
zero_doses <- zero_doses[,c("Serotype", "prop")]
zero_doses$Serotype <- factor(zero_doses$Serotype, levels = c("PCV-13", "Additional PCV-15", "Additional PCV-20", "Non-vaccine"))
zero_dose_pie <- ggplot(zero_doses, aes(x = "", (y = prop), fill = Serotype)) +
  theme_void() +
  geom_col(color = "black") +
  geom_text(aes(x = 1.65, label = scales::percent(prop, accuracy = 1)), position = position_stack(vjust = .5), size = 9)+
  coord_polar(theta = "y") +
  scale_fill_manual(values = c("#e99d96", "#f7f1dd", "#f1c3b2", "#8e99d0"))
zero_pie_no_legend <- zero_dose_pie + theme(legend.position="none", panel.grid.major = element_blank(),
                                            panel.grid.minor = element_blank(),
                                            axis.text.x=element_blank(),
                                            plot.margin = unit(c(-.75,-.75,-.75,-.75),"cm"))
zero_pie_no_legend

one_doses <- serotype_pies[which(serotype_pies$doses == 1),]
sum(one_doses$frequency) #sums to 44
one_doses$prop <- one_doses$frequency/44
sum(one_doses$prop)  #sums to 1
one_doses$percent = one_doses$prop*100
one_doses <- as.data.frame(one_doses)
colnames(one_doses)
one_doses <- one_doses[,c("Serotype", "prop")]
one_doses$Serotype <- factor(one_doses$Serotype, levels = c("PCV-13", "Additional PCV-15", "Additional PCV-20", "Non-vaccine"))
one_dose_pie <- ggplot(one_doses, aes(x = "", (y = prop), fill = Serotype)) +
  theme_void() +
  geom_col(color = "black") +
  geom_text(aes(x = 1.65, label = scales::percent(prop, accuracy = 1)), position = position_stack(vjust = .5), size =9)+
  coord_polar(theta = "y") +
  scale_fill_manual(values = c("#e99d96", "#f7f1dd", "#f1c3b2", "#8e99d0"))
one_pie_no_legend <- one_dose_pie + theme(legend.position="none", panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank(),
                                          axis.text.x=element_blank(),
                                          plot.margin = unit(c(-.75,-.75,-.75,-.75),"cm"))
one_pie_no_legend

two_doses <- serotype_pies[which(serotype_pies$doses == 2),]
sum(two_doses$frequency) #sums to 29
two_doses$prop <- two_doses$frequency/29
sum(two_doses$prop)  #sums to 1
two_doses$percent = two_doses$prop*100
two_doses <- as.data.frame(two_doses)
colnames(two_doses)
two_doses <- two_doses[,c("Serotype", "prop")]
two_doses$Serotype <- factor(two_doses$Serotype, levels = c("PCV-13", "Additional PCV-15", "Additional PCV-20", "Non-vaccine"))
two_dose_pie <- ggplot(two_doses, aes(x = "", (y = prop), fill = Serotype)) +
  theme_void() +
  geom_col(color = "black") +
  geom_text(aes(x = 1.65, label = scales::percent(prop, accuracy = 1)), position = position_stack(vjust = .5), size = 9)+
  coord_polar(theta = "y") +
  scale_fill_manual(values = c("#e99d96", "#f7f1dd", "#f1c3b2", "#8e99d0"))
two_pie_no_legend <- two_dose_pie + theme(legend.position="none", panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank(),
                                          axis.text.x=element_blank(),
                                          plot.margin = unit(c(-.75,-.75,-.75,-.75),"cm"))
two_pie_no_legend

three_doses <- serotype_pies[which(serotype_pies$doses == 3),]
sum(three_doses$frequency) #sums to 140
three_doses$prop <- three_doses$frequency/140
sum(three_doses$prop)  #sums to 1
three_doses$percent = three_doses$prop*100
three_doses <- as.data.frame(three_doses)
colnames(three_doses)
three_doses <- three_doses[,c("Serotype", "prop")]
three_doses$Serotype <- factor(three_doses$Serotype, levels = c("PCV-13", "Additional PCV-15", "Additional PCV-20", "Non-vaccine"))
three_dose_pie <- ggplot(three_doses, aes(x = "", (y = prop), fill = Serotype)) +
  theme_void() +
  geom_col(color = "black") +
  geom_text(aes(x = 1.65, label = scales::percent(prop, accuracy = 1)), position = position_stack(vjust = .5), size = 9)+
  coord_polar(theta = "y") +
  scale_fill_manual(values = c("#e99d96", "#f7f1dd", "#f1c3b2", "#8e99d0"))
three_pie_no_legend <- three_dose_pie + theme(legend.position="none", panel.grid.major = element_blank(),
                                              panel.grid.minor = element_blank(),
                                              axis.text.x=element_blank(),
                                              plot.margin = unit(c(-.75,-.75,-.75,-.75),"cm"))
three_pie_no_legend

#Cowplot assembly of the 4 pies
Pies_grid <- plot_grid(zero_pie_no_legend, NULL, one_pie_no_legend, NULL, NULL, NULL, two_pie_no_legend, NULL, three_pie_no_legend, labels = c("Zero doses", "", "One dose", "", "", "", "Two doses", "", "Three doses"), label_size = 24, label_fontface = "plain", label_x=0.5, label_y=0.95, hjust = 0.5, rel_widths = c(1,-0.4,1,1,-0.4,1,1,-0.4,1), rel_heights = c(1,0,1))
Pies_grid
#Add legend
pie_legend <- get_legend(zero_dose_pie + theme(legend.position = "bottom", legend.title = element_text(size=24, face="bold"), legend.text=element_text(size=24))+ 
                           theme(legend.margin = margin(t = 0, r=0, b=0, l=0)))
Fig1_pies_legend <- plot_grid(Pies_grid, NULL, pie_legend, ncol = 1, rel_heights = c(1, -0.05, 0.1), rel_widths = c(1, 0.2)) 
Fig1_pies_legend

#Cowplot assembly of the Fig1 bar plot and 4 pies, using the method where we place one plot_grid into another
left_side <- plot_grid(Fig1_bar, labels = c('A'), label_size = 20, label_fontface = "bold")
left_side
right_side <-plot_grid(Fig1_pies_legend, labels = c('B'), label_size = 20, label_fontface = "bold", label_x = 0.1)
right_side

Bar_pies <- plot_grid(left_side, NULL, right_side, nrow=1, rel_widths=c(1, -0.2, 1.4))
ggsave("Fig1_bar_pie.png", plot = Bar_pies, width = 25, height = 12, dpi = 800) 


#Aim2: Univariate analysis of all the factors in relation to colonization by a PCV13 serotype
#Use a chi-squared test for categorical variables
#Use Wilcoxon test for continuous variables (assuming non-normal distribution)
#Create a variable to indicate whether the detected strain is covered by PCV13
#Create a variable to indicate whether the detected strain is covered by PCV13

table(cohort$colonized_PCV13)
cohort %>% 
  group_by(colonized_PCV13) %>% 
  summarise( percent = 100 * n() / nrow( cohort ) )

#PCV-13 vs. non-PCV-13 Cochrane Armitage
year <- cohort$year
PCV13_colonized <- cohort$colonized_PCV13
aim2_PCV13_year <- table(PCV13_colonized, year)
aim2_PCV13_year
sum(aim2_PCV13_year)
prop.table(aim2_PCV13_year,
           margin = 2)
Test_aim2_PCV13 = chisq_test(aim2_PCV13_year,
                                   scores = list("year" = c(2016, 2017, 2018, 2019)))
Test_aim2_PCV13

#Univariate analysis of age at time of collection
tapply(cohort$age_days, cohort$colonized_PCV13, summary)
wilcox.test(age_days ~ colonized_PCV13, data = cohort, exact = FALSE)


#Univariate analysis of infant sex
table(cohort$sex, cohort$colonized_PCV13)
chisq.test(cohort$sex, cohort$colonized_PCV13, correct=FALSE)
prop.table(table(cohort$colonized_PCV13, cohort$sex), 1)

#Univariate analysis of maternal HIV status
table(cohort$mat_hiv, cohort$colonized_PCV13)
chisq.test(cohort$mat_hiv, cohort$colonized_PCV13, correct=FALSE)
prop.table(table(cohort$mat_hiv, cohort$colonized_PCV13), 1)

#Univariate analysis of breastfeeding status
table(cohort$colonized_PCV13, cohort$breastmilk)
chisq.test(cohort$colonized_PCV13, cohort$breastmilk, correct=FALSE)
prop.table(table(cohort$colonized_PCV13, cohort$breastmilk), 1)

#Univariate analysis of birthweight
tapply(cohort$bw, cohort$colonized_PCV13, summary)
wilcox.test(bw ~ colonized_PCV13, data = cohort, exact = FALSE)


#Univariate analysis of LBW status
table(cohort$colonized_PCV13, cohort$lbw)
chisq.test(cohort$colonized_PCV13, cohort$lbw, correct=FALSE)
prop.table(table(cohort$colonized_PCV13, cohort$lbw), 1)

#Univariate analysis of residence
table(cohort$colonized_PCV13, cohort$residence)
chisq.test(cohort$colonized_PCV13, cohort$residence, correct=FALSE)
prop.table(table(cohort$colonized_PCV13, cohort$residence), 1)

#Univariate analysis of sample year
table(cohort$colonized_PCV13, cohort$year)
prop.table(table(cohort$colonized_PCV13, cohort$year), 1)

#Create a dummy variable to make year continuous rather than categorical
cohort <- cohort %>% mutate(year_cont = case_when(
  (cohort$year == "2016") ~ 1,
  (cohort$year == "2017") ~ 2,
  (cohort$year == "2018") ~ 3,
  (cohort$year == "2019") ~4))
view(cohort)

tapply(cohort$year_cont, cohort$colonized_PCV13, summary)
wilcox.test(year_cont ~ colonized_PCV13, data = cohort)

#Univariate analysis of season
table(cohort$colonized_PCV13, cohort$season)
chisq.test(cohort$colonized_PCV13, cohort$season, correct=FALSE)
prop.table(table(cohort$colonized_PCV13, cohort$season), 1)

#Univariate analysis of antibiotics since prior visit
table(cohort$colonized_PCV13, cohort$inf_abx_any)
chisq.test(cohort$colonized_PCV13, cohort$inf_abx_any, correct=FALSE)
prop.table(table(cohort$colonized_PCV13, cohort$inf_abx_any), 1)

#Univariate analysis of amoxicillin since prior visit
table(cohort$colonized_PCV13, cohort$inf_abx_amox)
chisq.test(cohort$colonized_PCV13, cohort$inf_abx_amox, correct=FALSE)
prop.table(table(cohort$colonized_PCV13, cohort$inf_abx_amox), 1)

#Univariate analysis of cotrim since prior visit
table(cohort$colonized_PCV13, cohort$inf_abx_cotrim)
chisq.test(cohort$colonized_PCV13, cohort$inf_abx_cotrim, correct=FALSE)
prop.table(table(cohort$colonized_PCV13, cohort$inf_abx_cotrim), 1)

#Univariate analysis of metronidazole since prior visit
table(cohort$colonized_PCV13, cohort$inf_abx_metro)
chisq.test(cohort$colonized_PCV13, cohort$inf_abx_metro, correct=FALSE)
prop.table(table(cohort$colonized_PCV13, cohort$inf_abx_metro), 1)

#Univariate analysis of PCV13 doses prior to sample collection (continuous variable)
class(cohort$pcv_doses)
table(cohort$pcv_doses, cohort$colonized_PCV13)
prop.table(table(cohort$pcv_doses, cohort$colonized_PCV13), 1)
tapply(cohort$pcv_doses, cohort$colonized_PCV13, summary)
wilcox.test(pcv_doses ~ colonized_PCV13, data = cohort)

#Univariate analysis of respiratory virus diagnosis at time of collection
table(cohort$virus_any, cohort$colonized_PCV13)
chisq.test(cohort$virus_any, cohort$colonized_PCV13, correct=FALSE)
prop.table(table(cohort$colonized_PCV13, cohort$virus_any), 1)
unique(cohort$inf_rv)
table(cohort$virus_type, cohort$colonized_PCV13)
prop.table(table(cohort$virus_type, cohort$colonized_PCV13), 2)

#Among children with 3-doses of PCV-13, is there a difference in time since last dose among those colonized with PCV13?
#Subset to only look at kids with 3-doses at time of sampling
three_dose <- cohort[cohort$pcv_doses == 3,]
#Calculate the number of days between dose 3 receipt and sampling
three_dose$days_between_sample_dose3 <- three_dose$age_days-three_dose$pcv3
tapply(three_dose$days_between_sample_dose3, three_dose$colonized_PCV13, summary)
wilcox.test(days_between_sample_dose3 ~ colonized_PCV13, data = three_dose)
#Among children who have received any doses of vaccine, does the time since last dose differ among those colonized with PCV13 serotypes?
any_dose <- cohort[!cohort$pcv_doses == 0,]
tapply(any_dose$days_since_last_pcv, any_dose$colonized_PCV13, summary)
wilcox.test(days_since_last_pcv ~ colonized_PCV13, data = any_dose)


#Aim2 Multivariable Logistic regression - what factors make kids more likely to be colonized with a PVC13 vs. non-PCV13 serotype?
#Following the method at https://stats.oarc.ucla.edu/r/dae/logit-regression/#
#Show column names in the cohort dataframe
names(cohort)

#Step 1: Revise the data table to only include the predictor variables and outcome of interest
aim2MR <- select (cohort, c(SampleID, age_days, year_cont, season, sex, bw, lbw, mat_hiv, residence, breastmilk, current_uri, recent_uri, inf_abx_any, virus_any, pcv_doses, colonized_PCV13))

#Indicate which predictors should be treated as categorical variables
aim2MR$sex <- factor(aim2MR$sex) 
aim2MR$sex <- relevel(aim2MR$sex, ref = "F")
aim2MR$mat_hiv <- factor(aim2MR$mat_hiv)
aim2MR$breastmilk <- factor(aim2MR$breastmilk)
aim2MR$lbw <- factor(aim2MR$lbw)  
aim2MR$residence <- factor(aim2MR$residence)
aim2MR$season <- factor(aim2MR$season)
aim2MR$inf_abx_any <- factor(aim2MR$inf_abx_any)
aim2MR$virus_any <- factor(aim2MR$virus_any)

#Run the binomial glm model
mylogit_aim2 <- glm(colonized_PCV13 ~ age_days + year_cont + season + sex + bw + lbw + mat_hiv + residence + breastmilk + inf_abx_any + virus_any + pcv_doses, data = aim2MR, family = "binomial")
summary(mylogit_aim2)
confint(mylogit_aim2)
exp(coef(mylogit_aim2))
exp(cbind(OR = coef(mylogit_aim2), confint(mylogit_aim2)))


##Aim 3: Describe serotype epidemiology and antibiotic resistance##
##First step is to look the MIC50 for each antibiotic during the course of the study period
Abx_mic <- cohort[c("year", "AC", "PG", "CRO", "TS", "AZ")]

#We will calculate the MIC50 and IQR for each drug.
#We will use a Mann-Kendall Trend test to determine if MIC50 changes with time
#https://finnstats.com/index.php/2021/11/28/time-series-trend-analysis-in-r/#:~:text=Time%20series%20trend%20analysis%2C%20The,assumption%20about%20the%20data's%20normality.

AC_data <- table(Abx_mic$AC, Abx_mic$year)
AC_ts <- ts(AC_data,start = c(2016, 2017, 2018, 2019))
MannKendall(AC_ts)
AC_mic_data <- select(Abx_mic, c("AC","year"))
AC_mic_stats <- AC_mic_data %>%
  group_by(year) %>%
  summarise(across(starts_with('AC'), list(med = median, 
                                            first_quartile = ~quantile(., 0.25), 
                                            third_quartile = ~quantile(., 0.75))))
AC_mic_stats <- AC_mic_stats %>%
  mutate(Antibiotic = "Amoxicillin")
colnames(AC_mic_stats)
names(AC_mic_stats)[names(AC_mic_stats) == "AC_med"] <- "median"
names(AC_mic_stats)[names(AC_mic_stats) == "AC_first_quartile"] <- "1st quartile"
names(AC_mic_stats)[names(AC_mic_stats) == "AC_third_quartile"] <- "3rd quartile"
view(AC_mic_stats)

PG_data <- table(Abx_mic$PG, Abx_mic$year)
PG_ts <- ts(PG_data,start = c(2016, 2017, 2018, 2019))
MannKendall(PG_ts)
PG_data <- select(Abx_mic, c("PG","year"))
PG_mic_stats <- PG_data %>%
  group_by(year) %>%
  summarise(across(starts_with('PG'), list(med = median, 
                                           first_quartile = ~quantile(., 0.25), 
                                           third_quartile = ~quantile(., 0.75))))
PG_mic_stats <- PG_mic_stats %>%
  mutate(Antibiotic = "Penicillin")
colnames(PG_mic_stats)
names(PG_mic_stats)[names(PG_mic_stats) == "PG_med"] <- "median"
names(PG_mic_stats)[names(PG_mic_stats) == "PG_first_quartile"] <- "1st quartile"
names(PG_mic_stats)[names(PG_mic_stats) == "PG_third_quartile"] <- "3rd quartile"
view(PG_mic_stats)

CRO_data <- table(Abx_mic$CRO, Abx_mic$year)
CRO_ts <- ts(CRO_data,start = c(2016, 2017, 2018, 2019))
MannKendall(CRO_ts)
CRO_data <- select(Abx_mic, c("CRO","year"))
CRO_mic_stats <- CRO_data %>%
  group_by(year) %>%
  summarise(across(starts_with('CRO'), list(med = median, 
                                           first_quartile = ~quantile(., 0.25), 
                                           third_quartile = ~quantile(., 0.75))))
CRO_mic_stats <- CRO_mic_stats %>%
  mutate(Antibiotic = "Ceftriaxone")
colnames(CRO_mic_stats)
names(CRO_mic_stats)[names(CRO_mic_stats) == "CRO_med"] <- "median"
names(CRO_mic_stats)[names(CRO_mic_stats) == "CRO_first_quartile"] <- "1st quartile"
names(CRO_mic_stats)[names(CRO_mic_stats) == "CRO_third_quartile"] <- "3rd quartile"
view(CRO_mic_stats)

TS_data <- table(Abx_mic$TS, Abx_mic$year)
TS_ts <- ts(TS_data,start = c(2016, 2017, 2018, 2019))
MannKendall(TS_ts)
TS_data <- select(Abx_mic, c("TS","year"))
TS_mic_stats <- TS_data %>%
  group_by(year) %>%
  summarise(across(starts_with('TS'), list(med = median, 
                                            first_quartile = ~quantile(., 0.25), 
                                            third_quartile = ~quantile(., 0.75))))
TS_mic_stats <- TS_mic_stats %>%
  mutate(Antibiotic = "Trimethoprim-sulfamethoxazole")
colnames(TS_mic_stats)
names(TS_mic_stats)[names(TS_mic_stats) == "TS_med"] <- "median"
names(TS_mic_stats)[names(TS_mic_stats) == "TS_first_quartile"] <- "1st quartile"
names(TS_mic_stats)[names(TS_mic_stats) == "TS_third_quartile"] <- "3rd quartile"
view(TS_mic_stats)

AZ_data <- table(Abx_mic$AZ, Abx_mic$year)
AZ_ts <- ts(AZ_data,start = c(2016, 2017, 2018, 2019))
MannKendall(AZ_ts)
AZ_data <- select(Abx_mic, c("AZ","year"))
AZ_mic_stats <- AZ_data %>%
  group_by(year) %>%
  summarise(across(starts_with('AZ'), list(med = median, 
                                           first_quartile = ~quantile(., 0.25), 
                                           third_quartile = ~quantile(., 0.75))))
AZ_mic_stats <- AZ_mic_stats %>%
  mutate(Antibiotic = "Azithromycin")
colnames(AZ_mic_stats)
names(AZ_mic_stats)[names(AZ_mic_stats) == "AZ_med"] <- "median"
names(AZ_mic_stats)[names(AZ_mic_stats) == "AZ_first_quartile"] <- "1st quartile"
names(AZ_mic_stats)[names(AZ_mic_stats) == "AZ_third_quartile"] <- "3rd quartile"
view(AZ_mic_stats)

#Merge the dataframes by year
Abx_mic_stats <- rbind(AC_mic_stats, PG_mic_stats, CRO_mic_stats, TS_mic_stats, AZ_mic_stats)
view(Abx_mic_stats)

#Create dummy variable "antibiotic" so that we have a dataframe with year, mic, antibiotic
AC_data <- data.frame("Year"=c(cohort$year), "MIC"=c(cohort$AC))
AC_data <- AC_data %>%
  mutate(Antibiotic = "Amoxicillin")

PG_data <- data.frame("Year"=c(cohort$year), "MIC"=c(cohort$PG))
PG_data <- PG_data %>%
  mutate(Antibiotic = "Penicillin")

CRO_data <- data.frame("Year"=c(cohort$year), "MIC"=c(cohort$CRO))
CRO_data <- CRO_data %>%
  mutate(Antibiotic = "Ceftriaxone")

TS_data <- data.frame("Year"=c(cohort$year), "MIC"=c(cohort$TS))
TS_data <- TS_data %>%
  mutate(Antibiotic = "Trimethoprim-sulfamethoxazole")

AZ_data <- data.frame("Year"=c(cohort$year), "MIC"=c(cohort$AZ))
AZ_data <- AZ_data %>%
  mutate(Antibiotic = "Azithromycin")

Abx_data <- rbind(AC_data, PG_data, CRO_data, TS_data, AZ_data)
view(Abx_data)

#Create boxplots for susceptibility for each individual antibiotic
#Amoxicillin plot, S: ≤1.0 µg/mL, I: >1.0 & <4.0 µg/mL; R: ≥4.0 µg/mL 
#Recode values such that any values above the resistance cut off are coded as MIC=4.1
AC_data$MIC[AC_data$MIC > 4] <- 4.0
AC_data$Year = as.factor(AC_data$Year)
AC_ggbox <- ggplot(AC_data, aes(x=Year, y=MIC)) + 
  geom_boxplot(fill="#BE3428FF", color="black")+ 
  geom_jitter(shape=16, position=position_jitter(0.3))+
  theme_minimal()+
  labs(x=NULL)+
  ggtitle("Amoxicillin")+
  theme(plot.title = element_text(size = 18, face="bold", hjust=0.5))+
  theme(axis.text = element_text(size = 18))+
  theme(axis.title = element_text(size = 18))+
  theme(axis.text.x=element_blank())+
  geom_hline(yintercept=1, linetype="dotted", color = "black")+
  geom_hline(yintercept=4, linetype="dashed", color="black")
AC_ggbox

#Azithromycin plot, S : ≤0.5 µg/mL, I: >0.5 & <2.0 µg/mL; R: ≥2.0 µg/mL
#Recode values such that any values above the resistance cut off are coded as MIC=2.1
AZ_data$MIC[AZ_data$MIC > 2] <- 2.0
AZ_data$Year = as.factor(AZ_data$Year)
AZ_ggbox <- ggplot(AZ_data, aes(x=Year, y=MIC))+ 
  geom_boxplot(fill="#D95F30", color="black")+ 
  geom_jitter(shape=16, position=position_jitter(0.3))+
  theme_minimal()+
  labs(x=NULL, y=NULL)+
  ggtitle("Azithromycin")+
  theme(plot.title = element_text(size = 18, face="bold", hjust=0.5))+
  theme(axis.text = element_text(size = 18))+
  theme(axis.title = element_text(size = 18))+
  theme(axis.text.x=element_blank())+
  geom_hline(yintercept=0.5, linetype="dotted", color = "black")+
  geom_hline(yintercept=2, linetype="dashed", color="black")
AZ_ggbox

#Penicillin plot - meningitis
#Meningitis: S: ≤0.06 µg/mL; R: >0.06 µg/mL
#Recode values such that any values above the non-meningitis resistance cut off are coded as MIC=0.06
PGmen_data <- data.frame(PG_data)
PGmen_data$MIC[PG_data$MIC > 0.12] <- 0.12
PGmen_data$Year = as.factor(PGmen_data$Year)
PGmen_ggbox <- ggplot(PGmen_data, aes(x=Year, y=MIC)) + 
  geom_boxplot(fill="#8088b3", color="black")+ 
  ylim(0, 0.125)+
  geom_jitter(shape=16, position=position_jitter(0.3))+
  theme_minimal()+
  ggtitle("Penicillin (meningitis)")+
  theme(plot.title = element_text(size = 18, face="bold", hjust=0.5))+
  theme(axis.text = element_text(size = 18))+
  theme(axis.title = element_text(size = 18))+
  geom_hline(yintercept=0.06, linetype="dotted", color = "black")+
  geom_hline(yintercept=0.12, linetype="dashed", color = "black")
PGmen_ggbox

#Penicillin plot - non-meningitis
#Non-meningitis: S: ≤2.0 µg/mL; I: >4.0 µg/mL; R: ≥8.0 µg/mL
#Recode values such that any values above the non-meningitidis NS cut off are coded as MIC=8.1
PGnonmen_data <- data.frame(PG_data)
PGnonmen_data$MIC[PGnonmen_data$MIC > 8] <- 8.0
PGnonmen_data$Year = as.factor(PGnonmen_data$Year)
PGnonmen_ggbox <- ggplot(PGnonmen_data, aes(x=Year, y=MIC)) + 
  geom_boxplot(fill="#B0B5D0", color="black")+ 
  geom_jitter(shape=16, position=position_jitter(0.3))+
  labs(y=NULL) +
  theme_minimal()+
  ggtitle("Penicillin (non-meningitis)")+
  theme(plot.title = element_text(size = 18, face="bold", hjust=0.5))+
  theme(axis.text = element_text(size = 18))+
  theme(axis.title = element_text(size = 18))+
  geom_hline(yintercept=4, linetype="dotted", color = "black")
PGnonmen_ggbox

#Ceftriaxone plot: S: ≤0.5 µg/mL, I: >0.5 & <2.0 µg/mL; R: ≥2.0 µg/mL 
CRO_data$Year = as.factor(CRO_data$Year)
CRO_data$MIC[CRO_data$MIC > 2] <- 2.0
CRO_ggbox <- ggplot(CRO_data, aes(x=Year, y=MIC)) + 
  geom_boxplot(fill="#DABD61", color="black")+ 
  geom_jitter(shape=16, position=position_jitter(0.3))+
  labs(x=NULL, y=NULL)+
  theme_minimal()+
  ggtitle("Ceftriaxone")+
  theme(plot.title = element_text(size = 18, face = "bold", hjust=0.5))+
  theme(axis.text = element_text(size = 18))+
  theme(axis.title = element_text(size = 18))+
  geom_hline(yintercept=0.5, linetype="dotted", color = "black")+
  geom_hline(yintercept=2.0, linetype="dashed", color = "black")+
    theme(axis.text.x=element_blank())
CRO_ggbox

#Trimethoprim-sulfamethoxazole plot, S: ≤0.5 µg/mL, I: >0.5 & <4.0 µg/mL; R: ≥4.0 µg/mL
TS_data$MIC[TS_data$MIC > 4] <- 4.0
TS_data$Year = as.factor(TS_data$Year)
TS_ggbox <- ggplot(TS_data, aes(x=Year, y=MIC)) + 
    geom_boxplot(fill="#3A488A", color="black")+ 
    geom_jitter(shape=16, position=position_jitter(0.3))+
    theme_minimal()+
  labs(y=NULL)+
  ggtitle("TMP-SMX")+
  theme(plot.title = element_text(size = 18, face="bold", hjust=0.5))+
  theme(axis.text = element_text(size = 18))+
  theme(axis.title = element_text(size = 18))+
  geom_hline(yintercept=0.5, linetype="dotted", color = "black")+
  geom_hline(yintercept=4.0, linetype="dashed", color = "black")
TS_ggbox

#Combine the stacked barcharts above
library(patchwork)
layout <-"
ABC
DEF
"
Abx_box <- 
  AC_ggbox+
  AZ_ggbox+
  CRO_ggbox+ 
  PGmen_ggbox+ 
  PGnonmen_ggbox+ 
  TS_ggbox+
  plot_layout(design = layout)
Abx_box
ggsave("Abx_box.png", plot = Abx_box, width = 11, height = 8, dpi = 1200) 


#####SUPPLEMENTAL TABLE 2#####
#Select rows where isolates are resistant to amoxicillin
view(cohort)
AC_res <- cohort[cohort$AC_cat == 'R',]
nrow(AC_res) #3 isolates are resistant
AC_percent <- (3/264)*100
AC_percent
PG_nonmen <- cohort[cohort$PG_cat_nonmen == 'R',]
nrow(PG_nonmen) #0 isolates are resistant
PG_nonmen_percent <- (0/264)*100
PG_nonmen_percent
PG_men <- cohort[cohort$PG_cat_men == 'R',]
nrow(PG_men) #153 isolates are resistant
PG_men_percent <- (153/264)*100
PG_men_percent
CRO <- cohort[cohort$CRO_cat == 'R',]
nrow(CRO) #4 isolates are resistant
CRO_percent <- (4/264)*100
CRO_percent
TS <- cohort[cohort$TS_cat == 'R',]
nrow(TS) #48 isolates are resistant
TS_percent <- (48/264)*100
TS_percent
AZ <- cohort[cohort$AZ_cat == 'R',]
nrow(AZ) #15 isolates are resistant
AZ_percent <- (15/264)*100
AZ_percent

#Create a dummy variable sero_class so that it's easy to lump different serotypes, similar to what was done in the figures
Serotype_3 <- list("3")
Serotype_4 <- list("4")
Serotype_6A <- list("6A")
Serotype_6B <- list("6B")
Serotype_6C <- list("6C")
Serotype_14 <- list("14")
Serotype_18C <- list("18C")
Serotype_19A <- list("19A")
Serotype_19F <- list("19F")
Serotype_23F <- list("23F")
Other_PCV13 <- list ("1", "5", "7F", "9V")
Serotype_22F <- list("22F")
Serotype_33F <- list("33F")
Serotype_11A <- list("11A")
Serotype_15B <- list("15B")
Additional_PCV20 <- list("8", "10A", "12F")
Serotype_7C <- list("7C")
Serotype_9N <- list("9N")
Serotype_15A <- list("15A")
Serotype_15C <- list("15C")
Serotype_16F <- list("16F")
Serotype_17F <- list("17F")
Serotype_21 <- list("21")
Serotype_23A <- list("23A")
Serotype_23B <- list("23B")
Serotype_35B <- list("35B")
Serotype_35F <- list("35F")

cohort <- cohort %>% mutate(sero_class = case_when(
  (cohort$Serotype %in% Serotype_3) ~ "3",
  (cohort$Serotype %in% Serotype_4) ~ "4",
  (cohort$Serotype %in% Serotype_6A) ~ "6A",
  (cohort$Serotype %in% Serotype_6B) ~ "6B",
  (cohort$Serotype %in% Serotype_6C) ~ "6C",
  (cohort$Serotype %in% Serotype_14) ~ "14",
  (cohort$Serotype %in% Serotype_18C) ~ "18C",
  (cohort$Serotype %in% Serotype_19A) ~ "19A",
  (cohort$Serotype %in% Serotype_19F) ~ "19F",
  (cohort$Serotype %in% Serotype_23F) ~ "23F",
  (cohort$Serotype %in% Other_PCV13) ~ "Other PCV-13",
  (cohort$Serotype %in% Serotype_11A) ~ "11A",
  (cohort$Serotype %in% Serotype_15B) ~ "15B",
  (cohort$Serotype %in% Additional_PCV20) ~ "Additional PCV-20",
  (cohort$Serotype %in% Serotype_7C) ~ "7C",
  (cohort$Serotype %in% Serotype_9N) ~ "9N",
  (cohort$Serotype %in% Serotype_15A) ~ "15A",
  (cohort$Serotype %in% Serotype_15C) ~ "15C",
  (cohort$Serotype %in% Serotype_16F) ~ "16F",
  (cohort$Serotype %in% Serotype_17F) ~ "17F",
  (cohort$Serotype %in% Serotype_21) ~ "21",
  (cohort$Serotype %in% Serotype_23A) ~ "23A",
  (cohort$Serotype %in% Serotype_23B) ~ "23B",
  (cohort$Serotype %in% Serotype_35B) ~ "35B",
  (cohort$Serotype %in% Serotype_35F) ~ "35F",
  .default = "Non-vaccine"))
view(cohort)

table(cohort$spneumo_type)
table(cohort$sero_class)
AC_cat_num <- table(cohort$spneumo_type, cohort$AC_cat)
AC_cat_num
AC_cat_prop <- prop.table(table(cohort$spneumo_type, cohort$AC_cat), 1)
AC_cat_prop
AC_cat_num <- table(cohort$sero_class, cohort$AC_cat)
AC_cat_num
AC_cat_prop <- prop.table(table(cohort$sero_class, cohort$AC_cat), 1)
AC_cat_prop

PG_cat_nonmen_num <- table(cohort$spneumo_type, cohort$PG_cat_nonmen)
PG_cat_nonmen_num
PG_cat_nonmen_prop <- prop.table(table(cohort$spneumo_type, cohort$PG_cat_nonmen), 1)
PG_cat_nonmen_prop
PG_cat_nonmen_num <- table(cohort$sero_class, cohort$PG_cat_nonmen)
PG_cat_nonmen_num
colSums(PG_cat_nonmen_num)
PG_cat_nonmen_prop <- prop.table(table(cohort$sero_class, cohort$PG_cat_nonmen), 1)
PG_cat_nonmen_prop

PG_cat_men_num <- table(cohort$spneumo_type, cohort$PG_cat_men)
PG_cat_men_num
PG_cat_men_prop <- prop.table(table(cohort$spneumo_type, cohort$PG_cat_men), 1)
PG_cat_men_prop
PG_cat_men_num <- table(cohort$sero_class, cohort$PG_cat_men)
PG_cat_men_num
PG_cat_men_prop <- prop.table(table(cohort$sero_class, cohort$PG_cat_men), 1)
PG_cat_men_prop

CRO_cat_num <- table(cohort$spneumo_type, cohort$CRO_cat)
CRO_cat_num
CRO_cat_prop <- prop.table(table(cohort$spneumo_type, cohort$CRO_cat), 1)
CRO_cat_prop
CRO_cat_num <- table(cohort$sero_class, cohort$CRO_cat)
CRO_cat_num
CRO_cat_prop <- prop.table(table(cohort$sero_class, cohort$CRO_cat), 1)
CRO_cat_prop

TS_cat_num <- table(cohort$spneumo_type, cohort$TS_cat)
TS_cat_num
TS_cat_prop <- prop.table(table(cohort$spneumo_type, cohort$TS_cat), 1)
TS_cat_prop
TS_cat_num <- table(cohort$sero_class, cohort$TS_cat)
TS_cat_num
TS_cat_prop <- prop.table(table(cohort$sero_class, cohort$TS_cat), 1)
TS_cat_prop

AZ_cat_num <- table(cohort$spneumo_type, cohort$AZ_cat)
AZ_cat_num
AZ_cat_prop <- prop.table(table(cohort$spneumo_type, cohort$AZ_cat), 1)
AZ_cat_prop
AZ_cat_num <- table(cohort$sero_class, cohort$AZ_cat)
AZ_cat_num
AZ_cat_prop <- prop.table(table(cohort$sero_class, cohort$AZ_cat), 1)
AZ_cat_prop

#Identify class of antibiotic resistance for each of the isolates
PG_nonmen_resist <- filter(cohort, PG_cat_nonmen == "R")
PG_nonmen_resist <- PG_nonmen_resist %>%
  mutate(Abx_resist = "Penicillin, non-meningitis")
table(PG_nonmen_resist$Serotype)
table(PG_nonmen_resist$spneumo_type)
nrow(PG_nonmen_resist)
summary(PG_nonmen_resist$age_days)
PG_nonmen_resist %>%
  group_by(Serotype) %>%
  summarise(Percentage= 100 * n()/nrow(.))
PG_nonmen_resist %>%
  group_by(spneumo_type) %>%
  summarise(Percentage= 100 * n()/nrow(.))

PG_men_resist <- filter(cohort, PG_cat_men == "R")
PG_men_resist <- PG_men_resist %>%
  mutate(Abx_resist = "Penicillin, meningitis")
table(PG_men_resist$Serotype)
table(PG_men_resist$spneumo_type)
nrow(PG_men_resist)
summary(PG_men_resist$age_days)
PG_men_prop <- PG_men_resist %>%
  group_by(Serotype) %>%
  summarise(Percentage= 100 * n()/nrow(.))
PG_men_resist %>%
  group_by(spneumo_type) %>%
  summarise(Percentage= 100 * n()/nrow(.))

AC_resist <- filter(cohort, AC_cat == "R")
AC_resist <- AC_resist %>%
  mutate(Abx_resist = "Amoxicillin")
table(AC_resist$Serotype)
table(AC_resist$spneumo_type)
nrow(AC_resist)
summary(AC_resist$age_days)
AC_resist %>%
  group_by(Serotype) %>%
  summarise(Percentage= 100 * n()/nrow(.))
AC_resist %>%
  group_by(spneumo_type) %>%
  summarise(Percentage= 100 * n()/nrow(.))

AZ_resist <- filter(cohort, AZ_cat == "R")
AZ_resist <- AZ_resist %>%
  mutate(Abx_resist = "Azithromycin")
table(AZ_resist$Serotype)
table(AZ_resist$spneumo_type)
nrow(AZ_resist)
summary(AZ_resist$age_days)
AZ_resist %>%
  group_by(Serotype) %>%
  summarise(Percentage= 100 * n()/nrow(.))
AZ_resist %>%
  group_by(spneumo_type) %>%
  summarise(Percentage= 100 * n()/nrow(.))

CRO_resist <- filter(cohort, CRO_cat == "R")
CRO_resist <- CRO_resist %>%
  mutate(Abx_resist = "Ceftriaxone")
table(CRO_resist$Serotype)
table(CRO_resist$spneumo_type)
nrow(CRO_resist)
summary(CRO_resist$age_days)
CRO_resist %>%
  group_by(Serotype) %>%
  summarise(Percentage= 100 * n()/nrow(.))
CRO_resist %>%
  group_by(spneumo_type) %>%
  summarise(Percentage= 100 * n()/nrow(.))

TS_resist <- filter(cohort, TS_cat == "R")
TS_resist <- TS_resist %>%
  mutate(Abx_resist = "Trimethoprim Sulfamethoxazole")
table(TS_resist$Serotype)
table(TS_resist$spneumo_type)
nrow(TS_resist)
summary(TS_resist$age_days)
TS_resist %>% 
  group_by(Serotype) %>% 
  summarise(Percentage= 100 * n()/nrow(.))
TS_resist %>% 
  group_by(spneumo_type) %>% 
  summarise(Percentage= 100 * n()/nrow(.))

###Identify isolates that are resistant to more than one antibiotic###
Abx_resistance <- rbind(PG_nonmen_resist, AC_resist, AZ_resist, CRO_resist, TS_resist)
summary(Abx_resistance$age_days)
table(Abx_resistance$spneumo_type)
table(Abx_resistance$Serotype)
nrow(Abx_resistance)
#Remove the sample IDs to a separate section and de-duplicate
length(unique(Abx_resistance$SampleID)) #There are 163 samples that are resistant to antibiotics
#Proportion of isolates resistant to any antibiotic
163/264
abx_res_prop <- Abx_resistance %>% 
  group_by(Serotype) %>% 
  summarise(Percentage= 100 * n()/nrow(.))
view(abx_res_prop)
Abx_resistance %>% 
  group_by(spneumo_type) %>% 
  summarise(Percentage= 100 * n()/nrow(.))
#There are 163 isolates that are resistant to any antibiotic, with some being resistant to more than one antibiotic - which ones co-occur?
#Create a new variable to lump together the types of antibiotic resistance by class
#Groupings: Penicillins (PG_cat_nonmen, AC); Cephalosporins (CRO), Sulfonamides (TS), Macrolides (AZ)#
Penicillins <-list("Penicillin, non-meningitis", "Amoxicillin")
Cephalosporins <- list("Ceftriaxone")
Sulfonamides <- list("Trimethoprim Sulfamethoxazole")
Macrolides <- list("Azithromycin")
Abx_resistance <- Abx_resistance %>% mutate(resist_type = case_when(
  (Abx_resistance$Abx_resist %in% Penicillins) ~ "Penicillins",
  (Abx_resistance$Abx_resist %in% Cephalosporins) ~ "Cephalosporins",
  (Abx_resistance$Abx_resist %in% Sulfonamides) ~ "Sulfonamides",
  (Abx_resistance$Abx_resist %in% Macrolides) ~ "Macrolides",
  .default = "NA"))
view(Abx_resistance)

#Drop the NA rows
Abx_resistance = subset(Abx_resistance, select = -c(Abx_resist))
Abx_resistance <- Abx_resistance[Abx_resistance$resist_type != "NA", ]
Abx_resistance <- Abx_resistance %>% distinct()
view(Abx_resistance) 

#Create a resist_any variable and slim table that can be joined to the main cohort dataframe
Abx_resistance$resist_any <- Abx_resistance
  
#Select sample IDs that appear at least 3 times
MDR <-Abx_resistance %>% group_by(SampleID) %>% filter(n()>2) 
view(MDR)
table(MDR$SampleID)
table(MDR$SampleID, MDR$Serotype)

#Create a table summarizing S, R, and I for each antibiotic by year
#Amoxicillin susceptibility by year
AC_cat_overall <- table(cohort$AC_cat)
AC_cat_overall
prop.table(AC_cat_overall)
AC_cat_tab <- table(cohort$AC_cat, cohort$year)
AC_cat_tab
rowSums(AC_cat_tab)
AC_cat_prop <- prop.table((AC_cat_tab), 2)
AC_cat_prop
rowSums(AC_cat_prop)
AC_cat_prop <- as.data.frame(AC_cat_prop)
AC_cat_prop <- AC_cat_prop %>%
  mutate(Antibiotic = "Amoxicillin")
colnames(AC_cat_prop)
names(AC_cat_prop)[names(AC_cat_prop) == "Var1"] <- "Susceptibility"
names(AC_cat_prop)[names(AC_cat_prop) == "Var2"] <- "Year"
names(AC_cat_prop)[names(AC_cat_prop) == "Freq"] <- "Proportion"
#Amoxicillin stacked barchart
Fig3_colors <- c("#BE3428FF", "#DABD61", "#3A488A")
xlab <- c("2016", "2017", "2018", "2019")
Fig3_AC <- AC_cat_prop %>%
  mutate(Susceptibility = fct_relevel(Susceptibility, "R", "I", "S")) %>%
  ggplot( aes(fill=Susceptibility, x=Year, y=Proportion)) + geom_bar(stat= "identity", position = "fill", width = 0.9) +
  scale_fill_manual(values = Fig3_colors, labels=c("Resistant", "Intermediate", "Susceptible"))+
  scale_y_continuous(breaks=seq(0,1,by=0.2))+
  scale_x_discrete(labels=xlab) +
  labs(title="Amoxicillin", x="", y="Proportion of isolates", fill = "Susceptibility")+
  theme_classic() +
  theme(
    plot.title = element_text(size = 24, hjust = 0.5, face="bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 24),
    legend.position = "right",
    axis.text = element_text(size = 24),
    axis.title = element_text(size = 24))+
  theme(axis.text.x=element_blank())
Fig3_AC

#Penicillin meningitis susceptibility by year
PG_cat_men_overall <- table(cohort$PG_cat_men)
PG_cat_men_overall
prop.table(PG_cat_men_overall)
PG_cat_men_tab <- table(cohort$PG_cat_men, cohort$year)
PG_cat_men_tab
rowSums(PG_cat_men_tab)
PG_cat_men_prop <- prop.table((PG_cat_men_tab), 2)
PG_cat_men_prop
rowSums(PG_cat_men_prop)
PG_cat_men_prop <- as.data.frame(PG_cat_men_prop)
PG_cat_men_prop <- PG_cat_men_prop %>%
  mutate(Antibiotic = "Penicilin-meningitis")
colnames(PG_cat_men_prop)
names(PG_cat_men_prop)[names(PG_cat_men_prop) == "Var1"] <- "Susceptibility"
names(PG_cat_men_prop)[names(PG_cat_men_prop) == "Var2"] <- "Year"
names(PG_cat_men_prop)[names(PG_cat_men_prop) == "Freq"] <- "Proportion"
view(PG_cat_men_prop)
#Penicillin-meningitis stacked barchart
Fig3B_colors <- c("#BE3428FF", "#3A488A")
xlab <- c("2016", "2017", "2018", "2019")
Fig3_PG_men <- PG_cat_men_prop %>%
  mutate(Susceptibility = fct_relevel(Susceptibility, "R", "S")) %>%
  ggplot( aes(fill=Susceptibility, x=Year, y=Proportion)) + geom_bar(stat= "identity", position = "fill", width = 0.9) +
  scale_fill_manual(values = Fig3B_colors, labels=c("Resistant", "Susceptible")) +
  scale_y_continuous(breaks=seq(0,1,by=0.2))+  
  scale_x_discrete(labels=xlab) +
  labs(title="Penicillin (meningitis)", x="", y="Proportion of isolates", fill = "Susceptibility")+
  theme_classic() +
  theme(
    plot.title = element_text(size = 24, hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 24),
    legend.text = element_text(size = 24),
    legend.position = "right",
    axis.text = element_text(size = 24),
    axis.title = element_text(size = 24))
Fig3_PG_men

#Penicillin non-meningitis susceptibility by year
PG_cat_nonmen_overall <- table(cohort$PG_cat_nonmen)
PG_cat_nonmen_overall
prop.table(PG_cat_nonmen_overall)
PG_cat_nonmen_tab <- table(cohort$PG_cat_nonmen, cohort$year)
PG_cat_nonmen_tab
rowSums(PG_cat_nonmen_tab)
PG_cat_nonmen_prop <- prop.table((PG_cat_nonmen_tab), 2)
PG_cat_nonmen_prop
rowSums(PG_cat_nonmen_prop)
PG_cat_nonmen_prop <- as.data.frame(PG_cat_nonmen_prop)
PG_cat_nonmen_prop <- PG_cat_nonmen_prop %>%
  mutate(Antibiotic = "Penicilin - nonmeningitis")
colnames(PG_cat_nonmen_prop)
names(PG_cat_nonmen_prop)[names(PG_cat_nonmen_prop) == "Var1"] <- "Susceptibility"
names(PG_cat_nonmen_prop)[names(PG_cat_nonmen_prop) == "Var2"] <- "Year"
names(PG_cat_nonmen_prop)[names(PG_cat_nonmen_prop) == "Freq"] <- "Proportion"
#PG_cat_nonmen stacked barchart
Fig3C_colors <- c("#DABD61", "#3A488A")
xlab <- c("2016", "2017", "2018", "2019")
Fig3_PG_nonmen <- PG_cat_nonmen_prop %>%
  mutate(Susceptibility = fct_relevel(Susceptibility, "I", "S")) %>%
  ggplot( aes(fill=Susceptibility, x=Year, y=Proportion)) + geom_bar(stat= "identity", position = "fill", width = 0.9) +
  scale_fill_manual(values = Fig3C_colors, labels=c("Intermediate", "Susceptible")) +
  scale_y_continuous(breaks=seq(0,1,by=0.2))+   
  scale_x_discrete(labels=xlab) +
  labs(title="Penicillin (non-meningitis)", x="", y="", fill = "Susceptibility")+
  theme_classic() +
  theme(
    plot.title = element_text(size = 24, hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 24),
    legend.text = element_text(size = 24),
    legend.position = "right")+
  theme(axis.text.x = element_text(size = 24))+
  theme(axis.text.y=element_blank())
Fig3_PG_nonmen

#Ceftriaxone susceptibility by year
CRO_cat_overall <- table(cohort$CRO_cat)
CRO_cat_overall
prop.table(CRO_cat_overall)
CRO_cat_tab <- table(cohort$CRO_cat, cohort$year)
CRO_cat_tab
rowSums(CRO_cat_tab)
CRO_cat_prop <- prop.table((CRO_cat_tab), 2)
CRO_cat_prop
rowSums(CRO_cat_prop)
CRO_cat_prop <- as.data.frame(CRO_cat_prop)
CRO_cat_prop <- CRO_cat_prop %>%
  mutate(Antibiotic = "Ceftriaxone")
colnames(CRO_cat_prop)
names(CRO_cat_prop)[names(CRO_cat_prop) == "Var1"] <- "Susceptibility"
names(CRO_cat_prop)[names(CRO_cat_prop) == "Var2"] <- "Year"
names(CRO_cat_prop)[names(CRO_cat_prop) == "Freq"] <- "Proportion"
view(CRO_cat_prop)
#CRO_cat stacked barchart
Fig3_colors <- c("#BE3428FF", "#DABD61", "#3A488A")
xlab <- c("2016", "2017", "2018", "2019")
Fig3_CRO <- CRO_cat_prop %>%
  mutate(Susceptibility = fct_relevel(Susceptibility, "R", "I", "S")) %>%
  ggplot( aes(fill=Susceptibility, x=Year, y=Proportion)) + geom_bar(stat= "identity", position = "fill", width = 0.9) +
  scale_fill_manual(values = Fig3_colors, labels=c("Resistant", "Intermediate", "Susceptible")) +
  scale_y_continuous(breaks=seq(0,1,by=0.2))+  
  scale_x_discrete(labels=xlab) +
  labs(title="Ceftriaxone", x="", y="", fill = "Susceptibility")+
  theme_classic() +
  theme(
    plot.title = element_text(size = 24, hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 24),
    legend.text = element_text(size = 24),
    legend.position = "right")+
  theme(axis.text.x=element_blank(), axis.text.y=element_blank())
Fig3_CRO

#Trimethoprim sulfamethoxazole susceptibility by year
TS_cat_overall <- table(cohort$TS_cat)
TS_cat_overall
prop.table(TS_cat_overall)
TS_cat_tab <- table(cohort$TS_cat, cohort$year)
TS_cat_tab
rowSums(TS_cat_tab)
TS_cat_prop <- prop.table((TS_cat_tab), 2)
TS_cat_prop
rowSums(TS_cat_prop)
TS_cat_prop <- as.data.frame(TS_cat_prop)
TS_cat_prop <- TS_cat_prop %>%
  mutate(Antibiotic = "Trimethoprim-sulfamethoxazole")
colnames(TS_cat_prop)
names(TS_cat_prop)[names(TS_cat_prop) == "Var1"] <- "Susceptibility"
names(TS_cat_prop)[names(TS_cat_prop) == "Var2"] <- "Year"
names(TS_cat_prop)[names(TS_cat_prop) == "Freq"] <- "Proportion"
#TS_cat stacked barchart
Fig3_colors <- c("#BE3428FF", "#DABD61", "#3A488A")
xlab <- c("2016", "2017", "2018", "2019")
Fig3_TS <- TS_cat_prop %>%
  mutate(Susceptibility = fct_relevel(Susceptibility, "R", "I", "S")) %>%
  ggplot( aes(fill=Susceptibility, x=Year, y=Proportion)) + geom_bar(stat= "identity", position = "fill", width = 0.9) +
  scale_fill_manual(values = Fig3_colors, labels=c("Resistant", "Intermediate", "Susceptible")) +
  scale_y_continuous(breaks=seq(0,1,by=0.2))+  
  scale_x_discrete(labels=xlab) +
  labs(title="Trimethoprim-sulfamethoxazole", x="", y="")+
  theme_classic() +
  theme(
    plot.title = element_text(size = 24, hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 24),
    legend.text = element_text(size = 24),
    legend.position = "right")+
  theme(axis.text.x = element_text(size = 24))+
  theme(axis.text.y=element_blank())
Fig3_TS

#Azithromycin susceptibility by year
AZ_cat_overall <- table(cohort$AZ_cat)
AZ_cat_overall
prop.table(AZ_cat_overall)
AZ_cat_tab <- table(cohort$AZ_cat, cohort$year)
AZ_cat_tab
rowSums(AZ_cat_tab)
AZ_cat_prop <- prop.table((AZ_cat_tab), 2)
AZ_cat_prop
rowSums(AZ_cat_prop)
AZ_cat_prop <- as.data.frame(AZ_cat_prop)
AZ_cat_prop <- AZ_cat_prop %>%
  mutate(Antibiotic = "Azithromycin")
colnames(AZ_cat_prop)
names(AZ_cat_prop)[names(AZ_cat_prop) == "Var1"] <- "Susceptibility"
names(AZ_cat_prop)[names(AZ_cat_prop) == "Var2"] <- "Year"
names(AZ_cat_prop)[names(AZ_cat_prop) == "Freq"] <- "Proportion"
#AZ_cat stacked barchart
Fig3_colors <- c("#BE3428FF", "#DABD61", "#3A488A")
xlab <- c("2016", "2017", "2018", "2019")
Fig3_AZ <- AZ_cat_prop %>%
  mutate(Susceptibility = fct_relevel(Susceptibility, "R", "I", "S")) %>%
  ggplot( aes(fill=Susceptibility, x=Year, y=Proportion)) + geom_bar(stat= "identity", position = "fill", width = 0.9) +
  scale_fill_manual(values = Fig3_colors, labels=c("Resistant", "Intermediate", "Susceptible")) +
  scale_y_continuous(breaks=seq(0,1,by=0.2))+  
  scale_x_discrete(labels=xlab) +
  labs(title="Azithromycin", x="", y="", fill = "Susceptibility")+
  theme_classic() +
  theme(
    plot.title = element_text(size = 24, hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 24),
    legend.text = element_text(size = 24),
    legend.position = "right")+
  theme(axis.text.x=element_blank(), axis.text.y=element_blank())
Fig3_AZ

##FIGURE 3##
#Combine the stacked barcharts above
library("cowplot")
# extract the legend from one of the plots
legend_3 <- get_legend(Fig3_AC + guides(color = guide_legend(nrow = 1)) +
  theme(legend.position = "bottom"))

library(patchwork)
layout <-"
ABC
DEF
"
Fig3_combined <- 
  (Fig3_AC + theme(legend.position="none"))+
  (Fig3_AZ + theme(legend.position="none"))+
  (Fig3_CRO + theme(legend.position="none"))+ 
  (Fig3_PG_men + theme(legend.position="none"))+ 
  (Fig3_PG_nonmen + theme(legend.position="none"))+ 
  (Fig3_TS + theme(legend.position="none"))+
  plot_layout(design = layout)

Fig3_combined

Fig3_combined_legend <- plot_grid(Fig3_combined, legend_3, ncol = 1, rel_heights = c(1, 0.05)) 
Fig3_combined_legend
ggsave("Fig3_combined.png", plot = Fig3_combined_legend, width = 22, height = 16, dpi = 800)

#Create dummy variables for resistant and not susceptible to evaluate for each of the antibiotics and all antibiotics
#We are not counting the Penicillin meningitis cut-off
cohort<- cohort %>%
  mutate(resist_any = case_when(
   PG_cat_nonmen == "R" | AC_cat == "R" | AZ_cat == "R" |CRO_cat == "R" | TS_cat == "R" ~ 1, 
    TRUE ~ 0))
view(cohort)

cohort<- cohort %>%
  mutate(not_susceptible_any = case_when(
    PG_cat_nonmen == "I" | AC_cat == "I" | AZ_cat == "I" |CRO_cat == "I" | TS_cat == "I" ~ 1, PG_cat_nonmen == "R" | AC_cat == "R" | AZ_cat == "R" |CRO_cat == "R" | TS_cat == "R" ~ 1,
    TRUE ~ 0))
view(cohort)


#Summarize antibiotic susceptibility, stratified by year
aim3_NS_year_tab <- do.call(rbind , by(cohort$not_susceptible_any, cohort$year, table)) 
aim3_NS_year_tab
rowSums(aim3_NS_year_tab)
colSums(aim3_NS_year_tab)
67/264
197/264
aim3_NS_year_prop <- prop.table(aim3_NS_year_tab,1)
aim3_NS_year_prop
rowSums(aim3_NS_year_prop)


#Summarize antibiotic resistance, stratified by year
aim3_resist_year_tab <- do.call(rbind , by(cohort$resist_any, cohort$year, table)) 
aim3_resist_year_tab
rowSums(aim3_resist_year_tab)
colSums(aim3_resist_year_tab)
59/264
205/264
aim3_resist_year_prop <- prop.table(aim3_resist_year_tab,1)
aim3_resist_year_prop

#Run a Cochran-Armitage trend test to determine if there is a significant increase in non-susceptible strains during the study period
#Create a contingency table
year <- cohort$year
not_susceptible <- cohort$not_susceptible_any
aim3_NS_year <- table(not_susceptible, year)
aim3_NS_year
sum(aim3_NS_year)
prop.table(aim3_NS_year,
           margin = 2)
Test_aim3_NS = chisq_test(aim3_NS_year,
                              scores = list("year" = c(2016, 2017, 2018, 2019)))
Test_aim3_NS


#Run a Cochran-Armitage trend test to determine if there is a significant increase in resistant strains during the study period
#Create a contingency table
class(cohort$year)
year <- cohort$year
resist_any <- cohort$resist_any
aim3_resist_year <- table(resist_any, year)
aim3_resist_year
sum(aim3_resist_year)
  prop.table(aim3_resist_year,
           margin = 2)
library(coin)
Test_aim3_resist = chisq_test(aim3_resist_year,
                  scores = list("year" = c(2016, 2017, 2018, 2019)))
Test_aim3_resist


#Run a Cochrane Armitage test to see if there are changes in the proportions of PCV-13 vs. nonvaccine strains that are not susceptible
PCV13_isolates <- cohort[cohort$colonized_PCV13 == 1,]
year <- PCV13_isolates$year
resist_any <- PCV13_isolates$resist_any
PCV13_resist_year <- table(resist_any, year)
PCV13_resist_year
sum(PCV13_resist_year)
prop.table(PCV13_resist_year,
           margin = 2)
library(coin)
Test_PCV13_resist = chisq_test(PCV13_resist_year,
                              scores = list("year" = c(2016, 2017, 2018, 2019)))
Test_PCV13_resist

Non_vax_isolates <- cohort[cohort$colonized_PCV13 == 0,]
year <- Non_vax_isolates$year
resist_any <- Non_vax_isolates$resist_any
Nonvax_resist_year <- table(resist_any, year)
Nonvax_resist_year
sum(Nonvax_resist_year)
prop.table(Nonvax_resist_year,
           margin = 2)
library(coin)
Test_nonvax_resist = chisq_test(Nonvax_resist_year,
                               scores = list("year" = c(2016, 2017, 2018, 2019)))
Test_nonvax_resist


#Determine whether the proportion of PCV13 or nonPCV13 isolates that are not susceptible to an antibiotic changes over time
year <- PCV13_isolates$year
not_susceptible <- PCV13_isolates$not_susceptible_any
PCV13_NS_year <- table(not_susceptible, year)
PCV13_NS_year
sum(PCV13_NS_year)
prop.table(PCV13_NS_year,
           margin = 2)
Test_PCV13_NS = chisq_test(PCV13_NS_year,
                          scores = list("year" = c(2016, 2017, 2018, 2019)))
Test_PCV13_NS

Non_vax_isolates <- cohort[cohort$colonized_PCV13 == 0,]
year <- Non_vax_isolates$year
not_susceptible_any <- Non_vax_isolates$not_susceptible_any
Nonvax_NS_year <- table(not_susceptible_any, year)
Nonvax_NS_year
sum(Nonvax_NS_year)
prop.table(Nonvax_NS_year,
           margin = 2)
library(coin)
Test_nonvax_NS = chisq_test(Nonvax_NS_year,
                                scores = list("year" = c(2016, 2017, 2018, 2019)))
Test_nonvax_NS

#Create dummy variables for amoxicillin non-susceptibility or resistance
res_int <- list("R", "I")
cohort <- cohort %>% mutate(AC_not_susceptible = case_when(
  (cohort$AC_cat %in% res_int) ~ 1,
  .default = 0))
resistant <- list("R")
cohort <- cohort %>% mutate(AC_resistant = case_when(
  (cohort$AC_cat %in% resistant) ~ 1,
  .default = 0))
view(cohort)
#Summarize amoxicillin non-susceptibility, stratified by year
AC_NS_year_tab <- do.call(rbind , by(cohort$AC_not_susceptible, cohort$year, table)) 
AC_NS_year_tab
AC_NS_year_prop <- prop.table(AC_NS_year_tab,1)
AC_NS_year_prop
#Summarize amoxicillin resistance, stratified by year
AC_resist_year_tab <- do.call(rbind , by(cohort$AC_resistant, cohort$year, table)) 
AC_resist_year_tab
AC_resist_year_prop <- prop.table(AC_resist_year_tab,1)
AC_resist_year_prop
#Run a Cochran-Armitage trend test to determine if there is a significant increase in Amoxicillin non-susceptibility
class(cohort$year)
year <- cohort$year
AC_NS <- cohort$AC_not_susceptible
AC_NS_year <- table(AC_NS, year)
AC_NS_year
sum(AC_NS_year)
prop.table(AC_NS_year,
           margin = 2)
Test_AC_NS = chisq_test(AC_NS_year,
                              scores = list("year" = c(2016, 2017, 2018, 2019)))
Test_AC_NS
#Let's repeat this and look at PCV13 and non-PCV13 isolates separately
PCV13_isolates <- cohort[cohort$colonized_PCV13 == 1,]
year <- PCV13_isolates$year
AC_NS <- PCV13_isolates$AC_not_susceptible
AC_NS_year <- table(AC_NS, year)
AC_NS_year
sum(AC_NS_year)
prop.table(AC_NS_year,
           margin = 2)
Test_AC_NS = chisq_test(AC_NS_year,
                        scores = list("year" = c(2016, 2017, 2018, 2019)))
Test_AC_NS

Nonvax_isolates <- cohort[cohort$colonized_PCV13 == 0,]
year <- Nonvax_isolates$year
AC_NS <- Nonvax_isolates$AC_not_susceptible
AC_NS_year <- table(AC_NS, year)
AC_NS_year
sum(AC_NS_year)
prop.table(AC_NS_year,
           margin = 2)
Test_AC_NS = chisq_test(AC_NS_year,
                        scores = list("year" = c(2016, 2017, 2018, 2019)))
Test_AC_NS

#Run a Cochran-Armitage trend test to determine if there is a significant increase in Amoxicillin resistance
class(cohort$year)
year <- cohort$year
AC_resist <- cohort$AC_resistant
AC_resist_year <- table(AC_resist, year)
AC_resist_year
sum(AC_resist_year)
prop.table(AC_resist_year,
           margin = 2)
Test_AC_resist = chisq_test(AC_resist_year,
                        scores = list("year" = c(2016, 2017, 2018, 2019)))
Test_AC_resist
#Look at PCV13 and non-PCV13 separately
PCV13_isolates <- cohort[cohort$colonized_PCV13 == 1,]
year <- PCV13_isolates$year
AC_resist <- PCV13_isolates$AC_resistant
AC_resist_year <- table(AC_resist, year)
AC_resist_year
sum(AC_resist_year)
prop.table(AC_resist_year,
           margin = 2)
Test_AC_resist = chisq_test(AC_resist_year,
                            scores = list("year" = c(2016, 2017, 2018, 2019)))
Test_AC_resist

Non_vax_isolates <- cohort[cohort$colonized_PCV13 == 1,]
year <- Non_vax_isolates$year
AC_resist <- Non_vax_isolates$AC_resistant
AC_resist_year <- table(AC_resist, year)
AC_resist_year
sum(AC_resist_year)
prop.table(AC_resist_year,
           margin = 2)
Test_AC_resist = chisq_test(AC_resist_year,
                            scores = list("year" = c(2016, 2017, 2018, 2019)))
Test_AC_resist



#Create dummy variables for Penicillin meningitis non-susceptibility or resistance
res_int <- list("R", "I")
cohort <- cohort %>% mutate(PG_men_not_susceptible = case_when(
  (cohort$PG_cat_men %in% res_int) ~ 1,
  .default = 0))
resistant <- list("R")
cohort <- cohort %>% mutate(PG_men_resistant = case_when(
  (cohort$PG_cat_men %in% resistant) ~ 1,
  .default = 0))
#Summarize Penicillin meningitis non-susceptibility, stratified by year
PG_men_NS_year_tab <- do.call(rbind , by(cohort$PG_men_not_susceptible, cohort$year, table)) 
PG_men_NS_year_tab
PG_men_NS_year_prop <- prop.table(PG_men_NS_year_tab,1)
PG_men_NS_year_prop
#Summarize penicillin meningitis resistance, stratified by year
PG_men_resist_year_tab <- do.call(rbind , by(cohort$PG_men_resistant, cohort$year, table)) 
PG_men_resist_year_tab
PG_men_resist_year_prop <- prop.table(PG_men_resist_year_tab,1)
PG_men_resist_year_prop
#Run a Cochran-Armitage trend test to determine if there is a significant increase in Penicillin meningitis non-susceptibility
class(cohort$year)
year <- cohort$year
PG_men_NS <- cohort$PG_men_not_susceptible
PG_men_NS_year <- table(PG_men_NS, year)
PG_men_NS_year
sum(PG_men_NS_year)
prop.table(PG_men_NS_year,
           margin = 2)
Test_PG_men_NS = chisq_test(PG_men_NS_year,
                        scores = list("year" = c(2016, 2017, 2018, 2019)))
Test_PG_men_NS
#Look at PCV13 and non-PCV13 separately
PCV13_isolates <- cohort[cohort$colonized_PCV13 == 1,]
year <- PCV13_isolates$year
PG_men_NS <- PCV13_isolates$PG_men_not_susceptible
PG_men_NS_year <- table(PG_men_NS, year)
PG_men_NS_year
sum(PG_men_NS_year)
prop.table(PG_men_NS_year,
           margin = 2)
Test_PG_men_NS = chisq_test(PG_men_NS_year,
                            scores = list("year" = c(2016, 2017, 2018, 2019)))
Test_PG_men_NS

Non_vax_isolates <- cohort[cohort$colonized_PCV13 == 0,]
year <- Non_vax_isolates$year
PG_men_NS <- Non_vax_isolates$PG_men_not_susceptible
PG_men_NS_year <- table(PG_men_NS, year)
PG_men_NS_year
sum(PG_men_NS_year)
prop.table(PG_men_NS_year,
           margin = 2)
Test_PG_men_NS = chisq_test(PG_men_NS_year,
                            scores = list("year" = c(2016, 2017, 2018, 2019)))
Test_PG_men_NS

#Run a Cochran-Armitage trend test to determine if there is a significant increase in Penicillin meningitis resistance
class(cohort$year)
year <- cohort$year
PG_men_resist <- cohort$PG_men_resistant
PG_men_resist_year <- table(PG_men_resist, year)
PG_men_resist_year
sum(PG_men_resist_year)
prop.table(PG_men_resist_year,
           margin = 2)
Test_PG_men_resist = chisq_test(PG_men_resist_year,
                            scores = list("year" = c(2016, 2017, 2018, 2019)))
Test_PG_men_resist
#Look at PCV13 and non-PCV13 isolates separately
PCV13_isolates <- cohort[cohort$colonized_PCV13 == 1,]
year <- PCV13_isolates$year
PG_men_resist <- PCV13_isolates$PG_men_resistant
PG_men_resist_year <- table(PG_men_resist, year)
PG_men_resist_year
sum(PG_men_resist_year)
prop.table(PG_men_resist_year,
           margin = 2)
Test_PG_men_resist = chisq_test(PG_men_resist_year,
                                scores = list("year" = c(2016, 2017, 2018, 2019)))
Test_PG_men_resist

Nonvax_isolates <- cohort[cohort$colonized_PCV13 == 0,]
year <- Nonvax_isolates$year
PG_men_resist <- Nonvax_isolates$PG_men_resistant
PG_men_resist_year <- table(PG_men_resist, year)
PG_men_resist_year
sum(PG_men_resist_year)
prop.table(PG_men_resist_year,
           margin = 2)
Test_PG_men_resist = chisq_test(PG_men_resist_year,
                                scores = list("year" = c(2016, 2017, 2018, 2019)))
Test_PG_men_resist


#Create dummy variables for Penicillin NON-meningitis non-susceptibility or resistance
res_int <- list("R", "I")
cohort <- cohort %>% mutate(PG_nonmen_not_susceptible = case_when(
  (cohort$PG_cat_nonmen %in% res_int) ~ 1,
  .default = 0))
resistant <- list("R")
cohort <- cohort %>% mutate(PG_nonmen_resistant = case_when(
  (cohort$PG_cat_nonmen %in% resistant) ~ 1,
  .default = 0))
view(cohort)
#Summarize Penicillin NON-meningitis non-susceptibility, stratified by year
PG_nonmen_NS_year_tab <- do.call(rbind , by(cohort$PG_nonmen_not_susceptible, cohort$year, table)) 
PG_nonmen_NS_year_tab
PG_nonmen_NS_year_prop <- prop.table(PG_nonmen_NS_year_tab,1)
PG_nonmen_NS_year_prop
#Summarize penicillin NON-meningitis resistance, stratified by year
PG_nonmen_resist_year_tab <- do.call(rbind , by(cohort$PG_nonmen_resistant, cohort$year, table)) 
PG_nonmen_resist_year_tab
PG_nonmen_resist_year_prop <- prop.table(PG_nonmen_resist_year_tab,1)
PG_nonmen_resist_year_prop
#Run a Cochran-Armitage trend test to determine if there is a significant increase in Penicillin NON-meningitis non-susceptibility
class(cohort$year)
year <- cohort$year
PG_nonmen_NS <- cohort$PG_nonmen_not_susceptible
PG_nonmen_NS_year <- table(PG_nonmen_NS, year)
PG_nonmen_NS_year
sum(PG_nonmen_NS_year)
prop.table(PG_nonmen_NS_year,
           margin = 2)
Test_PG_nonmen_NS = chisq_test(PG_nonmen_NS_year,
                            scores = list("year" = c(2016, 2017, 2018, 2019)))
Test_PG_nonmen_NS
#See if there are differences in PCV13 and non-vaccine serotypes
PCV13_isolates <- cohort[cohort$colonized_PCV13 == 1,]
year <- PCV13_isolates$year
PG_nonmen_NS <- PCV13_isolates$PG_nonmen_not_susceptible
PG_nonmen_NS_year <- table(PG_nonmen_NS, year)
PG_nonmen_NS_year
sum(PG_nonmen_NS_year)
prop.table(PG_nonmen_NS_year,
           margin = 2)
Test_PG_nonmen_NS = chisq_test(PG_nonmen_NS_year,
                               scores = list("year" = c(2016, 2017, 2018, 2019)))
Test_PG_nonmen_NS

Non_vax_isolates <- cohort[cohort$colonized_PCV13 == 0,]
year <- Non_vax_isolates$year
PG_nonmen_NS <- Non_vax_isolates$PG_nonmen_not_susceptible
PG_nonmen_NS_year <- table(PG_nonmen_NS, year)
PG_nonmen_NS_year
sum(PG_nonmen_NS_year)
prop.table(PG_nonmen_NS_year,
           margin = 2)
Test_PG_nonmen_NS = chisq_test(PG_nonmen_NS_year,
                               scores = list("year" = c(2016, 2017, 2018, 2019)))
Test_PG_nonmen_NS

#Run a Cochran-Armitage trend test to determine if there is a significant increase in Penicillin NON-meningitis resistance
class(cohort$year)
year <- cohort$year
PG_nonmen_resist <- cohort$PG_nonmen_resistant
PG_nonmen_resist_year <- table(PG_nonmen_resist, year)
PG_nonmen_resist_year
sum(PG_nonmen_resist_year)
prop.table(PG_nonmen_resist_year,
           margin = 2)
Test_PG_nonmen_resist = chisq_test(PG_nonmen_resist_year,
                                scores = list("year" = c(2016, 2017, 2018, 2019)))
Test_PG_nonmen_resist
#Check PCV13 and non-PCV13 isolates
PCV13_isolates <- cohort[cohort$colonized_PCV13 == 1,]
year <- PCV13_isolates$year
PG_nonmen_resist <- PCV13_isolates$PG_nonmen_resistant
PG_nonmen_resist_year <- table(PG_nonmen_resist, year)
PG_nonmen_resist_year
sum(PG_nonmen_resist_year)
prop.table(PG_nonmen_resist_year,
           margin = 2)
Test_PG_nonmen_resist = chisq_test(PG_nonmen_resist_year,
                                   scores = list("year" = c(2016, 2017, 2018, 2019)))
Test_PG_nonmen_resist

Non_vax_isolates <- cohort[cohort$colonized_PCV13 == 0,]
year <- Non_vax_isolates$year
PG_nonmen_resist <- Non_vax_isolates$PG_nonmen_resistant
PG_nonmen_resist_year <- table(PG_nonmen_resist, year)
PG_nonmen_resist_year
sum(PG_nonmen_resist_year)
prop.table(PG_nonmen_resist_year,
           margin = 2)
Test_PG_nonmen_resist = chisq_test(PG_nonmen_resist_year,
                                   scores = list("year" = c(2016, 2017, 2018, 2019)))
Test_PG_nonmen_resist
#There are no penicillin resistant isolates at the non-meningitis cut-off

#Create dummy variables for Ceftriaxone non-susceptibility or resistance
res_int <- list("R", "I")
cohort <- cohort %>% mutate(CRO_not_susceptible = case_when(
  (cohort$CRO_cat %in% res_int) ~ 1,
  .default = 0))
resistant <- list("R")
cohort <- cohort %>% mutate(CRO_resistant = case_when(
  (cohort$CRO_cat %in% resistant) ~ 1,
  .default = 0))
view(cohort)
#Summarize Ceftriaxone non-susceptibility, stratified by year
CRO_NS_year_tab <- do.call(rbind , by(cohort$CRO_not_susceptible, cohort$year, table)) 
CRO_NS_year_tab
CRO_NS_year_prop <- prop.table(CRO_NS_year_tab,1)
CRO_NS_year_prop
#Summarize Ceftriaxone resistance, stratified by year
CRO_resist_year_tab <- do.call(rbind , by(cohort$CRO_resistant, cohort$year, table)) 
CRO_resist_year_tab
CRO_resist_year_prop <- prop.table(CRO_resist_year_tab,1)
CRO_resist_year_prop
#Run a Cochran-Armitage trend test to determine if there is a significant increase in Ceftriaxone non-susceptibility
class(cohort$year)
year <- cohort$year
CRO_NS <- cohort$CRO_not_susceptible
CRO_NS_year <- table(CRO_NS, year)
CRO_NS_year
sum(CRO_NS_year)
prop.table(CRO_NS_year,
           margin = 2)
Test_CRO_NS = chisq_test(CRO_NS_year,
                               scores = list("year" = c(2016, 2017, 2018, 2019)))
Test_CRO_NS
#Check for differences in PCV13 and non-PCV13 isolates
PCV13_isolates <- cohort[cohort$colonized_PCV13 == 1,]
year <- PCV13_isolates$year
CRO_NS <- PCV13_isolates$CRO_not_susceptible
CRO_NS_year <- table(CRO_NS, year)
CRO_NS_year
sum(CRO_NS_year)
prop.table(CRO_NS_year,
           margin = 2)
Test_CRO_NS = chisq_test(CRO_NS_year,
                         scores = list("year" = c(2016, 2017, 2018, 2019)))
Test_CRO_NS

Non_vax_isolates <- cohort[cohort$colonized_PCV13 == 0,]
year <- Non_vax_isolates$year
CRO_NS <- Non_vax_isolates$CRO_not_susceptible
CRO_NS_year <- table(CRO_NS, year)
CRO_NS_year
sum(CRO_NS_year)
prop.table(CRO_NS_year,
           margin = 2)
Test_CRO_NS = chisq_test(CRO_NS_year,
                         scores = list("year" = c(2016, 2017, 2018, 2019)))
Test_CRO_NS


#Run a Cochran-Armitage trend test to determine if there is a significant increase in Ceftriaxone resistance
class(cohort$year)
year <- cohort$year
CRO_resist <- cohort$CRO_resistant
CRO_resist_year <- table(CRO_resist, year)
CRO_resist_year
sum(CRO_resist_year)
prop.table(CRO_resist_year,
           margin = 2)
Test_CRO_resist = chisq_test(CRO_resist_year,
                                   scores = list("year" = c(2016, 2017, 2018, 2019)))
Test_CRO_resist
#Check for differences in PCV13 and non-PCV13 isolates
PCV13_isolates <- cohort[cohort$colonized_PCV13 == 1,]
year <- PCV13_isolates$year
CRO_resist <- PCV13_isolates$CRO_resistant
CRO_resist_year <- table(CRO_resist, year)
CRO_resist_year
sum(CRO_resist_year)
prop.table(CRO_resist_year,
           margin = 2)
Test_CRO_resist = chisq_test(CRO_resist_year,
                             scores = list("year" = c(2016, 2017, 2018, 2019)))
Test_CRO_resist

Non_vax_isolates <- cohort[cohort$colonized_PCV13 == 0,]
year <- Non_vax_isolates$year
CRO_resist <- Non_vax_isolates$CRO_resistant
CRO_resist_year <- table(CRO_resist, year)
CRO_resist_year
sum(CRO_resist_year)
prop.table(CRO_resist_year,
           margin = 2)
Test_CRO_resist = chisq_test(CRO_resist_year,
                             scores = list("year" = c(2016, 2017, 2018, 2019)))
Test_CRO_resist


#Create dummy variables for TS non-susceptibility or resistance
res_int <- list("R", "I")
cohort <- cohort %>% mutate(TS_not_susceptible = case_when(
  (cohort$TS_cat %in% res_int) ~ 1,
  .default = 0))
resistant <- list("R")
cohort <- cohort %>% mutate(TS_resistant = case_when(
  (cohort$TS_cat %in% resistant) ~ 1,
  .default = 0))
view(cohort)
#Summarize TS non-susceptibility, stratified by year
TS_NS_year_tab <- do.call(rbind , by(cohort$TS_not_susceptible, cohort$year, table)) 
TS_NS_year_tab
TS_NS_year_prop <- prop.table(TS_NS_year_tab,1)
TS_NS_year_prop
#Summarize TS resistance, stratified by year
TS_resist_year_tab <- do.call(rbind , by(cohort$TS_resistant, cohort$year, table)) 
TS_resist_year_tab
TS_resist_year_prop <- prop.table(TS_resist_year_tab,1)
TS_resist_year_prop
#Run a Cochran-Armitage trend test to determine if there is a significant increase in TS non-susceptibility
class(cohort$year)
year <- cohort$year
TS_NS <- cohort$TS_not_susceptible
TS_NS_year <- table(TS_NS, year)
TS_NS_year
sum(TS_NS_year)
prop.table(TS_NS_year,
           margin = 2)
Test_TS_NS = chisq_test(TS_NS_year,
                         scores = list("year" = c(2016, 2017, 2018, 2019)))
Test_TS_NS
#Check for differences in PCV13 vs. non-vax isolates
PCV13_isolates <- cohort[cohort$colonized_PCV13 == 1,]
year <- PCV13_isolates$year
TS_NS <- PCV13_isolates$TS_not_susceptible
TS_NS_year <- table(TS_NS, year)
TS_NS_year
sum(TS_NS_year)
prop.table(TS_NS_year,
           margin = 2)
Test_TS_NS = chisq_test(TS_NS_year,
                        scores = list("year" = c(2016, 2017, 2018, 2019)))
Test_TS_NS

Non_vax_isolates <- cohort[cohort$colonized_PCV13 == 0,]
year <- Non_vax_isolates$year
TS_NS <- Non_vax_isolates$TS_not_susceptible
TS_NS_year <- table(TS_NS, year)
TS_NS_year
sum(TS_NS_year)
prop.table(TS_NS_year,
           margin = 2)
Test_TS_NS = chisq_test(TS_NS_year,
                        scores = list("year" = c(2016, 2017, 2018, 2019)))
Test_TS_NS

#Run a Cochran-Armitage trend test to determine if there is a significant increase in TS resistance
class(cohort$year)
year <- cohort$year
TS_resist <- cohort$TS_resistant
TS_resist_year <- table(TS_resist, year)
TS_resist_year
sum(TS_resist_year)
prop.table(TS_resist_year,
           margin = 2)
Test_TS_resist = chisq_test(TS_resist_year,
                             scores = list("year" = c(2016, 2017, 2018, 2019)))
Test_TS_resist
#Check to see if there are differences in PCV13 vs. nonPCV13 isolates
PCV13_isolates <- cohort[cohort$colonized_PCV13 == 1,]
year <- PCV13_isolates$year
TS_resist <- PCV13_isolates$TS_resistant
TS_resist_year <- table(TS_resist, year)
TS_resist_year
sum(TS_resist_year)
prop.table(TS_resist_year,
           margin = 2)
Test_TS_resist = chisq_test(TS_resist_year,
                            scores = list("year" = c(2016, 2017, 2018, 2019)))
Test_TS_resist

Non_vax_isolates <- cohort[cohort$colonized_PCV13 == 0,]
year <- Non_vax_isolates$year
TS_resist <- Non_vax_isolates$TS_resistant
TS_resist_year <- table(TS_resist, year)
TS_resist_year
sum(TS_resist_year)
prop.table(TS_resist_year,
           margin = 2)
Test_TS_resist = chisq_test(TS_resist_year,
                            scores = list("year" = c(2016, 2017, 2018, 2019)))
Test_TS_resist


#Create dummy variables for Azithromycin non-susceptibility or resistance
res_int <- list("R", "I")
cohort <- cohort %>% mutate(AZ_not_susceptible = case_when(
  (cohort$AZ_cat %in% res_int) ~ 1,
  .default = 0))
resistant <- list("R")
cohort <- cohort %>% mutate(AZ_resistant = case_when(
  (cohort$AZ_cat %in% resistant) ~ 1,
  .default = 0))
view(cohort)
#Summarize AZ non-susceptibility, stratified by year
AZ_NS_year_tab <- do.call(rbind , by(cohort$AZ_not_susceptible, cohort$year, table)) 
AZ_NS_year_tab
AZ_NS_year_prop <- prop.table(AZ_NS_year_tab,1)
AZ_NS_year_prop
#Summarize AZ resistance, stratified by year
AZ_resist_year_tab <- do.call(rbind , by(cohort$AZ_resistant, cohort$year, table)) 
AZ_resist_year_tab
AZ_resist_year_prop <- prop.table(AZ_resist_year_tab,1)
AZ_resist_year_prop
#Run a Cochran-Armitage trend test to determine if there is a significant increase in AZ non-susceptibility
class(cohort$year)
year <- cohort$year
AZ_NS <- cohort$AZ_not_susceptible
AZ_NS_year <- table(AZ_NS, year)
AZ_NS_year
sum(AZ_NS_year)
prop.table(AZ_NS_year,
           margin = 2)
Test_AZ_NS = chisq_test(AZ_NS_year,
                        scores = list("year" = c(2016, 2017, 2018, 2019)))
Test_AZ_NS
#Check for differences in PCV13 vs. non-vax isolates
PCV13_isolates <- cohort[cohort$colonized_PCV13 == 1,]
year <- PCV13_isolates$year
AZ_NS <- PCV13_isolates$AZ_not_susceptible
AZ_NS_year <- table(AZ_NS, year)
AZ_NS_year
sum(AZ_NS_year)
prop.table(AZ_NS_year,
           margin = 2)
Test_AZ_NS = chisq_test(AZ_NS_year,
                        scores = list("year" = c(2016, 2017, 2018, 2019)))
Test_AZ_NS

Non_vax_isolates <- cohort[cohort$colonized_PCV13 == 0,]
year <- Non_vax_isolates$year
AZ_NS <- Non_vax_isolates$AZ_not_susceptible
AZ_NS_year <- table(AZ_NS, year)
AZ_NS_year
sum(AZ_NS_year)
prop.table(AZ_NS_year,
           margin = 2)
Test_AZ_NS = chisq_test(AZ_NS_year,
                        scores = list("year" = c(2016, 2017, 2018, 2019)))
Test_AZ_NS

#Run a Cochran-Armitage trend test to determine if there is a significant increase in AZ resistance
class(cohort$year)
year <- cohort$year
AZ_resist <- cohort$AZ_resistant
AZ_resist_year <- table(AZ_resist, year)
AZ_resist_year
sum(AZ_resist_year)
prop.table(AZ_resist_year,
           margin = 2)
Test_AZ_resist = chisq_test(AZ_resist_year,
                            scores = list("year" = c(2016, 2017, 2018, 2019)))
Test_AZ_resist
#Check to see if there differences in PCV13 and non-vax isolates
PCV13_isolates <- cohort[cohort$colonized_PCV13 == 1,]
year <- PCV13_isolates$year
AZ_resist <- PCV13_isolates$AZ_resistant
AZ_resist_year <- table(AZ_resist, year)
AZ_resist_year
sum(AZ_resist_year)
prop.table(AZ_resist_year,
           margin = 2)
Test_AZ_resist = chisq_test(AZ_resist_year,
                            scores = list("year" = c(2016, 2017, 2018, 2019)))
Test_AZ_resist

Non_vax_isolates <- cohort[cohort$colonized_PCV13 == 0,]
year <- Non_vax_isolates$year
AZ_resist <- Non_vax_isolates$AZ_resistant
AZ_resist_year <- table(AZ_resist, year)
AZ_resist_year
sum(AZ_resist_year)
prop.table(AZ_resist_year,
           margin = 2)
Test_AZ_resist = chisq_test(AZ_resist_year,
                            scores = list("year" = c(2016, 2017, 2018, 2019)))
Test_AZ_resist


#Evaluate non-susceptible isolates for each antibiotic
colnames(cohort)
#Amoxicillin not susceptible
table(cohort$spneumo_type, cohort$AC_not_susceptible)
prop.table(table(cohort$spneumo_type, cohort$AC_not_susceptible), 1)
table(cohort$sero_class, cohort$AC_not_susceptible)
prop.table(table(cohort$sero_class, cohort$AC_not_susceptible), 1)

#Azithromycin not susceptible
table(cohort$spneumo_type, cohort$AZ_not_susceptible)
prop.table(table(cohort$spneumo_type, cohort$AZ_not_susceptible), 1)
table(cohort$sero_class, cohort$AZ_not_susceptible)
prop.table(table(cohort$sero_class, cohort$AZ_not_susceptible), 1)

#Ceftriaxone not susceptible
table(cohort$spneumo_type, cohort$CRO_not_susceptible)
prop.table(table(cohort$spneumo_type, cohort$CRO_not_susceptible), 1)
table(cohort$sero_class, cohort$CRO_not_susceptible)
prop.table(table(cohort$sero_class, cohort$CRO_not_susceptible), 1)

#Penicillin meningitis not susceptible
table(cohort$spneumo_type, cohort$PG_men_not_susceptible)
prop.table(table(cohort$spneumo_type, cohort$PG_men_not_susceptible), 1)
table(cohort$sero_class, cohort$PG_men_not_susceptible)
prop.table(table(cohort$sero_class, cohort$PG_men_not_susceptible), 1)

#Penicillin NON-meningitis not susceptible
table(cohort$spneumo_type, cohort$PG_nonmen_not_susceptible)
prop.table(table(cohort$spneumo_type, cohort$PG_nonmen_not_susceptible), 1)
table(cohort$sero_class, cohort$PG_nonmen_not_susceptible)
prop.table(table(cohort$sero_class, cohort$PG_nonmen_not_susceptible), 1)

#TMP-SMX not susceptible
table(cohort$spneumo_type, cohort$TS_not_susceptible)
prop.table(table(cohort$spneumo_type, cohort$TS_not_susceptible), 1)
table(cohort$sero_class, cohort$TS_not_susceptible)
prop.table(table(cohort$sero_class, cohort$TS_not_susceptible), 1)


#Aim 3 Identify infant characteristics associated with colonization by antibiotic non-susceptible and resistant pneumococcal strains
#Univariate analysis of infant characteristics associated with colonization by non-susceptible strains
view(cohort)
#We'll use not-susceptible_any and resistant_any as our outcome variables
table (cohort$not_susceptible_any)

#Infant age in days vs. non-susceptible isolate
tapply(cohort$age_days, cohort$not_susceptible_any, summary)
wilcox.test(age_days ~ not_susceptible_any, data = cohort)

#Infant sex vs. non-susceptible isolate
table(cohort$not_susceptible_any, cohort$sex)
chisq.test(cohort$sex, cohort$not_susceptible_any, correct=FALSE) 
prop.table(table(cohort$sex, cohort$not_susceptible_any), 2)

#Maternal HIV status vs. non-susceptible isolate
table(cohort$not_susceptible_any, cohort$mat_hiv)
chisq.test(cohort$mat_hiv, cohort$not_susceptible_any, correct=FALSE) 
prop.table(table(cohort$mat_hiv, cohort$not_susceptible_any), 2)

#Breastfeeding status vs. non-susceptible isolate
table(cohort$not_susceptible_any, cohort$breastmilk)
chisq.test(cohort$breastmilk, cohort$not_susceptible_any, correct=FALSE) 
prop.table(table(cohort$breastmilk, cohort$not_susceptible_any), 2)

#Infant birthweight vs. non-susceptible isolate
tapply(cohort$bw, cohort$not_susceptible_any, summary)
wilcox.test(bw ~ not_susceptible_any, data = cohort)

#LBW status vs. non-susceptible isolate
table(cohort$not_susceptible_any, cohort$lbw)
chisq.test(cohort$lbw, cohort$not_susceptible_any, correct=FALSE) 
prop.table(table(cohort$lbw, cohort$not_susceptible_any), 2)

#Location of residence vs. non-susceptible isolate
table(cohort$not_susceptible_any, cohort$residence)
chisq.test(cohort$residence, cohort$not_susceptible_any, correct=FALSE) 
prop.table(table(cohort$residence, cohort$not_susceptible_any), 2)

#Year of collection - continuous vs. non-susceptible isolate
tapply(cohort$year_cont, cohort$not_susceptible_any, summary)
wilcox.test(year_cont ~ not_susceptible_any, data = cohort)
#year of collection as a categorical vs. non-susceptible isolate
table(cohort$not_susceptible_any, cohort$year)
chisq.test(cohort$year, cohort$not_susceptible_any, correct=FALSE) 
prop.table(table(cohort$year, cohort$not_susceptible_any), 1)

#Season of collection vs. non-susceptible isolate
table(cohort$not_susceptible_any, cohort$season)
chisq.test(cohort$season, cohort$not_susceptible_any, correct=FALSE) 
prop.table(table(cohort$season, cohort$not_susceptible_any), 2)

#Antibiotic exposure since prior visit vs. non-susceptible isolate
table(cohort$not_susceptible_any, cohort$inf_abx_any)
chisq.test(cohort$inf_abx_any, cohort$not_susceptible_any, correct=FALSE) 
prop.table(table(cohort$inf_abx_any, cohort$not_susceptible_any), 2)
#Amoxicillin since prior visit vs. non-susceptible isolate
table(cohort$not_susceptible_any, cohort$inf_abx_amox)
chisq.test(cohort$inf_abx_amox, cohort$not_susceptible_any, correct=FALSE) 
prop.table(table(cohort$inf_abx_amox, cohort$not_susceptible_any), 2)
#TS since prior visit vs. non-susceptible isolate
table(cohort$not_susceptible_any, cohort$inf_abx_cotrim)
chisq.test(cohort$inf_abx_cotrim, cohort$not_susceptible_any, correct=FALSE) 
prop.table(table(cohort$inf_abx_cotrim, cohort$not_susceptible_any), 2)
#Metronidazole since prior visit vs. non-susceptible isolate
table(cohort$not_susceptible_any, cohort$inf_abx_metro)
chisq.test(cohort$inf_abx_metro, cohort$not_susceptible_any, correct=FALSE) 
prop.table(table(cohort$inf_abx_metro, cohort$not_susceptible_any), 2)

#PCV doses - continuous vs. non-susceptible isolate
tapply(cohort$pcv_doses, cohort$not_susceptible_any, summary)
wilcox.test(pcv_doses ~ not_susceptible_any, data = cohort)
#PCV doses as a categorical vs. non-susceptible isolate
table(cohort$not_susceptible_any, cohort$pcv_doses)
chisq.test(cohort$pcv_doses, cohort$not_susceptible_any, correct=FALSE) 
prop.table(table(cohort$pcv_doses, cohort$not_susceptible_any), 2)

#Respiratory virus dx at visit vs. non-susceptible isolate
table(cohort$not_susceptible_any, cohort$virus_any)
chisq.test(cohort$virus_any, cohort$not_susceptible_any, correct=FALSE) 
prop.table(table(cohort$virus_any, cohort$not_susceptible_any), 2)

#Univariate analysis of respiratory virus diagnosis at time of collection
table(cohort$virus_type, cohort$not_susceptible_any)
chisq.test(cohort$virus_type, cohort$not_susceptible_any, correct=FALSE)
prop.table(table(cohort$not_susceptible_any, cohort$virus_type), 1)

#Colonized with PCV13-covered serotype vs. non-susceptible isolate
table(cohort$not_susceptible_any, cohort$colonized_PCV13)
chisq.test(cohort$colonized_PCV13, cohort$not_susceptible_any, correct=FALSE) 
prop.table(table(cohort$colonized_PCV13, cohort$not_susceptible_any), 2)

#Non-susceptibility Multivariable logistic regression
#Indicate which predictors should be treated as categorical variables
cohort$season <- factor(cohort$season)
cohort$sex <- factor(cohort$sex)
cohort$lbw <- factor(cohort$lbw)
cohort$mat_hiv <- factor(cohort$mat_hiv)
cohort$residence <- factor(cohort$residence)
cohort$breastmilk <- factor(cohort$breastmilk)
cohort$inf_abx_any <- factor(cohort$inf_abx_any)
cohort$virus_any <- factor(cohort$virus_any)
cohort$colonized_PCV13 <- factor(cohort$colonized_PCV13)
levels(cohort$sex)
cohort$sex <- relevel(cohort$sex, ref = "M")

#Run the binomial glm model
mylogit_NS <- glm(not_susceptible_any ~ age_days + year_cont + season + sex + bw + lbw + mat_hiv + residence + breastmilk + inf_abx_any + pcv_doses + virus_any + colonized_PCV13, data = cohort, family = "binomial")
summary(mylogit_NS)
confint(mylogit_NS)
exp(coef(mylogit_NS))
exp(cbind(OR = coef(mylogit_NS), confint(mylogit_NS)))
#None of the significant features survived the multivariable model


table (cohort$resist_any)

#Univariate analysis of infant characteristics associated with colonization by resistant strains
#Infant age in days vs. resistant isolate
tapply(cohort$age_days, cohort$resist_any, summary)
wilcox.test(age_days ~ resist_any, data = cohort)

#Infant sex vs. resistant isolate
table(cohort$resist_any, cohort$sex)
chisq.test(cohort$sex, cohort$resist_any, correct=FALSE) 
prop.table(table(cohort$sex, cohort$resist_any), 2)

#Maternal HIV status vs. resistant isolate
table(cohort$resist_any, cohort$mat_hiv)
chisq.test(cohort$mat_hiv, cohort$resist_any, correct=FALSE) 
prop.table(table(cohort$mat_hiv, cohort$resist_any), 2)

#Breastfeeding status vs. resistant isolate
table(cohort$resist_any, cohort$breastmilk)
chisq.test(cohort$breastmilk, cohort$resist_any, correct=FALSE) 
prop.table(table(cohort$breastmilk, cohort$resist_any), 2)

#Infant birthweight vs. resistant isolate
tapply(cohort$bw, cohort$resist_any, summary)
wilcox.test(bw ~ resist_any, data = cohort)

#LBW status vs. resistant isolate
table(cohort$resist_any, cohort$lbw)
chisq.test(cohort$lbw, cohort$resist_any, correct=FALSE) 
prop.table(table(cohort$lbw, cohort$resist_any), 2)

#Location of residence vs. resistant isolate
table(cohort$resist_any, cohort$residence)
chisq.test(cohort$residence, cohort$resist_any, correct=FALSE) 
prop.table(table(cohort$residence, cohort$resist_any), 2)

#Year of collection - continuous vs. resistant isolate
tapply(cohort$year_cont, cohort$resist_any, summary)
wilcox.test(year_cont ~ resist_any, data = cohort)
#year of collection as a categorical vs. non-susceptible isolate
table(cohort$resist_any, cohort$year)
chisq.test(cohort$year, cohort$resist_any, correct=FALSE) 
prop.table(table(cohort$year, cohort$resist_any), 1)

#Season of collection vs. resistant isolate
table(cohort$resist_any, cohort$season)
chisq.test(cohort$season, cohort$resist_any, correct=FALSE) 
prop.table(table(cohort$season, cohort$resist_any), 2)

#Antibiotic exposure since prior visit vs. resistant isolate
table(cohort$resist_any, cohort$inf_abx_any)
chisq.test(cohort$inf_abx_any, cohort$resist_any, correct=FALSE) 
prop.table(table(cohort$inf_abx_any, cohort$resist_any), 2)
#Amoxicillin since prior visit vs. resistant isolate
table(cohort$resist_any, cohort$inf_abx_amox)
chisq.test(cohort$inf_abx_amox, cohort$resist_any, correct=FALSE) 
prop.table(table(cohort$inf_abx_amox, cohort$resist_any), 2)
#TS since prior visit vs. resistant isolate
table(cohort$resist_any, cohort$inf_abx_cotrim)
chisq.test(cohort$inf_abx_cotrim, cohort$resist_any, correct=FALSE) 
prop.table(table(cohort$inf_abx_cotrim, cohort$resist_any), 2)
#Metronidazole since prior visit vs. resistant isolate
table(cohort$resist_any, cohort$inf_abx_metro)
chisq.test(cohort$inf_abx_metro, cohort$resist_any, correct=FALSE) 
prop.table(table(cohort$inf_abx_metro, cohort$resist_any), 2)

#PCV doses - continuous vs. resistant isolate
tapply(cohort$pcv_doses, cohort$resist_any, summary)
wilcox.test(pcv_doses ~ resist_any, data = cohort)
#PCV doses as a categorical vs. non-susceptible isolate
table(cohort$resist_any, cohort$pcv_doses)
chisq.test(cohort$pcv_doses, cohort$resist_any, correct=FALSE) 
prop.table(table(cohort$pcv_doses, cohort$resist_any), 2)

#Respiratory virus dx at visit vs. resistant isolate
table(cohort$resist_any, cohort$virus_any)
chisq.test(cohort$virus_any, cohort$resist_any, correct=FALSE) 
prop.table(table(cohort$virus_any, cohort$resist_any), 2)
table(cohort$resist_any, cohort$virus_type)
prop.table(table(cohort$virus_type, cohort$resist_any), 2)

#Colonized with PCV13-covered serotype vs. resistant isolate
table(cohort$resist_any, cohort$colonized_PCV13)
chisq.test(cohort$colonized_PCV13, cohort$resist_any, correct=FALSE) 
prop.table(table(cohort$colonized_PCV13, cohort$resist_any), 2)



view(cohort)
#Resistance Multivariable logistic regression
#Indicate which predictors should be treated as categorical variables
cohort$season <- factor(cohort$season)
cohort$sex <- factor(cohort$sex)
cohort$lbw <- factor(cohort$lbw)
cohort$mat_hiv <- factor(cohort$mat_hiv)
cohort$residence <- factor(cohort$residence)
cohort$breastmilk <- factor(cohort$breastmilk)
cohort$inf_abx_any <- factor(cohort$inf_abx_any)
cohort$virus_any <- factor(cohort$virus_any)
cohort$colonized_PCV13 <- factor(cohort$colonized_PCV13)
levels(cohort$sex)
cohort$sex <- relevel(cohort$sex, ref = "M")

#Run the binomial glm model
mylogit_resist <- glm(resist_any ~ age_days + year_cont + season + sex + bw + lbw + mat_hiv + residence + breastmilk + inf_abx_any + pcv_doses + virus_any + colonized_PCV13, data = cohort, family = "binomial")
summary(mylogit_resist)
confint(mylogit_resist)
exp(coef(mylogit_resist))
exp(cbind(OR = coef(mylogit_resist), confint(mylogit_resist)))
#No significant associations between infant characteristics and carriage of an antibiotic-resistant isolate.

