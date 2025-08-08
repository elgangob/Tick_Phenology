#Common Packages

#Tidyverse contains multiple important dataframe/graph packages
install.packages("tidyverse")
#Allows for multiple plots on a page
install.packages("gridExtra")
#Allows for publication-ready ggplots
install.packages("ggpubr")
#Better designed plots/graphs
install.packages("paletteer")
#Statistical linear fit models
install.packages("MASS")
#Statistical linear fit models
install.packages("glmmTMB")

library(tidyverse)
library(gridExtra)
library(ggpubr)
library(paletteer)
library(MASS)
library(glmmTMB)

#Importing qPCR CSV File
df <- read.csv("/Users/bharatelango/Coding/R_Studio_Code/Tick_Phenology/form-1__tick-pcr-results.csv")

#Removing unnecessary columns
df_imp <- df %>%
  subset(select = c(
    X1_Tick_ID,
    X4_Borrelia_burgdorfe,
    X5_Babesiosis,
    X6_Anaplasmosis)
  )

#Moving positive anaplasmosis results to top
df_imp <- df_imp %>%
  arrange(desc(X6_Anaplasmosis))

#Adding the Location Name Column
df_loc <- df_imp %>%
  mutate(
    Location = case_when(
      str_starts(df_imp[ , "X1_Tick_ID"], "NE12") ~ "Crow's Nest - Deep Forest",
      str_starts(df_imp[ , "X1_Tick_ID"], "NE13") ~ "Stroud South",
      str_starts(df_imp[ , "X1_Tick_ID"], "NE14") ~ "Stroud Electric",
      str_starts(df_imp[ , "X1_Tick_ID"], "NE15") ~ "Gwynedd",
      str_starts(df_imp[ , "X1_Tick_ID"], "NE16") ~ "Hildacy",
      str_starts(df_imp[ , "X1_Tick_ID"], "NE17") ~ "Cheslan South",
      str_starts(df_imp[ , "X1_Tick_ID"], "NE18") ~ "Cheslan West",
      str_starts(df_imp[ , "X1_Tick_ID"], "NE19") ~ "Crow's Nest - North",
      str_starts(df_imp[ , "X1_Tick_ID"], "MW") ~ "Midwest",
      .default = "NA"
      )
    )

#Making a count object of samples for all locations for all bacterium
loc_count_anaplasmosis <- df_loc %>%
  count(Location, name = "Total Anaplasmosis Samples") %>%
  left_join(
    df_loc %>%
      filter(X6_Anaplasmosis == "Positive") %>%
      count(Location, name = "Positive Anaplasmosis Samples"),
    by = "Location"
    ) %>%
  left_join(
    df_loc %>%
      filter(X6_Anaplasmosis == "Negative") %>%
      count(Location, name = "Negative Anaplasmosis Samples"),
    by = "Location"
  )
#Borrelia
loc_count_borrelia <- df_loc %>%
  count(Location, name = "Total Borrelia Samples") %>%
  left_join(
    df_loc %>%
      filter(X4_Borrelia_burgdorfe == "Positive") %>%
      count(Location, name = "Positive Borrelia Samples"),
    by = "Location"
  ) %>%
  left_join(
    df_loc %>%
      filter(X4_Borrelia_burgdorfe == "Negative") %>%
      count(Location, name = "Negative Borrelia Samples"),
    by = "Location"
  )
#Babesiosis
loc_count_babesiosis  <- df_loc %>%
  count(Location, name = "Total Babesiosis  Samples") %>%
  left_join(
    df_loc %>%
      filter(X5_Babesiosis == "Positive") %>%
      count(Location, name = "Positive Babesiosis  Samples"),
    by = "Location"
  ) %>%
  left_join(
    df_loc %>%
      filter(X5_Babesiosis == "Negative") %>%
      count(Location, name = "Negative Babesiosis  Samples"),
    by = "Location"
  )

#Replace NA with 0 in count df
loc_count_anaplasmosis[is.na(loc_count_anaplasmosis)] <- 0
#Borrelia
loc_count_borrelia[is.na(loc_count_borrelia)] <- 0
#Babesiosis
loc_count_babesiosis[is.na(loc_count_babesiosis)] <- 0

#Add percentage of samples column that are positive
loc_count_anaplasmosis<- loc_count_anaplasmosis %>%
  mutate(
    Prevalence = `Positive Anaplasmosis Samples`/`Total Anaplasmosis Samples`
  )
#Borrelia
loc_count_borrelia <- loc_count_borrelia %>%
  mutate(
    Prevalence = `Positive Borrelia Samples`/`Total Borrelia Samples`
  )
#Babesiosis
loc_count_babesiosis  <- loc_count_babesiosis  %>%
  mutate(
    Prevalence = `Positive Babesiosis  Samples`/`Total Babesiosis  Samples`
  )

#Sorted columns

#Anaplasma
loc_count_anaplasmosis <- loc_count_anaplasmosis %>%
  arrange(desc(Prevalence))
#Borrelia
loc_count_borrelia <- loc_count_borrelia %>%
  arrange(desc(Prevalence))
#Babesiosis
loc_count_babesiosis  <- loc_count_babesiosis  %>%
  arrange(desc(Prevalence))
  
  
#GGplots and Graphs for Count
g1 <- loc_count_anaplasmosis %>% ggplot(aes(x = Location, y = Prevalence)) +
  geom_bar(stat = "Identity", color = "black", fill = "pink")

g2 <- loc_count_babesiosis %>% ggplot(aes(x = Location, y = Prevalence)) +
  geom_bar(stat = "Identity", color = "black", fill = "pink")

g3 <- loc_count_borrelia %>% ggplot(aes(x = Location, y = Prevalence)) +
  geom_bar(stat = "Identity", color = "black", fill = "pink")

grid.arrange(g1, g2, g3, ncol = 1)

#ERI Calculations

#ERI = (# of Ticks Collected)/(Distance Dragged m^2) * (# of Positive Ana)/(# of Ticks Tested)
#Importing Field Phenology Data
df_field <- read.csv("/Users/bharatelango/Coding/R_Studio_Code/Tick_Phenology/Phenology_Field_Data - RawFlaggingData.csv")

#Important columns kept
df_field_imp <- df_field %>%
  subset(select = (c(Location, Habitat, Week.Num, Year, Species, Life.Stage, Pool.number, X..of.ticks, Sorted.Label)))

#Once distance column is filled run this!
df_field_imp <- df_field_imp %>% mutate(Distance_m2 = 400)

#I.scap Nymphs rows kept
df_field_nymph <- df_field_imp %>%
  filter(Species == 3 & Life.Stage == 2)

#Only I.scap nymphs that have been run through qPCR rows kept
df_field_nymph_curr <- df_field_nymph %>%
  filter(Sorted.Label %in% df_loc$X1_Tick_ID)


#Number of Ticks Collected


#Distance for Each Location


#Bernouli Regression Preperation
df_bern <- df_imp

#Changing Positives to 1 and Negatives to 0
for(i in 2:4) {
  df_bern[ , i] <- case_when(
    str_starts(df_bern[,i], "Pos") ~ 1,
    .default = 0
  )
}

#Creating a categorical column for regression calculation
df_bern <- df_bern %>%
  mutate(
    Location_ID = case_when(
      str_starts(df_bern[ , "X1_Tick_ID"], "NE12") ~ "NE12",
      str_starts(df_bern[ , "X1_Tick_ID"], "NE13") ~ "NE13",
      str_starts(df_bern[ , "X1_Tick_ID"], "NE14") ~ "NE14",
      str_starts(df_bern[ , "X1_Tick_ID"], "NE15") ~ "NE15",
      str_starts(df_bern[ , "X1_Tick_ID"], "NE16") ~ "NE16",
      str_starts(df_bern[ , "X1_Tick_ID"], "NE17") ~ "NE17",
      str_starts(df_bern[ , "X1_Tick_ID"], "NE18") ~ "NE18",
      str_starts(df_bern[ , "X1_Tick_ID"], "NE19") ~ "NE19",
      str_starts(df_bern[ , "X1_Tick_ID"], "MW") ~ "MW0",
      .default = NA
    )
  )

df_bern$Location_ID <- as.character(df_bern$Location_ID)

#Running the binomial GLM for each pathogen
bern_borrelia <- df_bern %>% 
  glm(X4_Borrelia_burgdorfe ~ Location_ID, ., family = binomial)

bern_babesiosis <- df_bern %>% 
  glm(X5_Babesiosis ~ Location_ID, ., family = binomial)

bern_anaplasmosis <- df_bern %>% 
  glm(X6_Anaplasmosis ~ Location_ID, ., family = binomial)

#Mixed Model Preparation
df_mix <- df_bern %>% 
  left_join(df_field_imp %>% dplyr::select(Habitat, Sorted.Label),
            by = c("X1_Tick_ID" = "Sorted.Label")) %>% 
  rename(tick_ID2 = Habitat)

#Fixing some mistakes in database
df_mix$tick_ID2 <- df_mix$tick_ID2 %>% 
  replace_na(4)

#Ordering the dataframe based on location ID then normal ID
df_mix <- df_mix %>% 
  arrange(tick_ID2, X1_Tick_ID)

df_mix$tick_ID2 <- as.character(df_mix$tick_ID2)

#Creating tick_ID2 as categorical column
df_mix <- df_mix %>% 
  left_join(df_field_imp %>% dplyr::select(Location, Sorted.Label),
            by = c("X1_Tick_ID" = "Sorted.Label"))

df_mix <- df_mix %>%
  unite("tick_ID2", Location:tick_ID2, remove = TRUE)




#Poisson
#y = total # of ticks with anaplasma
#x = total # of ticks per preserve
#dummy() to make the locations a dummy variable 


df_pois <- loc_count_anaplasmosis %>%
  mutate(
    Location_ID = case_when(
      str_starts(loc_count_anaplasmosis[ , "Location"], "Crow's Nest - Deep Forest") ~ "NE12",
      str_starts(loc_count_anaplasmosis[ , "Location"], "Stroud South") ~ "NE13",
      str_starts(loc_count_anaplasmosis[ , "Location"], "Stroud Electric") ~ "NE14",
      str_starts(loc_count_anaplasmosis[ , "Location"], "Gwynedd") ~ "NE15",
      str_starts(loc_count_anaplasmosis[ , "Location"], "Hildacy") ~ "NE16",
      str_starts(loc_count_anaplasmosis[ , "Location"], "Cheslan South") ~ "NE17",
      str_starts(loc_count_anaplasmosis[ , "Location"], "Cheslan West") ~ "NE18",
      str_starts(loc_count_anaplasmosis[ , "Location"], "Crow's Nest - North") ~ "NE19",
      str_starts(loc_count_anaplasmosis[ , "Location"], "MW") ~ "MW0",
      .default = NA
    )
  )

# pois <- df_pois %>% 
#   glm(`Positive Anaplasmosis Samples` ~ Location_ID + offset(`Total Anaplasmosis Samples`), ., family=poisson(link="log"))

# install.packages("fastDummies")
# library(fastDummies)
# 
# example <- dummy_cols(loc_count_anaplasmosis, select_columns = "Location")
