# --------------------------------------------------------------------------------------------- #
# --------------------------------------- ALDH2 overview -------------------------------------- #
# programmer: Zhengting He
# date: Feburary 4th, 2021
# evidence strength assessment
# --------------------------------------------------------------------------------------------- #

require(tidyverse)

setwd("F:/Box Sync/Duke Kunshan University Intern/6 ALDH2 SR/data synthesis/meta-meta-analysis/data")
ma_data <- read.csv("data extraction table (meta-analysis)_original authors' data.csv")

# clean p-value
ma_data <- ma_data %>%
             mutate(p_value_sym = ifelse(str_detect(ma_data$p_value, "<"), "less than", NA)) %>%
               select(c(Group:ES_measurement, p_value_sym, p_value, model:calculate_method))

ma_data <- ma_data %>%
             mutate(p_value = recode(p_value, "N/A" = "NA")) %>%
               mutate(p_value = ifelse(str_detect(p_value, "<"), str_replace(p_value, "<", ""), p_value)) %>%
                 mutate(p_value = ifelse(str_detect(p_value, "-"), eval(parse(text = p_value)), p_value))   