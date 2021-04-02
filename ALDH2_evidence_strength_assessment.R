# --------------------------------------------------------------------------------------------- #
# --------------------------------------- ALDH2 overview -------------------------------------- #
# programmer: Johnathan He
# date: April 1-2, 2021
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
                 mutate(p_value = ifelse(str_detect(p_value, "-"), eval(parse(text = p_value)), p_value)) %>%
                   mutate(p_value = as.numeric(p_value)) %>%
                     mutate(assess = NA)

# assessment
# class 5: all associations with P > 0.05
# class 4: all other risk factors with P < 0.05
# class 3: significant summary associations (P < 0.001) per random-effects calculation & cases > 1000
# class 2: significant summary associaions (P < 0.000001) per random-effects calculation & cases > 1000 &
#          the largest study with 95% CI excluding the null
# class 1: cases > 1000 & significant summary associations (P < 0.000001) per random-effects calculation &
#          no evidence of small-study effects & no evidence of excess of significance bias &
#          prediction intervals not including the null & largest study normally significant (P < 0.05)
#          heterogeneity not large (I2 < 50%)

ma_data <- ma_data %>%
             mutate(assess = case_when(is.na(p_value) == TRUE ~ "N/A",
                                       is.na(p_value_sym) == TRUE & p_value < 10^(-6) & N > 1000 & I2 < 50 ~ "may be class 1",
                                       p_value_sym == "less than" & p_value <= 10^(-6) & N > 1000 & I2 < 50 ~ "may be class 1",
                                       is.na(p_value_sym) == TRUE & p_value < 10^(-6) & N > 1000 & I2 >= 50 ~ "may be class 2",
                                       p_value_sym == "less than" & p_value <= 10^(-6) & N > 1000 & I2 >= 50 ~ "may be class 2",
                                       is.na(p_value_sym) == TRUE & p_value >= 10^(-6) & p_value < 10^(-3) & N > 1000 ~ "class 3",
                                       p_value_sym == "less than" & p_value > 10^(-6) & p_value <= 10^(-3) & N > 1000 ~ "class 3",
                                       is.na(p_value_sym) == TRUE & p_value >= 10^(-3) & p_value < 0.05 ~ "class 4",
                                       p_value_sym == "less than" & p_value > 10^(-3) & p_value <= 0.05 ~ "class 4",
                                       p_value > 0.05 ~ "class 5"
                                       )
             )

table(ma_data$assess)
write.csv(ma_data, file = "F:/Box Sync/Duke Kunshan University Intern/6 ALDH2 SR/data synthesis/evidence strength assessment/output/evidence strength assessment.csv")