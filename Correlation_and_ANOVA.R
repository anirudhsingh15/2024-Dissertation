library(ggstatsplot) 
library(ggplot2)
library(rstatix)

# set working directory
setwd("~/Documents/s4/Project_Data/stat")

# read the clinical data
data <- read.csv("consolidated.csv",sep = ",")
data <- data[,-c(1,3,12,31:53)]
data <- data[complete.cases(data),]

# check correlation
# serum calcium with eIS - CACTI
# Fit linear regression model
lm_model <- lm(IS_CACTI_exA ~ Serum_Ca_mg_dL, data = data)
summary(lm_model)  # Display model summary to get R-squared and p-value

# Calculate R-squared value
r_squared <- summary(lm_model)$r.squared

# Extract p-value
p_value <- coef(summary(lm_model))["Serum_Ca_mg_dL", "Pr(>|t|)"]

# Create scatter plot with trend line
ggplot(data, aes(x = Serum_Ca_mg_dL, y = IS_CACTI_exA)) +
  geom_point(color = "#ffa295") +  # Add points
  geom_smooth(method = "lm", se = FALSE, color = "black") +  # Add trend line
  
  labs(title = "Serum Calcium against eIS", subtitle = paste("r² =", round(r_squared, 7), "     " , "p-value =", format.pval(p_value)), 
       x = "Serum Calcium mg/dL", y = "eIS - CACTI") +  # Add title and axis labels
  theme_pubr() + theme(plot.title = element_text(face = "bold"))# Set the theme

# vitamin D with eIS - CACTI
# Fit linear regression model
lm_model2 <- lm(IS_CACTI_exA ~ Vit_D_ng_mL, data = data)
summary(lm_model2)  # Display model summary to get R-squared and p-value

# Calculate R-squared value
r_squared2 <- summary(lm_model2)$r.squared

# Extract p-value
p_value2 <- coef(summary(lm_model2))["Vit_D_ng_mL", "Pr(>|t|)"]

# Create scatter plot with trend line
ggplot(data, aes(x = Vit_D_ng_mL, y = IS_CACTI_exA)) +
  geom_point(color = "#ffa295") +  # Add points
  geom_smooth(method = "lm", se = FALSE, color = "black") +  # Add trend line
  
  labs(title = "Vitamin D against eIS", subtitle = paste("r² =", round(r_squared2, 7), "     " , "p-value =", format.pval(p_value2)), 
       x = "Vitamin D ng/mL", y = "eIS - CACTI") +  # Add title and axis labels
  theme_pubr() + theme(plot.title = element_text(face = "bold"))# Set the theme

##age - eIS
lmage <- lm(IS_CACTI_exA ~ Age, data = data)
summary(lmage)
r2_age <- summary(lmage)$r.squared
p_age <- coef(summary(lmage))["Age", "Pr(>|t|)"]
ggplot(data, aes(x = Age, y = IS_CACTI_exA)) +
  geom_point(color = "#00AFBB") + 
  geom_smooth(method = "lm", se = FALSE, color = "black") +  # Add trend line
  labs(title = "Distribution of Age and eIS - CACTI", 
       subtitle = paste("r² =", round(r2_age, 7), "     " , "p-value =", format.pval(p_age)),
       x = "Age", y = "eIS - CACTI") + 
  theme_pubr() + theme(plot.title = element_text(face = "bold"))# Set the theme
##bmi - eIS
lmbmi <- lm(IS_CACTI_exA ~ BMI, data = data)
summary(lmbmi)
r2_bmi <- summary(lmbmi)$r.squared
p_bmi <- coef(summary(lmbmi))["BMI", "Pr(>|t|)"]
ggplot(data, aes(x = BMI, y = IS_CACTI_exA)) +
  geom_point(color = "#E7B800") + 
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  labs(title = "Distribution of BMI and eIS - CACTI",  
       subtitle = paste("r² =", round(r2_bmi, 7), "     " , "p-value =", format.pval(p_bmi)),
       x = "BMI", y = "eIS - CACTI") + theme_pubr() + theme(plot.title = element_text(face = "bold"))# Set the theme
##hba1c - eIS
lma1c <- lm(IS_CACTI_exA ~ HbA1C, data = data)
summary(lma1c)
r2_a1c <- summary(lma1c)$r.squared
p_a1c <- coef(summary(lma1c))["HbA1C", "Pr(>|t|)"]
ggplot(data, aes(x = HbA1C, y = IS_CACTI_exA)) +
  geom_point(color = "#90D26D") + 
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  labs(title = "Distribution of HbA1C and eIS - CACTI", 
       subtitle = paste("r² =", round(r2_a1c, 7), "     " , "p-value =", format.pval(p_a1c)),
       x = "HbA1C %", y = "eIS - CACTI") + theme_pubr()  + theme(plot.title = element_text(face = "bold"))

#divide the data into terciles
#Age terciles
dfa <- data %>%
  mutate(terciles = ntile(data$Age, 3)) %>%
  mutate(Age_terciles = if_else(terciles == 1, 'Low', if_else(terciles == 2, 'Medium', 'High'))) %>%
  arrange(data$Age)
dfa

#BMI terciles
dfb <- data %>%
  mutate(terciles = ntile(data$BMI, 3)) %>%
  mutate(BMI_terciles = if_else(terciles == 1, 'Low', if_else(terciles == 2, 'Medium', 'High'))) %>%
  arrange(data$BMI)
dfb

#HbA1C terciles
dfc <- data %>%
  mutate(terciles = ntile(data$HbA1C, 3)) %>%
  mutate(HbA1C_terciles = if_else(terciles == 1, 'Low', if_else(terciles == 2, 'Medium', 'High'))) %>%
  arrange(data$HbA1C)
dfc

#vit d terciles
dfd <- data %>%
  mutate(terciles = ntile(data$Vit_D_ng_mL, 3)) %>%
  mutate(vitD_terciles = if_else(terciles == 1, 'Low', if_else(terciles == 2, 'Medium', 'High'))) %>%
  arrange(data$Vit_D_ng_mL)
dfd

#perform ANOVA for the divided terciles
#ANOVA Vitamin D terciles
##vitD and cacti
d_is_aov <- aov(IS_CACTI_exA ~ vitD_terciles, data = dfd)
summary(d_is_aov)
report(d_is_aov)
p_valuex <- summary(d_is_aov)[[1]]$`Pr(>F)`[1]

px <- ggboxplot(dfd, x = "vitD_terciles", y = "IS_CACTI_exA", color = "vitD_terciles", 
                palette = c("#00AFBB", "#E7B800","#90D26D")) +
  labs(title = "ANOVA Boxplots - Vitamin D terciles vs eIS - CACTI", subtitle = paste("p =", round(p_valuex, 3)), 
       x = "Vitamin D terciles", y = "eIS - CACTI") +
  theme_pubr() + theme(plot.title = element_text(face = "bold"))
print(px)

#ANOVA Age terciles
##age and vit d
age_d_aov <- aov(Vit_D_ng_mL ~ Age_terciles, data = dfa)
summary(age_d_aov)
report(age_d_aov)
p_value1 <- summary(age_d_aov)[[1]]$`Pr(>F)`[1]

p1 <- ggboxplot(dfa, x = "Age_terciles", y = "Vit_D_ng_mL", color = "Age_terciles", 
                palette = c("#00AFBB", "#E7B800","#90D26D")) +
  labs(title = "ANOVA Boxplots - Age terciles vs Vitamin D", subtitle = paste("p =", round(p_value1, 3)), 
       x = "Age terciles", y = "Vitamin D conc - ng/mL") +
  theme_pubr() + theme(plot.title = element_text(face = "bold"))
print(p1)

##age and ca
age_ca_aov <- aov(Serum_Ca_mg_dL ~ Age_terciles, data = dfa)
summary(age_ca_aov)
report(age_ca_aov)
p_value2 <- summary(age_ca_aov)[[1]]$`Pr(>F)`[1]

p2 <- ggboxplot(dfa, x = "Age_terciles", y = "Serum_Ca_mg_dL", color = "Age_terciles", 
                palette = c("#00AFBB", "#E7B800","#90D26D")) +
  labs(title = "ANOVA Boxplots - Age terciles vs Serum Calcium", subtitle = paste("p =", round(p_value2, 3)), 
       x = "Age terciles", y = "Serum Calcium - mg/dL") +
  theme_pubr() + theme(plot.title = element_text(face = "bold"))
print(p2)

##age and cacti
age_is_aov <- aov(IS_CACTI_exA ~ Age_terciles, data = dfa)
summary(age_is_aov)
report(age_is_aov)
p_value3 <- summary(age_is_aov)[[1]]$`Pr(>F)`[1]

p3 <- ggboxplot(dfa, x = "Age_terciles", y = "IS_CACTI_exA", color = "Age_terciles", 
                palette = c("#00AFBB", "#E7B800","#90D26D")) +
  labs(title = "ANOVA Boxplots - Age terciles vs eIS - CACTI", subtitle = paste("p =", sprintf("%.2e", p_value3)), 
       x = "Age terciles", y = "eIS - CACTI") +
  theme_pubr() + theme(plot.title = element_text(face = "bold"))
print(p3)


#ANOVA BMI terciles
##bmi and vit d
bmi_d_aov <- aov(Vit_D_ng_mL ~ BMI_terciles, data = dfb)
summary(bmi_d_aov)
report(bmi_d_aov)
p_value4 <- summary(bmi_d_aov)[[1]]$`Pr(>F)`[1]

p4 <- ggboxplot(dfb, x = "BMI_terciles", y = "Vit_D_ng_mL", color = "BMI_terciles", 
                palette = c("#00AFBB", "#E7B800","#90D26D")) +
  labs(title = "ANOVA Boxplots - BMI terciles vs Vitamin D", subtitle = paste("p =", round(p_value4, 3)), 
       x = "BMI terciles", y = "Vitamin D conc - ng/mL") +
  theme_pubr() + theme(plot.title = element_text(face = "bold"))
print(p4)

##bmi and ca
bmi_ca_aov <- aov(Serum_Ca_mg_dL ~ BMI_terciles, data = dfb)
summary(bmi_ca_aov)
report(bmi_ca_aov)
p_value5 <- summary(bmi_ca_aov)[[1]]$`Pr(>F)`[1]

p5 <- ggboxplot(dfb, x = "BMI_terciles", y = "Serum_Ca_mg_dL", color = "BMI_terciles", 
                palette = c("#00AFBB", "#E7B800","#90D26D")) +
  labs(title = "ANOVA Boxplots - BMI terciles vs Serum Calcium", subtitle = paste("p =", round(p_value5, 3)), 
       x = "BMI terciles", y = "Serum Calcium - mg/dL") +
  theme_pubr() + theme(plot.title = element_text(face = "bold"))
print(p5)

##bmi and cacti
bmi_is_aov <- aov(IS_CACTI_exA ~ BMI_terciles, data = dfb)
summary(bmi_is_aov)
report(bmi_is_aov)
p_value6 <- summary(bmi_is_aov)[[1]]$`Pr(>F)`[1]

p6 <- ggboxplot(dfb, x = "BMI_terciles", y = "IS_CACTI_exA", color = "BMI_terciles", 
                palette = c("#00AFBB", "#E7B800","#90D26D")) +
  labs(title = "ANOVA Boxplots - BMI terciles vs eIS - CACTI", subtitle = paste("p =", sprintf("%.2e", p_value6)), 
       x = "BMI terciles", y = "eIS - CACTI") +
  theme_pubr() + theme(plot.title = element_text(face = "bold"))
print(p6)


#ANOVA HbA1C terciles
##hba1c and vit d
a1c_d_aov <- aov(Vit_D_ng_mL ~ HbA1C_terciles, data = dfc)
summary(a1c_d_aov)
report(a1c_d_aov)
p_value7 <- summary(a1c_d_aov)[[1]]$`Pr(>F)`[1]

p7 <- ggboxplot(dfc, x = "HbA1C_terciles", y = "Vit_D_ng_mL", color = "HbA1C_terciles", 
                palette = c("#00AFBB", "#E7B800","#90D26D")) +
  labs(title = "ANOVA Boxplots - HbA1C terciles vs Vitamin D", subtitle = paste("p =", round(p_value7, 3)), 
       x = "HbA1C terciles", y = "Vitamin D conc - ng/mL") +
  theme_pubr() + theme(plot.title = element_text(face = "bold"))
print(p7)

##hba1c and ca
a1c_ca_aov <- aov(Serum_Ca_mg_dL ~ HbA1C_terciles, data = dfc)
summary(a1c_ca_aov)
report(a1c_ca_aov)
p_value8 <- summary(a1c_ca_aov)[[1]]$`Pr(>F)`[1]

p8 <- ggboxplot(dfc, x = "HbA1C_terciles", y = "Serum_Ca_mg_dL", color = "HbA1C_terciles", 
                palette = c("#00AFBB", "#E7B800","#90D26D")) +
  labs(title = "ANOVA Boxplots - HbA1C terciles vs Serum Calcium", subtitle = paste("p =", round(p_value8, 3)), 
       x = "HbA1C terciles", y = "Serum Calcium - mg/dL") +
  theme_pubr() + theme(plot.title = element_text(face = "bold"))
print(p8)

##hba1c and cacti
a1c_is_aov <- aov(IS_CACTI_exA ~ HbA1C_terciles, data = dfc)
summary(a1c_is_aov)
report(a1c_is_aov)
p_value9 <- summary(a1c_is_aov)[[1]]$`Pr(>F)`[1]

p9 <- ggboxplot(dfc, x = "HbA1C_terciles", y = "IS_CACTI_exA", color = "HbA1C_terciles", 
                palette = c("#00AFBB", "#E7B800","#90D26D")) +
  labs(title = "ANOVA Boxplots - HbA1C terciles vs eIS - CACTI", subtitle = paste("p =", round(p_value9, 3)), 
       x = "HbA1C terciles", y = "eIS - CACTI") +
  theme_pubr() + theme(plot.title = element_text(face = "bold"))
print(p9)


#correlation matrix and correlogram
corr <- cor(dir, method = "spearman") 
summary(corr)

c.pval.mtx <- cor_pmat(
  dir,
  vars = NULL,
  method = "spearman",
  conf.level = 0.95
)
cor_plot(
  corr,
  method = "square",
  type = "full",
  palette = NULL,
  p.mat = c.pval.mtx,
  significant.level = 0.05,
  insignificant = "blank",
  label = FALSE,
  font.label = list(size = 0.05, color = "black")
)