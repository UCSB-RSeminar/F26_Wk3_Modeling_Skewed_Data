###
#1. Load in the data and packages
###

install.packages("ggpubr")
install.packages("MASS")

library(readxl)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(moments)
library(lme4)
library(MASS)

post_clean_div_pred <- read_excel("post_clean_div_pred.xlsx")
View(post_clean_div_pred)

###
#2. Examine the skew of the data, let's look at total abundance ! #
###

# Create histograms for total abundance.This will tell us if our data is postive or negative skewed.

untransformed_total_a <- ggplot(post_clean_div_pred,aes(x=total_abundance)) +
  geom_histogram(binwidth = 2,
                 fill = 'cornflowerblue',
                 color = 'black') +
  labs(title = "Histogram of Total Abundance",
       x= "total abundance",
       y = 'Count')+
  theme_minimal()

print(untransformed_total_a)

# lets impose what a normal distribution as a red dashed line, looks like over our data, a blue solid line. 

ggarrange(
  ggdensity(
    post_clean_div_pred$total_abundance, 
    color = "blue",
    title = "untransformed total abundance") +
    stat_overlay_normal_density(
      color= "red", 
      linetype = "dashed"))



# Let's create QQ plots for total abundance.

qqnorm(post_clean_div_pred$total_abundance)
qqline(post_clean_div_pred$total_abundance)

# Let's get the skewnewss coeffient for all the continous response variables.


skewness(post_clean_div_pred$total_abundance, na.rm = TRUE)
kurtosis(post_clean_div_pred$total_abundance,na.rm = TRUE)

# Best practice is to n run a shapiro-wilk test for normality. If p>0.05, we reject the null hypothesis that the data is normally distributed.  W values close to 1 suggest normality. 
shapiro.test(post_clean_div_pred$total_abundance)

# Oh no! My data is skewed, but don't worry, we are not screwed!

#So we have determined via our histogram that the data is right skewed.
# QQ plots confirmed our data diverges from normality, skewness coefficient told us that our data is positively skewed,
# and the shapiro-wilk test confirmed our data is NOT normal.

####
####2b) Identify and remove any duplicates. Deal w/NA values...
###

# You may have noticed that in the code chunk above I added  '.na.rm = TRUE' when I called the skewness function. 
# Try removing the na.rm = True chunk and rerun the skewness function.

# Identify and remove duplicates in dataset

duplicated(post_clean_div_pred)
sum(duplicated(post_clean_div_pred))

# Deal w/NA values 
sum(is.na(post_clean_div_pred))
na_positions <- which(is.na(post_clean_div_pred), arr.ind= TRUE)
print(na_positions)

# Oopsies! There are two rows with NA values. For now let's omit them so it doesn't make our models messy. 
post_clean_div_pred <- post_clean_div_pred[-c(3,37), ]

print(post_clean_div_pred)

sum(is.na(post_clean_div_pred))


####
#### 3) What happens if we run a Gaussian linear model on this data?
####

# Predictor will be LiDAR DSM 60 ( absolute value of shrub volume loss calculated w/lidar data) and 
# response will be total abundance.


str(post_clean_div_pred)

post_clean_div_pred <- post_clean_div_pred %>%
  mutate(lidar_dsm_60 = as.numeric(as.character(lidar_dsm_60)))



model_untransformed_total_a <- lm(total_abundance ~ lidar_dsm_60, 
                                  data = post_clean_div_pred,
                                  na.action = na.omit)
summary(model_untransformed_total_a)

# So we didn't get a significant p-value, but we do see a pretty strong negative association between lidar_dsm_60 and
# total abundance.

###
### 4) Non-normal data transformation
###

#Ok, so our data is non-normal, what do we do now?

#Is it positive and moderately rigfht right-skewed? Then try ln() transformation.

#Is it positive and skewed? Try boxcox.

#Is it negative or full of zeroes? Try Yeo-Johnson transformation.

#Nothing working?!?! Switch to a non-parametric GLM, but still need to understand your data distribution...

### 4a) Because it is positive and right-skewed I am going to pick the log10 tranformation

transformed_total_a <- log10(post_clean_div_pred$total_abundance)

#Let's plot our transformed data set.


# histogram of tranformed 

transformed_total_a_hist <- ggplot(
  post_clean_div_pred,
  aes(x = log10(total_abundance))
) +
  geom_histogram(
    binwidth = 0.2,
    fill = "cornflowerblue",
    color = "black"
  ) +
  labs(
    title = "Histogram of Transformed Total Abundance",
    x = "Log10-transformed total abundance",
    y = "Count"
  ) +
  theme_minimal()

print(transformed_total_a_hist)



# let's overlay the normal distrubtution onto our tranformed data.
ggarrange(
  ggdensity(transformed_total_a, fill = "cornflowerblue", title = "log10 of total abundance")+
    stat_overlay_normal_density(color= "red", linetype = "dashed")
)



# Let's create QQ plots for all of the continuous response variables. 

qqnorm(log10(post_clean_div_pred$total_abundance))
qqline(log10(post_clean_div_pred$total_abundance))

# Let's get the skewnewss coeffient for total abundance.

skewness(transformed_total_a)
kurtosis(transformed_total_a)

# Let's run a shapiro-wilk test for normality. If p>0.05, we fail to reject H0 which states the data is non-normally distributed.
shapiro.test(transformed_total_a)

# Oof. Our data looks a lot more normal, but shapiro-wilk is telling us that it is still non-normal.


#### 4b)Let's do a boxcox transformation to find the ideal lambda.

#Note how I never ln transformed the actual column in post_clean_div_pred, I created a new object or applied the function within the function.
# Let's try a boxcox Maximum likelihood estimation

# First we fit a linear model with UNTRANSFORMED data
library(MASS)

bx_model <- lm(total_abundance ~ lidar_dsm_60, data= post_clean_div_pred)

# Then we create a plot that shows the log-likelihood curve with vertical dashed lines marking the 
# 95% CI for lambda.
# If the CI includes a 1, transformation not likely needed. 
# If the CI includes 0.5 > sqrt
# If the CI includes a 0, log transform is justified.
bxcx <- boxcox(bx_model, lambda = seq(-1, 1, by = 0.01), plotit = TRUE)

lambda_best <- bxcx$x[which.max(bxcx$y)]
lambda_best

# Our best lambda is -0.09, which we wound round to 0, justifying a log transformation.
# Let's get the confidence interval which shows  the range where 95% of lambda values will fall into.
ci_cut <- max(bxcx$y) - qchisq(0.95, 1) / 2
range(bxcx$x[bxcx$y > ci_cut])

# CI [-0.34,0.17]
#Apply the Transformation
post_clean_div_pred <- post_clean_div_pred |>
  mutate(abundance_bx = if (abs(lambda_best) < 0.01) log(total_abundance)
         else (total_abundance^lambda_best - 1) / lambda_best)

bx_model_transformed <- lm(abundance_bx ~ lidar_dsm_60, data = post_clean_div_pred)
summary(bx_model_transformed)
shapiro.test(residuals(bx_model_transformed))

# not rejecting the null hypothesis 

qqnorm(residuals(bx_model_transformed)); qqline(residuals(bx_model_transformed))

car::qqp(residuals(bx_model_transformed))




# When examining the output of the lm we get a significant negative association between total abundance and lidar dsm 60.

# The Q-Q plot looks better, but there is still some deviation from the normal distribution



###
### 5. Moving to non-parametric linear models
###

# So poission and negative binomial models work well for count data, but how do we know which one to choose?
# if the variance of your data is greater than the mean > negative binomial
# If the variance of your data is roughly equal to the mean > poisson

#5a) Fit a poisson regression model & a NB regression model to our data.

p_model <- glm(total_abundance ~ lidar_dsm_60, family = "poisson", data = post_clean_div_pred)

nb_model <- glm.nb(total_abundance ~ lidar_dsm_60, data = post_clean_div_pred)



summary(nb_model)
summary(p_model)
### 5b Gather evidence as to which statistical model is best.


#  Dharma test
  library(DHARMa)
  simulateResiduals(p_model, plot = TRUE)
  simulateResiduals(nb_model, plot = TRUE)
  
# Dispersion Test
install.packages("AER")
  library(AER)
dispersiontest(p_model)

# The dispersion tests gives us a z>1, indicating the data is likely overdispersed.
# Given our two different lines of evidence, I will choose the NB model

summary(nb_model)

#



plot <- ggplot(post_clean_div_pred, 
               aes(x = lidar_dsm_60, y = total_abundance)) +
  geom_point() +
  labs(
    title = "Relationship between Lidar Estimated Shrub Volume Loss and Total Post-fire Seedbank Abundance",
    x = "Lidar DSM 60",
    y = " Untransformed Total Seedling Abundance"
  ) +
  theme_minimal()

print(plot)

plot2 <- ggplot(post_clean_div_pred, 
               aes(x =lidar_dsm_60, y = abundance_bx)) +
  geom_point() +
  labs(
    title = "Relationship between Lidar Estimated Shrub Volume Loss and Total Post-fire Seedbank Abundance",
    x = "Lidar DSM 60",
    y = "Tranformed Total Seedling Abundance"
  ) +
  theme_minimal()

  
print(plot2)


#


