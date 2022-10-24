# This is the step by step way to analyze data collected in the mitigation bank#

#### PACKAGES #####
library(tidyverse)
library(nlme)
library(vegan)
library(emmeans)
library(BBmisc)
library(plsdepot)
library(metafolio)
library(stats)
library(factoextra)

###############################################################################
# Each subplot is an experimental unit, so when there are multiple data 
# collections within a subplot (i.e. 2 transects in each), data needs to be 
# averaged so that there is only one value per subplot to avoid pseudoreplication.
# This type of analysis needs to be done for each variable before it can be 
# analyzed. 

# Below is an example of code that can be used to average data for plant data, 
# called 

plantdata <- read.csv("https://raw.githubusercontent.com/michaelaw18/mitigation_bank/main/2022_mitbank_plantdata.csv")
head(plantdata)

# In this data file we have
# 1. sample: the subplot, transect, and sampling location
# 2. location: subplot
# 3. richness through bgcov - different measurements of plant community
# 4. date: the julian date of data collection

# Now we want to average the plant community measurements per location (subplot)
# since that is the experimental unit at the mitigation bank.

plantdata_fin <- plantdata %>% #pulling from the plantdata csv we imported
  group_by(location) %>%  #location is what we want to average by
  summarise_at(vars(richness, grassrich, grasscov, forbrich, #summarise for each plant community metric
                    forbcov,invasrich,invascov,woodyrich,woodycov,bgcov,date),
               list(name=mean))

for ( col in 2:ncol(plantdata_fin)){
  colnames(plantdata_fin)[col] <-  sub("_name", "", colnames(plantdata_fin)[col])
}

head(plantdata_fin)

# Now that there is one value per experimental unit (subplot) this can be added 
# to the master data sheet to be compared with other variables such as 
# soil nutrients and enzymes. You can export this to excel as a csv by using
# the write.csv() function.

################################################################################
################################################################################
# THE MASTER DATA SHEET
###############################################################################

mbdf <- read.csv("https://raw.githubusercontent.com/michaelaw18/mitigation_bank/main/woods22_master_mitbank_file.csv")
head(mbdf)

# This data sheet is organized by year_subplot since each experimental unit is
# based on the subplot, and we've repeated collection by year. This datasheet
# has data from 2019, 2020, 2021 and 2022, collected by Michaela Woods.

# Variables in the dataset 
# 1. year_subplot: one column per subplot in each year
# 2. year: year of data collection
# 3. treatment: a combination of soil and seed treatments
#      h5 - high diversity, low legume
#      h20 - high diversity, high legume
#      l5 - low diversity, low legume
#      l20 - low diversity, high legume
#      ws - whole soil
#      lc - leaf compost
#      wslc - whole soil and leaf compost
#      ctrl - no soil treatment
# 4. subplot: subplot location without year
# 5. tx: soil treatment
#      A: mulch and whole soil
#      B: mulch
#      C: whole soil
#      D: control - no treatment
# 6. seed: the seed mix applied
#      h5 - high diversity, low legume
#      h20 - high diversity, high legume
#      l5 - low diversity, low legume
#      l20 - low diversity, high legume
# 7-18. plant community metrics
# 19. Julian date of collection (doesn't include year)
# 20-22: soil enzymes
#      po - phenol oxidase
#      bg - beta glucosidase
#      per - peroxidase
# 23-41: soil chemistry and nutrients
# 42. tree_ht: the growth rate of tree height (cm) from May 2022 - September 2022
# 43. tree_bd: the growth rate of tree basal diameter (mm) from May 2022 - September 2022

# The main independent variables are year and treatment (soil tx and seed mix)
# soil and plant metrics can be compared or used as independent variables as well
# depending on the question being asked 


##################################################################################
### DATA ANALYSIS WITHIN THIS EXPERIMENTAL DESIGN
#################################################################################

# COMPARING YEARS

# To assess a dependent variable as a function of year, we have to consider
# whether treatment should be included within the model. Below are several 
# questions we can ask about time using peroxidase activity as a response 
# variable

# Question 1: How does peroxidase activity change over time within all treatments?

mbdf$year <- as.factor(mbdf$year) #make sure year is read as a categorical variable

summary(aov(per~year, data=mbdf)) #use ANOVA to determine if differences in groups occur
TukeyHSD(aov(per~year, data=mbdf)) #use Tukey test to determine where differences occur

#Question 2: How does peroxidase activity change over time and with treatment?

summary(aov(per~year*treatment, data=mbdf)) 
TukeyHSD(aov(per~year*treatment, data=mbdf))
    # This creates a huge output list since there are a ton of combinations of 
    # treatments and years. It may make more sense to subset certain years or
    # only consider seed mix OR soil treatment

# Question 3: How does peroxidase activity change over time with soil treatment
#   regardless of seed mix used?

    # Here we are going to use seed mix as a random effect, this means that we
    # will use a linear mixed effects model that eliminates the variation 
    # in peroxidase activity by seed mix so that we are only seeing variation
    # that would be attributed to either year or soil. 

summary(lme(per~year*tx, random=~1|seed, data=mbdf, na.action=na.omit))
lsmeans(lme(per~year*tx, random=~1|seed, data=mbdf, na.action=na.omit), 
        pairwise~year*tx, adjust= "Tukey")

    # This still generates a lot of combinations but it is more digestible.
    # You can swap the placement of tx and seed to consider the effect of 
    # seed with soil treatment as a random effect

###############################################################################

# COMPARING TREATMENTS WITHIN ONE YEAR

# First, we pick a year and subset the data

mbdf22 <- filter(mbdf, year=="2022")

# Question 1: Did peroxidase depend on soil and seed mix treatments?

summary(aov(per~treatment, data=mbdf22)) 
TukeyHSD(aov(per~treatment, data=mbdf22))

# Question 2: Did peroxidase depend on soil treatment alone, without variation
#     from seed mix?

summary(lme(per~tx, random=~1|seed, data=mbdf22, na.action=na.omit))
lsmeans(lme(per~tx, random=~1|seed, data=mbdf22, na.action=na.omit), 
        pairwise~tx, adjust= "Tukey")

################################################################################
# MULTIVARIATE ANALYSIS
# If we want to look at all of the enzyme activities by year, we can use a principal
# component analysis

# Step 1: normalize data
ppo<-normalize(mbdf$po, method ="range", range = c(0, 1), margin = 1L, on.constant = "quiet")
pper<-normalize(mbdf$per, method ="range", range = c(0, 1), margin = 1L, on.constant = "quiet")
pbg<-normalize(mbdf$bg, method ="range", range = c(0, 1), margin = 1L, on.constant = "quiet")

#Create dataframe with normalized variables and grouping variable (year)
mbdf.enz<-with(mbdf, data.frame(year,ppo, pper, pbg))
mbdf.enz<-na.omit(mbdf.enz)

#Using ML, group enzyme activities based on similarities
enzi.cor<-nipals(mbdf.enz[,2:4])
plot(enzi.cor,main="Circle of Correlation", cex.main=1)
cor(mbdf.enz[,2:4])

#PLOTTING

#Create the PCA object
os.pca<-prcomp(~ppo+pper+pbg, mbdf.enz[,2:4], na.action=na.omit, center=TRUE,scale.=TRUE)

#Scree plot to see the percentage of explained variance of each dimension
fviz_eig(os.pca)

#create PCA plot
fviz_pca_biplot(os.pca, #pca made with prcomp
                col.ind=mbdf.enz$year, #grouping value
                geom=c("point"), #makes points unlabelled
                addEllipses = TRUE, #creates 95% confidence intervals
                legend.title = "Year", #gives title to group legend
                repel= TRUE) #reduces overlap

#Analysis of similarity (ANOSIM) to determine if grouping is statistically significant
enz<-(mbdf.enz[,2:4])
enz1<-na.omit(enz)
anosim(enz1, mbdf.enz$year, permutations = 999, distance = "euclidean", strata=NULL, parallel=1)

#################################################################################
# Use non-metric multidimensional scaling to compare plant communities over time

#make dataframe including 2020, 2021, and 2022 
mbdfveg <- filter(mbdf, year != "2019")
head(mbdfveg)

#Select community metrics and grouping category 
mbdfvegnmdsyear<-select(mbdfveg, year, grassrich, forbrich, woodyrich, invasrich)
head(mbdfvegnmdsyear)

#create nmds
yearvegnmds<-mbdfvegnmdsyear[2:5]
set.seed(2)
fyearvegnmds <- metaMDS(comm=yearvegnmds, distance="bray", trace =F, trymax=999)
fyearvegnmds

#anosim to see if grouping by year is significant
anosim(yearvegnmds, mbdfvegnmdsyear$year)

#plot

yearveghnmdsscores <- as.data.frame(scores(fyearvegnmds))
yearvegnmdsgraph <- data.frame(yearveghnmdsscores, NMDS1=yearveghnmdsscores$NMDS1, NMDS2=yearveghnmdsscores$NMDS2, year=mbdfvegnmdsyear$year)

ggplot(yearvegnmdsgraph, aes(NMDS1, NMDS2, colour=year))+
  stat_ellipse(data=yearvegnmdsgraph,aes(x=NMDS1, geom="path", position="identity", level=0.95, y=NMDS2, colour=year),size=1, linetype=1)+
  geom_point(aes(shape=year), size=3.5)
