rm(list=ls())
curdir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(curdir)

# Import Database information
database <- read.csv(file = "Database_Information.csv", sep = ';')
names(database)
# First, keep datasets with confidence ratings
table(database$Confidence_type)
database <- subset(database,substr(tolower(database$Confidence_type),1,17) == "confidence rating")

# Second, let's keep only datasets with confidence ratings after the decision
table(database$Conf_simultaneous_with_decision)
database <- subset(database,Conf_simultaneous_with_decision == 'no')


table(database$Confidence_scale)
