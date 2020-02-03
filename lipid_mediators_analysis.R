setwd("C:/Users/Abrar/Dropbox/UNC_OneDrive/Shaikh Lab/EPA mansucript/Lipid Mediators Analysis/")
library(car)
library(ggpubr)

dt = read.csv("lipid_mediators_heart.csv", header = TRUE, fill = TRUE)
columns <- colnames(dt[,2:27]) #take all columns except for diet group (2-28 liver/adipose) & 2-27 (heart)

datalist = list()
i = 1
for(col in columns) #only loop through the column names specified above
{
  #print(col) #printing column name
  formula = as.formula(paste(col, "~ Diet")) #creating the comparison formula with the column
  
  if(shapiro.test(dt[,col])$p.value > 0.05){ #if data is normal
    if(bartlett.test(formula, dt)$p.value > 0.05){ #if the data has equal variances
      datalist[[i]] = print(compare_means(formula,  data = dt, method = "anova")) 
      res.aov <- aov(formula, data = dt)
      i = i+1
      datalist[[i]] <- TukeyHSD(res.aov) #Tukey adjusts the p-values (TukeyHSD can only be done after an ANOVA)
    } else { #if the data is normal but does NOT have equal variances
        #run the Welch ANOVA
        datalist[[i]] = print(oneway.test(formula, data = dt))
        i = i+1
        #run a Pairwise t-test with no assumption of equal variances
        datalist[[i]] = print(pairwise.t.test(dt[,col], dt$Diet, p.adjust.method = "bonferroni", pool.sd = FALSE))
      }
  } else { #if the data is not normal
      datalist[[i]] = print(compare_means(formula,  data = dt, method = "wilcox.test"))
  }
  i = i+1
}

#output results to a data file
sink('lipid_mediators_heart_statistics.txt')
print(datalist)
sink()
