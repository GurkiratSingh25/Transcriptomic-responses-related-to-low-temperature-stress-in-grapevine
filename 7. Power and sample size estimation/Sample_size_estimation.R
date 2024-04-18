if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("RNASeqPower")


library(RNASeqPower)



#Previously, effect size per gene comparison was calculated by Cohens_d formula and power statistic per gene was calcuated using pwr.t.test()


#GREM4
#Average effect size for all DEGs
#GREM_Chill_V_Control = 5.85
#GREM_Freeze_v_Control = 5.44


# What would the power be for a study with 4 per group, to detect a
#  2-fold change, given deep (12.98x) coverage?
rnapower(12.98, cv=0.5, effect=2, n=4, alpha=.05) #0.40
rnapower(12.98, cv=0.5, effect=3, n=4, alpha=.05) #0.78
rnapower(12.98, cv=0.5, effect=4, n=4, alpha=.05) #0.93
rnapower(12.98, cv=0.5, effect=5, n=4, alpha=.05) #0.98

# Compute a sample size for several combinations of parameters
temp <- rnapower(12.98, cv=.5, effect=c(1.5, 2, 4, 5), alpha=c(.05, .01, .001),
                 power=c(.8, .9))
round(temp,1)
# Result is an array with dimensions of effect size, alpha, and power
# which contains the sample size per group for each combination.
#round(temp,1)
#, , 0.8

#0.05 0.01 0.001
#1.5 31.2 46.5  67.9
#2   10.7 15.9  23.2
#4    2.7  4.0   5.8
#5    2.0  2.9   4.3

#, , 0.9
#
#0.05 0.01 0.001
#1.5 41.8 59.2  83.2
#2   14.3 20.2  28.5
#4    3.6  5.1   7.1
#5    2.7  3.8   5.3





#CabSav

#Average effect size for all DEGs
#Cab_Chill_V_Control = 5.81
#Cab_Freeze_v_Control = 4.91

# What would the power be for a study with 4 per group, to detect a
#  2-fold change, given deep (11.21x) coverage?
rnapower(11.21, cv=0.5, effect=2, n=4, alpha=.05) #0.39
rnapower(11.21, cv=0.5, effect=3, n=4, alpha=.05) #0.76
rnapower(11.21, cv=0.5, effect=4, n=4, alpha=.05) #0.92
rnapower(11.21, cv=0.5, effect=5, n=4, alpha=.05) #0.97

# Compute a sample size for several combinations of parameters
temp <- rnapower(11.21, cv=.5, effect=c(1.5, 2, 4, 5), alpha=c(.05, .01, .001),
                 power=c(.8, .9))
round(temp,1)
# Result is an array with dimensions of effect size, alpha, and power
# which contains the sample size per group for each combination.
#round(temp,1)

#0, , 0.8

#0.05 0.01 0.001
#1.5 32.4 48.2  70.5
#2   11.1 16.5  24.1
#4    2.8  4.1   6.0
#5    2.1  3.1   4.5

#, , 0.9

#0.05 0.01 0.001
#1.5 43.4 61.4  86.3
#2   14.8 21.0  29.5
#4    3.7  5.3   7.4
#5    2.8  3.9   5.5




