#The 3-factor model (can follow the 3-factor notebook for more details)

library(stats)
library(psych)
library(polycor)
library(lavaan)

####### NOTES 1 THE CORRELATION MATRIX ############################

#read in the coded undergraduate response data and take the items (questions) of interest
all_student_data <- read.csv('path_to_your_data_here.csv')
#student_data <- all_student_data[,c(6,7,9,15,16,17,19,21,23,24)] 
student_data <- all_student_data[,c(6,7,15,16,17,19,21,23,24)] #Q8 dropped for KMO test

#make sure the incoming data is interpreted as ordinal, not as interval - you can check that without this line, hetcor produces Pearson correlations
student_data_cat <- sapply(student_data, as.factor)

#apply the hetcor function which produces the polychoric correlation matrix and associated information
#check for 0 cell values in contingency tables
corr_info = hetcor(student_data_cat)
PCM <- corr_info$cor #the correlation matrix
PCM_type <- corr_info$type #the type of correlation computed for each pair - all should be polychoric
PCM_assumption_test <- corr_info$tests #the p-value for each test of bivariate normality (uses pearson, want < .1)

#Assessing correlation matrix as input to FA
# eyeball correlations, want some but not too large
# check Bartlett's test of sphericity - tests the null hypothesis that PCM is an identity matrix (reject with small p)
BToS_pvalue <- cortest.bartlett(PCM, n = 82)$p.value #n is the number of samples
# check for multicolinearity and singularity by checking the determinant
detMtest = det(PCM)
# check eigenvalues - should all be positive and not too close to zero (smoothing)
eigsPCM = eigen(PCM, only.values = TRUE)$values
# Kaiser-Meyer-Olkin (KMO) test of sampling adequacy 
KMO_summary <- KMO(PCM)


########## NOTES 2 EFA ##################

######### PART 2.1 DETERMINING NUMBER OF FACTORS TO RETAIN #######

#take a look at the eigenvalues of the correlation matrix and plot them  
#run the parallel analysis - which also will show you the scree plot and eigs greater than 1
fa.parallel(PCM, n.obs = 82, fm="pa", fa ="both") #note: this analysis runs 1-factor PA for many randomly simulated corr matices, it's not unusual that some warnings will be produced 

######## PART 2.2 RUNNING FACTOR ANALYSIS #######
fa_3 <- fa(r = PCM, nfactors = 3, fm = "uls", rotate = "promax")


######## PART 3 CFA ########

### 3-factor model ###

#define the model

model_3 <- 'F3_A =~ Q3_6 + Q3_14 + Q3_15 + Q3_18 + Q3_23
            F3_B =~ Q3_5 + Q3_22 + Q3_23
            F3_C =~ Q3_5 + Q3_16 + Q3_20'

#we fit our model providing the entire dataset
#ordinal data must be labelled as ordered; doing so will use the WLSMV estimator for (robust) test statistics 

fit_3 <- cfa(model_3, data = student_data, ordered = c("Q3_5","Q3_6", "Q3_14", "Q3_15", "Q3_16", "Q3_18", "Q3_20","Q3_22", "Q3_23"))

#look at some output

#summary(fit_1, fit.measures=TRUE)
sanity_3 <- standardizedSolution(fit_3) #compare factor loadings and uniquenesses to EFA soln (double check output is sensible)

#check the modifaction indices
mi <- modindices(fit_3)
#print(mi)


####### PART 4: RELIABILITY ##########
# check ordinal alpha for the item measure, and each factor 
#entire survey (since we use a polychoric corr matrix, we don't have to do anything except call alpha on the correct matrix)
overall_alpha <- alpha(PCM, check.keys=TRUE)
#find the sub-correlation matrices associated with each factor 
corr_FA <- hetcor(sapply(all_student_data[,c(7,15,16,19,24)],as.factor))$cor
corr_FB <- hetcor(sapply(all_student_data[,c(6,23,24)],as.factor))$cor
corr_FC <- hetcor(sapply(all_student_data[,c(6,17,21)],as.factor))$cor
#check their alphas
alpha_A <- alpha(corr_FA)
alpha_B <- alpha(corr_FB)
alpha_C <- alpha(corr_FC)


