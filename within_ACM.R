library(plyr)
library(ggplot2)


# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



ACM <- read.csv("DATOS_copy.csv", header=T, sep=',')


#===== Tidying Up

ACM$SEX[ACM$SEX == 1] <- "Male" 
ACM$SEX[ACM$SEX == 2] <- "Female" 

ACM<-cbind(ACM,TTN_Status =ACM$VarGroup )
ACM$TTN_Status <- as.character(ACM$TTN_Status)
ACM$TTN_Status[ACM$TTN_Status == "Genotype_Negative"] <- "TTN_Negative"
ACM$TTN_Status[ACM$TTN_Status == "otherDCMgenes_Positive"] <- "TTN_Negative"

#===== Counts

count(ACM, 'VarGroup')
                # VarGroup freq
# 1      Genotype_Negative  121
# 2 otherDCMgenes_Positive    6
# 3           TTN_Positive   14

count(ACM, 'SEX')
     # SEX freq
# 1 Female    3
# 2   Male  138


#===== Density Plots

#p1 <-ggplot(ACM, aes(DIASTOLIC_DIAMETER)) + geom_density()
#p2 <-ggplot(ACM, aes(INITIAL_EJECTION_FRACTION)) + geom_density()
#p3 <-ggplot(ACM, aes(AGE_AT_INITIAL_CLINICAL_ASSESSMENT)) + geom_density()
#p4 <-ggplot(ACM, aes(AGE_AT_ONSET_OF_SYMPTOMS)) + geom_density()

p1<-ggplot(ACM) + geom_density(aes(x = DIASTOLIC_DIAMETER,colour = VarGroup)) + ggtitle("DD")
p2<-ggplot(ACM) + geom_density(aes(x = INITIAL_EJECTION_FRACTION,colour = VarGroup)) + ggtitle("EF")
p3<-ggplot(ACM) + geom_density(aes(x = AGE_AT_INITIAL_CLINICAL_ASSESSMENT,colour = VarGroup)) + ggtitle("AGE - Consultation")
p4<-ggplot(ACM) + geom_density(aes(x = AGE_AT_ONSET_OF_SYMPTOMS,colour = VarGroup)) + ggtitle("Age - Symptoms")

png(file = "ACM_QC_density.png", bg = "transparent", height=8, width=12, units = 'in', res = 300)
multiplot(p1, p2, p3, p4, cols=2)
dev.off()	

p1<-ggplot(ACM) + geom_density(aes(x = DIASTOLIC_DIAMETER,colour = TTN_Status)) + ggtitle("DD")
p2<-ggplot(ACM) + geom_density(aes(x = INITIAL_EJECTION_FRACTION,colour = TTN_Status)) + ggtitle("EF")
p3<-ggplot(ACM) + geom_density(aes(x = AGE_AT_INITIAL_CLINICAL_ASSESSMENT,colour = TTN_Status)) + ggtitle("AGE - Consultation")
p4<-ggplot(ACM) + geom_density(aes(x = AGE_AT_ONSET_OF_SYMPTOMS,colour = TTN_Status)) + ggtitle("Age - Symptoms")

png(file = "ACM_QC_density_TTNYesNo.png", bg = "transparent", height=8, width=12, units = 'in', res = 300)
multiplot(p1, p2, p3, p4, cols=2)
dev.off()

#===== Box Plots

p1 <- ggplot(ACM, aes(factor(VarGroup), DIASTOLIC_DIAMETER)) + geom_boxplot()
p2 <- ggplot(ACM, aes(factor(VarGroup), INITIAL_EJECTION_FRACTION)) + geom_boxplot()
p3 <- ggplot(ACM, aes(factor(VarGroup), AGE_AT_INITIAL_CLINICAL_ASSESSMENT)) + geom_boxplot()
p4 <- ggplot(ACM, aes(factor(VarGroup), AGE_AT_ONSET_OF_SYMPTOMS)) + geom_boxplot()

png(file = "ACM_QC_boxPlots.png", bg = "transparent", height=8, width=12, units = 'in', res = 300)
multiplot(p1, p2, p3, p4, cols=2)
dev.off()

#===== Sub Populations

ACM_TTN_nonTTN <- ACM[ which(ACM$VarGroup!='Genotype_Negative'), ]
ACM_TTN_others <- ACM[ which(ACM$VarGroup!='otherDCMgenes_Positive'), ]
ACM_nonTTN_others <- ACM[ which(ACM$VarGroup!='TTN_Positive'), ]
	



#===== AGE - onset on symptom 

summary(ACM$AGE_AT_ONSET_OF_SYMPTOMS)
sd(ACM$AGE_AT_ONSET_OF_SYMPTOMS, na.rm = TRUE)
summary(ACM[ which(ACM$VarGroup=='TTN_Positive'), ]$AGE_AT_ONSET_OF_SYMPTOMS)
sd(ACM[ which(ACM$VarGroup=='TTN_Positive'), ]$AGE_AT_ONSET_OF_SYMPTOMS,na.rm = TRUE)
summary(ACM[ which(ACM$VarGroup=='otherDCMgenes_Positive'), ]$AGE_AT_ONSET_OF_SYMPTOMS)
sd(ACM[ which(ACM$VarGroup=='otherDCMgenes_Positive'), ]$AGE_AT_ONSET_OF_SYMPTOMS,na.rm = TRUE)
summary(ACM[ which(ACM$VarGroup=='Genotype_Negative'), ]$AGE_AT_ONSET_OF_SYMPTOMS)
sd(ACM[ which(ACM$VarGroup=='Genotype_Negative'), ]$AGE_AT_ONSET_OF_SYMPTOMS,na.rm = TRUE)
wilcox.test(AGE_AT_ONSET_OF_SYMPTOMS ~ VarGroup, data=ACM_TTN_nonTTN) 
wilcox.test(AGE_AT_ONSET_OF_SYMPTOMS ~ VarGroup, data=ACM_TTN_others) 
wilcox.test(AGE_AT_ONSET_OF_SYMPTOMS ~ VarGroup, data=ACM_nonTTN_others) 



#===== AGE - AT INITIAL CLINICAL ASSESSMENT
summary(ACM$AGE_AT_INITIAL_CLINICAL_ASSESSMENT)
sd(ACM$AGE_AT_INITIAL_CLINICAL_ASSESSMENT, na.rm = TRUE)
summary(ACM[ which(ACM$VarGroup=='TTN_Positive'), ]$AGE_AT_INITIAL_CLINICAL_ASSESSMENT)
sd(ACM[ which(ACM$VarGroup=='TTN_Positive'), ]$AGE_AT_INITIAL_CLINICAL_ASSESSMENT,na.rm = TRUE)
summary(ACM[ which(ACM$VarGroup=='otherDCMgenes_Positive'), ]$AGE_AT_INITIAL_CLINICAL_ASSESSMENT)
sd(ACM[ which(ACM$VarGroup=='otherDCMgenes_Positive'), ]$AGE_AT_INITIAL_CLINICAL_ASSESSMENT,na.rm = TRUE)
summary(ACM[ which(ACM$VarGroup=='Genotype_Negative'), ]$AGE_AT_INITIAL_CLINICAL_ASSESSMENT)
sd(ACM[ which(ACM$VarGroup=='Genotype_Negative'), ]$AGE_AT_INITIAL_CLINICAL_ASSESSMENT,na.rm = TRUE)
wilcox.test(AGE_AT_INITIAL_CLINICAL_ASSESSMENT ~ VarGroup, data=ACM_TTN_nonTTN) 
#Ties : Above gives a warning because there are less than 50 values in my TTN group, and there are ties (3 patients have same age [59])
wilcox.test(AGE_AT_INITIAL_CLINICAL_ASSESSMENT ~ VarGroup, data=ACM_TTN_others) 
wilcox.test(AGE_AT_INITIAL_CLINICAL_ASSESSMENT ~ VarGroup, data=ACM_nonTTN_others) 



#===== SEX 

count(ACM, 'SEX')
count(ACM[ which(ACM$VarGroup=='TTN_Positive'), ],'SEX')
count(ACM[ which(ACM$VarGroup=='otherDCMgenes_Positive'), ],'SEX')
count(ACM[ which(ACM$VarGroup=='Genotype_Negative'), ],'SEX')
 
ACM_TTN_nonTTN_sex_data <- matrix(c(13, 1, 6, 0), nrow = 2)
fisher.test(ACM_TTN_nonTTN_sex_data)

ACM_TTN_others_sex_data <- matrix(c(13, 1, 119, 2), nrow = 2)
fisher.test(ACM_TTN_others_sex_data)
	
ACM_nonTTN_others_sex_data <- matrix(c(6, 0, 119, 2), nrow = 2)
fisher.test(ACM_nonTTN_others_sex_data)
	
	 
	   
#===== DD 


summary(ACM$DIASTOLIC_DIAMETER)
sd(ACM$DIASTOLIC_DIAMETER, na.rm = TRUE)
summary(ACM[ which(ACM$VarGroup=='TTN_Positive'), ]$DIASTOLIC_DIAMETER)
sd(ACM[ which(ACM$VarGroup=='TTN_Positive'), ]$DIASTOLIC_DIAMETER,na.rm = TRUE)
summary(ACM[ which(ACM$VarGroup=='otherDCMgenes_Positive'), ]$DIASTOLIC_DIAMETER)
sd(ACM[ which(ACM$VarGroup=='otherDCMgenes_Positive'), ]$DIASTOLIC_DIAMETER,na.rm = TRUE)
summary(ACM[ which(ACM$VarGroup=='Genotype_Negative'), ]$DIASTOLIC_DIAMETER)
sd(ACM[ which(ACM$VarGroup=='Genotype_Negative'), ]$DIASTOLIC_DIAMETER,na.rm = TRUE)
wilcox.test(DIASTOLIC_DIAMETER ~ VarGroup, data=ACM_TTN_nonTTN) 
wilcox.test(DIASTOLIC_DIAMETER ~ VarGroup, data=ACM_TTN_others) 
wilcox.test(DIASTOLIC_DIAMETER ~ VarGroup, data=ACM_nonTTN_others) 




#===== EF 
 
 summary(ACM$INITIAL_EJECTION_FRACTION)
sd(ACM$INITIAL_EJECTION_FRACTION, na.rm = TRUE)
summary(ACM[ which(ACM$VarGroup=='TTN_Positive'), ]$INITIAL_EJECTION_FRACTION)
sd(ACM[ which(ACM$VarGroup=='TTN_Positive'), ]$INITIAL_EJECTION_FRACTION,na.rm = TRUE)
summary(ACM[ which(ACM$VarGroup=='otherDCMgenes_Positive'), ]$INITIAL_EJECTION_FRACTION)
sd(ACM[ which(ACM$VarGroup=='otherDCMgenes_Positive'), ]$INITIAL_EJECTION_FRACTION,na.rm = TRUE)
summary(ACM[ which(ACM$VarGroup=='Genotype_Negative'), ]$INITIAL_EJECTION_FRACTION)
sd(ACM[ which(ACM$VarGroup=='Genotype_Negative'), ]$INITIAL_EJECTION_FRACTION,na.rm = TRUE)
wilcox.test(INITIAL_EJECTION_FRACTION ~ VarGroup, data=ACM_TTN_nonTTN) 
wilcox.test(INITIAL_EJECTION_FRACTION ~ VarGroup, data=ACM_TTN_others) 
wilcox.test(INITIAL_EJECTION_FRACTION ~ VarGroup, data=ACM_nonTTN_others) 

