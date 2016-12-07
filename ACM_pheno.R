library(plyr)
library(ggplot2)
library(rmarkdown)


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

ACM_raw <- read.csv("DATOS_copy.csv", header=T, sep=',')
DCM_raw <- read.csv("Paz_DCM_pheno.csv", header=T, sep=',')
HVOL_raw <- read.csv("Paz_hvol_pheno.csv", header=T, sep=',')
postPCA_FinalDataset<- read.table("Filtered_SampleIDs.txt", header=T, sep='\t')


dim(ACM_raw)
dim(DCM_raw)
dim(HVOL_raw)
dim(postPCA_FinalDataset)

### DCM - post filtering
DCM_filtered <- subset(DCM_raw, DCM_raw$Bru.Number %in% postPCA_FinalDataset$IID , )
dim(DCM_filtered)
#366  10

### HVOL - post filtering
HVOL_filtered <- subset(HVOL_raw, HVOL_raw$Bru.Number %in% postPCA_FinalDataset$IID , )
dim(HVOL_filtered)
#445   7


# ACM - no filtering
dim(ACM_raw)
#141  59


#===== Tidying Up

ACM_raw$SEX[ACM_raw$SEX == 1] <- "M" 
ACM_raw$SEX[ACM_raw$SEX == 2] <- "F" 
ACM_raw$SEX <- as.factor(ACM_raw$SEX)

#Only work with selected columns 
ACM <- subset(ACM_raw, select = c(1,7,14,32,30))
DCM <- subset(DCM_filtered, select = c(1,6,3,7,9,5))
HVOL <- subset(HVOL_filtered, select = c(1,5,3,6,7, 4))


#Add TTN Status for ACM
ACM<-cbind(group = 'ACM', ACM,Ethnicity ='C',TTN_Status =ACM_raw$VarGroup )

#Add default TTN status for DCM & HVOL - will update TTN positives later 
DCM<-cbind(group = 'DCM', DCM,TTN_Status ="DCM_TTN_Negative" )
HVOL<-cbind(group = 'HVOL', HVOL,TTN_Status ="HVOL_TTN_Status_unknown" )

head(ACM)
head(DCM)
head(HVOL)

dim(ACM)
dim(DCM)
dim(HVOL)

colnames(ACM) <- c("group", "BRU_ID", "SEX", "AGE_AT_scan", "EF","DD","Ethnicity","TTN_Status")
colnames(DCM) <- c("group", "BRU_ID", "SEX", "AGE_AT_scan", "EF","DD","Ethnicity","TTN_Status")
colnames(HVOL) <- c("group", "BRU_ID", "SEX", "AGE_AT_scan", "EF","DD","Ethnicity","TTN_Status")

#Update TTN status  for ACM 
ACM$TTN_Status <- as.character(ACM$TTN_Status)
ACM$TTN_Status[ACM$TTN_Status == "Genotype_Negative"] <- "ACM_TTN_Negative"
ACM$TTN_Status[ACM$TTN_Status == "otherDCMgenes_Positive"] <- "ACM_TTN_Negative"
ACM$TTN_Status[ACM$TTN_Status == "TTN_Positive"] <- "ACM_TTN_Positive"
ACM$TTN_Status <- as.factor(ACM$TTN_Status)

#Update TTN status  for DCM [from Nicky - Need to get Liz's list]
DCM$TTN_Status <- as.character(DCM$TTN_Status)
DCM$TTN_Status[DCM$BRU_ID == '10AC00360'] <- "DCM_TTN_Positive"
DCM$TTN_Status[DCM$BRU_ID == '10AG01687'] <- "DCM_TTN_Positive"
DCM$TTN_Status[DCM$BRU_ID == '10AH00506'] <- "DCM_TTN_Positive"
DCM$TTN_Status[DCM$BRU_ID == '10AL00811'] <- "DCM_TTN_Positive"
DCM$TTN_Status[DCM$BRU_ID == '10AP01669'] <- "DCM_TTN_Positive"
DCM$TTN_Status[DCM$BRU_ID == '10BK01803'] <- "DCM_TTN_Positive"
DCM$TTN_Status[DCM$BRU_ID == '10CH00929'] <- "DCM_TTN_Positive"
DCM$TTN_Status[DCM$BRU_ID == '10CP00605'] <- "DCM_TTN_Positive"
DCM$TTN_Status[DCM$BRU_ID == '10CP01917'] <- "DCM_TTN_Positive"
DCM$TTN_Status[DCM$BRU_ID == '10CS01472'] <- "DCM_TTN_Positive"
DCM$TTN_Status[DCM$BRU_ID == '10DH00882'] <- "DCM_TTN_Positive"
DCM$TTN_Status[DCM$BRU_ID == '10DR00069'] <- "DCM_TTN_Positive"
DCM$TTN_Status[DCM$BRU_ID == '10DW00512'] <- "DCM_TTN_Positive"
DCM$TTN_Status[DCM$BRU_ID == '10JC00262'] <- "DCM_TTN_Positive"
DCM$TTN_Status[DCM$BRU_ID == '10JF01881'] <- "DCM_TTN_Positive"
DCM$TTN_Status[DCM$BRU_ID == '10JL01453'] <- "DCM_TTN_Positive"
DCM$TTN_Status[DCM$BRU_ID == '10JM01592'] <- "DCM_TTN_Positive"
DCM$TTN_Status[DCM$BRU_ID == '10JS00403'] <- "DCM_TTN_Positive"
DCM$TTN_Status[DCM$BRU_ID == '10KF00073'] <- "DCM_TTN_Positive"
DCM$TTN_Status[DCM$BRU_ID == '10LB00067'] <- "DCM_TTN_Positive"
DCM$TTN_Status[DCM$BRU_ID == '10MK03495'] <- "DCM_TTN_Positive"
DCM$TTN_Status[DCM$BRU_ID == '10MV00215'] <- "DCM_TTN_Positive"
DCM$TTN_Status[DCM$BRU_ID == '10PC00581'] <- "DCM_TTN_Positive"
DCM$TTN_Status[DCM$BRU_ID == '10PP00413'] <- "DCM_TTN_Positive"
DCM$TTN_Status[DCM$BRU_ID == '10PP00465'] <- "DCM_TTN_Positive"
DCM$TTN_Status[DCM$BRU_ID == '10RC00179'] <- "DCM_TTN_Positive"
DCM$TTN_Status[DCM$BRU_ID == '10RN00513'] <- "DCM_TTN_Positive"
DCM$TTN_Status[DCM$BRU_ID == '10RW00786'] <- "DCM_TTN_Positive"
DCM$TTN_Status[DCM$BRU_ID == '10SR02386'] <- "DCM_TTN_Positive"
DCM$TTN_Status[DCM$BRU_ID == '10TM00933'] <- "DCM_TTN_Positive"
DCM$TTN_Status[DCM$BRU_ID == '10TM01563'] <- "DCM_TTN_Positive"
DCM$TTN_Status[DCM$BRU_ID == '10TM04444'] <- "DCM_TTN_Positive"
DCM$TTN_Status[DCM$BRU_ID == '10VC02667'] <- "DCM_TTN_Positive"
DCM$TTN_Status[DCM$BRU_ID == '10VC02667'] <- "DCM_TTN_Positive"
DCM$TTN_Status[DCM$BRU_ID == '10WS00448'] <- "DCM_TTN_Positive"
DCM$TTN_Status[DCM$BRU_ID == '11SM00049'] <- "DCM_TTN_Positive"
DCM$TTN_Status[DCM$BRU_ID == '12AM00400'] <- "DCM_TTN_Positive"
DCM$TTN_Status[DCM$BRU_ID == '12GB01392'] <- "DCM_TTN_Positive"
DCM$TTN_Status[DCM$BRU_ID == '12JL00046'] <- "DCM_TTN_Positive"
DCM$TTN_Status[DCM$BRU_ID == '12MH00114'] <- "DCM_TTN_Positive"
DCM$TTN_Status[DCM$BRU_ID == '12ML01299'] <- "DCM_TTN_Positive"
DCM$TTN_Status[DCM$BRU_ID == '12NG01428'] <- "DCM_TTN_Positive"
DCM$TTN_Status[DCM$BRU_ID == '12NP00167'] <- "DCM_TTN_Positive"
DCM$TTN_Status[DCM$BRU_ID == '12PF00041'] <- "DCM_TTN_Positive"
DCM$TTN_Status[DCM$BRU_ID == '12RH00086'] <- "DCM_TTN_Positive"
DCM$TTN_Status[DCM$BRU_ID == '12SR01222'] <- "DCM_TTN_Positive"
DCM$TTN_Status <- as.factor(DCM$TTN_Status)

#update data missing from Paz's pheno file 
DCM$AGE_AT_scan[DCM$BRU_ID == '10CH00929'] <- 71
DCM$AGE_AT_scan[DCM$BRU_ID == '10PO01022'] <- 44
DCM$AGE_AT_scan[DCM$BRU_ID == '10RM00840'] <- 66


summary(ACM)
summary(DCM)
summary(HVOL)

All<-rbind(ACM,DCM,HVOL)
dim(All)

#===== Counts
count(All, 'group')
  # group freq
# 1   ACM  141
# 2   DCM  366
# 3  HVOL  445

summary(ACM)
summary(DCM)
summary(HVOL)

All[ which(is.na(All$AGE_AT_scan) ), ]
All[ which(is.na(All$EF) ), ]
All[ which(is.na(All$DD) ), ]

#===== Sub Populations

ACM_DCM <- All[ which(All$group!='HVOL'), ]
ACM_HVOL <- All[ which(All$group!='DCM'), ]
DCM_HVOL <- All[ which(All$group!='ACM'), ]

ACM_DCM_TTN_Only <- ACM_DCM[ which(ACM_DCM$TTN_Status!='ACM_TTN_Negative'), ]
ACM_DCM_TTN_Only <- ACM_DCM_TTN_Only[ which(ACM_DCM_TTN_Only$TTN_Status!='DCM_TTN_Negative'), ]

summary(ACM_DCM)
summary(ACM_HVOL)
summary(DCM_HVOL)


#===== Exploratory Data Analysis : Box Plots & Density Plots 

Age_Box <- ggplot(All, aes(factor(group), AGE_AT_scan)) + geom_boxplot()	 #3 missing values 
Age_Density <- ggplot(ACM_DCM) + geom_density(aes(x = AGE_AT_scan,colour = group)) + ggtitle("Age")
Age_Density_TTN <- ggplot(ACM_DCM_TTN_Only) + geom_density(aes(x = AGE_AT_scan,colour = TTN_Status)) + ggtitle("Age - TTN Status")

EF_Box <- ggplot(All, aes(factor(group), EF)) + geom_boxplot()	 #12 missing values
EF_Density <- ggplot(ACM_DCM) + geom_density(aes(x = EF,colour = group)) + ggtitle("EF")
Gender <- ggplot(All, aes(x=group, fill=SEX)) + geom_bar(position="dodge", stat="count")

png(file = "ExploratoryPlots.png", bg = "transparent", height=8, width=12, units = 'in', res = 300)
multiplot(Age_Box, Age_Density,Age_Density_TTN, Gender, cols=2)
dev.off()




#===== AGE - AGE_AT_scan

summary(ACM_DCM_TTN_Only$AGE_AT_scan)
wilcox.test(AGE_AT_scan ~ TTN_Status, data=ACM_DCM_TTN_Only) 



summary(All$AGE_AT_scan)
sd(All$AGE_AT_scan, na.rm = TRUE)
summary(All[ which(All$group=='ACM'), ]$AGE_AT_scan)
sd(All[ which(All$group=='ACM'), ]$AGE_AT_scan,na.rm = TRUE)
summary(All[ which(All$group=='DCM'), ]$AGE_AT_scan)
sd(All[ which(All$group=='DCM'), ]$AGE_AT_scan,na.rm = TRUE)
summary(All[ which(All$group=='HVOL'), ]$AGE_AT_scan)
sd(All[ which(All$group=='HVOL'), ]$AGE_AT_scan,na.rm = TRUE)
wilcox.test(AGE_AT_scan ~ group, data=ACM_DCM) 
wilcox.test(AGE_AT_scan ~ group, data=ACM_HVOL) 
wilcox.test(AGE_AT_scan ~ group, data=DCM_HVOL) 



#===== SEX 

count(All, 'SEX')
count(All[ which(All$group=='ACM'), ],'SEX')
count(All[ which(All$group=='DCM'), ],'SEX')
count(All[ which(All$group=='HVOL'), ],'SEX')
ACM_vs_DCM <- matrix(c(138, 3, 255, 111), nrow = 2)
fisher.test(ACM_vs_DCM)
ACM_vs_HVOL <- matrix(c(138, 3, 201, 244), nrow = 2)
fisher.test(ACM_vs_HVOL)
DCM_vs_HVOL <- matrix(c(255, 111, 201, 244), nrow = 2)
fisher.test(DCM_vs_HVOL)
	


#===== Ethnicity
count(All, 'Ethnicity')
count(All[ which(All$group=='ACM'), ],'Ethnicity')
count(All[ which(All$group=='DCM'), ],'Ethnicity')
count(All[ which(All$group=='HVOL'), ],'Ethnicity')
ACM_vs_DCM <- matrix(c(141, 0, 366, 0), nrow = 2)
fisher.test(ACM_vs_DCM)
ACM_vs_HVOL <- matrix(c(141, 0, 445, 0), nrow = 2)
fisher.test(ACM_vs_HVOL)
DCM_vs_HVOL <- matrix(c(366, 0, 445, 0), nrow = 2)
fisher.test(DCM_vs_HVOL)


	
		 
#===== EF

summary(All$EF)
sd(All$EF, na.rm = TRUE)
summary(All[ which(All$group=='ACM'), ]$EF)
sd(All[ which(All$group=='ACM'), ]$EF,na.rm = TRUE)
summary(All[ which(All$group=='DCM'), ]$EF)
sd(All[ which(All$group=='DCM'), ]$EF,na.rm = TRUE)
summary(All[ which(All$group=='HVOL'), ]$EF)
sd(All[ which(All$group=='HVOL'), ]$EF,na.rm = TRUE)
wilcox.test(EF ~ group, data=ACM_DCM) 
wilcox.test(EF ~ group, data=ACM_HVOL) 
wilcox.test(EF ~ group, data=DCM_HVOL) 


#===== DD

summary(All$DD)
sd(All$DD, na.rm = TRUE)
summary(All[ which(All$group=='ACM'), ]$DD)
sd(All[ which(All$group=='ACM'), ]$DD,na.rm = TRUE)
summary(All[ which(All$group=='DCM'), ]$DD)
sd(All[ which(All$group=='DCM'), ]$DD,na.rm = TRUE)
summary(All[ which(All$group=='HVOL'), ]$DD)
sd(All[ which(All$group=='HVOL'), ]$DD,na.rm = TRUE)
wilcox.test(DD ~ group, data=ACM_DCM) 
wilcox.test(DD ~ group, data=ACM_HVOL) 
wilcox.test(DD ~ group, data=DCM_HVOL) 

