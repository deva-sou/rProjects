##########################################################################
#  EX1 # PAM50 intrinsic breast cancer subtypes # Projet Biostatistiques #
#  Authors : Sam MAXWELL, Enzo GOMES, Déva SOU
##########################################################################

# Objectif = Analyze the patient data according to the 50-gene expression.

data = read.csv("Data_pam50.csv", header=T, sep=",", row.names = 1)
View(data)
newdata=scale(data[,2:51]) # nomralisation des data pour le boxplot
View(newdata) 
boxplot(newdata, col="mediumslateblue") # boxplot 

# importation des librairies : 
library(FactoMineR) 
library(ade4) 
library(factoextra)

newdata2=data[2:51] # on selectionne les données à analyser

cor(newdata2) # correlation entre les variables 
# pairs(newdata2) # Marche pas 
# pcaData = dudi.pca(df = newdata2) # on regarde l'histogramme pour déterminer le nombre de dimention
pcaData = dudi.pca(df = newdata2, scannf = FALSE, nf = 2) # faire la PCA 
fviz_eig(pcaData) # tracer histogramme (plus propre avec factoextra)
histo=inertia.dudi(pcaData) 
hist(histo$tot.inertia$inertia)
inertia.dudi(pcaData, col.inertia=TRUE)

## problème Q5 : que doit on faire "Analyze the quality of the representation of the data : variables and individuals"

# cercle de correlation avec ade4
s.corcircle(pcaData$co)
pcaData$co

# graph des individus et valeurs de cos2
fviz_pca_ind(pcaData,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

# tracer le cercle de correlation propre :
fviz_pca_var(pcaData,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

# Biplot
res <- scatter(pcaData, clab.row = 0, posieig = "none")
s.class(pcaData$li, 
        data$subtype,
        col = c("#00AFBB",  "#FC4E07", "aquamarine3", "blueviolet", "darkgoldenrod1"),
        add.plot = TRUE,         # Add onto the scatter plot
        cstar = 0,               # Remove stars
        cellipse = 0             # Remove ellipses
)

grp <- as.factor(data$subtype)
fviz_mca_ind(pcaData,  
             habillage = grp,
             addEllipses = TRUE, 
             palette = c("#00AFBB",  "#FC4E07", "aquamarine3", "blueviolet", "darkgoldenrod1"), 
             label = "none"
             )

# Clustering : 

mat <- dist(newdata, method="euclidean") # Ces 2 etapes permettent de baser la CAH sur les distances du chi2
CAH <- hclust(mat, method="average")
plot(CAH)
grpecluster<-cutree(CAH,k=5)

grpecluster

res <- scatter(pcaData, clab.row = 0, posieig = "none")
s.class(pcaData$li, 
        fac = as.factor(grpecluster),
        col = c("darkgoldenrod1", "#00AFBB", "blueviolet","#FC4E07","aquamarine3"),
        add.plot = TRUE,         # Add onto the scatter plot
)


###################
# EX2 # expre.txt #
###################

data = read.table("DATA-20201123/expre.txt", header = T, sep = ",")

boxplot(data$Expression~data$Subtype,
        ylab="Expression", 
        xlab="Subtypes", 
        col=c("aquamarine3", "mediumpurple1", "mediumorchid1")
)

GrpA = data[data$Subtype=="A", "Expression"]
GrpB = data[data$Subtype=="B", "Expression"]
GrpC = data[data$Subtype=="C", "Expression"]

moyA = by(data$Subtype=="A", data$Expression)

moyA = mean(GrpA)
moyB = mean(GrpB)
moyC = mean(GrpC)
moy = c(moyA, moyB, moyC)

## Test ANOVA : 

AOV = aov(data$Expression~data$Subtype)
summary(AOV)

# D'après les résultat de l'analyse de variance, on observe une différence sugnificative entre les 3 groupes
# (p=0.0196 ; F=3.991 ; DF=2)

#  test post-hoc : 
TukeyHSD(AOV)

# Interpretation : 
##################
# L'analyse de variance montre une différence significative de l'expression des gènes entre les 3 groupes de patients
# (p=0.0196 ; F=3.991 ; DF=2 ; Analysis of Variance Model)
# Les comparaisons avec le test de Tukey montrent que deux groupent ce distinguent de par l'expression des gènes : 
# le groupe A et C  (C-A ; p=0.0258808 ; Tukey test). 
# Concernant les autres comparaisons de l'expression des gènes, le groupe B n'as pas de différence significative avec le groupe A 
# (B-A ; p=0.1486211 ; Tukey test), ainsi qu'avec le groupe C (C-B ; p=0.7363837 ; Tukey test). 
# On peut donc conclure que le groupe de patient A présente une expression des gène moins importante que le groupe de patient C. 
# En revanche on ne peut conclure entre les groupe B et A et B et C. 

c(mean(data[data$Subtype=="A", ""]))

## Test One Way : 
OW=oneway.test(data$Expression~data$Subtype)
print(OW)

# result : (p=0.0293 ; F=3.6357 ; df=2.00 ; Test for Equal Means in a One-Way Layout) 

# Kruskal-Wallis : 
KT = kruskal.test(data$Expression~data$Subtype)
print(KT)

# result : (p=0.08463 ; chi-suared=4.9389 ; df=2)

# Two Way ANOVA : 

AOV2 = aov(data$Expression~data$Subtype * data$Gender)
summary(AOV2)

# Test post-hoc : 
TukeyHSD(AOV2)

# Interpretation : 
##################
# L'analyse de variance à deux facteur met en avant une différence significative de l'expression des gènes entre les groupes de patient et le genre (masculin ou feminin) 
# (p=0.0452 ; F=3.131 ; DF=2 ; ANOVA à deux facteurs). 
# Le teste post-hoc montre que l'expression des gènes d'intéret n'est significativement différent qu'entre les patients feminins du groupe C et les patients masculins du groupe A.  
# (A:m-C:f ; p=0.0107037). 

# Version EX3 clustering : 

# groupe de gènes : 

data = read.csv("Patient.csv", header=T, sep=",")
plot(data)

K2 = kmeans(data, centers = 2)
K3 = kmeans(data, centers = 3)
K4 = kmeans(data, centers = 4)


library(factoextra)
library(ggplot2)
library(gridExtra)


plot1 <- fviz_cluster(K2, data=data) + ggtitle("k = 2")
plot2 <- fviz_cluster(K3, data=data) + ggtitle("k = 3")
plot3 <- fviz_cluster(K4, data=data) + ggtitle("k = 4")

grid.arrange(plot1, plot2, plot3, nrow = 1)

# groupe de patients : 

data2 = data.frame(t(data))

K2 = kmeans(data2, centers = 2)
K3 = kmeans(data2, centers = 3)
K4 = kmeans(data2, centers = 4)


plot1 <- fviz_cluster(K2, data=data2) + ggtitle("k = 2")
plot2 <- fviz_cluster(K3, data=data2) + ggtitle("k = 3")
# plot3 <- fviz_cluster(K4, data=data2) + ggtitle("k = 4")

grid.arrange(plot1, plot2, nrow = 1)










