setwd("/Macro")
library(geomorph)
library(RRPP)
library(rgl)
library(ape)
library(geiger)
library(phytools)
library(png)

######################MEAN SHAPE and SIZE####################
#mean shape calculation

Toothedwhales_phy<-read.morphologika("TW_phy_Rep.txt")
CV_Toothedwhales<-read.csv("TW_phy_Rep.csv",sep=",")
head(CV_Toothedwhales)

ToothedWhales_phy.gpa<-gpagen(Toothedwhales_phy)
Toothedwhales_gdf<- geomorph.data.frame(Shape=ToothedWhales_phy.gpa$coords,Csize=ToothedWhales_phy.gpa$Csize,Species=CV_Toothedwhales$Species2)
Shape<-two.d.array(ToothedWhales_phy.gpa$coords)
shape.means<-rowsum(Shape,Toothedwhales_gdf$Species)/as.vector(table(Toothedwhales_gdf$Species))

#mean CS calculation

Cs<-Toothedwhales_gdf$Csize
Cs.mean<-rowsum(Cs,Toothedwhales_gdf$Species)/as.vector(table(Toothedwhales_gdf$Species))
logCs<-log10(Cs.mean)

#########################TREE:DROP.TIP################

##Drop branches in MCGowen et al 2009
setwd("/Online Resource 1")
McGowen2009<-read.nexus("McGowen2009.nex")
plot(McGowen2009,cex=0.5)


#taxa I need
Myclade=c("Be_arnuxii",
          "Ce_commersonii",
          "Ce_eutropia",
          "Ce_heavisidii",
          "De_leucas",
          "De_capensis",
          "De_delphis",
          "De_tropicalis",
          "Fe_attenuata",
          "Gl_macrorhynchus",
          "Gl_melas",
          "Gr_griseus",
          "Hy_ampullatus",
          "Hy_planifrons",
          "In_pacificus",
          "In_geoffrensis",
          "Ko_breviceps",
          "Ko_simus",
          "La_hosei",
          "La_albirostris",
          "La_australis",
          "La_obscurus",
          "La_acutus",
          "Li_vexillifer",
          "Li_borealis",
          "Li_peronii",
          "Me_bidens",
          "Me_bowdoini",
          "Me_europaeus",
          "Me_ginkgodens",
          "Me_hectori",
          "Me_mirus",
          "Mo_monoceros",
          "Ne_phocaenoides",
          "Orcaella_brevirostris",
          "Orcaella_heinsohni",
          "Or_orca",
          "Pe_electra",
          "Au_dioptrica",
          "Ph_phocoena",
          "Ph_sinus",
          "Ph_spinipinnis",
          "Ph_dalli",
          "Pl_gangetica",
          "Po_blainvillei",
          "Ps_crassidens",
          "La_cruciger",
          "La_obliquidens",
          "Sousa_chinensis",
          "Sousa_plumbea",
          "Sousa_teuszii",
          "So_fluviatilis",
          "St_attenuata",
          "St_coeruleoalba",
          "St_frontalis",
          "St_longirostris",
          "St_bredanensis",
          "Tu_aduncus",
          "Tu_truncatus",
          "Zi_cavirostris")

#tip.to:
seq(1:Ntip(McGowen2009))[-which(McGowen2009$tip.label%in%Myclade)]->tip.to

#"tip.to": displays taxa I do not need
tip.to

#function drop tip to emilinate taxa I do not need
phy<-drop.tip(McGowen2009, tip.to)
plot(phy,cex=0.5)

#change tip labels name in new phy as my coords file
phy$tip.label<-c("Delphinus tropicalis",
                 "Delphinus delphis",
                 "Delphinus capensis",
                 "Tursiops truncatus",
                 "Tursiops aduncus",
                 "Stenella coeruleoalba",
                 "Stenella frontalis",
                 "Stenella longirostris",
                 "Lagenodelphis hosei",
                 "Stenella attenuata",
                 "Sousa teuszii",
                 "Sousa chinensis",
                 "Sousa plumbea",
                 "Sotalia fluviatilis",
                 "Steno bredanensis",
                 "Feresa attenuata",
                 "Globicephala melas",
                 "Globicephala macrorhynchus",
                 "Peponocephala electra",
                 "Pseudorca crassidens",
                 "Grampus griseus",
                 "Orcaella brevirostris",
                 "Orcaella heinsohni",
                 "Lissodelphis peronii",
                 "Lissodelphis borealis",
                 "Lagenorhynchus obscurus",
                 "Lagenorhynchus obliquidens",
                 "Lagenorhynchus cruciger",
                 "Lagenorhynchus australis",
                 "Cephalorhynchus heavisidii",
                 "Cephalorhynchus eutropia",
                 "Cephalorhynchus commersoni",
                 "Lagenorhynchus albirostris",
                 "Lagenorhynchus acutus",
                 "Orcinus orca",
                 "Phocoena sinus",
                 "Phocoena spinipinnis",
                 "Phocoena dioptrica",
                 "Phocoena phocoena",
                 "Phocoenoides dalli",
                 "Neophocaena phocaenoides",
                 "Delphinapterus leucas",
                 "Monodon monoceros",
                 "Inia geoffrensis",
                 "Pontoporia blanvillei",
                 "Lipotes vexillifer",
                 "Mesoplodon bowdoini",
                 "Mesoplodon ginkgondens",
                 "Mesoplodon europaeus",
                 "Mesoplodon mirus",
                 "Mesoplodon hectori",
                 "Mesoplodon bidens",
                 "Indopacetus pacificus",
                 "Hyperodont planifrons",
                 "Hyperodont ampullatus",
                 "Ziphius cavirostris",
                 "Berardius arnuxii",
                 "Platanista gangetica",
                 "Kogia breviceps",
                 "Kogia sima")
plot(phy,cex = 0.5)


#export new pruned tree
setwd("/Online Resource 1")
write.tree(phy,"McGowen2009_pruned.txt")

#########################TREE SHAPE MATCH############
#matching tree and shape data
#name check

td<-treedata(phy,shape.means,sort=TRUE)
shape.data=td$data#this is the one that needs to be saved
tree.matched=td$phy

name.check(phy, shape.means)#OK

#export mean shape 
TW_ms<-arrayspecs(shape.data,28,3)
dimnames(TW_ms)
dim(TW_ms)#28 3 60

# write.morphologika
# authors:
#   João Coelho
#   David Navega
#   Department of Life Sciences
#   University of Coimbra

setwd("/Online Resource 1")
write.morphologika(TW_ms,"TW_ms.txt")

##files created(McGowen2009_pruned.txt, TW_ms.txt) 
##to be used in next analyses

####################################PHYSIGNAL####################


# Calculate measures of phylogenetic signal(Kmult)

PS.shape.Toothedwhales <- physignal(shape.means,phy=phy,print.progress = T)
PS.shape.Toothedwhales$phy.signal
PS.shape.Toothedwhales$pvalue

#kmult= k (Blomberg,2003)
PS.Cs.Toothedwhales <- physignal(logCs,phy=phy, print.progress = T)
PS.Cs.Toothedwhales$phy.signal
PS.Cs.Toothedwhales$pvalue


#########################PCA##############################

Y.gpa<-gpagen(TW_ms)
ref<-mshape(Y.gpa$coords)
findMeanSpec(Y.gpa$coords)

#import 3D model
setwd("/Online Resource 1")
mesh<-read.ply("Tt_1984-1757_NHM_ASCII-dec.ply")

#coords of the specimen in the Replicate coords file
mesh.coord<-Toothedwhales_phy[,,217]
refmesh<-warpRefMesh(mesh,mesh.coord,ref, color = NULL)

#import variables
TW_phy.csv<-read.csv("TW_parameters_phy.csv",sep=",")

#name row as tre and coords data
rownames(TW_phy.csv)<-TW_phy.csv[,2]
head(TW_phy.csv)

#match variables with tree and coords data
TW_phy.csv[match(rownames(TW_phy.csv),tree.matched$tip.label),]

TW_gdf<-geomorph.data.frame(Shape=TW_ms,CS=TW_phy.csv$Cs,Species=TW_phy.csv$Species2,Family=TW_phy.csv$Family,
                            Biosonar=TW_phy.csv$Biosonar..Jensen.et.al.2018.,Diet=TW_phy.csv$Diet,Ecology=TW_phy.csv$ECOLOGY,
                            SST=TW_phy.csv$SST,L=TW_phy.csv$Log10.Length.cm..,BM=TW_phy.csv$Log10.Body.mass.g..)

#PCA on whole mean shape:Figure 3

TW_gp<-as.factor(paste(TW_gdf$Family))
TW_gp2<- as.factor(paste(TW_gdf$Species))
TW_col.gp<-rainbow(length(levels(TW_gp))) 
names(TW_col.gp) <-levels(TW_gp)
TW_col.gp2<-TW_col.gp[match(TW_gp, names(TW_col.gp))]

label<-c("Dt",
         "Dd",
         "Dc",
         "Tt",
         "Ta",
         "Sco",
         "Sfr",
         "Sl",
         "Lh",
         "Sa",
         "St",
         "Sch",
         "Sp",
         "Sf",
         "Sb",
         "Fa",
         "Gm",
         "Gma",
         "Pe",
         "Pc",
         "Gg",
         "Ob",
         "Oh",
         "Lip",
         "Lib",
         "Laobs",
         "Laobl",
         "Lac",
         "Laau",
         "Ch",
         "Ce",
         "Cc",
         "Laal",
         "Laac",
         "Oo",
         "Ps",
         "Psp",
         "Pdi",
         "Pp",
         "Pda",
         "Np",
         "Dl",
         "Monmon",
         "Ig",
         "Pb",
         "Lv",
         "Mbo",
         "Mg",
         "Me",
         "Mm",
         "Mh",
         "Mbi",
         "Ip",
         "Hp",
         "Ha",
         "Zc",
         "Ba",
         "Pg",
         "Kb",
         "Ks")

TW.pca<-plotTangentSpace(Y.gpa$coords,groups=TW_col.gp2,pch=length(levels(TW_gp)),label=label,legend = F,lwd=1)
#PC2*-1 to reverse PC2

##3D warping

#shape change along PC1
Pc1min<-plotRefToTarget(M1=ref,M2=TW.pca$pc.shapes$PC1min,mesh=refmesh,method="surface")
PC1max<-plotRefToTarget(M1=ref,M2=TW.pca$pc.shapes$PC1max,mesh=refmesh,method="surface")
##shape change along PC2
PC2min<-plotRefToTarget(M1=ref,M2=TW.pca$pc.shapes$PC2min,mesh=refmesh,method="surface")
PC2max<-plotRefToTarget(M1=ref,M2=TW.pca$pc.shapes$PC2max,mesh=refmesh,method="surface")
#plot
plotspec(Pc1min,digitspec =TW.pca$pc.shapes$PC1min,box=F,axes=F,type="shade")
plotspec(PC1max,digitspec =TW.pca$pc.shapes$PC1max,box=F,axes=F,type="shade")
plotspec(PC2min,digitspec =TW.pca$pc.shapes$PC2min,box=F,axes=F,type="shade")
plotspec(PC2max,digitspec =TW.pca$pc.shapes$PC2max,box=F,axes=F,type="shade")
PCA_summary<-TW.pca$pc.summary$importance

setwd("/Online Resource 1")
#write.table(PCA_summary,"PCsummary.xls",sep ="\t")

##PCA: Figure 4
#plot three graphs
par(mfrow=c(3,1))

xlab <-paste("PC1","(",round(TW.pca$pc.summary$importance[2,1]*100, 1),"%)", sep="")
ylab <- paste("PC2","(", round(TW.pca$pc.summary$importance[2,2]*100, 1), "%)", sep="")

#PCA Diet A
TW_gpDiet<-as.factor(paste(TW_gdf$Diet))
TW_gp2Diet<- as.factor(paste(TW_gdf$Species))
TW_col.gpDiet<-rainbow(length(levels(TW_gpDiet))) 
names(TW_col.gpDiet) <-levels(TW_gpDiet)
TW_col.gp2Diet<-TW_col.gpDiet[match(TW_gpDiet,names(TW_col.gpDiet))]

#PC2*-1 to reverse PC2
plot(TW.pca$pc.scores[,1],TW.pca$pc.scores[,2],cex=2,pch=21,bg=TW_col.gp2Diet,lwd=1,xlab=xlab, ylab=ylab, asp=T,main="Diet",
     text(TW.pca$pc.scores[,1:2],labels=label,pos=4))

# Add convex hull polygons to the PCA plot:
colour = c("red", "green","blue") # colour for the three groups
for(j in 1:nlevels(TW_gpDiet)) {
  # Get edge points (used to plot convex hull):
  edge_points <- rownames(TW.pca$pc.scores[which(TW_gpDiet == levels(TW_gpDiet)[j]),])[
    chull(TW.pca$pc.scores[which(TW_gpDiet == levels(TW_gpDiet)[j]), c(1,2)])]
  # Plot convex hull as polygon:
  polygon(TW.pca$pc.scores[edge_points, c(1,2)], col = adjustcolor(colour[j],
                                                                   alpha.f = 0.3) , border = colour[j])
} # alpha gives the degree of transparency of the polygon

legend("topright",pch=21,legend=c("Fish/Mammals", "Squid","Fish"),cex = 0.7,col=c("#00FF00FF", "#0000FFFF","#FF0000FF"),lwd=4)

#PCA Ecology B
TW_gpEcology<-as.factor(paste(TW_gdf$Ecology))
TW_gp2Ecology<- as.factor(paste(TW_gdf$Species))
TW_col.gpEcology<-rainbow(length(levels(TW_gpEcology))) 
names(TW_col.gpEcology) <-levels(TW_gpEcology)
TW_col.gp2Ecology<-TW_col.gpEcology[match(TW_gpEcology,names(TW_col.gpEcology))]

#PC2*-1 to reverse PC2
plot(TW.pca$pc.scores[,1],TW.pca$pc.scores[,2],cex=2,pch=21,bg=TW_col.gp2Ecology,lwd=1,xlab=xlab, ylab=ylab, asp=T,main="Ecology",
     text(TW.pca$pc.scores[,1:2],labels=label,pos=4))

# Add convex hull polygons to the PCA plot:
for(j in 1:nlevels(TW_gpEcology)) {
  # Get edge points (used to plot convex hull):
  edge_points <- rownames(TW.pca$pc.scores[which(TW_gpEcology == levels(TW_gpEcology)[j]),])[
    chull(TW.pca$pc.scores[which(TW_gpEcology == levels(TW_gpEcology)[j]), c(1,2)])]
  # Plot convex hull as polygon:
  polygon(TW.pca$pc.scores[edge_points, c(1,2)], col = adjustcolor(colour[j],
                                                                   alpha.f = 0.3) , border = colour[j])
} # alpha gives the degree of transparency of the polygon

#PCA Biosonar C
TW_gpBio<-as.factor(paste(TW_gdf$Biosonar))
TW_gp2Bio<- as.factor(paste(TW_gdf$Species))
TW_col.gpBio<-rainbow(length(levels(TW_gpBio))) 
names(TW_col.gpBio) <-levels(TW_gpBio)
TW_col.gp2Bio<-TW_col.gpBio[match(TW_gpBio,names(TW_col.gpBio))]

#PC2*-1 to reverse PC2
plot(TW.pca$pc.scores[,1],TW.pca$pc.scores[,2],cex=2,pch=21,bg=TW_col.gp2Bio,lwd=1,xlab=xlab, ylab=ylab, asp=T,main="Biosonar",
     text(TW.pca$pc.scores[,1:2],labels=label,pos=4))

# Add convex hull polygons to the PCA plot:
for(j in 1:nlevels(TW_gpBio)) {
  # Get edge points (used to plot convex hull):
  edge_points <- rownames(TW.pca$pc.scores[which(TW_gpBio == levels(TW_gpBio)[j]),])[
    chull(TW.pca$pc.scores[which(TW_gpBio == levels(TW_gpBio)[j]), c(1,2)])]
  # Plot convex hull as polygon:
  polygon(TW.pca$pc.scores[edge_points, c(1,2)], col = adjustcolor(colour[j],
                                                                   alpha.f = 0.3) , border = colour[j])
} # alpha gives the degree of transparency of the polygon

##################################Figure 2#####
###Figure 2 

#GroupSize
CentroidSize<-as.matrix(TW_phy.csv)[,18]
mode(CentroidSize)<-'numeric'#should be labeled as numeric
#str(TW_phy.csv$Cs)#continuous variable should be labeled as numeric
#CentroidSize<-as.numeric(CentroidSize)
#CentroidSize<-CentroidSize[tree.matched$tip.label]#assigning cranium size as a tree label
obj<-contMap(tree.matched,CentroidSize,outline=TRUE,lwd=2,fsize=.6);axisPhylo()

par(mfrow=c(1,1))
plot(obj,direction="rightwards",type="fan", fsize=c(0.5,1),lwd=2,outline=T,palette="gray",leg.txt="CS")


#adding images to our phylogeny
#create a vector with all the png names that need to be imported in the tree


spp2<-c("Delphinus delphis",
        "Tursiops truncatus",
        "Stenella coeruleoalba",
        "Stenella longirostris",
        "Lagenodelphis hosei",
        "Sousa plumbea",
        "Steno bredanensis",
        "Feresa attenuata",
        "Globicephala melas",
        "Peponocephala electra",
        "Pseudorca crassidens",
        "Grampus griseus",
        "Orcaella heinsohni",
        "Lissodelphis peronii",
        "Lagenorhynchus obliquidens",
        "Cephalorhynchus commersoni",
        "Cephalorhynchus heavisidii",
        "Lagenorhynchus albirostris",
         "Orcinus orca",
        "Phocoena spinipinnis",
        "Phocoena phocoena",
        "Delphinapterus leucas",
        "Monodon monoceros",
        "Inia geoffrensis",
        "Pontoporia blanvillei",
        "Lipotes vexillifer",
        "Mesoplodon europaeus",
        "Indopacetus pacificus",
        "Hyperodont planifrons",
        "Ziphius cavirostris",
        "Berardius arnuxii",
        "Platanista gangetica",
          "Kogia sima")

##add them at each tip. tree tip must have same name of the png image.
#to display arrow change col 
for(i in 1:length(spp2)){
  xy<-add.arrow(obj,spp2[i],col="transparent",arrl=7,lwd=3,hedl=10)
  arrl2<-if(spp2[i]%in%c("Orcinus orca","Tursiops truncatus"))101
  else if(spp2[i]=="Globicephala melas") 0.1302674 else 89
  add.arrow(obj,spp2[i],col="transparent",arrl=arrl2,lwd=3,hedl=3)
  img2<-readPNG(source=paste(spp2[i],".png",sep=""))
  asp<-dim(img2)[1]/dim(img2)[2]
  rasterImage(img2,xy$x0-w/2,xy$y0-w/2*asp,xy$x0+w/2,xy$y0+w/2*asp)
}

##################################################################


######ALLOMETRY:Figure 5

#Are species' shape differences just a manifestation of shape allometry?
##Allometry

allometry<-procD.lm(Shape~CS,data=TW_gdf)
allometry$aov.table

par(mfrow=c(1,1))
all.plot<-plot(allometry, type = "regression", predictor = TW_phy.csv$Cs,reg.type = "RegScore",pch = 21, bg=TW_col.gp2, xlab ="log[Cs]",cex=2)
all.plot$RegScore
#text(all.plot$RegScore[,1],all.plot$GM$fitted[,,1:60],labels=label)
##text(allometry$X[,-1],allometry$Y,labels=label,pos=1)


########procD.lm####################

TW_phy.csv #paramenters previously imported
TW_gdf #geomorph.dataframe previously created

###procD.lm variables

#Biosonar
Shape.Bios.lm<-procD.lm(Shape ~ Biosonar,data=TW_gdf)
Shape.Bios.lm$aov.table
Cs.Bios.lm<-procD.lm(CS ~ Biosonar,data=TW_gdf)
Cs.Bios.lm$aov.table

#Diet
Shape.Diet.lm<-procD.lm(Shape ~ Diet,data=TW_gdf)
Shape.Diet.lm$aov.table
Cs.Diet.lm<-procD.lm(CS ~ Diet,data=TW_gdf)
Cs.Diet.lm$aov.table

#Diving ecology
Shape.Ecology.lm<-procD.lm(Shape ~ Ecology,data=TW_gdf)
Shape.Ecology.lm$aov.table
Cs.Ecology.lm<-procD.lm(CS ~ Ecology,data=TW_gdf)
Cs.Ecology.lm$aov.table

#SST
Shape.SST.lm<-procD.lm(Shape ~ SST,data=TW_gdf)
Shape.SST.lm$aov.table
Cs.SST.lm<-procD.lm(CS ~ SST, data=TW_gdf)
Cs.SST.lm$aov.table

#Length
Shape.L.lm<-procD.lm(Shape ~ L,data=TW_gdf)
Shape.L.lm$aov.table
Cs.L.lm<-procD.lm(CS ~ L,data=TW_gdf)
Cs.L.lm$aov.table

#BodyMass
Shape.BodyM.lm<-procD.lm(Shape ~ BM,data=TW_gdf)
Shape.BodyM.lm$aov.table
Cs.BodyM.lm<-procD.lm(CS ~ BM,data=TW_gdf)
Cs.BodyM.lm$aov.table

BodyM.Diet.lm<-procD.lm(BM ~ Diet,data=TW_gdf)
BodyM.Diet.lm$aov.table

###############procD.lm with residuals
allometry.res<-allometry$residuals

#procD.lm with residuals
Cs.Bios.res<-procD.lm(allometry.res ~ Biosonar,data=TW_gdf)
Cs.Bios.res$aov.table

Cs.Diet.res<-procD.lm(allometry.res ~ Diet,data=TW_gdf)
Cs.Diet.res$aov.table

Cs.Ecology.res<-procD.lm(allometry.res ~ Ecology,data=TW_gdf)
Cs.Ecology.res$aov.table

Cs.SST.res<-procD.lm(allometry.res ~ SST, data=TW_gdf)
Cs.SST.res$aov.table

Cs.L.lm.res<-procD.lm(allometry.res ~ L,data=TW_gdf)
Cs.L.lm.res$aov.table

Cs.BM.lm.res<-procD.lm(allometry.res ~ BM,data=TW_gdf)
Cs.BM.lm.res$aov.table

#################PGLS############

#evolutionary allometry
allometry.phy<-procD.pgls(f1=Shape ~ CS, phy=tree.matched,data=TW_gdf)
summary(allometry.phy)


#####procD.pgls variables

#Biosonar
Shape.Bios.pgls<-procD.pgls(Shape ~ Biosonar,tree.matched,data= TW_gdf,iter=999)
Shape.Bios.pgls$aov.table
Cs.Bios.pgls<-procD.pgls(CS ~ Biosonar,tree.matched,data= TW_gdf,iter=999)
Cs.Bios.pgls$aov.table

#Diet
Shape.Diet.pgls<-procD.pgls(Shape ~ Diet,phy=tree.matched,data=TW_gdf)
Shape.Diet.pgls$aov.table
Cs.Diet.pgls<-procD.pgls(CS ~ Diet,phy=tree.matched,data=TW_gdf)
Cs.Diet.pgls$aov.table

#Ecology
Shape.Ecology.pgls<-procD.pgls(Shape ~ Ecology,phy=tree.matched,data=TW_gdf)
Shape.Ecology.pgls$aov.table
Cs.Ecology.pgls<-procD.pgls(CS ~ Ecology,phy=tree.matched,data=TW_gdf)
Cs.Ecology.pgls$aov.table

#SST
Shape.SST.pgls<-procD.pgls(Shape ~ SST,phy=tree.matched,data=TW_gdf)
Shape.SST.pgls$aov.table
Cs.SST.pgls<-procD.pgls(CS ~ SST, tree.matched,data=TW_gdf)
Cs.SST.pgls$aov.table

#L
Shape.L.pgls<-procD.pgls(Shape ~ L,phy=tree.matched,data=TW_gdf)
Shape.L.pgls$aov.table
Cs.L.pgls<-procD.pgls(CS ~ L,phy=tree.matched,data=TW_gdf)
Cs.L.pgls$aov.table

#BodyMass
Shape.BodyM.pgls<-procD.pgls(Shape ~ BM,phy=tree.matched,data=TW_gdf)
Shape.BodyM.pgls$aov.table
Cs.BodyM.pgls<-procD.pgls(CS ~ BM,tree.matched,data=TW_gdf)
Cs.BodyM.pgls$aov.table

BodyM.Diet.pgls<-procD.pgls(BM ~ Diet,tree.matched,data=TW_gdf)
BodyM.Diet.pgls$aov.table



#######procD.pgls with residuals 
allometry.phy.res<-allometry.phy$pgls.residuals


#procD.pgls with residuals

Bios.pgls.res<-procD.pgls(allometry.phy.res~ Biosonar,phy=tree.matched,data=TW_gdf)
Bios.pgls.res$aov.table

Diet.pgls.res<-procD.pgls(allometry.phy.res~ Diet,phy=tree.matched,data=TW_gdf)
Diet.pgls.res$aov.table

Ecology.pgls.res<-procD.pgls(allometry.phy.res ~ Ecology,phy=tree.matched,data=TW_gdf)
Ecology.pgls.res$aov.table

SST.pgls.res<-procD.pgls(allometry.phy.res~ SST,phy=tree.matched,data=TW_gdf)
SST.pgls.res$aov.table

L.pgls.res<-procD.pgls(allometry.phy.res ~ L,phy=tree.matched,data=TW_gdf)
L.pgls.res$aov.table

BodyM.pgls.res<-procD.pgls(allometry.phy.res ~ BM,phy=tree.matched,data=TW_gdf)
BodyM.pgls.res$aov.table


#############################################################
####################56 species########

#use na.omit to create the b vector
Dt_56spp<-na.omit(TW_phy.csv$Prey.Mean)
b56spp<-c(1,46,47,48)
Dt_56spp<-TW_phy.csv[-c(b56spp),]
head(Dt_56spp)

TW_ms_56spp<-TW_ms[,,-b56spp]
dimnames(TW_ms_56spp)


Myclade56spp<-c("Delphinus delphis" ,  "Delphinus capensis", "Tursiops truncatus",  "Tursiops aduncus",  "Stenella coeruleoalba" ,  "Stenella frontalis",  
                "Stenella longirostris" ,  "Lagenodelphis hosei" ,  "Stenella attenuata" ,  "Sousa teuszii" ,  "Sousa chinensis" ,  "Sousa plumbea" , 
                "Sotalia fluviatilis" ,  "Steno bredanensis" ,  "Feresa attenuata" , "Globicephala melas",  "Globicephala macrorhynchus" , "Peponocephala electra" , 
                "Pseudorca crassidens" , "Grampus griseus" ,  "Orcaella brevirostris", "Orcaella heinsohni", "Lissodelphis peronii",  "Lissodelphis borealis" , 
                "Lagenorhynchus obscurus" ,  "Lagenorhynchus obliquidens",  "Lagenorhynchus cruciger", "Lagenorhynchus australis", "Cephalorhynchus heavisidii" ,
                "Cephalorhynchus eutropia",  "Cephalorhynchus commersoni",  "Lagenorhynchus albirostris",  "Lagenorhynchus acutus", "Orcinus orca", "Phocoena sinus" , 
                "Phocoena spinipinnis",  "Phocoena dioptrica",  "Phocoena phocoena" ,  "Phocoenoides dalli" , "Neophocaena phocaenoides", "Delphinapterus leucas", 
                "Monodon monoceros" ,  "Inia geoffrensis" , "Pontoporia blanvillei" ,  "Mesoplodon europaeus",  "Mesoplodon mirus",  "Mesoplodon hectori" , 
                "Mesoplodon bidens" , "Indopacetus pacificus", "Hyperodont planifrons" ,  "Hyperodont ampullatus" ,  "Ziphius cavirostris" ,  "Berardius arnuxii",
                "Platanista gangetica",  "Kogia breviceps", "Kogia sima")


seq(1:Ntip(phy))[-which(phy$tip.label%in%Myclade56spp)]->tip.to_56spp


phy_56spp<-drop.tip(phy, tip.to_56spp)
phy_56spp$tip.label#56

Dt_56spp[match(rownames(Dt_56spp),phy_56spp$tip.label),]

Dt_56spp_gdf<-geomorph.data.frame(Shape=TW_ms_56spp,CS=Dt_56spp$Cs,
                                  Preymean=Dt_56spp$Prey.Mean,Preymax=Dt_56spp$Prey.Max,Preymin=Dt_56spp$Prey.min)


#######56species procD.lm

#Prey mean

Dt_56spp.Shape.Preymean.lm<-procD.lm(Shape ~ Preymean,data=Dt_56spp_gdf)
Dt_56spp.Shape.Preymean.lm$aov.table
Dt_56spp.Cs.Preymean.lm<-procD.lm(CS ~ Preymean,data=Dt_56spp_gdf)
Dt_56spp.Cs.Preymean.lm$aov.table

#Prey max

Dt_56spp.Shape.Preymax.lm<-procD.lm(Shape ~ Preymax,data=Dt_56spp_gdf)
Dt_56spp.Shape.Preymax.lm$aov.table
Dt_56spp.Cs.Preymax.lm<-procD.lm(CS ~ Preymax,data=Dt_56spp_gdf)
Dt_56spp.Cs.Preymax.lm$aov.table

#Prey min

Dt_56spp.Shape.Preymin.lm<-procD.lm(Shape ~ Preymin,data=Dt_56spp_gdf)
Dt_56spp.Shape.Preymin.lm$aov.table
Dt_56spp.Cs.Preymin.lm<-procD.lm(CS ~ Preymin,data=Dt_56spp_gdf)
Dt_56spp.Cs.Preymin.lm$aov.table

###############procD.lm with residuals####56 spp
##Allometry 56 spp

allometry_56spp<-procD.lm(Shape~CS,data=Dt_56spp_gdf)
allometry_56spp$aov.table

par(mfrow=c(1,1))
all.plot_56spp<-plot(allometry_56spp, type = "regression", predictor = Dt_56spp_gdf$CS,reg.type = "RegScore",pch = 21, xlab ="log[Cs]",cex=2)
all.plot_56spp$RegScore


allometry.res_56spp<-allometry_56spp$residuals

#procD.lm with residuals
Cs.Preymean.res_56spp<-procD.lm(allometry.res_56spp ~ Preymean,data=Dt_56spp_gdf)
Cs.Preymean.res_56spp$aov.table

Cs.Preymax.res_56spp<-procD.lm(allometry.res_56spp ~ Preymax,data=Dt_56spp_gdf)
Cs.Preymax.res_56spp$aov.table

Cs.Preymin.res_56spp<-procD.lm(allometry.res_56spp ~ Preymin,data=Dt_56spp_gdf)
Cs.Preymin.res_56spp$aov.table

#######56species procD.pgls


#Prey mean

Dt_56spp.Shape.Preymean.pgls<-procD.pgls(Shape ~ Preymean,phy=phy_56spp,data=Dt_56spp_gdf)
Dt_56spp.Shape.Preymean.pgls$aov.table
Dt_56spp.Cs.Preymean.pgls<-procD.pgls(CS ~ Preymean,phy=phy_56spp,data=Dt_56spp_gdf)
Dt_56spp.Cs.Preymean.pgls$aov.table

#Prey max

Dt_56spp.Shape.Preymax.pgls<-procD.pgls(Shape ~ Preymax,phy=phy_56spp,data=Dt_56spp_gdf)
Dt_56spp.Shape.Preymax.pgls$aov.table
Dt_56spp.Cs.Preymax.pgls<-procD.pgls(CS ~ Preymax,phy=phy_56spp,data=Dt_56spp_gdf)
Dt_56spp.Cs.Preymax.pgls$aov.table


#Prey min

Dt_56spp.Shape.Preymin.pgls<-procD.pgls(Shape ~ Preymin,phy=phy_56spp,data=Dt_56spp_gdf)
Dt_56spp.Shape.Preymin.pgls$aov.table
Dt_56spp.Cs.Preymin.pgls<-procD.pgls(CS ~ Preymin,phy=phy_56spp,data=Dt_56spp_gdf)
Dt_56spp.Cs.Preymin.pgls$aov.table


###############procD.pgls with residuals####56 spp
##Allometry 56 spp

allometry.phy_56spp<-procD.pgls(Shape~CS,phy=phy_56spp,data=Dt_56spp_gdf)
allometry.phy_56spp$aov.table


allometry.phy_56spp.res<-allometry.phy_56spp$pgls.residuals

#procD.pgls with residuals
Cs.Preymean.pgls.res_56spp<-procD.pgls(allometry.phy_56spp.res ~ Preymean,phy=phy_56spp,data=Dt_56spp_gdf)
Cs.Preymean.pgls.res_56spp$aov.table

Cs.Preymax.pgls.res_56spp<-procD.pgls(allometry.phy_56spp.res ~ Preymax,phy=phy_56spp,data=Dt_56spp_gdf)
Cs.Preymax.pgls.res_56spp$aov.table

Cs.Preymin.pgls.res_56spp<-procD.pgls(allometry.phy_56spp.res ~ Preymin,phy=phy_56spp,data=Dt_56spp_gdf)
Cs.Preymin.pgls.res_56spp$aov.table



#############################################
############################31 species#########################


##31 species procD.lm
######EQ
#use na.omit to create the b vector
b31spp<-c(1,3,5,7,9:13,16,19,22:24,28:31,34,36,38,47,48,51:55,57)
Dt_31spp<-TW_phy.csv[-c(b31spp),]



TW_ms_31spp<-TW_ms[,,-b31spp]
dimnames(TW_ms_31spp)



Myclade31spp<-c("Delphinus delphis",  "Tursiops truncatus" ,  "Stenella coeruleoalba" ,  "Stenella longirostris"  ,  "Sotalia fluviatilis" ,  "Steno bredanensis", 
                "Globicephala melas",  "Globicephala macrorhynchus",  "Pseudorca crassidens",  "Grampus griseus"  , "Lissodelphis borealis", "Lagenorhynchus obscurus",
                "Lagenorhynchus obliquidens" , "Cephalorhynchus commersoni",  "Lagenorhynchus albirostris" , "Orcinus orca" ,  "Phocoena spinipinnis", "Phocoena phocoena",
                "Phocoenoides dalli" ,  "Neophocaena phocaenoides"  ,  "Delphinapterus leucas"  ,  "Monodon monoceros",  "Inia geoffrensis" ,  "Pontoporia blanvillei"  , 
                "Lipotes vexillifer", "Mesoplodon europaeus",  "Mesoplodon mirus" ,  "Ziphius cavirostris", "Platanista gangetica", "Kogia breviceps"  , "Kogia sima")

seq(1:Ntip(phy))[-which(phy$tip.label%in%Myclade31spp)]->tip.to_31spp

phy_31spp<-drop.tip(phy, tip.to_31spp)
phy_31spp$tip.label

Dt_31spp[match(rownames(Dt_31spp),phy_31spp$tip.label),]

Dt_31spp_gdf<-geomorph.data.frame(Shape=TW_ms_31spp,CS=Dt_31spp$Cs,
                                   EQ=Dt_31spp$Lo10.EQ.,BrainM=Dt_31spp$Log10.Brain.mass..g..)

#EQ
Dt_31spp.Shape.EQ.lm<-procD.lm(Shape ~ EQ,data=Dt_31spp_gdf)
Dt_31spp.Shape.EQ.lm$aov.table
Dt_31spp.Cs.EQ.lm<-procD.lm(CS ~ EQ,data=Dt_31spp_gdf)
Dt_31spp.Cs.EQ.lm$aov.table

#Brain Mass
Dt_31spp.Shape.BrainM.lm<-procD.lm(Shape ~ BrainM,data=Dt_31spp_gdf)
Dt_31spp.Shape.BrainM.lm$aov.table
Dt_31spp.CS.BrainM.lm<-procD.lm(CS ~ BrainM,data=Dt_31spp_gdf)
Dt_31spp.CS.BrainM.lm$aov.table

###############procD.lm with residuals####31 spp
##Allometry 31 spp

allometry_31spp<-procD.lm(Shape~CS,data=Dt_31spp_gdf)
allometry_31spp$aov.table

par(mfrow=c(1,1))
all.plot_31spp<-plot(allometry_31spp, type = "regression", predictor = Dt_31spp_gdf$CS,reg.type = "RegScore",pch = 21, xlab ="log[Cs]",cex=2)
all.plot_31spp$RegScore


allometry.res_31spp<-allometry_31spp$residuals

#procD.lm with residuals
Cs.EQ.res_31spp<-procD.lm(allometry.res_31spp ~ EQ,data=Dt_31spp_gdf)
Cs.EQ.res_31spp$aov.table

Cs.BrainM.res_31spp<-procD.lm(allometry.res_31spp ~ BrainM,data=Dt_31spp_gdf)
Cs.BrainM.res_31spp$aov.table



#########31 species procD.pgls

#EQ
Dt_31spp.Shape.EQ.pgls<-procD.pgls(Shape ~ EQ,phy=phy_31spp,data=Dt_31spp_gdf)
Dt_31spp.Shape.EQ.pgls$aov.table
Dt_31spp.Cs.EQ.pgls<-procD.pgls(CS ~ EQ,phy=phy_31spp,data=Dt_31spp_gdf)
Dt_31spp.Cs.EQ.pgls$aov.table

#Brain Mass

Dt_31spp.Shape.BrainM.pgls<-procD.pgls(Shape ~ BrainM,phy=phy_31spp,data=Dt_31spp_gdf)
Dt_31spp.Shape.BrainM.pgls$aov.table
Dt_31spp.Cs.BrainM.pgls<-procD.pgls(CS ~ BrainM,phy=phy_31spp,data=Dt_31spp_gdf)
Dt_31spp.Cs.BrainM.pgls$aov.table

###############procD.pgls with residuals####31 spp
##Allometry 31 spp

allometry.phy_31spp<-procD.pgls(Shape~CS,phy=phy_31spp,data=Dt_31spp_gdf)
allometry.phy_31spp$aov.table


allometry.phy_31spp.res<-allometry.phy_31spp$pgls.residuals

#procD.pgls with residuals
Cs.EQ.pgls.res_31spp<-procD.pgls(allometry.phy_31spp.res ~ EQ,phy=phy_31spp,data=Dt_31spp_gdf)
Cs.EQ.pgls.res_31spp$aov.table

Cs.BrainM.pgls.res_31spp<-procD.pgls(allometry.phy_31spp.res ~ BrainM,phy=phy_31spp,data=Dt_31spp_gdf)
Cs.BrainM.pgls.res_31spp$aov.table




###############################################
#############################26 species#########################


#na.omit(TW_phy.csv$KhZmin)
b_kHz<-c(1,  2 , 3,  4,  6 , 7,  8,  9 ,10, 11 ,13, 14, 15, 18, 19, 23, 24,
         25, 27, 31, 33, 34, 36, 37, 38, 45, 46, 47 ,48, 50, 51, 54, 57, 60) 


Dt_26spp<-TW_phy.csv[-(b_kHz),]

TW_ms_26spp<-TW_ms[,,-b_kHz]
dimnames(TW_ms_26spp)

Myclade26spp<-c("Tursiops aduncus" , "Sousa chinensis", "Feresa attenuata", 
                "Globicephala melas", "Pseudorca crassidens",  "Grampus griseus",
                "Orcaella brevirostris" , "Lagenorhynchus obscurus" , "Lagenorhynchus cruciger", "Lagenorhynchus australis",
                "Cephalorhynchus heavisidii", "Cephalorhynchus commersoni", "Orcinus orca" , "Phocoena phocoena", "Phocoenoides dalli",
                "Neophocaena phocaenoides"  , "Delphinapterus leucas" , "Monodon monoceros"   , "Inia geoffrensis", "Mesoplodon europaeus",
                "Mesoplodon bidens", "Indopacetus pacificus", "Hyperodont ampullatus", "Ziphius cavirostris" , "Platanista gangetica","Kogia breviceps")

seq(1:Ntip(phy))[-which(phy$tip.label%in%Myclade26spp)]->tip.to_26spp

phy_26spp<-drop.tip(phy, tip.to_26spp)
phy_26spp$tip.label

Dt_26spp[match(rownames(Dt_26spp),phy_26spp$tip.label),]

Dt_26spp_gdf<-geomorph.data.frame(Shape=TW_ms_26spp,CS=Dt_26spp$Cs,
                                  kHzmin=Dt_26spp$KhZmin,kHzmax=Dt_26spp$KhZmax,
                                  Biosonar=Dt_26spp$Biosonar..Jensen.et.al.2018.)

#kHzmin
Dt_26spp.Shape.kHzmin.lm<-procD.lm(Shape ~ kHzmin,data=Dt_26spp_gdf)
Dt_26spp.Shape.kHzmin.lm$aov.table
Dt_26spp.Cs.kHzmin.lm<-procD.lm(CS ~ kHzmin,data=Dt_26spp_gdf)
Dt_26spp.Cs.kHzmin.lm$aov.table

#kHzmax
Dt_26spp.Shape.kHzmax.lm<-procD.lm(Shape ~ kHzmax,data=Dt_26spp_gdf)
Dt_26spp.Shape.kHzmax.lm$aov.table
Dt_26spp.CS.kHzmax.lm<-procD.lm(CS ~ kHzmax,data=Dt_26spp_gdf)
Dt_26spp.CS.kHzmax.lm$aov.table

#Biosonar
Dt_26spp.Shape.Biosonar.lm<-procD.lm(Shape ~ Biosonar,data=Dt_26spp_gdf)
Dt_26spp.Shape.Biosonar.lm$aov.table
Dt_26spp.Cs.Biosonar.lm<-procD.lm(CS ~ Biosonar,data=Dt_26spp_gdf)
Dt_26spp.Cs.Biosonar.lm$aov.table


###############procD.lm with residuals####26 spp
##Allometry 26 spp

allometry_26spp<-procD.lm(Shape~CS,data=Dt_26spp_gdf)
allometry_26spp$aov.table

par(mfrow=c(1,1))
all.plot_26spp<-plot(allometry_26spp, type = "regression", predictor = Dt_26spp_gdf$CS,reg.type = "RegScore",pch = 21, xlab ="log[Cs]",cex=2)
all.plot_26spp$RegScore


allometry.res_26spp<-allometry_26spp$residuals

#procD.lm with residuals
Cs.kHzmin.res_26spp<-procD.lm(allometry.res_26spp ~ kHzmin,data=Dt_26spp_gdf)
Cs.kHzmin.res_26spp$aov.table

Cs.kHzmax.res_26spp<-procD.lm(allometry.res_26spp ~ kHzmax,data=Dt_26spp_gdf)
Cs.kHzmax.res_26spp$aov.table

Cs.Biosonar.res_26spp<-procD.lm(allometry.res_26spp ~ Biosonar,data=Dt_26spp_gdf)
Cs.Biosonar.res_26spp$aov.table

#########26 species procD.pgls

#kHzmin
Dt_26spp.Shape.kHzmin.pgls<-procD.pgls(Shape ~ kHzmin,phy=phy_26spp,data=Dt_26spp_gdf)
Dt_26spp.Shape.kHzmin.pgls$aov.table
Dt_26spp.Cs.kHzmin.pgls<-procD.pgls(CS ~ kHzmin,phy=phy_26spp,data=Dt_26spp_gdf)
Dt_26spp.Cs.kHzmin.pgls$aov.table

#kHzmax

Dt_26spp.Shape.kHzmax.pgls<-procD.pgls(Shape ~ kHzmax,phy=phy_26spp,data=Dt_26spp_gdf)
Dt_26spp.Shape.kHzmax.pgls$aov.table
Dt_26spp.Cs.kHzmax.pgls<-procD.pgls(CS ~ kHzmax,phy=phy_26spp,data=Dt_26spp_gdf)
Dt_26spp.Cs.kHzmax.pgls$aov.table

#Biosonar
Dt_26spp.Shape.Biosonar.pgls<-procD.pgls(Shape ~ Biosonar,phy=phy_26spp,data=Dt_26spp_gdf)
Dt_26spp.Shape.Biosonar.pgls$aov.table
Dt_26spp.Cs.Biosonar.pgls<-procD.pgls(CS ~ Biosonar,phy=phy_26spp,data=Dt_26spp_gdf)
Dt_26spp.Cs.Biosonar.pgls$aov.table

###############procD.pgls with residuals####26 spp
##Allometry 26 spp

allometry.phy_26spp<-procD.pgls(Shape~CS,phy=phy_26spp,data=Dt_26spp_gdf)
allometry.phy_26spp$aov.table


allometry.phy_26spp.res<-allometry.phy_26spp$pgls.residuals

#procD.pgls with residuals
Cs.kHzmin.pgls.res_26spp<-procD.pgls(allometry.phy_26spp.res ~ kHzmin,phy=phy_26spp,data=Dt_26spp_gdf)
Cs.kHzmin.pgls.res_26spp$aov.table

Cs.kHzmax.pgls.res_26spp<-procD.pgls(allometry.phy_26spp.res ~ kHzmax,phy=phy_26spp,data=Dt_26spp_gdf)
Cs.kHzmax.pgls.res_26spp$aov.table

Cs.Biosonar.pgls.res_26spp<-procD.pgls(allometry.phy_26spp.res ~ Biosonar,phy=phy_26spp,data=Dt_26spp_gdf)
Cs.Biosonar.pgls.res_26spp$aov.table



######Plot 4 phy together

par(mfrow=c(1,4))

plot(phy,cex=0.5)
plot(phy_56spp,cex=0.5)
plot(phy_31spp,cex=0.5)
plot(phy_26spp,cex=0.5)


save.image("/Online Resource 1/Macro_script.RData")
