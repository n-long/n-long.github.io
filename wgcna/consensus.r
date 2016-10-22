library(WGCNA)
library(DESeq2)
library(cluster)
options(stringsAsFactors = FALSE);
setwd(dir = '/home/nlong/Public/clc')
norm=read.csv("clcnorm", sep='\t', header=T)
still=read.csv("clcstill", sep='\t', header=T)
datExprnorm = as.data.frame(norm[, -c(1)]);
datExprstill = as.data.frame(still[, -c(1)]);
rownames(datExprnorm) = norm$gene
rownames(datExprstill) = still$gene
datExprnorm=t(datExprnorm)
datExprstill=t(datExprstill)

datExprnorm<-as.matrix(t(datExprnorm))
datExprstill<-as.matrix(t(datExprstill))
condition<-(c("p1", "p3", "a1", "a3", "p1", "p3", "a1", "a3"))

coldata<-data.frame(row.names=colnames(datExprnorm), condition)
ddsnorm=DESeqDataSetFromMatrix(datExprnorm, colData=coldata, design=~condition)
ddsnorm2=DESeq(ddsnorm,betaPrior=F)
vsdnorm=getVarianceStabilizedData(ddsnorm2)
vsd2norm=t(vsdnorm)

coldatastill<-data.frame(row.names=colnames(datExprstill), condition)
ddsstill=DESeqDataSetFromMatrix(datExprstill, colData=coldatastill, design=~condition)
ddsstill2=DESeq(ddsstill,betaPrior=F)
vsdstill=getVarianceStabilizedData(ddsstill2)
#vsd2still=t(vsdstill)
nSets=2
setLabels= c("norm", "still")
multiExpr=vector(mode="list", length=nSets)
multiExpr[[1]]=list(data=as.data.frame(t(vsdnorm)));
names(multiExpr[[1]]$data)=norm$gene
rownames(multiExpr[[1]]$data)=names(norm)[-c(1)]
multiExpr[[2]]=list(data=as.data.frame(t(vsdstill)));
names(multiExpr[[2]]$data)=still$gene;
rownames(multiExpr[[2]]$data)=names(still)[-c(1)]
exprSize=checkSets(multiExpr)
multiExpr2=t(multiExpr)

gsg = goodSamplesGenesMS(multiExpr, verbose = 3);
gsg$allOK

if (!gsg$allOK)
{
  # Print information about the removed genes:
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes],
                                              collapse = ", ")))
  for (set in 1:exprSize$nSets)
  {
    if (sum(!gsg$goodSamples[[set]]))
      printFlush(paste("In set", setLabels[set], "removing samples",
                       paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
    # Remove the offending genes and samples
    multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
  }
  # Update exprSize
  exprSize = checkSets(multiExpr)
}          

gsg = goodSamplesGenes(vsdnorm, verbose = 3);
gsg$allOK
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(vsdnorm)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(vsdnorm)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  vsdnorm = vsdnorm[gsg$goodSamples, gsg$goodGenes]
}


net = blockwiseConsensusModules(
  multiExpr, power = 18, minModuleSize = 30, deepSplit = 2, 
 maxBlockSize=15000, 
 randomSeed=12345, 
 
 saveIndividualTOMs=TRUE, individualTOMFileNames = "individualTOM-Set%s-Block%b.RData",
  checkPower=TRUE, checkMissingData = TRUE,
  pamRespectsDendro = FALSE, method = "hybrid", networkType = "signed", TOMType = "signed",
  saveConsensusTOMs=TRUE, consensusTOMFileNames= "consensusTOM-block.%b.Rdata",
  mergeCutHeight = 0.25, numericLabels = TRUE,
  minKMEtoStay = 0.2, checkMinModuleSize=TRUE, detectCutHeight=0.8,
  saveTOMs = TRUE, verbose = 5)

softPower = 14;
# Initialize an appropriate array to hold the adjacencies
adjacencies = array(0, dim = c(nSets, nGenes, nGenes));
# Calculate adjacencies in each individual data set
for (set in 1:nSets)
  adjacencies[set, , ] = abs(cor(multiExpr[[set]]$data, use = "p"))^softPower;

####NORM ONLY AUTO BLOCK########
norm2net = blockwiseModules(vsdnorm, power = 10,
                       TOMType = "signed", minModuleSize = 10, deepSplit =2, maxBlockSize=5000,
                       reassignThreshold = 0, mergeCutHeight = 0.15, cutHeight="hybrid", networkType="signed",
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE, minKMEtoStay=0.3,
                       saveTOMFileBase = "normnetTOM",
                       verbose = 3)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "normnet-02-networkConstruction-auto.RData")

# Isolate the module labels in the order they appear in ordered module eigengenes
normModuleLabels = substring(names(normMEs), 3)
consModuleLabels = substring(names(consMEs[[1]]$data), 3)
# Convert the numeric module labels to color labels
normModules = labels2colors(as.numeric(normModuleLabels))
consModules = labels2colors(as.numeric(consModuleLabels))
# Numbers of female and consensus modules
nnormMods = length(normModules)
nConsMods = length(consModules)
# Initialize tables of p-values and of the corresponding counts
pTable = matrix(0, nrow = nnormMods, ncol = nConsMods);
CountTbl = matrix(0, nrow = nnormMods, ncol = nConsMods);
# Execute all pairwaise comparisons
for (fmod in 1:nnormMods)
  for (cmod in 1:nConsMods)
  {
    normMembers = (normaleColors == normModules[fmod]);
    consMembers = (moduleColors == consModules[cmod]);
    pTable[fmod, cmod] = -log10(fisher.test(normMembers, consMembers, alternative = "greater")$p.value);
    CountTbl[fmod, cmod] = sum(normaleColors == normModules[fmod] & moduleColors ==
                                 consModules[cmod])
  }
    
    
    
femModuleLabels = substring(names(femaleMEs), 3)
consModuleLabels = substring(names(consMEs[[1]]$data), 3)
# Convert the numeric module labels to color labels
femModules = labels2colors(as.numeric(femModuleLabels))
consModules = labels2colors(as.numeric(consModuleLabels))
# Numbers of female and consensus modules
nFemMods = length(femModules)
nConsMods = length(consModules)
# Initialize tables of p-values and of the corresponding counts
pTable = matrix(0, nrow = nFemMods, ncol = nConsMods);
CountTbl = matrix(0, nrow = nFemMods, ncol = nConsMods);
# Execute all pairwaise comparisons
pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)]);
pTable[pTable>50 ] = 50 ;
# Marginal counts (really module sizes)
femModTotals = apply(CountTbl, 1, sum)
consModTotals = apply(CountTbl, 2, sum)
# Actual plotting
sizeGrWindow(10,7 );
par(mfrow=c(1,1));
par(cex = 1.0);
par(mar=c(8, 10.4, 2.7, 1)+0.3);
# Use function labeledHeatmap to produce the color-coded table with all the trimmings
labeledHeatmap(Matrix = pTable,
               2xLabels = paste(" ", consModules),
               yLabels = paste(" ", femModules),
               colorLabels = TRUE,
               xSymbols = paste("Cons ", consModules, ": ", consModTotals, sep=""),
               ySymbols = paste("Fem ", femModules, ": ", femModTotals, sep=""),
               textMatrix = CountTbl,
               colors = greenWhiteRed(100)[50:100],
               main = "Correspondence of Female set-specific and Female-Male consensus modules",
               cex.text = 1.0, cex.lab = 1.0, setStdMargins = FALSE);
  

consMEs.unord = multiSetMEs(multiExpr, universalColors = moduleLabels, excludeGrey = TRUE)
#GS = list();
kME = list();
for (set in 1:nSets)
{
  #GS[[set]] = corAndPvalue(multiExpr[[set]]$data, Traits[[set]]$data);
  kME[[set]] = corAndPvalue(multiExpr[[set]]$data, consMEs.unord[[set]]$data);
}


femModuleLabels = substring(names(femaleMEs), 3)
consModuleLabels = substring(names(consMEs[[1]]$data), 3)
# Convert the numeric module labels to color labels
femModules = labels2colors(as.numeric(femModuleLabels))
consModules = labels2colors(as.numeric(consModuleLabels))
# Numbers of female and consensus modules
nFemMods = length(femModules)
nConsMods = length(consModules)
# Initialize tables of p-values and of the corresponding counts
pTable = matrix(0, nrow = nFemMods, ncol = nConsMods);
CountTbl = matrix(0, nrow = nFemMods, ncol = nConsMods);
# Execute all pairwaise comparisons
for (fmod in 1:nFemMods)
  for (cmod in 1:nConsMods)
  {
    femMembers2 = (femaleColors2 == femModules[fmod]);
    consMembers = (moduleColors == consModules[cmod]);
    pTable[fmod, cmod] = -log10(fisher.test(femMembers2, consMembers, alternative = "greater")$p.value);
    CountTbl[fmod, cmod] = sum(femaleColors2 == femModules[fmod] & moduleColors ==
                                 consModules[cmod])
  }






pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)]);
pTable[pTable>50 ] = 50 ;
# Marginal counts (really module sizes)
femModTotals = apply(CountTbl, 1, sum)
consModTotals = apply(CountTbl, 2, sum)
# Actual plotting
sizeGrWindow(10,7 );
par(mfrow=c(1,1));
par(cex = 1.0);
par(mar=c(8, 10.4, 2.7, 1)+0.3);
# Use function labeledHeatmap to produce the color-coded table with all the trimmings
labeledHeatmap(Matrix = pTable,
               xLabels = paste(" ", consModules),
               yLabels = paste(" ", femModules),
               colorLabels = TRUE,
               xSymbols = paste("Cons ", consModules, ": ", consModTotals, sep=""),
               ySymbols = paste("Fem ", femModules, ": ", femModTotals, sep=""),
               textMatrix = CountTbl,
               colors = greenWhiteRed(100)[50:100],
               main = "Correspondence of Female set-specific and Female-Male consensus modules",
               cex.text = 1.0, cex.lab = 1.0, setStdMargins = FALSE);

gene.names=colnames(multiExpr[[1]]$data)

module_colors= setdiff(unique(consModules), "grey")
for (color in module_colors){
  module=gene.names[which(consModules==color)]
  write.table(module, paste("cons_",color, ".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
}
