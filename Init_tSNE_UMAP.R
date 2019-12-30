Path <- "/home/nikolay/WABI/K_Pietras/Easy_scRNAseq_tSNE_Cluster/Ultra_Large_DataSets/"
setwd(Path)
files<-list.files(pattern="*.txt")

library("Rtsne")
library("dimRed")
library("umap")
library("pals")

for(k in 3:3) # 2:length(files)
{

#READ DATA, FILTER AND SUMMARIZE IT
file <- files[k]
print(paste0("START WORKING WITH FILE ", file))
name_cap <- toupper(matrix(unlist(strsplit(file,"\\.")),ncol=2,byrow=TRUE)[1])
name <- matrix(unlist(strsplit(file,"\\.")),ncol=2,byrow=TRUE)[1]
pdf(paste0("Easy_Init/PLOTS_",name,".pdf"), paper="a4r", width = 210, height = 297)

expr<-read.delim(file,header=TRUE,row.names=1,check.names=FALSE,sep="\t")
expr<-expr[rowMeans(as.matrix(expr))>=1,]
N_cells<-dim(expr)[2]
print(expr[1:5,1:5])
print(paste0("DATASET CONTAINS ",dim(expr)[1]," GENES AND ",N_cells," CELLS"))
print(paste0("PERCENT OF MISSING VALUES = ",round(100*sum(expr==0)/(dim(expr)[1]*dim(expr)[2]),0)," %"))

#READ CLUSTER ASSIGNMENTS FOR ALL CELLS
cluster<-read.delim(paste0("Easy_Output/CLUSTER_ASSIGNMENT_",file),header=TRUE,sep="\t")
print(head(cluster))
N_clust<-length(sort(unique(cluster$CLUSTER)))
if(N_clust <= 25){colors <- cols25(N_clust)}else{colors <- rainbow(N_clust)}
names(colors) <- sort(unique(cluster$CLUSTER))
my_col <- colors[as.character(cluster$CLUSTER)]

#READ CLACULATED PREVIOUSLY OPTIMAL PERPLEXITY AND NUMBER OF SIGNIFICANT PRINCIPAL COMPONENTS
log_perp<-system(paste0("grep 'OPTIMAL PERPLEXITY = ' ",
                        Path,paste0("Easy_Output/LOG_",name,".txt")),intern=TRUE)
optPerp<-as.numeric(gsub("\"","",matrix(unlist(strsplit(log_perp," = ")),ncol=2,byrow=TRUE)[,2]))
print(paste0("OPTIMAL PERPLEXITY = ",optPerp))

log_pc<-system(paste0("grep 'PRINCIPAL COMPONENTS = ' ",
                      Path,paste0("Easy_Output/LOG_",name,".txt")),intern=TRUE)
optPC<-as.numeric(matrix(unlist(strsplit(gsub("\"","",matrix(unlist(strsplit(log_pc," = ")),
                                           ncol=2,byrow=TRUE)[,2]),",")),ncol=2,byrow=TRUE)[1])
print(paste0("OPTIMAL NUMBER OF PRINCIPAL COMPONENTS = ",optPC))

#COMPUTE PCA
PC <- prcomp(t(log10(expr + 1)), center=TRUE, scale=FALSE)
plot(PC$x[,1:2], col=my_col,main=paste0("PCA PLOT: ",name_cap),pch=19)

#COMPUTE LAPLACIAN EIGENMAPS
leim <- LaplacianEigenmaps()
emb <- leim@fun(as(t(log10(expr + 1)),"dimRedData"), leim@stdpars)
plot(emb@data@data,col=my_col,main=paste0("LAPLACIAN EIGENMAPS PLOT: ",name_cap),pch=19)


############### ###################### tSNE ###############################################################

print("COMPUTING tSNE PLOTS")

#INITIALIZATION RANDOM
tsne_opt_perp <- Rtsne(t(log10(expr+1)),initial_dims=optPC,verbose=FALSE,check_duplicates=FALSE,
                       perplexity=optPerp,dims=2,max_iter=10000)
plot(tsne_opt_perp$Y,col=my_col,xlab="tSNE1",ylab="tSNE2",cex=0.8,
     main=paste0("tSNE PLOT ",name_cap,": RANDOM INITIALIZATION"),pch=19)

#INITIALIZATION PCA
tsne_init_pca <- Rtsne(t(log10(expr+1)),initial_dims=optPC,verbose=FALSE,check_duplicates=FALSE,
                       perplexity=optPerp,dims=2,max_iter=10000,Y_init=PC$x[,1:2])
#tsne_pca_init <- fftRtsne(t(log10(expr+20)),check_duplicates=FALSE,perplexity=350,verbose=TRUE,
#                                 dims=2,max_iter=10000,initialization=expr[1:2,])
plot(tsne_init_pca$Y,col=my_col,xlab="tSNE1",ylab="tSNE2",cex=0.8,
        main=paste0("tSNE PLOT ",name_cap,": PCA INITIALIZATION"),pch=19)

#INITIALIZATION LAPLACIAN EIGENMAPS
tsne_init_le <- Rtsne(t(log10(expr+1)),initial_dims=optPC,verbose=FALSE,check_duplicates=FALSE,
                   perplexity=optPerp,dims=2,max_iter=10000,Y_init=emb@data@data)
plot(tsne_init_le$Y,col=my_col,xlab="tSNE1",ylab="tSNE2",cex=0.8,
     main=paste0("tSNE PLOT ",name_cap,": LAPLACIAN EIGENMAPS INITIALIZATION"),pch=19)


##################################### UMAP ##############################################################

print("COMPUTING UMAP PLOTS")
custom.settings = umap.defaults
custom.settings$n_neighbors = optPerp
custom.settings$min_dist = 0.1
n_epochs = 10000
method = "umap-learn"
config = custom.settings

#INITIALIZATION RANDOM
umap_init_rand <- umap(PC$x[,1:optPC], config=config, method=method, init="random", n_epochs=n_epochs)
plot(umap_init_rand$layout,col=my_col,xlab="UMAP1",ylab="UMAP2",cex=0.8,
     main=paste0("UMAP PLOT ",name_cap,": RANDOM INITIALIZATION"),pch=19)

#INITIALIZATION PCA
umap_init_pca <- umap(PC$x[,1:optPC], config = config, method=method, init=PC$x[,1:2], n_epochs=n_epochs)
plot(umap_init_pca$layout,col=my_col,xlab="UMAP1",ylab="UMAP2",cex=0.8,
     main=paste0("UMAP PLOT ",name_cap,": PCA INITIALIZATION"),pch=19)

#INITIALIZATION LAPLACIAN EIGNEMAPS
umap_init_le <- umap(PC$x[,1:optPC], config=config, method=method, init=emb@data@data,n_epochs=n_epochs)
plot(umap_init_le$layout,col=my_col,xlab="UMAP1",ylab="UMAP2",cex=0.8,
     main=paste0("UMAP PLOT ",name_cap,": LAPLACIAN EIGENMAPS INITIALIZATION"),pch=19)

print("*********************************************************************************************")
dev.off()

}
        