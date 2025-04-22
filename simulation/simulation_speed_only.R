
#####################################################################################################
################# 
## SIMULATION from the figure

library(dplyr)

load(file = "~/Documents/Research/github/SpaTopic_benchmarking/simulation/rerun/data/LIBDsubsample.RData")
### subsample

ARI_lda<-list()
ARI_knn<-list()
fig_data<-list()
fig_true<-list()
fig_lda<-list()
fig_knn<-list()
time_lda<-list()
time_knn<-list()


for(i in 1:5){

  set.seed(NULL)
    
  n.celltypes<-20
  perc<-matrix(0,nrow = 20, ncol = 4)
  perc[1:5,1]<-rgamma(5,1,1)
  perc[6:10,2]<-rgamma(5,1,1)
  perc[11:15,3]<-rgamma(5,1,1)
  perc[16:20,4]<-rgamma(5,1,1)
  perc<-t(perc)/apply(perc,2,sum)
  perc<-t(perc)

  subsample$niche<-subsample$label

  simulation_MIF<-function(topic,dist,n.cells,n.celltypes){
    
    C<-NULL
    ### multinomial distribution 
    for (i in 1:n.cells){
      C[i]<- sample(n.celltypes, 1, FALSE, dist[,topic[i]])
    }  
    
    return(C)
  }

  n.celltypes<-nrow(perc)

  set.seed(123)
  subsample$type<-simulation_MIF(subsample$niche,perc,nrow(subsample),n.celltypes)


  ### niches
  fig_true[[i]]<-ggplot(subsample,aes(x = x_coor,y = y_coor, color = as.factor(label)))+geom_point(size = 0.8)

  ### celltype
  fig_data[[i]]<-ggplot(subsample,aes(x = x_coor,y = y_coor, color = as.factor(type)))+geom_point(size = 0.8)


  #### running results with SpaTopic----------------------------------------------------------

  dataset_simu<-subsample %>% 
    mutate(X = x_coor, Y = y_coor, image = 1,type = type) %>%
    select(image,X,Y,type)

  dataset_simu<-list(dataset_simu) 

  library(SpaTopic)
  time_lda[[i]]<-system.time(gibbs.res<-SpaTopic_inference(dataset_simu, sigma = 50, ntopics = 4, region_radius = 60))

  #### visualziation
  prob<-as.matrix(gibbs.res$Z.trace)
  topic<-as.factor(apply(prob,1,which.max))
  subsample$topic<-topic
  fig_lda[[i]]<-ggplot(subsample,aes(x = x_coor,y = y_coor, color = as.factor(topic)))+geom_point(size = 0.8)

  ARI_lda[[i]]<-adjustedRandIndex(subsample$topic,subsample$niche)
  #fig_lda[[i]]
  #ARI_lda[[i]]

  #####running result with KNN-kmeans in Seurat v5----------------------------------------------------------------------


  ### Create Seurat object
  library(Seurat)
  seurat.obj<-CreateSeuratObject(counts = matrix(0,nrow = 10, ncol = nrow(subsample)))
  seurat.obj<-AddMetaData(seurat.obj,metadata = subsample)

  coords <- CreateFOV(
    data.frame(X = seurat.obj$x_coor,Y = seurat.obj$y_coor),
    type = c("centroids"),assay = "RNA"
  )
  seurat.obj[["simu"]] <-coords

  time_knn[[i]]<-system.time(nano.obj <- BuildNicheAssay(object = seurat.obj, fov = "simu", group.by = "type",
                                          niches.k = 4))


  fig_knn[[i]]<-ggplot(subsample,aes(x = x_coor,y = y_coor, color = as.factor(nano.obj$niches)))+geom_point(size = 0.8)
  ARI_knn[[i]]<-adjustedRandIndex(nano.obj$niches,subsample$niche)

  #fig_knn[[i]]
  #ARI_knn[[i]]
  gc()

}

saveRDS(list(ARI_lda = ARI_lda,ARI_knn = ARI_knn,fig_data = fig_data,fig_true = fig_true, fig_lda = fig_lda,
             fig_knn = fig_knn,time_knn= time_knn,time_lda = time_lda),
        file = "~/Documents/Research/github/SpaTopic_benchmarking/simulation/rerun/Simulation_figure_s50_r60_nint10_times_new.rds")

library(ggpubr)
result<-as.data.frame(cbind(ARI_lda,ARI_knn))
result %>% tidyr::gather("method","ARI") %>% mutate(ARI = unlist(ARI)) %>%
  ggboxplot(., "method","ARI")


