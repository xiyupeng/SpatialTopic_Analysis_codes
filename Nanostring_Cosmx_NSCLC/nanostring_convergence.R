### Use Lung5-1 sample of nanostring NSCLC dataset as benchmarking

## We use Seurat v5 package to visualize the results.
## If you still use Seurat v4, you will have the error
library(Seurat, quietly = TRUE);packageVersion("Seurat")
## Load the Seurat object for the image
load("~/github/SpaTopic_test/data/nanostring_example.rdata")
## for large dataset
options(future.globals.maxSize = 1e9)

library(ggplot2)
library(SpaTopic)
library(sf)

## Prepare input from Seurat Object
dataset<-Seurat5obj_to_SpaTopic(object = nano.obj, group.by = "predicted.annotation.l1",image = "image1")
head(dataset)

# Gibbs sampling for SpaTopic (trace the Gibbs Sampling line)
system.time(gibbs.res<-SpaTopic_inference(dataset, ntopics = 7, sigma = 50, 
                                          region_radius = 400,burnin = 0,trace = TRUE,thin = 20, niter = 200))

### plot the log likelihood of the method
trace_result<-data.frame(loglike = gibbs.res$loglike.trace, iteration = seq(0, 4000, 20))
ggplot(trace_result, aes(x = iteration, y = loglike))+ geom_line()+theme_minimal()


ggsave("~/flow_cytometry/spatial_topic/manuscript/images/convergence.png",width = 6,height = 3, units = "in")
