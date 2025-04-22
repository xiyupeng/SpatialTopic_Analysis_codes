library("grid")
library("gridExtra")

## modified from https://github.com/shangll123/SpatialPCA_analysis_codes with original author Lulu Shang

library(jpeg)
mandrill = readJPEG("/Users/xiyupeng/Documents/Research/github/SpaTopic_benchmarking/simulation/LIBDsimu.jpg", native = FALSE)
# > dim(mandrill)
# [1] 1100  984    3
# copy the image three times
mandrill.R = mandrill
mandrill.G = mandrill
mandrill.B = mandrill
# zero out the non-contributing channels for each image copy
mandrill.R[,,2:3] = 0
mandrill.G[,,1]=0
mandrill.B[,,1:2]=0

df = data.frame(
  red = matrix(mandrill[,,1], ncol=1),
  green = matrix(mandrill[,,2], ncol=1),
  blue = matrix(mandrill[,,3], ncol=1)
)

### compute the k-means clustering, to obtain regions from image by color
K = kmeans(df,5,nstart = 100)

df$label = K$cluster

table(df$label)

### Replace the color of each pixel in the image with the mean 
### R,G, and B values of the cluster in which the pixel resides:
# get the coloring
colors = data.frame(
  label = 1:nrow(K$centers), 
  R = K$centers[,"red"],
  G = K$centers[,"green"],
  B = K$centers[,"blue"]
)
# merge color codes on to df
# IMPORTANT: we must maintain the original order of the df after the merge!
df$order = 1:nrow(df)
df = merge(df, colors)
df = df[order(df$order),]
df$order = NULL

# get mean color channel values for each row of the df.
R = matrix(df$R, nrow=dim(mandrill)[1])
G = matrix(df$G, nrow=dim(mandrill)[1])
B = matrix(df$B, nrow=dim(mandrill)[1])

# reconstitute the segmented image in the same shape as the input image
mandrill.segmented = array(dim=dim(mandrill))
mandrill.segmented[,,1] = R
mandrill.segmented[,,2] = G
mandrill.segmented[,,3] = B

RGB_label = R*B*G

# View the result
pdf("simulation/segmented.pdf")
grid.raster(mandrill.segmented)
dev.off()

#save(df, mandrill.segmented, R,G,B, file = "From_image_region.RData")

#> df[1:4,]
#       label red green blue         R         G         B
#200066     2   1     1    1 0.9997333 0.9994392 0.9992982
#200067     2   1     1    1 0.9997333 0.9994392 0.9992982
#200068     2   1     1    1 0.9997333 0.9994392 0.9992982
#200069     2   1     1    1 0.9997333 0.9994392 0.9992982

pixel_ind = c(1:dim(df)[1])
x_coor = rep(1:dim(mandrill.segmented)[2], each=dim(mandrill.segmented)[1])
y_coor = -rep(1:dim(mandrill.segmented)[1], times=dim(mandrill.segmented)[2])

data_groundtruth = data.frame(pixel_ind, x_coor,y_coor )
data_groundtruth$label = as.integer(as.factor(c(RGB_label)))

# remove background 
data_groundtruth_use = data_groundtruth[-which(data_groundtruth$label==5),]

times<-6

set.seed(1234)
subsample = data_groundtruth_use[sample(1:dim(data_groundtruth_use)[1],10000*times*times,replace=F),]
# > min(subsample$y_coor)
# [1] -1591
subsample$y_coor = subsample$y_coor+1591 # make coordinates >=0
subsample$x_coor = subsample$x_coor *times
subsample$y_coor = subsample$y_coor *times
save(data_groundtruth, data_groundtruth_use,subsample, file = "simulation/simulate_spatial_cell_label_6times.RData")

save(subsample, file = "~/Documents/Research/github/SpaTopic_benchmarking/simulation/rerun/data/LIBDsubsample_6times.RData") 


