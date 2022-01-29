## Hierarchical Clustering
## Katherine Grisanzio
## March - 2019


#Load packages
library(factoextra)
library(cluster)
library(vegan)

#Create fake dataset for demonstration purposes
set.seed(333)
df <- data.frame(replicate(3, rnorm(n = 100, mean = 0, sd = 3)))


# DETERMINE NUMBER OF CLUSTERS

#Gap Statistic - Tibshirani (2001) recommendations
  #The gap statistic compares the change in within-cluster dispersion with that expected under a 
  #reference null distribution (a distribution with no obvious clustering)
  #The optimal number of clusters indicated by this metric is the solution that yields the largest 
  #gap statistic, signifying the clustering solution is far from a null distribution 
  #B = the number of Monte Carlo bootstrapping samples (“B” copies of the reference data sets)

gap_stat <- clusGap(df, FUN=hcut, K.max = 50, B = 50, d.power = 2, verbose = TRUE, method = "Tibs2001SEmax", SE.factor=1)

fviz_gap_stat(gap_stat, linecolor = "black")


#Elbow plot
  #The elbow method, a graphical method that shows the percentage of variance explained as a function of 
  #number of clusters.
  #Based on this method, the number of clusters should be chosen in such a way that adding an additional 
  #cluster doesn’t significantly improve the modeling of the data (or percentage of variance explained). 
  #When plotted, this point can be identified by locating the “elbow” in the plot.

fviz_nbclust(df, hcut, method="wss", linecolor = "black")


#Calinski-Harabasz
  #Identifies the number of clusters based on the ratio of between-cluster variance to within-cluster variance.
  #Larger scores denote more optimal clustering solutions since it indicates both a large separation between 
  #clusters and low separation within clusters.

#Script from an online source to generate calinski and wss graphs from 1-20 clusters
Distance <- function(cluster)
{
  # the center of the cluster, mean of all the points
  center <- colMeans(cluster)
  
  # calculate the summed squared error between every point and 
  # the center of that cluster 
  distance <- apply( cluster, 1, function(row)
  {
    sum( ( row - center )^2 )
  }) %>% sum()
  
  return(distance)
}

# calculate the within sum squared error manually for hierarchical clustering 
# [WSS] : pass in the dataset, and the resulting groups(cluster)

WSS <- function( data, groups )
{
  k <- max(groups)
  
  # loop through each groups (clusters) and obtain its 
  # within sum squared error 
  total <- lapply( 1:k, function(k)
  {
    # extract the data point within the cluster
    cluster <- subset( data, groups == k )
    
    distance <- Distance(cluster)
    return(distance)
  }) %>% unlist()
  
  return( sum(total) )
}


CHCriterion <- function( data, kmax, clustermethod, ...  )
{
  if( !clustermethod %in% c( "kmeanspp", "hclust" ) )
    stop( "method must be one of 'kmeanspp' or 'hclust'" )
  
  # total sum squared error (independent with the number of cluster k)
  tss <- Distance( cluster = data )
  
  # initialize a numeric vector storing the score
  wss <- numeric(kmax)
  
  # k starts from 2, cluster 1 is meaningless
  if( clustermethod == "kmeanspp" )
  {
    for( k in 2:kmax )
    {
      results <- Kmeanspp( data, k, ... )
      wss[k]  <- results$tot.withinss 
    }		
  }else # "hclust"
  {
    for( k in 2:kmax )
    {
      d <- dist( data, method = "euclidean" )
      clustering <- hclust( d, ... )
      groups <- cutree( clustering, k )
      wss[k] <- WSS( data = data, groups =  groups )
    }
  }		
  
  # between sum of square
  bss <- tss - wss[-1]
  
  # cluster count start from 2! 
  numerator <- bss / ( 1:(kmax-1) )
  denominator <- wss[-1] / ( nrow(data) - 2:kmax )
  
  criteria <- data.frame( k = 2:kmax,
                          CHIndex = numerator / denominator,
                          wss = wss[-1] )
  
  # convert to long format for plotting 
  criteria_long <- gather( criteria, "index", "value", -1 )
  
  plot <- ggplot( criteria_long, aes( k, value, color = index ) ) + 
    geom_line() + geom_point( aes( shape = index ), size = 3 ) +
    facet_wrap( ~ index, scale = "free_y" ) + 
    guides( color = FALSE, shape = FALSE )
  
  return( list( data = criteria, 
                plot = plot ) )
}

ClusterMethod <- function( data, k, noise.cut = 0, clustermethod, ... )
{
  if( !clustermethod %in% c( "kmeanspp", "hclust" ) )
    stop( "method must be one of 'kmeanspp' or 'hclust'" )
  
  # hierarchical clustering 
  if( clustermethod == "hclust" )
  {
    cluster   <- hclust( dist(data), ... )
    partition <- cutree( cluster, k )
    
  }else # kmeanspp
  {
    cluster   <- Kmeanspp( data = data, k = k, ... )
    partition <- cluster$cluster
  }	
  
  # equivalent to k
  cluster_num <- max(partition) 
  
  # calculate each cluster's size 
  cluster_size <- numeric(cluster_num)
  for( i in 1:cluster_num )
    cluster_size[i] <- sum( partition == i )
  
  # if there're cluster size smaller than the specified noise.cut, do :
  not_noise_num <- sum( cluster_size > noise.cut )
  
  if( cluster_num > not_noise_num )
  {
    # extract the cluster whose size is larger than noise.cut
    cluster_new <- (1:cluster_num)[ cluster_size > noise.cut ]
    
    # all the data points whose original cluster is smaller than the noise.cut
    # will be assigned to the same new cluster
    cluster_num <- not_noise_num + 1
    
    # new clustering number, assign the noise cluster's number first
    # then adjust the original cluster's number
    new <- rep( cluster_num, nrow(data) )
    
    for( i in 1:not_noise_num )
      new[ ( partition == cluster_new[i] ) ] <- i
    
    partition <- new
  }
  
  # boolean vector indicating which data point belongs to which cluster
  cluster_list <- lapply( 1:cluster_num, function(x)
  {
    return( partition == x )
  })
  
  cluster_result <- list( result      = cluster,	                        
                          partition   = partition,
                          clusternum  = cluster_num,
                          clusterlist = cluster_list )
  return(cluster_result)
}


criteria <- CHCriterion(data = df, kmax = 20, 
                        clustermethod = "hclust", method = "ward.D2")

criteria$plot



#Silhouette
  #A silhouette analysis is a commonly used metric of cluster cohesion (how similar an object 
  #is to its own cluster) compared to separation (how far it is to neighboring clusters)

fviz_nbclust(df, hcut, method="silhouette", linecolor = "black")



## RUNNING CLUSTERING USING HCLUST

#Dendrogram
#A tree diagram that shows the relative similarity between cases, and are organized into branches 
#that represent the clusters.

d <- dist(df, method = "euclidean")
fit <- stats::hclust(d, method = "ward.D2", members = NULL) ### HAD TO ADD THE PACKAGE

dend <- as.dendrogram(fit) 
plot(dend)

#Let's say our fake data shows evidence for 5 clusters - 

five_clusters <- cutree(fit, h=16) #Assign clusters by cutting tree at 5 clusters
df$five_clusters <- five_clusters #Save cluster assignment to dataframe

#Can plot the silhouette in more detail now that we have cluster number determined
sil = silhouette (five_clusters,d)
plot(sil)




