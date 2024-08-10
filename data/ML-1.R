data <- readRDS("data dump.rds")
tr1 <- as.data.frame(table(data$clinical.data$Subtype_Selected))
tr1 <- tr1[tr1$Freq >= 10,] #removes 6 tumor subtypes and 44 samples#

sel <- subset(data$clinical.data, data$clinical.data$Subtype_Selected %in% tr1$Var1)
sel <- sel[,c(1,3,35)]
m <- subset(data$matrisome.data, data$matrisome.data$sample %in% sel$sample)
m <- merge(sel,m,by="sample")
m[,1] <- NULL
m[is.na(m)] <- 0

library(Rtsne)

m2 <- normalize_input(as.matrix(m[,3:ncol(m)]))

#visualization with t-SNE#

set.seed(123)
rt <- Rtsne(m2, pca_center = F, perplexity = 663, eta = 553) 
d <- data.frame(labels=m[,2],tum.type=m[,1],x=rt$Y[,1],y=rt$Y[,2])

d$histo.type <- ifelse(d$tum.type == "OV","Gyn",
                       ifelse(d$tum.type == "UCEC", "Gyn",
                              ifelse(d$tum.type == "UCS", "Gyn",
                                     ifelse(d$tum.type == "GBM", "Neuroendocrine",
                                            ifelse(d$tum.type == "LGG", "Neuroendocrine",
                                                   ifelse(d$tum.type == "PCPG","Neuroendocrine",
                                                          ifelse(d$tum.type == "ACC", "Neuroendocrine",
                                                                 ifelse(d$tum.type == "COAD", "Gi",
                                                                        ifelse(d$tum.type == "READ", "Gi",
                                                                               ifelse(d$tum.type == "STAD", "Gi-associated",
                                                                                      ifelse(d$tum.type == "LIHC", "Gi-associated",
                                                                                             ifelse(d$tum.type == "ESCA", "Squamous",
                                                                                                    ifelse(d$tum.type == "HNSC", "Squamous",
                                                                                                           ifelse(d$tum.type == "LUSC", "Squamous",
                                                                                                                  ifelse(d$tum.type == "LAML", "Blood",
                                                                                                                         ifelse(d$tum.type == "KICH", "Renal",
                                                                                                                                ifelse(d$tum.type == "KIRC", "Renal",
                                                                                                                                       ifelse(d$tum.type == "KIRP", "Renal",
                                                                                                                                              "other"))))))))))))))))))

library(ggplot2)
#global graph - option is given to add a shape code for the cell/tissue-of-orign patterns#
ggplot(data = d, aes(x,y)) + geom_point(aes(colour = labels #,
                                            #shape = histo.type
), size = 1) + 
  #scale_shape_manual(values=as.factor(d$histo.type)) +
  theme_void() #+ facet_wrap(~labels)

#major tumor types#
ggplot(data = d, aes(x,y)) + geom_point(aes(colour = tum.type), size = 1) + 
  theme_void() #+ facet_wrap(~tum.type)

#cell/tissue-of-origin patterns#
ggplot(data = d, aes(x,y) ) + geom_point(aes(colour = histo.type), size = 1) + 
  theme_bw() #+ facet_wrap(~tum.type)

#correlation table and plot#

res <- aggregate(m[,c(3:ncol(m))], by=list(m$Subtype_Selected), mean)
rownames(res) <- res[,1]
res[,1] <- NULL
res <- t(res)
res <- as.matrix(res)
res <- cor(res, method = "spearman")
library(pheatmap)
out <- pheatmap(res, fontsize = 5, cutree_rows = 5, cutree_cols = 5)

#unsupervised clustering by gaussian mixture modelling#

library(mclust)
m2 <- as.data.frame(m2)
res <- list()
set.seed(123)
for(i in 1:10){
  z <- m2[sample(nrow(m2), 660), ]
  bic <- mclustBIC(z, G=c(2:20))
  mod1 <- Mclust(z, x = bic)
  res[[i]] <- data.frame(iteration=i, N.of.components=mod1$G)
}
res <- bind_rows(res)
n <- round(mean(res[,2]),0)
k <- kmeans(m2, n)
df <- data.frame(tum.type=m$Subtype_Selected, cl=k$cluster)

library(corrplot)
library(RColorBrewer)
corrplot(cor(t(table(df$tum.type,df$cl))), tl.cex = 0.3, tl.col = "black", col = c("blue","white","red"))

#density-based clustering, using UMAP to increase performance of DBSCAN#

library(umap)
library(dbscan)

clustembed <- umap(m2,
                   n_neighbors=30,
                   min_dist=0.01,
                   n_components=10,
                   random_state=42)

plot(density(table(df$tum.type))) #density of tumor subtype sizes
abline(v=30, col=c("red"))

dbs <- hdbscan(clustembed$layout, minPts = 60) #at least twice the theoretical size of the average tumor type. Smaller values, closer to real distribution, give N clusters ) N tumor types#

df <- data.frame(tum.type=m$Subtype_Selected, cl=dbs$cluster)

corrplot(cor(t(table(df$tum.type,df$cl))), tl.cex = 0.3, tl.col = "black", col = c("blue","white","red"))