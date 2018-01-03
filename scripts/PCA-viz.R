####
#
#
## PCA of KEGGs using only the KEGGs found from literature (or best of lists)
##
#
#
##

setwd("~/Documents/KEGG-IEEE-BigData/")

df.patients = read.csv("./data/table-kegg-clean-transpose.csv", header = TRUE, stringsAsFactors = FALSE)
numeric.df = df.patients[,-c(1, ncol(df.patients))]



pca = prcomp(log10(1e-9 + as.matrix(numeric.df)))
pca.res = cbind(df.patients$subject.type, df.patients$X, as.data.frame(pca$x))
names(pca.res)[2] = "subject.ID"
names(pca.res)[1] = "subject.type"
ggplot(pca.res, aes(x = PC1, y = PC2, colour = subject.type, label = subject.ID)) + 
  geom_point(size = 3) + 
  ggtitle("Original PCA") + 
  theme(legend.title=element_blank())-> p
print(p)


ggplotly(p)


### 
##

literature.keggs = read.csv("~/Google Drive/KEGGs-dec-2017/data/ibd_uc_cd_keggs.txt", header = FALSE, stringsAsFactors = FALSE)

col.names = names(numeric.df)
col.names.keggs = sapply(col.names, FUN = function(x) strsplit(x, split = "\\.")[[1]][1])
names(numeric.df) = col.names.keggs

pca = prcomp(log10(1e-9 + as.matrix(numeric.df[,intersect(literature.keggs$V1, col.names.keggs)])))
pca.res = cbind(df.patients$subject.type, df.patients$X, as.data.frame(pca$x))
names(pca.res)[2] = "subject.ID"
names(pca.res)[1] = "subject.type"
ggplot(pca.res, aes(x = PC1, y = PC2, colour = subject.type, label = subject.ID)) + 
  ggtitle("Using only KEGGs found from literature (22 out of 89 found in our data)") + 
  geom_point(size = 3) + 
  theme(legend.title=element_blank())-> p
print(p)


ggplotly(p)



#####

extreme.keggs = read.csv("~/Google Drive/KEGGs-dec-2017/results/KEGG-extreme-HE-ratios-19-12-2017.csv", header = TRUE, stringsAsFactors = FALSE)
extreme.keggs.clean = sapply(extreme.keggs$KEGG_names, FUN = function(x) strsplit(x, split = "\\(")[[1]][1])

pca = prcomp(log10(1e-9 + as.matrix(numeric.df[,intersect(extreme.keggs.clean, col.names.keggs)])))
pca.res = cbind(df.patients$subject.type, df.patients$X, as.data.frame(pca$x))
names(pca.res)[2] = "subject.ID"
names(pca.res)[1] = "subject.type"
ggplot(pca.res, aes(x = PC1, y = PC2, colour = subject.type, label = subject.ID)) + 
  geom_point(size = 3) + 
  ggtitle("Using extreme 67 KEGGs") +
  theme(legend.title=element_blank())-> p
print(p)


ggplotly(p)



###

line.keggs = read.csv("~/Google Drive/KEGGs-dec-2017/results/discriminative-KEGGs-median-HE-ratios-15-12-2017.csv", header = TRUE, stringsAsFactors = FALSE)
line.keggs.clean = sapply(line.keggs$KEGG_names, FUN = function(x) strsplit(x, split = "\\(")[[1]][1])

pca = prcomp(log10(1e-9 + as.matrix(numeric.df[,intersect(line.keggs.clean, col.names.keggs)])))
pca.res = cbind(df.patients$subject.type, df.patients$X, as.data.frame(pca$x))
names(pca.res)[2] = "subject.ID"
names(pca.res)[1] = "subject.type"
ggplot(pca.res, aes(x = PC1, y = PC2, colour = subject.type, label = subject.ID)) + 
  geom_point(size = 3) + 
  ggtitle("Using line 1513 KEGGs") +
  theme(legend.title=element_blank())-> p
print(p)


ggplotly(p)



####
# ayasdi group 4
####
group4 = read.csv("~/Google Drive/KEGGs-dec-2017/results/ayasdi/ayasdi-group-4_21-12-2017.csv", header = TRUE, stringsAsFactors = FALSE)
group4.keggs.clean = sapply(group4$KEGG_names, FUN = function(x) strsplit(x, split = "\\(")[[1]][1])

pca = prcomp(log10(1e-9 + as.matrix(numeric.df[,intersect(group4.keggs.clean, col.names.keggs)])))
pca.res = cbind(df.patients$subject.type, df.patients$X, as.data.frame(pca$x))
names(pca.res)[2] = "subject.ID"
names(pca.res)[1] = "subject.type"
ggplot(pca.res, aes(x = PC1, y = PC2, colour = subject.type, label = subject.ID)) + 
  geom_point(size = 3) + 
  ggtitle("Ayasdi Group 4 KEGGs") +
  theme(legend.title=element_blank())-> p
print(p)


ggplotly(p)


####
# ayasdi group 6
####
group6 = read.csv("~/Google Drive/KEGGs-dec-2017/results/ayasdi/ayasdi-group-6_21-12-2017.csv", header = TRUE, stringsAsFactors = FALSE)
group6.keggs.clean = sapply(group6$KEGG_names, FUN = function(x) strsplit(x, split = "\\(")[[1]][1])

pca = prcomp(log10(1e-9 + as.matrix(numeric.df[,intersect(group6.keggs.clean, col.names.keggs)])))
pca.res = cbind(df.patients$subject.type, df.patients$X, as.data.frame(pca$x))
names(pca.res)[2] = "subject.ID"
names(pca.res)[1] = "subject.type"
ggplot(pca.res, aes(x = PC1, y = PC2, colour = subject.type, label = subject.ID)) + 
  geom_point(size = 3) + 
  ggtitle("Ayasdi Group 6 KEGGs") +
  theme(legend.title=element_blank())-> p
print(p)


ggplotly(p)

