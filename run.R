# setup
# clear the environment
rm(list=ls())

DATA_DIR <- './data'
IMAGES_DIR <- './images'
OUTPUT_DIR <- './output'

make_dir <- function(d) {
  if (file.exists(d)) unlink(d, recursive=TRUE, force=TRUE)
  dir.create(d)
}
lapply(c(IMAGES_DIR, OUTPUT_DIR),make_dir)


## function that concatenates strings (useful for directory paths)
concat <- function(x1,x2) {
  result <- paste(x1,x2,sep="")
  return(result)
}

## function that checks to see if a package is installed and,if not,installs it
## portions of this code came from http://stackoverflow.com/questions/9341635/how-can-i-check-for-installed-r-packages-before-running-install-packages
load_package <- function(x) {
  if (x %in% rownames(installed.packages())) { 
    print(concat("package already installed: ", x))
  }
  else { 
    install.packages(x) 
  }
  library(x, character.only=TRUE)
}

###################################################
# Kaiser-Meyer-Olkin’s sampling adequacy criteria #
# R practice: Factor analysis                     #
# Minato Nakazawa (minato-nakazawa@umin.net)      #
# 27 June 2011                                    #
###################################################
# Tests whether there are a significant number of factors in the dataset. Technically, 
# tests the ratio of item-correlations to partial item correlations. If the partials 
# are similar to the raw correlations, it means the item doesn’t share much variance 
# with other items. The range of KMO is from 0.0 to 1.0 and desired values 
# are > 0.5. Variables with MSA being below 0.5 indicate that item does not belong
# to a group and may be removed form the factor analysis.
kmo <- function(x)
{
  # Prof. Shigenobu Aoki
  # http://aoki2.si.gunma-u.ac.jp/R/kmo.html
  x <- subset(x, complete.cases(x))  # Omit missing values
  r <- cor(x)                        # Correlation matrix
  r2 <- r^2                          # Squared correlation coefficients
  i <- solve(r)                      # Inverse matrix of correlation matrix
  d <- diag(i)                       # Diagonal elements of inverse matrix
  p2 <- (-i/sqrt(outer(d, d)))^2     # Squared partial correlation coefficients
  diag(r2) <- diag(p2) <- 0          # Delete diagonal elements
  KMO <- sum(r2)/(sum(r2)+sum(p2))
  MSA <- colSums(r2)/(colSums(r2)+colSums(p2))
  return(list(KMO=KMO, MSA=MSA))
}

###################################################
# Bartlett’s sphericity test                      #
# R practice: Factor analysis                     #
# Minato Nakazawa (minato-nakazawa@umin.net)      #
# 27 June 2011                                    #
###################################################
# Tests the hypothesis that correlations between variables are greater than would be 
# expected by chance: Technically, tests if the matrix is an identity matrix. 
# The p-value should be significant: i.e., the null hypothesis that all off-diagonal 
# correlations are zero is falsified.
Bartlett.sphericity.test <- function(x)
{
  # Prof. Shigenobu Aoki
  # http://aoki2.si.gunma-u.ac.jp/R/Bartlett.sphericity.test.html
  method <- "Bartlett’s test of sphericity"
  data.name <- deparse(substitute(x))
  x <- subset(x, complete.cases(x)) # Omit missing values
  n <- nrow(x)
  p <- ncol(x)
  chisq <- (1-n+(2*p+5)/6)*log(det(cor(x)))
  df <- p*(p-1)/2
  p.value <- pchisq(chisq, df, lower.tail=FALSE)
  names(chisq) <- "X-squared"
  names(df) <- "df"
  return(structure(list(statistic=chisq, parameter=df, p.value=p.value,
                        method=method, data.name=data.name), class="htest"))
}

# get some data
load_package("foreign")  # used to load a .sav SPSS file
y <- read.spss("http://www.subjectpool.com/ed_teach/y3method/factorexdata05.sav")
CURRENT_DATE <- Sys.time() # US CST (Chicago)
x <- as.data.frame(y)
# how many are n/a?
sum(is.na(x))
head(which(is.na(x)))
x[14647,]
# how many are NULL?
sum(is.null(x))
# how many are blank?
length(which(x == ""))
# 999 is an indication that the data is NA
for (i in 1:length(x)) { x[,i] <- ifelse(x[,i]==999,NA,x[,i]) }
# how many are n/a?
sum(is.na(x))
head(which(is.na(x)))
x[541,]
# The data 'x' consists of 538 cases with 102 variables.
# it can be saved as "factorexdata05.txt" by the following line
# write.table(x,"factorexdata05.txt",quote=FALSE,sep="\t",row.names=FALSE)
# if so, the data can be read by:
# x <- read.delim("factorexdata05.txt")
write.table(x,concat(DATA_DIR, "/factorexdata05.txt"),quote=FALSE,sep="\t",row.names=FALSE)

Ps <- x[,4:43] # Extract variables p1-p40
Ps <- subset(Ps, complete.cases(Ps)) # Omit missings (511 cases remain)
load_package("rela")
res <- paf(as.matrix(Ps))
summary(res)  # Automatically calculate KMO with MSA, determine the number of factors,
# calculate chi-square of Bartlett’s sphericity test, communalities and
# factor loadings.  Communalities are 1 minus uniquenesses.
(kmo(x))
(Bartlett.sphericity.test(x))
barplot(res$Eigenvalues[,1]) # First column of eigenvalues.
plot(res$Eigenvalues[,1],type="b")
abline(a=1,b=0,col="red")
(resv <- varimax(res$Factor.Loadings)) # Varimax rotation is possible later.
barplot(sort(colSums(loadings(resv)^2),decreasing=TRUE)) # screeplot using rotated SS loadings.
scores <- as.matrix(Ps) %*% as.matrix(resv$loadings) # Get factor scores in a simple manner.
load_package("psych")
cortest.bartlett(Ps) # Bartlett’s sphericity test.
res2 <- fa.parallel(Ps)
load_package("GPArotation")
(res3 <- fa(Ps, fm="minres", nfactors=8, rotate="oblimin")) # Factor loadings as $loadings

# Maximum Likelihood Factor Analysis
# entering raw data and extracting 9 factors, 
# with varimax rotation 
fit <- factanal(Ps, 9, rotation="varimax")
print(fit, digits=2, cutoff=.3, sort=TRUE)
# plot factor 1 by factor 2 
load <- fit$loadings[,1:2] 
plot(load,type="n") # set up plot 
text(load,labels=names(Ps),cex=.7) # add variable names

# Principal Axis Factor Analysis
fit <- fa(Ps, fm="pa", nfactors=9, rotate="varimax", scores="Bartlett")
fit   # print results

# Determine Number of Factors to Extract
load_package("nFactors")
ev <- eigen(cor(Ps))  # get eigenvalues
ap <- parallel(subject=nrow(Ps),var=ncol(Ps),
               rep=100,cent=.05)
nS <- nScree(x=ev$values, aparallel=ap$eigen$qevpea)
plotnScree(nS)

# PCA Variable Factor Map 
load_package("FactoMineR")
par(mfrow=c(1,2))
res.pca <- PCA(Ps)
dev.off()
