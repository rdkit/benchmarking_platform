## 
# $Id$
##
## File "Friedman_Test.Target_x_FP_Data.R"
## Feb 12, 2013
## hpg


# load required packages
# nonparametric and multple comparison testing
library(coin)
library(multcomp)
library(colorspace)

# lattice graphics
library(lattice)

# color palette for heatmap
library(RColorBrewer)

# function
friedman.test.with.post.hoc <- function(formu, data, to.print.friedman = T, 
                                        to.post.hoc.if.signif = T, to.plot.parallel = T, 
                                        to.plot.boxplot = T, signif.P = .05, 
                                        color.blocks.in.cor.plot = T, 
                                        jitter.Y.in.cor.plot =F) {
	# formu is a formula of the shape: 	Y ~ X | block
	# data is a long data.frame with three columns:    [[ Y (numeric), X (factor), block (factor) ]]
	
	# Note: This function doesn't handle NA's! In case of NA in Y in one of the blocks, 
        # then that entire block should be removed.
	
	# get the names out of the formula
	formu.names <- all.vars(formu)
	Y.name <- formu.names[1]
	X.name <- formu.names[2]
	block.name <- formu.names[3]
	
	if(dim(data)[2] >3) data <- data[,c(Y.name,X.name,block.name)]	
        # In case we have a "data" data frame with more then the three columns we need. 
        # This code will clean it from them...

	# Note: the function doesn't handle NA's. In case of NA in one of the block T outcomes, 
        # that entire block should be removed.

	# stopping in case there is NA in the Y vector
	if(sum(is.na(data[,Y.name])) > 0) stop("Function stopped: This function doesn't handle NA's. In case of NA in Y in one of the blocks, then that entire block should be removed.")
	
	# make sure that the number of factors goes with the actual values present in the data:
	data[,X.name ] <- factor(data[,X.name ])
	data[,block.name ] <- factor(data[,block.name ])
	number.of.X.levels <- length(levels(data[,X.name ]))
	if(number.of.X.levels == 2) { warning(paste("'",X.name,"'", 
                                      "has only two levels. Consider using paired wilcox.test instead of friedman test"))}
	
	# making the object that will hold the friedman test and the other.
	the.sym.test <- symmetry_test(formu, data = data,	### all pairwise comparisons	
                                      teststat = "max",
				      xtrafo = function(Y.data) { trafo( Y.data, factor_trafo = function(x) { model.matrix(~ x - 1) %*% t(contrMat(table(x), "Tukey")) } ) },
				      ytrafo = function(Y.data){ trafo(Y.data, numeric_trafo = rank, block = data[,block.name] ) }
						)
	if(to.print.friedman) { print(the.sym.test) }
	
	if(to.post.hoc.if.signif) {
            print(pvalue(the.sym.test))
	    if(pvalue(the.sym.test) < signif.P) {
		# the post hoc test
		The.post.hoc.P.values <- pvalue(the.sym.test, method = "single-step")	
                # this is the post hoc of the friedman test
									
		# plotting
		if(to.plot.parallel & to.plot.boxplot)	par(mfrow = c(1,2)) 
                # if we are plotting two plots, let's make sure we'll be able to see both
									
		if(to.plot.parallel) {
		    X.names <- levels(data[, X.name])
		    X.for.plot <- seq_along(X.names)
		    plot.xlim <- c(.7 , length(X.for.plot)+.3)	# adding some spacing from both sides of the plot
		    if(color.blocks.in.cor.plot) {
		        blocks.col <- rainbow_hcl(length(levels(data[,block.name])))
		    } else {
		        blocks.col <- 1 # black
		    }					
		    data2 <- data
		    if(jitter.Y.in.cor.plot) {
		        data2[,Y.name] <- jitter(data2[,Y.name])
	        	par.cor.plot.text <- "Parallel coordinates plot (with Jitter)"
                    } else {
			par.cor.plot.text <- "Parallel coordinates plot"
		    }				
		    # adding a Parallel coordinates plot
		    matplot(as.matrix(reshape(data2,  idvar=X.name, timevar=block.name,
					      direction="wide")[,-1])  , 
					      type = "l",  lty = 1, axes = FALSE, ylab = Y.name, 
					      xlim = plot.xlim,
					      col = blocks.col,
					      main = par.cor.plot.text)
		    axis(1, at = X.for.plot , labels = X.names) # plot X axis
		    axis(2) # plot Y axis
		    points(tapply(data[,Y.name], data[,X.name], median) ~ X.for.plot, col = "red",pch = 4, cex = 2, lwd = 5)
		}
				
		if (to.plot.boxplot) {
		    # first we create a function to create a new Y, by substracting different combinations of X levels from each other.
		    subtract.a.from.b <- function(a.b , the.data) { the.data[,a.b[2]] - the.data[,a.b[1]]}			
		    temp.wide <- reshape(data,  idvar=X.name, timevar=block.name, direction="wide") 	#[,-1]
		    wide.data <- as.matrix(t(temp.wide[,-1]))
		    colnames(wide.data) <- temp.wide[,1]				
		    Y.b.minus.a.combos <- apply(with(data,combn(levels(data[,X.name]), 2)), 
                                                2, subtract.a.from.b, the.data =wide.data)
		    names.b.minus.a.combos <- apply(with(data,combn(levels(data[,X.name]), 2)), 
                                                2, function(a.b) {paste(a.b[2],a.b[1],sep=" - ")})			
		    the.ylim <- range(Y.b.minus.a.combos)
		    the.ylim[2] <- the.ylim[2] + max(sd(Y.b.minus.a.combos))	# adding some space for the labels
		    is.signif.color <- ifelse(The.post.hoc.P.values < .05 , "green", "grey")			
		    boxplot(Y.b.minus.a.combos, names = names.b.minus.a.combos ,
				col = is.signif.color,
				main = "Boxplots (of the differences)",
				ylim = the.ylim )
		    legend("topright", legend = paste(names.b.minus.a.combos, 
                           rep(" ; PostHoc P.value:", number.of.X.levels),
                           round(The.post.hoc.P.values,5)) , 
                           fill =  is.signif.color )
		    abline(h = 0, col = "blue")			
		}
				
		list.to.return <- list(Friedman.Test = the.sym.test, PostHoc.Test = The.post.hoc.P.values)
		if(to.print.friedman) {print(list.to.return)}				
		return(list.to.return)
				
	    } else {
		print("The results where not significant, There is no need for a post hoc test")
		return(the.sym.test)
	    }					
	}

# Original credit (for linking online, to the package that performs the post hoc test) goes to "David Winsemius", see:
# http://tolstoy.newcastle.edu.au/R/e8/help/09/10/1416.html
} # friedman.test.with.post.hoc

##########################
## Start of main script ##
##########################

# get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# read data file into data frame, 'FP_Ranking.csv' is the original data file from Sereina
d.f <- data.frame(read.table(args[1], sep=",", header=T))

# create keys for categorical summaries
by1 <- unique(d.f[,1])
by2 <- unique(d.f[,2])

a.mean.by.fp <- aggregate(x = as.numeric(d.f[,3]), by = list(d.f[,1]), FUN = "mean")

a.mean.by.target <- aggregate(x = as.numeric(d.f[,3]), by = list(d.f[,2]), FUN = "mean")

a.mean.by.fp.and.target  <- aggregate(x = as.numeric(d.f[,3]), by = list(d.f[,1],d.f[,2]), FUN = "mean")

# look at 'mean ranks' of fingerprints and their distribution
ds <- a.mean.by.fp.and.target

# global friedman test		
friedman_test(x ~ Group.1 | Group.2, data = ds)		

# pairwise fp testing and p-value adjustment ... ~30 seconds per post hoc testing run (91 pairs)
res <- friedman.test.with.post.hoc(x ~ Group.1 | Group.2, 
                                   data=ds, to.print.friedman = F, 
                                   to.post.hoc.if.signif = T,  
                                   to.plot.parallel = F, 
                                   to.plot.boxplot = F, 
                                   signif.P = .05, 
                                   color.blocks.in.cor.plot = F, 
                                   jitter.Y.in.cor.plot = F)

# check dimension of post-hoc p-values in return list structure			 
#print(dim(res$PostHoc.Test))
		 
# now disect the returned results and assemble a data frame (test.res) for csv file creation
# this is a string array of the fp contrasts, separator is " - "
tt <- dimnames(res$PostHoc.Test)[[1]]
n.tt <- length(tt)
fp.1 <- rep("", n.tt)
fp.2 <- rep("", n.tt)
mean.rank.1 <- rep(NA, n.tt)
mean.rank.2 <- rep(NA, n.tt)
mean.rank.diff <- rep(NA, n.tt)
signif.test <- rep("", n.tt)

for (i in 1:n.tt) {
    tmp <- unlist(strsplit(tt[i]," - ",fixed=T))
    fp.1[i] <- tmp[1]
    fp.2[i] <- tmp[2]
    idx.1 <- which(a.mean.by.fp[,1] == fp.1[i])
    idx.2 <- which(a.mean.by.fp[,1] == fp.2[i])
    mean.rank.1[i] <- a.mean.by.fp[idx.1, 2]
    mean.rank.2[i] <- a.mean.by.fp[idx.2, 2]
    mean.rank.diff[i] <- mean.rank.1[i] - mean.rank.2[i]
}

# this is the p value array
p.val <- as.numeric(res$PostHoc.Test[,1])

# define the significance markers ... don't change order of lines!
signif.test[which(p.val >= 0.05)] <- "."	
signif.test[which(p.val > 0.1)]   <- ""	
signif.test[which(p.val < 0.05)]  <- "*"
signif.test[which(p.val < 0.01)]  <- "**"
signif.test[which(p.val < 0.001)] <- "***"

test.res <-data.frame(tt,fp.1,fp.2, mean.rank.1, mean.rank.2, mean.rank.diff, p.val, signif.test)
names(test.res) <- c("fp.contrast", "fp.1", "fp.2", "mean.rank.1", "mean.rank.2", "mean.rank.diff", "adjusted.p-value", "sig.diff")

# write the test results into a csv.file
write.csv(test.res, file=args[2])
		 

# now let's do some Bootstrap resampling ...

# d.f is original data set which we have read from file
# ds is	mean per fp and target 'cell'

fprint <- unique(factor(ds$Group.1))
target <- unique(factor(ds$Group.2))

n.f <- length(fprint)
n.t <- length(target)

# can resample in all (fp) x (target) blocks (repeats 1 .. 50) of original d.f
n.res = 100

# prepare result matrix
p.val.resampled <- matrix(NA, length(p.val), n.res)

for (i in 1:n.res) {
    print(paste("resampling run #",i))
	
    # create a new d.f.resampled
    d.f.resampled <- d.f
	
    for (j in 1:n.f) {
        idx.f <- which(as.character(d.f[,1]) == as.character(fprint[j]))
        for (k in 1:n.t) {
	    idx.t <- which(as.character(d.f[,2]) == as.character(target[k]))	
	    idx   <- intersect(idx.f, idx.t)
	  
	    d.f.resampled[idx, 3] <- sample(d.f[idx, 3], 50, replace = TRUE, prob = NULL)
        }
    }
		
    ds.resampled  <- aggregate(x = as.numeric(d.f.resampled[,3]), 
                               by = list(d.f.resampled[,1],d.f.resampled[,2]), 
                               FUN = "mean")
    res.resampled <- friedman.test.with.post.hoc(x ~ Group.1 | Group.2, data=ds.resampled, 
                                                 to.print.friedman = F, 
                                                 to.post.hoc.if.signif = T,  
						 to.plot.parallel = F, 
                                                 to.plot.boxplot = F, 
                                                 signif.P = .05, 
                                                 color.blocks.in.cor.plot = F, 
                                                 jitter.Y.in.cor.plot =F)	

    # vector of p values
    p.val.res <- as.numeric(res.resampled$PostHoc.Test[,1])

    # replace 0 by 1e-16	
    p.val.res[which(p.val.res < 1.0E-16)] <- 1.0E-16
    p.val.resampled[,i] <- p.val.res
	
} # n.res, end resampling loop
	
#pval.res <- data.frame(p.val.resampled)	
pval.res <- data.frame(test.res[,2], test.res[,3], test.res[,4], test.res[,5], p.val.resampled)	
write.csv(pval.res, file=args[3])

# done ....

	
# EOF	
