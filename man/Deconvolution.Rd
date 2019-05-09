\name{Deconvolution}
\alias{Deconvolution}
\title{Deconvolution function}
\description{
  Function to estimate the deconvoluted signals - deconvolve the mixture signals to cell-type specific signals for each gene.
}
\usage{
Deconvolution(data, weight, method = "LM", l = 0, u = 2^34)
}
\arguments{
	\item{data}{data matrix of the mixture signals, with genes in row, cell type in column, in anti-log scale.}
	
	\item{weight}{weight matrix, with cell type in row, tissue types in column.}
	
	\item{method}{methods used in estimating true signals for each tissue type.  Default to 'LM' for linear regression.  
					Other methods included 'LG' (logistic regression),'QP_LM' or 'QP_LG' (quadradic programming with constraint on the estimated parameter on linear/logistic regression)
					For methods 'LG' and 'QP_LG', input data is transformed into log-scaled and returned values are anti-logged.}

	\item{l,u}{values for the lower- (l) and upper- (u) bound used in setting the vector for values of b0 (bvec for solve.QP) in solving for quadratic programming. Defaults to 0 and 2^34 respectively.}
}
\value{
	\item{paraM}{a matrix of deconvoluted signals (anti-log scale) with genes in row, tissue types in columns}
}
\seealso{
	\code{\link{EstimateWeight}}, 
	\code{\link{GSM_QP}}, 
	\code{solve.QP}, 
}
\examples{
	## load package DSA
	library(DSA)
	
	## load sample data
	data('mix.signals')
	data('cell.gene')

	# In this data, thr first three samples contain signal from only one cell	
	pure <- mix[, 1:3]
	mix <- mix[, 4:14]

	# Generate a list object set of marker genes for each cell type
	unique_cell_type <- unique(cell.gene[,2]) 
	gene_list <- list()
	for( i in 1 : length(unique_cell_type)){
		gene_list[[i]] <- cell.gene[cell.gene[,2] == unique_cell_type[i],1]
	} 
	names(gene_list) <- unique_cell_type

	# Estimate weight matrix
	estimated_weight <- EstimateWeight(2^mix, gene_list, method="LM")
	
	# Estimate deconvoluted signals
	deconv <- as.matrix(Deconvolution(2^mix, t(estimated_weight$weight), method="QP_LM", l=min(2^mix), u=max(2^mix)))
	rownames(deconv) <- rownames(mix)
	colnames(deconv) <- unique_cell_type
	
	# check the estimated cell-type specific signals
	for(i in 1:3){
		print(cor(2^(pure[,i]),deconv[,i]))
	}
}