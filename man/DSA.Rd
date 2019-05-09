\name{DSA}
\alias{DSA}
\title{Interface function to implement DSA}
\description{
  This is the main interface function to implement digital sorting algorithm for cell type specific gene expression profiles.
}
\usage{
DSA(mix, cell.gene, weight = NULL, method = "LM", out.cell.file = NULL, 
	out.weight.file = NULL, log2 = TRUE, l = 0, u = 2^34, meq = 0)
}
\arguments{
	\item{mix}{matrix of mixture signals (genes in row, cell type in column).}

	\item{cell.gene}{the cell-gene key table: mapping between cell type and cell specific genes. This table needs to be only two columns: first column with gene symbols, second column the cell type.}

	\item{weight}{weight matrix (samples in row, cell types in column).   
			Default to 'NULL', where no weight matrix is provided by user, and the weight will be estimated by the function.}

	\item{method}{methods used in estimating weight and true signals for each tissue type.  Default to 'LM' for linear regression.
								Other methods included 'LG' (logistic regression),'QP_LM' or 'QP_LG' (quadradic programming with constraint on the estimated parameter on linear/logistic regression).}

	\item{out.cell.file}{file name to store the deconvoluted signals.  default to 'NULL' where no file will be created.}

	\item{out.weight.file}{file name to store the estimated weights.  default to 'NULL' where no file will be created.}

	\item{log2}{flag indicating if the input mixture signals are in log2 scale.  Default to 'TRUE'.}

	\item{l,u}{values for the lower- (l) and upper- (u) bound used in setting the vector for values of b0 (bvec for solve.QP) in solving for quadratic programming. Defaults to 0 and 2^34 respectively.}
	
	\item{meq}{default to zero (used to set meq for solve.QP)}

}
\details{
	This is the main function (an interface) to implement digital sorting algorithm for cell type specific gene expression profiles, by proving minimum of an data matrix for mixture signal and cell-gene key tables.  
	You may provide specific parameters for weight estimation and deconvolution.  
}
\value{
A list object containing two elements:
	\item{est.weight}{a list of estimated weight ('estimated_weight') and the model's mean square error ('mse'). estimated_weight is a matrix of cell types in row and samples in columns. mse defaults to 'NULL' if method is not 'LM'.}

	\item{deconv}{a matrix of deconvoluted signals: genes in row, cell types in columns.}
}
\seealso{
	\code{\link{EstimateWeight}}, 
	\code{\link{Deconvolution}}, 
	\code{\link{GSM_QP}}, 
	\code{solve.QP}, 
}
\examples{
	## load package DSA
	library(DSA)
	
	## load sample data
	data('mix.signals')
	data('cell.gene')
	data('weight')
	
	# In this data, thr first three samples contain signal from only one cell	
	pure <- mix[, 1:3]
	mix <- mix[, 4:14]

	weight.orig <- weight
	weight <- weight[4:14, ]

	## Exp 1: Obtain the cell-type specific signals estimated from DSA with default parameters: 
	##			estimated weights and linear model in deconvolution
	decov.LM <- DSA(mix, cell.gene)
	
	## Exp 2: Obtain the cell-type specific signals estimated from DSA with estimated 
	##			weight and quadratic programming
	dec.QLM <- DSA(mix, cell.gene, method="QP_LM")
	
	# check cell type proportions
	t.weight <- as.numeric(weight)/100
	e.weight <- as.numeric(t(dec.QLM$est.weight))
	limits <- range(t.weight, e.weight)

	plot(t.weight, e.weight, xlim=limits, pch=19, ylim=limits)
	cor(t.weight, e.weight)
	abline(a=0, b=1)

	# check the estimated cell-type specific signals
	par(mfrow=c(1,3))
	for(i in 1:3){
		print(cor(2^(pure[,i]),2^(dec.QLM$deconv[,i])))
	}
}
