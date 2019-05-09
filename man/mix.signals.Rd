\name{mix.signals}
\alias{mix.signals.rda}
\docType{data}
\title{mixture expression signals data matrix}
\description{
  matrix of mixture signals (genes in row, cell type in column) derived from the public microarray data GSE19830.  
  The triplicates were averaged and log2-transformed. 
}
\format{
	A matrix of 31099 rows (probes) and 14 columns (samples/tissues)
}
\usage{
data('mix.signals')
}
\references{
    Shen-Orr SS, Tibshirani R, Khatri P, Bodian DL et al. 
	Cell type-specific gene expression differences in complex tissues. Nat Methods 2010 Apr;7(4):287-9. PMID: 20208531
}
