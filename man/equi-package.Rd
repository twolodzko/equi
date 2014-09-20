\name{equi-package}
\encoding{utf8}
\alias{equi-package}
\docType{package}
\title{equi: A package for equipercentile equating}
\description{
Package for equipercentile equating in \emph{Equivalent Groups} (EG), \emph{Single Group} (SG) and \emph{Nonequivalent groups with Anchor Test} with \emph{Chained Equating} (NEAT-CE) designs.

The process of observed-score equating consists of several steps: (1) presmoothing, (2) estimating the score probabilities, (3) continuization (or postsmoothing), (4) computing the equating function, and (5) evaluating the equating results and computing accuracy measures (von Davier, 2011; von Davier, Holland & Thayer, 2004). The steps (1)-(3) are implemented with function \code{\link{smoothtab}}, step (4) with function \code{\link{equi}}, and step (5) with function \code{\link{see}} for bootstrap-based \emph{Standard Error of Equating} estimation.

Equipercentile equting function of scale \emph{X} to scores of scale \emph{Y} can be defined as follows:

\deqn{equi_Y(x) = F^{-1}_Y[F^{-1}_X(x)]}

where \emph{Fx} ia a cumulative distribution function of \emph{X}, and \emph{Fy} is a cumulative distribution function of \emph{Y}. For the equating function to be symmetric, \emph{Fx} and \emph{Fy} need to be continuized so their distribution is transformed from discrete into continous cumulative distribution function. In this package, \emph{Gaussian Kernel} postsmoothing could be used for this pourpose (von Davier, Holland & Thayer, 2004).

This package was desined during a research project by \enc{Wołodźko}{Wolodzko}, Kondratek & Szaleniec (2014).
}
\author{
Maintainer: Tymoteusz \enc{Wołodźko}{Wolodzko} <\email{twolodzko@gmail.com}>
}
\references{
Holland, P.W. & Thayer, D.T. (2000). Univariate and Bivariate Loglinear Models for Discrete Test Score Distributions. \emph{Journal of Educational and Behavioral Statistics, 25}(2), 133-183.

Kolen, M.J. & Brennan, R.J. (2004). \emph{Test Equating, Scaling, and Linking: Methods and Practices.} New York: Springer-Verlag.

von Davier, A.A. (Ed.) (2011). \emph{Statistical Models for Test Equating, Scaling, and Linking.} New York: Springer-Verlang.

von Davier, A.A., Holland, P.W. & Thayer, D.T. (2004). \emph{The Kernel Method of Test Equating.} New York: Springer-Verlag.

Wand, M.P. & Jones, M.C. (1995). \emph{Kernel Smoothing.} London: Chapman & Hall/CRC. 

\enc{Wołodźko}{Wolodzko}, T., Kondratek, B., & Szaleniec, H. (2014). \enc{Zrównanie}{Zrownywanie} ekwicentylowe na tle innych metod \enc{zrównywania}{zrownywania} na \enc{przykładzie}{przykladzie} sprawdzianu i egzaminu gimnazjalnego. \emph{Edukacja, 3}(128), 131-149.
}
\keyword{equating}
\keyword{psychometrics}
