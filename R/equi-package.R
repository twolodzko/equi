

#' equi: A package for equipercentile equating
#' 
#' Package for equipercentile equating in \emph{Equivalent Groups} (EG),
#' \emph{Single Group} (SG) and \emph{Nonequivalent groups with Anchor Test}
#' with \emph{Chained Equating} (NEAT-CE) designs.
#' 
#' The process of observed-score equating consists of several steps: (1) presmoothing,
#' (2) estimating the score probabilities, (3) continuization (or postsmoothing),
#' (4) computing the equating function, and (5) evaluating the equating results and
#' computing accuracy measures (von Davier, 2011; von Davier, Holland & Thayer, 2004).
#' The steps (1)-(3) are implemented with function \code{\link{smoothtab}}, step
#' (4) with function \code{\link{equi}}, and step (5) with function \code{\link{see}}
#' for bootstrap-based \emph{Standard Error of Equating} estimation.
#' 
#' Equipercentile equting function of scale \emph{X} to scores of scale \emph{Y}
#' can be defined as follows:
#' 
#' \deqn{equi_Y(x) = F^{-1}_Y[F_X(x)]}
#' 
#' where \emph{Fx} ia a cumulative distribution function of \emph{X}, and \emph{Fy}
#' is a cumulative distribution function of \emph{Y}. For the equating function to be
#' symmetric, \emph{Fx} and \emph{Fy} need to be continuized so their distribution is
#' transformed from discrete into continous cumulative distribution function. In this
#' package, \emph{Gaussian Kernel} postsmoothing could be used for this pourpose
#' (von Davier, Holland & Thayer, 2004).
#' 
#' This package was desined during a research project by \enc{Wołodźko}{Wolodzko},
#' Kondratek & Szaleniec (2014).
#' 
#' 
#' @docType package
#' @name equi-package
NULL
 


#' Lower-secondary school students exam scores
#' 
#' This dataset is an extract taken from a larger study of examination results
#' conducted by the Student Performance Analysis Section of the Educational Research
#' Institute (Szaleniec et al., 2013) on random sample of Polish schools. The dataset
#' consists of two samples out of eleven given in the oryginal research.
#' 
#' @format A data frame with 1705 observations on the following 5 variables.
#' 
#' \itemize{
#'    \item x,y,z  numeric vectors of exam scores.
#'    \item Sample factor vector inidcating one of the two samples.
#'    \item School factor vector indicating one of the randomly sampled Schools. 
#' }
#' 
#' @references
#' 
#' Szaleniec, H., Grudniewska, M., Kondratek, B., Kulon, F., Pokropek, A.,
#' \enc{Stożek}{Stozek}, E. & \enc{Żółtak}{Zoltak}, M. (2013).
#' \href{http://eduentuzjasci.pl/pl/publikacje-ee-lista/raporty/181-raport-z-badania/analiza-porownawcza-wynikow-egzaminow-zewnetrznych/946-analiza-porownawcza-wynikow-egzaminow-zewnetrznych-sprawdzian-w-szostej-klasie-szkoly-podstawowej-i-egzamin-gimnazjalny.html}{Analiza \enc{Porównawcza Wyników Egzaminów Zewnętrznych - Sprawdzian w Szóstej Klasie Szkoły Podstawowej i Egzamin Gimnazjalny}{Porownawcza Wynikow Egzaminow Zewnetrznych - Sprawdzian w Szostej Klasie Szkoly Podstawowej i Egzamin Gimnazjalny}}.
#' Educational Research Institute. \url{http://ibe.edu.pl/}
#' 
#' @name Tests
NULL


 



