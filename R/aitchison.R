#' @name X2Xclrv
#' @title Centered Log-Ratio Velocity Transformation
#'
#' @description
#' Transform phase components using centered log-ratio velocity (clrv)
#' transformation.
#'
#' @details
#' Let \eqn{X(t)} be a phase component defined on \eqn{t \in [0,1]}. centered log-ratio velocity (clrv)
#' transformation is defined as
#' \deqn{
#' \texttt{clrv}(X)(t) := \log(DX(t)) - \int_0^1 \log(DX(s))\mathrm{d}s,
#' \quad t \in [0, 1], }
#' where \eqn{D} is a differential operator.
#'
#' @param X A matrix representing phase components, or an object of class
#'   \code{syncftn}. If a matrix is provided, it must be a
#'   \ifelse{html}{\out{<i>T</i>}\eqn{\times}\out{<i>n</i>}}{\eqn{T\times n}} matrix of <i>n</i> phase components evaluated
#'   at the time points \code{t}.
#' @param t A numeric vector of length \ifelse{html}{\out{<i>T</i>}}{\eqn{T}} giving the time points at which
#'   the phase components are evaluated. This vector must start at 0 and end at
#'   1.
#'
#' @return \code{X2Xclrv()} returns an object of class \code{Xclrv}, which is a
#'   \ifelse{html}{(\out{<i>T</i>}--1)\eqn{\times}\out{<i>n</i>}}{\eqn{(T-1)\times n}} matrix of clrv transformed phase components
#'   evaluated over the intervals defined by the time points \code{t}.
#'
#' @references Kang S. and Oh H.-S. (2026) \dQuote{Multiview representation and clustering of
#' functional data,} \emph{Unpublished Manuscript}.
#'
#' @seealso [auc_sync()] and [fr_sync()] for time-synchronizing mappings.
#'
#' @examples
#' t <- seq(0, 1, length.out = 100)
#' sync <- auc_sync(seoul_bike$Ytilde[,1:10], seoul_bike$x, t)
#' plot(X2Xclrv(sync))
#'
#' @export
X2Xclrv <- function(X, t) UseMethod("X2Xclrv")

#' @name X2Xclrv
#' @exportS3Method X2Xclrv matrix
X2Xclrv.matrix <- function(X, t){
  if(t[1] != 0L || t[length(t)] != 1L) stop("t must be an increasing sequence that starts at zero and end at one")
  DX <- diff(X) / diff(t)
  logDX <- log(DX)
  res <- t(t(logDX) - colSums(diff(t) * logDX))
  colnames(res) <- colnames(X)
  return(structure(
    res,
    t = t,
    class = "Xclrv"
  ))
}

#' @name X2Xclrv
#' @exportS3Method X2Xclrv syncftn
X2Xclrv.syncftn <- function(X, t = attr(X, "t")){
  if(t[1] != 0L || t[length(t)] != 1L) stop("t must be an increasing sequence that starts at zero and end at one")
  DX <- diff(X$X) / diff(t)
  logDX <- log(DX)
  res <- t(t(logDX) - colSums(diff(t) * logDX))
  colnames(res) <- colnames(X$X)
  return(structure(
    res,
    t = t,
    class = "Xclrv"
  ))
}

#' @name X2Xclrv
#' @aliases plot.Xclrv
#'
#' @param x A \code{Xclrv} object, obtained as a result of the function
#'   \code{X2Xclrv()}.
#' @param ... Further graphical parameters supplied to the internal
#'   [graphics::matplot()] function.
#'
#' @return \code{plot.Xclrv()} plots phase components after centered log-ratio
#'   velocity transformation.
#' @export
plot.Xclrv <- function(x, ...){
  args <- list(...)
  t <- attr(x, "t")
  if(is.null(args$ylab)) args$ylab <- ""
  if(is.null(args$type)) args$type <- "l"
  if(is.null(args$xlab)) args$xlab <- "t"
  if(is.null(args$main)) args$main <- "Phase components (clrv)"
  do.call(graphics::matplot, c(list(x = rep(t, each = 2)[-c(1,2*(length(t)))],
                                    y = matrix(rep(x, each = 2), ncol = ncol(x))), args))
}


#' @name Xclrv2X
#' @title Inverse Centered Log-Ratio Velocity Transformation
#'
#' @description
#' Inverse of the centered log-ratio velocity (clrv) transformation.
#'
#' @details
#' Let \eqn{X(t)} be a phase component defined on \eqn{t \in [0,1]}. centered log-ratio velocity (clrv)
#' transformation is defined as
#' \deqn{
#' \texttt{clrv}(X)(t) := \log(DX(t)) - \int_0^1 \log(DX(s))\mathrm{d}s,
#' \quad t \in [0, 1], }
#' where \eqn{D} is a differential operator.
#'
#' If \eqn{X: [0,1]\to[T_1,T_2]} for \eqn{T_1 < T_2}, the inverse of the clrv transformation is defined as
#' \deqn{
#' T_1 + (T_2 - T_1)\frac{\int_0^t \exp\{\texttt{clrv}(X)(s)\}\mathrm{d}s}{\int_0^1 \exp\{\texttt{clrv}(X)(s)\}\mathrm{d}s},\quad t\in[0,1].
#' }
#'
#' @param Xclrv A \ifelse{html}{(\out{<i>T</i>}--1)\eqn{\times}\out{<i>n</i>}}{\eqn{(T-1)\times n}} matrix of clrv transformed phase components
#'   evaluated over the intervals defined by the time points \code{t}.
#' @param t A numeric vector of length \ifelse{html}{\out{<i>T</i>}}{\eqn{T}} giving the time points at which
#'   the phase components are evaluated. This vector must start at 0 and end at
#'   1.
#' @param x A numeric vector of length \ifelse{html}{\out{<i>m</i>}}{\eqn{m}}
#'   giving the observed time points.
#'
#' @return A \ifelse{html}{\out{<i>T</i>}\eqn{\times}\out{<i>n</i>}}{\eqn{T\times n}} matrix of inverse-clrv transformed phase components
#'   evaluated at the time points \code{t}.
#'
#' @references Kang S. and Oh H.-S. (2026) \dQuote{Multiview representation and clustering of
#' functional data,} \emph{Unpublished Manuscript}.
#'
#' @seealso [X2Xclrv()] for clrv transformation.
#'
#' @examples
#' t <- seq(0, 1, length.out = 100)
#' sync <- auc_sync(seoul_bike$Ytilde[,1:10], seoul_bike$x, t)
#' xclrv <- X2Xclrv(sync)
#' range(Xclrv2X(xclrv, t, seoul_bike$x) - sync$X)
#'
#' @export
Xclrv2X <- function(Xclrv, t, x){
  if(t[1] != 0L || t[length(t)] != 1L) stop("t must be an increasing sequence that starts at zero and end at one")
  res <- t(t(apply(exp(Xclrv) * diff(t), 2, cumsum)) / colSums(exp(Xclrv) * diff(t)))
  rbind(0, res) * diff(range(x))  + min(x)
}
