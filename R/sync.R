#' @name syncftn
#' @aliases auc_sync fr_sync
#' @title Time-Synchronizing Mappings
#'
#' @description
#' Area-under-the-curve (AUC) time-synchronizing mapping and Fisher--Rao (FR)
#' time-synchronizing mapping.
#'
#' @details
#' Let \eqn{\tilde{Y}(x)} be an observed function defined on \eqn{x \in [T_1,
#' T_2]}. The AUC time-synchronizing mapping is defined as
#' \deqn{
#' \varphi_{\tilde{Y}}(x)
#' := \left( \frac{\int_{T_1}^x |\tilde{Y}(s)|^p \, \mathrm{d}s}
#'      {\int_{T_1}^{T_2} |\tilde{Y}(s)|^p \, \mathrm{d}s} \right)^{1/p},
#' \quad x \in [T_1, T_2]. }
#' For further details, see Liu and Müller (2004) and Kang and Oh (2026).
#'
#' The FR time-synchronizing mapping is defined as
#' \deqn{
#' \varphi_{\tilde{Y}} := \operatorname*{argmin}_{\varphi} d(\tilde{Y}\circ \varphi^{-1}, Y_0),}
#' where
#' \deqn{
#' d(Y_1,Y_2):=\left\|\operatorname{sgn}(DY_1)\sqrt{\lvert DY_1\rvert} - \operatorname{sgn}(DY_2)\sqrt{\lvert DY_2\rvert} \right\|_2,}
#' with \eqn{\operatorname{sgn}(u)=1} if \eqn{u\ge 0} and \eqn{-1} otherwise.
#' The phase component is then defined as \eqn{X(t) :=
#' \varphi_{\tilde{Y}}^{-1}(t)}, and the amplitude
#' component is defined as \eqn{Y(t) := (\tilde{Y} \circ
#' \varphi_{\tilde{Y}}^{-1})(t)} for \eqn{t\in[0,1]}. For further details, see
#' Srivastava et al. (2011) and Kang and Oh (2026).
#'
#' @param Ytilde A \ifelse{html}{\out{<i>m</i>}\eqn{\times}\out{<i>n</i>}}{\eqn{m\times n}}
#'   matrix whose \ifelse{html}{\out{<i>i</i>}}{\eqn{i}}th column
#'   contains the values of the \ifelse{html}{\out{<i>i</i>}}{\eqn{i}}th observed function evaluated at the
#'   \ifelse{html}{\out{<i>m</i>}}{\eqn{m}} time points \code{x}.
#' @param x A numeric vector or a \code{syncftn} object, depending on the function.
#' \itemize{
#'  \item{\code{auc_sync()} and \code{fr_sync()}}: A numeric vector of length
#'  \ifelse{html}{\out{<i>m</i>}}{\eqn{m}} giving the observed time points
#'  corresponding to \code{Ytilde}. \item{\code{plot.syncftn()}}: A
#'  \code{syncftn} object, obtained as a result of the function
#'  \code{auc_sync()} or \code{fr_sync()}.
#' }
#' @param t A numeric vector of length \ifelse{html}{\out{<i>T</i>}}{\eqn{T}} giving the time points at which
#'   the phase and amplitude components are evaluated. This vector must start at
#'   0 and end at 1.
#' @param p A numeric value specifying the power parameter of the AUC
#'   time-synchronizing mapping. The default is \code{p = 1}.
#'
#' @return \code{auc_sync()} and \code{fr_sync()} returns a object of class
#'   \code{syncftn}, which is a list containing the following components,
#'   obtained from AUC synchronization and FR synchronization, respectively:
#'   \item{\code{X}}{A \ifelse{html}{\out{<i>T</i>}\eqn{\times}\out{<i>n</i>}}{\eqn{T\times n}} matrix of phase components evaluated
#'     at the time points \code{t}.}
#'   \item{\code{Y}}{A \ifelse{html}{\out{<i>T</i>}\eqn{\times}\out{<i>n</i>}}{\eqn{T\times n}} matrix of amplitude components
#'     evaluated at the time points \code{t}.}
#'
#' @references Kang S. and Oh H.-S. (2026) \dQuote{Multiview functional
#'   clustering using latent representations of phase and amplitude components,}
#'   \emph{Unpublished Manuscript}.
#'
#' Liu X. and Müller H.-G. (2004).\dQuote{Functional convex averaging and
#' synchronization for time-warped random curve,} \emph{Journal of the American
#' Statistical Association}, \strong{99}(467), 687--699.
#'
#' Srivastava A., Wu W., Kurtek S., Klassen E., and Marron J. S. (2011)
#' \dQuote{Registration of functional data using Fisher--Rao metric,} \emph{arXiv preprint
#' arXiv:1103.3817}.
#'
#' @seealso [fdasrvf::pair_align_functions()] for FR synchronization method.
#'   [X2Xclrv()] for centered log-ratio velocity transformation.  [fkmeans()]
#'   and [fkmedians()] for \eqn{k}-centroids clustering using phase and
#'   amplitude components.
#'
#' @examples
#' t <- seq(0, 1, length.out = 100)
#' sync <- auc_sync(seoul_bike$Ytilde[,1:10], seoul_bike$x, t)
#' plot(sync)
#' par(mfrow = c(1,2))
#' plot(sync, col = 1)
#'
#' template <- 5 * dnorm(t, 0.2, 0.1) + 5 * dnorm(t, 0.8, 0.1)
#' sync <- fr_sync(seoul_bike$Ytilde[,1:10], seoul_bike$x, t, template)
#' plot(sync, col = 1)
#' lines(t, template, col = 2)
#'
#' @export
auc_sync <- function(Ytilde, x, t, p = 1){
  if(t[1] != 0L || t[length(t)] != 1L) stop("t must be an increasing sequence that starts at zero and end at one")
  integrate.y <- apply(Ytilde, 2, function(v) sum(diff(x) * (abs(v[-1])^p + abs(v[-nrow(Ytilde)])^p) / 2))
  sync.y <- apply(Ytilde, 2, function(v) cumsum(diff(x) * (abs(v[-1])^p + abs(v[-nrow(Ytilde)])^p) / 2))
  Xinv <- rbind(rep(0, ncol(Ytilde)), (t(t(sync.y) / integrate.y))^(1/p))
  X <- apply(Xinv, 2, function(v) {
    stats::approx(x = v, y = x, xout = t)$y
  })
  Y <- sapply(1:ncol(Ytilde), function(i) {
    stats::approx(x = x, y = Ytilde[,i], xout = X[,i])$y
  })
  colnames(X) <- colnames(Y) <- colnames(Ytilde)
  return(structure(list(X = X, Y = Y), t = t, class = "syncftn"))
}

#' @name syncftn
#'
#' @param template A numeric vector of length \ifelse{html}{\out{<i>T</i>}}{\eqn{T}} giving the template
#'   function value evaluated at the time points \code{t}.
#'
#' @export
fr_sync <- function(Ytilde, x, t, template){
  if(t[1] != 0L || t[length(t)] != 1L) stop("t must be an increasing sequence that starts at zero and end at one")
  if(length(template) != length(t)) stop("length of t must be equivalent to length(template)")
  X <- matrix(nrow = length(t), ncol = ncol(Ytilde))
  Y <- matrix(nrow = length(t), ncol = ncol(Ytilde))
  template.x <- stats::approx(x = t / diff(range(t)) * diff(range(x)) - t[1] + x[1], y = template, xout = x)$y
  for(i in 1:ncol(Y)){
    f.warp <- fdasrvf::pair_align_functions(template.x, Ytilde[,i], time = x)
    X[,i] <- stats::approx((x - x[1])/(x[length(x)] - x[1]), (x[length(x)] - x[1]) * f.warp$gam + x[1], xout = t)$y
    Y[,i] <- stats::approx((x - x[1])/(x[length(x)] - x[1]), f.warp$f2tilde, xout = t)$y
  }
  colnames(X) <- colnames(Y) <- colnames(Ytilde)
  return(structure(list(X = X, Y = Y), t = t, class = "syncftn"))
}

#' @name syncftn
#' @aliases plot.syncftn
#'
#' @param phase_mode A character string.
#'  * `raw` (the default): Plots phase components in their original form.
#'  * `clrv`: Plots phase components after centered log-ratio velocity transformation. Refer to [X2Xclrv()].
#' @param ... Further graphical parameters supplied to the internal
#'   [graphics::matplot()] function. \code{main} and \code{ylab} arguments are ignored.
#'
#' @return \code{plot.syncftn()} plots phase and amplitude components of each
#'   observed function.
#' @export
plot.syncftn <- function(x, phase_mode = c("raw", "clrv"), ...){
  mode <- match.arg(tolower(phase_mode), choices = c("raw", "clrv"))
  args <- list(...)
  args$ylab <- ""
  if(is.null(args$type)) args$type <- "l"
  if(is.null(args$xlab)) args$xlab <- "t"
  withr::local_par(ask = TRUE)
  if(mode == "raw"){
    args$main <- "Phase components"
    do.call(graphics::matplot, c(list(attr(x,"t"), x$X), args))
  }else{
    args$main <- "Phase components (clrv)"
    do.call(plot.Xclrv, c(list(X2Xclrv.syncftn(x)), args))
  }
  args$main <- "Amplitude components"
  do.call(graphics::matplot, c(list(attr(x,"t"), x$Y), args))
}

