#' @name fkmeans
#' @aliases fkmeans_pre print.fkmeans fitted.fkmeans
#' @title Functional \eqn{k}-Means Clustering Using Phase and Amplitude
#'   Components
#'
#' @description
#' Conducts \ifelse{html}{\out{<i>k</i>}}{\eqn{k}}-means clustering by jointly
#' considering phase and amplitude variation. The relative importance of the two
#' components can be explicitly controlled by the user via the multiview
#' parameter \eqn{\alpha}. Optionally,
#' \ifelse{html}{\out{<i>k</i>}}{\eqn{k}}-means clustering can be performed
#' directly on the observed curves, rather than on their phase and amplitude
#' components. See Details below.
#'
#' @details
#' The distance between two observed functions is defined in terms of their
#' phase and amplitude components. For two functions with components
#' \eqn{(X_1, Y_1)} and \eqn{(X_2, Y_2)}, the distance is given by
#' \deqn{
#' \left\{
#'   \alpha \|\texttt{clrv}(X_1) - \texttt{clrv}(X_2)\|_2^2
#'   + \|Y_1 - Y_2\|_2^2
#' \right\}^{1/2},
#' }
#' where \eqn{\|\cdot\|_2} denotes the usual \eqn{\mathbb{L}^2} norm and
#' \eqn{\alpha \ge 0} is the multiview parameter. For the \code{clrv}
#' transformation, refer to [X2Xclrv()].
#'
#' Based on this distance,
#' \ifelse{html}{\out{<i>k</i>}}{\eqn{k}}-means clustering is performed.
#'
#' A reference value \eqn{\alpha_0}, which serves as a baseline around which
#' \eqn{\alpha} may be varied, is selected as follows. The value \eqn{\alpha_0}
#' is defined as the ratio of the total sum of squares of the amplitude
#' components to that of the phase components.
#'
#' See the documentation for [stats::kmeans()] and Kang and Oh (2026) for
#' further details.
#'
#' @param Xclrv A \ifelse{html}{(\out{<i>T</i>}--1)\eqn{\times}\out{<i>n</i>}}{\eqn{(T-1)\times n}} matrix of centered log-ratio
#'   velocity transformed phase components evaluated over the intervals defined
#'   by the time points \code{t}. Refer to [X2Xclrv()].
#' @param Y A \ifelse{html}{\out{<i>T</i>}\eqn{\times}\out{<i>n</i>}}{\eqn{T\times n}} matrix of amplitude components
#'     evaluated at the time points \code{t}.
#' @param t A numeric vector of length \ifelse{html}{\out{<i>T</i>}}{\eqn{T}} giving the time points at which
#'   the phase and amplitude components are evaluated. This vector must start at
#'   0 and end at 1.
#' @param alpha_scale A numeric indicating the value of multiview parameter. The
#'   multiview parameter \eqn{\alpha} is set to \eqn{\alpha = \alpha_0 \times
#'   }\code{alpha_scale}. See the details below. By default, set to 1.
#' @param k A numeric indicating the number of clusters.
#' @param itermax A numeric indicating the maximum number of iterations allowed
#'   in the \ifelse{html}{\out{<i>k</i>}}{\eqn{k}}-means algorithm. By default, set to 10.
#' @param nstart A numeric indicating the number of initial sets. By default, set to 1.
#' @param algorithm A character string indicating the algorithm for \ifelse{html}{\out{<i>k</i>}}{\eqn{k}}-means
#'   clustering. \code{"Hartigan-Wong"} algorithm is set as default.
#' @param trace A boolean. If \code{TRUE} and \code{algorithm ==
#'   "Hartigan-Wong"}, it prints tracing information on the console.
#'
#' @return \code{fkmeans_pre()} and \code{fkmeans()} return an object of class
#'   \code{fkmeans}, which is a list containing the following components:
#'   \item{\code{cluster}}{A vector of integers (from \code{1:k}) indicating the cluster to which each function is allocated.}
#'   \item{\code{centers.Xclrv}}{A \ifelse{html}{(\out{<i>T</i>}--1)\eqn{\times}\out{<i>k</i>}}{\eqn{(T-1)\times k}} matrix of phase components' cluster centers (centered log-ratio velocity transformed). This component is not returned when \code{sync_map == "none"}.}
#'   \item{\code{centers.Y}}{A \ifelse{html}{\out{<i>T</i>}\eqn{\times}\out{<i>k</i>}}{\eqn{T\times k}} matrix of amplitude components' cluster centers. This component is not returned when \code{sync_map == "none"}.}
#'   \item{\code{centers.Ytilde}}{A \ifelse{html}{\out{<i>T</i>}\eqn{\times}\out{<i>k</i>}}{\eqn{T\times k}} matrix of raw functions' cluster centers. This component is only returned when \code{sync_map == "none"}.}
#'   \item{\code{totss}}{The total sum of squares.}
#'   \item{\code{withinss}}{A vector of within-cluster sum of squares, one component per cluster.}
#'   \item{\code{tot.withinss}}{Total within-cluster sum of squares, i.e., \code{sum(withinss)}.}
#'   \item{\code{betweenss}}{The between-cluster sum of squares, i.e., \code{totss - tot.withinss}.}
#'   \item{\code{size}}{The number of functions in each cluster.}
#'   \item{\code{iter}}{The number of (outer) iterations.}
#'   \item{\code{ifault}}{Indicator of a possible algorithm problem; refer to [stats::kmeans()].}
#'   \item{\code{alpha0}}{The reference value \eqn{\alpha_0}.}
#'
#' @references Kang S. and Oh H.-S. (2026) \dQuote{Multiview representation and clustering of
#' functional data,} \emph{Unpublished Manuscript}.
#'
#' @seealso [stats::kmeans()] for multivariate \ifelse{html}{\out{<i>k</i>}}{\eqn{k}}-means clustering.
#'   [fkmedians_pre()] and [fkmedians()] for robust functional \ifelse{html}{\out{<i>k</i>}}{\eqn{k}}-medians
#'   clustering. [auc_sync()] and [fr_sync()] for time-synchronizing mappings.
#'   [X2Xclrv()] for centered log-ratio velocity transformation.
#'
#' @examples
#' t <- seq(0, 1, length.out = 100)
#' sync <- auc_sync(seoul_bike$Ytilde[,1:10], seoul_bike$x, t)
#' fkmeans_pre(X2Xclrv(sync), sync$Y, t, alpha_scale = 1, k = 2)
#' fkmeans(seoul_bike$Ytilde[,1:10], seoul_bike$x, t, sync_map = "auc",
#'   sync_args = 1, alpha_scale = 1, k = 2)
#'
#' @export
fkmeans_pre <- function(Xclrv, Y, t, alpha_scale = 1,
                        k, itermax = 10, nstart = 1,
                        algorithm = c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen"),
                        trace = FALSE){
  if(length(t) <= 1) stop("length of t must be >= 2")
  if(t[1] != 0L || t[length(t)] != 1L) stop("t must be an increasing sequence that starts at zero and end at one")
  if(length(t) != nrow(Y)) stop("length of t must be equivalent to nrow(Y)")
  if(length(t) - 1 != nrow(Xclrv)) stop("length of t - 1 must be equivalent to nrow(X)")
  if(length(alpha_scale) != 1) stop("length of alpha_scale must be 1")
  if(length(k) != 1) warning("only the first element of k is used")
  k <- k[1]
  w.Y <- c((t[2] - t[1]) / 2,
           diff(t, lag = 2) / 2,
           (t[length(t)] - t[length(t) - 1]) / 2)
  sst.X <- sum((Xclrv - rowMeans(Xclrv))^2 * diff(t))
  sst.Y <- sum((Y - rowMeans(Y))^2 * w.Y)
  alpha0 <- sst.Y / sst.X
  DATAMAT <- t(rbind(Xclrv * (alpha0 * alpha_scale)^(1/2) * diff(t),
                     Y * w.Y))
  algorithm <- match.arg(algorithm, c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen"))
  res <- stats::kmeans(x = DATAMAT, centers = k, iter.max = itermax,
                nstart = nstart, algorithm = algorithm, trace = trace)
  centers.Xclrv <- t(res$centers[,seq_along(diff(t))]) / ((alpha0 * alpha_scale)^(1/2) * diff(t))
  centers.Y <- t(res$centers[,length(diff(t)) + seq_along(t)]) / w.Y
  dimnames(centers.Xclrv) <- list(dimnames(Xclrv)[[1L]],1L:k)
  dimnames(centers.Y) <- list(dimnames(Y)[[1L]],1L:k)

  structure(list(cluster = res$cluster, centers.Xclrv = centers.Xclrv, centers.Y = centers.Y,
                 totss = res$totss, withinss = res$withinss, tot.withinss = res$tot.withinss,
                 betweenss = res$betweenss, size = res$size, iter = res$iter, ifault = res$ifault, alpha0 = alpha0),
            class = "fkmeans")
}

fkmeans_raw <- function(Ytilde, x,
                        k, itermax = 10, nstart = 1,
                        algorithm = c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen"),
                        trace = FALSE){
  if(length(x) <= 1) stop("length of x must be >= 2")
  if(length(x) != nrow(Ytilde)) stop("length of x must be equivalent to nrow(Ytilde)")
  if(length(k) != 1) warning("only the first element of k is used")
  k <- k[1]
  w.Y <- c((x[2] - x[1]) / 2,
           diff(x, lag = 2) / 2,
           (x[length(x)] - x[length(x) - 1]) / 2)
  DATAMAT <- t(Ytilde * w.Y)
  algorithm <- match.arg(algorithm, c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen"))
  res <- stats::kmeans(x = DATAMAT, centers = k, iter.max = itermax,
                       nstart = nstart, algorithm = algorithm, trace = trace)
  centers.Ytilde <- t(res$centers) / w.Y
  dimnames(centers.Ytilde) <- list(dimnames(Ytilde)[[1L]],1L:k)

  structure(list(cluster = res$cluster, centers.Ytilde = centers.Ytilde,
                 totss = res$totss, withinss = res$withinss, tot.withinss = res$tot.withinss,
                 betweenss = res$betweenss, size = res$size, iter = res$iter, ifault = res$ifault, alpha0 = NA),
            class = "fkmeans")
}

#' @name fkmeans
#'
#' @inheritParams syncftn
#' @param x A numeric vector of length \ifelse{html}{\out{<i>m</i>}}{\eqn{m}} giving the observed time points
#'   corresponding to \code{Ytilde}.
#' @param sync_map A character string. If \code{"auc"} (the default), AUC
#'   time-synchronizing mapping is used. If \code{"fr"}, FR time-synchronizing
#'   mapping is used. Refer to [auc_sync()] and [fr_sync()]. If \code{"none"},
#'   time-synchronizing mapping is not used, and the functional clustering is conducted on the
#'   observed curves \code{Ytilde}. Hence, for \code{"none"}, the arguments
#'   \code{t}, \code{sync_args}, and \code{alpha_scale} are ignored.
#' @param sync_args If \code{sync_map == "auc"} it represents a
#'   numeric indicating the parameter \ifelse{html}{\out{<i>p</i>}}{\eqn{p}}
#'   used in AUC time-synchronizing mapping. If \code{sync_map == "fr"} it
#'   represent the template function used in FR time-synchronizing mapping.
#' @export
fkmeans <- function(Ytilde, x, t, sync_map = c("auc", "fr", "none"), sync_args,
                    alpha_scale = 1, k, itermax = 10, nstart = 1,
                    algorithm = c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen"),
                    trace = FALSE){
  sync_map <- match.arg(sync_map, c("auc", "fr", "none"))
  if(sync_map == "none"){
    fkmeans_raw(Ytilde = Ytilde, x = x, k = k, itermax = itermax, nstart = nstart,
                algorithm = algorithm, trace = trace)
  }else{
    if(t[1] != 0L || t[length(t)] != 1L) stop("t must be an increasing sequence that starts at zero and end at one")
    if(sync_map == "auc"){
      if(!is.numeric(sync_args) || length(sync_args) != 1) stop("sync_args must be a single numeric; see auc_sync() documentation")
      sync <- auc_sync(Ytilde, x, t, sync_args)
    }else{
      if(length(sync_args) != length(t)) stop("sync_args must be a vector of length = length(t) representing a template function")
      sync <- fr_sync(Ytilde, x, t, sync_args)
    }
    algorithm <- match.arg(algorithm, c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen"))
    fkmeans_pre(Xclrv = X2Xclrv(sync$X, t), Y = sync$Y, t = t, alpha_scale = alpha_scale,
                k, itermax = itermax, nstart = nstart,
                algorithm = algorithm,
                trace = trace)
  }
}

#' @name fkmeans
#' @aliases print.fkmeans
#' @usage NULL
#' @export
print.fkmeans <- function(x, ...){
  cat(
    ngettext(length(x$size),
             paste("Functional k-means clustering with", length(x$size), "cluster of sizes", x$size[1], "\n\n"),
             paste("Functional k-means clustering with", length(x$size), "clusters of sizes",
                   paste(x$size, collapse = ", "), "\n\n"))
  )
  cat("Clustering vector: \n")
  print.default(x$cluster)
  cat("\nWithin cluster sum of squares by cluster: \n")
  print.default(x$withinss)
  cat(" (between_SS / total_SS =", round(x$betweenss / x$totss * 100, 1) ,"%)\n\n")
  cat("Available components: \n")
  print.default(attributes(x)$names)
}

#' @name fkmeans
#' @aliases fitted.fkmeans
#'
#' @param object A \code{fkmeans} object, obtained as a result of the function
#'   \code{fkmeans_pre()} or \code{fkmeans()}.
#' @param method A character string.
#'  - \code{"centers"}: Returns cluster centers for each curve.
#'  - \code{"classes"}: Returns a vector of class assignments.
#' @param ... Not used.
#' @return \code{print()} and \code{fitted()} methods are supported for the
#'   object of class \code{fkmeans}. \code{fitted.fkmeans()} with \code{method =
#'   "centers"} returns cluster centers (one for each input point) and
#'   \code{method = "classes"} returns a vector of class assignments.
#' @export
fitted.fkmeans <- function(object, method = c("centers", "classes"), ...){
  method <- match.arg(method, c("centers", "classes"))
  if (method == "centers")
    list(centers.X = object$centers.X.clrv[, object$cluster, drop = FALSE],
         centers.Y = object$centers.Y[, object$cluster, drop = FALSE])
  else
    object$cluster
}
