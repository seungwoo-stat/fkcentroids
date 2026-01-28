#' @name fkmedians
#' @aliases fkmedians_pre print.fkmedians fitted.fkmedians
#' @title Functional \eqn{k}-Medians Clustering Using Phase and Amplitude
#'   Components
#'
#' @description
#' Conducts \ifelse{html}{\out{<i>k</i>}}{\eqn{k}}-medians clustering by jointly considering phase and amplitude
#' variation. The relative importance of the two components can be explicitly
#' controlled by the user via the multiview parameter \eqn{\alpha}. Optionally,
#' \ifelse{html}{\out{<i>k</i>}}{\eqn{k}}-medians clustering can be performed
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
#' \ifelse{html}{\out{<i>k</i>}}{\eqn{k}}-medians clustering is performed. In
#' particular, Weiszfeld algorithm is used to find the geometric median function
#' for each cluster, implemented in [Kmedians::Kmedians()].
#'
#' A reference value \eqn{\alpha_0}, which serves as a baseline around which
#' \eqn{\alpha} may be varied, is selected as follows. The value \eqn{\alpha_0}
#' is defined as the ratio of the total sum of squares of the amplitude
#' components to that of the phase components.
#'
#' See the documentation for [Kmedians::Kmedians()] (Godichon-Baggioni and
#' Surendran, 2024) and Kang and Oh (2026) for further details.
#'
#' @inheritParams fkmeans
#' @param niter A numeric indicating the number of iterations for
#'   the \ifelse{html}{\out{<i>k</i>}}{\eqn{k}}-medians algorithm. By default, set to 20.
#'
#' @return \code{fkmedians_pre()} and \code{fkmedians()} return an object of class
#'   \code{fkmedians}, which is a list containing the following components:
#'   \item{\code{cluster}}{A vector of integers (from \code{1:k}) indicating the cluster to which each function is allocated.}
#'   \item{\code{centers.Xclrv}}{A \ifelse{html}{(\out{<i>T</i>}--1)\eqn{\times}\out{<i>k</i>}}{\eqn{(T-1)\times k}} matrix of phase components' cluster centers (centered log-ratio velocity transformed). This component is not returned when \code{sync_map == "none"}.}
#'   \item{\code{centers.Y}}{A \ifelse{html}{\out{<i>T</i>}\eqn{\times}\out{<i>k</i>}}{\eqn{T\times k}} matrix of amplitude components' cluster centers. This component is not returned when \code{sync_map == "none"}.}
#'   \item{\code{centers.Ytilde}}{A \ifelse{html}{\out{<i>T</i>}\eqn{\times}\out{<i>k</i>}}{\eqn{T\times k}} matrix of raw functions' cluster centers. This component is only returned when \code{sync_map == "none"}.}
#'   \item{\code{withinsrs}}{A vector of within-cluster sum of residuals, one component per cluster.}
#'   \item{\code{tot.withinsrs}}{Total within-cluster sum of residuals, i.e., \code{sum(withinsrs)}.}
#'   \item{\code{size}}{The number of functions in each cluster.}
#'   \item{\code{iter}}{The number of (outer) iterations.}
#'   \item{\code{alpha0}}{The reference value \eqn{\alpha_0}. This component is not returned when \code{sync_map == "none"}.}
#'
#' @references Godichon-Baggioni A. and Surendran S. (2024) \dQuote{A
#'   penalized criterion for selecting the number of clusters for K-medians,}
#'   \emph{Journal of Computational and Graphical Statistics}, \strong{33}(4),
#'   1298--1309.
#'
#' Kang S. and Oh H.-S. (2026) \dQuote{Multiview representation and clustering of
#' functional data,} \emph{Unpublished Manuscript}.
#'
#' @seealso [Kmedians::Kmedians()] for multivariate \ifelse{html}{\out{<i>k</i>}}{\eqn{k}}-medians clustering.
#'   [fkmeans_pre()] and [fkmeans()] for functional \ifelse{html}{\out{<i>k</i>}}{\eqn{k}}-means
#'   clustering. [auc_sync()] and [fr_sync()] for time-synchronizing mappings.
#'   [X2Xclrv()] for centered log-ratio velocity transformation.
#'
#' @examples
#' t <- seq(0, 1, length.out = 100)
#' sync <- auc_sync(seoul_bike$Ytilde[,1:10], seoul_bike$x, t)
#' fkmedians_pre(X2Xclrv(sync), sync$Y, t, alpha_scale = 1, k = 2, nstart = 10)
#' fkmedians(seoul_bike$Ytilde[,1:10], seoul_bike$x, t, sync_map = "auc",
#'   sync_args = 1, alpha_scale = 1, k = 2, nstart = 10)
#'
#' @export
fkmedians_pre <- function(Xclrv, Y, t, alpha_scale = 1,
                          k, niter = 20, nstart = 1){
  if(length(t) <= 1) stop("length of t must be >= 2")
  if(t[1] != 0L || t[length(t)] != 1L) stop("t must be an increasing sequence that starts at zero and end at one")
  if(length(t) != nrow(Y)) stop("length of t must be equal to nrow(Y)")
  if(length(t) - 1 != nrow(Xclrv)) stop("length of t - 1 must be equal to nrow(X)")
  if(length(alpha_scale) != 1) stop("length of alpha_scale must be 1")
  if(length(k) != 1) warning("only the first element of k is used")
  k <- k[1]
  w.Y <- c((t[2] - t[1]) / 2,
           diff(t, lag = 2) / 2,
           (t[length(t)] - t[length(t) - 1]) / 2)
  sst.X <- sum((Xclrv - rowMeans(Xclrv))^2 * diff(t))
  sst.Y <- sum((Y - rowMeans(Y))^2 * w.Y)
  alpha0 <- sst.Y / sst.X
  DATAMAT <- t(rbind(Xclrv * (alpha0 * alpha_scale)^(1/2) * sqrt(diff(t)),
                     Y * sqrt(w.Y)))
  res <- Kmedians::Kmedians(X = DATAMAT, nclust = k, ninit = nstart, niter = niter, method = "Offline", init = TRUE, par = FALSE)
  resb <- res$bestresult

  withinsrs <- sapply(1:k, function(k){
    sum(sqrt(colSums((t(DATAMAT)[,resb$cluster == k,drop = FALSE] - res$bestresult$centers[k,])^2)))
  })
  centers.Xclrv <- t(resb$centers[,seq_along(diff(t))]) / ((alpha0 * alpha_scale)^(1/2) * sqrt(diff(t)))
  centers.Y <- t(resb$centers[,length(diff(t)) + seq_along(t)]) / sqrt(w.Y)
  dimnames(centers.Xclrv) <- list(dimnames(Xclrv)[[1L]],1L:k)
  dimnames(centers.Y) <- list(dimnames(Y)[[1L]],1L:k)
  cluster <- resb$cluster
  if(!is.null(rn <- colnames(Y))) names(cluster) <- rn
  size <- as.numeric(table(factor(cluster, levels = 1:k)))

  structure(list(cluster = cluster, centers.Xclrv = centers.Xclrv, centers.Y = centers.Y,
                 withinsrs = withinsrs, tot.withinsrs = sum(withinsrs),
                 size = size, iter = niter, alpha0 = alpha0),
            class = "fkmedians")
}

fkmedians_raw <- function(Ytilde, x, k, niter = 20, nstart = 1){
  if(length(x) <= 1) stop("length of x must be >= 2")
  if(length(x) != nrow(Ytilde)) stop("length of x must be equivalent to nrow(Ytilde)")
  if(length(k) != 1) warning("only the first element of k is used")
  k <- k[1]
  w.Y <- c((x[2] - x[1]) / 2,
           diff(x, lag = 2) / 2,
           (x[length(x)] - x[length(x) - 1]) / 2)
  DATAMAT <- t(Ytilde * sqrt(w.Y))
  res <- Kmedians::Kmedians(X = DATAMAT, nclust = k, ninit = nstart, niter = niter, method = "Offline", init = TRUE, par = FALSE)
  resb <- res$bestresult

  withinsrs <- sapply(1:k, function(k){
    sum(sqrt(colSums((t(DATAMAT)[,resb$cluster == k,drop = FALSE] - res$bestresult$centers[k,])^2)))
  })
  centers.Ytilde <- t(resb$centers) / sqrt(w.Y)
  dimnames(centers.Ytilde) <- list(dimnames(Ytilde)[[1L]],1L:k)
  cluster <- resb$cluster
  if(!is.null(rn <- colnames(Ytilde))) names(cluster) <- rn
  size <- as.numeric(table(factor(cluster, levels = 1:k)))

  structure(list(cluster = cluster, centers.Ytilde = centers.Ytilde,
                 withinsrs = withinsrs, tot.withinsrs = sum(withinsrs),
                 size = size, iter = niter),
            class = "fkmedians")
}

#' @name fkmedians
#'
#' @inheritParams fkmeans
#' @export
fkmedians <- function(Ytilde, x, t, sync_map = c("auc", "fr", "none"), sync_args = NULL,
                      alpha_scale = 1, k, niter = 20, nstart = 1){
  sync_map <- match.arg(sync_map, c("auc", "fr", "none"))
  if(sync_map == "none"){
    fkmedians_raw(Ytilde = Ytilde, x = x, k = k, niter = niter, nstart = nstart)
  }else{
    if(t[1] != 0L || t[length(t)] != 1L) stop("t must be an increasing sequence that starts at zero and end at one")
    if(sync_map == "auc"){
      if(!is.numeric(sync_args) || length(sync_args) != 1) stop("sync_args must be a single numeric; see auc_sync() documentation")
      sync <- auc_sync(Ytilde, x, t, sync_args)
    }else{
      if(length(sync_args) != length(t)) stop("sync_args must be a vector of length = length(t) representing a template function")
      sync <- fr_sync(Ytilde, x, t, sync_args)
    }
    fkmedians_pre(Xclrv = X2Xclrv(sync$X, t), Y = sync$Y, t = t, alpha_scale = alpha_scale,
                  k, niter = niter, nstart = nstart)
  }
}

#' @name fkmedians
#' @aliases print.fkmedians
#' @usage NULL
#' @export
print.fkmedians <- function(x, ...){
  cat(
    ngettext(length(x$size),
             paste("Functional k-medians clustering with", length(x$size), "cluster of sizes", x$size[1], "\n\n"),
             paste("Functional k-medians clustering with", length(x$size), "clusters of sizes",
                   paste(x$size, collapse = ", "), "\n\n"))
  )
  cat("Clustering vector: \n")
  print.default(x$cluster)
  cat("\nWithin cluster sum of residuals by cluster: \n")
  print.default(x$withinsrs)
  # cat(" (between_SS / total_SS =", round(x$betweenss / x$totss * 100, 1) ,"%)\n\n")
  cat("\nAvailable components: \n")
  print.default(attributes(x)$names)
}



#' @name fkmedians
#' @aliases fitted.fkmedians
#'
#' @param object A \code{fkmedians} object, obtained as a result of the function
#'   \code{fkmedians_pre()} or \code{fkmedians()}.
#' @return \code{print()} and \code{fitted()} methods are supported for the
#'   object of class \code{fkmedians}. \code{fitted.fkmedians()} with \code{method =
#'   "centers"} returns cluster centers (one for each input point) and
#'   \code{method = "classes"} returns a vector of class assignments.
#' @export
fitted.fkmedians <- function(object, method = c("centers", "classes"), ...){
  method <- match.arg(method, c("centers", "classes"))
  if (method == "centers")
    list(centers.X = object$centers.Xclrv[, object$cluster, drop = FALSE],
         centers.Y = object$centers.Y[, object$cluster, drop = FALSE])
  else
    object$cluster
}
