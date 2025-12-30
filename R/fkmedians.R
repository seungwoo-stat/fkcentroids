#' @name fkmedians
#' @aliases fkmedians_pre print.fkmedians fitted.fkmedians
#'
#' @title title
#' @export
fkmedians_pre <- function(X.clrv, Y, t, alpha.scale = 1,
                          k, niter = 20, nstart = 1){
  if(length(t) <= 1) stop("length of t must be >= 2")
  if(length(t) != nrow(Y)) stop("length of t must be equal to nrow(Y)")
  if(length(t) - 1 != nrow(X.clrv)) stop("length of t - 1 must be equal to nrow(X)")
  if(length(alpha.scale) != 1) stop("length of alpha.scale must be 1")
  if(length(k) != 1) warning("only the first element of k is used")
  k <- k[1]
  w.Y <- c((t[2] - t[1]) / 2,
           diff(t, lag = 2) / 2,
           (t[length(t)] - t[length(t) - 1]) / 2)
  sst.X <- sum((X.clrv - rowMeans(X.clrv))^2 * diff(t))
  sst.Y <- sum((Y - rowMeans(Y))^2 * w.Y)
  alpha0 <- sst.Y / sst.X
  DATAMAT <- t(rbind(X.clrv * (alpha0 * alpha.scale)^(1/2) * diff(t),
                     Y * w.Y))
  res <- Kmedians::Kmedians(X = DATAMAT, nclust = k, ninit = nstart, niter = niter, method = "Offline", init = TRUE, par = FALSE)
  resb <- res$bestresult

  withinsrs <- sapply(1:k, function(k){
    sum(sqrt(colSums((t(DATAMAT)[,resb$cluster == k,drop = FALSE] - res$bestresult$centers[k,])^2)))
  })
  centers.X.clrv <- t(resb$centers[,seq_along(diff(t))]) / ((alpha0 * alpha.scale)^(1/2) * diff(t))
  centers.Y <- t(resb$centers[,length(diff(t)) + seq_along(t)]) / w.Y
  dimnames(centers.X.clrv) <- list(dimnames(X.clrv)[[1L]],1L:k)
  dimnames(centers.Y) <- list(dimnames(Y)[[1L]],1L:k)
  cluster <- resb$cluster
  if(!is.null(rn <- colnames(Y))) names(cluster) <- rn
  size <- as.numeric(table(factor(cluster, level = 1:k)))

  structure(list(cluster = cluster, centers.X.clrv = centers.X.clrv, centers.Y = centers.Y,
                 withinsrs = withinsrs, tot.withinsrs = sum(withinsrs),
                 size = size, iter = niter, alpha0 = alpha0),
            class = "fkmedians")
}

#' @name fkmedians
#' @title title
#' @export
fkmedians <- function(){

}

#' @name fkmedians
#'
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
#'
#' @param object A \code{fkmedians} object, obtained as a result of the function
#'   \code{fkmedians()}.
#' @param method A character.
#'  - centers: returns cluster centers for each curve.
#'  - classes: returns a vector of class assignments.
#' @return
#' @export
fitted.fkmedians <- function(object, method = c("centers", "classes"), ...){
  method <- match.arg(method)
  if (method == "centers")
    list(centers.X = object$centers.X.clrv[, object$cluster, drop = FALSE],
         centers.Y = object$centers.Y[, object$cluster, drop = FALSE])
  else
    object$cluster
}
