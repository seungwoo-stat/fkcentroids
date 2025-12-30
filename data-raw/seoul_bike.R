bike <- readRDS("data-raw/bike250401.rds")
bike.pattern <- by(bike$rent.time, bike$대여.대여소번호, identity)
bike.pattern <- bike.pattern[sapply(bike.pattern, length) >= 24]

x <- seq(0,24,length.out = 73) # every 20 minutes

bw <- 1
bike.pattern.kernel <- bike.pattern
for(i in seq_along(bike.pattern.kernel)){
  bike.pattern.kernel[[i]] <- sapply(x, \(tt) dnorm(tt, bike.pattern[[i]], bw) |> sum())
}

Y.tilde <- matrix(unlist(bike.pattern.kernel), nrow = length(x))
matplot(x, Y.tilde, type = "l", col = 1)
colnames(Y.tilde) <- names(bike.pattern)

seoul_bike <- list(Ytilde = Y.tilde, x = x)
usethis::use_data(seoul_bike, overwrite = TRUE)
