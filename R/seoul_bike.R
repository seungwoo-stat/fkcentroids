#' Seoul Public Bike Rental Records
#'
#' @description
#' Functions generated from Seoul's public bike rental records collected
#' between 00:00 and 24:00 on April 1st, 2025 (Tuesday). Each function represents
#' a single rental station that has at least 24 rentals over the day.
#'
#' Each function is constructed through a two-step process:
#' (1) each rental record is converted to time in hours, that is, a numeric
#' value between 0 and 24; and
#' (2) Gaussian convolution with bandwidth 1 is applied to the rental times
#' to generate a smooth function for each station.
#'
#' @name seoul_bike
#' @docType data
#' @keywords datasets
#' @usage data(seoul_bike)
#' @format A list of length 2 containing the following components:
#' \itemize{
#'   \item{\code{Ytilde}}: A \ifelse{html}{73\eqn{\times}1784}{{\eqn{73\times 1784}}}
#'   matrix whose <i>i</i>th column
#'   contains the values of the <i>i</i>th observed function (bike station)
#'   evaluated at the time points \code{x}.
#'   \item{\code{x}}: A length 73 numeric vector representing time in hours at
#'   20-minute intervals. For example, 18.333 corresponds to 6:20 p.m., i.e.,
#'   \code{seq(0, 1, length.out = 73)}.
#' }
#'
#' @source Seoul Open Data Plaza :
#' https://data.seoul.go.kr/dataList/OA-15182/F/1/datasetView.do
#'
#' @references Kang S. and Oh H.-S. (2026) "Multiview representation and clustering of
#' functional data," \emph{Unpublished Manuscript}.
NULL
