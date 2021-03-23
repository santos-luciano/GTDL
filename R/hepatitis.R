#' Hepatitis data
#'
#' @description A randomized clinical trial was carried out to investigate the effect of therapy with steroid in the treatment of acute
#' viral hepatitis (Gregory et al., 1975). Twenty and nine patients with this disease were randomized to receive a placebo or the steroid
#' treatment. Each patient was followed up for 16 weeks or until death( event of interest) or even loss of follow-up.
#'
#' @docType data
#' @keywords dataSets
#' @name hepatitis
#' @usage data(hepatitis)
#' @format This data frame contains the following columns:
#' \itemize{
#' \item {t:}{survival time in weeks}
#' \item {censured:} {censored = 0, dead = 1}
#' \item {group:} {control = 0, steroid = 1}
#' }
#' @references
#' 
#'Colosimo, E. A and  Giolo, S. R. Análise de sobrevivência aplicada.  Edgard Blucher: São Paulo. 2006.
#'
#' @examples
#'
#' data(hepatitis)
#' head(hepatitis)
#'