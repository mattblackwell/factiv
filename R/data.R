#' New Haven Voter Modelization Field Experiment
#'
#' A data set that contains a subset of the voter mobilzation field
#' experiment analyzed in Gerber and Green (2000).
#'
#' Data was cleaned and used in Bowers and Hansen (2009). It was a 2x2
#' factorial design with noncompliance on both factors. This is the
#' subset of the subjects who lived in a single-person household and
#' were not randomized to receive get-out-the-vote mailers. Blackwell
#' (2017) analyzed the full data adjusting for noncompliance.
#'
#' @docType data
#' @name newhaven
#' @format  A data frame with 7865 observations and 6 variables:
#' \describe{
#'  \item{ward}{ward of the residence.}
#'  \item{turnout_98}{indicator variable for voting in the 1998
#'   general election.}
#'  \item{inperson}{indicator variable for if the respondent received
#'   in-person canvassing.}
#'  \item{phone}{indicator variable for if the respondent received
#'   phone canvassing.}
#'  \item{inperson_rand}{indicator variable for if the respondent was
#'   randomized to receive in-person canvassing.}
#'  \item{phone_rand}{indicator variable for if the respondent was
#'   randomized to receive phone canvassing.}
#'  \item{age}{respondent age.}
#'  \item{maj_party}{indicator variable for if the respondent is a
#' registered Democrat or Republican.}
#'  \item{turnout_96}{indicator variable for whether or not the
#' respondent voted in the 1996 general election.}
#' }
#' @source https://doi.org/10.7910/DVN/Q6CID7
#' @references Gerber, Alan S., and Donald P. Green. “The Effects of
#'   Canvassing, Telephone Calls, and Direct Mail on Voter Turnout: A
#'   Field Experiment.” The American Political Science Review 94, no.
#'   3 (2000): 653. https://doi.org/10.2307/2585837.
#' 
#' Hansen, Ben B., and Jake Bowers. “Attributing Effects to a
#'   Cluster-Randomized Get-Out-the-Vote Campaign.” Journal of the
#'   American Statistical Association 104, no. 487 (2009):
#'   873–85. https://doi.org/10.1198/jasa.2009.ap06589.
#'
#' Blackwell, Matthew. “Instrumental Variable Methods for Conditional
#'   Effects and Causal Interaction in Voter Mobilization
#'   Experiments.” Journal of the American Statistical Association
#'   112, no. 518 (2017): 590–99.
#'   https://doi.org/10.1080/01621459.2016.1246363.
#' 
NULL
