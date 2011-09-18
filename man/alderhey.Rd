\name{alderhey}
\alias{alderhey}
\docType{data}
\title{Meningococcal disease data}
\description{
  This dataset consists of 6 list objects. The \code{cases} dataframe consists of 864 rows and 16 columns 
  of data on all 0-14 year olds diagnosed with meningococcal disease
  in the Merseyside region 1981 - 2007. The \code{control} data has been obtained from the Office of National 
  Statistics for the same time period and region, and provides data at the Output Area (OA) level over time. \code{cent}
  is a dataframe containing 4586 rows and 3 columns which report the coordinates of the centroids of each OA.
}
\usage{data(alderhey)}
\format{
  This data object contains a list of 6 objects:

   \code{$cases}: A dataframe comprised of 864 observations and  16 variables:
    \tabular{lll}{
    \tab \code{OA} \tab Code of the OA within which the case resides \cr
    \tab \code{id} \tab Patient ID \cr
    \tab \code{x}  \tab x-coordinate of the residence of the patient \cr
    \tab \code{y} \tab y-coordinate of the residence of the patient \cr
    \tab \code{age} \tab age of the patient to two decimal places \cr
    \tab \code{sex} \tab indicator variable for sex, with 1=Male, 0=Female \cr
    \tab \code{townsend} \tab Townsend score of the OA of residence during the year of diagnosis \cr
    \tab \code{day}  \tab Day of diagnosis \cr
    \tab \code{month} \tab Month of diagnosis \cr
    \tab \code{year} \tab Year of diagnosis \cr
    \tab \code{LA} \tab Local Authority within which the patient resides, with 1=Knowsley, \cr
    \tab           \tab     2=Liverpool,3=St Helens, 4=Sefton,5=Wirral\cr
    \tab \code{ind1} \tab Indicator function = 1 if patient resides in Knowsley, 0 otherwise \cr
    \tab \code{ind2} \tab Indicator function = 1 if patient resides in Liverpool, 0 otherwise \cr
    \tab \code{ind3} \tab Indicator function = 1 if patient resides in St Helens, 0 otherwise \cr
    \tab \code{ind4} \tab Indicator function = 1 if patient resides in Sefton, 0 otherwise \cr
    \tab \code{ind5} \tab Indicator function = 1 if patient resides in the Wirral, 0 otherwise
    }
 \code{$control}: A list of 11 objects
  \tabular{lll}{
  \tab \code{OAsum} \tab The population of each OA summed over time \cr
  \tab \code{OAid} \tab ID variable for each OA \cr
  \tab \code{sumall} \tab A 4586 by 27 matrix of population counts \cr
  \tab \code{countLA1} \tab A 4586 by 27 matrix of Knowsley population counts \cr
  \tab \code{countLA2} \tab A 4586 by 27 matrix of Liverpool population counts \cr
  \tab \code{countLA3} \tab A 4586 by 27 matrix of St Helens population counts \cr
  \tab \code{countLA4} \tab A 4586 by 27 matrix of Sefton population counts \cr
  \tab \code{countLA5} \tab A 4586 by 27 matrix of Wirral population counts \cr
  \tab \code{sumage} \tab A 4586 by 27 matrix of the sum of the ages of the 0-14 year old population \cr
  \tab \code{sumgender} \tab A 4586 by 27 matrix of the number of males for each OA at each year \cr
  \tab \code{Townsend} \tab A 4586 by 27 matrix of the Townsend score for each OA at each year
  }

 \code{$cent}: A data frame with 4586 rows and 3 columns:
  \tabular{lll}{
  \tab \code{OAid} \tab ID variable for each OA \cr 
  \tab \code{xcoord} \tab x-coordinate of the centroid of the OA \cr
  \tab \code{ycoord} \tab y-coordinate of the centroid of the OA
  }
}

\details{
  Note that in the \code{control} list the LA-specific population counts equal \code{0} if the OA is not part of the respective LA}

\source{
  Population estimates for intercensus years at the Local Authority level were obtained from \url{http://www.statistics.gov.uk/StatBase/Product.asp?vlnk=601}.

  Population counts for census years at the Output Area level were obtained from \url{http://casweb.mimas.ac.uk/} 
}
\references{
  Data relating to cases of meningococcal disease in Merseyside were obtained from Alder Hey Children's Hospital, Liverpool
  and the Health Protection Agency. See \cr \cr
  Stanton, M.C., Taylor-Robinson, D., Harris, D., Paize, F., Makwana, N., Hackett, S.J., Baines, P.B., Riordan, F.A.I., Marzouk, O., Thomson, A.P.J., Diggle, P., Hart, C.A and Carrol, E.A. (2011) Meningococcal disease in children in Merseyside, England: a 31 year descriptive study. \emph{Accepted, PLoS One}
}

\keyword{datasets}
