#############################################################################
#    Copyright (c) 2013-2018 Alexandra Kuznetsova, Per Bruun Brockhoff, and
#    Rune Haubo Bojesen Christensen
#
#    This file is part of the lmerTest package for R (*lmerTest*)
#
#    *lmerTest* is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#
#    *lmerTest* is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    A copy of the GNU General Public License is available at
#    <https://www.r-project.org/Licenses/> and/or
#    <http://www.gnu.org/licenses/>.
#############################################################################
#
# data_documentation.R - roxygen2 documentation for datasets.

# Datasets documented in this file:
#
# - carrots
# - ham
# - TVbo

##############################################
######## carrots
##############################################
#' Consumer Preference Mapping of Carrots
#'
#' In a consumer study 103 consumers scored their preference of 12 danish
#' carrot types on a scale from 1 to 7. Moreover the consumers scored the
#' degree of sweetness, bitterness and crispiness in the products.
#'
#' The carrots were harvested in autumn 1996 and tested in march 1997. In
#' addition to the consumer survey, the carrot products were evaluated by
#' a trained panel of tasters, the sensory panel, with respect to a
#' number of sensory (taste, odour and texture) properties. Since usually
#' a high number of (correlated) properties (variables) are used, in this
#' case 14, it is a common procedure to use a few, often 2, combined
#' variables that contain as much of the information in the sensory
#' variables as possible. This is achieved by extracting the first two
#' principal components in a principal components analysis (PCA) on the
#' product-by-property panel average data matrix. In this data set the
#' variables for the first two principal components are named
#' (\code{sens1} and \code{sens2}).
#'
#' @docType data
#'
#' @usage data(carrots)
#'
#' @format
#' \describe{
#' \item{Consumer}{factor with 103 levels: numbering identifying consumers.}
#' \item{Frequency}{factor with 5 levels; "How often do you eat carrots?"
#' 1: once a week or more, 2: once
#' every two weeks, 3: once every three weeks, 4: at least once month,
#' 5: less than once a month.}
#' \item{Gender}{factor with 2 levels. 1: male, 2:female.}
#' \item{Age}{factor with 4 levels. 1: less than 25 years, 2: 26-40 years,
#' 3: 41-60 years, 4 more than 61 years.}
#' \item{Homesize}{factor with two levels. Number of persons in the household.
#' 1: 1 or 2 persons, 2: 3 or more persons.}
#' \item{Work}{factor with 7 levels. different types of employment.
#' 1: unskilled worker(no education),
#' 2: skilled worker(with education), 3: office worker, 4: housewife (or man),
#' 5: independent
#' businessman/ self-employment, 6: student, 7: retired}
#' \item{Income}{factor with 4 levels. 1: <150000, 2: 150000-300000,
#' 3: 300000-500000, 4: >500000}
#' \item{Preference}{consumer score on a seven-point scale.}
#' \item{Sweetness}{consumer score on a seven-point scale.}
#' \item{Bitterness}{consumer score on a seven-point scale.}
#' \item{Crispness}{consumer score on a seven-point scale.}
#' \item{sens1}{first sensory variable derived from a PCA.}
#' \item{sens2}{second sensory variable derived from a PCA.}
#' \item{Product}{factor on 12 levels.}
#' }
#'
#' @keywords datasets
#' @source Per Bruun Brockhoff, The Royal Veterinary and Agricultural University,
#' Denmark.
#'
#' @examples
#'
#' fm <- lmer(Preference ~ sens2 + Homesize + (1 + sens2 | Consumer), data=carrots)
#' anova(fm)
#'
"carrots"


##############################################
######## ham
##############################################
#' Conjoint Study of Dry Cured Ham
#'
#' One of the purposes of the study was to investigate the effect of
#' information given to the consumers measured in hedonic liking for the
#' hams. Two of the hams were Spanish and two were Norwegian, each origin
#' representing different salt levels and different aging time. The
#' information about origin was given in such way that both true and
#' false information was given. Essentially a 4x2 design with 4 samples
#' and 2 information levels. A total of 81 Consumers participated in the
#' study.
#'
#' @docType data
#'
#' @usage data(ham)
#'
#' @format
#' \describe{
#' \item{Consumer}{factor with 81 levels: numbering identifying consumers.}
#' \item{Product}{factor with four levels.}
#' \item{Informed.liking}{numeric: hedonic liking for the products.}
#' \item{Information}{factor with two levels.}
#' \item{Gender}{factor with two levels.}
#' \item{Age}{numeric: age of Consumer.}
#' }
#'
#' @keywords datasets
#'
#' @references
#' T. Næs, V. Lengard, S. Bølling Johansen, M. Hersleth (2010)
#' Alternative methods for combining design variables and consumer preference
#' with information about attitudes and demographics in conjoint analysis,
#' \emph{Food Quality and Preference}, 10-4, 368-378, ISSN 0950-3293,
#' \url{https://doi.org/10.1016/j.foodqual.2009.09.004}.
#'
#' @examples
#'
#' # Simple model for the ham data:
#' fm <- lmer(Informed.liking ~ Product*Information + (1|Consumer) , data=ham)
#'
#' # Anova table for the fixed effects:
#' anova(fm)
#'
#' \dontrun{
#' # Fit 'big' model:
#' fm <- lmer(Informed.liking ~ Product*Information*Gender*Age +
#'              + (1|Consumer) + (1|Consumer:Product) +
#'              (1|Consumer:Information),
#'            data=ham)
#' step_fm <- step(fm)
#' step_fm # Display elimination results
#' final_fm <- get_model(step_fm)
#' }
#'
"ham"


##############################################
######## TVbo
##############################################
#' Sensory Assesment of B&O TVs
#'
#' The TVbo dataset has kindly been made available by the Danish high-end
#' consumer electronics company
#' \href{https://www.bang-olufsen.com}{Bang & Olufsen}.
#' The main purpose was to assess 12 different TV sets (products) specified by
#' the two attributes Picture and TVset.
#' 15 different response variables (characteristics of the
#' product) were assessed by a trained panel with 8 assessors.
#'
#' @format
#' \describe{
#' \item{Assessor}{factor with 8 levels assessors.}
#' \item{TVset}{product factor with 3 levels.}
#' \item{Picture}{product factor with 4 levels.}
#' }
#' In addition the following 15 numeric (response) variables are the
#' characteristics on which the TV sets (products) are assessed:
#'
#' Coloursaturation, Colourbalance, Noise, Depth, Sharpness, Lightlevel,
#' Contrast, Sharpnessofmovement, Flickeringstationary, Flickeringmovement,
#' Distortion, Dimglasseffect, Cutting, Flossyedges, Elasticeffect.
#'
#' @docType data
#'
#' @usage data(TVbo)
#'
#' @examples
#'
#' fm <- lmer(Coloursaturation ~ TVset + Picture + (1|Assessor:TVset) +
#'              (1|Assessor), data=TVbo)
#' ranova(fm)
#' anova(fm)
#'
"TVbo"
