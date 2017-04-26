
<!-- README.md is generated from README.Rmd. Please edit that file -->
Inspect which packages use licenses recommended for data and see if they are, indeed, enriched for data packages, whatever that means.

<https://choosealicense.com/non-software/>

``` r
library(jsonlite)
library(tidyverse)
#> + ggplot2 2.2.1             Date: 2017-04-25
#> + tibble  1.3.0                R: 3.3.2
#> + tidyr   0.6.1.9000          OS: OS X El Capitan 10.11.6
#> + readr   1.1.0              GUI: X11
#> + purrr   0.2.2.9000      Locale: en_CA.UTF-8
#> + dplyr   0.5.0.9004          TZ: America/Vancouver
#> + stringr 1.2.0           
#> + forcats 0.2.0
#> Conflicts -----------------------------------------------------------------
#> * filter(),  from dplyr, masks stats::filter()
#> * lag(),     from dplyr, masks stats::lag()
library(forcats)
library(httr)

# dump <- GET("http://crandb.r-pkg.org/-/latest")
# ct <- content(dump)
# 
# ct %>%
#   toJSON() %>%
#   write_file("cran-packages.json")

ct <- fromJSON("cran-packages.json")

str(ct[[1]])
#> List of 18
#>  $ Package         : chr "A3"
#>  $ Type            : chr "Package"
#>  $ Title           : chr "Accurate, Adaptable, and Accessible Error Metrics for Predictive\nModels"
#>  $ Version         : chr "1.0.0"
#>  $ Date            : chr "2015-08-15"
#>  $ Author          : chr "Scott Fortmann-Roe"
#>  $ Maintainer      : chr "Scott Fortmann-Roe <scottfr@berkeley.edu>"
#>  $ Description     : chr "Supplies tools for tabulating and analyzing the results of predictive models. The methods employed are applicable to virtually "| __truncated__
#>  $ License         : chr "GPL (>= 2)"
#>  $ Depends         :List of 3
#>   ..$ R      : chr ">= 2.15.0"
#>   ..$ xtable : chr "*"
#>   ..$ pbapply: chr "*"
#>  $ Suggests        :List of 2
#>   ..$ randomForest: chr "*"
#>   ..$ e1071       : chr "*"
#>  $ NeedsCompilation: chr "no"
#>  $ Packaged        : chr "2015-08-16 14:17:33 UTC; scott"
#>  $ Repository      : chr "CRAN"
#>  $ Date/Publication: chr "2015-08-16 23:05:52"
#>  $ crandb_file_date: chr "2015-08-16 17:08:22"
#>  $ date            : chr "2015-08-16T23:05:52+00:00"
#>  $ releases        : list()
length(ct)
#> [1] 10476

## CC licenses interest me and I'm curious about 'Unlimited""
## assuming that a data package unlikely to need compilation
dp <- keep(
  ct,
  ~ grepl("Unlimited|CC", .x$License) && .x$NeedsCompilation == "no"
)
length(dp)
#> [1] 150

df <- tibble(
  license = map_chr(dp, "License"),
  name = map_chr(dp, "Package"),
  title = map_chr(dp, "Title"),
  desc = map_chr(dp, "Description"),
  url = map(dp, "URL"),
  pkg = dp
) %>% 
  mutate(license = fct_infreq(license)) %>% 
  arrange(license, name)
      
df %>% 
  count(license, sort = TRUE)
#> # A tibble: 11 × 2
#>                           license     n
#>                            <fctr> <int>
#> 1                             CC0    92
#> 2                       Unlimited    32
#> 3                 CC BY-NC-SA 4.0    11
#> 4                    CC BY-SA 4.0     7
#> 5                       CC BY 4.0     2
#> 6                    CC BY-NC 4.0     1
#> 7                 CC BY-NC-SA 3.0     1
#> 8     CC BY-SA 2.0 + file LICENSE     1
#> 9     CC BY-SA 4.0 + file LICENSE     1
#> 10 MIT + file LICENSE | Unlimited     1
#> 11       Unlimited | file LICENSE     1
```

``` r
f <- function(..., sep = "  \n") cat(paste(..., sep = sep), "\n\n")
df %>% 
  select(-pkg) %>% 
  mutate(license = as.character(license)) %>% 
  pwalk(f)
```

CC0
ameco
European Commission Annual Macro-Economic (AMECO) Database
Annual macro-economic database provided by the European Commission.
<http://github.com/expersso/ameco>

CC0
AntWeb
programmatic interface to the AntWeb
A complete programmatic interface to the AntWeb database from the<U+000a>California Academy of Sciences.
<https://github.com/ropensci/AntWeb>

CC0
aop
Adverse Outcome Pathway Analysis
Provides tools for analyzing adverse outcome pathways (AOPs) for pharmacological and toxicological research. Functionality includes the ability to perform causal network analysis of networks developed in and exported from Cytoscape or existing as R graph objects, and identifying the point of departure/screening/risk value from concentration- response data.

CC0
asaur
Data Sets for "Applied Survival Analysis Using R""
Data sets are referred to in the text "Applied Survival Analysis Using R" by Dirk F. Moore, Springer, 2016, ISBN: 978-3-319-31243-9, <DOI:10.1007/978-3-319-31245-3>.

CC0
averisk
Calculation of Average Population Attributable Fractions and Confidence Intervals
Average population attributable fractions are calculated for a set of risk factors (either binary or ordinal valued) for both prospective and case- control designs. Confidence intervals are found by Monte Carlo simulation. The method can be applied to either prospective or case control designs, provided an estimate of disease prevalence is provided. In addition to an exact calculation of AF, an approximate calculation, based on randomly sampling permutations has been implemented to ensure the calculation is computationally tractable when the number of risk factors is large.

CC0
babynames
US Baby Names 1880-2015
US baby names provided by the SSA. This package contains all names used for at least 5 children of either sex.
<http://github.com/hadley/babynames>

CC0
banxicoR
Download Data from the Bank of Mexico
Provides functions to scrape IQY calls to Bank of Mexico, downloading and ordering the data conveniently.

CC0
bea.R
Bureau of Economic Analysis API
Provides an R interface for the Bureau of Economic Analysis (BEA) API (see <http://www.bea.gov/API/bea_web_service_api_user_guide.htm> for more information) that serves two core purposes - 1. To Extract/Transform/Load data \[beaGet()\] from the BEA API as R-friendly formats in the user's work space \[transformation done by default in beaGet() can be modified using optional parameters; see, too, bea2List(), bea2Tab()\]. 2. To enable the search of descriptive meta data \[beaSearch()\]. Other features of the library exist mainly as intermediate methods or are in early stages of development. Important Note - You must have an API key to use this library. Register for a key at <http://www.bea.gov/API/signup/index.cfm> .
<https://github.com/us-bea/beaR>

CC0
bgmfiles
Example BGM Files for the Atlantis Ecosystem Model
A collection of box-geometry model (BGM) files for the Atlantis ecosystem model. Atlantis is a deterministic, biogeochemical, whole-of-ecosystem model (see <http://atlantis.cmar.csiro.au/> for more information).
<https://github.com/AustralianAntarcticDivision/bgmfiles/>

CC0
bikeshare14
Bay Area Bike Share Trips in 2014
Anonymised Bay Area bike share trip data for the year 2014. Also contains additional metadata on stations and weather.
<http://github.com/arunsrinivasan/bikeshare14>

CC0
bodenmiller
Profilling of Peripheral Blood Mononuclear Cells using CyTOF
This data package contains a subset of the Bodenmiller et al, Nat Biotech 2012 dataset for testing single cell, high dimensional analysis and visualization methods.
<https://github.com/yannabraham/bodenmiller>

CC0
chromer
Interface to Chromosome Counts Database API
A programmatic interface to the Chromosome Counts Database (<http://ccdb.tau.ac.il/>). This package is part of the rOpenSci suite (<http://ropensci.org>)
<http://www.github.com/ropensci/chromer>

CC0
corlink
Record Linkage, Incorporating Imputation for Missing Agreement Patterns, and Modeling Correlation Patterns Between Fields
A matrix of agreement patterns and counts for record pairs is the input for the procedure. An EM algorithm is used to impute plausible values for missing record pairs. A second EM algorithm, incorporating possible correlations between per-field agreement, is used to estimate posterior probabilities that each pair is a true match - i.e. constitutes the same individual.

CC0
csp
Correlates of State Policy Data Set in R
Provides the Correlates of State Policy data set for easy use in R.
<http://ippsr.msu.edu/public-policy/correlates-state-policy> <https://github.com/expersso/csp>

CC0
dam
Data Analysis Metabolomics
A collection of functions which aim to assist common computational workflow for analysis of matabolomic data..

CC0
dataRetrieval
Retrieval Functions for USGS and EPA Hydrologic and Water Quality Data
Collection of functions to help retrieve U.S. Geological Survey (USGS) and U.S. Environmental Protection Agency (EPA) water quality and hydrology data from web services. USGS web services are discovered from National Water Information System (NWIS) tools. Both EPA and USGS water quality data are obtained from the Water Quality Portal <https://www.waterqualitydata.us/>.
<https://pubs.er.usgs.gov/publication/tm4A10>

CC0
datastepr
An Implementation of a SAS-Style Data Step
Based on a SAS data step. This allows for row-wise dynamic building of data, iteratively importing slices of existing dataframes, conducting analyses, and exporting to a results frame. This is particularly useful for differential or time-series analyses, which are often not well suited to vector- based operations.
<https://github.com/bramtayl/datastepr>

CC0
EBrank
Empirical Bayes Ranking
Empirical Bayes ranking applicable to parallel-estimation settings where the estimated parameters are asymptotically unbiased and normal, with known standard errors. A mixture normal prior for each parameter is estimated using Empirical Bayes methods, subsequentially ranks for each parameter are simulated from the resulting joint posterior over all parameters (The marginal posterior densities for each parameter are assumed independent). Finally, experiments are ordered by expected posterior rank, although computations minimizing other plausible rank-loss functions are also given.

CC0
ecb
Programmatic Access to the European Central Bank's Statistical Data Warehouse (SDW)
Provides an interface to the European Central Bank's Statistical Data Warehouse API, allowing for programmatic retrieval of a vast quantity of statistical data.
<https://sdw.ecb.europa.eu/> <https://www.ecb.int>

CC0
ecospace
Simulating Community Assembly and Ecological Diversification Using Ecospace Frameworks
Implements stochastic simulations of community assembly (ecological diversification) using customizable ecospace frameworks (functional trait spaces). Provides a wrapper to calculate common ecological disparity and functional ecology statistical dynamics as a function of species richness. Functions are written so they will work in a parallel-computing environment.
<https://github.com/pnovack-gottshall/ecospace>, <http://www.ben.edu/faculty/pnovack-gottshall/index.html>

CC0
EGRET
Exploration and Graphics for RivEr Trends (EGRET)
Statistics and graphics for streamflow history, water quality trends, and the statistical modeling algorithm: Weighted Regressions on Time, Discharge, and Season (WRTDS).
<http://pubs.usgs.gov/tm/04/a10/>, <https://github.com/USGS-R/EGRET/wiki>

CC0
EGRETci
Exploration and Graphics for RivEr Trends (EGRET) Confidence Intervals
Collection of functions to evaluate uncertainty of results from water quality analysis using the Weighted Regressions on Time Discharge and Season (WRTDS) method. This package is an add-on to the EGRET package that performs the WRTDS analysis.
<https://github.com/USGS-R/EGRETci>

CC0
elevatr
Access Elevation Data from Various APIs
Several web services are available that provide access to elevation data. This package provides access to several of those services and returns elevation data either as a SpatialPointsDataFrame from point elevation services or as a raster object from raster elevation services. Currently, the package supports access to the Mapzen Elevation Service <https://mapzen.com/documentation/elevation/elevation-service/>, Mapzen Terrain Service <https://mapzen.com/documentation/terrain-tiles/>, Amazon Web Services Terrain Tiles <https://aws.amazon.com/public-datasets/terrain/> and the USGS Elevation Point Query Service <http://ned.usgs.gov/epqs/>.
<https://www.github.com/usepa/elevatr>

CC0
etl
Extract-Transform-Load Framework for Medium Data
A framework for loading medium-sized data from the Internet to a local or remote relational database management system. This package itself doesn't do much more than provide a toy example and set up the method structure. Packages that depend on this package will facilitate the construction and maintenance of their respective databases.
<http://github.com/beanumber/etl>

CC0
europop
Historical Populations of European Cities, 1500-1800
This dataset contains population estimates of all European cities with at least 10,000 inhabitants during the period 1500-1800. These data are adapted from Jan De Vries, "European Urbanization, 1500-1800" (1984).
<https://github.com/mdlincoln/europop>

CC0
fakeR
Simulates Data from a Data Frame of Different Variable Types
Generates fake data from a dataset of different variable types. The package contains the functions simulate\_dataset and simulate\_dataset\_ts to simulate time-independent and time-dependent data. It randomly samples character and factor variables from contingency tables and numeric and ordered factors from a multivariate normal distribution. It currently supports the simulation of stationary and zero-inflated count time series.

CC0
fancycut
A Fancy Version of 'base::cut'
Provides the function fancycut() which is like cut() except you can mix left open and right open intervals with point values, intervals that are closed on both ends and intervals that are open on both ends.

CC0
favnums
A Dataset of Favourite Numbers
A dataset of favourite numbers, selected from an online poll of over 30,000 people by Alex Bellos (<http://pages.bloomsbury.com/favouritenumber>).

CC0
fermicatsR
Fermi Large Area Telescope Catalogs
Data from various catalogs of astrophysical gamma-ray sources detected by NASA's Large Area Telescope (The Astrophysical Journal, 697, 1071, 2009 June 1), on board the Fermi gamma-ray satellite. More information on Fermi and its data products is available from the Fermi Science Support Center (<http://fermi.gsfc.nasa.gov/ssc/>).
<https://github.com/sazpark/fermicatsR.git>

CC0
FFTrees
Generate, Visualise, and Compare Fast and Frugal Decision Trees
Create, visualise, and test fast and frugal decision trees (FFTrees). FFTrees are very simple decision trees for classifying cases (i.e.; breast cancer patients) into one of two classes (e.g.; no cancer vs. true cancer) based on a small number of cues (e.g.; test results). FFTrees can be preferable to more complex algorithms because they are easy to communicate, require very little information, and are robust against overfitting.

CC0
fueleconomy
EPA fuel economy data
Fuel economy data from the EPA, 1985-2015, conveniently packaged for consumption by R users.
<http://github.com/hadley/fueleconomy>

CC0
gaiah
Genetic and Isotopic Assignment Accounting for Habitat Suitability
Tools for using genetic markers, stable isotope data, and habitat suitability data to calculate posterior probabilities of breeding origin of migrating birds.

CC0
gapminder
Data from Gapminder
An excerpt of the data available at Gapminder.org. For each of 142 countries, the package provides values for life expectancy, GDP per capita, and population, every five years, from 1952 to 2007.
<https://github.com/jennybc/gapminder>, <http://www.gapminder.org/data/>

CC0
geoknife
Web-Processing of Large Gridded Datasets
Processes gridded datasets found on the U.S. Geological Survey Geo Data Portal web application or elsewhere, using a web-enabled workflow that eliminates the need to download and store large datasets that are reliably hosted on the Internet. The package provides access to several data subset and summarization algorithms that are available on remote web processing servers.
<https://github.com/USGS-R/geoknife>

CC0
gesis
R Client for GESIS Data Catalogue (DBK)
Provides programmatic access to the GESIS - Leibniz-Institute for the Social Sciences Data Catalogue/Datenbestandkatalog (DBK), which maintains a large repository of data sets related to the social sciences. See <http://www.gesis.org> for more information.

CC0
hdr
Interface to the UNDR Human Development Report API
Provides a complete interface to the United Nations Development Programme Human Development Report API (<http://hdr.undp.org>). The API includes a large amount of human development data, including all the series used to compute the Human Development Index (HDI), as well as the HDI itself.
<http://hdr.undp.org> ; <https://github.com/expersso/hdr>

CC0
hflights
Flights that departed Houston in 2011
A data only package containing commercial domestic flights that<U+000a>departed Houston (IAH and HOU) in 2011.

CC0
housingData
U.S. Housing Data from 2008 to 2016
Monthly median home listing, sale price per square foot, and number of units sold for 2984 counties in the contiguous United States From 2008 to January 2016. Additional data sets containing geographical information and links to Wikipedia are also included.
<http://github.com/hafen/housingData>

CC0
ie2misc
Irucka Embry's Miscellaneous USGS Functions
A collection of Irucka Embry's miscellaneous USGS functions (processing .exp and .psf files, statistical error functions, "+" dyadic operator for use with NA, creating ADAPS and QW spreadsheet files, calculating saturated enthalpy). Irucka created these functions while a Cherokee Nation Technology Solutions (CNTS) United States Geological Survey (USGS) Contractor and/or USGS employee.
<https://github.com/iembry-USGS/ie2misc>

CC0
ig.vancouver.2014.topcolour
Instagram 2014 Vancouver Top Colour Dataset
A dataset of the top colours of photos from Instagram taken in 2014 in the city of Vancouver, British Columbia, Canada. It consists of: top colour and counts data. This data was obtained using the Instagram API. Instagram is a web photo sharing service. It can be found at: <https://instagram.com>. The Instagram API is documented at: <https://instagram.com/developer/>.

CC0
imputeTestbench
Test Bench for the Comparison of Imputation Methods
Provides a test bench for the comparison of missing data imputation methods in univariate time series. Imputation methods are compared using RMSE, MAE or MAPE error metrics. Proposed imputation methods and alternative error metrics can be used.
<http://www.neerajbokde.com/cran/imputetestbench>

CC0
inegiR
Integrate INEGI’s (Mexican Stats Office) API with R
Provides functions to download and parse information from INEGI (Official Mexican statistics agency).

CC0
inlmisc
Miscellaneous Functions for the USGS INL Project Office
A collection of functions for creating high-level graphics, performing raster-based analysis, processing MODFLOW-based models, and overlaying multi-polygon objects. Used to support packages and scripts written by researchers at the United States Geological Survey (USGS) Idaho National Laboratory (INL) Project Office.
<https://github.com/USGS-R/inlmisc>

CC0
KScorrect
Lilliefors-Corrected Kolmogorov-Smirnoff Goodness-of-Fit Tests
Implements the Lilliefors-corrected Kolmogorov-Smirnoff test for use in goodness-of-fit tests, suitable when population parameters are unknown and must be estimated by sample statistics. P-values are estimated by simulation. Can be used with a variety of continuous distributions, including normal, lognormal, univariate mixtures of normals, uniform, loguniform, exponential, gamma, and Weibull distributions. Functions to generate random numbers and calculate density, distribution, and quantile functions are provided for use with the log uniform and mixture distributions.
<https://github.com/pnovack-gottshall/KScorrect>

CC0
lakemorpho
Lake Morphometry Metrics in R
Lake morphometry metrics are used by limnologists to understand, among other things, the ecological processes in a lake. Traditionally, these metrics are calculated by hand, with planimeters, and increasingly with commercial GIS products. All of these methods work; however, they are either outdated, difficult to reproduce, or require expensive licenses to use. The lakemorpho package provides the tools to calculate a typical suite of these metrics from an input elevation model and lake polygon. The metrics currently supported are: fetch, major axis, minor axis, maximum length, maximum width, mean width,maximum depth, mean depth, shoreline development, shoreline length, surface area, and volume.
<http://www.github.com/USEPA/lakemorpho>

CC0
laketemps
Lake Temperatures Collected by Situ and Satellite Methods from 1985-2009
Lake temperature records, metadata, and climate drivers for 291 global lakes during the time period 1985-2009. Temperature observations were collected using satellite and in situ methods. Climatic drivers and geomorphometric characteristics were also compiled and are included for each lake. Data are part of the associated publication from the Global Lake Temperature Collaboration project (<http://www.laketemperature.org>). See citation('laketemps') for dataset attribution.

CC0
macleish
Retrieve Data from MacLeish Field Station
Download data from the Ada and Archibald MacLeish Field Station in Whately, MA. The Ada and Archibald MacLeish Field Station is a 260-acre patchwork of forest and farmland located in West Whately, MA that provides opportunities for faculty and students to pursue environmental research, outdoor education, and low-impact recreation (see <http://www.smith.edu/ceeds/macleish.php> for more information). This package contains weather data over several years, and spatial data on various man-made and natural structures.

CC0
maddison
Maddison Project Database
Contains the Maddison Project database, which provides estimates of GDP per capita for all countries in the world between AD 1 and 2010. See <http://www.ggdc.net/maddison> for more information.
<http://www.ggdc.net/maddison> <https://github.com/expersso/maddison>

CC0
maxmatching
Maximum Matching for General Weighted Graph
Computes the maximum matching for unweighted graph and maximum matching for (un)weighted bipartite graph efficiently.

CC0
mded
Measuring the Difference Between Two Empirical Distributions
Provides a function for measuring the difference between two independent or non-independent empirical distributions and returning a significance level of the difference.

CC0
mdsr
Complement to 'Modern Data Science with R'
A complement to *Modern Data Science with R* (ISBN: 978-1498724487, publisher URL: <https://www.crcpress.com/Modern-Data-Science-with-R/Baumer-Kaplan-Horton/p/book/9781498724487>). This package contains all of the data and code necessary to complete exercises and reproduce examples from the text. It also facilitates connections to the SQL database server used in the book.
<http://github.com/beanumber/mdsr>

CC0
MIDN
Nearly Exact Sample Size Calculation for Exact Powerful Nonrandomized Tests for Differences Between Binomial Proportions
Implementation of the mid-n algorithms presented in Wellek S (2015) <DOI:10.1111/stan.12063> Statistica Neerlandica 69, 358-373 for exact sample size calculation for superiority trials with binary outcome.

CC0
modeval
Evaluation of Classification Model Options
Designed to assist novice to intermediate analysts in choosing an optimal classification model, particularly for working with relatively small data sets. It provides cross-validated results comparing several different models at once using a consistent set of performance metrics, so users can hone in on the most promising approach rather than attempting single model fittings at a time. The package predefined 12 most common classification models, although users are free to select from the 200+ other options available in caret package.

CC0
mycor
Automatic Correlation and Regression Test in a Data Frame
Perform correlation and linear regression test among the numeric columns in a data frame automatically and make plots using pairs or lattice::parallelplot.

CC0
nasadata
Interface to Various NASA API's
Provides functions to access NASA's Earth Imagery and Assets API and the Earth Observatory Natural Event Tracker (EONET) webservice.

CC0
NeuralNetTools
Visualization and Analysis Tools for Neural Networks
Visualization and analysis tools to aid in the interpretation of neural network models. Functions are available for plotting, quantifying variable importance, conducting a sensitivity analysis, and obtaining a simple list of model weights.

CC0
NHLData
Scores for Every Season Since the Founding of the NHL in 1917
Each dataset contains scores for every game during a specific season of the NHL.

CC0
nycflights13
Flights that Departed NYC in 2013
Airline on-time data for all flights departing NYC in 2013. Also includes useful 'metadata' on airlines, airports, weather, and planes.
<http://github.com/hadley/nycflights13>

CC0
nzpullover
Driving Offences in New Zealand Between 2009 and 2016
Datasets of driving offences and fines in New Zealand between 2009 and 2016. Originally published by the New Zealand Police at <http://www.police.govt.nz/about-us/publication/road-policing-driver-offence-data-january-2009-september-2016>.
<https://github.com/nacnudus/nzpullover>

CC0
OECD
Search and Extract Data from the OECD
Search and extract data from the OECD.
<https://www.github.com/expersso/OECD>

CC0
okcupiddata
OkCupid Profile Data for Introductory Statistics and Data Science Courses
Cleaned profile data from "OkCupid Profile Data for Introductory Statistics and Data Science Courses" (Journal of Statistics Education 2015 <http://www.amstat.org/publications/jse/v23n2/kim.pdf>).
<http://github.com/rudeboybert/okcupiddata>

CC0
paleotree
Paleontological and Phylogenetic Analyses of Evolution
Provides tools for transforming, a posteriori time-scaling, and modifying phylogenies containing extinct (i.e. fossil) lineages. In particular, most users are interested in the functions timePaleoPhy(), bin\_timePaleoPhy(), cal3TimePaleoPhy() and bin\_cal3TimePaleoPhy(), which a posteriori time-scale cladograms of fossil taxa into dated phylogenies. This package also contains a large number of likelihood functions for estimating sampling and diversification rates from different types of data available from the fossil record (e.g. range data, occurrence data, etc). paleotree users can also simulate diversification and sampling in the fossil record using the function simFossilRecord(), which is a detailed simulator for branching birth-death-sampling processes composed of discrete taxonomic units arranged in ancestor-descendant relationships. Users can use simFossilRecord() to simulate diversification in incompletely sampled fossil records, under various models of morphological differentiation (i.e. the various patterns by which morphotaxa originate from one another), and with time-dependent, longevity-dependent and/or diversity-dependent rates of diversification, extinction and sampling. Additional functions allow users to translate simulated ancestor-descendant data from simFossilRecord() into standard time-scaled phylogenies or unscaled cladograms that reflect the relationships among taxon units.
<https://github.com/dwbapst/paleotree> , <http://webpages.sdsmt.edu/~dbapst/>

CC0
PARSE
Model-Based Clustering with Regularization Methods for High-Dimensional Data
Model-based clustering and identifying informative features based on regularization methods. The package includes three regularization methods - PAirwise Reciprocal fuSE (PARSE) penalty proposed by Wang, Zhou and Hoeting (2016), the adaptive L1 penalty (APL1) and the adaptive pairwise fusion penalty (APFP). Heatmaps are included to shown the identification of informative features.

CC0
pdftables
Programmatic Conversion of PDF Tables
Allows the user to convert PDF tables to formats more amenable to analysis ('.csv', '.xml', or '.xlsx') by wrapping the PDFTables API. In order to use the package, the user needs to sign up for an API account on the PDFTables website (<https://pdftables.com/pdf-to-excel-api>). The package works by taking a PDF file as input, uploading it to PDFTables, and returning a file with the extracted data.
<https://www.github.com/expersso/pdftables> , <https://pdftables.com>

CC0
pmc
Phylogenetic Monte Carlo
Monte Carlo based model choice for applied phylogenetics of continuous traits. Method described in Carl Boettiger, Graham Coop, Peter Ralph (2012) Is your phylogeny informative? Measuring the power of comparative methods, Evolution 66 (7) 2240-51. <doi:10.1111/j.1558-5646.2011.01574.x>.
<https://github.com/cboettig/pmc>

CC0
poliscidata
Datasets and Functions Featured in Pollock and Edwards R Companion to Essential of Political Analysis
Bundles the datasets and functions used in the book by Philip Pollock and Barry Edwards, An R Companion to Essentials of Political Analysis, Second Edition.

CC0
quickmapr
Quickly Map and Explore Spatial Data
While analyzing geospatial data, easy visualization is often needed that allows for quick plotting, and simple, but easy interactivity. Additionally, visualizing geospatial data in projected coordinates is also desirable. The 'quickmapr' package provides a simple method to visualize 'sp' and 'raster' objects, allows for basic zooming, panning, identifying, labeling, selecting, and measuring spatial objects. Importantly, it does not require that the data be in geographic coordinates.
<https://www.github.com/jhollist/quickmapr>

CC0
randomcoloR
Generate Attractive Random Colors
Simple methods to generate attractive random colors. The random colors are from a wrapper of 'randomColor.js' <https://github.com/davidmerfield/randomColor>. In addition, it also generates optimally distinct colors based on k-means (inspired by 'IWantHue' <https://github.com/medialab/iwanthue>).

CC0
raw
R Actuarial Workshops
In order to facilitate R instruction for actuaries, we have organized several sets of publicly available data of interest to non-life actuaries. In addition, we suggest a set of packages, which most practicing actuaries will use routinely. Finally, there is an R markdown skeleton for basic reserve analysis.
<http://pirategrunt.com/raw/>

CC0
rcorpora
A Collection of Small Text Corpora of Interesting Data
A collection of small text corpora of interesting data. It contains all data sets from <https://github.com/dariusk/corpora>. Some examples: names of animals: birds, dinosaurs, dogs; foods: beer categories, pizza toppings; geography: English towns, rivers, oceans; humans: authors, US presidents, occupations; science: elements, planets; words: adjectives, verbs, proverbs, US president quotes.
<https://github.com/gaborcsardi/rcorpora>

CC0
resampledata
Data Sets for Mathematical Statistics with Resampling in R
Package of data sets from "Mathematical Statistics with Resampling in R" (2011) by Laura Chihara and Tim Hesterberg.
<https://github.com/rudeboybert/resampledata>

CC0
rfigshare
An R Interface to 'figshare'
An interface to 'figshare' (<http://figshare.com>), a scientific repository to archive and assign 'DOIs' to data, software, figures, and more.
<https://github.com/ropensci/rfigshare>

CC0
rfishbase
R Interface to 'FishBase'
A programmatic interface to <http://www.fishbase.org>, re-written based on an accompanying 'RESTful' API. Access tables describing over 30,000 species of fish, their biology, ecology, morphology, and more. This package also supports experimental access to <http://www.sealifebase.org> data, which contains nearly 200,000 species records for all types of aquatic life not covered by 'FishBase.'
<https://github.com/ropensci/rfishbase>

CC0
rnaturalearthdata
World Vector Map Data from Natural Earth Used in 'rnaturalearth'
Vector map data from <http://www.naturalearthdata.com/>. Access functions are provided in the accompanying package 'rnaturalearth'.
<https://github.com/ropenscilabs/rnaturalearthdata>

CC0
rpdo
Pacific Decadal Oscillation Index Data
Monthly Pacific Decadal Oscillation (PDO) index values from January 1900 to February 2017. Includes download\_pdo() to scrape the latest values from <http://research.jisao.washington.edu/pdo/PDO.latest>.
<https://github.com/poissonconsulting/rpdo>

CC0
RSCABS
Rao-Scott Cochran-Armitage by Slices Trend Test
Performs the Rao-Scott Cochran-Armitage by Slices trend test (RSCABS) used in analysis of histopathological endpoints. It has functions for both command line operations along with a built in GUI.
<https://CRAN.R-project.org/package=RSCABS>

CC0
RSurvey
Geographic Information System Application
A geographic information system (GIS) graphical user interface (GUI) that provides data viewing, management, and analysis tools.
<https://github.com/USGS-R/RSurvey>

CC0
sbtools
USGS ScienceBase Tools
Tools for interacting with U.S. Geological Survey ScienceBase <https://www.sciencebase.gov> interfaces. ScienceBase is a data cataloging and collaborative data management platform. Functions included for querying ScienceBase, and creating and fetching datasets.
<https://github.com/USGS-R/sbtools>

CC0
spTest
Nonparametric Hypothesis Tests of Isotropy and Symmetry
Implements nonparametric hypothesis tests to check isotropy and symmetry properties for two-dimensional spatial data.

CC0
stocc
Fit a Spatial Occupancy Model via Gibbs Sampling
Fit a spatial-temporal occupancy models using a probit formulation instead of a traditional logit model.

CC0
superheat
A Graphical Tool for Exploring Complex Datasets Using Heatmaps
A system for generating extendable and customizable heatmaps for exploring complex datasets, including big data and data with multiple data types.

CC0
SWMPr
Retrieving, Organizing, and Analyzing Estuary Monitoring Data
Tools for retrieving, organizing, and analyzing environmental data from the System Wide Monitoring Program of the National Estuarine Research Reserve System <http://cdmo.baruch.sc.edu/>. These tools address common challenges associated with continuous time series data for environmental decision making.

CC0
titanic
Titanic Passenger Survival Data Set
This data set provides information on the fate of passengers on the fatal maiden voyage of the ocean liner "Titanic", summarized according to economic status (class), sex, age and survival. Whereas the base R Titanic data found by calling data("Titanic") is an array resulting from cross-tabulating 2201 observations, these data sets are the individual non-aggregated observations and formatted in a machine learning context with a training sample, a testing sample, and two additional data sets that can be used for deeper machine learning analysis. These data sets are also the data sets downloaded from the Kaggle competition and thus lowers the barrier to entry for users new to R or machine learing.
<https://github.com/paulhendricks/titanic>

CC0
treebase
Discovery, Access and Manipulation of 'TreeBASE' Phylogenies
Interface to the API for 'TreeBASE' <http://treebase.org> from 'R.' 'TreeBASE' is a repository of user-submitted phylogenetic trees (of species, population, or genes) and the data used to create them.
<https://github.com/ropensci/treebase>

CC0
ttbbeer
US Beer Statistics from TTB
U.S. Department of the Treasury, Alcohol and Tobacco Tax and Trade Bureau (TTB) collects data and reports on monthly beer industry production and operations. This data package includes a collection of 10 years (2006 - 2015) worth of data on materials used at U.S. breweries in pounds reported by the Brewer's Report of Operations and the Quarterly Brewer's Report of Operations forms, ready for data analysis. This package also includes historical tax rates on distilled spirits, wine, beer, champagne, and tobacco products as individual data sets.
<https://github.com/jasdumas/ttbbeer>

CC0
USGSstates2k
Replaced by 'states2k' -- United States of America Map with the NAD 1983 Albers Projection
A map of the USA from the United States Geological Survey (USGS). Irucka worked with this data set while a Cherokee Nation Technology Solutions (CNTS) USGS Contractor and/or USGS employee. It is replaced by 'states2k'.

CC0
vaersNDvax
Non-Domestic Vaccine Adverse Event Reporting System (VAERS) Vaccine Data for Present
Non-Domestic VAERS vaccine data for 01/01/2016 - 06/14/2016. If you want to explore the full VAERS data for 1990 - Present (data, symptoms, and vaccines), then check out the 'vaersND' package from the URL below. The URL and BugReports below correspond to the 'vaersND' package, of which 'vaersNDvax' is a small subset (2016 only). 'vaersND' is not hosted on CRAN due to the large size of the data set. To install the Suggested 'vaers' and 'vaersND' packages, use the following R code: 'devtools::install\_git("<https://gitlab.com/iembry/vaers.git>", build\_vignettes = TRUE)' and 'devtools::install\_git("<https://gitlab.com/iembry/vaersND.git>", build\_vignettes = TRUE)'. "VAERS is a national vaccine safety surveillance program co-sponsored by the US Centers for Disease Control and Prevention (CDC) and the US Food and Drug Administration (FDA). VAERS is a post-marketing safety surveillance program, collecting information about adverse events (possible side effects) that occur after the administration of vaccines licensed for use in the United States." For more information about the data, visit <https://vaers.hhs.gov/index>. For information about vaccination/immunization hazards, visit <http://www.questionuniverse.com/rethink.html/#vaccine>.
<https://gitlab.com/iembry/vaersND>

CC0
vaersvax
US Vaccine Adverse Event Reporting System (VAERS) Vaccine Data for Present
US VAERS vaccine data for 01/01/2016 - 06/14/2016. If you want to explore the full VAERS data for 1990 - Present (data, symptoms, and vaccines), then check out the 'vaers' package from the URL below. The URL and BugReports below correspond to the 'vaers' package, of which 'vaersvax' is a small subset (2016 only). 'vaers' is not hosted on CRAN due to the large size of the data set. To install the Suggested 'vaers' and 'vaersND' packages, use the following R code: 'devtools::install\_git("<https://gitlab.com/iembry/vaers.git>", build\_vignettes = TRUE)' and 'devtools::install\_git("<https://gitlab.com/iembry/vaersND.git>", build\_vignettes = TRUE)'. "VAERS is a national vaccine safety surveillance program co-sponsored by the US Centers for Disease Control and Prevention (CDC) and the US Food and Drug Administration (FDA). VAERS is a post-marketing safety surveillance program, collecting information about adverse events (possible side effects) that occur after the administration of vaccines licensed for use in the United States." For more information about the data, visit <https://vaers.hhs.gov/index>. For information about vaccination/immunization hazards, visit <http://www.questionuniverse.com/rethink.html/#vaccine>.
<https://gitlab.com/iembry/vaers>

CC0
valottery
Results from the Virginia Lottery Draw Games
Historical results for the state of Virginia lottery draw games. Data were downloaded from <https://www.valottery.com/>.

CC0
WHO
R Client for the World Health Organization API
Provides programmatic access to the World Health Organization API.
<https://www.github.com/expersso/WHO> , <http://www.who.int>

CC0
WRTDStidal
Weighted Regression for Water Quality Evaluation in Tidal Waters
An adaptation for estuaries (tidal waters) of weighted regression on time, discharge, and season to evaluate trends in water quality time series.

CC0
xtractomatic
Accessing Environmental Data from ERD's ERDDAP Server
Contains three functions that access environmental data from ERD's ERDDAP service <http://coastwatch.pfeg.noaa.gov/erddap>. The xtracto() function extracts data along a trajectory for a given "radius" around the point. The xtracto\_3D() function extracts data in a box. The xtractogon() function extracts data in a polygon. There are also two helper functions to obtain information about available data.
<https://github.com/rmendels/xtractomatic>

Unlimited
antitrust
Tools for Antitrust Practitioners
A collection of tools for antitrust practitioners, including the ability to calibrate different consumer demand systems and simulate the effects mergers under different competitive regimes.

Unlimited
aqfig
Functions to help display air quality model output and<U+000a>monitoring data
This package contains functions to help display air quality model output and monitoring data, such as creating color scatterplots, color legends, etc.

Unlimited
boot
Bootstrap Functions (Originally by Angelo Canty for S)
Functions and datasets for bootstrapping from the book "Bootstrap Methods and Their Application" by A. C. Davison and D. V. Hinkley (1997, CUP), originally written by Angelo Canty for S.

Unlimited
BurStFin
Burns Statistics Financial
A suite of functions for finance, including the estimation<U+000a>of variance matrices via a statistical factor model or<U+000a>Ledoit-Wolf shrinkage.
<http://www.burns-stat.com/>

Unlimited
BurStMisc
Burns Statistics Miscellaneous
Script search, corner, genetic optimization, permutation tests, write expect test.

Unlimited
CIAAWconsensus
Isotope Ratio Meta-Analysis
Calculation of consensus values for atomic weights, isotope amount ratios, and isotopic abundances with the associated uncertainties using multivariate meta-regression approach for consensus building.

Unlimited
DAAGxtras
Data Sets and Functions, supplementary to DAAG
various data sets used in additional exercises for<U+000a>the book Maindonald, J.H. and Braun, W.J. (3rd edn 2010)<U+000a>"Data Analysis and Graphics Using R", and for a<U+000a>'Data Mining' course. Note that a number of datasets<U+000a>that were in earlier versions of this package have been<U+000a>transferred to the DAAG package.
<http://www.maths.anu.edu.au/~johnm>

Unlimited
dhglm
Double Hierarchical Generalized Linear Models
Double hierarchical generalized linear models in which the mean, dispersion parameters for variance of random effects, and residual variance (overdispersion) can be further modeled as random-effect models.

Unlimited
hierarchicalDS
Functions For Performing Hierarchical Analysis of Distance Sampling Data
Functions for performing hierarchical analysis of distance sampling data, with ability to use an areal spatial ICAR model on top of user supplied covariates to get at variation in abundance intensity. The detection model can be specified as a function of observer and individual covariates, where a parametric model is supposed for the population level distribution of covariate values. The model uses data augmentation and a reversible jump MCMC algorithm to sample animals that were never observed. Also included is the ability to include point independence (increasing correlation multiple observer's observations as a function of distance, with independence assumed for distance=0 or first distance bin), as well as the ability to model species misclassification rates using a multinomial logit formulation on data from double observers. New in version 2.1 is the ability to include zero inflation, but this is only recommended for cases where sample sizes and spatial coverage of the survey are high.

Unlimited
HIV.LifeTables
HIV calibrated model life tables for countries with generalized HIV epidemics
The functions in this package produce a complete set of mortality rates as a function of a combination of HIV prevalence and either life expectancy at birth (e0), child mortality (5q0), or child mortality with adult mortality (45q15)

Unlimited
learningr
Data and functions to accompany the book "Learning R".
Crabs in the English channel, deer skulls, English monarchs, half-caste Manga characters, Jamaican cities, Shakespeare's The Tempest, drugged up cyclists and sexually transmitted diseases.

Unlimited
LSD
Lots of Superior Depictions
Create lots of colorful plots in a plethora of variations (try the LSD demotour() )

Unlimited
mdhglm
Multivariate Double Hierarchical Generalized Linear Models
Allows various models for multivariate response variables where each response is assumed to follow double hierarchical generalized linear models. In double hierarchical generalized linear models, the mean, dispersion parameters for variance of random effects, and residual variance can be further modeled as random-effect models.

Unlimited
ModelMap
Modeling and Map Production using Random Forest and Stochastic Gradient Boosting
Creates sophisticated models of training data and validates the models with an independent test set, cross validation, or in the case of Random Forest Models, with Out Of Bag (OOB) predictions on the training data. Create graphs and tables of the model validation results. Applies these models to GIS .img files of predictors to create detailed prediction surfaces. Handles large predictor files for map making, by reading in the .img files in chunks, and output to the .txt file the prediction for each data chunk, before reading the next chunk of data.

Unlimited
MPV
Data Sets from Montgomery, Peck and Vining's Book
Most of this package consists of data sets from the textbook Introduction to Linear Regression Analysis (3rd ed), by Montgomery et al. Some additional data sets and functions useful in an undergraduate regression course are included.

Unlimited
muckrock
Data on Freedom of Information Act Requests
A data package containing public domain information on requests made by the 'MuckRock' (<https://www.muckrock.com/>) project under the United States Freedom of Information Act.
<https://github.com/Ironholds/muckrock/>

Unlimited
nonmem2R
Loading NONMEM Output Files and Simulate with Parameter Uncertainty
Loading NONMEM (NONlinear Mixed-Effect Modeling, <http://www.iconplc.com/innovation/solutions/nonmem/> ) output files and simulate with parameter uncertainty.

Unlimited
pathological
Path Manipulation Utilities
Utilities for paths, files and directories.
<https://github.com/richierocks/pathological>

Unlimited
qfasar
Quantitative Fatty Acid Signature Analysis in R
An implementation of Quantitative Fatty Acid Signature Analysis (QFASA) in R. QFASA is a method of estimating the diet composition of predators. The fundamental unit of information in QFASA is a fatty acid signature (signature), which is a vector of proportions describing the composition of fatty acids within lipids. Signature data from at least one predator and from samples of all potential prey types are required. Calibration coefficients, which adjust for the differential metabolism of individual fatty acids by predators, are also required. Given those data inputs, a predator signature is modeled as a mixture of prey signatures and its diet estimate is obtained as the mixture that minimizes a measure of distance between the observed and modeled signatures. A variety of estimation options and simulation capabilities are implemented. Please refer to the vignette for additional details and references.

Unlimited
readOffice
Read Text Out of Modern Office Files
Reads in text from 'unstructured' modern Microsoft Office files (XML based files) such as Word and PowerPoint. This does not read in structured data (from Excel or Access) as there are many other great packages to that do so already.

Unlimited
rebus
Build Regular Expressions in a Human Readable Way
Build regular expressions piece by piece using human readable code. This package is designed for interactive use. For package development, use the rebus.\* dependencies.

Unlimited
rebus.base
Core Functionality for the 'rebus' Package
Build regular expressions piece by piece using human readable code. This package contains core functionality, and is primarily intended to be used by package developers.
<https://github.com/richierocks/rebus.base>

Unlimited
rebus.datetimes
Date and Time Extensions for the 'rebus' Package
Build regular expressions piece by piece using human readable code. This package contains date and time functionality, and is primarily intended to be used by package developers.

Unlimited
rebus.numbers
Numeric Extensions for the 'rebus' Package
Build regular expressions piece by piece using human readable code. This package contains number-related functionality, and is primarily intended to be used by package developers.

Unlimited
rebus.unicode
Unicode Extensions for the 'rebus' Package
Build regular expressions piece by piece using human readable code. This package contains Unicode functionality, and is primarily intended to be used by package developers.

Unlimited
runittotestthat
Convert 'RUnit' Test Functions into 'testthat' Tests
Automatically convert a file or package worth of 'RUnit' test functions into 'testthat' tests.

Unlimited
seawaveQ
U.S. Geological Survey seawaveQ model
A model and utilities for analyzing trends in chemical concentrations in streams with a seasonal wave (seawave) and adjustment for streamflow (Q) and other ancillary variables
<http://dx.doi.org/10.3133/ofr20131255>

Unlimited
setter
Mutators that Work with Pipes
Mutators to set attributes of variables, that work well in a pipe (much like stats::setNames()).
<https://bitbucket.org/richierocks/setter>

Unlimited
sig
Print Function Signatures
Print function signatures and find overly complicated code.

Unlimited
Simile
Interact with Simile Models
Allows a Simile model saved as a compiled binary to be loaded, parameterized, executed and interrogated. This version works with Simile v5.97 on.

Unlimited
TaoTeProgramming
Illustrations from Tao Te Programming
Art-like behavior based on randomness

Unlimited
TotalCopheneticIndex
Total Cophenetic Index
Quantifies how balanced a phylogenetic tree is, using the Total Cophenetic Index - per A. Mir, F. Rossello, L. A. Rotger (2013), A new balance index for phylogenetic trees. Math. Biosci. 241, 125-136 <DOI:10.1016/j.mbs.2012.10.005>.
<https://github.com/ms609/tci>

CC BY-NC-SA 4.0
BaTFLED3D
Bayesian Tensor Factorization Linked to External Data
BaTFLED is a machine learning algorithm designed to make predictions and determine interactions in data that varies along three independent modes. For example BaTFLED was developed to predict the growth of cell lines when treated with drugs at different doses. The first mode corresponds to cell lines and incorporates predictors such as cell line genomics and growth conditions. The second mode corresponds to drugs and incorporates predictors indicating known targets and structural features. The third mode corresponds to dose and there are no dose-specific predictors (although the algorithm is capable of including predictors for the third mode if present). See 'BaTFLED3D\_vignette.rmd' for a simulated example.

CC BY-NC-SA 4.0
D3GB
Interactive Genome Browser with R
Creates interactive genome browser with 'R'. It joins the data analysis power of R and the visualization libraries of JavaScript in one package.
<http://d3gb.usal.es>

CC BY-NC-SA 4.0
IDmining
Intrinsic Dimension for Data Mining
Contains techniques for mining large high-dimensional data sets by using the concept of Intrinsic Dimension (ID). Here the ID is not necessarily integer. It is extended to fractal dimensions. And the Morisita estimator is used for the ID estimation, but other tools are included as well.
&lt; <https://www.sites.google.com/site/jeangolayresearch/> &gt;

CC BY-NC-SA 4.0
modes
Find the Modes and Assess the Modality of Complex and Mixture Distributions, Especially with Big Datasets
Designed with a dual purpose of accurately estimating the mode (or modes) as well as characterizing the modality of data. The specific application area includes complex or mixture distributions particularly in a big data environment. The heterogeneous nature of (big) data may require deep introspective statistical and machine learning techniques, but these statistical tools often fail when applied without first understanding the data. In small datasets, this often isn't a big issue, but when dealing with large scale data analysis or big data thoroughly inspecting each dimension typically yields an O(n^n-1) problem. As such, dealing with big data require an alternative toolkit. This package not only identifies the mode or modes for various data types, it also provides a programmatic way of understanding the modality (i.e. unimodal, bimodal, etc.) of a dataset (whether it's big data or not). See <http://www.sdeevi.com/modes_package> for examples and discussion.
<http://www.sdeevi.com/modes_package> <https://github.com/sathish-deevi/modes-Package/>

CC BY-NC-SA 4.0
netCoin
Interactive Networks with R
Create interactive networked coincidences. It joins the data analysis power of R to study coincidences and the visualization libraries of JavaScript in one package.

CC BY-NC-SA 4.0
nettools
A Network Comparison Framework
A collection of network inference methods for co-expression networks, quantitative network distances and a novel framework for network stability analysis.

CC BY-NC-SA 4.0
PepSAVIms
PepSAVI-MS Data Analysis
An implementation of the data processing and data analysis portion of a pipeline named the PepSAVI-MS which is currently under development by the Hicks laboratory at the University of North Carolina. The statistical analysis package presented herein provides a collection of software tools used to facilitate the prioritization of putative bioactive peptides from a complex biological matrix. Tools are provided to deconvolute mass spectrometry features into a single representation for each peptide charge state, filter compounds to include only those possibly contributing to the observed bioactivity, and prioritize these remaining compounds for those most likely contributing to each bioactivity data set.
<https://github.com/dpritchLibre/PepSAVIms>

CC BY-NC-SA 4.0
Radviz
Project Multidimensional Data in 2D Space
An implementation of the radviz projection in R. It enables the visualization of multidimensional data while maintaining the relation to the original dimensions. This package provides functions to create and plot radviz projections, and a number of summary plots that enable comparison and analysis. For reference see Ankerst et al. (1996) <http://citeseer.ist.psu.edu/viewdoc/summary?doi=10.1.1.68.1811> for original implementation, see Di Caro et al. (2010) <DOI:10.1007/978-3-642-13672-6_13> for the original method for dimensional anchor arrangements.
<http://github.com/yannabraham/Radviz>

CC BY-NC-SA 4.0
RchivalTag
Analyzing Archival Tagging Data
A set of functions to generate, access and analyze standard data products from archival tagging data.

CC BY-NC-SA 4.0
tigger
R Tools for Inferring New Immunoglobulin Alleles from Rep-Seq Data
Infers the V genotype of an individual from immunoglobulin (Ig) repertoire-sequencing (Rep-Seq) data, including detection of any novel alleles. This information is then used to correct existing V allele calls from among the sample sequences.
<http://tigger.readthedocs.io>

CC BY-NC-SA 4.0
USAboundaries
Historical and Contemporary Boundaries of the United States of America
The boundaries for geographical units in the United States of America contained in this package include state, county, congressional district, and zip code tabulation area. Contemporary boundaries are provided by the U.S. Census Bureau (public domain). Historical boundaries for the years from 1629 to 2000 are provided form the Newberry Library's 'Atlas of Historical County Boundaries' (licensed CC BY-NC-SA). Additional high resolution data is provided in the 'USAboundariesData' package; this package provides an interface to access that data.
<https://github.com/ropensci/USAboundaries>

CC BY-SA 4.0
apercu
Apercu is Giving you a Quick Look at your Data
The goal is to print an "aperçu", a short view of a vector, a matrix, a data.frame, a list or an array. By default, it prints the first 5 elements of each dimension. By default, the number of columns is equal to the number of lines. If you want to control the selection of the elements, you can pass a list, with each element being a vector giving the selection for each dimension.

CC BY-SA 4.0
arnie
"Arnie" box office records 1982-2014
Arnold Schwarzenegger movie weekend box office records from<U+000a>1982-2014
<https://github.com/imanuelcostigan/arnie>

CC BY-SA 4.0
diverse
Diversity Measures for Complex Systems
Computes the most common diversity measures used in social and other sciences, and includes new measures from interdisciplinary research.
<https://github.com/mguevara/diverse>

CC BY-SA 4.0
inspectr
Perform Basic Checks of Dataframes
Check one column or multiple columns of a dataframe using the preset basic checks or your own functions. Enables checks without knowledge of lapply() or sapply().

CC BY-SA 4.0
install.load
Check, Install and Load CRAN & USGS GRAN Packages
The function `install_load` checks the local R library(ies) to see if the required package(s) is/are installed or not. If the package(s) is/are not installed, then the package(s) will be installed along with the required dependency(ies). This function pulls source or binary packages from the Rstudio-sponsored CRAN mirror and/or the USGS GRAN Repository. Lastly, the chosen package(s) is/are loaded. The function `load_package` simply loads the provided packages.
<https://gitlab.com/iembry/install.load>

CC BY-SA 4.0
shazam
Immunoglobulin Somatic Hypermutation Analysis
Provides a computational framework for Bayesian estimation of antigen-driven selection in immunoglobulin (Ig) sequences, providing an intuitive means of analyzing selection by quantifying the degree of selective pressure. Also provides tools to profile mutations in Ig sequences, build models of somatic hypermutation (SHM) in Ig sequences, and make model-dependent distance comparisons of Ig repertoires.
<http://shazam.readthedocs.io>

CC BY-SA 4.0
stackoverflow
Stack Overflow's Greatest Hits
Consists of helper functions collected from StackOverflow.com, a question and answer site for professional and enthusiast programmers.
<http://stackoverflow.com>

CC BY 4.0
pkmon
Least-Squares Estimator under k-Monotony Constraint for Discrete Functions
We implement two least-squares estimators under k-monotony constraint using a method based on the Support Reduction Algorithm from Groeneboom et al (2008) <DOI:10.1111/j.1467-9469.2007.00588.x>. The first one is a projection estimator on the set of k-monotone discrete functions. The second one is a projection on the set of k-monotone discrete probabilities. This package provides functions to generate samples from the spline basis from Lefevre and Loisel (2013) <DOI:10.1239/jap/1378401239>, and from mixtures of splines.

CC BY 4.0
reproducer
Reproduce Statistical Analyses and Meta-Analyses
Includes data analysis functions (e.g., to calculate effect sizes and 95% Confidence Intervals (CI) on Standardised Effect Sizes (d) for ABBA cross-over repeated-measures experimental designs), data presentation functions (e.g., density curve overlaid on histogram), and the data sets analyzed in different research papers in software engineering (e.g., related to software defect prediction or multi-site experiment concerning the extent to which structured abstracts were clearer and more complete than conventional abstracts) to streamline reproducible research in software engineering.
<http://madeyski.e-informatyka.pl/reproducible-research/>

CC BY-NC 4.0
sdcTarget
Statistical Disclosure Control Substitution Matrix Calculator
Classes and methods to calculate and evaluate target matrices for statistical disclosure control.
<http://statistics.lazaridis.eu>

CC BY-NC-SA 3.0
DATforDCEMRI
Deconvolution Analysis Tool for Dynamic Contrast Enhanced MRI
This package performs voxel-wise deconvolution analysis of<U+000a>DCE-MRI contrast agent concentration versus time data and<U+000a>generates the Impulse Response Function, which can be used to<U+000a>approximate commonly utilized kinetic parameters such as Ktrans<U+000a>and ve. An interactive advanced voxel diagnosis tool (AVDT) is<U+000a>also provided to facilitate easy navigation of voxel-wise data.

CC BY-SA 2.0 + file LICENSE
zipcode
U.S. ZIP Code database for geocoding
This package contains a database of city, state, latitude,<U+000a>and longitude information for U.S. ZIP codes from the<U+000a>CivicSpace Database (August 2004) augmented by Daniel Coven's<U+000a>federalgovernmentzipcodes.us web site (updated January 22,<U+000a>2012). Previous versions of this package (before 1.0) were<U+000a>based solely on the CivicSpace data, so an original version of<U+000a>the CivicSpace database is also included.

CC BY-SA 4.0 + file LICENSE
igraphdata
A Collection of Network Data Sets for the 'igraph' Package
A small collection of various network data sets, to use with the 'igraph' package: the Enron email network, various food webs, interactions in the immunoglobulin protein, the karate club network, Koenigsberg's bridges, visuotactile brain areas of the macaque monkey, UK faculty friendship network, domestic US flights network, etc.
<http://igraph.org>

MIT + file LICENSE | Unlimited
labeling
Axis Labeling
Provides a range of axis labeling algorithms

Unlimited | file LICENSE
waterData
Retrieval, Analysis, and Anomaly Calculation of Daily Hydrologic Time Series Data
Imports U.S. Geological Survey (USGS) daily hydrologic data from USGS web services (see <https://waterservices.usgs.gov/> for more information), plots the data, addresses some common data problems, and calculates and plots anomalies.
<https://pubs.usgs.gov/of/2012/1168/>

``` r
#x <- df[2, ]
#cat(paste(x$license, x$name, x$title, x$desc, x$url, sep = "  \n"))
```
