#' Solar array characteristics in the Chesapeake Bay Watershed.
#'
#' A GeoJson dataset containing polygons representing ground mounted solar arrays and associated geospatial 
#'  attributes.
#'
#' @format A data frame with 5548 rows and 24 variables:
#' \describe{
#'   \item{X}{price, in US dollars}
#'   \item{cultivated}{percentage of cultivated pixels within the polygon using the annual cropland data layers}
#'   \item{impervious}{mean percentage impervious surface per pixel within polygon boundaries from NLCD}
#'   \item{percent_tree_cover}{mean percentage of treecover per pixel within polygon boundaires from NLCD}
#'   \item{slope}{mean slope of pixels within polygon boundaries from USGS 3DEP 10m DEM}
#'   \item{year}{integer year in which polygon first appeared. 2030 for simulated polygons}
#'   \item{geometry}{GeoJSON format POLYGON coordinates in epsg:3857 proejction}
#'   \item{index}{unique identifier per feature}
#'   \item{mtfcc}{TIGER road type classifier for S1400 or S1200 road closest to polygon }
#'   \item{linearid}{Unique identifier for transmission lines}
#'   \item{road_dist}{distance (m) to nearest TIGER S1400 or S1200 road}
#'   \item{VOLTAGE}{voltage rating of nearest transmission line}
#'   \item{ID}{unique identifier for nearest transmission line}
#'   \item{line_dist}{distance (m) to nearest transmission line}
#'   \item{GAP_Sts}{GAP status code [1-4] from USGS PADUS data indicating level of biodiversity protection}
#'   \item{FIPS}{U.S. Census Bureau census tract identifier for 2020 census tract containing polygon centroids}
#'   \item{POPULATION}{2010 population of census tract containing polygon centroids}
#'   \item{income}{2010 median household income of census tract containing polygon centroids}
#'   \item{housing}{2010 housing density of census tract containing polygon centroids}
#'   \item{index_padus}{row identifier for USGS PADUS polygon containing polygon centroids}
#'   \item{Shape_Area}{Size (m2) of parcel containing polygon centroids}
#'   \item{county}{US Census Bureau FIPS code of county countaining polygon centroid}
#'   \item{index_parcel}{row identifier for parcel containing polygon centroids}
#' }
#' 
#' @source TIGER road data from ESRI Living Atlas \url{https://services.arcgis.com/P3ePLMYs2RVChkJx/arcgis/rest/services/TIGER_Roads_2021_view/FeatureServer}
#' @source NLCD Landcover data \url{https://www.sciencebase.gov/catalog/item/5f21cef582cef313ed940043}
#' @source USDA Cropland Data Layers \url{https://nassgeodata.gmu.edu/CropScape/}
#' @source Homeland Infrastructure Foundation-Level transmission line data from ESRI Living Atlas \url{https://services2.arcgis.com/FiaPA4ga0iQKduv3/arcgis/rest/services/US_Electric_Power_Transmission_Lines/FeatureServer}
#' @source USGS PADUS Fee Manager data from ESRI Living Atlas \url{https://services.arcgis.com/P3ePLMYs2RVChkJx/arcgis/rest/services/USA_Protected_Areas_Fee_Manager/FeatureServer}
#' 
#' 
#' @author Michael Evans \email{mevans@chesapeakeconservancy.org}
"solar_analysis_data"