packages<-c("cowplot", "dplyr", "geosphere", "ggplot2", "ggExtra", "maps", "maptools", "readxl", 
            "rgdal", "rgeos", "sf", "sp", "spatialreg", "spdep", "stringr","tidyr", "viridis")
sapply(packages, library, character.only=T)

#The Data
se.data <- read.csv("./Data/childpov18_southfull.csv", 
                   colClasses = c("character", "character", "character", 
                                  "numeric", "numeric", "numeric", 
                                  "numeric", "numeric", "numeric", 
                                  "numeric", "numeric", "numeric", 
                                  "numeric", "numeric", "numeric", 
                                  "numeric", "numeric", "numeric", 
                                  "numeric", "numeric", "numeric", 
                                  "numeric", "numeric", "numeric", 
                                  "numeric", "numeric", "numeric", 
                                  "numeric", "numeric", "numeric", 
                                  "numeric", "numeric", "numeric", 
                                  "numeric", "numeric"))

#Format Variables
names(se.data)[names(se.data)=="X2016.child.poverty"] <- "child.pov.2016"
se.data$FIPS <- str_pad(se.data$FIPS, 5, "left", pad = 0)

#create subsets
satl.data <- subset(se.data, State %in% c("DE", "DC", "FL", "GA", "MD", "NC", "SC", "VA", "WV"))
esc.data <- subset(se.data, State %in% c("AL", "KY", "MS", "TN"))
wsc.data <- subset(se.data, State %in% c("AR", "LA", "OK", "TX"))

#creating counties from shapefiles
se.shape<-readOGR(dsn="./Shapefiles",layer="SE_Counties_2016")
se.shape@data <- se.shape@data[-c(6:54,57)]
satl.shape <- subset(se.shape, State %in% c("DE", "DC", "FL", "GA", "MD", "NC", "SC", "VA", "WV"))
esc.shape <- subset(se.shape, State %in% c("AL", "KY", "MS", "TN"))
wsc.shape <- subset(se.shape, State %in% c("AR", "LA", "OK", "TX"))

#The FIPS 
fips <- county.fips

#The Polygons
world <- map_data("world")
states <- map_data("state")
counties <- map_data("county")

#The Join
counties$polyname <- paste(counties$region, counties$subregion, sep = ",")
counties <- counties %>% left_join(fips, by = c("polyname" = "polyname"))
counties$fips <- as.character(counties$fips)
counties <- counties %>% left_join(allpov, by = c("fips" = "FIPS"))

#Subsets
southern_states <- subset(states, region %in% 
                            c("texas", "arkansas", "louisiana", "mississippi", 
                              "alabama", "georgia", "florida", "north carolina",
                              "south carolina", "tennessee", "oklahoma", 
                              "kentucky", "west virginia", "virginia", 
                              "maryland", "delaware", "district of columbia"))

southern_counties <- subset(counties, region %in% 
                              c("texas", "arkansas", "louisiana", "mississippi", 
                                "alabama", "georgia", "florida", "north carolina",
                                "south carolina", "tennessee", "oklahoma", 
                                "kentucky", "west virginia", "virginia", 
                                "maryland", "delaware", "district of columbia"))
subset <- c("texas", "arkansas", "louisiana", "mississippi", 
            "alabama", "georgia", "florida", "north carolina",
            "south carolina", "tennessee", "oklahoma", 
            "kentucky", "west virginia", "virginia", 
            "maryland", "delaware", "district of columbia")

#Get FIPs
all_counties <- map(database = "county", region = subset, fill=T, plot=F)
all_IDs <- gsub(".*,","",all_counties$names)
fips.codes <- separate(data = fips, col = polyname, into = c("state", "county"), sep = ",")
county_IDs <- fips.codes$fips
county_names <- fips.codes$county
unique_IDs <- unique(county_IDs)

#create polygons
counties_sp <- map2SpatialPolygons(all_counties,
                                   all_IDs,
                                   CRS("+proj=longlat"),
                                   checkHoles = FALSE)
s#names(counties_sp@polygons) <- all_IDs

#Create neighbors
county.neighb.data<-poly2nb(counties_sp, queen=T)
names(neighb.data) <- names(counties_sp@polygons)

#Create list of neighbors
county.neighb.list<-nb2listw(county.neighb.data,style="W", zero.policy = TRUE) #, row.names(names(neighb.data))

#Create xy for all counties
county.xy<-centroid(counties_sp)
#? rownames(county.xy)<-county_IDs centroid id /= xy
colnames(all.xy)<-cbind("x","y")

#Create distance centroid
county.dist.k1 <-knn2nb(knearneigh(county.xy, k=1, longlat = TRUE))

#Determine max k distance value
county.max.k1 <-max(unlist(nbdists(county.dist.k1, county.xy, longlat=TRUE)))

#Calculate neighbors based on distance
county.sp.dist.k1 <-dnearneigh(county.xy, d1=0, d2=1 * county.max.k1, longlat = TRUE)

#Create neighbor list
county.dist.neighb.k1 <-nb2listw(county.sp.dist.k1,style="W", zero.policy = TRUE)

#spatial ERROR regression model based on distance
county.dist.err.k1 <-errorsarlm(child.pov.2016 ~ rural + urban + lnmanufacturing + lnag + lnretail + lnhealthss + lnconstruction + lnlesshs + lnunemployment + lnsinglemom + lnblack + lnhispanic + lnuninsured + lnincome_ratio + lnteenbirth + lnunmarried, data = allpov, listw = county.dist.neighb.k1)
summary(all.dist.err.k1)

all.r2 <- summary(all.dist.err.k1, correlation=TRUE, Nagelkerke = TRUE)


##############################################################
if (rgeosStatus()) {
  all_counties <- unionSpatialPolygons(all_counties, IDs = all_counties$names)
  }

counties %>% group_by(County_Name) %>% summarize(count=n())
