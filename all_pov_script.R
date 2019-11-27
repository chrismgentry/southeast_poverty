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
se.shape<-readOGR(dsn="./Shapefiles",layer="SE_Counties_2016", stringsAsFactors = FALSE)
se.shape@data <- se.shape@data[-c(6:54,57:63)]
names(se.shape@polygons) <- se.shape@data$FIPS
satl.shape <- subset(se.shape, STATE_NAME %in% c("Delaware", "District of Columbia", "Florida", 
                                                 "Georgia", "Maryland", "North Carolina", "South Carolina",
                                                 "Virginia", "West Virginia"))
esc.shape <- subset(se.shape, STATE_NAME %in% c("Alabama", "Kentucky", "Mississippi",
                                           "Tennessee"))
wsc.shape <- subset(se.shape, STATE_NAME %in% c("Arkansas", "Louisiana", "Oklahoma",
                                           "Texas"))

#The FIPS 
fips <- county.fips
fips$fips <- str_pad(fips$fips, 5, "left", pad = 0)

#The Polygons
world <- map_data("world")
states <- map_data("state")
counties <- map_data("county")

#The Join
counties$polyname <- paste(counties$region, counties$subregion, sep = ",")
counties <- counties %>% left_join(fips, by = c("polyname" = "polyname"))
counties <- counties %>% left_join(se.data, by = c("fips" = "FIPS"))

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

#Create neighbors
se.neighbors<-poly2nb(se.shape, queen=T, row.names = se.shape$FIPS)
names(se.neighbors) <- names(se.shape$FIPS)
esc.neighbors<-poly2nb(esc.shape, queen=T, row.names = esc.shape$FIPS)
names(esc.neighbors) <- names(esc.shape$FIPS)
satl.neighbors<-poly2nb(satl.shape, queen=T, row.names = satl.shape$FIPS)
names(satl.neighbors) <- names(satl.shape$FIPS)
wsc.neighbors<-poly2nb(wsc.shape, queen=T, row.names = wsc.shape$FIPS)
names(wsc.neighbors) <- names(wsc.shape$FIPS)

#Create list of neighbors
se.neighbors.list<-nb2listw(se.neighbors, style="W", zero.policy = TRUE, 
                            row.names(names(se.neighbors)))
esc.neighbors.list<-nb2listw(esc.neighbors, style="W", zero.policy = TRUE, 
                            row.names(names(esc.neighbors)))
satl.neighbors.list<-nb2listw(satl.neighbors, style="W", zero.policy = TRUE, 
                            row.names(names(satl.neighbors)))
wsc.neighbors.list<-nb2listw(wsc.neighbors, style="W", zero.policy = TRUE, 
                            row.names(names(wsc.neighbors)))
#Create xy for all counties
county.xy<-centroid(se.shape)
colnames(county.xy)<-cbind("x","y")
esc.xy<-centroid(esc.shape)
colnames(esc.xy)<-cbind("x","y")
satl.xy<-centroid(satl.shape)
colnames(satl.xy)<-cbind("x","y")
wsc.xy<-centroid(wsc.shape)
colnames(wsc.xy)<-cbind("x","y")

#Create distance centroid
county.k1 <-knn2nb(knearneigh(county.xy, k=1, longlat = TRUE))
esc.k5 <-knn2nb(knearneigh(esc.xy, k=5, longlat = TRUE))
satl.k2 <-knn2nb(knearneigh(satl.xy, k=2, longlat = TRUE))
wsc.k1 <-knn2nb(knearneigh(wsc.xy, k=1, longlat = TRUE))

#Determine max k distance value
county.max.k1 <-max(unlist(nbdists(county.k1, county.xy, longlat=TRUE)))
esc.max.k5 <-max(unlist(nbdists(esc.k5, esc.xy, longlat=TRUE)))
satl.max.k2 <-max(unlist(nbdists(satl.k2, satl.xy, longlat=TRUE)))
wsc.max.k1 <-max(unlist(nbdists(wsc.k1, wsc.xy, longlat=TRUE)))

#Calculate neighbors based on distance
county.dist.k1 <-dnearneigh(county.xy, d1=0, d2=1 * county.max.k1, longlat = TRUE)
esc.dist.k5 <-dnearneigh(esc.xy, d1=0, d2=1 * esc.max.k5, longlat = TRUE)
satl.dist.k2 <-dnearneigh(satl.xy, d1=0, d2=1 * satl.max.k2, longlat = TRUE)
wsc.dist.k1 <-dnearneigh(wsc.xy, d1=0, d2=1 * wsc.max.k1, longlat = TRUE)

#Create neighbor list
county.k1.neighbors <-nb2listw(county.dist.k1, style="W", zero.policy = TRUE)
esc.k5.neighbors <-nb2listw(esc.dist.k5, style="W", zero.policy = TRUE)
satl.k2.neighbors <-nb2listw(satl.dist.k2, style="W", zero.policy = TRUE)
wsc.k1.neighbors <-nb2listw(wsc.dist.k1, style="W", zero.policy = TRUE)

#Creating the equation
equation <- child.pov.2016 ~ rural + urban + lnmanufacturing + lnag + 
  lnretail + lnhealthss + lnconstruction + lnlesshs + lnunemployment + 
  lnsinglemom + lnblack + lnhispanic + lnuninsured + lnincome_ratio + 
  lnteenbirth + lnunmarried

options(scipen = 5)

#OLS
ols <- lm(equation, data=se.data)
summary(ols)

#Morans Test
# cont.morans <- lm.morantest(ols, se.neighbors.list)
# cont.morans, did not produce significant results
dist.morans <- lm.morantest(ols, county.k1.neighbors)
dist.morans

#LaGrange Multiplier
# cont.lm.tests <- lm.LMtests(ols, se.neighbors.list, test="all")
# cont.lm.tests, did not produce significant results
dist.lm.tests <- lm.LMtests(ols, county.k1.neighbors, test="all")
dist.lm.tests

#Distance Lag Model
dist.lag.model <- spatialreg::lagsarlm(equation, data=se.data, 
                                       county.k1.neighbors)
dist.lag.summary <- summary(dist.lag.model, Nagelkerke = TRUE)

dist.lag.data <- cbind.data.frame(se.data$FIPS,
                                  dist.lag.summary$fitted.values,
                                  dist.lag.summary$residual,
                                  se.data$child.pov.2016,
                                  se.data$urban,
                                  se.data$lnretail,
                                  se.data$lnhealthss,
                                  se.data$lnconstruction,
                                  se.data$lnlesshs,
                                  se.data$lnunemployment,
                                  se.data$lnsinglemom, 
                                  se.data$lnhispanic,
                                  se.data$lnuninsured, 
                                  se.data$lnincome_ratio,
                                  se.data$lnunmarried,
                                    stringsAsFactors = FALSE)

#Renaming columns
colnames(dist.lag.data) <- c("fips","fitted","resid","childpov", "urban","retail",
                               "healthcare","construction","less_hs","unemployed",
                               "single_mom","hispanic","uninsured","income_ratio","unmarried")

#quantiles
quantiles_fit <- dist.lag.data %>%
  pull(fitted) %>%
  quantile(probs = seq(0, 1, length.out = 4), na.rm = TRUE)

quantiles_pov <- dist.lag.data %>%
  pull(childpov) %>%
  quantile(probs = seq(0, 1, length.out = 4), na.rm = TRUE)

#ranks
fit_rank <- cut(dist.lag.data$fitted, 
               breaks= quantiles_fit, 
               labels=c("1", "2", "3"), 
               na.rm = TRUE, 
               include.lowest = TRUE)

pov_rank <- cut(dist.lag.data$childpov, 
                breaks= quantiles_pov, 
                labels=c("1", "2", "3"), 
                na.rm = TRUE,
                include.lowest = TRUE)

#Join ranks and combined column to dataset
dist.lag.data$fit_score <- as.numeric(fit_rank)
dist.lag.data$pov_score <- as.numeric(pov_rank)
dist.lag.data$fit_pov <- paste(as.numeric(dist.lag.data$fit_score), 
                                 "-", 
                                 as.numeric(dist.lag.data$pov_score))

#legend
legend_colors <- tibble(
  x = c(3,2,1,3,2,1,3,2,1),
  y = c(3,3,3,2,2,2,1,1,1),
  z = c("#574249", "#627f8c", "#64acbe", "#985356", "#ad9ea5", "#b0d5df", "#c85a5a", "#e4acac", "#e8e8e8"))

xlabel <- "Modeled,Low \u2192 High"
xlabel <- gsub(",", "\n", xlabel)
ylabel <- "Observed,Low \u2192 High"
ylabel <- gsub(",", "\n", ylabel)

legend <- ggplot(legend_colors, aes(x,y)) + 
  geom_tile(aes(fill=z)) + 
  theme_minimal() + theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  labs(x = xlabel, y = ylabel) + 
  scale_fill_identity() +
  ggtitle("Legend") +
  theme(axis.title.y = element_text(face = "italic", hjust = 0.5, size = 8)) +
  theme(axis.title.x = element_text(face = "italic", hjust = 0.5, size = 8)) +
  theme(plot.title = element_text(face="bold", hjust = 0.5, size = 10))

#attach data to counties
county.data <- southern_counties %>% 
  left_join(dist.lag.data, by = c("fips" = "fips")) %>% fortify

#attach colors
bivariate_color_scale <- tibble(
  "3 - 3" = "#574249", 
  "2 - 3" = "#627f8c",
  "1 - 3" = "#64acbe",
  "3 - 2" = "#985356",
  "2 - 2" = "#ad9ea5",
  "1 - 2" = "#b0d5df",
  "3 - 1" = "#c85a5a",
  "2 - 1" = "#e4acac",
  "1 - 1" = "#e8e8e8") %>%
  gather("group", "fill")

county.data <- county.data %>% 
  left_join(bivariate_color_scale, by = c("fit_pov" = "group"))

#fit map
fit_pov_map <- ggplot() + 
  geom_polygon(data = world, aes(x=long,y=lat, group=group), fill = "gray95", color = "gray30") +
  geom_polygon(data = states, aes(x=long,y=lat, group=group), fill = "gray", color = "gray30") +
  geom_polygon(data = county.data, aes(x=long, y=lat, group=group, fill = fill)) + 
  geom_polygon(data = southern_counties, aes(x=long,y=lat, group=group), fill = NA, color = "black", size = 0.05) +
  geom_polygon(data = southern_states, aes(x=long,y=lat, group=group), fill = NA, color = "white") +
  coord_map("conic", lat0 = 30, xlim=c(-106,-77), ylim=c(24.5,40.5)) +
  scale_fill_identity() +
  theme_grey() + theme(legend.position="bottom") + theme(legend.title.align=0.5) +
  theme(panel.background = element_rect(fill = 'deepskyblue'),
        panel.grid.major = element_line(colour = NA)) +
  labs(x = "Longitude", y = "Latitude", fill = "Child Poverty", 
       title = "Bivariate Map of Observed vs Modeled Poverty Values") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

#final map
final_map <- ggdraw() +
  draw_plot(fit_pov_map, x = 0, y = 0, width = 1, height = 1) +
  draw_plot(legend, x = 0.55, y = 0.075, width = 0.15, height = 0.25) 
final_map

#Regional Analyses and other models
se.cont.lag.model <- spatialreg::lagsarlm(equation, data=se.data,
                                          se.neighbors.list)
se.cont.lag.summary <- summary(se.cont.lag.model, Nagelkerke = TRUE)
se.cont.lag.summary

se.cont.err.model <- spatialreg::errorsarlm(equation, data=se.data,
                                          se.neighbors.list)
se.cont.err.summary <- summary(se.cont.err.model, Nagelkerke = TRUE)
se.cont.err.summary

dist.err.model <- spatialreg::errorsarlm(equation, data=se.data, 
                                       county.k1.neighbors)
dist.err.summary <- summary(dist.err.model, Nagelkerke = TRUE)
dist.err.summary

#Original distance error models + lag
#ESC

esc.dist.lag.model <- spatialreg::lagsarlm(equation, data=esc.data, 
                                       esc.k5.neighbors)
esc.dist.lag.summary <- summary(esc.dist.lag.model, Nagelkerke = TRUE)
esc.dist.lag.summary
esc.dist.err.model <- spatialreg::errorsarlm(equation, data=esc.data, 
                                             esc.k5.neighbors)
esc.dist.err.summary <- summary(esc.dist.err.model, Nagelkerke = TRUE)
esc.dist.err.summary

#SATL
satl.dist.lag.model <- spatialreg::lagsarlm(equation, data=satl.data, 
                                       satl.k2.neighbors)
satl.dist.lag.summary <- summary(satl.dist.lag.model, Nagelkerke = TRUE)
satl.dist.lag.summary
satl.dist.err.model <- spatialreg::errorsarlm(equation, data=satl.data, 
                                         satl.k2.neighbors)
satl.dist.err.summary <- summary(satl.dist.err.model, Nagelkerke = TRUE)
satl.dist.err.summary

#WSC
wsc.dist.lag.model <- spatialreg::lagsarlm(equation, data=wsc.data, 
                                       wsc.k1.neighbors)
wsc.dist.lag.summary <- summary(wsc.dist.lag.model, Nagelkerke = TRUE)
wsc.dist.lag.summary
wsc.dist.err.model <- spatialreg::errorsarlm(equation, data=wsc.data, 
                                         wsc.k1.neighbors)
wsc.dist.err.summary <- summary(wsc.dist.err.model, Nagelkerke = TRUE)
wsc.dist.err.summary

#additional distance models

#Create distance centroid
county.k1 <-knn2nb(knearneigh(county.xy, k=1, longlat = TRUE))
county.k2 <-knn2nb(knearneigh(county.xy, k=2, longlat = TRUE))
county.k3 <-knn2nb(knearneigh(county.xy, k=3, longlat = TRUE))
county.k4 <-knn2nb(knearneigh(county.xy, k=4, longlat = TRUE))
county.k5 <-knn2nb(knearneigh(county.xy, k=5, longlat = TRUE))
esc.k1 <-knn2nb(knearneigh(esc.xy, k=1, longlat = TRUE))
esc.k2 <-knn2nb(knearneigh(esc.xy, k=2, longlat = TRUE))
esc.k3 <-knn2nb(knearneigh(esc.xy, k=3, longlat = TRUE))
esc.k4 <-knn2nb(knearneigh(esc.xy, k=4, longlat = TRUE))
esc.k5 <-knn2nb(knearneigh(esc.xy, k=5, longlat = TRUE))
satl.k1 <-knn2nb(knearneigh(satl.xy, k=1, longlat = TRUE))
satl.k2 <-knn2nb(knearneigh(satl.xy, k=2, longlat = TRUE))
satl.k3 <-knn2nb(knearneigh(satl.xy, k=3, longlat = TRUE))
satl.k4 <-knn2nb(knearneigh(satl.xy, k=4, longlat = TRUE))
satl.k5 <-knn2nb(knearneigh(satl.xy, k=5, longlat = TRUE))
wsc.k1 <-knn2nb(knearneigh(wsc.xy, k=1, longlat = TRUE))
wsc.k2 <-knn2nb(knearneigh(wsc.xy, k=2, longlat = TRUE))
wsc.k3 <-knn2nb(knearneigh(wsc.xy, k=3, longlat = TRUE))
wsc.k4 <-knn2nb(knearneigh(wsc.xy, k=4, longlat = TRUE))
wsc.k5 <-knn2nb(knearneigh(wsc.xy, k=5, longlat = TRUE))

#Determine max k distance value
county.max.k1 <-max(unlist(nbdists(county.k1, county.xy, longlat=TRUE)))
county.max.k2 <-max(unlist(nbdists(county.k2, county.xy, longlat=TRUE)))
county.max.k3 <-max(unlist(nbdists(county.k3, county.xy, longlat=TRUE)))
county.max.k4 <-max(unlist(nbdists(county.k4, county.xy, longlat=TRUE)))
county.max.k5 <-max(unlist(nbdists(county.k5, county.xy, longlat=TRUE)))
esc.max.k1 <-max(unlist(nbdists(esc.k1, esc.xy, longlat=TRUE)))
esc.max.k2 <-max(unlist(nbdists(esc.k2, esc.xy, longlat=TRUE)))
esc.max.k3 <-max(unlist(nbdists(esc.k3, esc.xy, longlat=TRUE)))
esc.max.k4 <-max(unlist(nbdists(esc.k4, esc.xy, longlat=TRUE)))
esc.max.k5 <-max(unlist(nbdists(esc.k5, esc.xy, longlat=TRUE)))
satl.max.k1 <-max(unlist(nbdists(satl.k1, satl.xy, longlat=TRUE)))
satl.max.k2 <-max(unlist(nbdists(satl.k2, satl.xy, longlat=TRUE)))
satl.max.k3 <-max(unlist(nbdists(satl.k3, satl.xy, longlat=TRUE)))
satl.max.k4 <-max(unlist(nbdists(satl.k4, satl.xy, longlat=TRUE)))
satl.max.k5 <-max(unlist(nbdists(satl.k5, satl.xy, longlat=TRUE)))
wsc.max.k1 <-max(unlist(nbdists(wsc.k1, wsc.xy, longlat=TRUE)))
wsc.max.k2 <-max(unlist(nbdists(wsc.k2, wsc.xy, longlat=TRUE)))
wsc.max.k3 <-max(unlist(nbdists(wsc.k3, wsc.xy, longlat=TRUE)))
wsc.max.k4 <-max(unlist(nbdists(wsc.k4, wsc.xy, longlat=TRUE)))
wsc.max.k5 <-max(unlist(nbdists(wsc.k5, wsc.xy, longlat=TRUE)))

#Calculate neighbors based on distance
county.dist.k1 <-dnearneigh(county.xy, d1=0, d2=1 * county.max.k1, longlat = TRUE)
county.dist.k2 <-dnearneigh(county.xy, d1=0, d2=1 * county.max.k2, longlat = TRUE)
county.dist.k3 <-dnearneigh(county.xy, d1=0, d2=1 * county.max.k3, longlat = TRUE)
county.dist.k4 <-dnearneigh(county.xy, d1=0, d2=1 * county.max.k4, longlat = TRUE)
county.dist.k5 <-dnearneigh(county.xy, d1=0, d2=1 * county.max.k5, longlat = TRUE)
esc.dist.k1 <-dnearneigh(esc.xy, d1=0, d2=1 * esc.max.k1, longlat = TRUE)
esc.dist.k2 <-dnearneigh(esc.xy, d1=0, d2=1 * esc.max.k2, longlat = TRUE)
esc.dist.k3 <-dnearneigh(esc.xy, d1=0, d2=1 * esc.max.k3, longlat = TRUE)
esc.dist.k4 <-dnearneigh(esc.xy, d1=0, d2=1 * esc.max.k4, longlat = TRUE)
esc.dist.k5 <-dnearneigh(esc.xy, d1=0, d2=1 * esc.max.k5, longlat = TRUE)
satl.dist.k1 <-dnearneigh(satl.xy, d1=0, d2=1 * satl.max.k1, longlat = TRUE)
satl.dist.k2 <-dnearneigh(satl.xy, d1=0, d2=1 * satl.max.k2, longlat = TRUE)
satl.dist.k3 <-dnearneigh(satl.xy, d1=0, d2=1 * satl.max.k3, longlat = TRUE)
satl.dist.k4 <-dnearneigh(satl.xy, d1=0, d2=1 * satl.max.k4, longlat = TRUE)
satl.dist.k5 <-dnearneigh(satl.xy, d1=0, d2=1 * satl.max.k5, longlat = TRUE)
wsc.dist.k1 <-dnearneigh(wsc.xy, d1=0, d2=1 * wsc.max.k1, longlat = TRUE)
wsc.dist.k2 <-dnearneigh(wsc.xy, d1=0, d2=1 * wsc.max.k2, longlat = TRUE)
wsc.dist.k3 <-dnearneigh(wsc.xy, d1=0, d2=1 * wsc.max.k3, longlat = TRUE)
wsc.dist.k4 <-dnearneigh(wsc.xy, d1=0, d2=1 * wsc.max.k4, longlat = TRUE)
wsc.dist.k5 <-dnearneigh(wsc.xy, d1=0, d2=1 * wsc.max.k5, longlat = TRUE)

#Create neighbor list
county.k1.neighbors <-nb2listw(county.dist.k1, style="W", zero.policy = TRUE)
county.k2.neighbors <-nb2listw(county.dist.k2, style="W", zero.policy = TRUE)
county.k3.neighbors <-nb2listw(county.dist.k3, style="W", zero.policy = TRUE)
county.k4.neighbors <-nb2listw(county.dist.k4, style="W", zero.policy = TRUE)
county.k5.neighbors <-nb2listw(county.dist.k5, style="W", zero.policy = TRUE)
esc.k1.neighbors <-nb2listw(esc.dist.k1, style="W", zero.policy = TRUE)
esc.k2.neighbors <-nb2listw(esc.dist.k2, style="W", zero.policy = TRUE)
esc.k3.neighbors <-nb2listw(esc.dist.k3, style="W", zero.policy = TRUE)
esc.k4.neighbors <-nb2listw(esc.dist.k4, style="W", zero.policy = TRUE)
esc.k5.neighbors <-nb2listw(esc.dist.k5, style="W", zero.policy = TRUE)
satl.k1.neighbors <-nb2listw(satl.dist.k1, style="W", zero.policy = TRUE)
satl.k2.neighbors <-nb2listw(satl.dist.k2, style="W", zero.policy = TRUE)
satl.k3.neighbors <-nb2listw(satl.dist.k3, style="W", zero.policy = TRUE)
satl.k4.neighbors <-nb2listw(satl.dist.k4, style="W", zero.policy = TRUE)
satl.k5.neighbors <-nb2listw(satl.dist.k5, style="W", zero.policy = TRUE)
wsc.k1.neighbors <-nb2listw(wsc.dist.k1, style="W", zero.policy = TRUE)
wsc.k2.neighbors <-nb2listw(wsc.dist.k2, style="W", zero.policy = TRUE)
wsc.k3.neighbors <-nb2listw(wsc.dist.k3, style="W", zero.policy = TRUE)
wsc.k4.neighbors <-nb2listw(wsc.dist.k4, style="W", zero.policy = TRUE)
wsc.k5.neighbors <-nb2listw(wsc.dist.k5, style="W", zero.policy = TRUE)

#Alter for each dist model
wsc.dist5.lag.model <- spatialreg::lagsarlm(equation, data=wsc.data, 
                                           wsc.k5.neighbors)
wsc.dist5.lag.summary <- summary(wsc.dist5.lag.model, Nagelkerke = TRUE)
wsc.dist5.lag.summary

wsc.dist5.err.model <- spatialreg::errorsarlm(equation, data=wsc.data, 
                                             wsc.k5.neighbors)
wsc.dist5.err.summary <- summary(wsc.dist5.err.model, Nagelkerke = TRUE)
wsc.dist5.err.summary

#Regional Cont Models
wsc.cont.lag.model <- spatialreg::lagsarlm(equation, data=wsc.data, 
                                       wsc.neighbors.list)
wsc.cont.lag.summary <- summary(wsc.cont.lag.model, Nagelkerke = TRUE)
wsc.cont.lag.summary

wsc.cont.err.model <- spatialreg::errorsarlm(equation, data=wsc.data, 
                                       wsc.neighbors.list)
wsc.cont.err.summary <- summary(wsc.cont.err.model, Nagelkerke = TRUE)
wsc.cont.err.summary

#Regular OLS
wsc.ols <- lm(equation, data=wsc.data)
summary(wsc.ols)
