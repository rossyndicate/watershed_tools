library(sf)
library(sp)
library(mapview)
library(elevatr)
library(rgdal)
library(tidyverse)
library(streamstats)
library(geoknife)
# devtools::install_github("giswqs/whiteboxR")
# library(whitebox)
# devtools::install_github('markwh/streamstats')

# 1. setup ####

setwd('~/Desktop/untracked/watershed_geohax/') #SET PATH HERE
unlink('scratch', recursive=TRUE)
unlink('out', recursive=TRUE)

choose_projection = function(lat, long){

    if(lat <= 15 && lat >= -15){ #equatorial
        PROJ4 = paste0('+proj=laea +lon_0=', long)
    } else { #temperate or polar
        PROJ4 = paste0('+proj=laea +lat_0=', lat, ' +lon_0=', long)
    }

    return(PROJ4)
}

sites = read.csv('site_data.csv')
site = sites[5,] #this example is for one site; could easily be loopified
rownames(site) = 1:nrow(site)

PROJ4 = choose_projection(mean(site$latitude, na.rm=TRUE),
    mean(site$longitude, na.rm=TRUE))

#convert to spatial object and project from WGS84
WGS84_crs = 4326 #EPSG code for coordinate reference system
site_sf = site %>%
    sf::st_as_sf(coords=c('longitude', 'latitude'), crs=WGS84_crs) %>%
    sf::st_transform(PROJ4)

#get DEM, 1 is lowest res, broadest area; 14 is highest res
dem = elevatr::get_elev_raster(site_sf, z=10, expand=20000)
mapview(dem) + mapview(site_sf)

# 2. watershed delineation: legit method ####

#delineate watershed with streamstats package
z = streamstats::delineateWatershed(site$longitude, site$latitude,
    crs=WGS84_crs)

#save watershed boundary as shapefile
#streamstats::toSp and streamstats::writeShapefile are broken;
#the bodies of those functions are extracted and modified below:
tpf = tempfile(fileext='.geojson')
streamstats::writeGeoJSON(z, file=tpf, what='boundary')
spatialdf = rgdal::readOGR(tpf)
unlink(tpf)
dir.create('out', showWarnings=FALSE)
rgdal::writeOGR(spatialdf, dsn='out', layer='wsboundary',
    driver='ESRI Shapefile')

mapview(dem) + mapview(site_sf) + mapview(spatialdf)
streamstats::leafletWatershed(z)


# 3. summarize raster data across watershed area using geoknife package ####

#locate a spatial dataset
webdatasets = query('webdata')
title(webdatasets)
set_num = grep('evapotranspiration', abstract(webdatasets))
set_num = grep('evapotranspiration', title(webdatasets))
title(webdatasets[set_num])
spset = webdata(webdatasets[set_num[2]])

#set extents of dataset that will be retrieved
query(spset, 'variables')
variables(spset) = 'et'
query(spset, 'times')
times(spset) = as.POSIXct(c("2014-01-01", "2015-01-01"))

#summarize by watershed area delineated above
wsboundary = rgdal::readOGR(dsn='out', layer='wsboundary')

pgon = sp::SpatialPolygons(wsboundary@polygons,
        # proj4string=sp::CRS(PROJ4)) %>%
        #only one CRS is accepted? should probs reproject wsboundary tp WGS84
        proj4string=sp::CRS('+proj=longlat +datum=WGS84')) %>%
    simplegeom()
spset_job = geoknife(stencil=pgon, fabric=spset, wait=TRUE)
spset_cutout = result(spset_job, with.units=TRUE)
check(spset_job)
head(spset_cutout)

#summarize by point
station = as.data.frame(t(site[,c('longitude', 'latitude')]))
colnames(station) = paste(site$region, site$sitecode, sep='_')
station = simplegeom(station)

spset_job = geoknife(stencil=station, fabric=spset, wait=TRUE)
spset_summ = result(spset_job, with.units=TRUE)
check(spset_job)
head(spset_summ)


# 0. rubbish: formerly great delineation methods that no longer work ####

# #bounding box must encompass full watershed (you'll know for sure later)
# box = st_bbox(dem, crs=WGS84_crs) %>% sf::st_transform(PROJ4) %>% st_as_sfc()
# box = st_bbox(dem) %>% st_as_sfc()
# mapview(dem) + mapview(box) +  mapview(site_sf)

#Save files so that whitebox can access them
dir.create('scratch', showWarnings=FALSE)
writeRaster(dem, filename='scratch/dem.tif', overwrite=TRUE)
st_write(site_sf, 'scratch/sites.shp', delete_layer=TRUE)

#Fill single cell pits (for hydrologic correctness)
wbt_fill_single_cell_pits('scratch/dem.tif', 'scratch/breach1.tif')

#uses lindasy's algorithm (better than depression filling; BROKEN)
wbt_breach_depressions('scratch/breach1.tif', 'scratch/breach2.tif',
# wbt_breach_depressions('scratch/dem.tif', 'scratch/breach2.tif',
    flat_increment=1)

#Get flow direction and accumulation
wbt_d8_pointer('scratch/breach2.tif', 'scratch/d8_pntr.tif')
wbt_d8_flow_accumulation('scratch/d8_pntr.tif', 'scratch/d8_flow.tif',
    out_type='catchment area')
# wbt_d_inf_pointer('scratch/breach2.tif', 'scratch/d_inf_pntr.tif')
# wbt_d_inf_flow_accumulation('scratch/breach2.tif', 'scratch/d_inf_flow.tif',
#     out_type='catchment area')

#snap_points to talweg
wbt_snap_pour_points('scratch/sites.shp', 'scratch/d8_flow.tif',
    'scratch/snapped_sites.shp', snap_dist=5000)
snapped = rgdal::readOGR('scratch/snapped_sites.shp')
mapview(dem) + mapview(site_sf) + mapview(snapped)

#Watershed delineation as "whole watersheds'
dir.create('scratch/basins')

# wbt_watershed('scratch/d8_pntr.tif', 'scratch/snapped_sites.shp', 'scratch/wshed.tif')
# ws = raster('scratch/wshed.tif')
# mapview(ws)
wbt_unnest_basins('scratch/d8_pntr.tif', 'scratch/snapped_sites.shp',
    'scratch/basins/sheds.tif')

#Get a list of the watersheds created by `unnest_basins`
sheds = list.files('scratch/basins', full.names=TRUE)
sh1 = raster(sheds[1])
mapview(sh1)

#use stars package to transform raster watershed outlines into shapefiles
shed_stacker = function(x){
    read_stars(sheds[x]) %>%
        st_as_sf(merge=TRUE, use_integer=TRUE) %>%
        rename(id=1) %>%
        group_by(id) %>%
        summarize()
}

# Use purrr::map to apply the raster-shapefile transformation to all
# rasters to a list of shapefiles (map_dfr doesn't play nice with sf for
# unknown reasons)
s = purrr::map(1:length(sheds), shed_stacker)
mapview(s)

#Use do.call to bind these sf objects into a single one
shape_sheds = do.call('rbind', s) %>% arrange(id)

#read in and subset flow accumulation raster
# facc = raster('scratch/d_inf_flow.tif')
facc = raster('scratch/d8_flow.tif') %>%
    crop(shape_sheds)
mapview(facc)

#crop and mask elevatr DEM
only = dem %>%
    crop(., shape_sheds) %>%
    mask(., shape_sheds)

#something is wonky

# mapview(dem) + mapview(site_sf) + mapview(s) + mapview(only)
mapview(only) + mapview(site_sf)
mapview(site_sf)
mapview(box)
mapview(s)
mapview(shape_sheds)
mapview(facc_sub)
mapview(only)

# #Convert to matrix so rayshader is happy
# pmat <- matrix(raster::extract(only,raster::extent(only),buffer=300),
#     nrow=ncol(only),ncol=nrow(only))
#
# #Generate a hillshade
# raymat = ray_shade(pmat,sunangle=330)
#
# #use rayshader commands to generate map
# #rglwidget embeds output in html
# pmat %>%
#     sphere_shade(texture='desert') %>%
#     add_shadow(raymat) %>%
#     plot_3d(pmat,zscale=10,fov=0,theta=135,zoom=0.75,phi=45,
#         windowsize=c(750,750))
# rglwidget()
