library(sf)#
library(mapview)###
# library(mapedit)
# library(rayshader)
library(elevatr)#
library(raster)#
# devtools::install_github("giswqs/whiteboxR")
# library(whitebox) #
library(stars) #
# library(rgl)
# library(rgdal)
library(RMariaDB) #
library(stringr) #
library(dplyr) #
# knitr::knit_hooks$set(webgl = hook_webgl)

# setwd('~/git/streampulse/other_projects/watershed_data/')
setwd('~/Desktop/untracked/sp_watershed_delin')

#read in mysql pw
conf = readLines('/home/mike/git/streampulse/server_copy/sp/config.py')
# conf = readLines('/home/aaron/sp/config.py')
ind = which(lapply(conf, function(x) grepl('MYSQL_PW', x)) == TRUE)
pw = str_match(conf[ind], '.*\\"(.*)\\"')[2]

#read in site and results tables from mysql
con = dbConnect(RMariaDB::MariaDB(), dbname='sp', username='root', password=pw)
sites = as_tibble(dbReadTable(con, "site")) %>%
    arrange(latitude, longitude) %>%
    mutate(regionsite=paste(region, site, sep='_'))

# projection_table = list('equatorial'='+proj=laea +lon_0=-125.859375',
#     'temperate'='+proj=laea +lat_0=37.71859032558816 +lon_0=-86.484375',
#     'polar'='+proj=laea +lat_0=90.0 +lon_0=-117.0703125')

# transforms raster watershed outlines into shapefiles using stars
shed_stacker = function(x){
    read_stars(sheds[x]) %>%
        st_as_sf(merge=TRUE, use_integer=TRUE) %>%
        rename(id=1) %>%
        group_by(id) %>%
        summarize()
}

regions = unique(sites$region)
for(i in 1:length(regions)){

    subset = sites[sites$region == regions[i],
        c('regionsite', 'longitude', 'latitude')]

    #choose reasonable projection
    meanlat = mean(subset$latitude, na.rm=TRUE)
    meanlon = mean(subset$longitude, na.rm=TRUE)
    if(meanlat <= 15 && meanlat >= -15){ #equatorial
        PROJ4 = paste0('+proj=laea +lon_0=', meanlon)
    # } else if(meanlat >= 75 || meanlat <= -75){ #polar
    } else { #temperate or polar
        PROJ4 = paste0('+proj=laea +lat_0=', meanlat, ' +lon_0=', meanlon)
    }

    #convert to spatial object and project from WGS 84
    subset = subset %>%
        st_as_sf(coords=c('longitude','latitude'), crs=4326) %>%
        st_transform(PROJ4)

    #get DEM, 14 is highest res, broadest area; 1 is converse
    dem = get_elev_raster(subset, z=8)
    mapview(dem) + mapview(subset)

    #generate a box and check topo basemap for full watershed capture
    box = st_bbox(dem) %>% st_as_sfc()

    # Double check that we captured the whole watershed
    mapview(box) +  mapview(subset)

    #Save files so that whitebox can call the data
    writeRaster(dem, filename='dem.tif', overwrite=TRUE)
    st_write(subset, 'sites.shp', delete_layer=TRUE)

    #breach concavities using lindsay's algorithm (better than depression filling)
    # breach_depressions('./breach.tif', './filled.tif', fill_pits=TRUE)
    whitebox::breach_single_cell_pits('./dem.tif', './breach1.tif')
    whitebox::breach_depressions('./breach1.tif', './breached.tif')

    #accumulate flow and partial (inf) flow
    whitebox::d8_pointer('./breached.tif', './d8_pntr.tif')
    whitebox::d8_flow_accumulation('./breached.tif', './d8_flow.tif',
        out_type='catchment area')
    whitebox::d_inf_flow_accumulation('./breached.tif', './d_inf_flow.tif',
        out_type='catchment area')

    #snap_points, isolate stream grid cells
    whitebox::snap_pour_points('./sites.shp', './d8_flow.tif',
        './snapped_sites.shp', 100)
    whitebox::extract_streams('./d8_flow.tif', output='./streams.tif',
        threshold=11)

    #Watershed delineation as "whole watersheds'
    whitebox::unnest_basins('./d8_pntr.tif', './snapped_sites.shp',
        './basins/sheds.tif')

    #Read in flow accumulation algorithm
    fac = raster('./d8_flow.tif')

    #Get a list of the watershed created by `unnest_basins`
    sheds = list.files('./basins', full.names=TRUE)

    #convert all raster watershed outlines to shapefiles
    s = purrr::map(1:length(sheds), shed_stacker)

    #Use do.call to bind these sf objects into a single one
    shape_sheds = do.call('rbind', s) %>% arrange(id)

    #subset flow accumulation
    fac_sub = crop(fac, shape_sheds)

    mapview(fac_sub) + mapview(shape_sheds) + mapview(subset)

    #crop DEM and add watershed shapes to it?
    dem_crop = dem %>%
        crop(., shape_sheds) %>%
        mask(., shape_sheds)
}

#NLCD via FedData package ####
library(FedData)
nlcd = get_nlcd(dem, 'test', year=2011, raw.dir='.')

nlcd_sub <- nlcd %>%
    crop(., st_transform(shape_sheds, st_crs(nlcd))) %>%
    mask(., st_transform(shape_sheds, st_crs(nlcd)))

nlcd_key <- tibble(value=c(11,12,21,22,23,24,31,41,42,43,51,52,71,72,73,74,81,82,90,95),
    label=c('Open Water','Ice/Snow','Dev','Dev','Dev3','Dev4',
        'Barren','Dec.Forest','Evergreen.Forest','Mixed.Forest',
        'DwarfScrub','Shrub/Scrub','Grassland','Sedge','Lichens',
        'Moss','Pasture','Crops','WoodyWetlands','Wetlands'),
    colors=c('lightblue','White','red','red','red2','red3','brown4',
        'green4','green4','green4','brown','tan','green2',
        'darkolivegreen','green','blue','yellow3','orange3',
        'blue','blue'))

land_data <- raster::as.data.frame(nlcd_sub, xy=TRUE) %>%
    rename(cover=3) %>%
    filter(!is.na(cover)) %>%
    mutate(cover=as.character(cover))

order <- getValues(nlcd_sub) %>% unique()

sub_key <- nlcd_key %>%
    filter(value %in% order) %>%
    arrange(value)

mapview(nlcd_sub)

raster::image(nlcd_sub, col=sub_key$colors, axes = FALSE, xlab='', ylab='')

ggplot() +
    geom_raster(data=land_data, aes(x=x, y=y, fill=cover)) +
    coord_quickmap() +
    scale_fill_manual(labels=sub_key$label, values=sub_key$colors)

#other options ####
# devtools::install_github("USGS-R/nhdplusTools")
library(nhdplusTools)
