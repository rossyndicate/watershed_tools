library(sf)#
library(mapview)###
# library(elevatr)#
library(raster)#
# devtools::install_github("giswqs/whiteboxR")
# library(stars) #
library(RMariaDB) #
library(stringr) #
library(dplyr) #
# devtools::install_github("USGS-R/nhdplusTools")
library(nhdplusTools)
# devtools::install_github("jsta/nhdR")
library(nhdR)
# knitr::knit_hooks$set(webgl = hook_webgl)

setwd('~/Desktop/untracked/sp_watershed_delin') #this can be any empty directory

# #read in mysql pw
# conf = readLines('/home/mike/git/streampulse/server_copy/sp/config.py')
# # conf = readLines('/home/aaron/sp/config.py')
# ind = which(lapply(conf, function(x) grepl('MYSQL_PW', x)) == TRUE)
# pw = str_match(conf[ind], '.*\\"(.*)\\"')[2]
#
# #read in site table from mysql
# con = dbConnect(RMariaDB::MariaDB(), dbname='sp', username='root', password=pw)
# sites = as_tibble(dbReadTable(con, "site")) %>%
#     arrange(latitude, longitude) %>%
#     mutate(regionsite=paste(region, site, sep='_'))

# saveRDS(sites, '~/Desktop/site_data.rds')
sites = readRDS('site_data.rds')

#append projection string to sites df
sites$PROJ4 = paste0('+proj=laea +lat_0=', sites$latitude,
    ' +lon_0=', sites$latitude)
equatorial_sites = sites$latitude <= 15 & sites$latitude >= -15
sites$PROJ4[equatorial_sites] =
    paste0('+proj=laea +lon_0=', sites$latitude)[equatorial_sites]

WSG84 = 4326
NAD83 = 4269

#COMID is the NHDPlus identifier for any reach. This gets it from lat/long
comid_from_point = function(lat, long, crs) {
    pt = st_point(c(long, lat))
    ptc = st_sfc(pt, crs=crs)
    COMID = discover_nhdplus_id(ptc)
    if(! length(COMID)) COMID = NA
    return(COMID)
}
sites$COMID = unlist(mapply(comid_from_point, sites$latitude,
    sites$longitude, WSG84))

#VPU is an NHDPlus processing unit. needed for the next function to work
vpu_from_point = function(lat, long, crs) {
    pt = st_point(c(long, lat))
    ptc = st_sfc(pt, crs=crs)
    VPU = find_vpu(ptc)
    return(VPU)
}
sites$VPU = unlist(mapply(vpu_from_point, sites$latitude,
    sites$longitude, WSG84))
sites$VPU[is.na(sites$COMID)] = NA

#this calculates how far along a reach any given point falls. That way when we pull in
#watershed summary data for a reach, we can adjust it according to how much
#of the total upstream area actually contributes to the point in question.
calc_reach_props = function(df) {

    message(paste0('The nhdR package downloads NHDPlusV2 components to ',
        nhdR:::nhd_path(), '. Unfortunately this cannot be changed.',
        ' Fortunately, each component need only be downloaded once.'))

    out = c()
    for(i in 1:nrow(df)) {

        print(i)
        vpu = df$VPU[i]

        if(is.na(vpu)) {
            out = append(out, 'NA')
            message(paste0('Skipping row ', i, '. Not associated with NHD.'))
            next
        }

        fl = nhd_plus_load(vpu=vpu, component='NHDSnapshot',
            dsn='NHDFlowline', approve_all_dl=TRUE)
        fl_etc = nhd_plus_load(vpu=vpu, component='NHDPlusAttributes',
            dsn='PlusFlowlineVAA', approve_all_dl=TRUE)

        colnames(fl)[colnames(fl) == 'ComID'] = 'COMID'
        colnames(fl)[colnames(fl) == 'ReachCode'] = 'REACHCODE'
        # fl = rename_if(fl, 'COMID'='ComID', 'REACHCODE'='ReachCode')
        fl = fl[fl$COMID %in% df$COMID[i],]
        fl = left_join(fl, fl_etc[, c('ComID', 'ToMeas', 'FromMeas')],
            by=c('COMID'='ComID'))

        pt = st_point(c(df$longitude[i], df$latitude[i]))
        ptc = st_sfc(pt, crs=WSG84)
        ptct = st_transform(ptc, crs=NAD83)
        x = suppressWarnings(get_flowline_index(fl, points=ptc))
        out = append(out, x$REACH_meas) #100=upstream end; 0=downstream end
    }

    return(out)
}

sites$reach_proportion = calc_reach_props(sites)

subset = subset %>%
    st_as_sf(coords=c('longitude','latitude'), crs=4326) %>%
    st_transform(PROJ4)

#get DEM, 14 is highest res, smallest area; 1 is lowest res, broadest area
dem = get_elev_raster(subset, z=8)
mapview(dem) + mapview(subset)

devtools::install_github("jsta/nhdR")

#convert to spatial object and project from WGS 84
# subset = subset %>%
#     st_as_sf(coords=c('longitude','latitude'), crs=4326) %>%
#     st_transform(PROJ4)

