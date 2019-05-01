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

#append projection string to sites df
sites$PROJ4 = paste0('+proj=laea +lat_0=', sites$latitude,
    ' +lon_0=', sites$latitude)
equatorial_sites = sites$latitude <= 15 & sites$latitude >= -15
sites$PROJ4[equatorial_sites] =
    paste0('+proj=laea +lon_0=', sites$latitude)[equatorial_sites]

WSG84 = 4326
NAD83 = 4269

comid_from_point = function(lat, long, crs) {
    pt = st_point(c(long, lat))
    ptc = st_sfc(pt, crs=crs)
    COMID = discover_nhdplus_id(ptc)
    if(! length(COMID)) COMID = NA
    return(COMID)
}
sites$COMID = unlist(mapply(comid_from_point, sites$latitude,
    sites$longitude, WSG84))

vpu_from_point = function(lat, long, crs) {
    pt = st_point(c(long, lat))
    ptc = st_sfc(pt, crs=crs)
    VPU = find_vpu(ptc)
    return(VPU)
}
sites$VPU = unlist(mapply(vpu_from_point, sites$latitude,
    sites$longitude, WSG84))
sites$VPU[is.na(sites$COMID)] = NA

# nhd_plus_info(vpu = 4, component = "NHDSnapshot", dsn = "NHDWaterbody")
# nhd_plus_info(vpu = 1, component = "NHDPlusAttributes", dsn = "PlusFlow")
# nhd_plus_get(vpu = 1, component = "NHDPlusAttributes")
# nhd_plus_list(vpu = 1, component = "NHDSnapshot")
# nhd_plus_list(vpu = 1)

# z = nhd_plus_load(vpu=1, component='NHDPlusAttributes', dsn='PlusFlow',
function(df) {

    out = c()
    for(i in nrow(df)) {

        vpu = df$VPU[i]
        lat = df$latitude[i]
        long = df$longitude[i]
        comid = df$COMID[i]

        if(is.na(VPU)) {
            out = append(out, 'NA')
            message(paste0('Skipping row ', i, '. Not associated with NHD.'))
            next
        }

        message(paste0('Downloading NHDPlusV2 flowline components for VPU "',
            vpu, '". This only needs to happen once.'))
        fl = nhd_plus_load(vpu=vpu, component='NHDSnapshot',
            dsn='NHDFlowline', approve_all_dl=TRUE)
        fl_etc = nhd_plus_load(vpu=vpu, component='NHDPlusAttributes',
            dsn='PlusFlowlineVAA', approve_all_dl=TRUE)
        fl_etc = tryCatch(nhd_plus_load(vpu=vpu, component='NHDPlusAttributes',
            dsn='PlusFlowlineVAA', approve_all_dl=TRUE),
            error=function(e) {
                if(grepl('not found', e)) {
                    return('bug?')
                } else {
                    return(e)
                }
            })
        if('character' %in% class(fl_etc) && fl_etc == 'bug?') {
            nhdR_dir = paste0(Sys.getenv("XDG_DATA_HOME",  "~/.local/share"),
                '/nhdR/NHDPlus')
            nhdR_files = list.files(nhdR_dir, pattern='7z')
            # nhdR_dirs = list.dirs(nhdR_dir, full.names=FALSE, recursive=FALSE)
            bug = ! any(grepl(paste0(vpu, '_NHDPlusAttributes'), nhdR_files))
            if(bug) {
                message(paste('\nDisregard previous message.',
                    'Caught bug in nhdR package.'))
                vpu_num = str_match(vpu, '[0-9]+')[,1]
                to_mv = grepl(paste0(vpu_num, '[A-Z]?_NHDPlusAttributes'),
                    nhdR_files)
                file.remove(paste0(nhdR_dir, '/', nhdR_files[to_mv]))
                # file.rename(paste0(nhdR_dir, '/', nhdR_files[to_mv]),
                #     paste0(nhdR_dir, '/..'))
            }

        }
        # warning=function(w){print('aaa');print(w)})

        fl = fl[fl$ComID %in% comid,]
        fl = left_join(fl, fl_etc[, c('ComID', 'ToMeas', 'FromMeas')], by='ComID') %>%
            rename(COMID=ComID, REACHCODE=ReachCode)

        pt = st_point(c(-111.6690, 33.57110))
        pt = st_point(c(sites$longitude, sites$latitude))
        ptc = st_sfc(pt, crs=WSG84)
        ptct = st_transform(ptc, crs=NAD83)
        ptct = sf_project(ptc, crs=NAD83)
        q = get_flowline_index(z3, points=ptct) #REACH_meas: 100=upstream end; 0=down
    }

    mapview(z3) + mapview(ptct)

    subset = subset %>%
        st_as_sf(coords=c('longitude','latitude'), crs=4326) %>%
        st_transform(PROJ4)

    #get DEM, 14 is highest res, broadest area; 1 is converse
    dem = get_elev_raster(subset, z=8)
    mapview(dem) + mapview(subset)

    devtools::install_github("jsta/nhdR")

    #convert to spatial object and project from WGS 84
    # subset = subset %>%
    #     st_as_sf(coords=c('longitude','latitude'), crs=4326) %>%
    #     st_transform(PROJ4)

