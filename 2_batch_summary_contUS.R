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
library(RCurl)

setwd('~/Desktop/untracked/sp_watershed_delin') #this can be any empty directory

#read in mysql pw
conf = readLines('/home/mike/git/streampulse/server_copy/sp/config.py')
# conf = readLines('/home/aaron/sp/config.py')
ind = which(lapply(conf, function(x) grepl('MYSQL_PW', x)) == TRUE)
pw = str_match(conf[ind], '.*\\"(.*)\\"')[2]

#read in site table from mysql
con = dbConnect(RMariaDB::MariaDB(), dbname='sp', username='root', password=pw)
sites = as_tibble(dbReadTable(con, "site")) %>%
    arrange(latitude, longitude) %>%
    mutate(regionsite=paste(region, site, sep='_'))

# saveRDS(sites, '~/Desktop/site_data.rds')
# sites = readRDS('site_data.rds')

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
# A value of 0 means upstream end; 1 means downstream end.
calc_reach_props = function(df) {

    message(paste0('The nhdR package downloads NHDPlusV2 components to ',
        nhdR:::nhd_path(), '. Unfortunately this cannot be changed.',
        ' Fortunately, each component need only be downloaded once.'))

    out = c()
    for(i in 1:nrow(df)) {

        print(i)
        vpu = df$VPU[i]

        if(is.na(vpu)) {
            out = append(out, NA)
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
        out = append(out, 1 - x$REACH_meas / 100) #0=upstream end; 1=downstream end
    }

    return(out)
}
sites$reach_proportion = calc_reach_props(sites)

#acquire basic watershed and catchment characteristics from the nhdplusv2
merge_nhdplus_attr = function(df) {

    message(paste0('The nhdR package downloads NHDPlusV2 components to ',
        nhdR:::nhd_path(), '. Unfortunately this cannot be changed.',
        ' Fortunately, each component need only be downloaded once.'))

    out = data.frame()
    no_nhd_sites_yet = TRUE
    pre_nhd_site_counter = 0 #for prepending NA rows
    for(i in 1:nrow(df)) {

        print(i)
        vpu = df$VPU[i]

        if(is.na(vpu)) {
            if(no_nhd_sites_yet){
                pre_nhd_site_counter = pre_nhd_site_counter + 1
            } else {
                na_rows = as.data.frame(matrix(NA, nrow=1,
                    ncol=ncol(out), dimnames=list(NULL, colnames(out))))
                out = rbind(na_rows, out)
            }
            message(paste0('Skipping row ', i, '. Not associated with NHD.'))
            next
        } else {
            no_nhd_sites_yet = FALSE
        }

        fl = nhd_plus_load(vpu=vpu, component='NHDSnapshot',
            dsn='NHDFlowline', approve_all_dl=TRUE)
        fl_etc = nhd_plus_load(vpu=vpu, component='NHDPlusAttributes',
            dsn='PlusFlowlineVAA', approve_all_dl=TRUE)

        colnames(fl)[colnames(fl) == 'ComID'] = 'COMID'
        colnames(fl)[colnames(fl) == 'ReachCode'] = 'REACHCODE'
        fl = fl[fl$COMID %in% df$COMID[i],]
        # fl = left_join(fl, fl_etc[, c('ComID', 'ToMeas', 'FromMeas')],
        fl = left_join(fl[, c('COMID', 'REACHCODE')], fl_etc,
            by=c('COMID'='ComID'))
        out = rbind(out, as.data.frame(fl))
    }

    if(nrow(out) > 0){
        na_rows = as.data.frame(matrix(NA, nrow=pre_nhd_site_counter,
            ncol=ncol(out), dimnames=list(NULL, colnames(out))))
        out = rbind(na_rows, out)
        out = cbind(df, out)
    }

    return(out)
}
sites = merge_nhdplus_attr(sites)

#acquire streamcat catchment watershed characteristics
merge_streamcat_attr = function(df, storage_path){

    out = data.frame()
    no_nhd_sites_yet = TRUE
    pre_nhd_site_counter = 0 #for prepending NA rows
    for(i in 1:nrow(df)) {

        print(i)
        vpu = df$VPU[i]

        if(is.na(vpu)) {
            if(no_nhd_sites_yet){
                pre_nhd_site_counter = pre_nhd_site_counter + 1
            } else {
                na_rows = as.data.frame(matrix(NA, nrow=1,
                    ncol=ncol(out), dimnames=list(NULL, colnames(out))))
                out = rbind(na_rows, out)
            }
            message(paste0('Skipping row ', i, '. Not associated with NHD.'))
            next
        } else {
            no_nhd_sites_yet = FALSE
        }


        ####PULL STREAMCAT DATA BY STATE, STORE ON DISK
        #APPEND TO GIANT 1000-COLUMN CSV? OR KEEP SEPARATE?
        #PROBS THE LATTER
        USstate = df$region[i]
        ftpdir = paste0('ftp://newftp.epa.gov/EPADataCommons/ORD/',
            'NHDPlusLandscapeAttributes/StreamCat/States/')
        existing_files = list.files(storage_path)
        if(USstate %in% ){
        }
        # table = 'PredictedBioCondition'

        #get streamcat files associated with US state
        url_list = getURL(ftpdir, dirlistonly=TRUE)
        url_list = strsplit(url_list, split='\n')[[1]]
        url_list = url_list[grep(USstate, url_list)]

        #Loop through files on FTP, download, and append
        for(i in 1:length(url_list)){
            print(i)
            temp = tempfile()
            download.file(paste0(ftpdir, url_list[i]), temp)
            #replace .zip with .csv in file name
            csv_file = gsub('.zip', '.csv', url_list[i])
            #read in csv
            tmp_metric = read.csv(unz(temp, csv_file))
            tmp_metric = read.csv(temp, csv_file)
            #Append to final table
            if(i == 1){
                outdf = tmp_metric
            }else{
                outdf = rbind(outdf, tmp_metric)
            }
        }



        fl = nhd_plus_load(vpu=vpu, component='NHDSnapshot',
            dsn='NHDFlowline', approve_all_dl=TRUE)
        fl_etc = nhd_plus_load(vpu=vpu, component='NHDPlusAttributes',
            dsn='PlusFlowlineVAA', approve_all_dl=TRUE)

        colnames(fl)[colnames(fl) == 'ComID'] = 'COMID'
        colnames(fl)[colnames(fl) == 'ReachCode'] = 'REACHCODE'
        fl = fl[fl$COMID %in% df$COMID[i],]
        # fl = left_join(fl, fl_etc[, c('ComID', 'ToMeas', 'FromMeas')],
        fl = left_join(fl[, c('COMID', 'REACHCODE')], fl_etc,
            by=c('COMID'='ComID'))
        out = rbind(out, as.data.frame(fl))
    }

    if(nrow(out) > 0){
        na_rows = as.data.frame(matrix(NA, nrow=pre_nhd_site_counter,
            ncol=ncol(out), dimnames=list(NULL, colnames(out))))
        out = rbind(na_rows, out)
        out = cbind(df, out)
    }

    return(out)
}
# storage_path='/home/mike/.local/share/StreamCat'

#adjust catchment and watershed area for reach proportion
sites$AreaSqKM = sites$AreaSqKM * sites$reach_proportion
    #make a function to do this generically?




# subset = subset %>%
#     st_as_sf(coords=c('longitude', 'latitude'), crs=4326) %>%
#     st_transform(PROJ4)
#
# #get DEM, 14 is highest res, smallest area; 1 is lowest res, broadest area
# dem = get_elev_raster(subset, z=8)
# mapview(dem) + mapview(subset)
#
# devtools::install_github("jsta/nhdR")
#
# #convert to spatial object and project from WGS 84
# # subset = subset %>%
# #     st_as_sf(coords=c('longitude','latitude'), crs=4326) %>%
# #     st_transform(PROJ4)
#

