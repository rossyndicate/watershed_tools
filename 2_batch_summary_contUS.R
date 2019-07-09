library(sf)#
library(mapview)#
# library(elevatr)#
library(raster)#
# devtools::install_github("giswqs/whiteboxR")
# library(stars) #
library(RMariaDB) #
library(stringr) #
library(dplyr) #
library(plyr) #
# devtools::install_github("USGS-R/nhdplusTools")
library(nhdplusTools)
# devtools::install_github("jsta/nhdR")
library(nhdR)
# knitr::knit_hooks$set(webgl = hook_webgl)
library(RCurl)
library(MODISTools)

setwd('~/git/watershed_tools/') #this can be any empty directory

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
# sites = readRDS('site_data.rds')
sites = readRDS('site_table_for_ws.rds')
sites = sites[sites$by > 0,]
sites = sites[sites$site != 'PILG1',]

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

sites = sites[! is.na(sites$COMID),]
# write.csv(sites, '~/git/watershed_tools/sites_with_comid.csv', row.names=FALSE)

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

#does not need to be a separate function. just use nhdplusv2_from_comid
#and calculate reach proportion from the output
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

####NOTE: calc_reach_props fills up system memory, so use the following code to
####save sites dataframe, then close R, reopen it, and start from here.
# write.csv(sites, row.names=FALSE,
#     '/home/mike/Dropbox/streampulse/data/watershed_data/sites_with_reachprop.csv')
# sites = read.csv('/home/mike/Dropbox/streampulse/data/watershed_data/sites_with_reachprop.csv',
#     stringsAsFactors=FALSE)

#acquire basic watershed and catchment characteristics from the nhdplusv2
####TO-DO: figure out why PlusFlowAR contains no data
nhdplusv2_from_comid = function(VPU, COMID, component, DSN, silent=FALSE) {

    if(! silent){
        message(paste0('The nhdR package downloads NHDPlusV2 components to ',
            nhdR:::nhd_path(), '. Unfortunately this cannot be changed.',
            ' Fortunately, each component need only be downloaded once.'))
    }

    data = nhd_plus_load(vpu=VPU, component=component,
        dsn=DSN, approve_all_dl=TRUE)

    colnames(data)[colnames(data) == 'ComID'] = 'COMID'
    colnames(data)[colnames(data) == 'ReachCode'] = 'REACHCODE'
    data = data[data$COMID == COMID,]

    return(data)
}

nhdplus_sets = list(#'NHDSnapshot'='NHDFlowline',
    'NHDPlusAttributes'='PlusFlowlineVAA', #'NHDPlusAttributes'='PlusFlowAR',
    'NHDPlusAttributes'='ElevSlope')
nhdplus_data = data.frame()

for(j in 1:nrow(sites)){
    if(is.na(sites$COMID[j])) next

    for(i in 1:length(nhdplus_sets)){
        print(paste(j, nhdplus_sets[[i]]))

        if(i == 1 || initerr){
            row_base = try(nhdplusv2_from_comid(sites$VPU[j], sites$COMID[j],
                names(nhdplus_sets[i]), nhdplus_sets[[i]], silent=TRUE))
            if(class(row_base)[1] == 'try-error'){
                initerr = TRUE
                print('init error')
            } else {
                initerr = FALSE
            }
        } else {
            row_ext = try(nhdplusv2_from_comid(sites$VPU[j], sites$COMID[j],
                names(nhdplus_sets[i]), nhdplus_sets[[i]], silent=TRUE))
            if(class(row_ext)[1] == 'try-error'){
                print('error')
                next
            }
            row_base = left_join(row_base, row_ext)
        }

    }

    nhdplus_data = rbind.fill(nhdplus_data, row_base)
}

colnames(nhdplus_data) = toupper(colnames(nhdplus_data))
# lengthkm_dupe = which(colnames(nhdplus_data) == 'LENGTHKM')
# nhdplus_data = nhdplus_data[,-lengthkm_dupe[1]]
nhdplus_data = nhdplus_data[, ! duplicated(colnames(nhdplus_data))]
nhdplus_data = select(nhdplus_data, COMID, STREAMORDE, FROMMEAS, TOMEAS, SLOPE,
    REACHCODE, AREASQKM, TOTDASQKM, MAXELEVSMO, MINELEVSMO)
sites = left_join(sites, nhdplus_data, by='COMID')
sites = sites[! duplicated(sites$site),]
sites$AREASQKM_corr = round(sites$AREASQKM * sites$reach_proportion, 5)
sites$TOTDASQKM_corr = sites$TOTDASQKM - (sites$AREASQKM - sites$AREASQKM_corr)
sites$areal_corr_factor = sites$TOTDASQKM_corr / sites$TOTDASQKM

# write.csv(sites, row.names=FALSE,
#     '/home/mike/Dropbox/streampulse/data/watershed_data/sites_with_nhd.csv')
# sites = read.csv('/home/mike/Dropbox/streampulse/data/watershed_data/sites_with_nhd.csv',
#     stringsAsFactors=FALSE)

#find out which streamcat datasets are available
query_streamcat_datasets = function(keyword=NULL){

    ftpdir = paste0('ftp://newftp.epa.gov/EPADataCommons/ORD/',
        'NHDPlusLandscapeAttributes/StreamCat/States/')

    url_list = getURL(ftpdir, dirlistonly=TRUE)
    url_list = strsplit(url_list, split='\n')[[1]]

    if(! is.null(keyword)){
        url_list = url_list[grep(keyword, url_list, ignore.case=TRUE)]
    }

    return(url_list)
}
query_streamcat_datasets('ripbuf')

#acquire lots of diverse data from streamcat
streamcat_from_comid = function(USstate, COMID, dataset){

    ftpdir = paste0('ftp://newftp.epa.gov/EPADataCommons/ORD/',
        'NHDPlusLandscapeAttributes/StreamCat/States/')
    zip_name = paste0(dataset, '_', USstate, '.zip')

    csv_name = gsub('.zip', '.csv', zip_name)
    temp = tempfile()
    download.file(paste0(ftpdir, zip_name), temp)
    data = read.csv(unz(temp, csv_name), stringsAsFactors=FALSE)
    data = data[data$COMID == COMID,]

    return(data)
}

streamcat_sets = c('Elevation', 'PRISM_1981_2010', 'NLCD2011', 'Runoff',
    'STATSGO_Set2', 'NADP', 'GeoChemPhys1', 'GeoChemPhys2', 'BFI')
streamcat_data = data.frame()


for(j in 1:nrow(sites)){
    if(is.na(sites$COMID[j])) next
    initerr = FALSE

    for(i in 1:length(streamcat_sets)){
        print(paste(j, streamcat_sets[i]))

        if(i == 1 || initerr){
            # row_base = tryCatch(
            row_base = try(streamcat_from_comid(sites$region[j], sites$COMID[j],
                streamcat_sets[i]))
            if(class(row_base)[1] == 'try-error'){
                initerr = TRUE
                print('init error')
            } else {
                initerr = FALSE
            }
            #     error=function(e){print('fail')},
            #     warning=function(w){NULL},
            # )
        } else {
            row_ext = try(streamcat_from_comid(sites$region[j], sites$COMID[j],
                streamcat_sets[i]))
            if(class(row_ext)[1] == 'try-error'){
                print('error')
                next
            }
            row_base = left_join(row_base, row_ext)
        }

    }

    streamcat_data = rbind.fill(streamcat_data, row_base)
}

write.csv(streamcat_data, row.names=FALSE,
    '/home/mike/Dropbox/streampulse/data/watershed_data/streamcat.csv')
# streamcat_data = read.csv('/home/mike/Dropbox/streampulse/data/watershed_data/streamcat.csv',
#     stringsAsFactors=FALSE)
streamcat_data = select(streamcat_data, COMID, ElevWs, Precip8110Ws, Tmin8110Ws,
    Tmax8110Ws, Tmean8110Ws, RunoffWs, matches('^Pct[a-zA-z]+2011Ws$'),
    PermWs, RckDepWs, OmWs, WtDepWs, matches('^[a-zA-z0-9]_2008Ws$'), BFIWs,
    NWs, Al2O3Ws, CaOWs, Fe2O3Ws, K2OWs, MgOWs, Na2OWs, P2O5Ws, SWs, SiO2Ws)
streamcat_data = mutate(streamcat_data,
    precip_runoff_ratio=Precip8110Ws / RunoffWs)

sites = left_join(sites, streamcat_data, by='COMID')
sites = sites[! duplicated(sites$site),]


#acquire MODIS data (this section not yet begun)
# VNP13A1
mt_bands("MOD13Q1")
subset1 <- mt_subset(product = "MOD13Q1",
    lat = 40,
    lon = -110,
    band = "250m_16_days_NDVI",
    start = "2004-01-01",
    end = "2004-02-01",
    km_lr = 10,
    km_ab = 10,
    site_name = "testsite",
    internal = TRUE,
    progress = FALSE)

dfx <- data.frame("site_name" = paste("test",1:2))
dfx$lat <- 40
dfx$lon <- -110

# test batch download
subsets <- mt_batch_subset(dfx = dfx,
    product = "MOD11A2",
    band = "LST_Day_1km",
    internal = TRUE,
    start = "2004-01-01",
    end = "2004-02-01")


#adjust catchment and watershed area for reach proportion
sites$AreaSqKM = sites$AreaSqKM * sites$reach_proportion
    #make a function to do this generically?




#### SCRAPS ####

#this works, but i replaced it with a more rudimentary function that
#just acquires data and leaves postprocessing to the user
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

#acquire streamcat catchment watershed characteristics (batch; incomplete)
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

