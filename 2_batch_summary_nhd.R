#batch summarizing of watershed data for the continental USA
#Mike Vlah (vlahm13@gmail.com)

#use these tools to:
    #1. acquire NHDPlusV2 COMIDs based on lat/long
    #2. acquire NHDPlusV2 VPU IDs based on lat/long
    #3. determine the position of a site along its NHD reach and calculate reach
    #   proportion, then adjust linear and areal summary data accordingly
    #4. use COMID and VPU to acquire NHDPlusV2 data for your sites
    #5. use COMID to acquire StreamCat data for your sites

#see NHDPlusV2 docs (1) and StreamCat variable list (2) for help.
#1. ftp://ftp.horizon-systems.com/NHDplus/NHDPlusV21/Documentation/NHDPlusV2_User_Guide.pdf
#2. ftp://newftp.epa.gov/EPADataCommons/ORD/NHDPlusLandscapeAttributes/StreamCat/Documentation/VariableList-QuickReference.html

#NOTE: these tools necessarily download a lot of large datasets and store them
#in memory. keep an eye on your usage.

library(stringr)
library(dplyr)
library(plyr)
# remotes::install_github("USGS-R/nhdplusTools")
library(nhdplusTools)
# remotes::install_github("jsta/nhdR")
library(nhdR)
library(RCurl)
library(sf)
# library(MODISTools)

# 1. setup and helper functions ####

#set your working directory to the location of your site data file.
#   site data must include latitude and longitude columns (decimal degrees)
# setwd('~/Desktop/untracked/watershed_geohax/')
sites = read.csv('site_data.csv', stringsAsFactors=FALSE)
WGS84 = 4326 #EPSG code for coordinate reference system

comid_from_point = function(lat, long, CRS) {
    pt = sf::st_point(c(long, lat))
    ptc = sf::st_sfc(pt, crs=CRS)
    COMID = nhdplusTools::discover_nhdplus_id(ptc)
    if(! length(COMID)) COMID = NA
    return(COMID)
}

vpu_from_point = function(lat, long, CRS) {
    pt = sf::st_point(c(long, lat))
    ptc = sf::st_sfc(pt, crs=CRS)
    VPU = nhdR::find_vpu(ptc)
    return(VPU)
}

#this calculates how far along a reach any given point falls. That way when we pull in
#watershed summary data for a reach, we can adjust it according to how much
#of the total upstream area actually contributes to the point in question.
# A value of 0 means upstream end; 1 means downstream end.
calc_reach_prop = function(VPU, COMID, lat, long, CRS, quiet=FALSE){

    if(! quiet){
        message(paste0('The nhdR package downloads NHDPlusV2 components to ',
            nhdR:::nhd_path(), '. Unfortunately this cannot be changed.',
            ' Fortunately, each component need only be downloaded once.'))
    }

    fl = nhdR::nhd_plus_load(vpu=VPU, component='NHDSnapshot',
        dsn='NHDFlowline', approve_all_dl=TRUE)
    fl_etc = nhdR::nhd_plus_load(vpu=VPU, component='NHDPlusAttributes',
        dsn='PlusFlowlineVAA', approve_all_dl=TRUE)

    colnames(fl)[colnames(fl) == 'ComID'] = 'COMID'
    colnames(fl)[colnames(fl) == 'ReachCode'] = 'REACHCODE'
    fl = fl[fl$COMID == COMID,]
    fl = left_join(fl, fl_etc[, c('ComID', 'ToMeas', 'FromMeas')],
        by=c('COMID'='ComID'))

    pt = sf::st_point(c(long, lat))
    ptc = sf::st_sfc(pt, crs=CRS)
    ptct = sf::st_transform(ptc, crs=4269) #CRS for NAD 83
    x = suppressWarnings(nhdplusTools::get_flowline_index(fl, points=ptct))
    out = 1 - x$REACH_meas / 100 #0=upstream end; 1=downstream end

    return(out)
}

#this acquires nhdplusv2 data for a single site by COMID.
#it's just a thin wrapper around nhdR::nhd_plus_load
nhdplusv2_from_comid = function(VPU, COMID, component, DSN, quiet=FALSE) {

    if(! quiet){
        message(paste0('The nhdR package downloads NHDPlusV2 components to ',
            nhdR:::nhd_path(), '. Unfortunately this cannot be changed.',
            ' Fortunately, each component need only be downloaded once.'))
    }

    data = nhdR::nhd_plus_load(vpu=VPU, component=component,
        dsn=DSN, approve_all_dl=TRUE)

    colnames(data)[colnames(data) == 'ComID'] = 'COMID'
    colnames(data)[colnames(data) == 'ReachCode'] = 'REACHCODE'
    data = data[data$COMID == COMID,]

    return(data)
}

#this calls nhdplusv2_from_comid repeatedly to get data for all your sites.
#the dataframe must include COMID and VPU columns
nhdplusv2_bulk = function(site_df, nhdplusv2_sets, quiet=FALSE){

    nhdplus_data = data.frame()
    if(any(is.na(site_df$COMID))) stop('you have missing COMIDs')

    for(j in 1:nrow(site_df)){
        for(i in 1:length(setlist)){
            print(paste(j, nhdplusv2_sets[[i]]))

            if(i == 1 || initerr){
                row_base = try(nhdplusv2_from_comid(site_df$VPU[j],
                    site_df$COMID[j], names(setlist[i]), setlist[[i]],
                    quiet=quiet))
                if('try-error' %in% class(row_base) || nrow(row_base) > 1){
                    initerr = TRUE
                    row_base = data.frame(COMID=site_df$COMID[j])
                } else {
                    initerr = FALSE
                }
            } else {
                row_ext = try(nhdplusv2_from_comid(site_df$VPU[j],
                    site_df$COMID[j], names(setlist[i]), setlist[[i]],
                    quiet=quiet))
                if(! 'try-error' %in% class(row_ext) && nrow(row_ext) == 1){
                    row_base = left_join(row_base, row_ext)
                }
            }

        }

        if(nrow(row_base) > 1){
            row_base = data.frame(COMID=site_df$COMID[j])
        }

        nhdplus_data = rbind.fill(nhdplus_data, row_base)
    }

    return(nhdplus_data)
}

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

#this function acquires streamcat data for a single site by NHDPlusV2 COMID.
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

#this calls streamcat_from_comid repeatedly to get data for all your sites
#the dataframe must include COMID and region columns, where "region" refers to
#each state's 2-letter abbreviation.
streamcat_bulk = function(site_df, streamcat_sets){

    streamcat_data = data.frame()
    if(any(is.na(site_df$COMID))) stop('you have missing COMIDs')

    for(j in 1:nrow(site_df)){
        for(i in 1:length(streamcat_sets)){
            print(paste(j, streamcat_sets[i]))

            if(i == 1 || initerr){
                row_base = try(streamcat_from_comid(site_df$region[j],
                    site_df$COMID[j], streamcat_sets[i]))
                if('try-error' %in% class(row_base) || nrow(row_base) > 1){
                    initerr = TRUE
                    row_base = data.frame(COMID=site_df$COMID[j])
                } else {
                    initerr = FALSE
                }
            } else {
                row_ext = try(streamcat_from_comid(site_df$region[j],
                    site_df$COMID[j], streamcat_sets[i]))
                if(! 'try-error' %in% class(row_ext) && nrow(row_ext) == 1){
                    row_base = left_join(row_base, row_ext)
                }
            }

        }

        if(nrow(row_base) > 1){
            row_base = data.frame(COMID=site_df$COMID[j])
        }

        streamcat_data = rbind.fill(streamcat_data, row_base)
    }

    return(streamcat_data)
}


# 2. get NHDPlusV2 data ####

#COMID is the NHD identifier for any reach in the continental U.S.
#add COMIDs to your site table. If this doesn't work, update nhdplusTools
sites$COMID = unlist(mapply(comid_from_point, sites$latitude,
    sites$longitude, WGS84))
sites = sites[! is.na(sites$COMID),]

#VPU == NHD vector processing unit. NHDPlusV2 data are downloaded per VPU.
#add VPUs to your site table and determine reach proportions.
sites$VPU = unlist(mapply(vpu_from_point, sites$latitude,
    sites$longitude, WGS84))

sites$reach_proportion = unlist(mapply(calc_reach_prop, sites$VPU, sites$COMID,
    sites$latitude, sites$longitude, WGS84, quiet=TRUE))

#construct list of DSN=component pairs to acquire. see NHDPlus docs for more.
setlist = list('NHDPlusAttributes'='PlusFlowlineVAA',
    'NHDPlusAttributes'='ElevSlope')
    #'NHDSnapshot'='NHDFlowline', 'NHDPlusAttributes'='PlusFlowAR'

#retrieve NHDPlusV2 data
nhdplusv2_data = nhdplusv2_bulk(sites, setlist, quiet=TRUE)

#nhd variable names do not have consistent naming conventions. sometimes they're
#all caps; other times camel case. here's a crude way to deal with that.
colnames(nhdplusv2_data) = toupper(colnames(nhdplusv2_data))
nhdplusv2_data = nhdplusv2_data[, ! duplicated(colnames(nhdplusv2_data))]

#pick out the variables you want, then join them to your site data
nhdplusv2_data = select(nhdplusv2_data, COMID, STREAMORDE, FROMMEAS, TOMEAS, SLOPE,
    REACHCODE, AREASQKM, TOTDASQKM, MAXELEVSMO, MINELEVSMO)
sites = left_join(sites, nhdplusv2_data, by='COMID')
sites = sites[! duplicated(sites$sitecode),]

#correct catchment area (AREASQKM) based on where each site falls within its reach.
#use this to correct watershed area (TOTDASQKM) and to determine an areal
#correction factor that can be multiplied with any areal summary data.
sites$AREASQKM_corr = round(sites$AREASQKM * sites$reach_proportion, 5)
sites$TOTDASQKM_corr = sites$TOTDASQKM - (sites$AREASQKM - sites$AREASQKM_corr)
sites$areal_corr_factor = sites$TOTDASQKM_corr / sites$TOTDASQKM


# 3. get StreamCat data ####

#find out which streamcat datasets are available (case insensitive)
query_streamcat_datasets()
query_streamcat_datasets('ripbuf')

#construct vector of streamcat datasets to acquire (check variable list for deets)
setlist2 = c('Elevation', 'PRISM_1981_2010', 'NLCD2011', 'Runoff',
    'STATSGO_Set2', 'NADP', 'GeoChemPhys1', 'GeoChemPhys2', 'BFI')

streamcat_data = streamcat_bulk(sites, setlist2)

#pick out the variables you want, then join them to your site data
streamcat_data = select(streamcat_data, COMID, ElevWs, Precip8110Ws, Tmin8110Ws,
        Tmax8110Ws, Tmean8110Ws, RunoffWs, matches('^Pct[a-zA-z]+2011Ws$'),
        PermWs, RckDepWs, OmWs, WtDepWs, matches('^[a-zA-z0-9]_2008Ws$'), BFIWs,
        NWs, Al2O3Ws, CaOWs, Fe2O3Ws, K2OWs, MgOWs, Na2OWs, P2O5Ws, SWs, SiO2Ws) %>%
    mutate(precip_runoff_ratio=Precip8110Ws / RunoffWs)

sites = left_join(sites, streamcat_data, by='COMID')
sites = sites[! duplicated(sites$sitecode),]

#save yer data
sites = arrange(sites, region, sitecode)
write.csv(sites, 'watershed_summary_data.csv', row.names=FALSE)


# 4. get MODIS data (this section incomplete) ####
# VNP13A1
mt_bands("MOD13Q1")
subset1 = mt_subset(product = "MOD13Q1",
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

dfx = data.frame("site_name" = paste("test",1:2))
dfx$lat = 40
dfx$lon = -110

# test batch download
subsets = mt_batch_subset(dfx = dfx,
    product = "MOD11A2",
    band = "LST_Day_1km",
    internal = TRUE,
    start = "2004-01-01",
    end = "2004-02-01")


#### SCRAPS ####

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

