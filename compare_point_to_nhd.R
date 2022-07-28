library(tidyverse)
library(nhdplusTools)
library(sf)
library(mapview)
mv <- mapview::mapview

options(timeout = 5000)

#setup ####

#establish this path to local nhd highres data
nhd_hr_dir <- '~/git/macrosheds/data_acquisition/data/general/nhd_hr'
#and this one if you want to save interactive map outputs (also activate commented chunk in loop)
# mapview_save_dir <- '~/git/macrosheds/data_acquisition/output/sites_vs_NHD'

buf <- function(site, buf_dist){

    site_buf <- sf::st_buffer(x = site,
                              dist = buf_dist)
    site_box <- st_bbox(site_buf)

    subset_file <- tempfile(fileext = '.gpkg')

    subset <- try({
        nhdplusTools::subset_nhdplus(bbox = site_box,
                                     output_file = subset_file,
                                     nhdplus_data = 'download',
                                     return_data = TRUE,
                                     overwrite = FALSE,
                                     out_prj = 4326) %>%
            suppressMessages()
    }, silent = TRUE)

    if(inherits(subset, 'try-error') || ! length(subset)){# || nrow(subset[[1]]) < 2){
        print('incrementing buffer distance by 500')
        buf_dist <- buf_dist + 500
        if(buf_dist > 5000) stop()
        buf(site = site, buf_dist)
    } else {
        return(list(subset = subset, box = site_box))
    }
}

# site_csv <- read_csv('~/git/macrosheds/data_acquisition/data/general/site_data.csv')
# sites = tribble(
#     ~site_code, ~domain, ~latitude, ~longitude, ~epsg,
#     'a', 'b', 41.98732, -72.60537, 4269,
#     'a', 'b', 38.59614, -77.05603, 4269,
# )

sites = read_csv('~/Downloads/duplicate_COMIDS.csv') %>%
    select(-site_name) %>%
    rename(latitude = lat, longitude = lon, site_code = nhdplus_id, site_name = long_name) %>%
    mutate(epsg = 4269, NHD_COMID = NULL)
# sites=sites[4:nrow(sites),]

sites$NHD_COMID <- '?'
total_len <- nrow(sites)

# loop 1: NHDPlusV2 or NHD-HR ####

#this loop is for identifying whether a point is on the NHDPlusV2 or the NHD-HR,
#   or neither. for NHM seg_ids, see the next loop
prev_huc4 <- 'none'
for(i in 1:total_len){

    print('---')
    print(i)
    site <- sites[i, ]
    site_code <- site$site_code
    epsg <- site$epsg

    print(paste('site:', site_code))
    print(paste('name:', sites$site_name[i]))

    site <- sf::st_point(c(site$longitude, site$latitude)) %>%
                         sf::st_sfc(crs = site$epsg)

    nextt <- FALSE
    comid <- tryCatch({
        nhdplusTools::discover_nhdplus_id(site)
    }, error = function(e) nextt <<- TRUE)
    if(nextt) {
        print('couldnt find comid')
        next
    }

    out <- suppressWarnings(buf(site = site, 1000))
    subset <- out$subset
    site_box <- out$box

    huc12 <- get_huc12(site)
    # print(huc12$huc12)
    huc4 <- substr(huc12$huc12, 1, 4)[1]
    nhdplusTools::download_nhdplushr(nhd_hr_dir,
                                     hu_list = huc4) %>%
        invisible()

    if(huc4 != prev_huc4){
        HRflowlines <- nhdplusTools::get_nhdplushr(file.path(nhd_hr_dir,
                                                             substr(huc4, 1, 2)),
                                                   file.path(nhd_hr_dir,
                                                             paste0(huc4, '.gpkg')),
                                                   layers = 'NHDFlowline',
                                                   proj = epsg)$NHDFlowline
    } else {
        print('using previous NHD HR HUC')
    }

    prev_huc4 <- huc4

    NHD_HR <- suppressWarnings(sf::st_crop(HRflowlines, site_box))

    NHDPlus <- subset$NHDFlowline_Network %>%
        st_transform(crs = epsg)

    dist_to_nearest_NHDPlus_flowline <- min(st_distance(NHDPlus, site))
    print(paste('Dist to NHDPlus:', round(dist_to_nearest_NHDPlus_flowline, 2), 'm'))
    dist_to_nearest_NHDHR_flowline <- min(st_distance(NHD_HR, site))
    print(paste('Dist to NHD_HR:', round(dist_to_nearest_NHDHR_flowline, 2), 'm'))

    options(scipen = 100)
    xx = mv(NHD_HR, color = 'darkslategray3') + mv(NHDPlus, color='deepskyblue4') + mv(site, color='red')
    print(xx)

    # mapview_save_path <- file.path(mapview_save_dir,
    #                                paste0(dmn, '_', site_code, '.html'))
    # mapview::mapshot(xx,
    #                  url = mapview_save_path)
    # print(paste('map saved to', mapview_save_path))
    # print(xx)

    system('spd-say "chili chili chili"')
    x <- readline(cat('This point is on: [A] an NHDPlus flowline, [B] an NHD_HR flowline, or [C] neither >\n'))

    if(x == 'A'){
        sites[i, 'NHD_COMID'] <- as.character(comid)
        print(comid)
    } else if(x == 'B'){
        sites[i, 'NHD_COMID'] <- 'HR only'
    } else if(x == 'C'){
        sites[i, 'NHD_COMID'] <- 'too small'
    } else {
        stop(paste("'A', 'B', or 'C'"))
    }
}
