#documentation and helper functions are included in the function definition below.
#if you're in Rstudio, collapse the function with Alt+o. there's an example included at the bottom.
#you'll need these packages: dplyr, glue, sf, elevatr, raster, whitebox, mapview

delineate_watershed_from_point <- function(lat,
                                           long,
                                           crs,
                                           machine_status = 'n00b',
                                           write_dir,
                                           write_name,
                                           verbose = TRUE){

    #lat: numeric representing latitude of the pour point in decimal degrees
    #   (negative indicates southern hemisphere)
    #long: numeric representing longitude of the pour point in decimal degrees
    #   (negative indicates west of prime meridian)
    #crs: numeric representing the coordinate reference system (e.g. 4326 for WSG84)
    #machine_status: either '1337', indicating that your machine has >= 16 GB
    #   RAM, or 'n00b', indicating < 16 GB RAM. DEM resolution is chosen accordingly
    #write_dir: character. the directory in which to write output shapefile
    #write_name: character. the basename of the shapefile components to be written. i.e.
    #   <write_name>.shp, <write_name>.shx, <write_name>.prj, <write_name>.dbf
    #verbose: logical. determines the amount of informative messaging during run

    #details: Output will have CRS 4326 (WGS 84), though processing is done
    #   on projected data. the projection specifications are determined
    #   automatically, based on pour point location. Note that for mega-huge
    #   watersheds, this could result in an inaccurate watershed area calculation.

    #returns: a list containing the following components:
    #   watershed_area_ha: the area of the delineated watershed in hectares
    #       (meters squared divided by 10,000)
    #   buffer_radius_m: the width (meters) around the site location that was used when
    #       requesting a DEM (digital elevation model)
    #   snap_distance_m: the search radius (meters) around the pour point that was used
    #       to choose a stream to snap the pour point to.
    #   snap_method: either "standard", which snaps the pour point to the cell
    #       within snap_distance_m with highest flow accumulation, or "jenson",
    #       which snaps to the nearest flow line
    #   dem_resolution: passed to elevatr::get_elev_raster (z parameter).
    #       depends on supplied machine_status

    #NOTE: in addition to the packages loaded below, you'll need mapview to
    #   visualize delineated watershed boundaries

    require(dplyr)
    require(glue)
    require(sf)
    require(elevatr)
    require(raster)
    require(whitebox)

    sm <- suppressMessages
    sw <- suppressWarnings

    #the functions below are all helpers for subroutines of
    #   delineate_watershed_from_point. they should all stand alone pretty
    #   well, and some are generally useful. they are defined locally here,
    #   just for convenient distribution.

    #moving shapefiles can be annoying, since they're actually represented by
    #   3-4 files
    move_shapefiles <- function(shp_files, from_dir, to_dir,
                                new_name_vec = NULL){

        #shp_files is a character vector of filenames with .shp extension
        #   (.shx, .prj, .dbf are handled internally and don't need to be listed)
        #from_dir and to_dir are strings representing the source and destination
        #   directories, respectively
        #new_name_vec is an optional character vector of new names for each shape file.
        #   these can end in ".shp", but don't need to

        if(any(! grepl('\\.shp$', shp_files))){
            stop('All components of shp_files must end in ".shp"')
        }

        if(length(shp_files) != length(new_name_vec)){
            stop('new_name_vec must have the same length as shp_files')
        }

        dir.create(to_dir,
                   showWarnings = FALSE,
                   recursive = TRUE)

        for(i in 1:length(shp_files)){

            shapefile_base <- strsplit(shp_files[i], '\\.shp')[[1]]

            files_to_move <- list.files(path = from_dir,
                                        pattern = shapefile_base)

            extensions <- stringr::str_match(files_to_move,
                                    paste0(shapefile_base, '(\\.[a-z]{3})'))[, 2]

            if(is.null(new_name_vec)){
                new_name_base <- rep(shapefile_base, length(files_to_move))
            } else {
                new_name_base <- strsplit(new_name_vec[i], '\\.shp$')[[1]]
                new_name_base <- rep(new_name_base, length(files_to_move))
            }

            mapply(function(x, nm, ext) file.rename(from = paste(from_dir,
                                                                 x,
                                                                 sep = '/'),
                                                    to = glue('{td}/{n}{ex}',
                                                              td = to_dir,
                                                              n = nm,
                                                              ex = ext)),
                   x = files_to_move,
                   nm = new_name_base,
                   ext = extensions)
        }

        return()
    }

    #prompt users for stuff, provide single-character responses,
    #   reprompt if they don't choose one of the expected responses
    get_response_1char <- function(msg, possible_chars,
                                   subsequent_prompt = FALSE){

        #msg: character. a message that will be used to prompt the user
        #possible_chars: character vector of acceptable single-character responses

        if(subsequent_prompt){
            cat(paste('Please choose one of:',
                      paste(possible_chars,
                            collapse = ', '),
                      '\n> '))
        } else {
            cat(msg)
        }

        ch <- as.character(readLines(con = stdin(), 1))

        if(length(ch) == 1 && ch %in% possible_chars){
            return(ch)
        } else {
            get_response_1char(msg, possible_chars, subsequent_prompt = TRUE)
        }
    }

    #chooose an appropriate projection, based on location
    choose_projection <- function(lat = NULL, long = NULL,
                                  unprojected = FALSE){

        if(unprojected){
            PROJ4 <- glue('+proj=longlat +datum=WGS84 +no_defs ',
                          '+ellps=WGS84 +towgs84=0,0,0')
            return(PROJ4)
        }

        if(is.null(lat) || is.null(long)){
            stop('If projecting, lat and long are required.')
        }

        if(lat <= 15 && lat >= -15){ #equatorial
            PROJ4 = glue('+proj=laea +lon_0=', long)
        } else { #temperate or polar
            PROJ4 = glue('+proj=laea +lat_0=', lat, ' +lon_0=', long)
        }

        return(PROJ4)
    }

    #for determining whether the DEM extent wasn't big enough to allow full
    #   delineation
    raster_intersection_summary <- function(wb, dem){

        #wb is a delineated watershed boundary as a rasterLayer
        #dem is a DEM rasterLayer

        summary_out <- list()

        #convert wb to sf object
        wb <- wb %>%
            raster::rasterToPolygons() %>%
            sf::st_as_sf()

        #get edge of DEM as sf object
        dem_edge <- raster::boundaries(dem) %>%
            raster::reclassify(matrix(c(0, NA),
                                      ncol = 2)) %>%
            raster::rasterToPolygons() %>%
            sf::st_as_sf()

        #tally raster cells
        summary_out$n_wb_cells <- length(wb$geometry)
        summary_out$n_dem_cells <- length(dem_edge$geometry)

        #tally intersections; calc percent of wb cells that overlap
        intersections <- sf::st_intersects(wb, dem_edge) %>%
            as.matrix() %>%
            apply(MARGIN = 2,
                  FUN = sum) %>%
            table()

        true_intersections <- sum(intersections[names(intersections) > 0])

        summary_out$n_intersections <- true_intersections
        summary_out$pct_wb_cells_intersect <- true_intersections /
            summary_out$n_wb_cells * 100

        return(summary_out)
    }

    #the workhorse
    delineate_watershed_apriori <- function(lat, long, crs,
                                            machine_status = 'n00b',
                                            verbose = FALSE){

        #lat: numeric representing latitude in decimal degrees
        #   (negative indicates southern hemisphere)
        #long: numeric representing longitude in decimal degrees
        #   (negative indicates west of prime meridian)
        #crs: numeric representing the coordinate reference system (e.g. WSG84)
        #machine_status: either '1337', indicating that your machine has >= 16 GB
        #   RAM, or 'n00b', indicating < 16 GB RAM. DEM resolution is chosen accordingly
        #verbose: logical. determines the amount of informative messaging during run

        #returns the location of candidate watershed boundary files

        tmp <- tempdir()
        inspection_dir <- glue(tmp, '/INSPECT_THESE')
        dem_f <- glue(tmp, '/dem.tif')
        point_f <- glue(tmp, '/point.shp')
        d8_f <- glue(tmp, '/d8_pntr.tif')
        flow_f <- glue(tmp, '/flow.tif')

        dir.create(path = inspection_dir,
                   showWarnings = FALSE)

        proj <- choose_projection(lat = lat,
                                  long = long)

        site <- tibble(x = lat,
                       y = long) %>%
            sf::st_as_sf(coords = c("y", "x"),
                         crs = crs) %>%
            sf::st_transform(proj)
        # sf::st_transform(4326) #WGS 84 (would be nice to do this unprojected)

        #prepare for delineation loops
        buffer_radius <- 100
        dem_coverage_insufficient <- FALSE
        while_loop_begin <- TRUE

        #snap site to flowlines 3 different ways. delineate watershed boundaries (wb)
        #for each unique snap. if the delineations get cut off, get more elevation data
        #and try again
        while(while_loop_begin || dem_coverage_insufficient){

            while_loop_begin <- FALSE

            if(machine_status == '1337'){
                dem_resolution <- case_when(
                    buffer_radius <= 1e4 ~ 12,
                    buffer_radius == 1e5 ~ 11,
                    buffer_radius == 1e6 ~ 10,
                    buffer_radius == 1e7 ~ 8,
                    buffer_radius == 1e8 ~ 6,
                    buffer_radius == 1e9 ~ 4,
                    buffer_radius >= 1e10 ~ 2)
            } else if(machine_status == 'n00b'){
                dem_resolution <- case_when(
                    buffer_radius <= 1e4 ~ 10,
                    buffer_radius == 1e5 ~ 8,
                    buffer_radius == 1e6 ~ 6,
                    buffer_radius == 1e7 ~ 4,
                    buffer_radius == 1e8 ~ 2,
                    buffer_radius >= 1e9 ~ 1)
            } else {
                stop('machine_status must be either "1337" or "n00b"')
            }

            site_buf <- sf::st_buffer(x = site,
                                      dist = buffer_radius)
            dem <- elevatr::get_elev_raster(locations = site_buf,
                                            z = dem_resolution,
                                            verbose = verbose)

            raster::writeRaster(x = dem,
                                filename = dem_f,
                                overwrite = TRUE)

            sf::st_write(obj = site,
                         dsn = point_f,
                         delete_layer = TRUE,
                         quiet = TRUE)

            whitebox::wbt_fill_single_cell_pits(dem = dem_f,
                                                output = dem_f)

            whitebox::wbt_breach_depressions(dem = dem_f,
                                             output = dem_f,
                                             flat_increment = 0.01)

            whitebox::wbt_d8_pointer(dem = dem_f,
                                     output = d8_f)

            whitebox::wbt_d8_flow_accumulation(input = dem_f,
                                               output = flow_f,
                                               out_type = 'catchment area')

            snap1_f <- glue(tmp, '/snap1_jenson_dist150.shp')
            whitebox::wbt_jenson_snap_pour_points(pour_pts = point_f,
                                                  streams = flow_f,
                                                  output = snap1_f,
                                                  snap_dist = 150)
            snap2_f <- glue(tmp, '/snap2_standard_dist50.shp')
            whitebox::wbt_snap_pour_points(pour_pts = point_f,
                                           flow_accum = flow_f,
                                           output = snap2_f,
                                           snap_dist = 50)
            snap3_f <- glue(tmp, '/snap3_standard_dist150.shp')
            whitebox::wbt_snap_pour_points(pour_pts = point_f,
                                           flow_accum = flow_f,
                                           output = snap3_f,
                                           snap_dist = 150)

            #the site has been snapped 3 different ways. identify unique snap locations.
            snap1 <- sf::st_read(snap1_f, quiet = TRUE)
            snap2 <- sf::st_read(snap2_f, quiet = TRUE)
            snap3 <- sf::st_read(snap3_f, quiet = TRUE)
            unique_snaps_f <- snap1_f
            if(! identical(snap1, snap2)) unique_snaps_f <- c(unique_snaps_f, snap2_f)
            if(! identical(snap1, snap3)) unique_snaps_f <- c(unique_snaps_f, snap3_f)

            #good for experimenting with snap specs:
            # delineate_watershed_test2(tmp, point_f, flow_f,
            #                           d8_f, 'standard', 1000)

            #delineate each unique location
            for(i in 1:length(unique_snaps_f)){

                rgx <- stringr::str_match(unique_snaps_f[i],
                                 '.*?_(standard|jenson)_dist([0-9]+)\\.shp$')
                snap_method <- rgx[, 2]
                snap_distance <- rgx[, 3]

                wb_f <- glue('{path}/wb{n}_buffer{b}_{typ}_dist{dst}.tif',
                             path = tmp,
                             n = i,
                             b = buffer_radius,
                             typ = snap_method,
                             dst = snap_distance)

                whitebox::wbt_watershed(d8_pntr = d8_f,
                                        pour_pts = unique_snaps_f[i],
                                        output = wb_f)

                wb <- raster::raster(wb_f)

                #check how many wb cells coincide with the edge of the DEM.
                #If > 0.1% or > 5, broader DEM needed
                smry <- raster_intersection_summary(wb = wb,
                                                    dem = dem)

                if(verbose){
                    print(glue('buffer radius: {br}; snap: {sn}/{tot}; ',
                               'n intersecting cells: {ni}; pct intersect: {pct}',
                               br = buffer_radius,
                               sn = i,
                               tot = length(unique_snaps_f),
                               ni = round(smry$n_intersections, 2),
                               pct = round(smry$pct_wb_cells_intersect, 2)))
                }

                if(smry$pct_wb_cells_intersect > 0.1 || smry$n_intersections > 5){
                    buffer_radius_new <- buffer_radius * 10
                    dem_coverage_insufficient <- TRUE
                } else {
                    dem_coverage_insufficient <- FALSE
                    buffer_radius_new <- buffer_radius

                    #write and record temp files for the technician to visually inspect
                    wb_sf <- wb %>%
                        raster::rasterToPolygons() %>%
                        sf::st_as_sf() %>%
                        sf::st_buffer(dist = 0.1) %>%
                        sf::st_union() %>%
                        sf::st_as_sf()#again? ugh.

                    wb_sf <- sf::st_transform(wb_sf, 4326) #EPSG for WGS84

                    wb_sf_f <- glue('{path}/wb{n}_BUF{b}{typ}DIST{dst}RES{res}.shp',
                                    path = inspection_dir,
                                    n = i,
                                    b = buffer_radius,
                                    typ = snap_method,
                                    dst = snap_distance,
                                    res = dem_resolution)

                    sf::st_write(obj = wb_sf,
                                 dsn = wb_sf_f,
                                 delete_dsn = TRUE,
                                 quiet = TRUE)
                }
            }

            buffer_radius <- buffer_radius_new
        } #end while loop

        if(verbose){
            message(glue('Candidate delineations are in: ', inspection_dir))
        }

        return(inspection_dir)
    }

    if(verbose){
        message('Beginning watershed delineation')
    }

    inspection_dir <- sw(delineate_watershed_apriori(
        lat = lat,
        long = long,
        crs = crs,
        machine_status = machine_status,
        verbose = verbose))

    files_to_inspect <- list.files(path = inspection_dir,
                                   pattern = '.shp')

    #if only one delineation, write it to write_dir
    if(length(files_to_inspect) == 1){

        selection <- files_to_inspect[1]

        move_shapefiles(shp_files = selection,
                        from_dir = inspection_dir,
                        to_dir = write_dir)

        if(verbose){
            message(glue('Delineation successful. Shapefile written to ',
                         write_dir))
        }

    #otherwise, inspect all delineations and choose one
    } else {

        nshapes <- length(files_to_inspect)

        wb_selections <- paste(paste0('[',
                                      c(1:nshapes, 'A'),
                                      ']'),
                               c(files_to_inspect, 'Abort delineation'),
                               sep = ': ',
                               collapse = '\n')

        helper_code <- glue('mapview::mapview(sf::st_read("{wd}/{f}"))',
                            wd = inspection_dir,
                            f = files_to_inspect) %>%
            paste(collapse = '\n\n')

        msg <- glue('Visually inspect the watershed boundary candidate shapefiles ',
                    'in {td}, then enter the number corresponding to the ',
                    'one that looks most legit. Here\'s some ',
                    'helper code you can paste into a separate R instance ',
                    ':\n\n{hc}\n\nChoices:\n{sel}\n\nEnter choice here > ',
                    hc = helper_code,
                    sel = wb_selections,
                    td = inspection_dir)

        resp <- get_response_1char(msg = msg,
                                   possible_chars = c(1:nshapes, 'A'))

        if(resp == 'A'){
            message('Delineation aborted')
            return()
        }

        selection <- files_to_inspect[as.numeric(resp)]

        move_shapefiles(shp_files = selection,
                        from_dir = inspection_dir,
                        to_dir = write_dir,
                        new_name_vec = write_name)

        message(glue('Selection {s}:\n\t{sel}\nwas written to:\n\t{sdr}',
                     s = resp,
                     sel = selection,
                     sdr = write_dir))
    }

    #calculate watershed area in hectares
    wb <- sf::st_read(glue('{d}/{s}.shp',
                           d = write_dir,
                           s = write_name),
                      quiet = TRUE)

    ws_area_ha <- as.numeric(sf::st_area(wb)) / 10000

    #return the specifications of the correctly delineated watershed, and some
    #   other goodies
    rgx <- stringr::str_match(selection,
                     paste0('^wb[0-9]+_BUF([0-9]+)(standard|jenson)',
                            'DIST([0-9]+)RES([0-9]+)\\.shp$'))

    deets <- list(name = write_name,
                  watershed_area_ha = ws_area_ha,
                  buffer_radius_m = as.numeric(rgx[, 2]),
                  snap_distance_m = as.numeric(rgx[, 4]),
                  snap_method = rgx[, 3],
                  dem_resolution = as.numeric(rgx[, 5]))

    return(deets)
}

deets <- delineate_watershed_from_point(lat = 44.21013,
                                        long = -122.2571,
                                        crs = 4326,
                                        write_dir = '~/some_path',
                                        write_name = 'ultimate_watershed')
