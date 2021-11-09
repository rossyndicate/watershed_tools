library(raster)
library(tidyverse)
library(sf)
library(rgee)

#setup ####

# setup for rgee is super annoying, and must be done anew on each machine
rgee::ee_Initialize(email = 'vlahm13@gmail.com',
                    drive = TRUE)

abp_wb = sf::st_read('~/Downloads/Shape Files/Al-Beuhler_Pond.shp')
cp_wb = sf::st_read('~/Downloads/Shape Files/Circuit_Pond.shp')
gpp_wb = sf::st_read('~/Downloads/Shape Files/GardenPond_Pond.shp')

get_nlcd_summary = function(watershed_boundary,
                            epochs = c(1992, 2001, 2004, 2006, 2008, 2011, 2013, 2016)){

    #summarizes landcover classes from NLCD for a spatial extent. Returns
    #   a two-element list. element 1 is a tibble containing cell tallies for
    #   each class and year (epoch).
    #   Element 2 is a named vector of impervious surface percentages for each
    #   epoch.

    #watershed_boundary is an sf object representing the outline of a watershed.
    #   NLCD categories will be summarized for that areal extent.
    #epochs is a numeric vector of years for which NLCD data will be collected.
    #   available epochs are 1992, 2001, 2004, 2006, 2008, 2011, 2013, 2016

    if(! all(epochs %in% c(1992, 2001, 2004, 2006, 2008, 2011, 2013, 2016))){
        stop('epochs must be one or more of 1992, 2001, 2004, 2006, 2008, 2011, 2013, 2016')
    }

    color_key = tribble(
        ~id, ~hexcode, ~category,
        11,"466b9f","Open water",
        12,"d1def8","Perennial ice/snow",
        21,"dec5c5","Developed, open space",
        22,"d99282","Developed, low intensity",
        23,"eb0000","Developed, medium intensity",
        24,"ab0000","Developed high intensity",
        31,"b3ac9f","Barren land (rock/sand/clay)",
        41,"68ab5f","Deciduous forest",
        42,"1c5f2c","Evergreen forest",
        43,"b5c58f","Mixed forest",
        51,"af963c","Dwarf scrub",
        52,"ccb879","Shrub/scrub",
        71,"dfdfc2","Grassland/herbaceous",
        72,"d1d182","Sedge/herbaceous",
        73,"a3cc51","Lichens",
        74,"82ba9e","Moss",
        81,"dcd939","Pasture/hay",
        82,"ab6c28","Cultivated crops",
        90,"b8d9eb","Woody wetlands",
        95,"6c9fb8","Emergent herbaceous wetlands",
    )

    watershed_boundary_ee = sf_as_ee(watershed_boundary)

    tmp = tempdir()
    tmp = paste0(tmp, '/nlcd')
    dir.create(tmp, showWarnings = FALSE)

    pct_impervious_surface = rep(NA, length(epochs))
    names(pct_impervious_surface) = epochs
    for(i in seq_along(epochs)){

        subset_id = paste0('NLCD', as.character(epochs[i]))

        #write landcover rasters

        img = ee$ImageCollection('USGS/NLCD')$
            select('landcover')$
            filter(ee$Filter$eq('system:index', subset_id))$
            first()$
            clip(watershed_boundary_ee)

        ee_as_raster(image = img,
                     region = watershed_boundary_ee$geometry(),
                     dsn = paste0(tmp, '/', subset_id, '.tif'))

        #get impervious surface %

        img = ee$ImageCollection('USGS/NLCD')$
            select('impervious')$
            filter(ee$Filter$eq('system:index', subset_id))$
            first()

        imp_surf = tryCatch({
            ee_extract(x = img,
                       y = watershed_boundary,
                       scale = 30,
                       fun = ee$Reducer$mean(),
                       sf = TRUE) %>%
                pull(impervious)
        }, error = function(e) NA)

        pct_impervious_surface[i] = imp_surf
    }

    nlcd_summary = color_key %>%
        select(id, category) %>%
        mutate(id = as.character(id))

    for(e in epochs){

        subset_id = paste0('NLCD', as.character(e))

        epoch_rst = raster::raster(paste0(tmp, '/', subset_id, '.tif'))
        tabulated_values = raster::values(epoch_rst) %>%
            table() %>%
            as_tibble() %>%
            dplyr::rename(id = '.',
                          !!sym(paste0('CellTally', as.character(e))) := 'n')

        nlcd_summary = full_join(nlcd_summary,
                                 tabulated_values,
                                 by = 'id')
    }

    return(list(landcover = nlcd_summary,
                pct_impervious_surface = pct_impervious_surface))
}

# summarize landcover and % impervious surfaces ####

abp_summary = get_nlcd_summary(watershed_boundary = abp_wb,
                               epochs = 2016)

cp_summary = get_nlcd_summary(watershed_boundary = cp_wb,
                              epochs = 2016)

gpp_summary = get_nlcd_summary(watershed_boundary = gpp_wb,
                               epochs = 2016)

saveRDS(abp_summary, '~/Downloads/al_beuhler_pond.rds')
saveRDS(cp_summary, '~/Downloads/circuit_pond.rds')
saveRDS(gpp_summary, '~/Downloads/gardenpond_pond.rds')
