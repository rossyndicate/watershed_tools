dem_unprojected <- raster::projectRaster(dem, crs = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
dem_bounds <- sf::st_bbox(dem_unprojected)
# highway_types <- c('motorway', 'trunk', 'primary')
highway_types <- c('motorway', 'trunk', 'primary', 'secondary', 'tertiary')
highway_types <- c(highway_types,
                   paste(highway_types, 'link', sep = '_'))
roads_query <- opq(dem_bounds) %>%
    add_osm_feature(key = 'highway',
                    value = highway_types)
roads <- osmdata_sf(roads_query)
plot(roads$osm_lines)

streams_query <- opq(dem_bounds) %>%
    add_osm_feature(key = 'waterway',
                    value = c('river', 'stream'))
streams <- osmdata_sf(streams_query)
streams$osm_lines$geometry
plot(streams$osm_lines)

streams_f <- glue(tmp, '/streams.shp')
roads_f <- glue(tmp, '/roads.shp')
sf::st_transform(streams$osm_lines$geometry,
                 crs = proj) %>%
    sf::st_write(dsn = tmp,
                 layer = 'streams.shp',
                 driver = 'ESRI Shapefile',
                 delete_layer = TRUE,
                 silent = TRUE)
sf::st_transform(roads$osm_lines$geometry,
                 crs = proj) %>%
    sf::st_write(dsn = tmp,
                 layer = 'roads.shp',
                 driver = 'ESRI Shapefile',
                 delete_layer = TRUE,
                 silent = TRUE)

whitebox::wbt_burn_streams_at_roads(dem = dem_f,
                                    streams = streams_f,
                                    roads = roads_f,
                                    output = dem_f,
                                    width = 50)

