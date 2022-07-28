# watershed_tools
master
development ground for general watershed tools (delineation, viz, remote sensing, etc...)


### `2_batch_summary_contUS.R`
 + Includes functions for:
   1. acquiring NHDPlusV2 COMIDs based on lat/long
   2. acquiring NHDPlusV2 VPU IDs based on lat/long
   3. using COMID and VPU to acquire NHDPlusV2 data for your sites
   4. using COMID to acquire StreamCat data for your sites
   5. determining the position of a site along its NHD reach and calculating reach proportion, so that you can adjust linear and areal summary data accordingly
 + These only work in the continental USA.
 + Broader data acquisition (MODIS, NEON, etc.) using watershed boundary shapefiles is coming soon.
 + A global equivalent of these tools is also hopefully coming soon.
 + These functions will eventually be packaged.

