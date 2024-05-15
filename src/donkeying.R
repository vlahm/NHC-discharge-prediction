library(macrosheds)
library(tidyverse)

whitebox::wbt_init(exe_path = '~/git/others_projects/whitebox-tools/target/release/whitebox_tools')

setwd('~/git/macrosheds/papers/nhc_q')
# options(whitebox.exe_path = "~/docker_whitebox.sh")
# options(whitebox.exe_path = "/home/mike/git/others_projects/whitebox-tools")

ss = readr::read_csv('cfg/sites.csv') %>%
    filter(! grepl('Down', site_code))

i = 6
macrosheds::ms_delineate_watershed(ss$lat[i], ss$lon[i],
                                   write_dir = 'in/ws_bounds',
                                   write_name = ss$site_code[i])


zz = list.files('in/ws_bounds/', pattern = '\\.shp', full.names = T)

wbs = list()
ind = 0
for(z in zz){
    ind = ind + 1
    qfg = sf::st_read(z, quiet = TRUE)
    wbs[[ind]] = qfg
}

combined_map <- mapview(wbs[[1]], col.regions = i, layer.name = 1)
for (i in 2:length(wbs)) {
    combined_map <- combined_map + mapview(wbs[[i]], col.regions = i, layer.name = i)
}
combined_map

mv = mapview::mapview

mv(wbs[[1]], col.regions = 'blue', layer.name= wbs[[1]]$site_code)+
    mv(wbs[[2]], col.regions = 'pink', layer.name= wbs[[2]]$site_code) +
    mv(wbs[[3]], col.regions = 'orange', layer.name= wbs[[3]]$site_code) +
    mv(wbs[[4]], col.regions = 'yellow', layer.name= wbs[[4]]$site_code) +
    mv(wbs[[5]], col.regions = 'green', layer.name= wbs[[5]]$site_code) +
    mv(wbs[[6]], col.regions = 'red', layer.name= wbs[[6]]$site_code)
