library(macrosheds)
library(tidyverse)
library(mapview)

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

combined_map <- mapview(wbs[[1]], col.regions = 1, layer.name = 1)
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


ff = read_csv('in/field_Q/erwin.csv',
              col_types = cols(.default = 'character'))
ff %>%
    mutate(time = str_pad(`T`, pad = '0', width = 5),
           date = gsub('\\.', '/', D),
           datetime = paste(D, time),) %>%
    mutate(datetime = mdy_hm(datetime, tz = 'US/Eastern'),
           datetime = with_tz(datetime, 'UTC'),
           discharge = as.numeric(discharge)) %>%
    select(datetime, discharge) %>%
    write_csv('in/field_Q/erwin.csv')

#plot modeled Q against rating Q, data, etc

site_ = 'erwin'
qq = read_csv('raw/nhc_discharge_detrended.csv')
qq$NHC_discharge_Ls = qq$NHC_discharge_m3s * 1000
plot(qq$DateTime_UTC, qq$NHC_discharge_Ls, ylim = c(0, 5000), col = 'gray50', type ='l',
     ylab = 'Discharge (L/s)', xlab = '2023-24')
legend('top', legend = c('rating-derived', 'measurements', 'predicted'),
       fill = c('gray50', 'blue', 'red'))

qp = read_csv(glue('out/lm_out/predictions/{site_}.csv'))
# qp = read_csv(glue('out/seasons_and_slugs/predictions/erwin.csv'))
# qp = read_csv(glue('out/no_season_no_slugs/predictions/erwin.csv'))
points(qp$datetime, qp$Q_predicted, col = 'red', pch = '.')
qm = read_csv(glue('in/field_Q/{site_}.csv'))
points(qm$datetime, qm$discharge, col = 'blue', pch = 42)

plot(qp$datetime, qp$Q_predicted, col = 'red', pch = '.')

jj = qq %>% filter(DateTime_UTC > as.Date('2023-12-25'))
plot(jj$DateTime_UTC, jj$NHC_discharge_Ls, ylim = c(0, 5000))
jj[sample(1:nrow(jj), 100, replace=F), ] %>%
    arrange(DateTime_UTC) %>%
    select(datetime = DateTime_UTC, discharge = NHC_discharge_Ls) %>%
    write_csv('in/field_Q/erwin.csv')

plot(qp$datetime, qp$Q_predicted, ylim = c(0, 5000), col = 'red', type ='l',
     ylab = 'Discharge (L/s)', xlab = '')
points(qm$datetime, qm$discharge, col = 'blue', pch = 42)


#eval stuff
group_by(qp, season) %>% summarize(zz = median(Q_predicted))
left_join(qm, select(qp, datetime, Q_predicted)) %>%
    ggplot2::ggplot(aes(x = discharge, y = Q_predicted)) +
    geom_point()
