z = read_csv('junk/mgas_salt-slug_summary.csv')
z = z %>%
    select(Date, Site, slug_in_time, discharge_m3s) %>%
    mutate(time = str_pad(slug_in_time, width = 8, pad = '0')) %>%
    mutate(datetime = mdy_hms(paste(Date, slug_in_time), tz = 'US/Eastern'),
           discharge = discharge_m3s * 1000) %>%
    select(-Date, -slug_in_time, -time, -discharge_m3s) %>%
    mutate(Site = case_when(Site == 'cb_downstream' ~ 'NHCcbriffle',
                            Site == 'wb_downstream' ~ 'NHCwbriffle',
                            Site == 'cb_upstream' ~ 'NHCcbpoolUp',
                            Site == 'wb_upstream' ~ 'NHCwbpoolUp',
                            TRUE ~ Site),
           datetime = with_tz(datetime, 'UTC')) %>%
    group_split(Site)

sitecodes = sapply(z, function(x) x$Site[1])
z = lapply(z, function(x) x = select(x, -Site))

mapply(function(a, b) write_csv(a, paste0('in/field_Q/', b, '.csv')),
       a = z, b = sitecodes)
