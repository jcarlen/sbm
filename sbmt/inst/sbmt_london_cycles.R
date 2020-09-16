# Run sbmt on london cycles data used in ppsbm paper. 
# Jane Carlen
# Created: 9-14-20
# Download data from: http://cmatias.perso.math.cnrs.fr/
#   linked  "http://cmatias.perso.math.cnrs.fr/Docs/ppsbm-files.zip (under "R code with datasets analyses")
#   with the entry for the  publication: "Catherine Matias, Tabea Rebafka & Fanny Villers, A semiparametric extension of the stochastic block model for longitudinal networks. Biometrika, 105(3): 665-680, 2018."
#
# NOTE: Set working directory to top of downloaded ppsbm-files folder
# -------------------------------------------------------------------------------------------------------------
# libraries ----
# devtools::install_github("jcarlen/sbm", subdir = "sbmt") 
library(sbmt) 
library(tidyverse)
library(cowplot)

# load data ----
load("data/London_cycles_dynsbm.Rdata") #for data.day1
locations <- read.table('data/London_cycles_locations_stations.txt',header=T,sep=";")

#   station info ----

station_ids <- c(data.day1$EndStation.Id,data.day1$StartStation.Id)  #see toto below
station_numbers <- setNames(station_ids, match(station_ids, unique(sort(station_ids))))

#   data to sbmt format ----

# hourly time slice
data.day1$hour = floor(data.day1$Start.Date/3600)
# any self-edges in data? no
sum(data.day1$StartStation.Number == data.day1$EndStation.Number)

# make edgelist
lbs_edgelist = lapply(0:max(data.day1$hour), function(i) {
  data.day1 %>% filter(hour == i) %>% 
    select(StartStation.Number, EndStation.Number) %>% 
    group_by(StartStation.Number, EndStation.Number) %>% 
    summarize(count = n()) %>% 
    ungroup() %>% 
    rename(from = StartStation.Number, to = EndStation.Number) %>% 
    as.matrix()
})

# run sbmt with and without degree correction ----

set.seed(1)
Sys.time()
lbs_output_dc0 = sbmt(lbs_edgelist, directed = TRUE, maxComms = 6, degreeCorrect = 0, klPerNetwork = 5)
Sys.time()
lbs_output_dc3 = sbmt(lbs_edgelist, directed = TRUE, maxComms = 6, degreeCorrect = 3, klPerNetwork = 5)
Sys.time()

# examine and compare output ----

#   plot found blocks by station id ----

found_comm_stations_dc0 = setNames(lbs_output_dc0$FoundComms, nm = as.vector(station_numbers[names(lbs_output_dc0$FoundComms)])); table(found_comm_stations_dc0)
found_comm_stations_dc3 = setNames(lbs_output_dc3$FoundComms, nm = as.vector(station_numbers[names(lbs_output_dc3$FoundComms)])); table(found_comm_stations_dc3)

par(mfrow = c(1,2))

# without correction
plot(locations[match(names(found_comm_stations_dc0), locations$Id), c("long","lat")], 
  col = found_comm_stations_dc0 + 1,
  pch = 16, main = "sbmt blocks no degree correction")

# with
plot(locations[match(names(found_comm_stations_dc3), locations$Id), c("long","lat")], 
  col = found_comm_stations_dc3 + 1,
  pch = 16, main = "sbmt blocks degree corrected")

#   plot estimated block-to-block parameters ----

# no degree correction
plot(lbs_output_dc0, mai = rep(.3, 4))
# degree correction
plot(lbs_output_dc3, mai = rep(.3, 4))

# compare highly active station block ----

# no degreee correction
station_output_info = cbind(locations[match(names(found_comm_stations_dc0), locations$Id),], block = found_comm_stations_dc0)
belgrove_block = station_output_info %>% filter(name == "Belgrove Street, Kings Cross") %>% pull(block)
station_output_info %>% filter(block == belgrove_block) #without degree correcton, a small group containing only Belgrove Street, Kings Cross; Waterloo Station 3, Waterloo; Waterloo Station 1, Waterloo 

plot_grid(nrow =2,
          data.day1 %>% filter(StartStation.Name %in% (station_output_info %>% filter(block == belgrove_block) %>% pull(name))) %>%
            group_by(StartStation.Name, hour) %>% summarize(count = n()) %>% 
            ggplot(aes(x = hour, y = count, group = StartStation.Name)) + geom_line() + ggtitle("out-traffic from block"),
          
          data.day1 %>% filter(EndStation.Name %in% (station_output_info %>% filter(block == belgrove_block) %>% pull(name))) %>% 
            group_by(EndStation.Name, hour) %>% summarize(count = n()) %>% 
            ggplot(aes(x = hour, y = count, group = EndStation.Name)) + geom_line() + ggtitle("in-traffic to block"))


# degreee correction

station_output_info = cbind(locations[match(names(found_comm_stations_dc3), locations$Id),], block = found_comm_stations_dc3)
belgrove_block = station_output_info %>% filter(name == "Belgrove Street, Kings Cross") %>% pull(block)
station_output_info %>% filter(block == belgrove_block) #without degree correcton, a small group containing only Belgrove Street, Kings Cross; Waterloo Station 3, Waterloo; Waterloo Station 1, Waterloo 

plot_grid(nrow =2,
          data.day1 %>% filter(StartStation.Name %in% (station_output_info %>% filter(block == belgrove_block) %>% pull(name))) %>%
            group_by(StartStation.Name, hour) %>% summarize(count = n()) %>% 
            ggplot(aes(x = hour, y = count, group = StartStation.Name)) +geom_line() + ggtitle("out-traffic from block"),
                   
          data.day1 %>% filter(EndStation.Name %in% (station_output_info %>% filter(block == belgrove_block) %>% pull(name))) %>% 
            group_by(EndStation.Name, hour) %>% summarize(count = n()) %>% 
            ggplot(aes(x = hour, y = count, group = EndStation.Name)) +geom_line() + ggtitle("in-traffic to block"))
