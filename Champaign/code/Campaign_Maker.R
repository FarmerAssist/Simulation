# This script creates campaign for calibration.
# Spatial boundaries of the campaign are determined by county limits
# Same crop is grown at each pixel for all the years

library(pSIMSSiteMaker)
library(pSIMCampaignManager)
setwd("/mnt/iccp_storage/Regional_Calibration/DNN_APSIM")

fname <-  file.path(getwd(), "Campaign.nc4")



Create_Empty_Campaign(lat=seq(37, 43.5, by=0.25),
                      lon=seq(-86, -92.75, by=-0.25),
                      num_scen=1,
                      filename =fname)

params <- list (
  pars = c("tt_emerg_to_endjuv","rue", "potKernelWt",
           "tt_flower_to_maturity", "row_spacing",
           "sowing_density","depth","head_grain_no_max",
           "water_fraction_full", "subsequent_amount","initial_amount"),
  Upper.limit = c(500,2.1,400,
                  1100,800,10,50,
                  1000,0.95,400,400),
  lower.limit = c(200,1.5,200,
                  600,500,6,10,
                  700,0.05,5,5),
  unit=c('C/day','C/day','C/day','C/day','C/day','C/day','C/day','C/day','C/day',
         'C/day','C/day')
)

#----------------------------------------- LHC
n.knot <- 100 ## Total ensembles

grid.lhc <- pmap_dfc(params[c(1,2,3)], function(pars, Upper.limit, lower.limit, probs){
  
  probs <- gen_latin_nD(t(matrix(0:1, ncol = 1, nrow = 2)), n.knot)
  
  eval(parse(text = paste0("qunif", "(p,", lower.limit, ",", Upper.limit, ")")), list(p = probs)) %>%
    as.data.frame()
}) %>%
  `colnames<-`(params$pars)

grid <- grid.lhc

#------------------------------------------

Add_Scenario(fname, nrow(grid)-1)
prop <- Inspect_Camp(fname)
num_scen <- Get_Camp_dim(fname)$Scen
count <- length(prop$Lat)*length(prop$Lon)
Inspect_Camp(fname)

for(param in params$pars) {
  print(param)
  new.values <-  seq_along(prop$Scen) %>%
    purrr::map(~matrix(  grid[[param]][.x], nrow = length(prop$Lat), ncol = length(prop$Lon)))
  
  #debugonce(AddVar_Campaign)
  AddVar_Campaign(fname,
                  Variable = list(Name=param,
                                  Unit=params$unit[which(params$pars==param)],
                                  missingValue=-99,
                                  value= new.values,
                                  longname="",
                                  prec="float"
                  ),
                  attr = list('long_name',"")
  )
}

#----------------Adding met
new.values <-seq_along(prop$Scen) %>%
  purrr::map(~matrix(sample(c(1:9), prop$Count,TRUE), nrow = length(prop$Lat), ncol = length(prop$Lon)))

AddVar_Campaign(fname,
                Variable = list(Name='file',
                                Unit='Mapping',
                                missingValue=-99,
                                value= new.values,
                                longname="",
                                prec="float"
                ),
                attr = list('long_name',"met00000.met,met00001.met,met00002.met,met00003.met,
                            met00004.met,met00005.met,met00006.met,met00007.met,met00008.met,
                            met00009.met"))

Edit_mapping_var (fname, 'file' , 'long_name', "met00000.met,met00001.met,met00002.met,met00003.met,met00004.met,met00005.met,met00006.met,met00007.met,met00008.met,met00009.met")



#-------------------------------- pdate
new.values <-seq_along(prop$Scen) %>%
  purrr::map(~matrix(runif(1, 90, 150) %>% floor(), nrow = length(prop$Lat), ncol = length(prop$Lon)))


AddVar_Campaign(fname,Variable = list(Name=paste0('pdate'),
                                      Unit='Julian day',
                                      missingValue=-99,
                                      prec='float',
                                      longname="",
                                      value= new.values))

######################################################################################################

remove_var_campaign(fname, varnames=c('myvar'))



