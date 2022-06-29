
#NOTE - this will ask you to load in the data file (IUCN_Data.rds)
source("IUCN_source_code_GH.R")

#############################################################################
##FIGURE 1
###################################################################################

g1 <- ggplot(datAll2) + geom_bar(aes(status, fill = IslandEndemic)) +
  theme_bw() + xlab("IUCN Category") + ylab("Number of species") +
   scale_fill_manual(values=cbbPalette[3:4]) +
  theme(legend.position = "top") + ggtitle("a)") +
  labs(fill='') 

###work out proportions of species in each category 
#work out proportions (relative to total species in specific dataset)
#remember totals have DD species removed already, but include extinct species
tt_end <- table(end$status) / nrow(end)
tt_noend <- table(noend$status) / nrow(noend)

ttd <- data.frame("Status" = c(names(tt_end), names(tt_noend)),
                  "Prop" = c(as.vector(tt_end), as.vector(tt_noend)),
                  "Type" = c(rep("Endemic", length(tt_end)),
                             rep("Non-endemic", length(tt_end))))

ttd$Status <- factor(ttd$Status, levels = c("EX", "CR", "EN", "VU", "NT", "LC"))

g2 <- ggplot(ttd) + geom_col(aes(Status, Prop, fill = Type), position = "dodge") +
  theme_bw() + xlab("IUCN Category") + ylab("Proportion") +
  labs(fill="Species type") + scale_fill_manual(values=cbbPalette[3:4])+
  theme(legend.position = "none") + ggtitle("b)")


gridExtra::grid.arrange(g1, g2, nrow = 1)

################################################################
##########Null Model: Figure S2############################################
##################################################################

#across five categories, calculate number of island endemics in each,
#then randomly sample (without replacement) the same number from the
#whole species set and calculate the proportion of island endemics in 
#each. Repeat 9999 times and compare with observed

n <- table(datAll$IslandEndemic)["Endemic"]

null_end <- table(end$status)
null_end <- null_end[-which(names(null_end) == "LC"| names(null_end) =="NT")]
null_end["Threatened"] <- sum(null_end[c("CR", "EN", "VU")])
null_end <- data.frame("Category" = names(null_end), 
                       "Value" = as.vector(null_end))

##RUN NULL MODEL FOR STATUS
#run null_mod 9999 times
null999 <- replicate(9999, null_mod(datAll, n), simplify = "matrix") %>%
  t() %>%
  as.data.frame()

#drop the LC and NT columns
null999$LC <- NULL
null999$NT <- NULL

#add all threatened sp colum
null999 <- mutate(null999, Threatened = c(CR + EN + VU))

#make long version of it for ggplot
null999_long <- tidyr::pivot_longer(null999, cols = EX:Threatened,
                    names_to = "Category")

null999_long$Category <- factor(null999_long$Category, 
                                levels = c("EX", "CR", "EN", "VU", "Threatened"))

##violin plots of each category's null distribution with observed as *
g3 <- ggplot(data = null999_long, aes(Category, value)) + 
  geom_violin(fill = cbbPalette[1]) +
  theme_bw() + xlab("IUCN Category") + ylab("Number of species") + 
  stat_summary(fun.data=mean_sdl, geom = "crossbar", width=0.05,
               fill=cbbPalette[2]) +
  geom_point(data = null_end, aes(Category, Value), size = 4,
             shape = 18, col = 666666)

g3



zP_res <- matrix(ncol = 2, nrow = ncol(null999))
colnames(zP_res) <- c("Z", "P")

for (i in 1:ncol(null999)){
  obs <- null_end$Value[i]
  zz <- zP(null999[,i], obs, type = "ES")
  zP_res[i,] <- c(round(zz[1], 2), round(zz[2], 4))
}
zP_res <- as_tibble(zP_res) %>% 
  mutate("Category" = c("EX", "CR", "EN", "VU", "TH"))

  
####################################################
##ORDERS: FIGURE 5
###################################################################

#change datALL to 'end' to run with just the island endemic pool
ord_pool <- datAll #either datAll or end

#Select Orders with more than 150 species
Ord_Names <- table(datAll$Order)[table(datAll$Order) > 150] %>%
  sort %>%
  names()

##get observed numbers
#number of threatened island endemics in each category
th_obs <- filter(datAll, IslandEndemic == "Endemic", 
                 status %in% c("EX", "CR", "EN", "VU"))$status %>% table()
th_obs <- th_obs[c("EX", "CR", "EN", "VU")]
th_obs <- c(th_obs, "TH" = sum(th_obs[c("CR", "EN", "VU")]))


pp1 <- vector("list", length = length(th_obs))
pp2 <- vector("list", length = length(th_obs))

for (i in 1:length(th_obs)){

  pp <- plot_null_order(ord_pool, n = th_obs[i], Ord_Names, type = names(th_obs)[i],
                        zp_type = "ES")
  pp1[[i]] <- pp[[1]]
  pp1[[i]]$Category <- names(th_obs)[i]
  pp2[[i]] <- pp[[2]]
  pp2[[i]]$Category <- names(th_obs)[i]
  pp2[[i]]$Significance <- pp[[3]]$sig
  pp2[[i]]$Z <- pp[[3]]$Z
  pp2[[i]]$P <- pp[[3]]$P

}

come_together <- do.call(rbind, pp1)
come_together$Category <- factor(come_together$Category, 
                                 levels = c("EX", "CR", "EN", "VU", "TH"))
right_now <- do.call(rbind, pp2)
right_now$Category <- factor(right_now$Category, 
                             levels = c("EX", "CR", "EN", "VU", "TH"))

##violin plots of each category's null distribution with observed as *
g4 <- ggplot(data = come_together, aes(Order, value)) + 
  theme_bw() + xlab("Order") + ylab("Number of species") + 
  stat_summary(fun.data=mean_sdl, geom = "crossbar", width=0.05,
               fill=cbbPalette[2])  +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
   geom_point(data =  right_now, aes(Order, Number, col = Significance), size = 3,
              shape = 18) +
  facet_wrap(~Category, scales = "free_y") +
   scale_color_manual(labels = c("Non-significant", "Significant"),
                      values=cbbPalette[1:2])


g4

##############################################################
#############Make Trait Figure (Fig 6 / S4)#######################
########################################################

#set the species: 'All' extant species, or 'Isl' endemic extant species
Pool = "All"
#set the point size
PS <- 5
#set the point colours (sig and NS)
PCs <- c(cbbPalette[3], 666666)

#table for null model results
null_tab <- matrix(ncol = 3, nrow = 20)
colnames(null_tab) <- c("Variable", "ES", "P")
null_tab[,1] <- c("Volancy_Extant", "Volancy_Extinct", "HWI",
                  "Mass_Extant", "Mass_Extinct", "Forest", "Grassland",
                  "Marine", "Wetland", "Carnivore", "Herbivore", "Omnivore",
                  "Aq_Predator", "Frugivore", "Granivore", "Invertivore",
                  "Nectarivore", "Omnivore", "Vertivore", "Beak_volume")

###VOLANCY PLOT##############################

##Run null_volancy function
vol_run <- null_volancy(datExtant, datAll, trait = "volancy",
                        pool = Pool)

##Violin Null distribution plot
vr2 <- data.frame("Null_values" = vol_run[[2]][[1]], "Type" = "Extant")
vr3 <- data.frame("Null_values" = vol_run[[3]][[1]], "Type" = "Extinct")
vr4 <- rbind(vr2, vr3)

vr4$Type <- factor(vr4$Type, levels = c("Extant", "Extinct"))
vrOB <- data.frame("Obs" = c(vol_run[[2]][[2]], vol_run[[3]][[2]]),
                   "Type" = c("Extant", "Extinct"))

vrOB$Type <- factor(vrOB$Type, levels = c("Extant", "Extinct"))

#set point colour based on significance of ES
vol_col1 <- ifelse(vol_run[[2]][[3]][2] > 0.975 | vol_run[[2]][[3]][2] < 0.025,
                   PCs[1], PCs[2])
vol_col2 <- ifelse(vol_run[[3]][[3]][2] > 0.975 | vol_run[[3]][[3]][2] < 0.025,
                   PCs[1], PCs[2])
vol_col <- c(vol_col1, vol_col2)

g5a <- ggplot(data = vr4, aes(Type, Null_values)) + 
  theme_bw() + xlab("Species group") + ylab("No. of flightless species") + 
  geom_violin(fill = cbbPalette[1],  bw = 1.1) +
  stat_summary(fun.data=mean_sdl, geom = "crossbar", width=0.05,
               fill=cbbPalette[5]) +
  geom_point(data =  vrOB, aes(Type, Obs), size = PS,
             shape = 18, col = vol_col)  + ggtitle("a) Volancy")

#null model results
null_tab[1, 2:3] <- vol_run[[2]][[3]] #extant
null_tab[2, 2:3] <- vol_run[[3]][[3]] #extinct

###HWI##############################
TER = FALSE #normal analyses
datExtant_HWI <- datExtant #normal analyses

###remove marine species to check
# datExtant_HWI <- filter(datExtant, Habitat != "Marine")
# TER = TRUE
#########

hwi_run <- null_volancy(datExtant_HWI, datAll, zp_type = "ES", 
                        trait = "HWI",
                        pool = Pool, terrestrial = TER)
#make plot
vr4_hwi <- data.frame("Null_values" = hwi_run[[2]][[1]], "Type" = "Extant")
vrOB_hwi <- data.frame("Obs" = hwi_run[[2]][[2]],
                       "Type" = "Extant")
#set point colour based on significance of ES
hwi_col <- ifelse(hwi_run[[2]][[3]][2] > 0.975 | hwi_run[[2]][[3]][2] < 0.025,
                  PCs[1], PCs[2])

g5w <- ggplot(data = vr4_hwi, aes(Type, Null_values)) + 
  theme_bw() + xlab("") + ylab("Median HWI") + 
  geom_violin(fill = cbbPalette[1]) +
  stat_summary(fun.data=mean_sdl, geom = "crossbar", width=0.05,
               fill=cbbPalette[5]) +
  geom_point(data =  vrOB_hwi, aes(Type, Obs), size = PS,
             shape = 18, col = hwi_col) + ggtitle("b) HWI")

#null model results
null_tab[3, 2:3] <- hwi_run[[2]][[3]]

###BODY MASS##############################
mass_run <- null_volancy(datExtant, datAll, zp_type = "ES", 
                         trait = "Mass",
                         pool = Pool)
#make plot
vr_mass <- data.frame("Null_values" = c(mass_run[[2]][[1]],mass_run[[4]][[1]]),
                      "Type" = c("Extant", "Extinct"))

vr_mass$Type <- factor(vr_mass$Type, levels = c("Extant", "Extinct"))

vrOB_mass <- data.frame("Obs" =  c(mass_run[[2]][[2]],mass_run[[4]][[2]]),
                        "Type" = c("Extant", "Extinct"))

vrOB_mass$Type <- factor(vrOB_mass$Type, levels = c("Extant", "Extinct"))

#set point colour based on significance of ES
mass_col1 <- ifelse(mass_run[[2]][[3]][2] > 0.975 | mass_run[[2]][[3]][2] < 0.025,
                    PCs[1], PCs[2])
mass_col2 <- ifelse(mass_run[[4]][[3]][2] > 0.975 | mass_run[[4]][[3]][2] < 0.025,
                    PCs[1], PCs[2])
mass_col <- c(mass_col1, mass_col2)

g5c <- ggplot(data = vr_mass, aes(Type, Null_values)) + 
  theme_bw() + xlab("Species group") + ylab("Median body mass (g)") + 
  geom_violin(fill = cbbPalette[1],  bw = 3) +
  stat_summary(fun.data=mean_sdl, geom = "crossbar", width=0.05,
               fill=cbbPalette[5]) +
  geom_point(data =  vrOB_mass, aes(Type, Obs), size = PS,
             shape = 18, col = mass_col) + ggtitle("d) Body mass")

#null model results
null_tab[4, 2:3] <- mass_run[[2]][[3]] #extant
null_tab[5, 2:3] <- mass_run[[4]][[3]] #extinct

###HABITAT##############################
#hists all normal
hab_run <- null_volancy(datExtant, datAll, zp_type = "ES", 
                        trait = "Habitat",
                        pool = Pool)

vr_hab <- hab_run[[1]][[1]]

vr_hab4 <- tidyr::pivot_longer(vr_hab, cols = Forest:Wetland,
                               names_to = "Habitat")

vrOB_hab <- data.frame("Obs" = as.vector(hab_run[[1]][[2]]),
                       "Habitat" = names(hab_run[[1]][[2]]))

#set point colour based on significance of ES
#order in the results table is Forest, Grass, Marine, Wet - which matches figure
hab_col <- sapply(hab_run[[1]][[3]][2,], function(x){
  ifelse(x > 0.975 | x < 0.025, PCs[1], PCs[2])
})

g5h <- ggplot(data = vr_hab4, aes(Habitat, value)) + 
  theme_bw() + xlab("Habitat") + ylab("Number of species") + 
  geom_violin(fill = cbbPalette[1],  bw = 5) +
  stat_summary(fun.data=mean_sdl, geom = "crossbar", width=0.05,
               fill=cbbPalette[5]) +
  geom_point(data =  vrOB_hab , aes(Habitat, Obs), size = PS,
             shape = 18, col = hab_col) + ggtitle("e) Habitat")

#null model results
null_tab[6:9, 2:3] <- hab_run[[1]][[3]] %>% t()


###TROPHIC LEVEL##############################
lev_run <- null_volancy(datExtant, datAll, zp_type = "ES", 
                        trait = "Trophic.Level",
                        pool = Pool)

#null model results
null_tab[10:12, 2:3] <- lev_run %>% t()

###TROPHIC NICHE##############################

niche_run <- null_volancy(datExtant, datAll, zp_type = "ES", 
                        trait = "Trophic.Niche",
                        pool = Pool)

vr_niche <- niche_run[[1]][[1]]

vr_niche4 <- tidyr::pivot_longer(vr_niche, cols = 'Aquatic predator':'Vertivore',
                               names_to = "Niche")

vrOB_niche <- data.frame("Obs" = as.vector(niche_run[[1]][[2]]),
                       "Niche" = names(niche_run[[1]][[2]]))

#set point colour based on significance of ES
#order in the results table is matches figure
tro_col <- sapply(niche_run[[1]][[3]][2,], function(x){
  ifelse(x > 0.975 | x < 0.025, PCs[1], PCs[2])
})


g5n <- ggplot(data = vr_niche4, aes(Niche, value)) + 
  theme_bw() + xlab("Niche") + ylab("Number of species") + 
  geom_violin(fill = cbbPalette[1],  bw = 5) +
  stat_summary(fun.data=mean_sdl, geom = "crossbar", width=0.05,
               fill=cbbPalette[5]) +
  geom_point(data =  vrOB_niche, aes(Niche, Obs), size = PS,
             shape = 18, col = tro_col) + ggtitle("f) Trophic niche")

#null model results
null_tab[13:19, 2:3] <- niche_run[[1]][[3]]  %>% t()

######beak hypervolume############################

#In the source code, we load in the null values for all and isl species pools
#TO GENERATE THESE FILES FROM SCRATCH, see the code below in the hypervolume section

if (Pool == "All"){
  beak_null_vals <- as.vector(unlist(beak_run_all[[1]]))
  beak_null_df <- data.frame("Type" = rep("Extant", length(beak_null_vals)),
                             "value" = beak_null_vals)
  beak_obs <- data.frame("Type" = "Extant", "value" = beak_run_all[[2]])
} else {
  beak_null_vals <- as.vector(unlist(beak_run_isl[[1]]))
  beak_null_df <- data.frame("Type" = rep("Extant", length(beak_null_vals)),
                             "value" = beak_null_vals)
  beak_obs <- data.frame("Type" = "Extant", "value" = beak_run_isl[[2]])
}

beak_ZP <- zP(beak_null_vals, beak_obs$value, type = "ES")

#set point colour based on significance of ES
beak_col <- ifelse(beak_ZP[2] > 0.975 | beak_ZP[2] < 0.025,
                  PCs[1], PCs[2])

g5B <- ggplot(data = beak_null_df, aes(Type, value)) + 
  theme_bw() + xlab("") + ylab("Volume of hypervolume") + 
  geom_violin(fill = cbbPalette[1],  bw = "nrd") +
  stat_summary(fun.data=mean_sdl, geom = "crossbar", width=0.05,
               fill=cbbPalette[5]) +
  geom_point(data =  beak_obs, aes(Type, value), size = PS,
             shape = 18, col = beak_col) + ggtitle("c) Beak variation")

#null model results
null_tab[20, 2:3] <- beak_ZP

###################################################
######FOUR-PART TRAIT FIGURE#################
##############################################

gl <- list(g5a, g5w, g5B, g5c, g5h, g5n)

lay <- rbind(c(1,1,2,3),
             c(4,4,5,5),
             c(6,6,6,6))

gridExtra::grid.arrange(grobs = gl, layout_matrix = lay)


#null model results table
null_tab <- as.data.frame(null_tab)

##########################################################
################THREAT ANALYSES: FIGURE 3#############################
#################################################################

##the data (ll) is loaded in in the source code

##split ll into a vector of species names,
#and a list of threats
llSN <- vapply(ll, function(x) x[[1]], FUN.VALUE = character(1))
ll <- lapply(ll, function(x) x[[2]])

#create filtered versions of ll and dat for just island endemics,
# and within these, CR, EN and VU and all threatened
wCR <- which(datExtant$IslandEndemic == "Endemic" & datExtant$status == "CR")
wEN <- which(datExtant$IslandEndemic == "Endemic" & datExtant$status == "EN")
wVU <- which(datExtant$IslandEndemic == "Endemic" & datExtant$status == "VU")
wTH <- c(wCR, wEN, wVU)
wList <- list("CR" = wCR, "EN" = wEN, "VU" = wVU, "TH" = wTH)

#timing filters
ti_all <- c("Future","Ongoing","Past, Likely to Return","Past, Unlikely to Return",
            "Unknown")
#severity filters
#NOTE - a few severity records are coded (by IUCN) as NA
sev_all <- c("Causing/Could cause fluctuations","Negligible declines","No decline",
             "Rapid Declines","Slow, Significant Declines","Unknown",
             "Very Rapid Declines", NA)

#timing filters
ti_sub <- c("Ongoing")
#severity filters
#NOTE - a few severity records are coded (by IUCN) as NA
sev_sub <- c("Rapid Declines","Very Rapid Declines")

###run threats on each of these filters
#NOTE it does the matching due to threats and main dataset not having same
#alphabetical order
aLL <- lapply(wList, function(x){
  wDat <- datExtant[x,]
  high_hopes <- match(wDat$species, llSN)
  wLL <- ll[high_hopes]
  if (!identical(llSN[high_hopes],wDat$species)) stop("High hopes")
  threats(wLL, wDat, ti = ti_all, sev = sev_all)
})

#with subsets of timing and severity
#set scope to TRUE to test if you also filter out only major and moderate scope threats
aLL2 <- lapply(wList, function(x){
  wDat <- datExtant[x,]
  high_hopes <- match(wDat$species, llSN)
  wLL <- ll[high_hopes]
  if (!identical(llSN[high_hopes],wDat$species)) stop("High hopes")
  threats(wLL, wDat, ti = ti_sub, sev = sev_sub, scope = FALSE)
})

#sums of threats across species within each category
fre <- vapply(aLL, function(y) colSums(y[ ,1:11]), FUN.VALUE = numeric(11))
#version for subsets
fre2 <- vapply(aLL2, function(y) colSums(y[ ,1:11]), FUN.VALUE = numeric(11))

#remove the "other" (threat 12) column
gv3 <- data.frame("Threat" = 1:11, fre)
gv32 <- data.frame("Threat" = 1:11, fre2)

gv3_long <- tidyr::pivot_longer(gv3, cols = CR:TH,
                                    names_to = "Category")
gv32_long <- tidyr::pivot_longer(gv32, cols = CR:TH,
                                names_to = "Category")
#merge full and subset long dfs
gv_long <- rbind(gv3_long, gv32_long)

gv_long$Type <- c(rep("All",  nrow(gv3_long)),rep("Subset",  nrow(gv32_long)))

gv_long$Category <- factor(gv_long$Category, levels = c("CR", "EN", "VU", "TH"))


##for new stacked version need to remove TH
gv_long <- filter(gv_long, Category != "TH")

#make horizontal grouped bar plot
g6 <- ggplot(gv_long) + geom_col(aes(Threat, value, fill = Category),
                            position = "stack") +
  theme_bw() + xlab("Threat") + ylab("Number of species") +
  labs(fill="IUCN Category") + scale_fill_manual(values=cbbPalette[c(4,3,1)]) +
  scale_x_discrete(limits = factor(1:11)) +
  facet_wrap(~Type, scales = "free_y")

g6

##X2 analysis for each category vs. threat

#create the contingency table
thrt2 <- gv3[,c("CR", "EN", "VU")] %>%
  as.data.frame() %>%
  apply(2, as.numeric)
rownames(thrt2) <- gv32$Threat
#run test
chisq_thrt <- chisq.test(thrt2)
##calculate contribution of each cell to overall X2 statistic
contrib_thrt <- 100*chisq_thrt$residuals^2/chisq_thrt$statistic

#plot of residuals and contribution: FIGURE S1
par(mfrow=c(1,2))
corrplot::corrplot(chisq_thrt$residuals, is.cor = FALSE, cl.pos = "n")
corrplot::corrplot(contrib_thrt, is.cor = FALSE, cl.pos = "n")

#####################################################################
##Get country information and make a spatial map: a) threatened, b) extinct
##FIGURE 4
################################################################

#data (thrt_map) loaded in with source code

##Make test spatial heat map
world <- ne_countries(scale = "medium", returnclass = "sf")

#need to turn thrt_map long, but only non NA rows
thrt_map2 <- thrt_map %>%
  dplyr::select(GeographicalOrigin:LON5) 

thrt_map3 <- thrt_map2 %>%
  dplyr::select(GeographicalOrigin:LON)

for (i in 1:nrow(thrt_map2)){
  
  if (!is.na(thrt_map2$GeographicalOrigin2[i])){
    elas1 <- thrt_map2[i, c("GeographicalOrigin2", "LAT2", "LON2")]
    names(elas1) <-  c("GeographicalOrigin", "LAT", "LON")
    thrt_map3 <- rbind(thrt_map3, elas1)
  } else{
    next
  }
  
  if (!is.na(thrt_map2$GeographicalOrigin3[i])){
    elas2 <- thrt_map2[i, c("GeographicalOrigin3", "LAT3", "LON3")]
    names(elas2) <-  c("GeographicalOrigin", "LAT", "LON")
    thrt_map3 <- rbind(thrt_map3, elas2)
  } else{
    next
  }
  
  if (!is.na(thrt_map2$GeographicalOrigin4[i])){
    elas3 <- thrt_map2[i, c("GeographicalOrigin4", "LAT4", "LON4")]
    names(elas3) <-  c("GeographicalOrigin", "LAT", "LON")
    thrt_map3 <- rbind(thrt_map3, elas3)
  } else{
    next
  }
  
  if (!is.na(thrt_map2$GeographicalOrigin5[i])){
    elas4 <- thrt_map2[i, c("GeographicalOrigin5", "LAT5", "LON5")]
    names(elas4) <-  c("GeographicalOrigin", "LAT", "LON")
    thrt_map3 <- rbind(thrt_map3, elas4)
  } else{
    next
  }
  
}#eo for


##work out number of species in each archipelago
arch_N <- thrt_map3 %>%
  group_by(GeographicalOrigin) %>%
  summarise("N" = n())

#add lat and lon to arch_N for each unique archipelago
arch_M <- match(arch_N$GeographicalOrigin, thrt_map3$GeographicalOrigin)
arch_N$LAT <- thrt_map3$LAT[arch_M]
arch_N$LON <- thrt_map3$LON[arch_M]

breaks <- c(1, 10, 30, 60)
labs <- c("1", "10", "30", "60")

g8a <- ggplot(data = world) +
  geom_sf(color = "black", fill = "lightgreen") +
  geom_point(data = arch_N, aes(x = LON, y = LAT, size = N),  
             shape = 21, fill = "royalblue3", alpha = 0.7) +
  labs(size="No. of species") +
  theme(legend.position = "top") + xlab("") + ylab("") + 
  ggtitle("a)") +
  theme(plot.title = element_text(vjust = -8))+
  scale_size_continuous(range = c(2,12),
                        limits = c(1,64),
                        labels = labs,
                        breaks = breaks) + theme_classic()

##Plot B: Extinct Species Hotspot 

ext_new_isl <- filter(ext_new, IslandEndemic == "Yes")

#Group by island names, then count number of each island grp using n()
#and keep lat and lon columns
ext_new_isl_grp <- ext_new_isl %>% 
  group_by(GeographicalOrigin_New) %>%
  summarise(N = n(),
            LAT = unique(LAT),
            LON = unique(LON))

breaks <- c(1, 10, 20, 30)
labs <- c("1", "10", "20", "30")
g8b <- ggplot(data = world) +
  geom_sf(color = "black", fill = "lightgreen") +
  geom_point(data = ext_new_isl_grp, aes(x = LON, y = LAT, size = N),  
             shape = 21, fill = "orangered1", alpha = 0.7) +
  labs(size="No. of species") +
  theme(legend.position = "top") + xlab("") + ylab("") + ggtitle("b)") +
  theme(plot.title = element_text(vjust = -8)) +
  scale_size_continuous(range = c(2,12),
                        limits = c(1,32),
                        labels = labs,
                        breaks = breaks) + theme_classic()


gridExtra::grid.arrange(g8a, g8b, nrow = 2)

########################################################################################
####################HYPERVOLUME ANALYSES################################
#################################################################################

##Version of hypervolume used: 2.0.12
##Version of BAT used: 2.8.1

#Note - this code calculate results from scratch, but the results from
#the paper have been loaded in already (h_list)

##create new version of datExtant
traitExtant <- datExtant
                      
##Now we need to remove the five Kiwi species
w3 <- which(traitExtant$species %in% c("Apteryx australis", "Apteryx haastii",
                                       "Apteryx mantelli",
                                       "Apteryx owenii",
                                       "Apteryx rowi"))
traitExtant <- traitExtant[-w3,] 
if (nrow(traitExtant) != 10943) stop("Atom Heart Mother")

#log transform traits
traitExtant[,c("Mass",
               "Beak.Length.culmen",	"Beak.Length.nares", "Beak.Width",	"Beak.Depth",
               "Tarsus.Length", "Wing.Length", "Secondary1", "Tail.Length")] <- 
  apply(traitExtant[,c("Mass",
                       "Beak.Length.culmen",	"Beak.Length.nares", "Beak.Width",	"Beak.Depth",
                       "Tarsus.Length", "Wing.Length", "Secondary1", "Tail.Length")], 2, 
        function(x) log(x))


##################################################################
#######BODY SIZE CORRECTED TRAIT APPROACH (ALTERNATIVE)################
###########################################################################

#function to return residuals from trait-mass regression
#resp = morpho trait to use as a response
# trait_residuals <- function(traitExtant, resp){
#   modBlur <- lm(traitExtant[,resp] ~ traitExtant$Mass)
#   residBlur <- residuals(modBlur)
#   return(residBlur)
# }

#get residuals, create col names and join with traitExtant
# residTraits <- vapply(c("Beak.Length.culmen",	"Beak.Length.nares", "Beak.Width",
#                         "Beak.Depth",
#                      "Tarsus.Length", "Wing.Length", "Secondary1", "Tail.Length"),
#       function(x) trait_residuals(traitExtant, resp = x),
# FUN.VALUE = numeric(nrow(traitExtant)))
#
# colnames(residTraits) <- paste0(colnames(residTraits), "_resid")
#
# traitExtant <- cbind(traitExtant, residTraits) %>%
#   as.data.frame()

#Do the PCA axes - use center and scale as data not scaled here
# PCA <- prcomp(traitExtant[,c("Mass",
#                "Beak.Length.culmen_resid",	"Beak.Length.nares_resid", "Beak.Width_resid",
#                "Beak.Depth_resid", "Tarsus.Length_resid", "Wing.Length_resid",
#                "Secondary1_resid", "Tail.Length_resid")], center = TRUE, scale = TRUE)

##################################################################################################

PCA <- prcomp(traitExtant[,c("Mass",
                             "Beak.Length.culmen",	"Beak.Length.nares", "Beak.Width",
                             "Beak.Depth", "Tarsus.Length", "Wing.Length",
                             "Secondary1", "Tail.Length")], center = TRUE, scale = TRUE)

PCA_axes <- PCA$x 
rownames(PCA_axes) <- traitExtant$species

#for main analyses
p2 <- PCA_axes %>% as.data.frame() 
#add in other vars
p2$IslandEndemic <- traitExtant$IslandEndemic
p2$Threatened <- traitExtant$Threatened

##create subsets of data
#all global species that are not threatened island endemics
p2b <- rbind(filter(p2, IslandEndemic != "Endemic"),
             filter(p2, IslandEndemic == "Endemic" & Threatened == 0))

#filter out just island endemics
p3 <- filter(p2, IslandEndemic == "Endemic")
#non-threatened island endemics
p3b <- filter(p3, Threatened == 0)

#just threatened island endemics 
p4 <- filter(p3, Threatened == 1)

if (nrow(p2b) + nrow(p4) != nrow(p2)) stop("LZ4")


##########################################################
###########Using BAT for everything rather than hypervolume
##############################################################

#create a PA matrix with ncol of all world species and
#then rows for the three hypervolumes we build
sp_p2 <- rownames(p2)
mat_p2 <- matrix(0, nrow = 3, ncol = length(sp_p2))
colnames(mat_p2) <- sp_p2
rownames(mat_p2) <- c("nonThrtEnd", "NOthr", "thr")

sp_p2b <- rownames(p2b)
sp_p3b <- rownames(p3b)
sp_p4 <- rownames(p4)

mat_p2["nonThrtEnd", match(sp_p2b, sp_p2)] = 1
mat_p2["NOthr", match(sp_p3b, sp_p2)] = 1
mat_p2["thr", match(sp_p4, sp_p2)] = 1

if (sum(mat_p2["thr",]) != 526) stop("Us and them") #should be 530 minus the 4 threatened kiwis

#use kernel.build to build the volumes
BAT_hypers <- BAT::kernel.build(mat_p2, p2[,1:5], 
                                method = "svm", svm.gamma = 0.8)

bb_134 <- BAT::kernel.beta(BAT_hypers, func = "jaccard", comp = TRUE)

#calculate the unique fractions
union_14 <- sum(as.matrix(bb_134$Shared)["thr", "nonThrtEnd"], 
                as.matrix(bb_134$Unique_to_Cols)["thr", "nonThrtEnd"],
                as.matrix(bb_134$Unique_to_Rows)["thr", "nonThrtEnd"])

uni_14 <- c("nonThrtEnd" = as.matrix(bb_134$Unique_to_Cols)["thr", "nonThrtEnd"] /
              union_14,
            "thr" = as.matrix(bb_134$Unique_to_Rows)["thr", "nonThrtEnd"] /
              union_14)

hv_set34 <- hypervolume::hypervolume_set(BAT_hypers[[2]], BAT_hypers[[3]],
                                         check.memory=FALSE) #non thr isl end vs. thrt isl end
alp_34 <- BAT::kernel.alpha(hv_set34)

uni_34 <- c("NOthr" = alp_34["Unique component of (NOthr) relative to (thr)"] /
                 alp_34["Union of (NOthr, thr)"],
               "thr" = alp_34["Unique component of (thr) relative to (NOthr)"] /
                 alp_34["Union of (NOthr, thr)"])

#1 = all world's birds, 3 = non-threatened island endemics, 4 = threatened endemics
#the results in the paper all 34 (i.e. 1 is not focused on)
part_all <- matrix(c("Total14" = as.matrix(bb_134$Btotal)[3,1], 
                     "Repl14" = as.matrix(bb_134$Brepl)[3,1],
                     "Rich14" = as.matrix(bb_134$Brich)[3,1], 
                     "Total34" = as.matrix(bb_134$Btotal)[3,2], 
                     "Repl34" = as.matrix(bb_134$Brepl)[3,2],
                     "Rich34" = as.matrix(bb_134$Brich)[3,2]), 
                   ncol = 2, nrow = 3)
rownames(part_all) <- c("Total", "Repl", "Rich")
colnames(part_all) <- c("14", "34")
part_all #beta-div results 
uni_34 #unique fractions


#############################################################
###Hypervolume null modelling: Fig 8#############################
##############################################################

##Run using parallel processing. Note, with 999 iterations, it takes
#a long time to run.
#Change pool to 'isl' to use the island endemics species pool

# cores = 9
# cl = makeCluster(cores); on.exit(stopCluster(cl))
# registerDoParallel(cl)
# i = 1 #Dummy line for RStudio warnings
# 
# ##main parallel for loop (need to specify null_model)
# h0 = foreach(i=seq(from=1, to=999, by=1))  %dopar% { 
#   library(hypervolume)
#   library(dplyr)
#   library(BAT)
#   Fits <- hyper_null(p2, no_axe = 5)
#   Fits
# }
# 
# h = matrix(unlist(h0), ncol = 9, byrow = TRUE)
# colnames(h) <- names(h0[[1]])


#results of this already loaded in (h_list) through source code
#h_list: A list with five elements: i) the null model values, ii) unused volume calcs (ignore),
#iii) the raw hypervolumes, iv) beta diversity values, v) unique fractions
h <- h_list[[1]]; bb_134 <- h_list[[4]]; uni_34 <- h_list[[5]]

hh <-  as.data.frame(h) %>%
  select(jaccard34,  Uniq_Frac_34_3, Uniq_Frac_34_4)

hh_long <- tidyr::pivot_longer(hh, cols = jaccard34:Uniq_Frac_34_4,
                               names_to = "Metric")

hh_long$Metric <- factor(hh_long$Metric, 
  levels = c("jaccard34", "Uniq_Frac_34_3", "Uniq_Frac_34_4"))

#add observed values in from above
h_obs <- data.frame("Metric" = c("jaccard34", "Uniq_Frac_34_3", "Uniq_Frac_34_4"),
                    "value" = c(as.matrix(bb_134$Btotal)[3,2], uni_34))

h_obs$Metric <- factor(h_obs$Metric, 
                         levels = c("jaccard34", "Uniq_Frac_34_3", "Uniq_Frac_34_4"))

#null model results
zP_res_hh <- matrix(ncol = 2, nrow = ncol(hh))
colnames(zP_res_hh) <- c("Z", "P")

for (i in 1:ncol(hh)){
  obs_hh <- h_obs$value[i]
  zz_hh <- zP(hh[,i], obs_hh, type = "ES")
  zP_res_hh[i,] <- c(round(zz_hh[1], 2), round(zz_hh[2], 4))
}

#set point colour based on significance of ES
PCs <- c(cbbPalette[3], 666666)
hyp_col <- sapply(zP_res_hh[,2], function(x){
  ifelse(x > 0.975 | x < 0.025, PCs[1], PCs[2])
})

##build plot
g10 <- ggplot(data = hh_long, aes(Metric, value)) + 
  theme_bw() + xlab("Metric") + ylab("Unique fraction / Btotal") + 
  geom_violin(fill = cbbPalette[1],  bw = 0.05)  +
  stat_summary(fun.data=mean_sdl, geom = "crossbar", width=0.05,
               fill=cbbPalette[5]) +
  geom_point(data =  h_obs, aes(Metric, value), col = hyp_col, size = 4,
             shape = 18) +
  scale_x_discrete(labels=c("jaccard34" = "Btotal",
                            "Uniq_Frac_34_3" = "Unique-NTEs", 
                            "Uniq_Frac_34_4" = "Unique-TEs"))

g10


#############################################################################
#######BEAK SHAPE NULL MODELLING###########################################
##############################################################################

#do PCA on the four beak traits (non-residual versions but still log transformed)

##create new version of datExtant;
#this time we include the Kiwis as not outliers in terms of beak
traitExtant2 <- datExtant

#log-transform traits
traitExtant2[,c("Beak.Length.culmen",	"Beak.Length.nares", "Beak.Width",	"Beak.Depth")] <- 
  apply(traitExtant2[,c("Beak.Length.culmen",	"Beak.Length.nares", "Beak.Width",	"Beak.Depth")], 2, 
        function(x) log(x))
#do PCA
PCA_beak <- prcomp(traitExtant2[,c("Beak.Length.culmen",	"Beak.Length.nares", 
                                  "Beak.Width",
                                  "Beak.Depth")], center = TRUE, scale = TRUE)
PCA_beak_axes <- PCA_beak$x 
rownames(PCA_beak_axes) <- traitExtant2$species

#add the IE and thrt data to PCA axes
pBeak <- PCA_beak_axes %>% as.data.frame() 
pBeak$IslandEndemic <- traitExtant2$IslandEndemic
pBeak$Threatened <- traitExtant2$Threatened

#filter out thrt IEs
BeakIET <- filter(pBeak, IslandEndemic == "Endemic", Threatened == 1)

#create a PA matrix with ncol of all thrt IE sp, and one row
sp_beak <- rownames(BeakIET)
mat_beak <- matrix(1, nrow = 1, ncol = length(sp_beak))
colnames(mat_beak) <- sp_beak
rownames(mat_beak) <- c("thrt")

#we use 1.2 for svm.gamma in beaks
#use kernel.build to build the volumes
BAT_beak <- BAT::kernel.build(mat_beak, BeakIET[,1:4], 
                              method = "svm", svm.gamma = 1.2)
#hypervolume volume
alpha_beak <- BAT::kernel.alpha(BAT_beak) %>% as.vector()

##NULL MODELLING
##Run across parallel processing
cores = 9
cl = makeCluster(cores); on.exit(stopCluster(cl))
registerDoParallel(cl)
i = 1 #Dummy line for RStudio warnings

##main parallel for loop (need to specify null_model)
beak_run_all = foreach(i=seq(from=1, to=9, by=1))  %dopar% { 
  library(hypervolume)
  library(dplyr)
  library(BAT)
  bra <- beak_null(pBeak, pool = "All")
  bra
}
##THIS CAN THEN BE USED IN THE CODE ABOVE TO GENERATE THE FIGURE 
beak_run_all <- list(beak_run_all, alpha_beak)

###########################################################
#######Individual Trait Analyses############################
########################################################

#check species names in datExtant match those in trait dataset used in PCA
identical(datExtant$species, traitExtant2$species)

#PCA_beak calculated during the previous analyis 
datExtant2 <- cbind(datExtant, PCA_beak$x[,1:2]) %>% as.data.frame()

IE <- filter(datExtant2, IslandEndemic == "Endemic")

IE %>% 
  group_by(Threatened) %>%
  summarise("avg" = median(Mass), "avg2" = median(HWI),
            "avg3" = median(PC1))

#Trophic Niche
#Filter out 7 most species rich niches
IE_filt <- filter(IE, Trophic.Niche %in% c("Aquatic predator", "Frugivore",
                                           "Granivore", "Invertivore",
                                           "Nectarivore", "Omnivore","Vertivore"))

csTN <- rbind(table(filter(IE_filt, Threatened == 1)$Trophic.Niche),
              table(filter(IE_filt, Threatened == 0)$Trophic.Niche))
rownames(csTN) <- c("Thrt", "Non")
chisq.test(csTN)

#Habitat
#remove rock and human modified
IE_filt2 <- filter(IE, Habitat %in% c("Forest", "Grassland", "Marine", "Wetland"))

csH <- rbind(table(filter(IE_filt2, Threatened == 1)$Habitat),
             table(filter(IE_filt2, Threatened == 0)$Habitat))
rownames(csH) <- c("Thrt", "Non")
chisq.test(csH)

#beak
wilcox.test(PC1 ~ as.character(Threatened), data = IE, alternative = "two.sided")
wilcox.test(PC2 ~ as.character(Threatened), data = IE, alternative = "two.sided")

#mass
wilcox.test(Mass ~ as.character(Threatened), data = IE, alternative = "two.sided")

#HWI
wilcox.test(HWI ~ as.character(Threatened), data = IE, alternative = "two.sided")

