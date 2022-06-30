#######################################################################################
################SOURCE CODE###############################################
################################################################################

#devtools::install_github("ianmoran11/mmtable2")
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(mmtable2)
library(rredlist)
library(corrplot)
library(sf)
library("rnaturalearth")
library("rnaturalearthdata")
library(grid)
library(gridExtra)
library(ggstance)
library(ggpubr)
library(patchwork)
library(RColorBrewer)
library(BAT)
library(VGAM)

library(foreach)
library(doParallel)
library(cluster)


#http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
cbbPalette <- c("#009E73", "#56B4E9","#0072B2", "#000000", "#E69F00",  "#F0E442", 
                "#D55E00", "#CC79A7")

###################################################################
##Data loading
######################################################

#IUCN_Data <- list(ext_new, dat, ll, thrt_map, beak_run_all, beak_run_isl, h_list)

IUCN_Data <- readRDS(file.choose()) #load IUCN_Data.rds

##extinct species dataset
ext_new <- IUCN_Data[[1]]

##extant species dataset
dat <- IUCN_Data[[2]]

##species threat data
#a list where each element containts a list with two elements,
#the first being species name and the second being the threat info from IUCN
ll <- IUCN_Data[[3]]

#map data
thrt_map <- IUCN_Data[[4]]

#beak morphospace (all pool) null model values
beak_run_all <- IUCN_Data[[5]]

#beak morphospace (isl pool) null model values
beak_run_all <- IUCN_Data[[6]]

#Hypervolume beta-div and null model results
#A list with five elements: i) the null model values (all species pool),
#ii) the null model values (island endemic species pool),
#iii) the raw hypervolumes, iv) beta diversity values, v) unique fractions
h_list <- IUCN_Data[[7]]
  
###########################################################################################


###################################################################
##Data Formatting
######################################################


#change the 5 extinct in the wild species to extinct
ext_new$status <- vapply(ext_new$status, function(x) ifelse(x == "EW", "EX", x), 
                         FUN.VALUE = character(1))

#check all species unique
if (length(unique(dat$species)) != nrow(dat)) stop("ATWA")

##for the purposes of this study, group habitat types.
#Forest and woodland = Forest
#Grassland and shrubland = Grassland
#Wetland, Coastal & Riverine = Wetland
dat$Habitat <- vapply(dat$Habitat, function(x){
  if (is.na(x)){
    "NA"
  } else if (x == "Coastal" | x == "Riverine"){
    "Wetland"
  } else if (x == "Woodland"){
    "Forest"
  } else if (x == "Shrubland"){
    "Grassland"
  } else{
    x
  }
}, FUN.VALUE = character(1)) %>%
  as.vector()

##remove data deficient species (DD) - NOTE This is before any analysis so reduces total richness
dat <- filter(dat, status != "DD")

#select same columns from both
datJR <- select(dat, Order, species, status, IslandEndemic, Volancy, 
                Mass, Habitat, Trophic.Level, Trophic.Niche, HWI, LAT,
                Beak.Length.culmen,	Beak.Length.nares,Beak.Width,	Beak.Depth,
                Tarsus.Length,Wing.Length,Secondary1, Tail.Length)

extJR <- select(ext_new, Order, species, status, IslandEndemic, Volancy,
                Mass, Habitat, Trophic.Level, Trophic.Niche, HWI,  LAT,
                Beak.Length.culmen,	Beak.Length.nares,Beak.Width,	Beak.Depth,
                Tarsus.Length,Wing.Length,Secondary1, Tail.Length)

if(!identical(colnames(datJR), colnames(extJR))) stop("Mother Nature Son")

#merge
datAll <- rbind(datJR, extJR) %>% as.data.frame()

#11 species classified as weak flyer (1 extant, 10 extinct), convert to Flightless 
wvE <- which(datAll$Volancy == "Weak flyer")
datAll$Volancy[wvE] <- "Flightless"

#create extinct column (1 for extinct, 0 for not)
datAll <- mutate(datAll, Extinct = ifelse(status == "EX", 1, 0))

#turn status into factor with correct order for bar plot
datAll$status <- factor(datAll$status,levels = c("EX", "CR", "EN", "VU", "NT", "LC"))

#change Island endemic column to say Endemic - Non-endemic
datAll$IslandEndemic <- vapply(datAll$IslandEndemic, function(x){
  switch(x,
         "Yes" = "Endemic",
         "No" = "Non-endemic")
}, FUN.VALUE = character(1))


#split full dataset into island endemics and not
end <- filter(datAll, IslandEndemic == "Endemic")
noend <- filter(datAll, IslandEndemic == "Non-endemic")

##create new dataset version without LC or NT
datAll2 <- filter(datAll, !status %in% c("LC", "NT"))


##create new dataset by filtering out just extant species
datExtant <- filter(datAll, status != "EX")

##create a column of just threatened
datExtant <- mutate(datExtant, 
                    Threatened = ifelse(status == "VU" | status == "EN" | status == "CR", 1, 0))


###########################################################
#######FUNCTIONS##########################################
#########################################################################


##Get Z-scores for each null distribution (and 2-tailed P-value)
#if null distributions roughly normal use standard P-value and SES,
#if not use the Lhotsky ES (type == "ES") and 2-sided p-value


zP <- function(dis, obs, type = "ES"){
  if (type == "SES"){
    z <- (obs - mean(dis)) / sd(dis)
    p <- 2*pnorm(-abs(z))
  } else if (type =="ES"){
    dis <- c(obs, dis)# don't forget to add the obs into the null values
    res <- (sum(dis < obs) + sum(dis == obs)/2) / length(dis)
    z <- VGAM::probitlink(res)
    p <- res
  }
  return(c(z,p))
}



#null model function to randomly sample the main dataset
null_mod <- function(d, n, type = "status", numer = FALSE){
  #d = main dataset, n = number of species to sample, type = to sample status,
  #order or volant status etc.
  #numer = is the column you are sampling numeric or categorical - if numeric,
  #don't create a table, just return raw vector
  #returns number of n species in each of the categories
  s <- sample(d[ ,type], n, replace = FALSE)
  if (!numer){
    td2 <- table(s)
  } else {
    td2 <- s
  }
  return(td2)
}


####################################################
##ORDERS
###################################################################

  #function to extract the number of a given species from every element in
  #the replicate list
  #li = replicate list; on = order name
extract_list <- function(li, on){
  vv <- vapply(li, function(y){
    w <- which(names(y) == on)
    if (length(w) == 0){
      0
    } else{
      as.vector(y[w])
    }
  }, FUN.VALUE = numeric(1))
  return(vv)
}

  ##function to run the null model 9999 times, and extract n species in each
  #go, where n = number of CR species, or EN species etc,
  #Note - rowsums in ord_mat wont necessarily = n, because some sampled orders
  #will not be in Ord_Names (i.e. they will be in orders with < 150 species), but
  #this is fine
  null_order <- function(dat, n, Ord_Names){
    
    null999_order <- replicate(9999, null_mod(dat, n, type = "Order")) 
    
    #run these names through the extract_list function to get the null distributions
    #for each
    ord_mat <- matrix(NA, ncol = length(Ord_Names), nrow = 9999)
    colnames(ord_mat) <- Ord_Names
    
    for (i in 1:length(Ord_Names)){
      ord_mat[,i] <- extract_list(null999_order, on = Ord_Names[i])
    }
    
    return(ord_mat)
    
  }
  
  ##function to plot the null distributions for each order, alongside
  #the observed values
  #Arguments: n = number of CR species, en species etc; type = IUCN
  #category; bonferroni for p-value correction
  plot_null_order <- function(dat, n, Ord_Names, type = "CR",
                              bonferroni = FALSE, zp_type = "ES"){
    
    ord_mat <- null_order(dat, n, Ord_Names)
    ord_df <- data.frame(ord_mat)
    cns <- colnames(ord_df)
    
    #get observed order numbers within type
    if (type != "TH"){
      datIE <- filter(dat, IslandEndemic == "Endemic" & status == type)
    } else {
      datIE <- filter(dat, IslandEndemic == "Endemic" & 
                        status %in% c("CR", "EN", "VU"))
    }
    if (nrow(datIE) != n) stop("BYOB")
    tab_type_order <- table(datIE$Order)
    #remove order names not in Ord_Names (i.e. not in most speciose orders)
    ord_rem <- which(!names(tab_type_order) %in% Ord_Names)
    tab_type_order <- tab_type_order[-ord_rem]
    #create df and match up order column with tab_type_order
    type_order_obs <- data.frame("Order" = cns, "Number" = 0)
    mm <- match(names(tab_type_order), type_order_obs$Order)
    type_order_obs$Number[mm] <- as.vector(tab_type_order)
    
    #make long version of ord_mat for ggplot
    null999_order_long <- tidyr::pivot_longer(ord_df, cols = cns[1]:cns[length(cns)],
                                              names_to = "Order")
    
    un_null <- unique(null999_order_long$Order)
    zp_null <- matrix(nrow = length(un_null), ncol = 2)
    rownames(zp_null) <- un_null
    colnames(zp_null) <- c("Z", "P")
    for (j in 1:length(un_null)){
      dums <- filter(null999_order_long, Order == un_null[j])
      obs_null <- filter(type_order_obs, Order == un_null[j])["Number"] %>%
        unlist() %>%
        as.vector()
      zp_null[j,] <- zP(dums$value, obs_null, zp_type) %>% round(3)
    }
    zp_null <- as.data.frame(zp_null)
    
    #do bonferroni correction on p-value from either SES
    #or ES approaches and determine whether a p-value is sig or not
    if (bonferroni & zp_type == "SES"){
      alpha <- 0.05 / nrow(zp_null)
      zp_null$sig <- ifelse(zp_null$P < alpha, "Sig.", "Non-sig.")
    } else if ((!bonferroni) & zp_type == "SES"){
      alpha <- 0.05 
      zp_null$sig <- ifelse(zp_null$P < alpha, "Sig.", "Non-sig.")
    } else if (bonferroni & zp_type == "ES"){
      alpha1 <- (0.05 / nrow(zp_null)) / 2
      alpha2 <- 1 - alpha1
      zp_null$sig <- ifelse(zp_null$P < alpha1 | zp_null$P > alpha2, "Sig.", "Non-sig.")
    } else if ((!bonferroni) & zp_type == "ES"){
      alpha1 <- 0.025
      alpha2 <- 1 - alpha1
      zp_null$sig <- ifelse(zp_null$P < alpha1 | zp_null$P > alpha2, "Sig.", "Non-sig.")
    }
    return(list(null999_order_long,  type_order_obs, zp_null))
  }
  
  
  
  #########################################################
  ####TRAIT NULL MODELS###########################
  ################################################################
  
  #function to run the null model, and return the null values, observed value
  #and z-score and p-value
  #dataset = datExtant, or datALL
  #n = number of threatened (or extinct) sp to sample in numm model
  #obs = observed no. of flightless species for that dataset
  #zp_type = use the SES or ES effect size for Z and P
  
  null_volancy_int <- function(dataset, n, obs, zp_type = "ES", trait = "volancy"){
    
    if (trait == "volancy"){
      null999_volancy <- replicate(9999, null_mod(dataset, n, type = "Volancy"), 
                                   simplify = FALSE) 
      #if no flightless species in a given sample, return 0   
      vv_volancy <- vapply(null999_volancy, function(y){
        w <- which(names(y) == "Flightless")
        if (length(w) == 0){
          0
        } else{
          as.vector(y["Flightless"])
        }
      }, FUN.VALUE = numeric(1))
      
      #z and p
      zp_null <- zP(vv_volancy, obs, zp_type) %>% round(3)
      
      #results list
      sun1 <- list(vv_volancy, obs, zp_null)
      return(sun1)
    }#eo if volancy
    
    if (trait == "Trophic.Level"){
      null999_level <- replicate(9999, null_mod(dataset, n, type = "Trophic.Level"), 
                                 simplify = FALSE) 
      #if no species in a given level, return 0   
      vv_level <- lapply(null999_level, function(y){
        if (!all(c("Carnivore", "Herbivore", "Omnivore", "Scavenger") %in% names(y))){
          dum_vec <- c("Carnivore" = 0, "Herbivore" = 0, "Omnivore" = 0, 
                       "Scavenger" = 0)
          dum_match <- match(names(y), names(dum_vec))
          dum_vec[dum_match] <- y
          dum_vec
        } else{
          y
        }
      })
      #check all colnames are in correct order
      if (!all(sapply(vv_level, 
                      function(x) identical(names(x), 
                                            c("Carnivore", "Herbivore", "Omnivore", "Scavenger"))))){
        stop("zeppelin")
      }
      #turn into dataframe
      vv_level_df <- do.call(rbind.data.frame, vv_level)
      colnames(vv_level_df) <- names(vv_level[[1]])
      #get rid of scavenger
      vv_level_df <- vv_level_df[,c("Carnivore", "Herbivore", "Omnivore")]
      #no scavengers in obs
      if (!identical(names(obs), colnames(vv_level_df))) stop("Led")
      #run null model (all look normal)
      zp_null <- sapply(1:3, function(x) zP(vv_level_df[,x], obs[x], zp_type)) %>%
        round(3)
      rownames(zp_null) <- c("Z", "P")
      colnames(zp_null) <- names(obs)
      return(zp_null)
    }#eo if trophic level
    
    if (trait == "Trophic.Niche"){
      null999_niche <- replicate(9999, null_mod(dataset, n, type = "Trophic.Niche"), 
                                 simplify = FALSE) 
      #if no species in a given level, return 0   
      vv_niche <- lapply(null999_niche, function(y){
        if (!all(c("Aquatic predator", "Frugivore", "Granivore", "Herbivore aquatic",
                   "Herbivore terrestrial", "Invertivore", "Nectarivore", 
                   "Omnivore", "Scavenger", "Vertivore") %in% names(y))){
          dum_vec <- c("Aquatic predator" = 0, "Frugivore" = 0, "Granivore" = 0, 
                       "Herbivore aquatic" = 0,
                       "Herbivore terrestrial" = 0, "Invertivore" = 0, 
                       "Nectarivore" = 0, 
                       "Omnivore" = 0, "Scavenger" = 0, "Vertivore" = 0)
          dum_match <- match(names(y), names(dum_vec))
          dum_vec[dum_match] <- y
          dum_vec
        } else{
          y
        }
      })
      #check all colnames are in correct order
      if (!all(sapply(vv_niche, 
                      function(x) identical(names(x), 
                                            c("Aquatic predator", "Frugivore", "Granivore", "Herbivore aquatic",
                                              "Herbivore terrestrial", "Invertivore", "Nectarivore", 
                                              "Omnivore", "Scavenger", "Vertivore"))))){
        stop("zeppelin")
      }
      #turn into dataframe
      vv_niche_df <- do.call(rbind.data.frame, vv_niche)
      colnames(vv_niche_df) <- names(vv_niche[[1]])
      #get rid of species poor niches
      vv_niche_df <- vv_niche_df[,c("Aquatic predator", "Frugivore", "Granivore",
                                    "Invertivore", "Nectarivore", "Omnivore", "Vertivore")]
      if (!identical(names(obs), colnames(vv_niche_df))) stop("Led")
      #run null model (all look normal)
      zp_null <- sapply(1:7, function(x) zP(vv_niche_df[,x], obs[x], zp_type)) %>%
        round(3)
      rownames(zp_null) <- c("Z", "P")
      colnames(zp_null) <- names(obs)
      sun1 <- list(vv_niche_df, obs, zp_null)
      return(sun1)
    }#eo if trophic niche
    if (trait == "Mass"  | trait == "HWI" | trait == "RAN"){
      null999_mass <- replicate(9999, null_mod(dataset, n, type = trait, numer = TRUE), 
                                simplify = FALSE) 
      
      dum_vec <- vapply(null999_mass, median, FUN.VALUE = numeric(1))
      #histogram skewed (for Mass) so use effect size
      zp_null <- zP(dum_vec, obs, zp_type) %>% round(3)
      sun1 <- list(dum_vec, obs, zp_null)
      return(sun1)
    }
    
    if (trait == "Habitat"){
      null999_habitat <- replicate(9999, null_mod(dataset, n, type = "Habitat"), 
                                   simplify = FALSE) 
      #if no species in a given level, return 0   
      vv_habitat <- lapply(null999_habitat, function(y){
        if (!all(c("Forest", "Grassland", "Marine", 
                   "Wetland") %in% names(y))){
          
          dum_y <- c("Forest" = 0, "Grassland" = 0, "Marine" = 0, 
                     "Wetland" = 0)
          
          whi_y <- match(names(y), names(dum_y))
          
          for (i in 1:length(whi_y)){
            if(!is.na(whi_y[i])){
              dum_y[whi_y[i]]  <- y[i]
            }
          }
        } else{
          dum_y <- y[c("Forest", "Grassland", "Marine", 
                       "Wetland")]
        }
        dum_y
      })
      #check all colnames are in correct order
      if (!all(sapply(vv_habitat, 
                      function(x) identical(names(x), 
                                            c("Forest", "Grassland", "Marine", 
                                              "Wetland"))))){
        stop("zeppelin")
      }
      #turn into dataframe
      vv_habitat_df <- do.call(rbind.data.frame, vv_habitat)
      colnames(vv_habitat_df) <- names(vv_habitat[[1]])
      if (!identical(names(obs), colnames(vv_habitat_df))) stop("Led")
      #run null model (all look normal)
      zp_null <- sapply(1:4, function(x) zP(vv_habitat_df[,x], obs[x], zp_type)) %>%
        round(3)
      rownames(zp_null) <- c("Z", "P")
      colnames(zp_null) <- names(obs)
      sun1 <- list(vv_habitat_df, obs, zp_null)
      return(sun1)
    }#eo if habitat
    
  }#eo function
  
  #pool = whether the source pool for the random draws should be "All" extant
  #species or just "Isl" endemic extant species
  #terrestrial = whether we are using the full datExtant (F) or removing marine species
  #(TRUE) for the HWI analysis
  
  null_volancy <- function(datExtant, datAll, zp_type = "ES", trait = "volancy",
                           pool = "All", terrestrial = FALSE){
    
    #all extant
    n1 <- table(datExtant$Threatened)["1"] %>% as.vector()
    
    ##then fun just for island endemics
    datIE <- filter(datExtant, IslandEndemic == "Endemic")
    datIE_thrt <- filter(datIE, Threatened == 1)#just threatened island endemics
    n2 <- nrow(datIE_thrt)
    if (n2 != 530 & (!terrestrial)) stop("n2 wrong") #should equal 530
    
    ##then for extinct species
    n3 <- filter(datAll, IslandEndemic == "Endemic", status == "EX") %>% nrow()
    if (n3 != 149) stop("n3 wrong") #should equal 149
    
    #set the main source pool datasets
    pool_dataset <- switch(pool,
                           "All" = datExtant,
                           "Isl" = datIE)
    
    pool_dataset_extinct <- switch(pool,
                                   "All" = datAll,
                                   "Isl" = filter(datAll, IslandEndemic == "Endemic"))
    
    #small function to re-run a certain null model exercise for 
    #CR, EN, VU individually
    #classif = 1 of CR, EN or VU
    int_fun <- function(datIE_thrt, trait = "Mass", datPool, zp_type,
                        classif = "CR"){
      
      dumd <- filter(datIE_thrt, status == classif)
      if (trait == "Mass" | trait == "HWI"){
        obsd <- median(dumd[,trait])
      } else{
        obsd <- table(dumd[,trait])
        if (trait == "Habitat"){
          if (length(obsd) != 4) stop("dazed and confused")
        } else if (trait == "Trophic.Niche"){
          obsd <-  obsd[c("Aquatic predator", "Frugivore", "Granivore",
                          "Invertivore", "Nectarivore", "Omnivore", "Vertivore")]
          if (length(obsd) != 7) stop("dazed and confused")
        }
      }
      nd <- nrow(dumd)
      null_volancy_int(datPool, nd, obsd, zp_type, trait = trait)
    }
    
    if (trait == "volancy"){
      
      #all: get observed: number of threatened flightless species
      obs_volancy <- filter(datExtant, Threatened == 1 & Volancy == "Flightless") %>%
        nrow()
      #run first time on full dataset
      sun1 <- null_volancy_int(datExtant, n1, obs_volancy, zp_type, trait = trait)
      
      #:isl end.: no. of threatened flightless island endemics
      obs_volancy_IE <- filter(datIE, Threatened == 1 & Volancy == "Flightless") %>%
        nrow()
      sun2 <- null_volancy_int(pool_dataset, n2, obs_volancy_IE, zp_type, trait = trait)
      
      #extinct
      obs_volancy_EX <- filter(datAll, IslandEndemic == "Endemic",
                               status == "EX" & Volancy == "Flightless") %>%
        nrow()
      sun3 <- null_volancy_int(pool_dataset_extinct, n3, obs_volancy_EX, zp_type, 
                               trait = trait)
      return(list(sun1, sun2, sun3))
      
    } else if (trait == "Trophic.Level"){#eo if volancy
      
      obs_level_IE <- table(datIE_thrt$Trophic.Level)#has no scavengers
      sun2 <- null_volancy_int(pool_dataset, n2, obs_level_IE, zp_type, trait = trait)
      return(sun2)
      
    } else if (trait == "Trophic.Niche"){
      
      obs_niche_IE <- table(datIE_thrt$Trophic.Niche)
      obs_niche_IE <-  obs_niche_IE[c("Aquatic predator", "Frugivore", "Granivore",
                                      "Invertivore", "Nectarivore", "Omnivore", "Vertivore")]
      sun2 <- null_volancy_int(pool_dataset, n2, obs_niche_IE, zp_type, trait = trait)
      #CR, EN and VU
      sun3 <- lapply(c("CR", "EN", "VU"), function(x){
        int_fun(datIE_thrt, trait = "Trophic.Niche", pool_dataset, zp_type,
                classif = x)
      })
      return(list(sun2, sun3))
      
    } else if (trait == "Habitat"){
      
      obs_habitat_IE <- table(datIE_thrt$Habitat)
      sun2 <- null_volancy_int(pool_dataset, n2, obs_habitat_IE, zp_type, trait = trait)
      #CR, EN and VU
      sun3 <- lapply(c("CR", "EN", "VU"), function(x){
        int_fun(datIE_thrt, trait = "Habitat", pool_dataset, zp_type,
                classif = x)
      })
      return(list(sun2, sun3))
      
    } else if (trait == "Mass"){
      #all thrt sp
      obs_mass_all <- median(filter(datExtant, Threatened == 1)$Mass)
      sun1 <- null_volancy_int(datExtant, n1, obs_mass_all, zp_type, trait = trait)
      #thrt isl end
      obs_mass_IE <- median(datIE_thrt$Mass)
      sun2 <- null_volancy_int(pool_dataset, n2, obs_mass_IE, zp_type, trait = trait)
      #CR, EN and VU
      sun3 <- lapply(c("CR", "EN", "VU"), function(x){
        int_fun(datIE_thrt, trait = "Mass", datPool = pool_dataset, zp_type,
                classif = x)
      })
      #extinct
      obs_mass_EX <- filter(datAll, IslandEndemic == "Endemic",
                            status == "EX")$Mass %>%
        median()
      sun4 <- null_volancy_int(pool_dataset_extinct, n3, obs_mass_EX, zp_type, trait = trait)
      return(list(sun1, sun2, sun3, sun4))
      
    } else if (trait == "HWI"){
      #all thrt sp
      obs_hwi_all <- median(filter(datExtant, Threatened == 1)$HWI)
      sun1 <- null_volancy_int(datExtant, n1, obs_hwi_all, zp_type, trait = trait)
      #thrt isl end
      obs_hwi_IE <-  median(datIE_thrt$HWI)
      sun2 <- null_volancy_int(pool_dataset, n2, obs_hwi_IE, zp_type, trait = trait)
      #CR, EN and VU
      sun3 <- lapply(c("CR", "EN", "VU"), function(x){
        int_fun(datIE_thrt, trait = "HWI", datPool = pool_dataset, zp_type,
                classif = x)
      })
      return(list(sun1, sun2, sun3))
    } else if (trait == "RAN"){
      obs_RAN_IE <-  median(datIE_thrt$RAN)
      sun2 <- null_volancy_int(pool_dataset, n2, obs_RAN_IE, zp_type, trait = trait)
      return(sun2)
    } else{
      stop("trait not recognised")
    }
  }#eo function
  
  
######################################################################
##################THREAT ANALYSIS##############################
  ###############################################################
  
  ##function to build threats matrix
  threats <- function(ll, dat, ti = ti_all, sev = sev_all, scope = FALSE){
    
    llM <- matrix(0, ncol = 12, nrow = length(ll))
    colnames(llM) <- 1:12
    rownames(llM) <- dat$species
    
    for (j in 1:length(ll)){
      
      dd <- ll[[j]]
      
      #if no threats skip
      if (length(dd) == 0) next

      #filter out timing and severity
      dd2 <- filter(dd, timing %in% ti)
      dd2 <- filter(dd2, severity %in% sev)
      
      if (scope){
        dd2 <- filter(dd2, scope %in% c("Majority (50-90%)", "Whole (>90%)"))
      }
      #if filtering has removed all threats
      if (nrow(dd2) == 0) next
      
      #extract code and turn into single number
      cd <- dd2$code
      cd2 <- sub("\\..*", "", cd) #use regex
      cd3 <- as.numeric(unique(cd2))
      llM[j, cd3] <- 1
    }#eo j
    return(llM)
  }#eo threats
  
  
  #############################################################################
  ##hypervolume null modelling
  ##########################################################################
  
  #3 hypervolumes: all non thrt-endemic species, all island endemic, all threatened isl endemic
  
  #Question: is volume and overlap of threatened island endemics bigger or smaller than
  #expected by chance
  
  #Method: Randomly partition the 1705 island endemic species from the full pool, then
  #randomly classify 530 of these as threatened. Then build the three hypervolumes,
  #extract volumes and overlaps.
  
  #Uses BAT package method
  
  #sv = gamma parameter of svm method
  #pool = all or isl (use all species pool or just island endemics);
  hyper_null <- function(p2, sv = 0.8, no_axe = 5, pool = "all"){
    
    if (pool == "isl") p2 <- filter(p2, IslandEndemic == "Endemic")
    
    dum <- p2
    
    #Randomly partition the 1707 island endemic species from the full pool
    #this just randomly samples island endemic status and renames those in dum
    newIE <- sample(dum$IslandEndemic, nrow(dum), replace = FALSE)
    dum$IslandEndemic <- newIE
    
    #filter out real island endemics from p2. Then randomly sample the 0s and 1s,
    #from the threatened column, and assign these to the threatened column for just
    #the island endemic rows in dum
    #So, randomly classify 526 of the island endemics in dum as threatened
    dumIE <- filter(p2, IslandEndemic == "Endemic") #just island endemics
    dumS <- sample(dumIE$Threatened, nrow(dumIE), replace = FALSE)
    dum[which(dum$IslandEndemic == "Endemic"), "Threatened"] <- dumS
    
    ##NOTE the above code does not keep the number of threatened species in the whole dataset
    #constant but it DOES keep the number of threatened sp within island endemics constant,
    #and that is all we are analysing now
    
    if (sum(dum[which(dum$IslandEndemic == "Endemic"), "Threatened"]) != 526 | 
        length(which(dum$IslandEndemic == "Endemic")) != 1702){
      stop ("Meddle")
    }
    
    ##Make the different datasets for hypervolume construction, using the random dum
    
    #all global species that are not threatened island endemics
    p2b <- rbind(filter(dum, IslandEndemic != "Endemic"),
                 filter(dum, IslandEndemic == "Endemic" & Threatened == 0))
    
    #filter out just island endemics
    p3 <- filter(dum, IslandEndemic == "Endemic")
    #non-threatened island endemics
    p3b <- filter(p3, Threatened == 0)
    
    #just threatened island endemics 
    p4 <- filter(p3, Threatened == 1)
    
    if (nrow(p2b) + nrow(p4) != nrow(dum)) stop("LZ4")
    
    #create a PA matrix with ncol of all world species and
    #then rows for the three hypervolumes we build
    sp_p2 <- rownames(dum)
    mat_p2 <- matrix(0, nrow = 3, ncol = length(sp_p2))
    colnames(mat_p2) <- sp_p2
    rownames(mat_p2) <- c("nonThrtEnd", "NOthr", "thr")
    
    sp_p2b <- rownames(p2b)
    sp_p3b <- rownames(p3b)
    sp_p4 <- rownames(p4)
    
    mat_p2["nonThrtEnd", match(sp_p2b, sp_p2)] = 1
    mat_p2["NOthr", match(sp_p3b, sp_p2)] = 1
    mat_p2["thr", match(sp_p4, sp_p2)] = 1
    
    #use kernel.build to build the volumes
    BAT_hypers <- BAT::kernel.build(mat_p2, dum[,1:no_axe], 
                                    method = "svm", svm.gamma = 0.8)
    
    
    bb_134 <- BAT::kernel.beta(BAT_hypers, func = "jaccard", comp = TRUE)
    
    
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
    
    part_all <- matrix(c("Total14" = as.matrix(bb_134$Btotal)[3,1], 
                         "Repl14" = as.matrix(bb_134$Brepl)[3,1],
                         "Rich14" = as.matrix(bb_134$Brich)[3,1], 
                         "Total34" = as.matrix(bb_134$Btotal)[3,2], 
                         "Repl34" = as.matrix(bb_134$Brepl)[3,2],
                         "Rich34" = as.matrix(bb_134$Brich)[3,2]), 
                       ncol = 2, nrow = 3)
    rownames(part_all) <- c("Total", "Repl", "Rich")
    colnames(part_all) <- c("14", "34")
    
    ##calculate hypervolume volume
    vol1 <-  BAT_hypers@HVList[[1]]@Volume #nonThrtEnd
    vol3 <-  BAT_hypers@HVList[[2]]@Volume #NOthr
    vol4 <-  BAT_hypers@HVList[[3]]@Volume #thr
    #check names are in right order
    if (!(BAT_hypers@HVList[[1]]@Name == "nonThrtEnd" &
          BAT_hypers@HVList[[2]]@Name == "NOthr" &
          BAT_hypers@HVList[[3]]@Name == "thr")) stop("SL2") 
    
    res <- c(as.matrix(bb_134$Btotal)[3,1], as.matrix(bb_134$Btotal)[3,2], 
             uni_14, uni_34,
             vol1, vol3, vol4)
    
    names(res) <- c("jaccard14", "jaccard34",
                    "Uniq_Frac_14_1", "Uniq_Frac_14_4",
                    "Uniq_Frac_34_3", "Uniq_Frac_34_4",
                    "Volume1", "Volume3", "Volume4")
    return(res)
  }
  
  #############################################################################
  ##Beak shape hypervolume null modelling
  ##########################################################################
  
  #pool = All or Isl
  beak_null <- function(pBeak, pool = "All"){
    
    dum <- pBeak
    
    if (pool == "All"){
      
      #Randomly partition the 1707 island endemic species from the full pool
      newIE <- sample(pBeak$IslandEndemic, nrow(pBeak), replace = FALSE)
      dum$IslandEndemic <- newIE
      
    } #else keep island endemic species as they are, but just randomize which are threatened
    
    #filter out real island endemics from pBeak. Then randomly sample the 0s and 1s,
    #from the threatened column, and assign these to the threatened column for just
    #the island endemic rows in dum
    #So, randomly classify 530 of the island endemics in dum as threatened
    dumIE <- filter(pBeak, IslandEndemic == "Endemic") #just island endemics
    dumS <- sample(dumIE$Threatened, nrow(dumIE), replace = FALSE)
    dum[which(dum$IslandEndemic == "Endemic"), "Threatened"] <- dumS
    
    #filter out thrt IEs
    BeakIET_dum <- filter(dum, IslandEndemic == "Endemic", Threatened == 1)
    
    #create a PA matrix with ncol of all thrt IE sp, and one row
    sp_beak_dum <- rownames(BeakIET_dum)
    mat_beak_dum <- matrix(1, nrow = 1, ncol = length(sp_beak_dum))
    colnames(mat_beak_dum) <- sp_beak_dum
    rownames(mat_beak_dum) <- c("thrt")
    
    #use kernel.build to build the volumes (takes ~ 40 seconds per run)
    BAT_beak <- BAT::kernel.build(mat_beak_dum, BeakIET_dum[,1:4], 
                                  method = "svm", svm.gamma = 1.2)
    
    #hypervolume volume
    alpha_beak_dum <- BAT::kernel.alpha(BAT_beak) %>% as.vector()
    return(alpha_beak_dum)
  }
  
  
  
  
  