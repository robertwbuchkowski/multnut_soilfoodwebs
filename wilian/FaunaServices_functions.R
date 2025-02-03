# Color scheme
# "Herbivory & Algivory", "Litter transformation","Bacterivory","Fungivory","Soil transformation","Predation"
FWPalette <- rev(c("chartreuse3","darkorange3","goldenrod1","goldenrod4","gray50","orangered2")) # "cadetblue2",


sum2 <- function(X) {
  if(length(X)==sum(is.na(X))) # if only NAs, return NA
  {
    return(sum(X, na.rm = F))
  } else # else, return the sum of non na's
    return(sum(X, na.rm = T))
}

dataprepare.coregroups <- function(Raw_Dataset,taxa_selection_included){

  # Use dplyr to select
  Raw_Dataset_selected <- Raw_Dataset %>%
    dplyr::select(biome,repeated,biomass,litter_fauna,extraction_method,sample_size,
                  datasetid,transect_id,latitude.y,longitude.y,country,lu_cover1,lu_cover2,transect_size,sample_year,sample_month,
                  temp,prec,elev,rad,evap,hii,climKG,biome_whitt,ecoregion_name,realm,
                  starts_with("d_"), starts_with("dbm_"), starts_with("fbm_")) # selecting all density and biomass

  Raw_Dataset_selected <- Raw_Dataset_selected %>%
    dplyr::select(-starts_with("d_se_"), -starts_with("dbm_se_"), -starts_with("fbm_se_")) # selecting all density and biomass


  ##remove columns with only 0s or NAs
  temp <- apply(Raw_Dataset_selected[,27:ncol(Raw_Dataset_selected)], 2, function(x) sum(is.na(x))!=nrow(Raw_Dataset_selected) | sum(x, na.rm = T )!=0)
  to.keep <- c(names(Raw_Dataset_selected)[1:26], names(temp)[which(temp==TRUE)])

  Raw_Dataset_selected <- Raw_Dataset_selected[,to.keep]
  Raw_Dataset_selected <- Raw_Dataset_selected[,-grep("unknown|others|other", names(Raw_Dataset_selected))]# removing others or unknown (because they will cause problems later and they are very unespecific)


  #remove lifestages that were not considered across all the studies
  Raw_Dataset_selected <- Raw_Dataset_selected[,-grep("pupae|eggs|cocoons", names(Raw_Dataset_selected))]

  # We do not take the information on adults or larvaes, just keep the total
  Raw_Dataset_selected <- Raw_Dataset_selected[,-grep("larvaes|adults", names(Raw_Dataset_selected))]

  #Remove heteroptera and homoptera as
  #Jerome already calculated the sum in d_hemiptera_total and fbm_hemiptera_total
  Raw_Dataset_selected<-subset(Raw_Dataset_selected,select=-c(d_heteroptera, d_homoptera))
  Raw_Dataset_selected<-subset(Raw_Dataset_selected,select=-c(fbm_heteroptera, fbm_homoptera))


  # sum of individual hymenoptera taxa (Jerome code)
  hymeno_sum_tmp <- apply(Raw_Dataset_selected[,c("d_ants","d_hymenoptera_not_ants","d_vespidae")],1,sum2)

  Raw_Dataset_selected$d_hymenoptera_total = if_else( is.na(hymeno_sum_tmp),Raw_Dataset_selected$d_hymenoptera,
                                                      if_else(
                                                        is.na(Raw_Dataset_selected$d_hymenoptera), hymeno_sum_tmp,
                                                        if_else(
                                                          Raw_Dataset_selected$d_hymenoptera >= hymeno_sum_tmp, Raw_Dataset_selected$d_hymenoptera, hymeno_sum_tmp
                                                        )
                                                      )
  )

  hymeno_sum_tmp <- apply(Raw_Dataset_selected[,c("fbm_ants","fbm_hymenoptera_not_ants")],1,sum2)
  Raw_Dataset_selected$fbm_hymenoptera_total = if_else( is.na(hymeno_sum_tmp),Raw_Dataset_selected$fbm_hymenoptera,
                                                        if_else(
                                                          is.na(Raw_Dataset_selected$fbm_hymenoptera), hymeno_sum_tmp,
                                                          if_else(
                                                            Raw_Dataset_selected$fbm_hymenoptera >= hymeno_sum_tmp, Raw_Dataset_selected$fbm_hymenoptera, hymeno_sum_tmp
                                                          )
                                                        )
  )

  ####Calculate arachnida_total as done as for hymenoptera_total
  # sum of individual hymenoptera taxa (Jerome code)
  arachnid_sum_tmp <- apply(Raw_Dataset_selected[,c("d_araneae","d_opiliones")],1,sum2)
  Raw_Dataset_selected$d_arachnida_total = if_else( is.na(arachnid_sum_tmp),Raw_Dataset_selected$d_arachnida,
                                                    if_else(
                                                      is.na(Raw_Dataset_selected$d_arachnida), arachnid_sum_tmp,
                                                      if_else(
                                                        Raw_Dataset_selected$d_arachnida >= arachnid_sum_tmp, Raw_Dataset_selected$d_arachnida, arachnid_sum_tmp
                                                      )
                                                    )
  )
  arachnid_sum_tmp <- apply(Raw_Dataset_selected[,c("fbm_araneae","fbm_opiliones")],1,sum2)
  Raw_Dataset_selected$fbm_arachnida_total = if_else( is.na(arachnid_sum_tmp),Raw_Dataset_selected$fbm_arachnida,
                                                      if_else(
                                                        is.na(Raw_Dataset_selected$fbm_arachnida), arachnid_sum_tmp,
                                                        if_else(
                                                          Raw_Dataset_selected$fbm_arachnida >= arachnid_sum_tmp, Raw_Dataset_selected$fbm_arachnida, arachnid_sum_tmp
                                                        )
                                                      )
  )

  # Transform to the long format
  Raw_Dataset_tr <- Raw_Dataset_selected %>%
    gather(-c(1:which(startsWith(colnames(Raw_Dataset_selected), "d_"))[1]-1), key="GroupFull", value="Value")

  Raw_Dataset_tr<-separate(Raw_Dataset_tr,GroupFull, c("Variable", "Group"), remove=FALSE,extra = "merge", fill = "right")

  ##remove NA cases
  Raw_Dataset_tr_f <- Raw_Dataset_tr[!is.na(Raw_Dataset_tr$Value),]

  Raw_Dataset_tr_f<-subset(Raw_Dataset_tr_f,select=-c(GroupFull))

  # convert to a wide format to have density and biomass in columns
  Raw_Dataset_wide <- Raw_Dataset_tr_f %>%
    pivot_wider(names_from="Variable", values_from="Value")

  # name the variables with units
  names(Raw_Dataset_wide)[names(Raw_Dataset_wide)=="d"] <- "Density_sqm"

  # convert biomass in grams (SI); biomass in mg
  Raw_Dataset_wide$fbm <- Raw_Dataset_wide$fbm/1000 # mg to g
  names(Raw_Dataset_wide)[names(Raw_Dataset_wide)=="fbm"] <- "Biomass_fresh_g_sqm"

  # and calculate average body mass
  Raw_Dataset_wide$BodyMass_g <- Raw_Dataset_wide$Biomass_fresh_g_sqm/Raw_Dataset_wide$Density_sqm

  # ISSUE: some monoliths with 0 density have biomass values
  Raw_Dataset_wide_NA <- subset(Raw_Dataset_wide, is.na(Density_sqm)==T)
  # No solution applied, the values are excluded
  Raw_Dataset_wide_noNA <- subset(Raw_Dataset_wide, is.na(Density_sqm)==F)
  Raw_Dataset_wide_noNA$BodyMass_g[is.nan(Raw_Dataset_wide_noNA$BodyMass_g)==T] <- NA
  Raw_Dataset_wide_noNA$BodyMass_g[Raw_Dataset_wide_noNA$BodyMass_g==Inf] <- NA

  Data_Transect <- Raw_Dataset_wide_noNA

  # Select transects that record ALL selected taxa (here only look at abundance, biomass is later)
  Check_Transect <- Data_Transect %>%
    group_by(transect_id,biome,lu_cover1) %>%
    dplyr::summarise(Key_taxa = all(taxa_selection_included %in% Group))

  # Select transects with all groups present
  Check_Transect_selection <- subset(Check_Transect,Key_taxa==T)

  # Select data from transects with all 'must' groups present
  Data_Transect_selection <- subset(Data_Transect, transect_id %in% Check_Transect_selection$transect_id)

  # Keep only data for our core group
  Data_Transect_selection <- subset(Data_Transect_selection, Group %in% taxa_selection_included)

  Data_Transect_selection
}


dataprepare.bodymasscore <- function(Raw_Dataset,taxa_selection_included){

  # Use dplyr to select
  Raw_Dataset_selected <- Raw_Dataset %>%
    dplyr::select(biome,repeated,biomass,litter_fauna,extraction_method,sample_size,
                  datasetid,transect_id,latitude.y,longitude.y,country,lu_cover1,lu_cover2,transect_size,sample_year,sample_month,
                  temp,prec,elev,rad,evap,hii,climKG,biome_whitt,ecoregion_name,realm,
                  starts_with("d_"), starts_with("dbm_"), starts_with("fbm_")) # selecting all density and biomass

  Raw_Dataset_selected <- Raw_Dataset_selected %>%
    dplyr::select(-starts_with("d_se_"), -starts_with("dbm_se_"), -starts_with("fbm_se_")) # selecting all density and biomass


  ##remove columns with only 0s or NAs
  temp <- apply(Raw_Dataset_selected[,27:ncol(Raw_Dataset_selected)], 2, function(x) sum(is.na(x))!=nrow(Raw_Dataset_selected) | sum(x, na.rm = T )!=0)
  to.keep <- c(names(Raw_Dataset_selected)[1:26], names(temp)[which(temp==TRUE)])

  Raw_Dataset_selected <- Raw_Dataset_selected[,to.keep]
  Raw_Dataset_selected <- Raw_Dataset_selected[,-grep("unknown|others|other", names(Raw_Dataset_selected))]# removing others or unknown (because they will cause problems later and they are very unespecific)


  #remove lifestages that were not considered across all the studies
  Raw_Dataset_selected <- Raw_Dataset_selected[,-grep("pupae|eggs|cocoons", names(Raw_Dataset_selected))]

  # We do not take the information on adults or larvaes, just keep the total
  Raw_Dataset_selected <- Raw_Dataset_selected[,-grep("larvaes|adults", names(Raw_Dataset_selected))]

  #Remove heteroptera and homoptera as
  #Jerome already calculated the sum in d_hemiptera_total and fbm_hemiptera_total
  Raw_Dataset_selected<-subset(Raw_Dataset_selected,select=-c(d_heteroptera, d_homoptera))
  Raw_Dataset_selected<-subset(Raw_Dataset_selected,select=-c(fbm_heteroptera, fbm_homoptera))


  # sum of individual hymenoptera taxa (Jerome code)
  hymeno_sum_tmp <- apply(Raw_Dataset_selected[,c("d_ants","d_hymenoptera_not_ants","d_vespidae")],1,sum2)

  Raw_Dataset_selected$d_hymenoptera_total = if_else( is.na(hymeno_sum_tmp),Raw_Dataset_selected$d_hymenoptera,
                                                      if_else(
                                                        is.na(Raw_Dataset_selected$d_hymenoptera), hymeno_sum_tmp,
                                                        if_else(
                                                          Raw_Dataset_selected$d_hymenoptera >= hymeno_sum_tmp, Raw_Dataset_selected$d_hymenoptera, hymeno_sum_tmp
                                                        )
                                                      )
  )

  hymeno_sum_tmp <- apply(Raw_Dataset_selected[,c("fbm_ants","fbm_hymenoptera_not_ants")],1,sum2)
  Raw_Dataset_selected$fbm_hymenoptera_total = if_else( is.na(hymeno_sum_tmp),Raw_Dataset_selected$fbm_hymenoptera,
                                                        if_else(
                                                          is.na(Raw_Dataset_selected$fbm_hymenoptera), hymeno_sum_tmp,
                                                          if_else(
                                                            Raw_Dataset_selected$fbm_hymenoptera >= hymeno_sum_tmp, Raw_Dataset_selected$fbm_hymenoptera, hymeno_sum_tmp
                                                          )
                                                        )
  )

  ####Calculate arachnida_total as done as for hymenoptera_total
  # sum of individual hymenoptera taxa (Jerome code)
  arachnid_sum_tmp <- apply(Raw_Dataset_selected[,c("d_araneae","d_opiliones")],1,sum2)
  Raw_Dataset_selected$d_arachnida_total = if_else( is.na(arachnid_sum_tmp),Raw_Dataset_selected$d_arachnida,
                                                    if_else(
                                                      is.na(Raw_Dataset_selected$d_arachnida), arachnid_sum_tmp,
                                                      if_else(
                                                        Raw_Dataset_selected$d_arachnida >= arachnid_sum_tmp, Raw_Dataset_selected$d_arachnida, arachnid_sum_tmp
                                                      )
                                                    )
  )
  arachnid_sum_tmp <- apply(Raw_Dataset_selected[,c("fbm_araneae","fbm_opiliones")],1,sum2)
  Raw_Dataset_selected$fbm_arachnida_total = if_else( is.na(arachnid_sum_tmp),Raw_Dataset_selected$fbm_arachnida,
                                                      if_else(
                                                        is.na(Raw_Dataset_selected$fbm_arachnida), arachnid_sum_tmp,
                                                        if_else(
                                                          Raw_Dataset_selected$fbm_arachnida >= arachnid_sum_tmp, Raw_Dataset_selected$fbm_arachnida, arachnid_sum_tmp
                                                        )
                                                      )
  )

  # Transform to the long format
  Raw_Dataset_tr <- Raw_Dataset_selected %>%
    gather(-c(1:which(startsWith(colnames(Raw_Dataset_selected), "d_"))[1]-1), key="GroupFull", value="Value")

  Raw_Dataset_tr<-separate(Raw_Dataset_tr,GroupFull, c("Variable", "Group"), remove=FALSE,extra = "merge", fill = "right")

  ##remove NA cases
  Raw_Dataset_tr_f <- Raw_Dataset_tr[!is.na(Raw_Dataset_tr$Value),]

  Raw_Dataset_tr_f<-subset(Raw_Dataset_tr_f,select=-c(GroupFull))

  # convert to a wide format to have density and biomass in columns
  Raw_Dataset_wide <- Raw_Dataset_tr_f %>%
    pivot_wider(names_from="Variable", values_from="Value")

  # name the variables with units
  names(Raw_Dataset_wide)[names(Raw_Dataset_wide)=="d"] <- "Density_sqm"

  # convert biomass in grams (SI); biomass in mg
  Raw_Dataset_wide$fbm <- Raw_Dataset_wide$fbm/1000 # mg to g
  names(Raw_Dataset_wide)[names(Raw_Dataset_wide)=="fbm"] <- "Biomass_fresh_g_sqm"

  # and calculate average body mass
  Raw_Dataset_wide$BodyMass_g <- Raw_Dataset_wide$Biomass_fresh_g_sqm/Raw_Dataset_wide$Density_sqm

  # ISSUE: some monoliths with 0 density have biomass values
  Raw_Dataset_wide_NA <- subset(Raw_Dataset_wide, is.na(Density_sqm)==T)
  # No solution applied, the values are excluded
  Raw_Dataset_wide_noNA <- subset(Raw_Dataset_wide, is.na(Density_sqm)==F)
  Raw_Dataset_wide_noNA$BodyMass_g[is.nan(Raw_Dataset_wide_noNA$BodyMass_g)==T] <- NA
  Raw_Dataset_wide_noNA$BodyMass_g[Raw_Dataset_wide_noNA$BodyMass_g==Inf] <- NA

  Data_Transect <- Raw_Dataset_wide_noNA

  #take all data with body mass for the groups
  Data_Transect_selectedGroup <- subset(Data_Transect, Group %in% taxa_selection_included)
  Data_Transect_selectedGroup <- subset(Data_Transect_selectedGroup,!is.na(BodyMass_g))

  Data_Transect_selectedGroup
}

bodymass.imputation<-function(Data_Transect_selection_noZeros,Data_Transect_noZeros_BodyMass,taxa_selection_included){

  listTransect_select<-unique(Data_Transect_selection_noZeros$transect_id)

  dataimput<-NULL
  for (i in taxa_selection_included){
    #select only data for a given group
    groupsub_bodymass<-subset(Data_Transect_noZeros_BodyMass,Group==i)
    groupsub_Tr_Selection<-subset(Data_Transect_selection_noZeros,Group==i)
    groupsub<-rbind(groupsub_bodymass,groupsub_Tr_Selection)
    groupsub<-groupsub %>% distinct()
    groupsub$biome<-as.factor(groupsub$biome)
    groupsub$climKG<-as.factor(groupsub$climKG)
    groupsub$lu_cover1<-as.factor(groupsub$lu_cover1)
    #data imputation
    ini <- mice(groupsub,m=2,maxit=0)
    pred <- ini$pred
    meth <- ini$meth
    print(meth)
    imp<-mice(groupsub,m=1,method = meth)
    datimput<-complete(imp,1)
    #select only data in transect_selection
    #(i.e. transects where all groups needed sampled)
    datimput <- subset(datimput, transect_id %in% listTransect_select)
    dataimput<-rbind(dataimput,datimput)
  }
  dataimput
}

# Function to build normal distributions for foodweb.preferences
min.f1f2 <- function(x, mu1, mu2, sd1, sd2) {
  f1 <- dnorm(x, mean=mu1, sd=sd1)
  f2 <- dnorm(x, mean=mu2, sd=sd2)
  pmin(f1, f2)
}

#Function to calcule trophic levels from food web energy fluxes
GetTL2 <- function(web){  # fonction pour calculer les niveaux trophiques
  ## transposer le réseau
  tweb <- t(web)

  ## La somme des lignes doit ?tre 1
  rs <- rowSums(tweb)
  for(i in 1:length(tweb[,1]))
    tweb[i,tweb[i,]>0] = tweb[i,tweb[i,]>0]/rs[i]

  nb.TL <- try(solve(diag(length(tweb[,1])) - tweb), T)

  if(class(nb.TL)[1]=="try-error")
    nbTL <- rep(NA, length(tweb[,1]))

  if(class(nb.TL)[1]!="try-error")
    nbTL <- rowSums(nb.TL)

  nbTL

}

foodweb.preferences<-function(Data,Raw_GuildsSoil, weight_by_biomass = TRUE){

  ### Assigning preferences and traits

  Data$Group[which(Data$Group=="coleoptera_total")]<-"coleoptera"
  Data$Group[which(Data$Group=="diptera_total")]<-"diptera"
  Data$Group[which(Data$Group=="lepidoptera_total")]<-"lepidoptera"
  Data$Group[which(Data$Group=="gastropoda_total")]<-"gastropoda"
  Data$Group[which(Data$Group=="hemiptera_total")]<-"hemiptera"
  Data$Group[which(Data$Group=="hymenoptera_total")]<-"ants"
  Data$Group[which(Data$Group=="arachnida_total")]<-"araneae"

  # Select guilds and adjjust traits
  Raw_GuildsSoil_selection <- subset(Raw_GuildsSoil, Group %in% unique(Data$Group))

  # Here we select trophic functions (resources) on which we focus:
  # Herbivory (=plants+algae), Detritivory (=litter+wood), Soil transformation (=soil), Microbivory (=bacteria+fungi)
  # In the case of macrofauna we just take proxies
  # plants = herbivory
  # litter = detritivory
  # soil = soil transformation
  # fungi = microbivory
  Raw_GuildsSoil_selection <- Raw_GuildsSoil_selection %>%
    select(-A, -W, -B) %>% ungroup()

  # Because we don't have microfauna, we set no predation in earthworms
  Raw_GuildsSoil_selection$Invt[Raw_GuildsSoil_selection$Group=="earthworms"] <- NA

  # Merging the traits and the data
  Data_traits <- merge(Data, Raw_GuildsSoil_selection, all.y = F)

  # Resource table is separate
  Data_Resources <- Raw_GuildsSoil[c(1,3,5,7),]
  Data_Resources <- Data_Resources %>%
    select(-A, -W, -B) %>% ungroup()

  # Adding arbitrary values for biomass (used only for illustration)
  Data_Resources$BodyMass_g_log10 <- 1
  Data_Resources$MassSD_log10 <- 1
  Data_Resources$PreyMean <- 1
  Data_Resources$PreySD <- 1


  ### Network reconstruction

  # To reconstruct the food web, we will generate a set of adjacency matrices based on different
  # parameters with food objects in rows and consumers in columns and multiply these matrices.

  # The following default assumptions were used for the exemplar reconstruction:
  #  (1) optimum PPMR was set to 100 (PPMRopt = log10(100)) (Brose et al., 2008)
  # body mass range of the optimum prey was set larger to the body mass range of the predator (PPMRwidth = 1.5)

  # Reconstruction coefficients, preparations
  # Fill data
  Data_traits$PPMRopt[is.na(Data_traits$PPMRopt)==T] <- 1
  Data_traits$PPMRwidth[is.na(Data_traits$PPMRwidth)==T] <- 1

  # Protection
  Data_traits$Predatory[is.na(Data_traits$Predatory)==T] <- 1
  Data_traits$Agility[is.na(Data_traits$Agility)==T] <- 1
  Data_traits$Metabolites[is.na(Data_traits$Metabolites)==T] <- 1
  Data_traits$Physical.protection[is.na(Data_traits$Physical.protection)==T] <- 1
  Data_traits$Protection <- Data_traits$Predatory*Data_traits$Agility*Data_traits$Metabolites*Data_traits$Physical.protection

  # we order everything alphabetically
  Data_traits <- Data_traits[order(Data_traits$transect_id,Data_traits$Group),]

  # we add log-transformed body mass to calculate size niche overlaps in log space later
  # mainly all taxa fall within 0.5 and 10 mg in fresh body mass
  Data_traits$BodyMass_g_log10 <- log10(Data_traits$BodyMass_g)

  # we add custom SD to the body size distribution (can be replaced with some empirical later)
  Data_traits$MassSD_log10 <- 1 # SD of +- order of magnitude in body mass

  ## Descriptive table

  # Raw data
  Data_All_save <- Data_traits
  Data_All_save[is.na(Data_All_save)] <- 0
  Data_All_save$PPMR[Data_All_save$PPMR=="0"] <- ""

  ## Food-web reconstruction

  # Function to reconstruct body size distribution using the mean and standard deviation of
  # living body mass and assuming it follows a log-normal pattern, which generally
  # observed both across species (Hutchinson & MacArthur, 1959; Blackburn & Gaston, 1994; Allen et al., 2006)
  # and within species (Stead et al., 2005; Potapov & et al., 2021).


  # Set reconstruction coefficients
  PPMRopt <- log10(100) # optimum predator-prey body mass ratio, default = 100
  PPMRwidth <- 1.5 # mutiplyer; breadth of the size spectrum of prey (standard deviation), default = sd of body mass of predator
  PPMRwidth_min <- 1 # minimum SD of body mass
  Protection <- 1 # mutiplyer; importance of protection mechanisms, from 0 to 1
  Omnivory <- 0.2 # mutiplyer; importance of auxiliary resources, "identity" for given, range from 0.1 to 0.9
  Predation <- 0.5 # mutiplyer; importance of predation for omnivores (% of total diet), range from 0.1 to 0.9
  SelfPredation <- 0.25 # mutiplyer; importance of self-predation for all nodes, range from 0 to 1
  SpatialMobility <- 0.25 # mutiplyer; importance of auxiliary microhabitats, range from 0.1 to 0.9


  # Applying coefficients
  # correct the PPMR coefficients accoording to the Reconstruction coefficients (see above)
  Data_traits$PPMRoptCorrected <- Data_traits$PPMRopt*PPMRopt
  Data_traits$PPMRwidthCorrected <- PPMRwidth

  # calculate PPMR distribution estimates for the optimum prey for each node
  Data_traits$PreyMean <- Data_traits$BodyMass_g_log10-Data_traits$PPMRoptCorrected
  Data_traits$PreySD <- Data_traits$MassSD_log10*Data_traits$PPMRwidthCorrected
  Data_traits$PreySD[Data_traits$PreySD<1] <- PPMRwidth_min # to avoid very narrow niche

  # replace 0.5 in spatial distribution (i.e. occasional occurrence) by the spatial coefficient
  Data_traits$above[Data_traits$above==0.5] <- SpatialMobility
  Data_traits$epi[Data_traits$epi==0.5] <- SpatialMobility
  Data_traits$hemi[Data_traits$hemi==0.5] <- SpatialMobility
  Data_traits$eu[Data_traits$eu==0.5] <- SpatialMobility

  ## Reconstruction loop

  # output table
  Foodweb_matrices <- list(NA)

  # loop
  for(l in 1:length(unique(Data_traits$transect_id))){

    # select a transect_id, add resources
    loop <- subset(Data_traits, transect_id==unique(Data_traits$transect_id)[l])

    if (length(which(is.na(loop$BodyMass_g)))>0){
      Foodweb_matrices[[l]] <-NA
    }
    if (length(which(is.na(loop$BodyMass_g)))==0) {
      loop <- bind_rows(loop, Data_Resources)

      # select groups
      consumers <- which(loop$BroadCategory!="Resource")
      resources <- which(loop$BroadCategory=="Resource")


      #FIRST: Produce resource preference matrix. For each node, feeding preferences are distributed
      # across resources according to pre-classified trophic guilds.

      # create an empty matrix
      matrix_Resources <- matrix(0,ncol=nrow(loop),nrow=nrow(loop))
      rownames(matrix_Resources) <- colnames(matrix_Resources) <- loop$Guild

      # set resource preferences as in the table
      matrix_Resources[resources,] <- t(loop[,c('P','L','S','F')])
      matrix_Resources[consumers,] <- t(loop[,rep("Invt",length(consumers))])

      # set ominvory coefficient
      matrix_Resources[matrix_Resources==0.5] <- Omnivory

      # fill missing
      matrix_Resources[is.na(matrix_Resources)==T] <- 0


      #SECOND: Produce size ditribution 'allometric' overlap matrix (to infer predator-prey interactions from PPMRs). Calculate overlaps of the size distributions correcting for PPMR for all potential predator-prey interactions. Overlap proportion [0;1] is used as the proxy of interaction strength. Optimum PPMRs were corrected for corresponding trait-based coefficients. Function to fill the matrix. Rows are preys, columns are predators.

      # create an empty matrix
      matrix_SizeOverlap <- matrix(0,ncol=nrow(loop),nrow=nrow(loop))
      rownames(matrix_SizeOverlap) <- colnames(matrix_SizeOverlap) <- loop$Guild

      # run a loop to calculate pairwise overlaps in size distributions (predator vs prey)
      for(r in 1:nrow(loop)){
        for(c in 1:nrow(loop)){
          overlap <- integrate(min.f1f2, -Inf, Inf, mu1=loop$BodyMass_g_log10[r], sd1=loop$MassSD_log[r],
                               mu2=loop$PreyMean[c], sd2=loop$PreySD[c])
          matrix_SizeOverlap[[r, c]] <- round(overlap$value,3)
        }
      }

      # fill the full overlap for resources feeding
      matrix_SizeOverlap[resources,] <- 1

      # Fill non-existing interactions with zero
      matrix_SizeOverlap[is.na(matrix_SizeOverlap)==T] <- 0


      #THIRD: Protection matrix

      # create an empty matrix
      matrix_Protection <- matrix(1,ncol=nrow(loop),nrow=nrow(loop))
      rownames(matrix_Protection) <- colnames(matrix_Protection) <- loop$Guild

      # run a loop to calculate isotopic distances (predator vs prey)
      # correct coefficients for protected fauna
      matrix_Protection <- matrix_Protection*loop$Protection # fill all rows with column protection
      matrix_Protection <- matrix_Protection*Protection # account for protection importance coefficient'

      # fill no protection for resources feeding
      matrix_Protection[resources,] <- 1


      #FOUR: Spatial overlaps

      # vertical distributions (Bray–Curtis dissimilarity)
      matrix_Spatial <- as.matrix(vegdist(loop[,c('above','epi','hemi','eu')],method="bray",na.rm = T))
      matrix_Spatial <- 1-matrix_Spatial
      rownames(matrix_Spatial) <- colnames(matrix_Spatial) <- loop$Guild

      # fill full overlap for resources feeding
      matrix_Spatial[resources,] <- 1

      # Fill non-existing interactions with zero
      matrix_Spatial[is.na(matrix_Spatial)==T] <- 0


      #FIFTH: We add biomass-based preferences directly here because we don't have biomasses of resources
      # and will have problems with fluxweb biomass preferences

      # fill the biomasses, scaled as proportions
      # We use square root biomass to avoid overrepresentation of dominant groups
      matrix_Biomass <- matrix(round(sqrt(loop$Biomass_fresh_g_sqm/sum(loop$Biomass_fresh_g_sqm,na.rm = T)), 5),
                               nrow=nrow(loop), ncol=nrow(loop), byrow = F)
      matrix_Biomass <- matrix_Biomass/colSums(matrix_Biomass, na.rm = T)[1] # Divide by coefficient to make sum = 1
      #colSums(matrix_Biomass, na.rm = T) # should be = 1 because of sqrt
      matrix_Biomass[resources,] <- 1 # fill resource biomasses to 1
      rownames(matrix_Biomass) <- colnames(matrix_Biomass) <- loop$Guild


      #FINALLY: The preference, allometric, protection, spatial and biomass matrices are multiplied to generate the final adjacency matrix that represents feeding preferences among food-web nodes.

      # merge all matrices except for resources
      if(weight_by_biomass){ # If true, include biomass in the weighted preference
        matrix_Adjacency <- matrix_SizeOverlap*matrix_Protection*matrix_Spatial*matrix_Biomass
      }else{ # If false, do not include biomass in the weighted preference
        matrix_Adjacency <- matrix_SizeOverlap*matrix_Protection*matrix_Spatial
      }
      # However, for predators the feeding preferences for prey should be distributed
      # Scale all fauna separately
      if (length(consumers)>1){
        matrix_Adjacency[consumers,] <-  scale(matrix_Adjacency[consumers,],
                                               center = FALSE,
                                               scale = colSums(matrix_Adjacency[consumers,],na.rm = T))
      } else {
        matrix_Adjacency[consumers,] <-  scale(matrix_Adjacency[consumers,],
                                               center = FALSE,
                                               scale = sum(matrix_Adjacency[consumers,],na.rm = T))
      }

      # Fill non-existing interactions with zero
      matrix_Adjacency[is.nan(matrix_Adjacency)==T] <- 0

      # Now scale predators to the proportion of fauna in the total diet and add resources separately
      # It is adjusting the feeding of ommnivores
      matrix_Adjacency[consumers,] <-
        t(t(matrix_Adjacency[consumers,])*matrix_Resources[consumers[1],]) # Scale ivrt
      matrix_Adjacency[resources,] <- matrix_Resources[resources,]

      # Check: The numbers here should be 0, 0.2 (omnivory coefficient), or 1
      # cbind(colSums(matrix_Adjacency[consumers,]))

      # Control the self-predation
      diag(matrix_Adjacency) <- diag(matrix_Adjacency)*SelfPredation # self-predation coefficient

      # Scale everything to 1
      matrix_Adjacency <-  scale(matrix_Adjacency,center = FALSE,scale = colSums(matrix_Adjacency))
      matrix_Adjacency[,resources] <- 0 # resources do not feed
      matrix_Adjacency[is.na(matrix_Adjacency)==T] <- 0

      # Should be all 1 or 0 for resources
      # colSums(matrix_Adjacency)

      # round the numbers
      matrix_Adjacency <- round(matrix_Adjacency,5)

      # data can be cleaned from improbable interactions, for example
      matrix_Adjacency[matrix_Adjacency<0.00001] <- 0


      # WRITE RESULTS

      Foodweb_matrices[[l]] <- matrix_Adjacency
    }


    print(paste(l,"completed"))
  }


  # Each matrix in the list is a transect_id, names are here:
  names(Foodweb_matrices) <- unique(Data_traits$transect_id)

  list(DataTraits=Data_traits,PreferenceMat=Foodweb_matrices,DataResources=Data_Resources)

}

foodweb.energyflux<-function(Data_traits,Data_Resources,Foodweb_matrices){

  ### Coefficients used for metabolism
  #Losses estimated by respiration /allometric relation for metabolism using Johnston and Sibly 2018
  #It should be noted that primary data from Johnston and Sibly is available if we want to go more precise than just macrofauna
  #For macrofauna (from Supp Table 1):
  # B0 = 20.90 ±0.47
  # a = 0.71 ±0.01
  # E = 0.64 ±0.01
  #with ln(B) = ln(B0) + a ln(M) – E/kT (1)
  #where B is standard metabolic rate (J hr-1), B0 is a normalisation constant, a represents the
  #allometric scaling exponent  M is individual bodymass (mg fresh mass),
  #E is the activation energy (eV), k is Boltzmann’s constant (8.62 × 10-5 eV K-1)
  #T is the experimental temperature (K)
  B0<- 20.90
  a<-0.71
  E<-0.64
  Kelvin<-273.15
  k<-8.62*10^(-5)

  #assimilation efficiencies following Lang et al. 2017:
  #herbivory link: e = 0.545 (for consumption "Resources: Living plant material")
  eh <- 0.545
  #carnivory link: e = 0.906
  ec<-0.906
  #detritivory link: e = 0.158 (for consumption "Resources: Leaf and root litter and Fungi")
  ed<-0.158
  #detritivory link for soil (addition not in Lang, for consumption "Resources: Soil organic matter")
  es<-0.05

  #conversion from J.hour^-1 to g C.year^-1
  #Citing Johnston and Sibly 2018: "Metabolism is then converted to respiration rate units (g C yr−1) by using
  #the conversion factors 37,490 J g C−1 = 20,100 J LO2−1 × (1/0.5363 g C LO2−1) and 8,760 h yr−1."
  Jper_gC<-37490
  hper_year<-8760

  ###Loop for calculating energy fluxes
  # output table
  Energyweb_matrices <- list(NA)

  #output food web properties
  outweb<-NULL

  # loop
  for(l in 1:length(unique(Data_traits$transect_id))){

    # select a transect_id, add resources
    loop <- subset(Data_traits, transect_id==unique(Data_traits$transect_id)[l])
    Temperature<-mean(loop$temp)
    loop <- bind_rows(loop, Data_Resources)

    # select groups
    consumers <- which(loop$BroadCategory!="Resource")
    resources <- which(loop$BroadCategory=="Resource")
    #correction of the paper of Johnston and Sibly: B0 instead of log(B0)
    Metabo_ln_perind<-a*log(loop$BodyMass_g[consumers]*10^3) + B0 - E/(k*(Temperature+Kelvin))
    Metabo_J_tot<-exp(Metabo_ln_perind)*loop$Density_sqm[consumers]

    matnumber<-which(names(Foodweb_matrices)==unique(Data_traits$transect_id)[l])
    matrix_Adjacency<-Foodweb_matrices[[matnumber]]
    eff<-rep(ec,length(colnames(matrix_Adjacency)))
    eff[which(colnames(matrix_Adjacency)=="Resources: Living plant material")]<-eh
    eff[which(colnames(matrix_Adjacency)=="Resources: Leaf and root litter")]<-ed
    eff[which(colnames(matrix_Adjacency)=="Resources: Soil organic matter")]<-es
    eff[which(colnames(matrix_Adjacency)=="Resources: Fungi")]<-ed

    Metabo_C_tot_peryear<-Metabo_J_tot/Jper_gC*hper_year

    #add resources to the loss vector
    loss<-c(Metabo_C_tot_peryear,rep(0,length(resources)))
    try.flux<-try(fluxing(matrix_Adjacency,losses=loss,efficiencies=eff,bioms.prefs = FALSE, bioms.losses=FALSE,ef.level="prey"),T)

    if (class(try.flux)[1]!="try-error"){
      fluxes<-try.flux
      #calculate now some metrics on the food web but a special function for this might be better
      TL<-GetTL2(fluxes)
      outweb<-rbind(c(transect_id=unique(Data_traits$transect_id)[l],
                      sumflux=sum(fluxes),maxTL=max(TL),Ngroups=length(consumers),sumOM=sum(fluxes["Resources: Soil organic matter",]),
                      sumplant=sum(fluxes["Resources: Living plant material",]),sumfungi=sum(fluxes["Resources: Fungi",]),
                      sumlitter=sum(fluxes["Resources: Leaf and root litter",])),outweb)
    }  else {
      fluxes<-NA
      outweb<-rbind(c(transect_id=unique(Data_traits$transect_id)[l],
                      sumflux=NA,maxTL=NA,Ngroups=length(consumers),sumOM=NA,
                      sumplant=NA,sumfungi=NA,sumlitter=NA),outweb)
    }
    # WRITE RESULTS

    Energyweb_matrices[[l]] <- fluxes
    print(paste(l,"completed"))
  }
  # Each matrix in the list is a transect_id, names are here:
  names(Energyweb_matrices) <- unique(Data_traits$transect_id)

  list(EnergyMat=Energyweb_matrices,FirstOutput=outweb)
}


foodweb.energyflux.plusdeath<-function(Data_traits,Data_Resources,Foodweb_matrices,deathrates){

  ### Coefficients used for metabolism
  #Losses estimated by respiration /allometric relation for metabolism using Johnston and Sibly 2018
  #It should be noted that primary data from Johnston and Sibly is available if we want to go more precise than just macrofauna
  #For macrofauna (from Supp Table 1):
  # B0 = 20.90 ±0.47
  # a = 0.71 ±0.01
  # E = 0.64 ±0.01
  #with ln(B) = ln(B0) + a ln(M) – E/kT (1)
  #where B is standard metabolic rate (J hr-1), B0 is a normalisation constant, a represents the
  #allometric scaling exponent  M is individual bodymass (mg fresh mass),
  #E is the activation energy (eV), k is Boltzmann’s constant (8.62 × 10-5 eV K-1)
  #T is the experimental temperature (K)
  B0<- 20.90
  a<-0.71
  E<-0.64
  Kelvin<-273.15
  k<-8.62*10^(-5)

  #assimilation efficiencies following Lang et al. 2017:
  #herbivory link: e = 0.545 (for consumption "Resources: Living plant material")
  eh <- 0.545
  #carnivory link: e = 0.906
  ec<-0.906
  #detritivory link: e = 0.158 (for consumption "Resources: Leaf and root litter and Fungi")
  ed<-0.158
  #detritivory link for soil (addition not in Lang, for consumption "Resources: Soil organic matter")
  es<-0.05

  #conversion from J.hour^-1 to g C.year^-1
  #Citing Johnston and Sibly 2018: "Metabolism is then converted to respiration rate units (g C yr−1) by using
  #the conversion factors 37,490 J g C−1 = 20,100 J LO2−1 × (1/0.5363 g C LO2−1) and 8,760 h yr−1."
  Jper_gC<-37490
  hper_year<-8760

  ###Loop for calculating energy fluxes
  # output table
  Energyweb_matrices <- list(NA)

  #output food web properties
  outweb<-NULL

  # loop
  for(l in 1:length(unique(Data_traits$transect_id))){

    # select a transect_id, add resources
    loop <- subset(Data_traits, transect_id==unique(Data_traits$transect_id)[l])
    Temperature<-mean(loop$temp)
    loop<-merge(loop,deathrates,by="Group")
    loop <- bind_rows(loop, Data_Resources)

    # select groups
    consumers <- which(loop$BroadCategory!="Resource")
    resources <- which(loop$BroadCategory=="Resource")
    #correction of the paper of Johnston and Sibly: B0 instead of log(B0)
    Metabo_ln_perind<-a*log(loop$BodyMass_g[consumers]*10^3) + B0 - E/(k*(Temperature+Kelvin))
    Metabo_J_tot<-exp(Metabo_ln_perind)*loop$Density_sqm[consumers]

    matnumber<-which(names(Foodweb_matrices)==unique(Data_traits$transect_id)[l])
    matrix_Adjacency<-Foodweb_matrices[[matnumber]]
    eff<-rep(ec,length(colnames(matrix_Adjacency)))
    eff[which(colnames(matrix_Adjacency)=="Resources: Living plant material")]<-eh
    eff[which(colnames(matrix_Adjacency)=="Resources: Leaf and root litter")]<-ed
    eff[which(colnames(matrix_Adjacency)=="Resources: Soil organic matter")]<-es
    eff[which(colnames(matrix_Adjacency)=="Resources: Fungi")]<-ed

    Metabo_C_tot_peryear<-Metabo_J_tot/Jper_gC*hper_year
    Death<-loop$death_rate[consumers]*loop$Biomass_fresh_g_sqm[consumers]*loop$Proportion.C[consumers]

    #add resources to the loss vector
    loss<-c(Metabo_C_tot_peryear+Death,rep(0,length(resources)))
    try.flux<-try(fluxing(matrix_Adjacency,losses=loss,efficiencies=eff,bioms.prefs = FALSE, bioms.losses=FALSE,ef.level="prey"),T)

    if (class(try.flux)[1]!="try-error"){
      fluxes<-try.flux
      #calculate now some metrics on the food web but a special function for this might be better
      TL<-GetTL2(fluxes)
      outweb<-rbind(c(transect_id=unique(Data_traits$transect_id)[l],
                      sumflux=sum(fluxes),maxTL=max(TL),Ngroups=length(consumers),sumOM=sum(fluxes["Resources: Soil organic matter",]),
                      sumplant=sum(fluxes["Resources: Living plant material",]),sumfungi=sum(fluxes["Resources: Fungi",]),
                      sumlitter=sum(fluxes["Resources: Leaf and root litter",])),outweb)
    }  else {
      fluxes<-NA
      outweb<-rbind(c(transect_id=unique(Data_traits$transect_id)[l],
                      sumflux=NA,maxTL=NA,Ngroups=length(consumers),sumOM=NA,
                      sumplant=NA,sumfungi=NA,sumlitter=NA),outweb)
    }
    # WRITE RESULTS

    Energyweb_matrices[[l]] <- fluxes
    print(paste(l,"completed"))
  }
  # Each matrix in the list is a transect_id, names are here:
  names(Energyweb_matrices) <- unique(Data_traits$transect_id)

  list(EnergyMat=Energyweb_matrices,FirstOutput=outweb)
}
