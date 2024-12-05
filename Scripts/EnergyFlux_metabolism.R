
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
    loop <- rbind.fill(loop, Data_Resources)
    
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