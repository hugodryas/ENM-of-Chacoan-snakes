# Adapted by Mario R. Moura
# Overprediction correction
MSDM_adapted <- function(records,    # Dataframe with presence records
                 SpNames,    # Vector with species names
                 RasterList, # List with rasters
                 Dispersion_Dist = NULL # Maximum dispersion distance
)
{
  
  # records = OccPresences
  # SpNames = names(ssp585_AvgGCM)
  # RasterList = ssp585_AvgGCM
  SpData <- records
  
  #Create a List
  Result <- list()
  dist.t <- list()
  
  # loop to process each species
  for (s in 1:length(SpNames)) {
    
    # Read the raster of the species:
    Adeq <- RasterList[[SpNames[s]]] 
    
    # Set output as NA if no presence is found among cells of the binary raster:
    if (length(Adeq[Adeq@data@values == TRUE]) ==0) {
      
      Result[[s]] <- Adeq
      dist.t[[s]] <- NA
    
    # If at least one presence exists among the binary raster cells, proceed:    
    }else{
      
      # Select the occurrence records of the target species:
      singleSpData <- SpData[SpData$sp == SpNames[s], ]
      
      # Transform occurrence records in a SpatialPoints object:
      pts1 <- singleSpData[ , c("x", "y")]
      coordinates(pts1) <- ~ x + y
      crs(pts1) <- crs(Adeq)
      
      # Change raster cell values from '0' to 'NA'. 
      AdeqBin <- Adeq 
      AdeqBin[AdeqBin[] == 0] <- NA
      
      # Detect contiguous patches formed by raster cells:
      AdeqBin <- raster::clump(AdeqBin)
      
      # Get patches as a data.frame object (x, y, clumps): 
      AdeqPoints_latlong <- data.frame(raster::rasterToPoints(AdeqBin))
      
      # Convert decimal degrees to geographical distance:
      AdeqPoints<-data.frame(
        (SoDA::geoXY(latitude=AdeqPoints_latlong$y, longitude=AdeqPoints_latlong$x, unit = 1000)),
        clumps=AdeqPoints_latlong$clumps)
      
      # Find the patches that contain presences records:
      polypoint <- unique(stats::na.omit(as.numeric(raster::extract(AdeqBin, pts1))))
      
      # Set to zero those raster cells within patches of suitable area (clumps) without species presence: 
      AdeqBin2 <- AdeqBin
      AdeqBin2[!AdeqBin2[] %in% polypoint] <- NA
      
      # Filter raster to include only patches with species presence:
      AdeqBin3 <- !is.na(AdeqBin2)
      
      # Filter patches that include species presence (occupied patches):
      CoordPathP <- as.data.frame(AdeqPoints[AdeqPoints[, 3] %in% polypoint,])
      
      # Filter patches that do not include species presence (patches unoccupied):
      CoordPathNP <- as.data.frame(AdeqPoints[!AdeqPoints[, 3] %in% polypoint,])
      
      # Check if there patches are either 'all occupied' or 'all unoccupied':
      if (nrow(CoordPathNP)==0 || nrow(CoordPathP)==0) {  # || stand for 'or'
        Mask2 <- Adeq
        Result[[s]] <- Mask2
        dist.t[[s]] <- NA
        
      # If number of 'occupied' and/or 'unoccupied' patches differ from zero, proceed:   
      }else{
        
        occupied_patches <- unique(CoordPathP[, 3])
        unoccupied_patches <- unique(CoordPathNP[, 3])
        
        # Create all possible combination between values of occupied_patches and unoccupied_patches:
        DistBetweenPoly0 <- expand.grid(occupied_patches, unoccupied_patches)
        
        # Create a empty column to register the edge-edge distance between patches:
        DistBetweenPoly0$Distance <- NA
        DistBetweenPoly0 <- as.matrix(DistBetweenPoly0)
        
        # Find Euclidean distance between patches with and without presences:
        for (i in 1:nrow(DistBetweenPoly0)) {
          
          # Select one combination of occupied_patches and unoccupied_patches:
          comb <- (DistBetweenPoly0[i, 1:2])
          
          # Get rasterPoints the selected combination of patches:
          A <- CoordPathP[CoordPathP[, 3] == comb[1], 1:2]
          B <- CoordPathNP[CoordPathNP[, 3] == comb[2], 1:2]
          
          # Get the nearest neighbour distance between cells of occupied and unoccupied selected patches:
          eucdist <- min(flexclust::dist2(A, B, method = 'euclidean', p = 2))
          
          # Store the minimum distance between pixels of occupied and unoccupied selected patches:
          DistBetweenPoly0[i, 3] <- eucdist
        }
        
        # If there is only one patch:
        if (nrow(DistBetweenPoly0) == 1){DistBetweenPoly0 <- DistBetweenPoly0}
          
        # If there are multiple patches: 
        else{
          
        # Order data.frame according to the unoccupied patch ID:
        DistBetweenPoly0 <- DistBetweenPoly0[order(DistBetweenPoly0[, 2]),]
          
        }
        
        # If the dispersion distance is not informed in the function argument:
        # OBR method------
        if(is.null(Dispersion_Dist)){
         
          # Find the maximum nearest neighbour distance between occurrence records:
          spraster <- rasterize(pts1, Adeq, field = 1)
          sps_latlong <- as.data.frame(as(spraster, 'SpatialPixels')@coords)
          sps <-SoDA::geoXY(latitude=sps_latlong$y, longitude=sps_latlong$x, unit = 1000)
          dist <- flexclust::dist2(sps, sps, method = 'euclidean', p = 2)
          dist[dist == 0] <- NA
          distmin <- apply(dist, 
                           MARGIN=1, # for each row
                           function(x) min(x, na.rm = TRUE))#
          CUT <- max(distmin) # km
          
          # Store the maximum nearest neighbour distance as the distance threshold:
          dist.t[[s]] <- data.frame(distdisp= CUT) 
          
          # Which unoccupied patches are at a distance <= CUT threshold from an occupied patch:
          Mask <- DistBetweenPoly0[DistBetweenPoly0[, 3] <= CUT, 2]
          
          # Get the ID of occupied patches and reachable unoccupied patches:
          FinalPatches<-c(Mask, occupied_patches)
          
          # Change the raster clump ID to zero if it does not match FinalPatches:
          Mask <- raster::match(AdeqBin, table=FinalPatches, nomatch=0)
          Mask <- Mask != 0 # set unmatched pixels to zero (this will include all raster extent)
          Mask[is.na(Adeq)] <- NA # set pixels outside the study area to NA
          Mask2 <- Adeq * Mask # mask unreachable patches
          
          # Save results as raster object
          Result[[s]] <- Mask2
          
        }
        
        # User-informed dispersal distance threshold:
        else{
          
          CUT <- as.numeric(Dispersion_Dist[[s]])
          
          # Which unoccupied patches are at a distance <= CUT threshold from an occupied patch:
          Mask <- DistBetweenPoly0[DistBetweenPoly0[, 3] <= CUT, 2]
          
          # Get the ID of occupied patches and reachable unoccupied patches:
          FinalPatches<-c(Mask, occupied_patches)
          
          # Change the raster clump ID to zero if it does not match FinalPatches:
          Mask <-raster::match(AdeqBin, table = FinalPatches, nomatch = 0)
          Mask <- Mask != 0 # set unmatched pixels to zero (this will include all raster extent)
          Mask[is.na(Adeq)] <- NA # set pixels outside the study area to NA
          Mask2 <- Adeq * Mask # mask unreachable patches
          
          # Save results as raster object
          Result[[s]] <- Mask2
        }
      } 
    }
  }
  
  
  if(is.null(Dispersion_Dist)){
    
    names(Result) <- SpNames
    names(dist.t) <- SpNames
    
    out <- list(MSDM = Result, Dispersion_t = dist.t)
    
  }else{
    
    names(Result) <- SpNames
    out <- c(Result)
  }
  
  return(out)
  
}

# Developed by José-Hidasi Neto
# https://rfunctions.blogspot.com/2016/10/calculating-temporal-beta-diversity-on.html
tempbetagrid<-function(oc1, oc2, index="jaccard", phylotree=phylo, phylobeta=F){
  tempturn<-numeric(nrow(oc1))
  tempnest<-numeric (nrow(oc1))
  tempbeta<-numeric(nrow(oc1))
  tempturnbeta<-numeric (nrow(oc1))
  for(i in 1:nrow(oc1) ){
    namesoc1<-names(oc1)[oc1[i,]==1]
    namesoc2<-names(oc2)[oc2[i,]==1]
    both<-namesoc1[namesoc1%in%namesoc2]
    bothmat<-rbind(rep(1,length(both)),rep(1,length(both)))
    colnames(bothmat)<-both
    namoc1<-namesoc1[namesoc1%in%namesoc2==FALSE]
    nam1mat<-rbind(rep(1,length(namoc1)),rep(0,length(namoc1)))
    colnames(nam1mat)<-namoc1
    namoc2<-namesoc2[namesoc2%in%namesoc1==FALSE]
    nam2mat<-rbind(rep(0,length(namoc2)),rep(1,length(namoc2)))
    colnames(nam2mat)<-namoc2
    matcomp<-cbind(bothmat,nam1mat,nam2mat)
    forprune<-t(data.frame(names(data.frame(matcomp))))
    colnames(forprune)<-forprune
    ifelse(phylobeta==T, betas<-phylo.beta.pair(matcomp, prune.sample(forprune, phylotree), index.family=index), betas<-beta.pair(matcomp, index.family=index) )
    tempturn[i]<-betas[[1]]
    tempnest[i]<-betas[[2]]
    tempbeta[i]<-betas[[3]]
    tempturnbeta[i]<-betas[[1]]/betas[[3]]}
  return(data.frame(turnover=tempturn,nestedness=tempnest,beta=tempbeta,turn.beta=tempturnbeta))}

