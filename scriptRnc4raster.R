#Télécharger données température

##########################################################################################################################################
library(ncdf4) #pour lire les nc4
library(raster) #pour lire les rasters
library(stringr) #
library(rgdal) #pour faire manipulations geospatialisés
library(gdalUtils) #pour faire manipulations geospatialisés
library(dplyr) 

wd <- "C:/IMERG/lst/"
setwd(wd)

files.hdf <- list.files(getwd(), pattern = "*.hdf$")

world <- readOGR("C:/IMERG/shp/world.shp")
poland <- subset(world, ADMIN == "Poland")

wgs84 = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"


for (i in 1:length(files.hdf)) {
  
  sds <- get_subdatasets(files.hdf[i]) #je lis le HDF avec toutes les variables
  #print sds
  name <- substr(files.hdf[i],1,nchar(files.hdf[i])-18) #me donne le nom "MYD21A1D.A2021001.h18v03.006"  = 001 = 01 janvier, si j'étais au 31 décembre je serais à 365. C'est cette information qui me permet d'identifier tuiles différentes mais du mm jours. 
  
  
  gdal_translate(sds[1], dst_dataset = paste0(name,".tif")) #prend uniquement la variable qui m'intéresse avec sds 1 . gls translate convertif HDF en tiff, et je lui dis de convertir uniquement la premiere bande du sds qui correspond à la variable que je veux 
  # (toutes les variables m'intéresse pas: view angle etc.. moi je m'interesse lst 1km)
  r <- raster(paste0(name,".tif"))
  r <- projectRaster(r, crs = wgs84, method = "ngb")
  writeRaster(r, paste0(name,".tif"), format = 'GTiff', overwrite = TRUE) #je suppose que ça me permet de convertir les HDF en tiff
  
  
}
#je vais me retrouver avec une tuile un tiff

# Je télécharge tous mes tiffs ? J'en ai 20, donc j'ai 5 jours. 4  tuiles par jour. On récupère information de la date, pour pouvoir aggréger les tuiles entres elles.
files.tif <- as.data.frame(list.files(getwd(), pattern = "*.tif$"))  #je créer un tableau avec identifiants des différents tiff
colnames(files.tif) <- "data"   #je renomme tableau
files.tif$day <- substr(files.tif$data,11,nchar(files.tif$data)-15) #j'extrait info  de la date et je créer  une colonne avec 
files.tif$day <- as.factor(files.tif$day) #transforme en facteur plutot quen valeur numerique
files.tif <- files.tif %>% mutate(select = as.integer((files.tif$day))) #je donne numero qui va me permettre d'iterer, 1 numéro pour un jour donnee

output <- NULL

for(i in unique(files.tif$select)) {
  
  select <- subset(files.tif, files.tif$select == i) #c'est ici que j'itère, va me permettre de subset tous mes files: selectionne tiff du mm jour (4 tuiles)
  
  rast.list <- list() #creer liste vide 
  j <- 1 
  
  #Boucle qui suit permet d'effectuer une fonction de raster sur chacun de mes identifiants.
  while (j < nrow(select)) {           #Intérêt while: si je fais une boucle for fonction raster va me refaire chaque fois le 1,2, 3 et le 4. Donc quand bcp données: ralentit. En gros permet de rasterisé une seule fois pour j=1 genre au début de la boucle puis après plus besoin car c deja fait
    #tant que j sera inférieur au nrow select il reprendra celui d'après. Une fois que c'est rasterisé, pas besoin de le faire à nouveau.
    r <- raster(select$data[1])  #function raster = methode pour construire un rasterlayer object
    #select.. nom du tiff qui nous intéresse, en gros permet de lire le premier raster typiquement la, donc une des tuiles des 4 tuiles du jour donnée à chaque fois. Il semble qu'il y ai une différence enrte raster et raster layer
    j <- j+1
    r2 <- resample(raster(select$data[j]),r,method='bilinear') #methode interpolation bilinéaire, on pourrait également utilisé méthode du nearest neighbor avec method='ngb'
    #Resample : Changes the spatial resolution of a raster dataset and sets rules for aggregating or interpolating values across the new pixel sizes. Use projectRaster if the target has a different coordinate reference system (projection).
    # resample select 2 tuiles du même jours et les resamples.
    #quand je veux clipper les tuiles ensembles: elles ont pas les mm formes, elles sont pas au même latitude , pas le mm extend. Pixel ne se superposent pas comme il faudrait. La sphéricité de la terre fait que en fonction de l'endroit ou on est: distance métrique, calcul déforme la réalité. 
    #pour parrer à ça: je resample = fait en sorte que tout les pixels de tous les raster est la mm taille ou tête? Il y aura légère deformation mais à l'échelle ou on travaille c'est pas grave.
    rast.list[[j]] <- r2
    
  }
  #rast.list correspond a une liste qui contient l'ensemble des rasters, qui ont tous le mm extend
  #Je vais enfin pouvoir les concaténer
  rast.list[[1]] <- r  # je sais pas pq mais le premier element de la liste ne contient rien avant de faire cela. Dois être une subtilité dans la boucle.
  
  rast.list$fun <- mean #jsp
  rast.mosaic <- do.call(mosaic,rast.list) #mosaic function qui permet de les concaténer. Do call permet exectuter sur ensemble des tiff dans rast.list
  
  r <- crop(rast.mosaic, extent(poland)) #découpe pour avoir shp que je veux 
  r <- mask(r, poland)
  
  #plot(r) : extend de la pologne , si i=1 pour premier janvier. On se rend compte qu'on a pas beaucoup de données pour certaine jour à cause des nuages. 
  #Conseil : quand je prendrais mes températures moyenne, prendre en considérations le nombre de pixel par exemple faire un ratio pour avoir nombre pixel ou y'a pas de données. 
  #faudra peut etre transformer les pixels ou y'a pas d'info (NA) en -99999 par exemple , et faire count pour compter le nombre de -99999 . Si je laisse en NA il pourra pas me compter 
  name <- substr(select$data[1],1,nchar(select$data[1])-15)
  
  date <- substr(name,11,nchar(name))
  date <- paste0(substr(date, 1, 5-1), "/", substr(date, 5, nchar(date)))
  
  #Stats#
  mean <- extract(r, poland, fun = mean, na.rm = TRUE, def = TRUE, weights = TRUE) # !!!! degrees kelvin !!!!
  med <- extract(r, poland, fun = median, def = TRUE, na.rm = TRUE)
  sd <- extract(r, poland, fun = sd, def = TRUE, na.rm = TRUE)
  max <- extract(r, poland, fun = max, def = TRUE, na.rm = TRUE)
  min <- extract(r, poland, fun = min, def = TRUE, na.rm = TRUE)
  
  results <- as.data.frame(cbind(date, mean, med, sd, max, min))
  colnames(results) <- c("date", "mean", "med", "sd", "max", "min")
  
  output <- rbind(output, results)
  
  writeRaster(r, paste0(wd,"/mosaic/",name,"_POLAND.tif"), format = 'GTiff', overwrite = TRUE)
  
  rm(mean)
  rm(med)
  rm(sd)
  rm(max)
  rm(min)
  
}

write.csv(output, file = "C:/IMERG/results/lst_stats_POLAND.csv", row.names=F)