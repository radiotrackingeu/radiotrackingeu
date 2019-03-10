######### --------- Global Variables --------- #########
globalVariables(c("j","i","timestamp","freq_tag","Station","receiver","Name","input"))
######### --------- Data --------- #########


#' Data set containing signal recordings of one bat on two stations.
#'
#' @author Dipl.-Phys. Ralf Zeidler, \email{ralf.zeidler@fridata.de}
#'
#' @name logger_data
#' @docType data
#' @author Ralf Zeidler \email{ralf.zeidler@fridata.de}
#' @references \url{www.radio-tracking.eu}
NULL

#' Data set containing signal recordings of one bat on two stations.
#'
#' It is a subst of logger_data to speed up process time while running the examples.
#'
#' @author Dipl.-Phys. Ralf Zeidler, \email{ralf.zeidler@fridata.de}
#'
#' @name ld_small
#' @docType data
#' @author Ralf Zeidler \email{ralf.zeidler@fridata.de}
#' @references \url{www.radio-tracking.eu}
NULL

#' Data set containing position and orientation of two tracking stations.
#'
#' @author Dipl.-Phys. Ralf Zeidler, \email{ralf.zeidler@fridata.de}
#'
#' @name receiver_data
#' @docType data
#' @author Ralf Zeidler \email{ralf.zeidler@fridata.de}
#' @references \url{www.radio-tracking.eu}
NULL

######### --------- Direction of arrival functions --------- #########

#' Calculates the smaller angle between two vectors.
#'
#' @author Dipl.-Phys. Ralf Zeidler, \email{ralf.zeidler@fridata.de}
#'
#' @param angle_a Angle of vector a.
#' @param angle_b Angle of vector b.
#' @return The smaller angle between both vectors - always less than or equal to 180 degrees.
#' @examples
#' angle_between(30,77)
#' #returns 47
#' angle_between(20,350)
#' #returns 30
#' @export
#'
angle_between <- function(angle_a,angle_b){
  ((((angle_b - angle_a) %% 360) + 540) %% 360) - 180
}

#' Calculates the DoA using the signal strengths at two antennas.
#'
#' @author Dipl.-Phys. Ralf Zeidler, \email{ralf.zeidler@fridata.de}
#'
#' @param sig_a Signal strength at antenna a.
#' @param sig_b Signal strength at antenna b.
#' @param angle_a Orientation of antenna a.
#' @param angle_b Orientation of antenna b.
#' @param dbLoss Signal drop at 90 degrees in dB.
#' @param option A string specifying the approximation method. "arccos" can only be used with antennas at 90 degree intervals, "linear" is more rough and "automatic" chooses depending on the angle between antennas.
#' @return Approximation of the DoA.
#' @examples
#' calc_angle(20,30,0,90,15,"arccos")
#' calc_angle(20,20,0,90,15,"arccos")
#' @export
#'
calc_angle <- function(sig_a, sig_b, angle_a, angle_b, dbLoss, option){
  #options: linear, arcos, lookup, automatic, old_linear
  #lookup still to be implemented
  if(option=="automatic"){
    if(abs(angle_between(angle_a,angle_b))==90){
      option<-"arccos"
    }else{
      option<-"linear"
    }
  }
  alpha<-angle_between(angle_a,angle_b)
  #if the left antenna is the left one
  if(alpha>0){
    sig_l<-sig_a
    sig_r<-sig_b
    angle_l<-angle_a
    angle_r<-angle_b
  }
  #if the left antenna is the right one
  if(alpha<0){
    sig_l<-sig_b
    sig_r<-sig_a
    angle_l<-angle_b
    angle_r<-angle_a
    alpha<-angle_between(angle_l,angle_r)
  }
  #if the the angle ist the same
  if(alpha==0) return(NA)
  delta_m<-(sig_l-sig_r)/dbLoss
  if(abs(delta_m)>1){
    return(NA)
  }
  switch(option,
         linear = {
           angle <- 1/2*(alpha-alpha*delta_m/(sin(pi*alpha/180)^2))+angle_l
         },
         arccos = {
           angle <- acos(delta_m)*90/pi+angle_l
         }
  )
  if(angle<0){
    angle<-angle+360
  }
  if(angle>360){
    angle<-angle-360
  }
  return(angle)
}

#' Matches timestamps by dividing into time slots.
#'
#' @author Dipl.-Phys. Ralf Zeidler, \email{ralf.zeidler@fridata.de}
#'
#' @param data Data frame with certain column names containing the signal logs.
#' @param station_time_error Maximum length of the time slots.
#' @param progress Whether this call is wrapped in a \code{shiny::withProgress} call.
#' @return Data frame with matched timestamps.
#' @examples
#' time_match_signals(ld_small)
#' @export
#'
time_match_signals <- function(data,station_time_error=0.3, progress=FALSE){
  matched_data<-NULL
  cnt_stats=0
  #for each station
  for(i in unique(data$Name)){
    if (progress) {
      setProgress(value=cnt_stats)
      incProgress(amount=0, detail = paste0("Station: ",i))
    }
    tmp_s<-subset(data,Name==i)
    num_tags=length(unique(tmp_s$freq_tag))
    #for each frequency tag
    for(l in unique(tmp_s$freq_tag)){
      tmp_sf<-subset(tmp_s,freq_tag==l)
      tmp_sf<-tmp_sf[order(tmp_sf$timestamp),]
      #calculate time difference between the loggings
      tmp_sf$td<-c(0,diff(tmp_sf$timestamp))
      tmp_s$ti <- as.POSIXct(NA, origin = "1970-01-01", tz="UTC")
      gc<-0
      tmp_sf$ti[1]<-tmp_sf$timestamp[1]
      for(i in 2:nrow(tmp_sf)){
        if(sum(tmp_sf$td[(i-gc):i])<=station_time_error){
          tmp_sf$ti[i]<- tmp_sf$timestamp[i-gc-1]
          if(any(duplicated(tmp_sf$receiver[(i-gc-1):i]))){
            tmp_sf$ti[i]<- tmp_sf$timestamp[i]
            gc<--1
          }
          gc<-gc+1
        }else{
          tmp_sf$ti[i]<- tmp_sf$timestamp[i]
          gc<-0
        }
      }
      matched_data<-rbind(matched_data,tmp_sf)
    }
    cnt_stats<-cnt_stats+1
  }

  matched_data$timestamp<-as.POSIXct(matched_data$ti, origin = "1970-01-01", tz ="UTC")
  return(matched_data)
}

#' Matches timestamps by fitting a smoothing spline.
#'
#' @author Dipl.-Phys. Ralf Zeidler, \email{ralf.zeidler@fridata.de}
#'
#' @param data Data frame with certain column names containing the signal logs.
#' @param spar_value Smoothing coefficient to be used with stats::smooth.spline.
#' @param progress Whether this call is wrapped in a \code{shiny::withProgress} call.
#' @return Data frame with matched timestamps.
#' @examples
#' smooth_to_time_match(ld_small,spar_value=0.01, progress=FALSE)
#' @export
#'
smooth_to_time_match <-function(data,spar_value=0.01, progress=FALSE){
  smoothed_data<-NULL
  cnt_recs=0
  #for each receiver
  for(i in unique(data$receiver)){
    if (progress) {
      setProgress(value=cnt_recs)
      incProgress(amount=0, detail = paste0("Receiver: ",i))
    }
    tmp_r<-subset(data,receiver==i)
    num_tags=length(unique(tmp_r$freq_tag))
    #for each frequency tag
    for(l in unique(tmp_r$freq_tag)){
      tmp_rf<-subset(tmp_r,freq_tag==l)
      if (nrow(tmp_rf)<5) {
        print(paste0('skipping freq "',l,'" on receiver "',i,'": not enough signals (',nrow(tmp_rf),')'))
        next
      }
      time_seq<-unique(c(round(tmp_rf$timestamp))) #,round(tmp_rf$timestamp)+1,round(tmp_rf$timestamp)-1)
      smoothed<-data.frame(max_signal=predict(
        smooth.spline(tmp_rf$timestamp,tmp_rf$max_signal,spar=spar_value),
        as.numeric(time_seq))$y,
        timestamp=time_seq,
        receiver=i,
        freq_tag=l,stringsAsFactors = FALSE)
      smoothed_data<-rbind(smoothed_data,smoothed)
      if(progress)
        incProgress(amount=1/num_tags)
    }
    cnt_recs<-cnt_recs+1
  }
  return(smoothed_data)
}

#' Wrapper function for doa_internal - requires time matched signals
#'
#' @author Dipl.-Phys. Ralf Zeidler, \email{ralf.zeidler@fridata.de}
#'
#' @export
#' @importFrom shiny withProgress
#'
#' @param signals Data frame with certain column names containing the signal logs.
#' @param receivers Data frame with certain column names containing receiver / antenna information (location, orientation).
#' @param dBLoss Signal drop at 90 degrees in dB - passed to calc_angle().
#' @param live_mode Information required by the rteu Logger shiny app.
#' @param live_update_interval Information required by the rteu Logger shiny app.
#' @param progress Whether this call is wrapped in a \code{shiny::withProgress} call.
#' @return Data frame containing the DoAs.
#' @examples
#' logger_data_tm <- time_match_signals(ld_small)
#' doa(logger_data_tm,receiver_data)
#'
doa <- function(signals, receivers,dBLoss=14, live_mode=FALSE, live_update_interval=15, progress = FALSE){
  data<-merge(signals,receivers,by.x="receiver",by.y="Name")
  # time_to_look_for<-NULL
  #for each timestamp of the smoothed data
  progress = FALSE
  if(!live_mode){
    time_to_look_for<-unique(data$timestamp)
  }else{
    if(length(unique(data$timestamp))<live_update_interval){
      end_point<-length(unique(data$timestamp))-1
    }else{
      end_point<-live_update_interval-1
    }
    time_to_look_for<-unique(data$timestamp)[order(unique(data$timestamp),decreasing = TRUE)][1:end_point]
  }
  if(progress)
    withProgress(min=0, max=length(time_to_look_for), value=0, expr={
      doa_internal(data, time_to_look_for, progress,dBLoss=input$dBLoss,doa_approx=input$doa_option_approximation)
    })
  else
    doa_internal(data, time_to_look_for,dBLoss=dBLoss, doa_approx="automatic", progress=FALSE)
}

#' Chooses which antennas to use to calculate the DoA.
#'
#' @author Dipl.-Phys. Ralf Zeidler, \email{ralf.zeidler@fridata.de}
#'
#' @param data Data frame with certain column names containing signal logs and receiver data.
#' @param time_to_look_for Vector of timestamps to be analyzed.
#' @param dBLoss Signal drop at 90 degrees in dB - passed to calc_angle().
#' @param doa_approx A string specifying DoA approximation method used in calc_angle(). One of "arccos", "linear", "automatic".
#' @param progress Whether this call is wrapped in a \code{shiny::withProgress} call.
#' @return dataframe containing the DoAs.
#'
doa_internal <- function(data, time_to_look_for, dBLoss=14, doa_approx="automatic", progress=FALSE) {
  cnt_timestamp=0
  split<-foreach(t=time_to_look_for,
                 .export=c("angle_between","calc_angle"),
                 .combine=rbind,
                 .inorder=FALSE) %dopar% {
                   if (progress)
                     setProgress(value=cnt_timestamp, message = "Computing Bearings... ")
                   #build subset for the timestamp
                   data_t<-subset(data,timestamp==t)
                   num_tags = length(data_t$freq_tag)
                   for(f in unique(data_t$freq_tag)) {
                     #build subset for the frequency
                     data_tf<-subset(data_t,freq_tag==f)
                     for(s in unique(data_tf$Station)) {
                       result<-NULL
                       output<-NULL
                       #build subset for the Station
                       data_tfs<-subset(data_tf, Station==s)
                       #sort using signal_strength
                       data_tfs<-unique(data_tfs[order(data_tfs$max_signal, decreasing = TRUE, na.last=NA),])
                       if(nrow(data_tfs)>1){
                         if(anyNA(data_tfs[1:2,]))
                           next
                         #check angle between strongest and second strongest and if it is smaller then 90 degree, calc it linearly
                         if(abs(angle_between(data_tfs[1,"Orientation"],data_tfs[2,"Orientation"]))<=120){
                           angle<-calc_angle(data_tfs[1,"max_signal"],data_tfs[2,"max_signal"],data_tfs[1,"Orientation"],data_tfs[2,"Orientation"],dBLoss,doa_approx)
                           result<-cbind.data.frame(timestamp=as.POSIXct(t,origin="1970-01-01",tz="UTC"),angle=angle,antennas=nrow(data_tfs),Station=s,freq_tag=f,strength=max(data_tfs$max_signal),stringsAsFactors=FALSE)
                         }else{
                           #back antenna plays a big role here
                           if(nrow(data_tfs)>2){
                             angle_1<-data_tfs[1,"Orientation"]
                             angle_2<-calc_angle(data_tfs[1,"max_signal"],data_tfs[3,"max_signal"],data_tfs[1,"Orientation"],data_tfs[3,"Orientation"],2*dBLoss,"linear")
                             angle<-angle_1+angle_between(angle_1,angle_2)/abs(angle_between(data_tfs[1,"Orientation"],data_tfs[2,"Orientation"]))*60
                             result<-cbind.data.frame(timestamp=as.POSIXct(t,origin="1970-01-01",tz="UTC"),angle=angle,antennas=nrow(data_tfs),Station=s,freq_tag=f,strength=max(data_tfs$max_signal),stringsAsFactors=FALSE)
                           }
                           if(nrow(data_tfs)==2){
                             angle<-data_tfs[1,"Orientation"]
                             result<-cbind.data.frame(timestamp=as.POSIXct(t,origin="1970-01-01",tz="UTC"),angle=angle,antennas=nrow(data_tfs),Station=s,freq_tag=f,strength=max(data_tfs$max_signal),stringsAsFactors=FALSE)
                           }
                         }
                       }
                       if(nrow(data_tfs)==1){
                         if(anyNA(data_tfs[1,]))
                           next
                         angle<-data_tfs[1,"Orientation"]
                         result<-cbind.data.frame(timestamp=as.POSIXct(t,origin="1970-01-01",tz="UTC"),angle=angle,antennas=nrow(data_tfs),Station=s,freq_tag=f,strength=max(data_tfs$max_signal),stringsAsFactors=FALSE)
                       }
                       cnt_timestamp<-cnt_timestamp+1
                       output<-rbind(output,result)
                     }
                   }
                   output
                 }
  return(split)
}

#' Converts longitude/latitude coordinates to UTM based WGS84 Datum.
#'
#' @author Dipl.-Phys. Ralf Zeidler, \email{ralf.zeidler@fridata.de}
#'
#' @export
#'
#' @importFrom sp coordinates
#' @importFrom sp proj4string
#' @importFrom sp CRS
#' @importFrom sp spTransform
#'
#' @param x Longitude
#' @param y Latitude
#' @return UTM coordinates in a data frame with columns X, Y and zone.
#'
wgstoutm<-function(x,y){
  tmp<-data.frame()
  for(i in 1:length(x)){
    zone<-(floor((x[i] + 180)/6) %% 60) + 1
    xy <- data.frame(cbind("X"=x[i],"Y"=y[i]))
    sp::coordinates(xy) <- c("X", "Y")
    sp::proj4string(xy) <- sp::CRS("+proj=longlat +datum=WGS84")  ## for example
    res <- sp::spTransform(xy, sp::CRS(paste("+proj=utm +zone=",zone," ellps=WGS84",sep='')))
    tmp<-rbind(tmp,cbind.data.frame(X=res$X,Y=res$Y,zone))
  }
  return(tmp)
}

#' Converts UTM coordinates to longitude/latitude based on WGS84 Datum.
#'
#' @author Dipl.-Phys. Ralf Zeidler, \email{ralf.zeidler@fridata.de}
#'
#' @importFrom sp coordinates
#' @importFrom sp proj4string
#' @importFrom sp CRS
#' @importFrom sp spTransform
#'
#' @param x Easting
#' @param y Northing
#' @param zone Zone
#' @return Data frame holding longitude (X) and latitude (Y) columns.
#' @export
#'
utmtowgs<-function(x,y,zone){
  tmp<-data.frame()
  for(i in 1:length(x)){
    xy <- data.frame(cbind("X"=x[i],"Y"=y[i]))
    sp::coordinates(xy) <- c("X", "Y")
    sp::proj4string(xy) <- sp::CRS(paste0("+proj=utm +zone=",zone[i]," +datum=WGS84"))  ## for example
    res <- sp::spTransform(xy, sp::CRS("+proj=longlat +datum=WGS84"))
    tmp<-rbind(tmp,as.data.frame(res))
  }
  return(tmp)
}

######### --------- Localization functions --------- #########

# progress: TRUE if function is wrapped in withProgress() call
#' Calculates the intersection of bearings of the same frequency and time interval.
#'
#' @author Dipl.-Phys. Ralf Zeidler, \email{ralf.zeidler@fridata.de}
#'
#' @export
#'
#' @importFrom sp coordinates
#' @importFrom utils combn
#'
#' @param receivers Data frame of all receivers with station names and positions.
#' @param bearings Data frame produced by doa().
#' @param only_one Whether a location should be estimated based on signal strength if only one bearing is available in a time slot.
#' @param time_error_inter_station Accepted time discrepancy between station clocks in seconds. Used in timematch_inter()
#' @param angles_allowed Vector of lower and upper limits for angles between bearings.
#' @param tri_option What to do if more than two bearings per time slot are available: "two_strongest" - just consider the two bearings with the strongest signal, "centroid" - calculate a centroid of all intersection points
#' @param tm_method Which method of time matching to use: "tm" - time slots, "spline" - spline smoothing See \code{\link{smooth_to_time_match}} and \code{\link{time_match_signals}}.
#' @param spar Smoothing factor for the spline to be used with \code{tm_method="spline"}.
#' @param progress Whether this call is wrapped in a \code{shiny::withProgress} call.
#' @return Data frame containing triangulated UTM coordinates.
#' @examples
#' logger_data_tm <- time_match_signals(ld_small)
#' #registerDoParallel(cores=2)
#' doa <- doa(logger_data_tm,receiver_data)
#' triangulate(receiver_data,doa)
#' #stopImplicitCluster()
#'
triangulate <- function(receivers, bearings, only_one=FALSE,time_error_inter_station=0.6,angles_allowed=c(30,150),tri_option="centroid",tm_method = "spline",spar=0.01, progress=FALSE) {
  progress=FALSE
  positions<-data.frame()
  #Calc UTM of Stations and add them
  stations<-na.omit(unique(receivers[,c("Station","Longitude","Latitude")]))
  stations<-stations[!duplicated(stations$Station),]
  stations_utm<-cbind(stations,utm=wgstoutm(stations[,"Longitude"],stations[,"Latitude"]))

  if(length(unique(stations_utm$utm.zone))>1){
    print("UTM Zone Problem")
  }
  #for each frequency tag
  freq_names=unique(bearings$freq_tag)
  num_freq_names=length(freq_names)
  cnt_freq_names=0
  result<-NULL
  for(i in freq_names){
    tmp_f <- subset(bearings,freq_tag==i)
    tmp_f <- switch(tm_method,
                    tm = timematch_inter(tmp_f,time_error_inter_station),
                    spline = smooth_to_time_match_bearings(tmp_f,receivers,spar))
    timestamps_unique<-unique(tmp_f$timestamp)
    num_timestamps_unique<-length(timestamps_unique)
    #for each times interval
    split<-foreach(j=timestamps_unique,
                   .export=c("tri_one","tri_two","tri_centroid","utmtowgs","coordinates","angle_between","triang"),
                   .packages=c("sp"),
                   .combine=rbind,
                   .inorder=FALSE) %dopar% {
                     if(progress)
                       setProgress(value=cnt_freq_names)
                     tmp_ft <- subset(tmp_f,timestamp==j)
                     tmp_fts <- merge(tmp_ft,stations_utm,by.x="Station",by.y="Station")
                     #calculate positions for two or more bearings in one slot
                     if(nrow(tmp_fts)==1&only_one){
                       positions<-cbind(timestamp=j,freq_tag=i,pos=tri_one(tmp_fts))
                     }
                     if(nrow(tmp_fts)>=2){
                       positions<-cbind(
                         timestamp=j,
                         freq_tag=i,
                         pos=switch(tri_option,
                                    centroid =  tri_centroid(tmp_fts,angles_allowed),
                                    two_strongest = tri_two(tmp_fts,angles_allowed)
                         )
                       )
                     }
                     positions
                   }
    result<-rbind(split,result)
    cnt_freq_names<-cnt_freq_names+1
  }
  if(nrow(result)>0){
    return(result[order(result$timestamp),])
  }
}

#' function to match times between two or more station
#'
#' @author Dipl.-Phys. Ralf Zeidler, \email{ralf.zeidler@fridata.de}
#'
#' @export
#' @importFrom stats complete.cases
#'
#' @param data bearing data
#' @param inter_error maximal error between stations
#' @return data frame with matched timestamps
#' @examples
#' #registerDoParallel(cores=2)
#' logger_data_tm <- time_match_signals(ld_small)
#' doa <- doa(logger_data_tm,receiver_data)
#' timematch_inter(doa)
#' #stopImplicitCluster()
#'
timematch_inter <- function(data,inter_error=0.6){
  if(nrow(data)==1) return(data)
  data<-data[complete.cases(data), ]
  tmp_s<-data[order(data$timestamp),]
  tmp_s$td <- c(0,diff(tmp_s$timestamp))
  tmp_s$ti <- NA
  gc<-0
  tmp_s$ti[1]<-tmp_s$timestamp[1]
  for(i in 2:nrow(tmp_s)){
    if(sum(tmp_s$td[(i-gc):i])<=inter_error){
      tmp_s$ti[i]<- tmp_s$timestamp[i-gc-1]
      if(any(duplicated(tmp_s$Station[(i-gc-1):i]))){
        tmp_s$ti[i]<- tmp_s$timestamp[i]
        gc<--1
      }
      gc<-gc+1
    }else{
      tmp_s$ti[i]<- tmp_s$timestamp[i]
      gc<-0
    }
  }
  tmp_s$timestamp<-as.POSIXct(tmp_s$ti,origin="1970-01-01")
  tmp_s$ti<-NULL
  tmp_s$td<-NULL
  return(tmp_s)
}

#' Does the actual triangulation. Assumes GK coordinates or UTM coordinates within the same zone.
#'
#' @author Dipl.-Phys. Ralf Zeidler, \email{ralf.zeidler@fridata.de}
#' @export
#'
#' @param x1 Easting of the first station.
#' @param y1 Northing of the first station.
#' @param alpha1 DoA of the first station.
#' @param x2 Easting of the second station.
#' @param y2 Northing of the second station.
#' @param alpha2 DoA of the second station.
#' @return Vector with Easting \code{x} and Northing \code{y} of triangulated position.
#' @examples
#' triang(0,0,45,10,0,315)
#'
triang <- function(x1,y1,alpha1,x2,y2,alpha2){
  # For Triangulation GK Coordinates are necessary!
  # First calculate tan keeping in mind that 0° in geo-coordinates are 90° in a x-y plane
  ta1 <- (alpha1%%360)/180*pi
  ta2 <- (alpha2%%360)/180*pi
  if(((alpha1-alpha2)%%180)==0){#print("No triangulation possible: all three points are on one line")
    return(c(NA,NA))}

  # Finding Intersection Using solver
  b<-c(x2-x1,y2-y1)
  a1<-c(sin(ta1),cos(ta1))
  a2<-c(-sin(ta2),-cos(ta2))
  a<-matrix(c(a1,a2),nrow=2)
  l<-solve(a,b)
  px<-x1+l[1]*sin(ta1)
  py<-y1+l[1]*cos(ta1)

  if(l[2]>0&l[1]>0)
  {
    return(c(px,py))
  }
  else{
    return(c(NA,NA))
  }
}

#' Calculates the centroid of intersection points of the given bearings and checks if intersection angles are in allowed range.
#'
#' @author Dipl.-Phys. Ralf Zeidler, \email{ralf.zeidler@fridata.de}
#'
#' @export
#'
#' @param tmp_fts Data frame of bearings. Expects columns \code{utm.X, utm.Y, utm.zone, angle}.
#' @param angles_allowed Vector of lower and upper limits for angles between bearings.
#' @return Data frame holding centroid position as lng/lat \code{X,Y} and UTM \code{utm.X, utm.Y}. Will hold \code{NA} if triangulation yields no result.
#'
tri_centroid <- function(tmp_fts,angles_allowed){
  tmp_positions<-data.frame()
  for(c in 1:dim(combn(nrow(tmp_fts),2))[2]){
    e<-combn(nrow(tmp_fts),2)[1,c]
    z<-combn(nrow(tmp_fts),2)[2,c]
    #order using singal strength and take first two
    #tmp_fts<-tmp_fts[order(tmp_fts$strength,decreasing = TRUE, na.last=NA),]
    if(anyNA(tmp_fts[c(e,z),]))
      next
    if(abs(angle_between(tmp_fts$angle[e],tmp_fts$angle[z]))<angles_allowed[1]|abs(angle_between(tmp_fts$angle[z],tmp_fts$angle[e]))>angles_allowed[2])
      next
    location<-triang(tmp_fts$utm.X[e],tmp_fts$utm.Y[e],tmp_fts$angle[e],tmp_fts$utm.X[z],tmp_fts$utm.Y[z],tmp_fts$angle[z])
    if(anyNA(location))
      next
    tmp_positions<-rbind(tmp_positions,cbind(utm.x=location[1],utm.y=location[2],utm.zone=tmp_fts$utm.zone[1]))
  }
  if(nrow(tmp_positions)>0){
    x<-mean(tmp_positions$utm.x,na.rm = T)
    y<-mean(tmp_positions$utm.y,na.rm = T)
    zone<-tmp_positions$utm.zone[1]
    location_wgs<-utmtowgs(x,y,zone)
    return(data.frame(location_wgs,utm.X=x,utm.Y=y))
  }else{
    return(data.frame(X=NA,Y=NA,utm.X=NA,utm.Y=NA))
  }
}

#' Approximates a position using one bearing and an estimation of the distance.
#'
#' @author Dipl.-Phys. Ralf Zeidler, \email{ralf.zeidler@fridata.de}
#'
#' @param tmp_fts Data frame of bearings. Expects columns \code{utm.X, utm.Y, utm.zone, angle, strength}.
#' @return Data frame holding approximated position as lng/lat \code{X,Y} and UTM \code{utm.X, utm.Y}.
#'
tri_one <- function(tmp_fts){
  str_mod<-50000
  location<-data.frame(utm.X=tmp_fts$utm.X+cospi((90-tmp_fts$angle)/180)/tmp_fts$strength*str_mod, utm.Y=tmp_fts$utm.Y+sinpi((90-tmp_fts$angle)/180)/tmp_fts$strength*str_mod)
  location_wgs<-utmtowgs(location$utm.X,location$utm.Y,tmp_fts$utm.zone)
  return(data.frame(location_wgs,location))
}

#' Calculates the position using only the bearings of two the strongest signals.
#'
#' @author Dipl.-Phys. Ralf Zeidler, \email{ralf.zeidler@fridata.de}
#'
#' @param tmp_fts Data frame of bearings. Expects columns \code{utm.X, utm.Y, utm.zone, angle}.
#' @param angles_allowed Vector of lower and upper limits for angles between bearings.
#' @return Data frame holding the calculated position as lng/lat \code{X,Y} and UTM \code{utm.X, utm.Y}. Will hold \code{NA} if triangulation yields no result.
#'
tri_two <- function(tmp_fts,angles_allowed){
  tmp_positions<-data.frame()
  #order using singal strength and take first two
  tmp_fts<-tmp_fts[order(tmp_fts$strength,decreasing = TRUE, na.last=NA),]
  if(anyNA(tmp_fts))
    return(data.frame(X=NA,Y=NA,utm.X=NA,utm.Y=NA))
  if(abs(angle_between(tmp_fts$angle[1],tmp_fts$angle[2]))<angles_allowed[1]|abs(angle_between(tmp_fts$angle[2],tmp_fts$angle[1]))>angles_allowed[2])
    return(data.frame(X=NA,Y=NA,utm.X=NA,utm.Y=NA))
  location<-triang(tmp_fts$utm.X[1],tmp_fts$utm.Y[1],tmp_fts$angle[1],tmp_fts$utm.X[2],tmp_fts$utm.Y[2],tmp_fts$angle[2])
  if(anyNA(location))
    return(data.frame(X=NA,Y=NA,utm.X=NA,utm.Y=NA))
  location_wgs<-utmtowgs(location[1],location[2],tmp_fts$utm.zone[1])
  if(nrow(location_wgs)>0){
    return(data.frame(location_wgs,utm.X=location[1],utm.Y=location[2]))
  }else{
    return(data.frame(X=NA,Y=NA,utm.X=NA,utm.Y=NA))
  }
}

#' uses a spline to time match the signal strength between several stations
#'
#' @author Dipl.-Phys. Ralf Zeidler, \email{ralf.zeidler@fridata.de}
#'
#' @export
#'
#' @importFrom shiny setProgress
#' @importFrom shiny incProgress
#' @importFrom stats na.omit
#' @importFrom stats predict
#' @importFrom stats smooth.spline
#'
#' @param data Data frame containing logger data.
#' @param receivers Data frame containing receiver information.
#' @param spar_value Smoothing factor of the spline.
#' @param progress Whether this call is wrapped in a \code{shiny::withProgress} call.
#'
#' @return return dataframe containg time matched signal strength using spline
#'
smooth_to_time_match_bearings <-function(data,receivers,spar_value=0.01, progress=FALSE){
  smoothed_data<-NULL
  cnt_recs=0
  #for each receiver
  for(i in unique(data$Station)){
    if (progress) {
      setProgress(value=cnt_recs)
      incProgress(amount=0, detail = paste0("Station: ",i))
    }
    tmp_r<-subset(data,Station==i)
    num_tags=length(unique(tmp_r$freq_tag))
    #for each frequency tag
    for(l in unique(tmp_r$freq_tag)){
      tmp_rf<-na.omit(subset(tmp_r,freq_tag==l))
      if (nrow(tmp_rf)<5) {
        print(paste0('skipping freq "',l,'" on receiver "',i,'": not enough signals (',nrow(tmp_rf),')'))
        next
      }
      time_seq<-unique(c(round(tmp_rf$timestamp))) #,round(tmp_rf$timestamp)+1,round(tmp_rf$timestamp)-1)
      smoothed<-data.frame(angle=predict(
        smooth.spline(tmp_rf$timestamp,tmp_rf$angle,spar=spar_value),
        as.numeric(time_seq))$y,
        timestamp=time_seq,
        Station=i,
        strength=predict(
          smooth.spline(tmp_rf$timestamp,tmp_rf$strength,spar=spar_value),
          as.numeric(time_seq))$y,
        antennas=0,
        freq_tag=l,stringsAsFactors = FALSE)
      smoothed_data<-rbind(smoothed_data,smoothed)
      if(progress)
        incProgress(amount=1/num_tags)
    }
    cnt_recs<-cnt_recs+1
  }
  return(smoothed_data)
}

#' Calculates an average position using a moving average for the x and y coordinates.
#'
#' @author Dipl.-Phys. Ralf Zeidler, \email{ralf.zeidler@fridata.de}
#'
#' @export
#' @importFrom foreach "%dopar%"
#' @importFrom foreach foreach
#' @importFrom stats median
#'
#' @param tri_data Data frame with triangulated positions. Expects columns \code{timestamp,freq_tag, pos.X, pos.utm.X, pos.utm.Y}.
#' @param time aggregated time
#' @param s_time Time window the average is applied to.
#' @param method What average to be used: "mean" or "median".
#' @return Data frame holding the moving average of the positions. Holds columns \code{timestamp, freq_tag, pos.utm.X, pos.utm.Y, utm.zone, X, Y}.
#'
centroid_fun <- function(tri_data,time,s_time,method="mean"){
  min_time<-min(tri_data$timestamp)
  max_time<-max(tri_data$timestamp)
  time_seq<-seq(round(min_time,"mins")-60,round(max_time,"mins")+60,by=60*time)
  utm<-foreach(i=time_seq,
               .combine=rbind,
               .inorder=FALSE) %dopar% {
                 tmp<-subset(tri_data,timestamp>=i&timestamp<i+60*s_time)
                 zone<-(floor((tmp$pos.X[1] + 180)/6) %% 60) + 1
                 if(nrow(tmp)>0){
                   data.frame(timestamp=i,
                              freq_tag=tmp$freq_tag[1],
                              pos.utm.X=switch(method,
                                               mean=mean(tmp$pos.utm.X),
                                               median=median(tmp$pos.utm.X)),
                              pos.utm.Y=switch(method,
                                               mean=mean(tmp$pos.utm.Y),
                                               median=median(tmp$pos.utm.Y)),
                              utm.zone=zone)
                 }

               }
  location_wgs<-utmtowgs(utm$pos.utm.X,utm$pos.utm.Y,utm$utm.zone)
  return(cbind(utm,pos=location_wgs))
}
