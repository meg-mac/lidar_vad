##R function to calculate wind speed and direction from raw lidar VAD data
     

lidar_VAD <- function(Nrays,Ngates,elevation,SNR_thresh,lat,lon,date,file_pattern,trailer_angle)
{

#Values You Dont Need to Change
Lgate = 30*sin(elevation*(pi/180)) # accounting for VAD step-stares not pointing directly vertically
x1 = seq(1,360,1)
x = matrix(x1,ncol = 12, nrow = 30) 
azimuth = matrix(data = seq(0,330,30))
library(rfgfs) #may need to dowload from github
declination=abs(magvar(lat, lon, h = 1, date = date))

#Retrieve a list of files in my data folder:
flist <- list.files(pattern = file_pattern, recursive = FALSE)

i = 1
for (i in 1:length(flist)) {

    unlisted <- sapply(seq(from=1, to=nchar(flist[i]), by=1), function(l)
                substr(flist[i]  ,l, l+5))
    
    #Saving portion of filenames to use later for times
    name2<- unlisted[18] 
    if (exists("name1")){
    name1 <-rbind(name1,name2) }
    if (!exists("name1")){
    (name1 <- name2) }

    #Dealing with file structure see file example on github
    raw_file <- read.table(file = flist[i], fileEncoding="ASCII", dec="\t", skip = 17, fill = TRUE, na.strings=c("","NA"))
    headers <- raw_file[complete.cases(raw_file),]
    n30 <-  c("Decimal_time_hours", "Azimuth_degrees",  "Elevation_degrees", "Pitch_degrees", "Roll_degrees")
    names(headers) <- paste(n30,sep="")
    
    values <- raw_file[!complete.cases(raw_file),]
    
    rep_header <- headers[rep(seq_len(nrow(headers)), each = 200), ]
    values$V5 <- rep_header$Azimuth_degrees
    n30_v <-  c("Range_Gate" , "radial_wind" , "Intensity_SNR_1" , "Beta_m_sr", "azimuth")
    names(values) <- paste(n30_v,sep="")
    
  
    #New Data Frame and Coercing into Numeric  
    radial_wind <- values
    
    radial_wind[] <- lapply(radial_wind, function(x) {
        if(is.factor(x)) as.numeric(as.character(x)) else x })

    
    a = radial_wind
    
    ## importing variables from the loaded .hpl file
    range_gate=a$Range_Gate  # first column of every .hpl file
    altitude=(range_gate+0.5)*Lgate 
    altitude = matrix(altitude,ncol = 12)  
    altitude= altitude[,1]
    
    #second column of every .hpl file    
    radial_velocity=a$radial_wind  
    radial_velocity=matrix(radial_velocity,nrow = Ngates, ncol = Nrays)
     
     # third column of every .hpl file   
    SNR=a$Intensity_SNR_1
    SNR=matrix(SNR,nrow = Ngates,ncol = Nrays) 
    
    #gives index for those range gates for which the average SNR across all Nrays exceeds
    #the imposed SNR threshold. Removes NA values  
    ind_good = which(rowMeans(SNR)>=SNR_thresh)
    ind_good = na.omit(as.numeric(ind_good))
                                                    
    #initialize    
    amp <- (rep(NA,length(ind_good)))
    phi <- rep(NA,length(ind_good))
    
  
    j = 1
    for (j in (1:length(ind_good))){
      # try({ # fitting the sine to only the high-SNR range gates
           
            y=(radial_velocity[ind_good[j],] )
     
            # fitting parameters (may fail during daytime convective conditions due to 
            #demanding 2*pi periodicity) which is why there is a try function
            #x0=c(max(y) ,2*pi/360 ,0.01 , mean(y)) 
            x=azimuth
            
            #nlsLM likes to work in data frames 
            df <- as.data.frame(y)
            df$x <- x
            
            #https://stats.stackexchange.com/questions/60500/how-to-find-a-good-fit-for-semi-sinusoidal-model-in-r
            #https://stats.stackexchange.com/questions/60994/fit-a-sinusoidal-term-to-data
            
            # FFT used to find angular frequency (omega)
            # Anything past the N/2 - 1 element due to the Nyquist-Shannon limit as well as the first value, 
            # which contains no frequency information, were removed. The importance of each frequency and omega 
            # corresponded to the largest absolute value of the FFT. 
            raw.fft = fft(y)
            truncated.fft = raw.fft[seq(1, length(y)/2 - 1)]
            truncated.fft[1] = 0
            omega = which.max(abs(truncated.fft)) * 2 * pi / length(y)
         
            amp[j] = NA
            phi[j] = NA
            
            #Occasionally nlsLM fails which is why there is a try function. If it fails, values remain NA.
            #If you want to see number of times nlsLM fails change silent = FALSE
            try({
              library(minpack.lm)
              fm1DNase2 <- nlsLM(df$y ~ (a*sin((b*x+c)*(pi/180))+d),  
                                      start=list(a = max(y), b= omega,
                                       c= 0.01 ,d = mean(y)) ,data=df, trace = FALSE)
                amp[j] = coef(fm1DNase2)[1]
                phi[j] = coef(fm1DNase2)[3] }, silent = TRUE)  
              
          
    } 
           
    #no scientific notation
    options(scipen=999)
    
    #redefining absolute value amp 
    ws= (abs(amp))
    
    # necessary to constrain occasional negative wind directions
    k = 1
    for (k in (1:length(phi))){
      if(is.na(phi[k])){
        phi[k]= NA}
      else if (phi[k] <= 0){
             phi[k] = -1*phi[k]-180
           }
    }
    
      
    #redefining wind direction 
    wdir=phi        
    
    #initalizing value
    ws_total <- (rep(NA,length(ind_good)))
    wdir_total <- rep(NA,length(ind_good))
    #accounting for VAD step-stares not pointing directly vertically    
    #correction to get the proper wind direction
    
    k = 1
    for (k in 1:length(ind_good)){
          ws_total[k] = ws[k]/cos(elevation*(pi/180))
          wdir_total[k] = wdir[k]+90}
   
          
    ws_result <- rep(NA, length(altitude))
    ws_result[ind_good] <- ws_total
    ws_total <- ws_result
    
    wdir_result <- rep(NA, length(altitude))
    wdir_result[ind_good] <- wdir_total
    wdir_total <- wdir_result
    
  
    #lowest three range gates are the HALO dead zone
    ws_total[1:3] =  NA
    wdir_total[1:3] = NA
    
    #minor QA to remove unreasonable wind speed and direction values (if any)
    ind1 = which(ws_total>20)
    ws_total[ind1] = NA
    wdir_total[ind1] = NA

    ind2 = which(wdir_total < 0 | wdir_total > 360)
    #ws_total[ind2] = NA
    wdir_total[ind2] = NA

    
    # correction (90 degrees stem from the HALO home position not aligned with the trailer azimuth)
    azimuth_trailer=90-trailer_angle+declination; 
    
    
    wdir_total_corr <- rep(NA, length(altitude))
    
    #multiple if statements here to account for the cyclic character of wind direction
    for (i in 1:200){
            wdir_total_corr[i] = wdir_total[i]-azimuth_trailer
            if (is.na(wdir_total_corr[i])){
              wdir_total_corr[i] = NA
            }else if (wdir_total_corr[i] > 360){
                wdir_total_corr[i]=wdir_total_corr[i]-360
            }else if (wdir_total_corr[i]<0){
                wdir_total_corr[i]=360+wdir_total_corr[i]
            }
    
    }



    #Saving Variables
    wsp2 <- ws_total
    if (exists("wsp1")){
    wsp1 <-cbind(wsp1,wsp2) }
    if (!exists("wsp1")){
    (wsp1 <- wsp2) }
    
    alt2 <- altitude
    if (exists("alt1")){
    alt1 <-cbind(alt1,alt2) }
    if (!exists("alt1")){
    (alt1 <- alt2) }
    
    wdir2 <- wdir_total_corr
    if (exists("wdir1")){
    wdir1 <-cbind(wdir1,wdir2) }
    if (!exists("wdir1")){
    (wdir1 <- wdir2) }
    
       

}

colnames(wsp1) <- name1
wind_speed <<- wsp1

altitude <<- alt2

colnames(wdir1) <- name1
wind_direction <<- wdir1
}
