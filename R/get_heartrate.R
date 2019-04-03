#' Preprocess and extract heart rate from smartphone video recordings.
#' 
#' A convenience wrapper for extracting heart rate features for each color
#' band from average pixel value per frame of video (processed hr) captured
#' using smartphone cameras.
#' 
#' The heartrate assay entails participants placing their finger over the 
#' camera for a period of time with the flash on.
#' 
#' @param heartrate_data A data frame with columns t, red, green and blue
#' @param window_length Length of the time window \emph{in seconds}, to be
#' considered while calculating the heart rate for each channel.
#' @param window_overlap Fraction in the interval [0, 1) specifying the amount of
#' window overlap.
#' @param method The algorithm used to calculate the heartrate, current methods
#' include ('acf','psd') which stand for autocorrelation function, and power
#' spectral density respectively. We will be adding support for peak picking 
#' algorithms, and wavelet methods later. 
#' 
#' @return list containing heart rate and confidence of the estimate for
#' each color (red, green, blue)
#' @seealso \code{\link{heartrate_data}}
#' @export
#' @author Meghasyam Tummalacherla, Phil Snyder 
#' @examples 
#' heartrate_data = heartrate_data[,c('t', 'red', 'green', 'blue')]
#' heartrate_ftrs = get_heartrate(heartrate_data)
#'  
get_heartrate <- function(heartrate_data,
                          window_length = 10,
                          window_overlap = 0.5,
                          method = 'acf') {
  ## We will throw away ~5s worth of data(180 samples) after filtering,
  ## keep this in mind
  
  heartrate_error_frame <- data.frame(red = NA, green = NA, blue = NA,
                                      error = NA, sampling_rate = NA)
  sampling_rate <- mhealthtools:::get_sampling_rate(heartrate_data)
  if (is.infinite(sampling_rate) || is.na(sampling_rate)) {
    heartrate_error_frame$error <- paste("Sampling Rate calculated from timestamp is Inf",
                                         "or NaN / timestamp not found in json")
    return(heartrate_error_frame)
  }
  
  # Convert window length from seconds to samples
  window_length <- round(sampling_rate * window_length)
  mean_filter_order <- 65
  if(sampling_rate <= 32){
    mean_filter_order <- 33
  }
  if(sampling_rate <= 18){
    mean_filter_order <- 19
  }
  if(sampling_rate <= 15){
    mean_filter_order <- 15
  }
  
  ##  Apply pre-processing filter to all heartrate data
  
  # We do so as we are using an IIR (running filter), so we do not need
  # to filter each window, as for a running filter the effects are local
  # ARMA filter's output depends only on the y_i's and x_i's used
  # It is also inline with our mean centering filter, since that is also
  # a running filter
  
  heartrate_data <- tryCatch({
    heartrate_data %>% 
      dplyr::select(red, green, blue) %>% 
      na.omit() %>% 
      lapply(get_filtered_signal,
             sampling_rate,
             mean_filter_order,
             method) %>% 
      as.data.frame()
  }, error = function(e){NA})
  if (all(is.na(heartrate_data))) {
    heartrate_error_frame$error <- "Error in filtering the signal"
    return(heartrate_error_frame)
  }
  
  
  # Split each color into segments based on window_length
  heartrate_data <- tryCatch({
    heartrate_data %>%
      dplyr::select(red, green, blue) %>%
      na.omit() %>%
      lapply(mhealthtools:::window_signal, window_length, window_overlap, 'rectangle')
  }, error = function(e) { NA })
  if (all(is.na(heartrate_data))) {
    heartrate_error_frame$error <- "red, green, blue cannot be read from JSON"
    return(heartrate_error_frame)
  }
  
  # Get HR for each filtered segment of each color
  heartrate_data <- heartrate_data %>%
    lapply(function(dfl) {
      dfl <- tryCatch({
        apply(dfl, 2, get_hr_from_time_series, sampling_rate, method)
      }, error = function(e) {c(hr= NA, confidence = NA) })
      dfl <- as.data.frame(t(dfl))
      colnames(dfl) <- c("hr", "confidence")
      return(dfl)
    })
  heartrate_data$error <- "none"
  
  if (sampling_rate < 55) {
    heartrate_data$error <- "Low sampling rate, at least 55FPS needed"
  }
  heartrate_data$sampling_rate <- sampling_rate
  return(heartrate_data)
}

#' Bandpass and sorted mean filter the given signal
#'
#' @param x A time series numeric data
#' @param sampling_rate The sampling rate (fs) of the signal
#' @param mean_filter_order The number of samples used in the sliding window
#' for the mean filtering function
#' @param method The algorithm used to estimate the heartrate, because the
#' preprocessing steps are different for each. method can be any of 
#' 'acf','psd' or 'peak' for algorithms based on autocorrelation, 
#' power spectral density and peak picking respectively
#' @return The filtered time series data
get_filtered_signal <- function(x,
                                sampling_rate,
                                mean_filter_order = 65,
                                method = 'acf') {
  
  # Defaults are set for 60Hz sampling rate
  x[is.na(x)] <- 0
  x <- x[round(3*sampling_rate):length(x)]
  # Ignore the first 3s
  
  sampling_rate_rounded <- round(sampling_rate)
  # Filter the signal based on fiters designed
  if(sampling_rate_rounded > 15){
    bf_low <- signal::butter(7, 5/(sampling_rate_rounded/2), type = 'low')
    bf_high <- signal::butter(7, 0.5/(sampling_rate_rounded/2), type = 'high')
  }else{
    bf_low <- signal::butter(7, 4/(sampling_rate_rounded/2), type = 'low')
    bf_high <- signal::butter(7, 0.5/(sampling_rate_rounded/2), type = 'high')
  }
  
  x <- signal::filter(bf_low, x) # lowpass
  x <- x[round(sampling_rate):length(x)] # 1s
  x <- signal::filter(bf_high, x) # highpass
  x <- x[round(sampling_rate):length(x)] # 1s @ 60Hz
  
  y <- x
  
  #################
  ## Mean centering filter design (For 60Hz Sampling Rate)
  ## The purpose of this is to make sure the waveform is uniform in range
  #################
  if(method == 'acf' || method == 'psd'|| method == 'peak'){
    y <- 0 * x
    sequence_limits <- seq((mean_filter_order + 1) / 2,
                           length(x) - (mean_filter_order - 1) / 2, 1)
    for (i in sequence_limits) {
      temp_sequence <- x[seq(i - (mean_filter_order - 1) / 2,
                             (i + (mean_filter_order - 1) / 2),1)]
      
      # y[i] <- (((x[i] - max(temp_sequence) - min(temp_sequence)) -
      #             (sum(temp_sequence) - max(temp_sequence)) / (mean_filter_order - 1)) /
      #            (max(temp_sequence) - min(temp_sequence) + 0.00001))
      y[i] <- (((x[i]) -
                  (sum(temp_sequence) - (max(temp_sequence)) + min(temp_sequence)) / (mean_filter_order - 2)) /
                 (max(temp_sequence) - min(temp_sequence) + 0.00001))
      
      # y[i] <- (x[i] - (max(temp_sequence) + min(temp_sequence))*0.5)/(max(temp_sequence) - min(temp_sequence))
      
      if(method == 'peak'){
        y[i] = (y[i]*(sign(y[i])+1)/2)
        y[i] = (y[i])^0.15
        y[i] = exp(y[i])^0.75
      }
    }
    y <- y[sequence_limits]
  }
  return(y)
}

#' Given a processed time series find its period using autocorrelation
#' and then convert it to heart rate (bpm)
#'
#' @param x A time series numeric data
#' @param sampling_rate The sampling rate (fs) of the time series data
#' @param method The algorithm used to estimate the heartrate, because the
#' preprocessing steps are different for each. method can be any of 
#' 'acf','psd' or 'peak' for algorithms based on autocorrelation, 
#' power spectral density and peak picking respectively
#' @param min_hr Minimum expected heart rate
#' @param max_hr Maximum expected heart rate
#' @return A named vector containing heart rate and the confidence of the result 
get_hr_from_time_series <- function(x, sampling_rate, method = 'acf', min_hr = 45, max_hr=210) {
  x[is.na(x)] <- 0
  
  if(method == 'acf'){
    
    max_lag = round(60 * sampling_rate / min_hr) # 4/3 fs is 45BPM
    min_lag = round(60 * sampling_rate / max_hr) # 1/3.5 fs is 210BPM
    
    x <- stats::acf(x, lag.max = max_lag, plot = F)$acf
    # x <- x-min(x)
    y <- 0 * x
    y[seq(min_lag, max_lag)] <- x[seq(min_lag, max_lag)]
    # y <- y-min(y)
    y_max_pos <- which.max(y)
    y_max <- max(y)
    y_min <- min(y)
    
    hr_initial_guess <- 60 * sampling_rate / (y_max_pos - 1)
    aliasedPeak <- getAliasingPeakLocation(hr = hr_initial_guess,
                                           actual_lag = y_max_pos,
                                           sampling_rate = sampling_rate,
                                           min_lag = min_lag,
                                           max_lag = max_lag)
    
    if(!is.na(aliasedPeak$earlier_peak[1])){
      peak_pos <- (y[aliasedPeak$earlier_peak]-y_min) > 0.7*(y_max-y_min)
      peak_pos <- aliasedPeak$earlier_peak[peak_pos]
      if(length(peak_pos>0)){
        hr_vec <- 60 * sampling_rate / (peak_pos - 1)
        # hr_pos <- which.min(abs(hr_vec - hr_initial_guess*(aliasedPeak$Npeaks+1)/(aliasedPeak$Npeaks)))
        hr <- mean(hr_vec,hr_initial_guess*(aliasedPeak$Npeaks+1)/(aliasedPeak$Npeaks))
        # confidence <- y[peak_pos[hr_pos]] / max(x)
        confidence <- mean(y[peak_pos]-y_min)/(max(x)-min(x))
      }else{
        hr <- hr_initial_guess
        confidence <- y_max/max(x)
      }
    }else{
      peak_magnitude_vec <- y[aliasedPeak$later_peak]
      status_flag <- (sum(peak_magnitude_vec > 0.7*y_max) == length(peak_magnitude_vec))
      if(!is.logical(status_flag) || is.na(status_flag)){
        status_flag <- F
      }
      if(status_flag){
        hr <- hr_initial_guess
        confidence <- y_max/max(x)
      }else{
        hr <- NA
      }
    }
  }
  
  if(method == 'psd'){
    x_spec <- mhealthtools:::get_spectrum(
      x, sampling_rate,nfreq = 2^round(log(length(x))/log(2))
    ) %>% dplyr::filter(freq>0.7, freq< 3.3)
    # 0.7Hz = 42BPM, 3.3HZ = 198BPM
    hr <- 60*x_spec$freq[which.max(x_spec$pdf)]
    confidence <- NA
  }
  
  if(method == 'peak'){
    
    x_max <- max(x)
    x_peaks <- pracma::findpeaks(x, minpeakdistance = 60 * sampling_rate / max_hr,
                                 minpeakheight = 0.85*x_max)
    
    # peak cleaning to be done
    
    # sort peaks by time when they occur not by amplitude
    x_peaks <- x_peaks[order(x_peaks[,2]),]
    peak_dist <- diff(x_peaks[,2])
    hr <- 60 * sampling_rate / (mean(peak_dist))
    confidence <- NA
  }
  
  
  # If hr or condidence is NaN, then return hr = 0 and confidence = 0
  if ((length(hr) == 0) || is.null(hr) || is.na(hr)) {
    confidence <- NA
    hr <- NA
  }
  
  return(c(hr, confidence))
}


getAliasingPeakLocation <- function(hr, actual_lag = NA,sampling_rate, min_lag, max_lag){
  # The following ranges are only valid if the minimum hr is 45bpm and
  # maximum hr is less than 240bpm, since for the acf of the ideal hr signal
  # Npeaks = floor(BPM/minHR) - floor(BPM/maxHR)
  # in the search ranges 60*fs/maxHR to 60*fs/minHR samples
  
  if(hr < 90){
    Npeaks = 1
    # 45BPM - 89BPM
  }else if(hr < 135){
    Npeaks = 2
    # 90BPM - 134BPM
  }else if(hr < 180){
    Npeaks = 3
    # 135BPM -179BPM
  }else if(hr < 225){
    Npeaks = 4
    # 180BPM - 225BPM
  }else if(hr <= 240){
    Npeaks = 5
  }
  
  if(is.na(actual_lag)){
    actual_lag = ceiling(sampling_rate*60/hr + 1)
  }
  
  if(actual_lag%%2 == 0){
    earlier_peak <- actual_lag/2
  }else{
    earlier_peak <- c(floor(actual_lag/2), ceiling(actual_lag/2))
  }
  
  
  if(Npeaks > 1){
    later_peak <- actual_lag*seq(2,Npeaks)
    later_peak[later_peak>max_lag] <- NA
  }else{
    later_peak <- NA
  }
  
  earlier_peak[earlier_peak < min_lag] <- NA
  
  return(list(Npeaks = Npeaks,
              earlier_peak = earlier_peak,
              later_peak = later_peak))
}
