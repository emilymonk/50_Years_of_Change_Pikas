#Code to calculate temperature statistics as in Krear's thesis, visualize, and summarize
#This code uses a metadata file that has all sensor names, temperature info, and time info
#All sensor files are read in from individual files, using the metadata file.

#Code by Emily Monk and Chris Ray
#January 2022

dev.off()
require(tidyverse)
require(mgcv) #for fitting a gam
require(lubridate) #for dealing with dates


#set working directory first

#This is the metadata file
files <- read.csv("Krear_Site_Sensor_Metadata.csv",stringsAsFactors=F)
filenames <- files[,1]
time.format <- files[,2]

#Note that the code in its current state will exclude "E_XIreplacement_2013_trim" and "PZ_2015_trim"
# as those sensors were close replacements or not in place long enough

###############################################################

#function f to c
get.temp.C <- function(x){
  temp.C <- x[,2] #* fills tempC with whatever's in the 2nd column right now
  if (names(x)[2]=="Temp.F"){
    temp.C <- (temp.C-32)*5/9 #* revises if needed
  }
  return(temp.C) #* returns the original or revised version
}


################################################################
################################################################
#Finds for a given month: 
# -"Days Above Freezing", 
# -"Mean Daily Max", 
# -"Mean Daily Min",
# -"Month Max"
# -"Month Min"
# -"Mean of Daily Max and Min"

#Based on a daily fit GAM
daily.fit.stats <- function(sensorfile="",xv,ys,month="NA",knots=6,t.min=-30,t.max=40){
  
  #* get a basic plot started by plotting all data points
  #* this allows us to visualize the model
  plot(xv,ys,col="purple",cex=0.75,bty="n",las=1,
       ylim=c(t.min,t.max),xlim=c(0,23),xaxt="n",main="",
       xlab="Hour",ylab=expression(paste("Temperature ( "^"o","C)")))
  axis(1,at=seq(0,23,by=2),labels=seq(0,23,by=2))
  
  md <- unique(xv$mday) #* these are the days: we will fit a gam through each days' temps
  daily.max <- rep(NA,length(md))#* will hold the fitted daily maxima
  abs.daily.max <- rep(NA,length(md))
  daily.min <- rep(NA,length(md))
  
  for (ii in md[1:length(md)]){ #* for every day in this month...
    
    
    jj <- which(xv$mday==ii) #* indexes the temperature records for this day
    xvv <- xv$hour[jj] #* grabs the hours associated w/temperature records for this day
    
    if (abs(max(ys[jj],na.rm=T)-min(ys[jj],na.rm=T))>2){ #* fit a curve only if temps vary
      
      mod.gam <- gam(ys[jj]~s(xvv,k=knots,bs="cc")) # fit, cyclic cubic spline
      
      xvv <- seq(min(xvv),(max(xvv)+1),by=1/60) # xvv holds every minute of the day
      
      yps <- predict(mod.gam,newdata=list(xvv=xvv)) # yps is the model
      
      lines(xvv,yps,col="purple",lwd=1) # show fitted y
      
      # plot the fitted max and min temps
      day.max <- max(yps) # this is the max of the model
      xx <- xvv[which(yps==day.max)] # the time corresponding to max temp
      points(xx[1],day.max,pch=1,col="orange",cex=1.5) # plotting a point for the max
      
      day.min <- min(yps) # minimum of the model
      xx <- xvv[which(yps==day.min)] # the time corresponding to min temp
      points(xx[1],day.min,pch=1,col="turquoise3",cex=1.5) # plot the min
      
      abs.daily.max[ii] <- day.max #* store this daily max
      daily.max[ii] <- day.max
      daily.min[ii] <- day.min
      
    } else {
      
      #when there is no variation in temps fitting a GAM doesn't work well, 
      #and there are no peaks/troughs to sample so the fit gets weird
      abs.daily.max[ii] <- max(ys) 
      daily.max[ii] <- max(ys)
      daily.min[ii] <- min(ys)
    }
  }
  
  ifelse(sensorfile=="",mysep<-"",mysep<-": ")
  
  text(0,t.max,paste(sensorfile,month,sep=mysep),cex=1.25,adj=0)
  
  legend("bottomright",bty="n",pch=c(NA,1,1),
         pt.cex=c(NA,1.5,1.5),col=c("purple","orange","turquoise3"),
         lty=c(1,NA,NA),lwd=c(1,NA,NA),
         legend=c("GAM fit",
                  "Daily maximum", # max for THIS DAY
                  "Daily minimum"))
  
  #Calculate days above freezing
  abf <- (sum(abs.daily.max>0,na.rm=T))
  #Calculate mean daily maximum
  mean.mon.max <- mean(daily.max, na.rm=T)
  #Calculate mean daily minimum
  mean.mon.min <- mean(daily.min, na.rm=T)
  #Calculate month max and min
  mon.max <- max(daily.max, na.rm = T)
  mon.min <- min(daily.min, na.rm = T)
  #Calculate the mean of daily max and min
  mean.max.min <- mean(c(daily.max, daily.min), na.rm = T)
  
  stat.names <- c("Days Above Freezing", "Mean Daily Max", "Mean Daily Min", "Month Maximum", "Month Minimum", "Mean of Daily Max and Min")
  
  value <- c(abf, mean.mon.max,mean.mon.min, mon.max, mon.min, mean.max.min)  #these are the actual values to go in the data frame
  #create the data frame with the monthly stats
  df <- data.frame(month=month,stat.names,value=round(value,1))
  return(df)
  
}



##########################################################

#This section produces a file with all monthly statistics for all sensor years.


#create empty data frame
mon_ab_df <- c("Aug", "Sep", "Oct", "Nov", "Dec","Jan", "Feb", "Mar", "Apr", "May",
               "Jun", "Jul","Aug", "Sep", "Oct")

header <- c("sensor.year","statistic", mon_ab_df) 
output <- matrix(NA,nrow=7,ncol=length(header))
df <- data.frame(output,stringsAsFactors=F)
colnames(df) <- header
df[,2] <- c("Month-Year",
            "Days above freezing", 
            "Mean daily max T", "Mean daily min T",
            "Month Maximum", "Month Minimum",
            "Mean of Daily Max and Min")
df

#* now loop through the files

for (i in 1:length(files$filename.new)){
  
  sensor<- files$filename.new[i]
  
  x <- read.csv(sensor,stringsAsFactors=F,header = TRUE)
  x[,2] <- get.temp.C(x)
  
  #re-name the first 2 columns
  colnames(x)[1:2] <- c("date.time","temp.C")
  
  if (time.format[i]=="12hr") form <- "%m/%d/%Y %I:%M:%S %p" 
  
  if (time.format[i]=="24hr") form <- "%m/%d/%Y %H:%M" 
  
  dt <- strptime(x$date.time,format=form) 
  #Note that there was some issues with the time format and which format included second and which didn't
  #The way the code is set up currently uses seconds for the 12 hour time format but not the 24 hour time 
  #There will likely be an error if future sensors don't follow this same comvention with the seconds.
  
  t.max <- max(x$temp.C)  #max for entire sensor year
  t.min <- min(x$temp.C)  #min for entire sensor year
  
  #df is the mostly empty data frame; put sensor name into rows under "sensor.year" col name
  df[1:7,c("sensor.year")] <- files$filename.new[i] # name the sensor-year
  
  
  # make sure all columns of monthly stats are set to "NA" before saving the stats unique to this sensor-year;
  # this is important for all sensor-years after the first one, so we don't have stats carrying over from one
  # sensor-year to the next, which happens when sensor-year t contains data for more months than sensor-year t+1
  df[1:7,3:17] <- NA
  
  #* now find the month-years in this sensor-year file
  mo.yr <- paste(month.abb[dt$mon+1],dt$year+1900,sep="")
  mo.yrs <- unique(mo.yr)
  mo.yrs
  
  #This is the part that should tell R where to start putting the data in the data frame.
  #* now we know we need to put month-year "Sep2009" under the first "Sep" column, etc;
  #* so we use grep() to figure out which column holds our very first month-year, and
  #* we further use length(grep()) to give us a true/false value for whether each of
  #* 4 possible months is the first month in our data for this sensor-year:
  ii <- 0 #* ii is the column index we're trying to find, set at zero initially
  if (length(grep("Aug",mo.yrs[1]))) ii <- 3 #* column = 3 only if mo.yrs[1] contains "Aug"
  if (length(grep("Sep",mo.yrs[1]))) ii <- 4 #* column = 4 only if mo.yrs[1] contains "Sep", etc.
  if (length(grep("Oct",mo.yrs[1]))) ii <- 5
  if (length(grep("Nov",mo.yrs[1]))) ii <- 6
  #* Sensors were never placed later than early Nov, so we don't have to worry about Dec, etc;
  
  #* now position all mo.yrs within the 1st row of df, starting at column ii
  df[1,ii:(ii+length(mo.yrs)-1)] <- mo.yrs
  #df #* looks to see if any bugs
  
  #* now loop through the temperature data by month-year
  
  for (my in mo.yrs){
    
    
    #j indexes all the positions in mo.yr that "my" currently represents (in the current iteration of the loop)
    j <- which(mo.yr==my)
    
    #* report whether this "my" has enough data to consider from this sensor-year
    data.days <- length(unique(dt$mday[j]))
    ii <- which(df[1,]==my); ii #* ii indexes the position of "my" in the top ROW (not header) of df
    if (data.days < 20) df[2:7,ii] <- "low data"
    
    #next means to skip the current iteration of a loop without terminating it.
    #So, if a month has less than 20 days in it, we are not going to use it in our analyses.
    if (data.days < 20) next
    
    xv <- dt[j]$hour #* xv is now the hour from each date-time within all data.days of this month-year
    
    ys <- x$temp.C[j] #* ys is now the temperature "
    #* for a 31-day month, both xv and ys will each be 31*k long, where k = # of temps logged daily
    
    k <- min(6,length(unique(xv))) #* k is # of temperatures logged daily
    
    stats <- daily.fit.stats(sensorfile=files$filename.new[i],xv=dt[j],ys=x$temp.C[j],month=my,knots=k)
    
    
    ii <- which(df[1,]==my); ii #* ii indexes the position of "my" in the top ROW (not header) of df
    df[2:7,ii] <- as.character(stats$value)
    df
    
  }
  suppressWarnings(write.csv(df,"revised_temp-stats-table_daily_fits.csv",TRUE))
  
}

#########################################################################################
#########################################################################################
#Now graph all values found above:

require(tidyverse) #tidy data and pipes used frequently
library(weathermetrics) #fahrenheit.to.celsius() function
library(naniar) #replace_with_na() function
library(readxl) #needed because Krear's data used here is an excel file
library(zoo) #as.yearmon is used

#This is the metadata file
files <- read.csv("Krear_Site_Sensor_Metadata.csv",stringsAsFactors=F)

#Column one has the names of files, name in file is filename.new
filenames <- files[,1]
#column with the info about if the times are in 24 or 12 hour format.
time.format <- files[,2]
#column with info about sensor position
s.position <- files[,3]

#This is the actual data with all stats to match Krear
data <- read.csv("revised_temp-stats-table_daily_fits.csv") 


#filter out rows
df <- filter(data, !(statistic %in% c("Month-Year", "statistic")))

#Make long
df <- pivot_longer(df, Aug:length(df), names_to = "Month", values_to = "values")

#Now load in Krear's data frame:
krear <- read_excel("Krear OG Data Formatted.xlsx", sheet = 3)
head(krear)


#Make long
krear <- pivot_longer(krear, Jul63:length(krear), names_to = "Month", values_to = "values")

#Make the month names match with our data for graphing purposes:
krear$Month <- gsub("6.","",krear$Month)

krear <- krear %>%
  mutate(values = as.numeric(values)) %>%
  mutate(Month = as.factor(Month))

#need to convert Krear's f to celcius, excluding the days above freezing count values
krear$values[krear$statistic != "Days Above Freezing"] <- fahrenheit.to.celsius(krear$values[krear$statistic != "Days Above Freezing"], 2)


graph_jul2Jun_f <- function(s.location, files, df, stats_to_graph, krear_stats_to_graph, graph_title ,stat_labels, y.lab, y.axis){
  krear_loc <- krear %>%
    filter(location == s.location)
  
  #Get the sensor year names from the files df for each location type
  loc_sensors <- files %>%
    filter(sensor.type %in% c(s.location)) %>%
    pull(filename.new)
  
  #now get only the data for the sensors at the location
  loc_recent_df <- df %>%
    filter(sensor.year %in% loc_sensors)
  
  #Deal with the duplicate months. 
  loc_recent_df$Month <- gsub("\\..*","",loc_recent_df$Month)
  
  
  loc_recent_df <- loc_recent_df %>% 
    replace_with_na(replace = list(values = "low data")) %>%
    mutate(values = as.numeric(values))
  
  graph.data <- subset(loc_recent_df, statistic %in% stats_to_graph)
  
  #Re-order the factor levels
  graph.data$Month <- factor(graph.data$Month, levels = c("Jul","Aug", "Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun"))
  
  #convert the months  to numeric
  graph.data$Month.num <- match(graph.data$Month,month.abb) 
  krear_loc$Month <- match(krear_loc$Month, month.abb)
  
  #now turn the months into date class, using pretend dates to get around the issue.
  graph.data$pretend_year <-  ifelse(graph.data$Month.num > 6, "2009", "2010")
  graph.data$Date <- as.yearmon(paste(graph.data$pretend_year, graph.data$Month.num), "%Y %m")
  graph.data$Date <- as.Date(graph.data$Date)
  
  krear_loc$pretend_year <-  ifelse(krear_loc$Month > 6, "2009", "2010")
  krear_loc$Date <- as.yearmon(paste(krear_loc$pretend_year, krear_loc$Month), "%Y %m")
  krear_loc$Date <- as.Date(krear_loc$Date)
  
  if(stats_to_graph == "Days above freezing") {
    
    g_poisson <- ggplot(graph.data, 
                        aes(x = Date, y = values)) +
      geom_point(aes(color = statistic, shape = statistic)) +
      theme_minimal() +
      #use this geom to do the gam calculation within the ggplot code. Plotted line is exactly the same as outside/normal calculation.
      stat_smooth(data = (subset(graph.data, statistic %in% stats_to_graph[1])), aes(x=Date, y = values), method="gam",  
                  method.args = list(family = "poisson", link="log"),
                  formula= y ~ s(x, bs = "cs"), 
                  inherit.aes = F, color = "#F8766D", fill = "#F8766D", alpha = .35) +
      
      labs(y = y.lab, title = graph_title) +
      #modify the color legend for the recent years stats
      scale_color_discrete(name = " Recent Years Statistic", labels = stat_labels)+ 
      #change the x axis labels from numeric to months
      scale_x_date(breaks = graph.data$Date, date_labels="%b") +
      ylim(y.axis) +
      theme(axis.title.x = element_blank(), axis.text = element_text(size = 12), axis.title.y = element_text(size = 14), 
            plot.title = element_text(size = 16), legend.title = element_text(size = 14), legend.text = element_text(size = 12))+
      
      #now add the Krear things to the graph:
      geom_point(data = subset(krear_loc, statistic %in% krear_stats_to_graph), 
                 aes(x = Date, y = values, shape = statistic)) +
      geom_line(data = subset(krear_loc, statistic %in% krear_stats_to_graph), 
                aes(x = Date, y = values, group = statistic))  +
      #modify the shape legend to only have historical labels in it.
      scale_shape_manual(values=c(16, 16), name = "Historical Data Statistic", breaks = stats_to_graph) 
    
    print(g_poisson)
    
  }
  
  else{
    
    
    if (length(stats_to_graph) > 1.1){
      
      
      graph <- ggplot(graph.data, 
                      aes(x = Date, y = values)) +
        geom_point(aes(color = statistic, shape = statistic)) +
        theme_minimal() +
        #use this geom to do the gam calculation within the ggplot code. Plotted line is exactly the same as outside/normal calculation.
        stat_smooth(data = (subset(graph.data, statistic %in% stats_to_graph[1])), aes(x=Date, y = values), method="gam",  formula= y ~ s(x, bs = "cs"), inherit.aes = F, color = "#F8766D", fill = "#F8766D", alpha = .35) +
        stat_smooth(data = (subset(graph.data, statistic %in% stats_to_graph[2])), aes(x=Date, y = values), method="gam",  formula= y ~ s(x, bs = "cs"), inherit.aes = F, color = "#00BFC4", fill = "#00BFC4", alpha = .35) +
        
        labs(y = y.lab, title = graph_title) +
        #modify the color legend for the recent years stats
        scale_color_discrete(name = " Recent Years Statistic", labels = stat_labels)+ 
        #change the x axis labels from numeric to months
        scale_x_date(breaks = graph.data$Date, date_labels="%b") +
        ylim(y.axis) +
        theme(axis.title.x = element_blank(), axis.text = element_text(size = 12), axis.title.y = element_text(size = 14), 
              plot.title = element_text(size = 16), legend.title = element_text(size = 14), legend.text = element_text(size = 12))+
        
        #now add the Krear things to the graph:
        geom_point(data = subset(krear_loc, statistic %in% krear_stats_to_graph), 
                   aes(x = Date, y = values, shape = statistic)) +
        geom_line(data = subset(krear_loc, statistic %in% krear_stats_to_graph), 
                  aes(x = Date, y = values, group = statistic))  +
        #modify the shape legend to only have historical labels in it.
        scale_shape_manual(values=c(16, 17, 16, 17), name = "Historical Data Statistic", breaks = stats_to_graph, labels = stat_labels) 
      
      
      print(graph)} else{
        
        
        
        graph <- ggplot(graph.data, 
                        aes(x = Date, y = values)) +
          geom_point(aes(color = statistic, shape = statistic)) +
          theme_minimal() +
          #use this geom to do the gam calculation within the ggplot code. Plotted line is exactly the same as outside/normal calculation.
          stat_smooth(data = (subset(graph.data, statistic %in% stats_to_graph[1])), aes(x=Date, y = values), method="gam",  formula= y ~ s(x, bs = "cs"), inherit.aes = F, color = "#F8766D", fill = "#F8766D", alpha = .35) +
          
          labs(y = y.lab, title = graph_title) +
          #modify the color legend for the recent years stats
          scale_color_discrete(name = " Recent Years Statistic", labels = stat_labels)+ 
          #change the x axis labels from numeric to months
          scale_x_date(breaks = graph.data$Date, date_labels="%b") +
          ylim(y.axis) +
          theme(axis.title.x = element_blank(), axis.text = element_text(size = 12), axis.title.y = element_text(size = 14), 
                plot.title = element_text(size = 16), legend.title = element_text(size = 14), legend.text = element_text(size = 12))+
          
          #now add the Krear things to the graph:
          geom_point(data = subset(krear_loc, statistic %in% krear_stats_to_graph), 
                     aes(x = Date, y = values, shape = statistic)) +
          geom_line(data = subset(krear_loc, statistic %in% krear_stats_to_graph), 
                    aes(x = Date, y = values, group = statistic))  +
          #modify the shape legend to only have historical labels in it.
          scale_shape_manual(values=c(16, 16), name = "Historical Data Statistic", breaks = stats_to_graph) 
        
        print(graph)
      }
  }
  
}
dev.off()
pdf("Recent vs Krear Comparison Plots, GAM Fits All.pdf")

#function input order:
#function(s.location, files, df, stats_to_graph, krear_stats_to_graph, graph_title ,stat_labels, y.lab)

#Max and Min stats
graph_deep_month_maxmin <- graph_jul2Jun_f("deep", files, df, c("Month Maximum", "Month Minimum"), c("Max.T", "Min.T"), "Deep Sensors, Monthly Max and Min",c("Monthly Maximum", "Monthly Minimum"), "Temperature (C)", c(-28,31))
graph_shallow_month_maxmin <- graph_jul2Jun_f("shallow", files, df, c("Month Maximum", "Month Minimum"), c("Max.T", "Min.T"), "Shallow Sensors, Monthly Max and Min",c("Monthly Maximum", "Monthly Minimum"), "Temperature (C)", c(-28,31))
graph_air_month_maxmin <- graph_jul2Jun_f("free air", files, df, c("Month Maximum", "Month Minimum"), c("Max.T", "Min.T"), "Free Air Sensors, Monthly Max and Min",c("Monthly Maximum", "Monthly Minimum"), "Temperature (C)", c(-28,31))

#Mean Daily Max and Min stats
graph_deep_daily2_maxmin <- graph_jul2Jun_f("deep", files, df, c("Mean daily max T", "Mean daily min T"), c("Mean.Daily.Max.T", "Mean.Daily.Min.T"), "Deep Sensors, Mean Daily Max and Min",c("Mean Daily Maximum", "Mean Daily Minimum"), "Temperature (C)", c(-28,31))
graph_shallow_daily2_maxmin <- graph_jul2Jun_f("shallow", files, df, c("Mean daily max T", "Mean daily min T"), c("Mean.Daily.Max.T", "Mean.Daily.Min.T"), "Shallow Sensors, Mean Daily Max and Min",c("Mean Daily Maximum", "Mean Daily Minimum"), "Temperature (C)", c(-28,31))
graph_air_daily2_maxmin <- graph_jul2Jun_f("free air", files, df, c("Mean daily max T", "Mean daily min T"), c("Mean.Daily.Max.T", "Mean.Daily.Min.T"), "Free Air Sensors, Mean Daily Max and Min",c("Mean Daily Maximum", "Mean Daily Minimum"), "Temperature (C)", c(-28,31))

#Mean of Daily Max and Min stats
graph_deep_mean_max_min <- graph_jul2Jun_f("deep", files, df, "Mean of Daily Max and Min", "Mean.Daily.MaxMin.T", "Deep Sensors, Mean of Daily Max and Min", "Mean of Daily Max and Min", "Temperature (C)", c(-28,31))
graph_shallow_mean_max_min <- graph_jul2Jun_f("shallow", files, df, "Mean of Daily Max and Min", "Mean.Daily.MaxMin.T", "Shallow Sensors, Mean of Daily Max and Min", "Mean of Daily Max and Min", "Temperature (C)", c(-28,31))
graph_air_mean_max_min <- graph_jul2Jun_f("free air", files, df, "Mean of Daily Max and Min", "Mean.Daily.MaxMin.T", "Free Air Sensors, Mean of Daily Max and Min", "Mean of Daily Max and Min", "Temperature (C)", c(-28,31))

#Days above Freezing stat
graph_deep_days_abv_0 <- graph_jul2Jun_f("deep", files, df, "Days above freezing", "Days Above Freezing", "Deep Sensors, Number of Days Above Freezing", "Days above Freezing", "Number of Days Above Freezing", c(-3,40))
graph_shallow_days_abv_0 <- graph_jul2Jun_f("shallow", files, df, "Days above freezing", "Days Above Freezing", "Shallow Sensors, Number of Days Above Freezing", "Days above Freezing", "Number of Days Above Freezing", c(-3,40))
graph_air_days_abv_0 <- graph_jul2Jun_f("free air", files, df, "Days above freezing", "Days Above Freezing", "Free Air Sensors, Number of Days Above Freezing", "Days above Freezing", "Number of Days Above Freezing", c(-3,40))


dev.off()

################################################################################################
################################################################################################
#Produce a table of summary statistics


#This is the actual data with all stats to match Krear
data <- read.csv("revised_temp-stats-table_daily_fits.csv") 


sum_stats_f <- function(data, files, s.location){
  
  #filter out rows
  df <- filter(data, !(statistic %in% c("Month-Year", "statistic")))
  #Make long
  df <- pivot_longer(df, Aug:length(df), names_to = "Month", values_to = "values")
  
  
  #Get the sensor year names from the files df for each location type
  loc_sensors <- files %>%
    filter(sensor.type %in% c(s.location)) %>%
    pull(filename.new)
  
  #now get only the data for the sensors at the location
  df <- df %>%
    filter(sensor.year %in% loc_sensors)
  
  
  #Deal with the duplicate months. 
  df$Month <- gsub("\\..*","",df$Month)
  df <- df %>% 
  replace_with_na(replace = list(values = "low data")) %>%
  mutate(values = as.numeric(values))

  #change to factor and order so output is in a better order
  df$Month <- factor(df$Month, levels = c("Jul","Aug", "Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun"))


  #Now summarize
  sum_stats <- df %>%
    group_by(statistic, Month) %>%
    summarise(mean=mean(values, na.rm=T), sd=sd(values, na.rm=T)) %>%
    ungroup()

  wide_sum_stats <- pivot_wider(sum_stats, names_from = statistic, values_from = c(mean, sd))
  
  return(wide_sum_stats)
  
}

shallow_sum_stats <- sum_stats_f(data, files, "shallow")
write.csv(shallow_sum_stats, "shallow_mean_sd.csv")

deep_sum_stats <- sum_stats_f(data, files, "deep")
write.csv(shallow_sum_stats, "deep_mean_sd.csv")

air_sum_stats <- sum_stats_f(data, files, "free air")
write.csv(shallow_sum_stats, "air_mean_sd.csv")
