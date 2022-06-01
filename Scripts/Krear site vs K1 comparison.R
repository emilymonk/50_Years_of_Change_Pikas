#Krear site vs K1 comparison


require(ggpubr)
require(tidyverse)
library(weathermetrics)
library(naniar)
library(forcats)
library(readxl)
library(lubridate)
require(mgcv)
library(zoo)

setwd("D:/Pika_Data")

#read in files. 
krear_site_metadata <- read.csv("Shallow_talus_Krear_site_vs_K1_2010-2020.csv")
head(krear_site_metadata)

krear_site_metadata <- subset(krear_site_metadata, select = -record)

# krear_site_metadata <- krear_site_metadata %>%
#   filter(!DeploymentID == c("2013-2014-Long3-Aidan-relic", "2014-2015-Long3-Aidan-relic"))

#-----------------------------------------------------------------------------------------------


#read in actual data

krear_site <- read.csv("Krear_site_recent_shallow_2010-2020.csv")

krear_site <-  krear_site %>%
  filter(DeploymentID %in% unique(krear_site_metadata$DeploymentID))
#   filter(!DeploymentID == c("2013-2014-Long3-Aidan-relic", "2014-2015-Long3-Aidan-relic"))

krear_site$date.time <- strptime(krear_site$date.time, "%m/%d/%Y %H:%M", tz = "MST") 

head(krear_site)
str(krear_site)

unique(krear_site$DeploymentID)

#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------


daily.fit.stats <- function(sensorfile="",xv,ys,month="NA",knots=6,t.min=-30,t.max=40, min_knots){
  
  #* for debugging start with the function below, run all of that with the debugging lines then go here.
  
  #* get the plot started by plotting all data points
  plot(xv,ys,col="purple",cex=0.75,bty="n",las=1,
       ylim=c(t.min,t.max),xlim=c(0,23),xaxt="n",main="",
       xlab="Hour",ylab=expression(paste("Temperature ( "^"o","C)")))
  axis(1,at=seq(0,23,by=2),labels=seq(0,23,by=2))
  
  md <- unique(xv$mday) #* these are the days: we will fit a gam through each days' temps
  daily.max <- rep(NA,length(md))#* will hold the fitted daily maxima
  abs.daily.max <- rep(NA,length(md))
  daily.min <- rep(NA,length(md))
  
  for (ii in md[1:length(md)]){ #* for every day in this month...
    
    #ii <- 8 #debugger
    jj <- which(xv$mday==ii) #* indexes the temperature records for this day
    
    if(length(jj)<min_knots) next
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

#This section is mostly copied and pasted from Krear script, 
# but I made a number of modifications to make it work with the format of the surrounding site data.

# Most significant difference is that this is all another function now, calling the above function.

sur_site_stats_ftn <- function(site_meta, site_data, csv_stats_title, min_knots){
  #create empty data frame
  mon_ab_df <- c("Jun", "Jul","Aug", "Sep", "Oct", "Nov", "Dec","Jan", "Feb", "Mar", "Apr", "May",
                 "Jun", "Jul","Aug", "Sep", "Oct")
  
  header <- c("DeploymentID","statistic", mon_ab_df) 
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
  
  
  for (i in 1:length(site_meta$DeploymentID)){
    #site_meta <- ml_meta # debugger
    #site_data <- ml #debugger
    #For mL this is the index of the file it seems to get stuck on
    #i <- 1 #* debugger
    #Get the specific sensor ID from the metadata
    sensor<- site_meta$DeploymentID[i]
    #Use the metadata senosr info to select only the values we want from the site data file
    x <- site_data %>%
      filter(DeploymentID %in% sensor)
    
    #The time format should be good to go
    
    t.max <- max(x$temp.C)  #max for entire sensor year
    t.min <- min(x$temp.C)  #min for entire sensor year
    
    #df is the mostly empty data frame; put sensor name into rows under "DeploymentID" col name
    df[1:7,c("DeploymentID")] <- sensor # name the sensor-year ID, we will work on getting the exact sensor year later if needed
    
    
    # make sure all columns of monthly stats are set to "NA" before saving the stats unique to this sensor-year;
    # this is important for all sensor-years after the first one, so we don't have stats carrying over from one
    # sensor-year to the next, which happens when sensor-year t contains data for more months than sensor-year t+1
    df[1:7,3:19] <- NA
    
    dt <- x$date.time
    
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
    if (length(grep("Jun",mo.yrs[1]))) ii <- 3
    if (length(grep("Jul",mo.yrs[1]))) ii <- 4
    if (length(grep("Aug",mo.yrs[1]))) ii <- 5 #* column = 3 only if mo.yrs[1] contains "Aug"
    if (length(grep("Sep",mo.yrs[1]))) ii <- 6 #* column = 4 only if mo.yrs[1] contains "Sep", etc.
    if (length(grep("Oct",mo.yrs[1]))) ii <- 7
    if (length(grep("Nov",mo.yrs[1]))) ii <- 8
    if (length(grep("Dec",mo.yrs[1]))) ii <- 9
    #* Sensors were never placed later than early Nov, so we don't have to worry about Dec, etc;
    #Except for the Aiden "relic" sensor haha
    
    #* now position all mo.yrs within the 1st row of df, starting at column ii
    df[1,ii:(ii+length(mo.yrs)-1)] <- mo.yrs
    #df #* looks to see if any bugs
    
    #* now loop through the temperature data by month-year
    
    for (my in mo.yrs){
      
      #my <- mo.yrs[2]; my #* debugger
      
      #j indexes all the positions in mo.yr that "my" currently represents (in the current iteration of the loop)
      j <- which(mo.yr==my)
      
      #* report whether this "my" has enough data to consider from this sensor-year
      data.days <- length(unique(dt$mday[j]))
      ii <- which(df[1,]==my); ii #* ii indexes the position of "my" in the top ROW (not header) of df
      if (data.days < 20) df[2:7,ii] <- "low data"
      
      #df #look to see if any bugs
      
      #next means to skip the current iteration of a loop without terminating it.
      #So, if a month has less than 20 days in it, we are not going to use it in our analyses.
      if (data.days < 20) next
      
      xv <- dt[j]$hour #* xv is now the hour from each date-time within all data.days of this month-year
      
      ys <- x$temp.C[j] #* ys is now the temperature "
      #* for a 31-day month, both xv and ys will each be 31*k long, where k = # of temps logged daily
      
      k <- min(6,length(unique(xv))) #* k is # of temperatures logged daily
      
      stats <- daily.fit.stats(sensorfile=sensor,xv=dt[j],ys=x$temp.C[j],month=my,knots=k, min_knots = min_knots)
      
      
      ii <- which(df[1,]==my); ii #* ii indexes the position of "my" in the top ROW (not header) of df
      df[2:7,ii] <- as.character(stats$value)
      df
      
    }
    suppressWarnings(write.csv(df,csv_stats_title, TRUE))
    
  }
}


##------------------------------------------------------------------------
##------------------------------------------------------------------------

#* Run function on Krear site data
sur_site_stats_ftn(krear_site_metadata, krear_site, "Krear_site_recent_shallow_temp-stats-table.csv", 6)
#With a minimum of 6 knots there are still a couple sensors that have fewer observations per day.
#For lack of a better solution, they were removed from the metadata file.
#Also, the data were trimmed so every day had a minimum of 6 observations 
## (meaning, if a sensor was pulled early in the morning before it could record 6 times that day was removed)


#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------

#Note that in this instance, Krear Site refers only to K1.
files <- read.csv("Krear_Site_Sensor_Metadata.csv",stringsAsFactors=F)

#Column one has the names of files, name in file is filename.new
filenames <- files[,1]
#column with the info about if the times are in 24 or 12 hour format.
time.format <- files[,2]
#column with info about sensor position
s.position <- files[,3]

#This is the actual data for K1 with all stats to match Krear
data <- read.csv("revised_temp-stats-table_daily_fits.csv") 


#filter out rows
df <- filter(data, !(statistic %in% c("Month-Year", "statistic")))

#Make long
df <- pivot_longer(df, Aug:length(df), names_to = "Month", values_to = "values")



#Now load in surrounding site data frames:
krear_sur_stats <- read.csv("Krear_site_recent_shallow_temp-stats-table.csv")
krear_sur_stats <- filter(krear_sur_stats, !(statistic %in% c("Month-Year", "statistic")))
krear_sur_stats <- pivot_longer(krear_sur_stats, Jun:length(krear_sur_stats), names_to = "Month", values_to = "values")
head(krear_sur_stats)


#------------------------------------------------------------------------------


graph_jul2Jun_f <- function(s.location, sur.location, files, df, sur_meta, sur_data, stats_to_graph, graph_title ,stat_labels, y.lab, y.axis){
  
  #Get the sensor year names from the files df for each location type
  loc_sensors <- files %>%
    filter(sensor.type %in% c(s.location)) %>%
    pull(filename.new)
  
  sur_sensors <- sur_meta %>%
    filter(Depth %in% c(sur.location)) %>%
    pull(DeploymentID)
  
  #now get only the data for the sensors at the location
  loc_recent_df <- df %>%
    filter(sensor.year %in% loc_sensors)
  
  sur_df <- sur_data %>%
    filter(DeploymentID %in% sur_sensors)
  
  
  #Deal with the duplicate months. 
  loc_recent_df$Month <- gsub("\\..*","",loc_recent_df$Month)
  sur_df$Month <- gsub("\\..*","",sur_df$Month)
  
  loc_recent_df <- loc_recent_df %>% 
    replace_with_na(replace = list(values = "low data")) %>%
    mutate(values = as.numeric(values))
  
  sur_df <- sur_df %>% 
    replace_with_na(replace = list(values = "low data")) %>%
    mutate(values = as.numeric(values))
  
  graph.data.krear <- subset(loc_recent_df, statistic %in% stats_to_graph)
  graph.data.sur <- subset(sur_df, statistic %in% stats_to_graph)
  
  #Re-order the factor levels
  graph.data.krear$Month <- factor(graph.data.krear$Month, levels = c("Jul","Aug", "Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun"))
  graph.data.sur$Month <- factor(graph.data.sur$Month, levels = c("Jul","Aug", "Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun"))
  
  #convert the months  to numeric
  graph.data.krear$Month.num <- match(graph.data.krear$Month,month.abb) 
  graph.data.sur$Month.num <- match(graph.data.sur$Month, month.abb)
  
  #now turn the months into date class, using pretend dates to get around the issue.
  graph.data.krear$pretend_year <-  ifelse(graph.data.krear$Month.num > 6, "2009", "2010")
  graph.data.krear$Date <- as.yearmon(paste(graph.data.krear$pretend_year, graph.data.krear$Month.num), "%Y %m")
  graph.data.krear$Date <- as.Date(graph.data.krear$Date)
  
  graph.data.sur$pretend_year <-  ifelse(graph.data.sur$Month.num > 6, "2009", "2010")
  graph.data.sur$Date <- as.yearmon(paste(graph.data.sur$pretend_year, graph.data.sur$Month.num), "%Y %m")
  graph.data.sur$Date <- as.Date(graph.data.sur$Date)
  
  if(stats_to_graph == "Days above freezing") {
    
    g_poisson <- ggplot(graph.data.sur, 
                        aes(x = Date, y = values)) +
      geom_point(aes(color = statistic, shape = statistic)) +
      theme_minimal() +
      #use this geom to do the gam calculation within the ggplot code. Plotted line is exactly the same as outside/normal calculation.
      stat_smooth(data = (subset(graph.data.sur, statistic %in% stats_to_graph[1])), aes(x=Date, y = values), 
                  method="gam",   method.args = list(family = "poisson", link="log"),
                  inherit.aes = F, color = "gray40", fill = "gray30", alpha = .35) +
      
      labs(y = y.lab, title = graph_title) +
      #modify the color legend for the recent years stats
      scale_color_discrete(name = "Krear Site Recent", labels = stat_labels)+ 
      #change the x axis labels from numeric to months
      scale_x_date(breaks = graph.data.sur$Date, date_labels="%b") +
      ylim(y.axis) +
      theme(axis.title.x = element_blank())+
      
      #now add the Krear things to the graph:
      geom_point(data = graph.data.krear, 
                 aes(x = Date, y = values, shape = statistic), color = "#F8766D") +
      stat_smooth(data = (subset(graph.data.krear, statistic %in% stats_to_graph[1])), aes(x=Date, y = values), 
                  method="gam",  method.args = list(family = "poisson", link="log"),
                  inherit.aes = F, color = "#F8766D", fill = "#F8766D", alpha = .35) +
      #modify the shape legend to only have historical labels in it.
      scale_shape_manual(values=c(16, 16), name = "K1 Recent", breaks = stats_to_graph, labels = stat_labels) +
      scale_color_manual(values=c("gray", "#F8766D"), name = "Krear Site", labels = stat_labels)
    
    
    print(g_poisson)
    
    
  }
  
  else{
    
    if (length(stats_to_graph) > 1.1){
      
      
      graph <- ggplot(graph.data.sur, 
                      aes(x = Date, y = values)) +
        geom_point(aes(color = statistic, shape = statistic), color = "gray40") +
        theme_minimal() +
        #use this geom to do the gam calculation within the ggplot code. Plotted line is exactly the same as outside/normal calculation.
        stat_smooth(data = (subset(graph.data.sur, statistic %in% stats_to_graph[1])), aes(x=Date, y = values), method="gam",  formula= y~s(x), inherit.aes = F, color = "gray40", fill = "gray30", alpha = .35) +
        stat_smooth(data = (subset(graph.data.sur, statistic %in% stats_to_graph[2])), aes(x=Date, y = values), method="gam",  formula= y~s(x), inherit.aes = F, color = "gray40", fill = "gray30", alpha = .35) +
        
        labs(y = y.lab, title = graph_title) +
        #modify the color legend for the recent years stats
        #scale_color_discrete(name = "Krear Site Recent", labels = stat_labels)+ 
        #change the x axis labels from numeric to months
        scale_x_date(breaks = graph.data.sur$Date, date_labels="%b") +
        ylim(y.axis) +
        theme(axis.title.x = element_blank())+
        
        #now add the Krear things to the graph:
        geom_point(data = graph.data.krear, 
                   aes(x = Date, y = values, shape = statistic, color = statistic)) +
        stat_smooth(data = (subset(graph.data.krear, statistic %in% stats_to_graph[1])), aes(x=Date, y = values), method="gam",  formula= y~s(x), inherit.aes = F, color = "#F8766D", fill = "#F8766D", alpha = .35) +
        stat_smooth(data = (subset(graph.data.krear, statistic %in% stats_to_graph[2])), aes(x=Date, y = values), method="gam",  formula= y~s(x), inherit.aes = F, color = "#00BFC4", fill = "#00BFC4", alpha = .35) +
        
        #modify the shape legend to only have historical labels in it.
        scale_shape_manual(values=c(16, 17, 16, 17), name = "K1 Recent", breaks = stats_to_graph, labels = stat_labels) +
        scale_color_manual(values = c("#F8766D", "#00BFC4", "gray","gray"), name =  "Krear Site", breaks = stats_to_graph, labels = stat_labels)
      
      
      print(graph)} else{
        
        
        
        graph <- ggplot(graph.data.sur, 
                        aes(x = Date, y = values)) +
          geom_point(aes(color = statistic, shape = statistic)) +
          theme_minimal() +
          #use this geom to do the gam calculation within the ggplot code. Plotted line is exactly the same as outside/normal calculation.
          stat_smooth(data = (subset(graph.data.sur, statistic %in% stats_to_graph[1])), aes(x=Date, y = values), method="gam",  formula= y~s(x), inherit.aes = F, color = "gray40", fill = "gray30", alpha = .35) +
          
          labs(y = y.lab, title = graph_title) +
          #modify the color legend for the recent years stats
          scale_color_discrete(name = "Krear Site Recent", labels = stat_labels)+ 
          #change the x axis labels from numeric to months
          scale_x_date(breaks = graph.data.sur$Date, date_labels="%b") +
          ylim(y.axis) +
          theme(axis.title.x = element_blank())+
          
          #now add the Krear things to the graph:
          geom_point(data = graph.data.krear, 
                     aes(x = Date, y = values, shape = statistic), color = "#F8766D") +
          stat_smooth(data = (subset(graph.data.krear, statistic %in% stats_to_graph[1])), aes(x=Date, y = values), method="gam",  formula= y~s(x), inherit.aes = F, color = "#F8766D", fill = "#F8766D", alpha = .35) +
          #modify the shape legend to only have historical labels in it.
          scale_shape_manual(values=c(16, 16), name = "K1 Recent", breaks = stats_to_graph, labels = stat_labels) +
          scale_color_manual(values=c("gray", "#F8766D"), name = "Krear Site", labels = stat_labels)
        
        print(graph)
      }
  }
  
}

#----------------------------------------------------------------------------------------------------------

dev.off()


pdf("revised_Krear Site Haypiles vs K1 Comparison Plots GAM All.pdf")

#function input order:
#function(s.location, sur.location, files, df, sur_meta, sur_data, stats_to_graph, graph_title ,stat_labels, y.lab, y.axis){

#Max and Min stats
graph_shallow_month_maxmin <- graph_jul2Jun_f("shallow", "Std", files, df, krear_site_metadata, krear_sur_stats, c("Month Maximum", "Month Minimum"),  "Krear Site vs K1,\nMonthly Max and Min",c("Monthly Maximum", "Monthly Minimum"), "Monthly temperature (C)", c(-28,31))

#Mean Daily Max and Min stats
graph_shallow_daily2_maxmin <- graph_jul2Jun_f("shallow", "Std", files, df, krear_site_metadata, krear_sur_stats, c("Mean daily max T", "Mean daily min T"), "Krear Site vs K1,\nMean Daily Max and Min",c("Mean Daily Maximum", "Mean Daily Minimum"), "Monthly mean of daily temperature (C)", c(-28,31))

#Mean of Daily Max and Min stats
graph_shallow_mean_max_min <- graph_jul2Jun_f("shallow", "Std", files, df, krear_site_metadata, krear_sur_stats, "Mean of Daily Max and Min",  "Krear Site vs K1,\nMean of Daily Max and Min", "Mean of Daily\nMax and Min", "Monthly mean of daily temperature (C)", c(-28,31))

#dev.off()

#pdf("revised_Days Abv F Only Krear Site Haypiles vs K1 Comparison Plots Fixed Colors.pdf", width=5.5, height=7)
#Days above Freezing stat
graph_air_days_abv_0 <- graph_jul2Jun_f("shallow", "Std", files, df, krear_site_metadata, krear_sur_stats, "Days above freezing",  "Krear Site vs K1,\nNumber of Days Above Freezing", "Days above Freezing", "Number of Days Above Freezing", c(-3,40))

dev.off()
