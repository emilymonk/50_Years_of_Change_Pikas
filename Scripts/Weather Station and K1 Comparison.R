#script to compare long form  data from surrounding sites and weather stations with krear metrics
#This basically modifies the Krear comparison code to work with this different data format

setwd("D:/Pika_Data")

library(tidyverse)
require(mgcv) #for fitting a gam
require(zoo) #for dealing with dates (as.Date)

#Lets start with the weather stations:

#read in the data:
c1 <- read.csv("c1_temp_1952-2020_draft.csv")
d1 <- read.csv("d1_temp_1952-2020_draft.csv")
head(c1)
str(c1)

#These files already have the daily max and mins listed.
#This code just calculates the monthly stats from the daily values


wthr_data_ftn <- function(data, hist_lower, hist_upper, recent_lower, recent_upper){

  wthr_data <- data
  
  #make date column into a date class
  wthr_data$date <- as.Date(wthr_data$date)
  
  #now subset the data by date range
  historical_subset <- wthr_data %>%
    select(local_site, year, date, max_temp, min_temp, mean_temp) %>%
    filter(date >= hist_lower & date <= hist_upper)
  
  recent_subset <- wthr_data %>%
    select(local_site, year, date, max_temp, min_temp, mean_temp) %>%
    filter(date >= recent_lower & date <= recent_upper)
  
  #get all the stats into a long form table
  recent_subset_stats <- recent_subset %>%
    group_by(month = month(date), year) %>%
    summarize(Days_Above_Freezing = sum(max_temp>0,na.rm=T), 
              Month_Max = max(max_temp),
              Month_Min = min(min_temp),
              Mean_Daily_Max = mean(max_temp),
              Mean_Daily_Min = mean(min_temp),
              Mean_of_Daily_Max_and_Min = mean(mean(max_temp), mean(min_temp))) %>%
    ungroup() %>%
    pivot_longer(cols = Days_Above_Freezing:Mean_of_Daily_Max_and_Min, names_to = "statistic", values_to = "values") %>%
    mutate(time = "Recent Years")
  
  historical_subset_stats <- historical_subset %>%
    group_by(month = month(date), year) %>%
    summarize(Days_Above_Freezing = sum(max_temp>0,na.rm=T), 
              Month_Max = max(max_temp),
              Month_Min = min(min_temp),
              Mean_Daily_Max = mean(max_temp),
              Mean_Daily_Min = mean(min_temp),
              Mean_of_Daily_Max_and_Min = mean(mean(max_temp), mean(min_temp))) %>%
    ungroup() %>%
    pivot_longer(cols = Days_Above_Freezing:Mean_of_Daily_Max_and_Min, names_to = "statistic", values_to = "values") %>%
    mutate(time = "Historical Years")
  
  wthr_stats <- rbind(recent_subset_stats, historical_subset_stats)
  head(wthr_stats)
  tail(wthr_stats)
  
  return(wthr_stats)
}
  
###################################################################################################

#This code also considers the data from the K1.

#This code is currently just copied and pasted from the Master Krear Comparison Script,
#For the sake of length I removed the comments and extra spaces. 
#Go to the above mentioned script for more readable code.

files <- read.csv("Krear_Site_Sensor_Metadata.csv",stringsAsFactors=F)
filenames <- files[,1]
time.format <- files[,2]
s.position <- files[,3]
data <- read.csv("temp-stats-table_daily_fits.csv") 
df <- filter(data, !(statistic %in% c("Month-Year", "statistic")))
df <- pivot_longer(df, Aug:length(df), names_to = "Month", values_to = "values")
krear <- read_excel("Krear OG Data Formatted.xlsx", sheet = 3)
krear <- pivot_longer(krear, Jul63:length(krear), names_to = "Month", values_to = "values")
krear$Month <- gsub("6.","",krear$Month)
krear <- krear %>%
  mutate(values = as.numeric(values)) %>%
  mutate(Month = as.factor(Month))
krear$values[krear$statistic != "Days Above Freezing"] <- fahrenheit.to.celsius(krear$values[krear$statistic != "Days Above Freezing"], 2)
krear_air <- krear %>%
  filter(location == "free air")%>%
  select(-location) %>%
  mutate(time = "Historical")


recent_air_sensors <- files %>%
    filter(sensor.type == "free air") %>%
    pull(filename.new)
recent_air_df <- df %>%
    filter(sensor.year %in% recent_air_sensors)
recent_air_df$Month <- gsub("\\..*","",recent_air_df$Month)
recent_air_df <- recent_air_df %>% 
    replace_with_na(replace = list(values = "low data")) %>%
    mutate(values = as.numeric(values)) %>%
  select(-X) %>%
  mutate(time = "Recent")
colnames(recent_air_df)
colnames(krear_air)

krear_site_data <- rbind(recent_air_df, krear_air)

###################################################################################################  
  
#now graph all the stats like I did for the Krear comparison, but add GAMs to the historical values
#The weather station only graphs have distinct colors from comparisons with the K1.
  
graph_jul2Jun_f <- function(wthr_stats, stats_to_graph, graph_title ,stat_labels, y.lab, y.axis, historical_leg, recent_leg){
    
    #subset data to the desired statistic
    graph.data <- subset(wthr_stats, statistic %in% stats_to_graph)
    
    
    #Now change the month column to addreviations then to factor
    graph.data$month <- month.abb[graph.data$month]
    
    #Re-order the factor levels
    graph.data$month <- factor(graph.data$month, levels = c("Jul","Aug", "Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun"))
    
    
    #convert the months  to numeric
    graph.data$month.num <- match(graph.data$month,month.abb) 
    
    #now turn the months into date class, using pretend dates to get around the issue.
    graph.data$pretend_year <-  ifelse(graph.data$month.num > 6, "2009", "2010")
    graph.data$pretend_Date <- as.yearmon(paste(graph.data$pretend_year, graph.data$month.num), "%Y %m")
    graph.data$pretend_Date <- as.Date(graph.data$pretend_Date)
    
    if(stats_to_graph == "Days_Above_Freezing") {
      g_poisson <- ggplot(data = subset(graph.data, time == "Recent Years"), 
             aes(x = pretend_Date, y = values)) +
        geom_point(aes(color = statistic, shape = statistic)) +
        theme_minimal() +
        #use this geom to do the gam calculation within the ggplot code. Plotted line is exactly the same as outside/normal calculation.
        stat_smooth(data = (subset(graph.data, statistic %in% stats_to_graph[1]  & time == "Recent Years")), aes(x=pretend_Date, y = values), 
                    method="gam",   method.args = list(family = "poisson", link="log"),
                    inherit.aes = F, color = "#7CAE00", fill = "#7CAE00", alpha = .35) +
        
        labs(y = y.lab, title = graph_title) +
        #modify the color legend for the recent years stats
        scale_color_discrete(name = recent_leg, labels = stat_labels)+ 
        #change the x axis labels from numeric to months
        scale_x_date(breaks = graph.data$pretend_Date, date_labels="%b") +
        ylim(y.axis) +
        theme(axis.title.x = element_blank(), axis.text = element_text(size = 12), axis.title.y = element_text(size = 14), 
              plot.title = element_text(size = 16), legend.title = element_text(size = 14), legend.text = element_text(size = 12))+
        
        #now add the Krear things to the graph:
        geom_point(data = subset(graph.data, time == "Historical Years"), 
                   aes(x = pretend_Date, y = values, shape = statistic), color = "gray25") +
        stat_smooth(data = (subset(graph.data, statistic %in% stats_to_graph[1] & time == "Historical Years")), aes(x=pretend_Date, y = values), 
                    method="gam",   method.args = list(family = "poisson", link="log"),
                    inherit.aes = F, color = "gray", fill = "gray", alpha = .3) +
        #modify the shape legend to only have historical labels in it.
        scale_shape_manual(values=c(16, 16), name = historical_leg, breaks = stats_to_graph, labels = stat_labels) +
        scale_color_manual(values=c("#7CAE00", "gray25"), name = recent_leg, breaks = stats_to_graph, labels = stat_labels)
      
      print(g_poisson)
     
      
    }
    
    else{
      
      if (length(stats_to_graph) > 1.1){
        
        
        graph <- ggplot(data = subset(graph.data, time == "Recent Years"), 
                        aes(x = pretend_Date, y = values)) +
          geom_point(aes(color = statistic, shape = statistic)) +
          theme_minimal() +
          #use this geom to do the gam calculation within the ggplot code. Plotted line is exactly the same as outside/normal calculation.
          stat_smooth(data = (subset(graph.data, statistic %in% stats_to_graph[1] & time == "Recent Years")), aes(x=pretend_Date, y = values), method="gam",  formula= y~s(x), inherit.aes = F, color = "#7CAE00", fill = "#7CAE00", alpha = .35) +
          stat_smooth(data = (subset(graph.data, statistic %in% stats_to_graph[2] & time == "Recent Years")), aes(x=pretend_Date, y = values), method="gam",  formula= y~s(x), inherit.aes = F, color = "#C77CFF", fill = "#C77CFF", alpha = .35) +
          
          labs(y = y.lab, title = graph_title) +
          #modify the color legend for the recent years stats
          scale_color_discrete(name = recent_leg, labels = stat_labels)+ 
          #change the x axis labels from numeric to months
          #scale_x_continuous(breaks = (1:12), labels = function(x) month.abb[x])+
          scale_x_date(breaks = graph.data$pretend_Date, date_labels="%b") +
          ylim(y.axis) +
          theme(axis.title.x = element_blank(), axis.text = element_text(size = 12), axis.title.y = element_text(size = 14), 
                plot.title = element_text(size = 16), legend.title = element_text(size = 14), legend.text = element_text(size = 12))+
          
          #now add the Krear things to the graph:
          geom_point(data = subset(graph.data, time == "Historical Years"), 
                     aes(x = pretend_Date, y = values, shape = statistic), color = "gray25") +
          stat_smooth(data = (subset(graph.data, statistic %in% stats_to_graph[1] & time == "Historical Years")), aes(x=pretend_Date, y = values), method="gam",  formula= y~s(x), inherit.aes = F, color = "gray", fill = "gray", alpha = .3) +
          stat_smooth(data = (subset(graph.data, statistic %in% stats_to_graph[2]& time == "Historical Years")), aes(x=pretend_Date, y = values), method="gam",  formula= y~s(x), inherit.aes = F, color = "gray", fill = "gray", alpha = .3) +
          
          #modify the shape legend to only have historical labels in it.
          scale_shape_manual(values=c(16, 17, 16, 17), name = historical_leg, breaks = stats_to_graph, labels = stat_labels) +
          scale_color_manual(values = c("#7CAE00", "#C77CFF", "#7CAE00",  "#C77CFF"), name = recent_leg, breaks = stats_to_graph, labels = stat_labels)
        
        
        print(graph)} else{
          
          
          
          graph <- ggplot(data = subset(graph.data, time == "Recent Years"), 
                          aes(x = pretend_Date, y = values)) +
            geom_point(aes(color = statistic, shape = statistic)) +
            theme_minimal() +
            #use this geom to do the gam calculation within the ggplot code. Plotted line is exactly the same as outside/normal calculation.
            stat_smooth(data = (subset(graph.data, statistic %in% stats_to_graph[1]  & time == "Recent Years")), aes(x=pretend_Date, y = values), method="gam",  formula= y~s(x), inherit.aes = F, color = "#7CAE00", fill = "#7CAE00", alpha = .35) +
            
            labs(y = y.lab, title = graph_title) +
            #modify the color legend for the recent years stats
            scale_color_discrete(name = recent_leg, labels = stat_labels)+ 
            #change the x axis labels from numeric to months
            scale_x_date(breaks = graph.data$pretend_Date, date_labels="%b") +
            ylim(y.axis) +
            theme(axis.title.x = element_blank(), axis.text = element_text(size = 12), axis.title.y = element_text(size = 14), 
                  plot.title = element_text(size = 16), legend.title = element_text(size = 14), legend.text = element_text(size = 12))+
            
            #now add the Krear things to the graph:
            geom_point(data = subset(graph.data, time == "Historical Years"), 
                       aes(x = pretend_Date, y = values, shape = statistic), color = "gray25") +
            stat_smooth(data = (subset(graph.data, statistic %in% stats_to_graph[1] & time == "Historical Years")), aes(x=pretend_Date, y = values), method="gam",  formula= y~s(x), inherit.aes = F, color = "gray", fill = "gray", alpha = .3) +
            #modify the shape legend to only have historical labels in it.
            scale_shape_manual(values=c(16, 16), name = historical_leg, breaks = stats_to_graph, labels = stat_labels) +
            scale_color_manual(values=c("#7CAE00", "gray25"), name = recent_leg, breaks = stats_to_graph, labels = stat_labels)
          
          print(graph)
        }
    }
    
}

c1_wthr_stats <- wthr_data_ftn(c1, "1958-07-01", "1969-06-30", "2009-07-01", "2020-06-30")
d1_wthr_stats <- wthr_data_ftn(d1, "1958-07-01", "1969-06-30", "2009-07-01", "2020-06-30")


dev.off()
pdf("NWT Ambient Temperature Data Comparison Plots GAM All.pdf")
  
  #function input order:
  #function(stats_to_graph, graph_title ,stat_labels, y.lab, y.axis, recent_leg, historical_leg)
  
  #Max and Min stats
  C1_month_maxmin <- graph_jul2Jun_f(c1_wthr_stats, c("Month_Max", "Month_Min"), "C1 Weather Station, Monthly Max and Min",c("Month Maximum", "Monthly Minimum"), "Monthly temperature (C)", c(-40,35), "C1, 1958-1969", "C1, 2009-2020")
  D1_month_maxmin <- graph_jul2Jun_f(d1_wthr_stats, c("Month_Max", "Month_Min"), "D1 Weather Station, Monthly Max and Min",c("Month Maximum", "Monthly Minimum"), "Monthly temperature (C)", c(-40,35), "D1, 1958-1969", "D1, 2009-2020")
  
  #Mean Daily Max and Min stats
  C1_daily_maxmin <- graph_jul2Jun_f(c1_wthr_stats, c("Mean_Daily_Max", "Mean_Daily_Min"), "C1 Weather Station, Mean Daily Max and Min",c("Mean Daily Maximum", "Mean Daily Minimum"), "Monthly temperature (C)", c(-40,35), "C1, 1958-1969", "C1, 2009-2020")
  D1_daily_maxmin <- graph_jul2Jun_f(d1_wthr_stats, c("Mean_Daily_Max", "Mean_Daily_Min"), "D1 Weather Station, Mean Daily Max and Min",c("Mean Daily Maximum", "Mean Daily Minimum"), "Monthly temperature (C)", c(-40,35), "D1, 1958-1969", "D1, 2009-2020")
  
  #Mean of Daily Max and Min stats 
  C1_mean_max_min <- graph_jul2Jun_f(c1_wthr_stats, "Mean_of_Daily_Max_and_Min",  "C1 Weather Station, Mean of Daily Max and Min", "Mean of Daily\nMax and Min", "Monthly mean of daily temperature (C)", c(-40,35), "C1, 1958-1969", "C1, 2009-2020")
  D1_mean_max_min <- graph_jul2Jun_f(d1_wthr_stats, "Mean_of_Daily_Max_and_Min",  "D1 Weather Station, Mean of Daily Max and Min", "Mean of Daily\nMax and Min", "Monthly mean of daily temperature (C)", c(-40,35), "D1, 1958-1969", "D1, 2009-2020")
  
  #Days above Freezing stat
  C1_days_abv_0 <- graph_jul2Jun_f(c1_wthr_stats, "Days_Above_Freezing", "C1 Weather Station, Number of Days Above Freezing",  "C1, Number of Days\nAbove Freezing", "Number of Days Above Freezing", c(-3,40), "NA", c("C1, 2009-2020", "C1, 1958-1969"))
  D1_days_abv_0 <- graph_jul2Jun_f(d1_wthr_stats, "Days_Above_Freezing", "D1 Weather Station, Number of Days Above Freezing", "D1, Number of Days\nAbove Freezing", "Number of Days Above Freezing", c(-3,40), "NA", c("D1, 2009-2020", "D1, 1958-1969"))
  #                             function(wthr_stats, stats_to_graph,      graph_title ,                                         stat_labels,                          y.lab,                          y.axis, historical_leg, recent_leg)
  
dev.off()

  #for reference here are the names of the stats for the weather data I made:
  # Days_Above_Freezing
  # Month_Max
  # Month_Min
  # Mean_Daily_Max
  # Mean_Daily_Min
  # Mean_of_Daily_Max_and_Min

#############################################################################################
  
#Now Write code to compare the weather data above to the Krear data:
# 2) Historical Krear (points) vs Historical C1 (GAM)
# 3) Historical Krear (points) vs Historical D1 (GAM)
# 4) Recent Krear (GAM) vs Recent C1 (GAM)
# 5) Recent Krear (GAM) vs Recent D1 (GAM)
  
#This code will end up being quite similar to the above, but different enough to warrant a new function

#for reference, here are the names of the statistics for the K1 (urgh, why did I make them all different?)
#Recent
# "Days above freezing"
# "Mean daily max T"
# "Mean daily min T"         
# "Month Maximum"
# "Month Minimum"
# "Mean of Daily Max and Min"

#Historical
# "Max.T"
# "Min.T"
# "Mean.Daily.Max.T"
# "Mean.Daily.Min.T"
# "Mean.Daily.MaxMin.T"
# "Days Above Freezing" 

#testing variables:
# stats_to_graph <- Mean_of_Daily_Max_and_Min
# stats_to_graph <- c("Month_Max", "Month_Min")
y.lab <- "number of days"
graph_title <- "D1 Historical"
 y.axis <- c(-4,45)
# stat_labels <- c("Month Maximum", "Month Minimum")
leg_shape <- "K1, 1963-1964"
leg_color <- "D1, 1958-1969"
time.var = "Historical Years"
wthr_stats <- d1_wthr_stats
stats_to_graph <- "Days_Above_Freezing"
stat_labels <- "Days Above Freezing"

krear_stats <- "Days Above Freezing"
krear_stats <- c("Max.T","Min.T")
krear_stats <- c("Month Maximum", "Month Minimum")



#function(wthr_stats, stats_to_graph, graph_title ,stat_labels, y.lab, y.axis, historical_leg, recent_leg){

krear_wthr_stn_f <- function(wthr_stats, time.var, krear_site_data, stats_to_graph, krear_stats, graph_title ,stat_labels, y.lab, y.axis, leg_shape, leg_color){

  krear_df <- subset(krear_site_data, statistic %in% krear_stats)

  #convert the months  to numeric
  krear_df$Month.num <- match(krear_df$Month,month.abb) 

  #now turn the months into date class, using pretend dates to get around the issue.
  krear_df$pretend_year <-  ifelse(krear_df$Month.num > 6, "2009", "2010")
  krear_df$pretend_Date <- as.yearmon(paste(krear_df$pretend_year, krear_df$Month.num), "%Y %m")
  krear_df$pretend_Date <- as.Date(krear_df$pretend_Date)

  #subset data to the desired statistic
  graph.data <- subset(wthr_stats, statistic %in% stats_to_graph)
  
  #Now change the month column to addreviations then to factor
  graph.data$month <- month.abb[graph.data$month]
  
  #convert the months  to numeric
  graph.data$month.num <- match(graph.data$month,month.abb) 
  
  
  #now turn the months into date class, using pretend dates to get around the issue.
  graph.data$pretend_year <-  ifelse(graph.data$month.num > 6, "2009", "2010")
  graph.data$pretend_Date <- as.yearmon(paste(graph.data$pretend_year, graph.data$month.num), "%Y %m")
  graph.data$pretend_Date <- as.Date(graph.data$pretend_Date)
  
  #Max and min graphs:
  if (length(stats_to_graph) > 1.1){

  graph <- ggplot(data = subset(graph.data, time == time.var), 
                  aes(x = pretend_Date, y = values)) +
    geom_point(aes(color = statistic, shape = statistic)) +
    theme_minimal() +
    
    labs(y = y.lab, title = graph_title) +
    #modify the color legend for the recent years stats
    #scale_color_discrete(name = leg_color, labels = stat_labels)+ 
    #change the x axis labels from numeric to months
    scale_x_date(breaks = graph.data$pretend_Date, date_labels="%b") +
    ylim(y.axis) +
    theme(axis.title.x = element_blank(), axis.text = element_text(size = 12), axis.title.y = element_text(size = 14), 
          plot.title = element_text(size = 16), legend.title = element_text(size = 14), legend.text = element_text(size = 12))
  
  if(krear_df$time[1] == "Historical"){
    graph <- graph +
      #use this geom to do the gam calculation within the ggplot code. Plotted line is exactly the same as outside/normal calculation.
      #this plots the historical weather station data
      stat_smooth(data = (subset(graph.data, statistic %in% stats_to_graph[1] & time == "Historical Years")), aes(x=pretend_Date, y = values), method="gam",  formula= y~s(x), inherit.aes = F, color = "gray", fill = "gray", alpha = .4) +
      stat_smooth(data = (subset(graph.data, statistic %in% stats_to_graph[2] & time == "Historical Years")), aes(x=pretend_Date, y = values), method="gam",  formula= y~s(x), inherit.aes = F, color = "gray", fill = "gray", alpha = .4) +
      
      #add K1 data
      geom_point(data = krear_df, 
                 aes(x = pretend_Date, y = values, shape = statistic, color = statistic)) +
      geom_line(data = krear_df, 
                aes(x = pretend_Date, y = values, group = statistic)) +
      scale_shape_manual(values=c(16, 17, 16, 17), name = leg_shape, breaks = stats_to_graph, labels = stat_labels) +
      scale_color_manual(values = c("black", "black", "gray", "gray"), name = leg_color, breaks = stats_to_graph, labels = stat_labels)
    print(graph)
    }
  else{
    #now add the Krear things to the graph:
   graph <- graph + geom_point(data = krear_df, 
               aes(x = pretend_Date, y = values, shape = statistic, color = statistic)) +
     
     #use this geom to do the gam calculation within the ggplot code. Plotted line is exactly the same as outside/normal calculation.
     #now plot the recent weather station data
     stat_smooth(data = (subset(graph.data, statistic %in% stats_to_graph[1] & time == "Recent Years")), aes(x=pretend_Date, y = values), method="gam",  formula= y~s(x), inherit.aes = F, color = "gray", fill = "gray", alpha = .4) +
     stat_smooth(data = (subset(graph.data, statistic %in% stats_to_graph[2] & time == "Recent Years")), aes(x=pretend_Date, y = values), method="gam",  formula= y~s(x), inherit.aes = F, color = "gray", fill = "gray", alpha = .4) +
     #add K1 data
      stat_smooth(data = (subset(krear_df, statistic %in% krear_stats[1])), aes(x=pretend_Date, y = values), method="gam",  formula= y~s(x), inherit.aes = F, color = "#F8766D", fill = "#F8766D", alpha = .35) +
      stat_smooth(data = (subset(krear_df, statistic %in% krear_stats[2])), aes(x=pretend_Date, y = values), method="gam",  formula= y~s(x), inherit.aes = F, color = "#00BFC4", fill = "#00BFC4", alpha = .35) +
    
      #modify the shape legend to only have historical labels in it.
     scale_shape_manual(values=c(16, 17, 16, 17), name = leg_shape, breaks = stats_to_graph, labels = stat_labels) +
     scale_color_manual(values = c("#F8766D", "#00BFC4", "gray","gray"), name =  leg_color, breaks = stats_to_graph, labels = stat_labels)
    print(graph)
  }
  
  }
  else{
    #days above freezing
    if(stats_to_graph == "Days_Above_Freezing"){
     
      graph_base <- ggplot(data = subset(graph.data, time == time.var), 
                      aes(x = pretend_Date, y = values)) +
        geom_point(aes(color = statistic, shape = statistic)) +
        theme_minimal() +
        
        labs(y = y.lab, title = graph_title) +
        #modify the color legend for the recent years stats
        scale_color_discrete(name = leg_color, labels = stat_labels)+ 
        #change the x axis labels from numeric to months
        scale_x_date(breaks = graph.data$pretend_Date, date_labels="%b") +
        ylim(y.axis) +
        theme(axis.title.x = element_blank(), axis.text = element_text(size = 12), axis.title.y = element_text(size = 14), 
              plot.title = element_text(size = 16), legend.title = element_text(size = 14), legend.text = element_text(size = 12))
      
      if(krear_df$time[1] == "Historical"){
        
        graph_hist <- graph_base +
          #add weather station historical GAM
          stat_smooth(data = (subset(graph.data, statistic %in% stats_to_graph[1] & time == "Historical Years")), aes(x=pretend_Date, y = values), 
                      method="gam", method.args = list(family = "poisson", link="log"),
                      inherit.aes = F, color = "gray", fill = "gray", alpha = .3) +
          
          geom_point(data = krear_df, 
                     aes(x = pretend_Date, y = values, shape = statistic)) +
          geom_line(data = krear_df, 
                    aes(x = pretend_Date, y = values, group = statistic)) +
          scale_shape_manual(values=c(16, 16), name = leg_shape, breaks = stats_to_graph, labels = stat_labels) +
          scale_color_manual(values = c("gray", "black"), name = leg_color, labels = stat_labels)
        print(graph_hist)
      }
      else{
        
        #now add the Krear things to the graph:
        graph_rec <- graph_base + geom_point(data = krear_df, 
                                    aes(x = pretend_Date, y = values, shape = statistic), color = "#F8766D") +
          #use this geom to do the gam calculation within the ggplot code. Plotted line is exactly the same as outside/normal calculation.
          stat_smooth(data = (subset(graph.data, statistic %in% stats_to_graph[1] & time == "Recent Years")), aes(x=pretend_Date, y = values), 
                      method="gam",  method.args = list(family = "poisson", link="log"),
                      inherit.aes = F, color = "gray", fill = "gray", alpha = .3) +
          
          
          stat_smooth(data = (subset(krear_df, statistic %in% krear_stats[1])), aes(x=pretend_Date, y = values), 
                      method="gam",  method.args = list(family = "poisson", link="log"),
                      inherit.aes = F, color = "#F8766D", fill = "#F8766D", alpha = .35) +
          
          #modify the shape legend to only have historical labels in it.
          scale_shape_manual(values=c(16, 16), name = leg_shape, breaks = stats_to_graph, labels = stat_labels)+
          scale_color_manual(values=c("gray", "black"), name = leg_color, labels = stat_labels) 
        print(graph_rec)
        }
        
      
    }
    else{
      #Monthly means
    graph <- ggplot(data = subset(graph.data, time == time.var), 
                    aes(x = pretend_Date, y = values)) +
      geom_point(aes(color = statistic, shape = statistic)) +
      theme_minimal() +
      
      labs(y = y.lab, title = graph_title) +
      #modify the color legend for the recent years stats
      scale_color_discrete(name = leg_color, labels = stat_labels)+ 
      #change the x axis labels from numeric to months
      scale_x_date(breaks = graph.data$pretend_Date, date_labels="%b") +
      ylim(y.axis) +
      theme(axis.title.x = element_blank(), axis.text = element_text(size = 12), axis.title.y = element_text(size = 14), 
            plot.title = element_text(size = 16), legend.title = element_text(size = 14), legend.text = element_text(size = 12))
    
    if(krear_df$time[1] == "Historical"){
      graph <- graph +
        #add weather station historical GAM
        stat_smooth(data = (subset(graph.data, statistic %in% stats_to_graph[1] & time == "Historical Years")), aes(x=pretend_Date, y = values), method="gam",  formula= y~s(x), inherit.aes = F, color = "gray", fill = "gray", alpha = .3) +
        
        geom_point(data = krear_df, 
                   aes(x = pretend_Date, y = values, shape = statistic)) +
        geom_line(data = krear_df, 
                  aes(x = pretend_Date, y = values, group = statistic)) +
        scale_shape_manual(values=c(16, 16), name = leg_shape, breaks = stats_to_graph, labels = stat_labels) +
        scale_color_manual(values = c("gray", "black"), name = leg_color, labels = stat_labels)
      print(graph)
    }
    else{
      #now add the Krear things to the graph:
      graph <- graph + geom_point(data = krear_df, 
                                  aes(x = pretend_Date, y = values, shape = statistic), color = "#F8766D") +
        #use this geom to do the gam calculation within the ggplot code. Plotted line is exactly the same as outside/normal calculation.
        stat_smooth(data = (subset(graph.data, statistic %in% stats_to_graph[1] & time == "Recent Years")), aes(x=pretend_Date, y = values), method="gam",  formula= y~s(x), inherit.aes = F, color = "gray", fill = "gray", alpha = .3) +
        
        
        stat_smooth(data = (subset(krear_df, statistic %in% krear_stats[1])), aes(x=pretend_Date, y = values), method="gam",  formula= y~s(x), inherit.aes = F, color = "#F8766D", fill = "#F8766D", alpha = .35) +
        
        #modify the shape legend to only have historical labels in it.
        scale_shape_manual(values=c(16, 16), name = leg_shape, breaks = stats_to_graph, labels = stat_labels)+
        scale_color_manual(values=c("gray", "black"), name = leg_color, labels = stat_labels) 
      print(graph)
    }
    }
    
  }
    
}


c1_wthr_stats <- wthr_data_ftn(c1, "1958-07-01", "1969-06-30", "2013-07-01", "2020-06-30")
d1_wthr_stats <- wthr_data_ftn(d1, "1958-07-01", "1969-06-30", "2013-07-01", "2020-06-30")


dev.off()
pdf("NWT Ambient Temperature and K1 Data Comparison Plots GAM All.pdf")

#                                       (wthr_stats,  time.var,       krear_site_data,  stats_to_graph,                 krear_stats,                            graph_title ,                                               stat_labels,                              y.lab,          y.axis,       leg_shape, leg_color){
#Monthly Max and Mins
C1_recent_month_maxmin <- krear_wthr_stn_f(c1_wthr_stats, "Recent Years", krear_site_data, c("Month_Max", "Month_Min"), c("Month Maximum", "Month Minimum"),"Recent Monthly Max and Min\nC1 Weather Station and K1",c("Monthly Maximum", "Monthly Minimum"), "Monthly temperature (C)", c(-40,35), "K1, 2013-2020", "C1, 2013-2020")
D1_recent_month_maxmin <- krear_wthr_stn_f(d1_wthr_stats, "Recent Years", krear_site_data, c("Month_Max", "Month_Min"), c("Month Maximum", "Month Minimum"),"Recent Monthly Max and Min\nD1 Weather Station and K1",c("Monthly Maximum", "Monthly Minimum"), "Monthly temperature (C)", c(-40,35), "K1, 2013-2020", "D1, 2013-2020")

C1_historical_month_maxmin <- krear_wthr_stn_f(c1_wthr_stats, "Historical Years", krear_site_data, c("Month_Max", "Month_Min"), c("Max.T","Min.T"),"Historical Monthly Max and Min\nC1 Weather Station and K1",c("Monthly Maximum", "Monthly Minimum"), "Monthly temperature (C)", c(-40,35), "K1, 1963-1964", "C1, 1958-1969")
D1_historical_month_maxmin <- krear_wthr_stn_f(d1_wthr_stats, "Historical Years", krear_site_data, c("Month_Max", "Month_Min"), c("Max.T","Min.T"),"Historical Monthly Max and Min\nD1 Weather Station and K1",c("Monthly Maximum", "Monthly Minimum"), "Monthly temperature (C)", c(-40,35), "K1, 1963-1964", "D1, 1958-1969")

#Mean Daily Max and Mins
C1_recent_month_maxmin <- krear_wthr_stn_f(c1_wthr_stats, "Recent Years", krear_site_data, c("Mean_Daily_Max", "Mean_Daily_Min"), c("Mean daily max T", "Mean daily min T"),"Recent Mean Daily Max and Mean of Daily Min\nC1 Weather Station and K1",c("Mean Daily Maximum", "Mean Daily Minimum"), "Monthly mean of daily temperature (C)", c(-40,35), "K1, 2013-2020", "C1, 2013-2020")
D1_recent_month_maxmin <- krear_wthr_stn_f(d1_wthr_stats, "Recent Years", krear_site_data, c("Mean_Daily_Max", "Mean_Daily_Min"), c("Mean daily max T", "Mean daily min T"),"Recent Mean Daily Max and Mean of Daily Min\nD1 Weather Station and K1",c("Mean Daily Maximum", "Mean Daily Minimum"), "Monthly mean of daily temperature (C)", c(-40,35), "K1, 2013-2020", "D1, 2013-2020")

C1_historical_month_maxmin <- krear_wthr_stn_f(c1_wthr_stats, "Historical Years", krear_site_data, c("Mean_Daily_Max", "Mean_Daily_Min"), c("Mean.Daily.Max.T","Mean.Daily.Min.T"),"Historical Mean Daily Max and Mean of Daily Min\nC1 Weather Station and K1",c("Mean Daily Maximum", "Mean Daily Minimum"), "Monthly mean of daily temperature (C)", c(-40,35), "K1, 1963-1964", "C1, 1958-1969")
D1_historical_month_maxmin <- krear_wthr_stn_f(d1_wthr_stats, "Historical Years", krear_site_data, c("Mean_Daily_Max", "Mean_Daily_Min"), c("Mean.Daily.Max.T","Mean.Daily.Min.T"),"Historical Mean Daily Max and Mean of Daily Min\nD1 Weather Station and K1",c("Mean Daily Maximum", "Mean Daily Minimum"), "Monthly mean of daily temperature (C)", c(-40,35), "K1, 1963-1964", "D1, 1958-1969")

#Mean of Daily Max and Min
C1_recent_mean_maxmin <- krear_wthr_stn_f(c1_wthr_stats, "Recent Years", krear_site_data, "Mean_of_Daily_Max_and_Min", "Mean of Daily Max and Min", "Recent Mean of Daily Max and Min\nC1 Weather Station and K1","Mean of Daily\nMaximum & Minimum", "Monthly mean of daily temperature (C)", c(-40,35), "K1, 2013-2020", "C1, 2013-2020")
D1_recent_mean_maxmin <- krear_wthr_stn_f(d1_wthr_stats, "Recent Years", krear_site_data, "Mean_of_Daily_Max_and_Min", "Mean of Daily Max and Min", "Recent Mean of Daily Max and Min\nD1 Weather Station and K1","Mean of Daily\nMaximum & Minimum", "Monthly mean of daily temperature (C)", c(-40,35), "K1, 2013-2020", "D1, 2013-2020")

C1_historical_mean_maxmin <- krear_wthr_stn_f(c1_wthr_stats, "Historical Years", krear_site_data, "Mean_of_Daily_Max_and_Min", "Mean.Daily.MaxMin.T", "Historical Mean Daily Max and Min\nC1 Weather Station and K1","Mean of Daily\nMaximum & Minimum", "Monthly mean of daily temperature (C)", c(-40,35), "K1, 1963-1964", "C1, 1958-1969")
D1_historical_mean_maxmin <- krear_wthr_stn_f(d1_wthr_stats, "Historical Years", krear_site_data, "Mean_of_Daily_Max_and_Min", "Mean.Daily.MaxMin.T", "Historical Mean Daily Max and Min\nD1 Weather Station and K1","Mean of Daily\nMaximum & Minimum", "Monthly mean of daily temperature (C)", c(-40,35), "K1, 1963-1964", "D1, 1958-1969")


#Days Above Freezing
C1_recent_mean_maxmin <- krear_wthr_stn_f(c1_wthr_stats, "Recent Years", krear_site_data, "Days_Above_Freezing", "Days above freezing", "Recent Number of Days Above Freezing\nC1 Weather Station and K1", "Days Above Freezing", "Number of Days", c(-3,40), "K1, 2013-2020", "C1, 2013-2020")
D1_recent_mean_maxmin <- krear_wthr_stn_f(d1_wthr_stats, "Recent Years", krear_site_data, "Days_Above_Freezing","Days above freezing",  "Recent Number of Days Above Freezing\nD1 Weather Station and K1", "Days Above Freezing", "Number of Days", c(-3,40), "K1, 2013-2020", "D1, 2013-2020")

C1_historical_mean_maxmin <- krear_wthr_stn_f(c1_wthr_stats, "Historical Years", krear_site_data, "Days_Above_Freezing", "Days Above Freezing", "Historical Number of Days Above Freezing\nC1 Weather Station and K1", "Days Above Freezing", "Number of Days", c(-3,40), "K1, 1963-1964", "C1, 1958-1969")
D1_historical_mean_maxmin <- krear_wthr_stn_f(d1_wthr_stats, "Historical Years", krear_site_data, "Days_Above_Freezing", "Days Above Freezing", "Historical Number of Days Above Freezing\nD1 Weather Station and K1", "Days Above Freezing", "Number of Days", c(-3,40), "K1, 1963-1964", "D1, 1958-1969")


dev.off()

