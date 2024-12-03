#import shark_data to bypass data compilation and modification (Section 1)
shark_data <- read.csv("shark_data.csv", header = TRUE, sep =",")

#Packages 
library(dplyr)
library(stringr)
library(rfishbase)
library(ggOceanMaps)
library(lme4) 
library(lmerTest) 
library(DHARMa) 
library(car)
library(ggeffects)
library(ggplot2)
library(sjPlot)
library(ggpubr)
library(patchwork)

#---------Section 1: Data compilation and modification
#---------Section 1.1: load, standardize and compile all csvs from sharkipedia

# Define the path to CSV files
path <- "/Users/joelholmes/Documents/test csv/csv"
#path for files from Marcus
path2 <- "/Users/joelholmes/Documents/Data"

# Get list of CSV files 
files <- list.files(path = path, pattern = "\\.csv$", full.names = TRUE)
files2 <- list.files(path = path2, pattern = "\\.csv$", full.names = TRUE)

#creating an empty list to store zscore data frames
standardized_data_list <- list()

#create standardize function to calc z-scores
standardize <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)}

file_index <- 1
#standardize each file for value column
for (file in files) {data <-
  read.csv(file)
data[, 2] <- as.numeric(data[, 2]) #to fix data type
data[, 2] <- ifelse(is.na(data[, 2]), 0, data[, 2]) #replace NA with 0
data[,2] <- standardize(data[,2])
data$ID <- sprintf("ID_%03d", file_index)
file_index <- file_index + 1
standardized_data_list[[length(standardized_data_list)+1]] <- data}

for (file in files2) {data <-
  read.csv(file)
data[, 2] <- as.numeric(data[, 2]) #to fix data type
data[, 2] <- ifelse(is.na(data[, 2]), 0, data[, 2]) #replace NA with 0
data[,2] <- standardize(data[,2])
data$ID <- sprintf("ID_%03d", file_index)
file_index <- file_index + 1
standardized_data_list[[length(standardized_data_list)+1]] <- data}

#combining all standardized csvs using dplyr package
combined_data <- bind_rows(standardized_data_list)

#change column names
colnames(combined_data)<- c("year","z_value","species","unit","ocean","lat","long","ID")

#---------Section 1.2: modifying data and quality control (some quality control done manually in excel)

#remove rows where multiple species have been combined
shark_data <- combined_data[!is.na(combined_data$species) & !str_detect(combined_data$species, ","), ]

#remove empty species 
shark_data <- shark_data %>% 
  filter(!is.na(species) & species != "")

#code to remove typos in ocean names using stringr package
shark_data <- shark_data %>%
  mutate(ocean = str_replace_all(ocean, c("atlatic" = "atlantic", "pacifc" = "pacific")))

#Add order and superorder by species (manually collected from Catelogue of Life)
species_order_super <- read.csv("species_order_super.csv", header = TRUE, sep =",")
shark_data <- left_join(shark_data, species_order_super, by = "species")

#add timeseries length range (max-min)
year_range <- shark_data %>% group_by(ID) %>% summarise(min_year = min(year), max_year = max(year))
year_range <- year_range %>% group_by(ID) %>% summarise(max_year-min_year)
colnames(year_range) <- c("ID", "tsrange")

shark_data <- left_join(shark_data,year_range, by="ID")

#Trend column -> regression slope
shark_data <- shark_data %>% group_by(ID) %>%
  summarise(slope = coef(lm(z_value ~ year))[2]) %>% left_join(shark_data, by = "ID")

#clean house 
rm(combined_data)
rm(species_order_super)
rm(year_range)

#---------Section 1.3: adding columns of data from rfishbase
# Remove leading/trailing spaces
shark_data$species <- trimws(shark_data$species) 

#Validate species names
species_data <- as.data.frame(shark_data$species)
species_data <- unique(species_data) #remove duplicates

species_data$valid_species <- validate_names(species_data$`shark_data$species`) 
species_input <- as.data.frame(unique(species_data$valid_species))

colnames(species_data) <-c("species","valid_species")
shark_data <- left_join(shark_data, species_data, by = "species")

#Extract data from rfishbase package
#add common names
common_species <- common_names(species_list = species_input$`unique(species_data$valid_species)`)
combined_common_names <- common_species %>% group_by(Species) %>%
  summarise(ComName = paste(unique(ComName), collapse = ", "))
colnames(combined_common_names) <- c("valid_species", "common_names")
shark_data <- left_join(shark_data, combined_common_names, by = "valid_species")

#length
length_data <- species(species_input$`unique(species_data$valid_species)`, fields = c("Species","Length"))
colnames(length_data) <- c("valid_species","MaxTotlength")
shark_data <- left_join(shark_data,length_data,by="valid_species")

#add dist2land using ggoceanmaps package
shark_data <- dist2land(shark_data, lon = ncol(9), lat = ncol(8))

#clean house 
rm(length_data)
rm(species_data)
rm(species_input)
rm(common_species)
rm(combined_common_names)

#---------Section 2: Analysis

#---------Section 2.1: binomial slope glmer testing using lme4, lmerTest, and DHARMa packages

#make slope binomial (negative = 0, positive= 1)
shark_data$binomial_slope <- ifelse(shark_data$slope>0, 1, 0)

#make subset grouping by ID and summarising for testing -> 1 slope per time series
test_subset <- shark_data%>% group_by(ID,species,ocean,tsrange,ldist,order,superorder,MaxTotlength,binomial_slope,slope, lat, long) %>% summarise()

#use glmer with binomial family
#---------Section 2.1.1: length 
glength_model <- glmer(binomial_slope~MaxTotlength+(1|species), data=test_subset, family = binomial("logit") )
summary(glength_model)
glength_modelx <- simulateResiduals(glength_model)
plot(glength_modelx)

#ggeffects - trial and error find x point where crosses 50% (length threshold)
v <- c(206,213,220)
predict_response(glength_model, terms = "MaxTotlength [v]")

#Test for an interaction between length and superorder
glength_super_model <- glmer(binomial_slope~MaxTotlength*superorder+(1|species), data=test_subset, family = binomial("logit"))
summary(glength_super_model)
#Uses car package
Anova(glength_super_model)

#subset
shark_subset <- subset(test_subset, superorder == "Selachimorpha")
ray_subset <-subset(test_subset, superorder == "Batoidea")
chimaera_subset <- subset(test_subset, superorder == "Holocephalimorpha")

#rerun models for all three superorders
glength_shark_model <- glmer(binomial_slope~MaxTotlength+(1|species), data=shark_subset, family = binomial("logit") )
summary(glength_shark_model)
#ggeffects - trial and error find x point where crosses 50% (length threshold)
v <- c(143,151,159)
predict_response(glength_shark_model, terms = "MaxTotlength [v]")

#How many species are longer than threshold
threshold_sharks <- shark_subset %>% filter(MaxTotlength>151)
threshold_sharks <- threshold_sharks %>% group_by(species) %>% summarise()

glength_ray_model <- glmer(binomial_slope~MaxTotlength+(1|species), data=ray_subset, family = binomial("logit") )
summary(glength_ray_model)

glength_chimaera_model <- glmer(binomial_slope~MaxTotlength+(1|species), data=chimaera_subset, family = binomial("logit") )
summary(glength_chimaera_model)

#clean house
rm(glength_modelx)
rm(chimaera_subset)
rm(glength_super_model)
rm(glength_chimaera_model)


#---------Section 2.1.2 dist2land
gdistance_model <- glmer(binomial_slope~ldist+(1|species), data=test_subset, family = binomial("logit") )
summary(gdistance_model)
gdistance_modelx <- simulateResiduals(gdistance_model)
plot(gdistance_modelx)

#ggeffects - trial and error find x point where crosses 50% (distance threshold)
v <- c(134,147, 160)
predict_response(gdistance_model, terms = "ldist [v]")
threshold_distance <- shark_subset %>% filter(ldist<22.2)
#clean house
rm(gdistance_modelx)

#---------Section 2.2: SLiding window analysis - looking for uplift in trend

# Create a data frame to store results
sliding_window_results <- data.frame(
  start_year = integer(),
  end_year = integer(),
  slope = numeric(),
  p_value = numeric(),
  r_squared = numeric(),
  std_error = numeric(),
  stringsAsFactors = FALSE
)

# Define the range of years in the dataset
start_year <- 1950
end_year <- 2020

# Loop through the dataset with a sliding window of 20 years, moving by 1 year
for (start in seq(start_year, end_year - 19)) {
  # Define the end of the 10-year period
  end <- start + 19
  
  # Subset the data to only include the current 20-year window
  shark_window <- shark_data %>% filter(year >= start & year <= end)
  
  # Fit a linear model on the window
  model <- lm(z_value ~ year, data = shark_window)
  
  # Extract relevant statistics (slope, p-value, R-squared, std error)
  slope <- coef(model)["year"]
  p_value <- summary(model)$coefficients[2, 4]
  r_squared <- summary(model)$r.squared
  std_error <- summary(model)$coefficients[2, 2]
  
  # Store the results in the results data frame
  sliding_window_results <- rbind(sliding_window_results, data.frame(
    start_year = start,
    end_year = end,
    slope = slope,
    p_value = p_value,
    r_squared = r_squared,
    std_error = std_error
  ))
}

#clean house
rm(shark_window)

#---------Section 2.3: abundance and abundance trend over time 

zscore_model <- lmer(z_value~year+(1|species), data = shark_data)
summary(zscore_model)

#clean house
rm(zscore_model)

slopeyear_model <- glmer(binomial_slope~year+(1|species), data=shark_data, family = binomial("logit") )
summary(slopeyear_model)

#get predicted probabilities from estimate
plogis(0.044161)
plogis(-87.9376)
#---------Section 2.4: species trends since 1990s which have uplift

#Make loop to look at each species 

uplift_subset <- shark_data %>% filter(year>1990)

uplift_results <- data.frame(
  species = character(),
  slope = numeric(),
  p_value = numeric(),
  r_squared = numeric(),
  std_error = numeric(),
  stringsAsFactors = FALSE
)


# Loop through the dataset 
for (species_name in unique(uplift_subset$species)) {
 
  uplift_species <- shark_data %>% filter(species==species_name)

  if (nrow(uplift_species) < 20) {
    next  
  }
  # Fit a linear model 
  model <- lm(z_value ~ year, data = uplift_species)
  
  # Extract relevant statistics
  slope <- coef(model)["year"]
  p_value <- summary(model)$coefficients[2, 4]
  r_squared <- summary(model)$r.squared
  std_error <- summary(model)$coefficients[2, 2]
  
  # Store the results in the results data frame
  uplift_results <- rbind(uplift_results, data.frame(
    species=species_name,
    slope = slope,
    p_value = p_value,
    r_squared = r_squared,
    std_error = std_error
  ))
}


write.csv(uplift_results, "uplift_results.csv", row.names = FALSE)

#Find number of species sig pos and neg (and not significant)
uplift_results$sig <- ifelse(uplift_results$p_value<0.05, "significant", "not")
uplift_results$trend <- ifelse(uplift_results$slope<0, "negative", "positive")
uplift_count <- uplift_results %>% group_by(sig,trend) %>% summarise(count=n())


#clean house
rm(uplift_subset)
rm(uplift_species)
rm(uplift_count)
rm(uplift_results)

#---------Section 3: Plots and Figures 
#---------Section 3.1: Figure 1 - data description

#Plot density of year freq.
dist_year <-ggplot(shark_data, aes(x=year))+ 
  geom_vline(aes(xintercept = median(year)), color= "red", linetype= "dashed", linewidth=1)+
  geom_density(alpha=.2, fill="darkblue")+
  theme_classic()+
  labs(x="Year", y="Proportion of data")+
  scale_x_continuous(breaks = seq(from=1900, to=2020, by= 20), expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0), limits = c(0,0.04))+
  theme(axis.text.x = element_text(size = 16, family = "Arial"),
        axis.text.y = element_text(size = 16, family = "Arial"),
        axis.title.x = element_text(size = 20, family = "Arial"),
        axis.title.y = element_text(size = 20, family = "Arial"),
        plot.margin = margin(3, 1, 1, 1, "cm"))

ggsave("dist_year.png", plot = dist_year, width = 10, height =6.5, dpi = 300)

#clean house
rm(dist_year)

#dist of ts length
ts_dist <- shark_data %>% group_by(ID,tsrange) %>% summarise()
dist_tsrange <- ggplot(ts_dist, aes(x=tsrange))+ 
  #geom_histogram(aes(y=..density..),binwidth = 10, center = 5, colour= "black",fill = "white") +
  geom_vline(aes(xintercept = median(tsrange)), color= "red", linetype= "dashed", linewidth=1)+
  geom_density(alpha=.2, fill="darkblue")+
  theme_classic()+
  labs(x="Length of time series (years)", y="Proportion of timeseries")+
  scale_x_continuous(breaks = seq(from=0, to=120, by= 20), expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0), , limits = c(0,0.04))+
  theme(axis.text.x = element_text(size = 16, family = "Arial"),
        axis.text.y = element_text(size = 16, family = "Arial"),
        axis.title.x = element_text(size = 20, family = "Arial"),
        axis.title.y = element_text(size = 20, family = "Arial"),
        plot.margin = margin(2, 1, 0, 0, "cm"))

ggsave("dist_tsrange.png", plot = dist_tsrange, width = 10, height = 6, dpi = 300)

#clean house
rm(ts_dist)
rm(dist_tsrange)

#ggoceanmap data distribution plot
lat_lon_counts <- shark_data %>% group_by( lat, long) %>% summarise(count = n_distinct(ID))
data_map <- basemap(limits = c(-180, 180, -70, 80), bathymetry = TRUE, bathy.style = "rcb") + 
  geom_point(data = lat_lon_counts, aes(x = long, y = lat, size = count), 
             color = "darkblue", alpha = 0.6) +
  scale_size_continuous(range = c(2,15)) + 
  labs(x = "", y = "", title = "", size = "Time series") +
  theme_classic()+
  theme(axis.text.x = element_text(size = 16, family = "Arial"),
        axis.text.y = element_text(size = 16, family = "Arial"),
        legend.text = element_text(size = 12, family = "Arial"),
        legend.title = element_text(size = 14, family = "Arial"),
        legend.key.size = unit(1.5, "lines"))+
  guides(fill= "none")

ggsave("dataglobe.png", plot = data_map, width = 10, height = 5, dpi = 300)

#clean house
rm(lat_lon_counts)
rm(data_map)

#Barplot number of ts per Chondricthyan order - add silhouettes from phylopic website
order_ts_data <- shark_data %>% 
  group_by(species, ID, order, superorder) %>% 
  summarise(n = n())

order_ts_count <- order_ts_data %>% 
  group_by(order, superorder) %>% 
  summarise(n = n())

orderbarplot <- ggplot(data = order_ts_count, aes(x = reorder(order, n), y=n, fill = superorder))+ 
  geom_bar(stat = "identity", width=0.6) +
  coord_flip() + theme_classic() +
  scale_fill_manual(values = c("Selachimorpha" = "darkblue", "Batoidea" = "turquoise", "Holocephalimorpha" = "violet"))+
  xlab("Taxonomic Order")+ylab("Timeseries count")+
  scale_y_continuous(limits=c(0,320),breaks = seq(from=0, to=300, by= 50), expand = c(0, 0))+
  labs(fill="Superorder")+
  theme(axis.text.x = element_text(size = 16, family = "Arial"),
        axis.text.y = element_text(size = 16, family = "Arial"),
        plot.margin = margin(1, 1, 1, 1, "cm"),
        axis.title.x = element_text(size = 20, family = "Arial"),
        axis.title.y = element_text(size = 20, family = "Arial"),
        legend.text = element_text(size = 12, family = "Arial"),
        legend.title = element_text(size = 14, family = "Arial")
        ,legend.position = "none")

ggsave("orderbar21.png", plot = orderbarplot, width = 12, height = 7, dpi = 300)

#clean house
rm(order_ts_data)
rm(order_ts_count)
rm(orderbarplot)

#---------Section 3.2: Figure 2 - bislope~Length for all, sharks, and rays (sj plots)

#sjplot glength_model with data from all 3 superorders 
glengthplot<- plot_model(glength_model, type= "pred", terms = c("MaxTotlength"), color= "darkblue")+  
  labs(
  title = "",
  x = "",
  y = "")+ 
  theme_classic()+
  geom_hline(yintercept = 0.5, color= "red", size=0.5)+
  ylim(0,1)+
  scale_x_continuous(breaks = seq(from=0, to=1800, by= 200),limits = c(0, 1800), expand = c(0, 0))+
  scale_y_continuous(breaks = seq(from=0, to=1, by= 0.25),limits = c(0, 1), expand = c(0, 0))+
  theme(axis.text.x = element_text(size = 16, family = "Arial"),
        axis.text.y = element_text(size = 16, family = "Arial"),
        plot.margin = margin(0, 0.7, 1, 1, "cm"))


all_length_dist <- ggplot(test_subset, aes(x=MaxTotlength))+ 
  geom_histogram(aes(y=..density..),binwidth = 50, center = 25, colour= "black",fill = "black") +
  geom_vline(aes(xintercept = median(MaxTotlength, na.rm= TRUE)), color= "red", linetype= "dashed", size=1)+
  theme_classic()+
  labs(x="", y="")+
  scale_x_continuous(breaks = seq(from=0, to=1800, by= 200),limits = c(0, 1800), expand = c(0, 0))+
  scale_y_continuous(breaks = seq(from=0, to=0.01, by= 0.0025),limits = c(0, 0.01), expand = c(0, 0))+
  theme(axis.title.x = element_blank(),  
        axis.text.x = element_blank(),   
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 16, family = "Arial"),
        plot.margin = margin(2, 0.7, 1, 0, "cm"))

combined_plot <- all_length_dist / glengthplot + plot_layout(heights = c(1, 4)) & theme(plot.margin = margin(1, 1, 0, 0, "cm"))

#shark plot
glength_shark_plot<- plot_model(glength_shark_model, type= "pred", terms = c("MaxTotlength"), color = "darkblue") +
  labs(
    title = "",
    x = "",
    y = "")+ 
  geom_hline(yintercept = 0.5, color= "red", size=0.5)+
  theme_classic()+
  ylim(0,1)+
  scale_x_continuous(breaks = seq(from=0, to=1800, by= 200),limits = c(0, 1800), expand = c(0, 0))+
  scale_y_continuous(breaks = seq(from=0, to=1, by= 0.25),limits = c(0, 1), expand = c(0, 0))+
  theme(axis.text.x = element_text(size = 16, family = "Arial"),
        axis.text.y = element_blank(),#element_text(size = 16, family = "Arial"),
        plot.margin = margin(0, 0.7, 0, 1, "cm"))

shark_length_dist <- ggplot(shark_subset, aes(x=MaxTotlength))+ 
  geom_histogram(aes(y=..density..),binwidth = 50, center = 25, colour= "black",fill = "black") +
  geom_vline(aes(xintercept = median(MaxTotlength, na.rm= TRUE)), color= "red", linetype= "dashed", size=1)+
  theme_classic()+
  labs(x="", y="")+
  scale_x_continuous(breaks = seq(from=0, to=1800, by= 200),limits = c(0, 1800), expand = c(0, 0))+
  scale_y_continuous(breaks = seq(from=0, to=0.01, by= 0.0025),limits = c(0, 0.01), expand = c(0, 0))+
  theme(axis.title.x = element_blank(),  # Remove x-axis title
        axis.text.x = element_blank(),   # Remove x-axis labels
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),#element_text(size = 16, family = "Arial"),
        plot.margin = margin(2, 0.7, 0, 0, "cm"))

combined_shark_plot <- shark_length_dist / glength_shark_plot + plot_layout(heights = c(1, 4)) & theme(plot.margin = margin(1, 1, 0, 0, "cm"))


#Ray plot
glength_ray_plot<- plot_model(glength_ray_model, type= "pred", terms = c("MaxTotlength"), color = "darkblue")+
  labs(
  title = "",
  x = "",
  y = "")+ 
  geom_hline(yintercept = 0.5, color= "red", size=0.5)+
  theme_classic()+
  ylim(0,1)+
  scale_x_continuous(breaks = seq(from=0, to=1800, by= 200),limits = c(0, 1800), expand = c(0, 0))+
  scale_y_continuous(breaks = seq(from=0, to=1, by= 0.25),limits = c(0, 1), expand = c(0, 0))+
  theme(axis.text.x = element_text(size = 16, family = "Arial"),
        axis.text.y = element_blank(),#element_text(size = 16, family = "Arial"),
        plot.margin = margin(0, 1, 0, 1, "cm"))


ray_length_dist <- ggplot(ray_subset, aes(x=MaxTotlength))+ 
  geom_histogram(aes(y=..density..),binwidth = 50, center = 25, colour= "black",fill = "black") +
  geom_vline(aes(xintercept = median(MaxTotlength, na.rm= TRUE)), color= "red", linetype= "dashed", size=1)+
  theme_classic()+
  labs(x="", y="")+
  scale_x_continuous(breaks = seq(from=0, to=1800, by= 200),limits = c(0, 1800), expand = c(0, 0))+
  scale_y_continuous(breaks = seq(from=0, to=0.01, by= 0.0025),limits = c(0, 0.01), expand = c(0, 0))+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),#element_text(size = 16, family = "Arial"),
        plot.margin = margin(2, 1, 0, 0, "cm"))

combined_ray_plot <- ray_length_dist / glength_ray_plot + plot_layout(heights = c(1, 4)) & theme(plot.margin = margin(1, 1, 0, 0, "cm"))

#Use patchwork package to arrange
figure2 <- combined_plot | combined_shark_plot | combined_ray_plot 

ggsave("figure2.png", plot = figure2, width= 20, height = 10, dpi = 300)


#clean house
rm(glength_model,glengthplot, all_length_dist, combined_plot)
rm(shark_subset,glength_shark_model, glength_shark_plot, shark_length_dist, combined_shark_plot)
rm(ray_subset, glength_ray_model, glength_ray_plot, ray_length_dist, combined_ray_plot)
rm(figure2)


#---------Section 3.3: Figure 3 - spatial and temporal chondrichthyan trends

#Map showing global trends over full temporal scale of the dataset
slope_summary <- test_subset %>% group_by(lat,long) %>% summarise(mean_slope = mean(slope))

slope_map <- basemap( limits = c(-180, 180, -70, 80), bathymetry = TRUE, bathy.style = "rcb") + 
  geom_point(data = slope_summary, aes(x = long,
                                     y = lat, size = abs(mean_slope), 
                                     color = ifelse(mean_slope < 0, "Negative", "Positive")), alpha = 0.5) +
  scale_size_continuous(range = c(2, 15), limits = c(0, 0.45)) +
  scale_color_manual(
    values = c("Negative" = "red3", "Positive" = "blue3"),
    name = "Trend direction")+
  labs(title = "", x = "", y = "") +
  theme_classic()+
  theme(axis.text.x = element_text(size = 16, family = "Arial"),
        axis.text.y = element_text(size = 16, family = "Arial"),
        legend.position = "none")

#Early slope map - 1920-1990
early_map <- shark_data %>% filter(between(year, 1920,1990))
early_map <- early_map %>% group_by(ID) %>%
  summarise(early_slope = coef(lm(z_value ~ year))[2]) %>% left_join(early_map, by = "ID")
early_map <- early_map %>% filter(!is.na(early_slope) & early_slope != 0)
early_summary <- early_map%>% group_by(ID, lat, long, early_slope) %>% summarise()
early_summary <- early_summary %>% group_by(lat,long) %>% summarise(meanearly_slope = mean(early_slope))

early_slope_map <- basemap(
  limits = c(-180, 180, -70, 80), , bathymetry = TRUE, bathy.style = "rcb") + 
  geom_point(data = early_summary, aes(
    x = long,
    y = lat,
    size = abs(meanearly_slope),  
    color = ifelse(meanearly_slope < 0, "Negative", "Positive")), 
    alpha = 0.5) +
  scale_size_continuous(range = c(2, 15), limits = c(0, 0.45), name = "Slope Value") +
  scale_color_manual(
    values = c("Negative" = "red3", "Positive" = "blue3"),
    name = "Trend direction")+
  labs(title = "", x = "", y = "") +
  theme_classic()+
  guides(fill= "none", size = guide_legend(override.aes = list(color = "red3", alpha = 1)), color = "none")+
  theme(axis.text.x = element_text(size = 16, family = "Arial"),
        axis.text.y = element_text(size = 16, family = "Arial"),
        legend.text = element_text(size = 12, family = "Arial"),
        legend.title = element_text(size = 14, family = "Arial"),
        legend.key.size = unit(1.5, "lines"))

#Late map 1990-2020

late_map <- shark_data %>% filter(between(year, 1990,2020))
late_map <- late_map %>% group_by(ID) %>%
  summarise(late_slope = coef(lm(z_value ~ year))[2]) %>% left_join(late_map, by = "ID")
late_map <- late_map %>% filter(!is.na(late_slope) & late_slope != 0)
late_summary <- late_map%>% group_by(ID, lat, long, late_slope) %>% summarise()
late_summary <- late_summary %>% group_by(lat,long) %>% summarise(meanlate_slope = mean(late_slope))

late_slope_map <- basemap(
  limits = c(-180, 180, -70, 80), , bathymetry = TRUE, bathy.style = "rcb") + 
  geom_point(data = late_summary, aes(
    x = long,
    y = lat,
    size = abs(meanlate_slope),  
    color = ifelse(meanlate_slope < 0, "Negative", "Positive")),
    alpha = 0.5) +
  scale_size_continuous(range = c(2, 15), limits = c(0, 0.45)) +
  scale_color_manual(
    values = c("Negative" = "red3", "Positive" = "blue3"),
    name = "Trend direction")+
  labs(title = "", x = "", y = "") +
  theme_classic()+
  guides(fill= "none", size = guide_legend(override.aes = list(color = "blue3", alpha = 1)), color = "none")+
  theme(axis.text.x = element_text(size = 16, family = "Arial"),  
        axis.text.y = element_text(size = 16, family = "Arial"),
        legend.text = element_text(size = 12, family = "Arial"),
        legend.title = element_text(size = 14, family = "Arial"),
        legend.key.size = unit(1.5, "lines"))

#Use patchwork package
slope_maps <- (slope_map/early_slope_map/late_slope_map) 

ggsave("slope_maps.png", plot = slope_maps, width = 10, height = 12, dpi = 300)

#clean house
rm(early_map)
rm(early_summary)
rm(early_slope_map)
rm(late_map)
rm(late_summary)
rm(late_slope_map)
rm(early_late)


#plot sliding window analysis results
sliding_window_results$start_year <- as.numeric(sliding_window_results$start_year)

sliding_window <- ggplot(sliding_window_results, aes(x = start_year, y = slope)) +
  geom_ribbon(aes(ymin=(slope-std_error), ymax= (slope+std_error)), fill= "darkblue", alpha = 0.2)+
  geom_line(color="darkblue")+
  geom_hline(yintercept = 0, color= "red", linewidth=0.5)+
  theme_classic() +
  scale_x_continuous(breaks = seq(from=1950, to=2000, by= 10),limits = c(1950, 2001))
  labs(title = "",
       x = "Year",
       y = "    Decreasing <- Trend in Abundance -> Increasing")+
  theme(axis.text.x = element_text(size = 16, family = "Arial"),
        axis.text.y = element_text(size = 16, family = "Arial"),
        axis.title.x = element_text(size = 20, family = "Arial"),
        axis.title.y = element_text(size = 20, family = "Arial"))

ggsave("sliding_window.png", plot = sliding_window, width = 10, height = 9, dpi = 300)


#clean house 
rm(sliding_window_results)
rm(sliding_window)

#plot slope~distance

gdistplot<- plot_model(gdistance_model, type= "pred", terms = c("ldist"), color= "darkblue")+  
  labs(
    title = "",
    x = "",
    y = "")+ 
  theme_classic()+
  geom_hline(yintercept = 0.5, color= "red", linewidth=0.5)+
  ylim(0,1)+
  scale_x_continuous(breaks = seq(from=0, to=2000, by= 250),limits = c(0, 2000), expand = c(0, 0))+
  scale_y_continuous(breaks = seq(from=0, to=1, by= 0.25),limits = c(0, 1), expand = c(0, 0))+
  theme(axis.text.x = element_text(size = 16, family = "Arial"),
        axis.text.y = element_text(size = 16, family = "Arial"),
        plot.margin = margin(1, 1, 1, 1, "cm"))

ggsave("distanceplot.png", plot = gdistplot, width = 10, height = 7.5, dpi = 300)