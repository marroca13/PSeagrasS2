library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(distributional)
library(forcats)
library(hrbrthemes)
library(viridis)

#-------
# Directory
#-------

wd <- getwd()
filename <- paste0(wd, '/Signature_sampledData_2024_OSW.csv')

#-------
# Import Data
#-------

# Read a comma-separated CSV file
data <- read.csv(filename, header = TRUE, stringsAsFactors = FALSE)
head(data)  # Preview the first few rows

names(data) <- sub("X", "", names(data))

# Create a new column 'class' based on the values in 'value'
data <- data %>%
  mutate(
    class = case_when(
      FONDO == "Halodule" ~ "Halodule",
      FONDO == "Ruppia" ~ "Ruppia"
    )
  )


####################### PlanetScope 2018
# Select relevant bands and class column
data <- data %>%
  select(PS_B1, PS_B2, PS_B3, PS_B4, class)

# Rename columns to their respective wavelengths
data <- data %>%
  rename(
    "485" = PS_B1,
    "565" = PS_B2,
    "665" = PS_B3,
    "865" = PS_B4
  )


# Create a new column 'class' based on the values in 'value' 2024
data <- data %>%
  mutate(
    sensor = "PS")

# PS data <- select(data, c("444", "490", "531", "565", "610", "666", "705", "865", "class"))
# data <- select(data, c("443", "492", "560", "666", "704", "740", "783", "834", "866", "class"))

# set.seed(123)
# data <- data.frame(
#   Band = rep(c("B1", "B2", "B3", "B4", "B5"), times = 3 * 500), # 5 bandas, 3 clases, 50 observaciones
#   Reflectance = c(rnorm(500, 0.01, 0.005), rnorm(500, 0.015, 0.005), rnorm(500, 0.02, 0.007)),
#   Class = rep(c("Seagrass", "ODW", "Sand/rock"), each = 500 * 5)
# )

data_long <- data %>%
  pivot_longer(cols = c("485", "565", "665", "865"), names_to = 'Band', values_to = 'Reflectance')
# I prefer controlling the list order of the class here. Sometimes mistakes happen in scale_fill_manual 
# and the wrong plots might be done, potentially causing misclassification
data_long$class <- factor(data_long$class, levels = c('Halodule', 'Ruppia'))


############## PlanetScope SuperDove
# Create a new column 'class' based on the values in 'value' 2024
data <- data %>%
  mutate(
    class = case_when(
      Bottom == "Hadolule" ~ "Halodule",
      Bottom == "Ruppia" ~ "Ruppia"
    )
  ) %>%
  filter(!is.na(class))  # Remove rows where 'class' is NA

# Select relevant bands and class column
data <- data %>%
  select(PS_B1, PS_B2, PS_B3, PS_B4, PS_B5, PS_B6, PS_B7, PS_B8, class)

# Rename columns to their respective wavelengths
data <- data %>%
  rename(
    "444" = PS_B1,
    "490" = PS_B2,
    "531" = PS_B3,
    "565" = PS_B4,
    "610" = PS_B5,
    "665" = PS_B6,
    "705" = PS_B7,
    "865" = PS_B8
  )

# Create a new column 'class' based on the values in 'value' 2024
data <- data %>%
  mutate(
    sensor = "PS")


data_long <- data %>%
  pivot_longer(cols = c("444", "490", "531", "565", "610", "665", "705", "865"), names_to = 'Band', values_to = 'Reflectance')
data_long$class <- factor(data_long$class, levels = c('Halodule', 'Ruppia'))


########################### Sentinel-2
# Select relevant bands and class column
data <- data %>%
  select(S2_B1, S2_B2, S2_B3, S2_B4, S2_B5, S2_B6, S2_B7, S2_B8, class)

# Rename columns to their respective wavelengths
data <- data %>%
  rename(
    "443" = S2_B1,
    "492" = S2_B2,
    "560" = S2_B3,
    "665" = S2_B4,
    "704" = S2_B5,
    "740" = S2_B6,
    "783" = S2_B7,
    "842" = S2_B8
  )

# Create a new column 'class' based on the values in 'value' 2024
data <- data %>%
  mutate(
    sensor = "S2")


data_long <- data %>%
  pivot_longer(cols = c("443", "492", "560", "665", "704", "740", "783", "842"), names_to = 'Band', values_to = 'Reflectance')
# I prefer controlling the list order of the class here. Sometimes mistakes happen in scale_fill_manual 
# and the wrong plots might be done, potentially causing misclassification
data_long$class <- factor(data_long$class, levels = c('Halodule', 'Ruppia'))





data_ribbon_PS <- data_long %>%
  group_by(class, Band) %>%
  summarise(mean_value = mean(Reflectance), sd_value = sd(Reflectance), count = length(Reflectance), quartile_first = quantile(Reflectance, c(0.25)), quartile_third = quantile(Reflectance, c(0.75)), class = class, Band = Band, sensor = sensor) %>%
  distinct() %>%
  ungroup()


data_ribbon_S2 <- data_long %>%
  group_by(class, Band) %>%
  summarise(mean_value = mean(Reflectance), sd_value = sd(Reflectance), count = length(Reflectance), quartile_first = quantile(Reflectance, c(0.25)), quartile_third = quantile(Reflectance, c(0.75)), class = class, Band = Band, sensor = sensor) %>%
  distinct() %>%
  ungroup()


## Merged data
data_ribbon <- full_join(data_ribbon_PS, data_ribbon_S2)

#-------
# Violin Plot
#-------

# Grouped
p_violin <- data_long %>%
  ggplot(aes(fill=class, y=Reflectance, x=Band)) + 
  # geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent") +
  geom_violin(position="dodge", alpha=0.5) +
  scale_fill_manual(values = c("#5C8944", "#FFD37F")) +
  theme_ipsum()  +
  xlab("Bands") +
  ylab("Reflectances") +
  ylim(0, 0.050)
  
print(p_violin)

#-------
# Ribbon plot
#-------
p_ribbon_se <- ggplot(data_ribbon, aes(x = Band, y = mean_value, fill = class, group = class)) +
  geom_ribbon(aes(ymin = mean_value-1.96*(sd_value/sqrt(count)), ymax = mean_value+1.96*(sd_value/sqrt(count))), alpha = 0.1) +
  geom_line(aes(y = mean_value, colour = class), linewidth = 1) +
  scale_fill_manual(values = c("#5C8944", "#FFD37F")) +
  scale_colour_manual(values = c("#5C8944", "#FFD37F")) +
  theme_ipsum()  +
  xlab("Bands") +
  ylab("Reflectances") +
  ylim(0, 0.090)

print(p_ribbon_se)


p_ribbon_percent <- ggplot(data_ribbon, aes(x = Band, y = mean_value, fill = class, group = class)) +
  geom_ribbon(aes(ymin = quartile_first, ymax = quartile_third), alpha = 0.1) +
  geom_line(aes(y = mean_value, colour = class), linewidth = 1) +
  geom_point(aes(colour = class), size = 3) +  # <-- This adds the dots
  scale_fill_manual(values = c("#94b739", "#cdcf6f")) +
  scale_colour_manual(values = c("#94b739", "#cdcf6f")) +
  theme_ipsum()  +
  xlab("Bands") +
  ylab("Reflectances") +
  ylim(0, 0.12) +
  theme(
    axis.title.x = element_text(hjust = 0.5, size = 16, family = "roboto"),
    axis.title.y = element_text(hjust = 0.5, size = 16, family = "roboto"),
    panel.grid.major = element_line(color = "grey80", size = 0.2),
    panel.grid.minor = element_blank()
  )



print(p_ribbon_percent)



------------------

# Create a combined label
data_ribbon <- data_ribbon %>%
  mutate(sensor = paste(class, sensor, sep = "-"))

# Convert Band column to numeric
data_ribbon$Band <- as.numeric(as.character(data_ribbon$Band))

p_ribbon_percent <- ggplot(data_ribbon, aes(x = Band, y = mean_value, fill = class, group = class)) +
  geom_ribbon(aes(ymin = quartile_first, ymax = quartile_third), alpha = 0.1) +
  geom_line(aes(y = mean_value, colour = class), linewidth = 1) +
  scale_fill_manual(values = c("#5C8944", "#FFD37F")) +
  scale_colour_manual(values = c("#5C8944", "#FFD37F")) +
  scale_x_continuous(name = "Wavelength (nm)", limits = c(400, 900), breaks = seq(400, 900, 100)) +
  scale_y_continuous(name = "Reflectance", limits = c(0, 0.09)) +
  theme_ipsum()

print(p_ribbon_percent)


# Plot with new grouping
p <- data_ribbon <- ggplot(data_ribbon, aes(x = Band, y = mean_value, fill = sensor, group = sensor)) +
  geom_ribbon(aes(ymin = quartile_first, ymax = quartile_third), alpha = 0.1) +
  geom_errorbar(aes(ymin = mean_value - sd_value/sqrt(count),
                    ymax = mean_value + sd_value/sqrt(count)),
                width = 0.2, linewidth = 0.7, alpha = 0.7, colour = "grey40")+
  geom_line(aes(y = mean_value, colour = sensor), linewidth = 1) +
  geom_point(aes(colour = sensor), size = 2) +
  scale_fill_manual(values = c(
    "Halodule-PS" = "#cdcf6f",
    "Halodule-S2" = "#cdcf6f",
    "Ruppia-PS" = "#94b739",
    "Ruppia-S2" = "#94b739"
  )) +
  scale_colour_manual(values = c(
    "Halodule-PS" = "#cdcf6f",
    "Halodule-S2" = "#cdcf6f",
    "Ruppia-PS" = "#94b739",
    "Ruppia-S2" = "#94b739"
  )) +
  theme_ipsum() +
  xlab("Wavelength (nm)") +
  ylab("Rrs (sr-1)") +
  ylim(0, 0.11) +
  theme(
    axis.title.x = element_text(hjust = 0.5, size = 14, family = "roboto"),
    axis.title.y = element_text(hjust = 0.5, size = 14, family = "roboto"),
    panel.grid.major = element_line(color = "grey80", size = 0.2),
    panel.grid.minor = element_blank()
  )

print(p)


## Export graph
ggsave(
  filename = "S2_PS_RH_18_errorBar_bar.tiff",
  plot = p,
  width = 800/100, height = 450/100,   # Convert pixels to inches (1 inch = 100px for 300dpi)
  dpi = 300
)

  
----------




################ New code

# Convert Band column to numeric
data_ribbon$Band <- as.numeric(as.character(data_ribbon$Band))

p_ribbon_percent <- ggplot(data_ribbon, aes(x = Band, y = mean_value, fill = class, group = class)) +
  geom_ribbon(aes(ymin = quartile_first, ymax = quartile_third), alpha = 0.1) +
  geom_line(aes(y = mean_value, colour = class), linewidth = 1) +
  scale_fill_manual(values = c("#5C8944", "#FFD37F")) +
  scale_colour_manual(values = c("#5C8944", "#FFD37F")) +
  scale_x_continuous(name = "Wavelength (nm)", limits = c(400, 900), breaks = seq(400, 900, 100)) +
  scale_y_continuous(name = "Reflectance", limits = c(0, 0.12)) +
  theme_ipsum()

print(p_ribbon_percent)

data_ribbon %>%
  filter(class == "nonSeagrass") %>%
  filter(is.na(quartile_first) | is.na(quartile_third))


library(zoo)

data_ribbon <- data_ribbon %>%
  group_by(class) %>%
  mutate(
    quartile_first = na.approx(quartile_first, na.rm = FALSE),
    quartile_third = na.approx(quartile_third, na.rm = FALSE)
  )

data_ribbon <- data_ribbon %>%
  filter(!is.na(quartile_first) & !is.na(quartile_third))

# Same plot code after data cleaning
print(p_ribbon_percent)

