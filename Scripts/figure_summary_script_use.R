## S2MET Registration
## Figure and Table Scripts
## 
## Author: Jeff Neyhart
## Last updated: March 21, 2018
## 

# Run the source script
repo_dir <- getwd()
source(file.path(repo_dir, "source_use.R"))


### Tables ###

## Filter out locations that were not observed
trial_info1 <- trial_info %>%
  filter(!grepl("SMT", environment), environment != "BCW17")


# Summarize the number of years, locations, environments that were planned
trial_info %>% 
  summarize_at(vars(trial:year), n_distinct)
# Summarize the number of years, locations, environments that were actually observed
trial_info1 %>% 
  summarize_at(vars(trial:year), n_distinct)


# Summarize the spanned area
trial_info %>% 
  select(latitude, longitude) %>% 
  gather(coordinate, value) %>% 
  group_by(coordinate) %>% 
  summarize_at(vars(value), funs(min, max), na.rm = T) %>% 
  mutate(range = max - min)


# Summarize the traits in each trial
n_distinct(S2_MET_tidy_use$trait)


## Create a table of locations, number of trials, years, traits, lat/lon
# First summarize n traits per location
location_traits <- S2_MET_tidy_use %>% 
  distinct(trial, environment, location, year, trait)

# Merge and create a table
trial_info_toprint <- trial_info1 %>% 
  # Round the lat/long values
  mutate_at(vars(latitude, longitude), ~round(., 2)) %>%
  left_join(., location_traits, by = c("trial", "environment", "year")) %>% 
  select(location = location.x, state, trial, year, latitude, longitude, population, trait) %>% 
  unite(location, location, state, sep = ", ") %>% 
  group_by(location) %>% 
  summarize(latitude = latitude[1], longitude = longitude[1], 
            trials = n_distinct(trial), years = n_distinct(year), 
            traits = n_distinct(trait)) %>% 
  rename_at(vars(trials:traits), funs(str_c("no. ", .)))

save_file <- file.path(fig_dir, "trial_description_table.csv")
write_csv(x = trial_info_toprint, path = save_file)





## Create a table of traits, units, and number of observations

# Change the names of the traits
trait_info <- S2_MET_tidy_use %>% 
  distinct(trait) %>% 
  left_join(., t3_traits, by = c("trait" = "Nickname")) %>% 
  select(trait, trait_name = Trait, unit = Unit) %>% 
  mutate(trait_name = str_to_title(trait_name))

# Add information of number of observations for each trait
trait_info1 <- S2_MET_tidy_use %>% 
  group_by(trait) %>% 
  summarize(observations = n(), environments = n_distinct(environment)) %>% 
  arrange(desc(observations)) %>%
  left_join(., trait_info)

# Write a table
trait_info_toprint <- trait_info1 %>% 
  select(trait = trait_name, unit, environments, observations) %>%
  rename_all(str_to_title)

save_file <- file.path(fig_dir, "trait_information_table.csv")
write_csv(x = trait_info_toprint, path = save_file)





### Figures ###



## Create a figure to show the difference in field designs
n_row <- 9
n_col <- 15
n_plots <- n_row * n_col

n_blks <- 9
n_check_plots <- 15
n_entries <- n_plots - n_check_plots
n_checks <- 3

# Create randomizations
library(FldTrial)

# Set the seed
set.seed(3900) # This one looks nice

# aibd example
aibd_example <- design.aibd(enviro = "example", exp.name = "example", nEntries = n_entries,
                            nChk2 = n_checks - 1, nFieldRows = n_row, nFieldCols = n_col, 
                            nBlks.min = n_blks)

# rcbd example
rcbd_example <- design.rcbd(enviro = "ex", exp.name = "ex", nBlks = 1, nEntries = n_entries,
                            nChks = n_checks, nChkReps = (n_check_plots / n_checks), 
                            nFieldRows = n_row, nFieldCols = n_col)

# Replacement matrix for line codes
line_code_replace <- c("0" = "Entry", "1" = "Check 1\n(Primary)", "2" = "Check 2", "3" = "Check 3")

# Create a color code for the line_code
line_color <- c("white", umn_palette(name = "Secondary_Tier1", n = 3)) %>% 
  set_names(line_code_replace)


# Function for altering scales
adj_scale <- function(x) set_names(x = seq(1.5, (max(x) - 1) + 0.5, 1), seq(1, max(x) - 1))


# Create a modification for each plot
g_mod <- list(
  geom_rect(color = "black", lwd = 0.5),
  geom_rect(aes(xmin = column_min, xmax = column_max + 1, ymin = row_min, ymax = row_max + 1), 
            fill = "white", alpha = 0, color = "black", lwd = 1),
  # scale_fill_manual(values = line_color),
  scale_fill_grey(start = 1, end = 0),
  scale_x_continuous(breaks = adj_scale, position = "top"),
  scale_y_continuous(breaks = adj_scale),
  xlab("Column"),
  ylab("Row"),
  theme_acs(),
  theme(legend.title = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.position = "bottom")
)


# Plot these
g_aibd <- aibd_example$aibd.dgn %>% 
  mutate(line_code = str_replace_all(line_code, line_code_replace),
         line_code = parse_factor(line_code, levels = line_code_replace)) %>%
  # Create the block dimensions
  group_by(blk) %>%
  mutate_at(vars(row, column), funs(min, max)) %>%
  ggplot(aes(xmin = column, xmax = column + 1, ymin = row, ymax = row + 1, fill = line_code)) + 
  g_mod


g_rcbd <- rcbd_example$rcbd.dgn %>% 
  mutate(line_code = str_replace_all(line_code, line_code_replace),
         line_code = parse_factor(line_code, levels = line_code_replace)) %>%
  mutate_at(vars(row, column), funs(min, max)) %>%
  ggplot(aes(xmin = column, xmax = column + 1, ymin = row, ymax = row + 1, fill = line_code)) + 
  g_mod

# Combine
g_plots <- plot_grid(g_aibd + theme(legend.position = "none"), 
                     g_rcbd + theme(legend.position = "none", axis.text.x = element_blank(),
                                    axis.ticks.x = element_blank(), axis.title.x = element_blank()), 
                     ncol = 1, labels = c("A", "B"), align = "hv")

g_plots1 <- plot_grid(g_plots, get_legend(g_aibd), ncol = 1, rel_heights = c(0.85, 0.15))


# Save the files
ggsave(filename = "field_design_examples.jpg", plot = g_plots1, path = fig_dir,
       width = 3.5, height = 6, dpi = 1000)








### Plot the trial locations and summarize trial sites by the number of years
### a location was included

# First summarize
trial_info_summ <- trial_info %>% 
  group_by(location) %>% 
  summarize(long = longitude[1], lat = latitude[1], 
            n_year = n_distinct(year)) %>%
  filter(!is.na(lat), ) %>%
  mutate(n_year = ifelse(n_year == 3, "3+", n_year))

# Get the map data for canada
canada <- map_data("world", "Canada")

# Download map data for US by county
usa_county <- map_data(map = "county")
# Download state data
usa_state <- map_data(map = "state")

# Adjust the groups in the states
usa_state <- usa_state %>%
  mutate(group = group + max(canada$group))

# Adjust the groups in the counties
usa_county <- usa_county %>%
  mutate(group = group + max(usa_state$group))

# Tidy and combine
north_america <- bind_rows(usa_state, usa_county, canada)

# Map
g_map <- ggplot(data = north_america, aes(x = long, y = lat, group = group)) +
  geom_polygon(fill = "grey85") +
  geom_polygon(data = canada, fill = NA, color = "grey50", lwd = 0.3) + # Add canada
  geom_polygon(data = usa_state, aes(x = long, y = lat, group = group), fill = NA, color = "grey50", lwd = 0.3) +
  geom_point(data = trial_info_summ, aes(x = long, y = lat, group = location, size = n_year)) +
  coord_fixed(ratio = 1.5, xlim = c(-125, -60), ylim = c(35, 50)) +
  scale_size_manual(name = "Years in\nExperiment", values = c(1, 2.5, 4)) +
  xlab("Longitude") +
  ylab("Latitude") + 
  theme_acs() + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        # legend.title = element_text(size = 10),
        # legend.text = element_text(size = 8),
        legend.position = "bottom")

# Add north mark and scale bar
g_map1 <- g_map + 
  ggsn::north(location = "bottomright", symbol = 12, x.min = -125, x.max = -60, 
              y.min = 37, y.max = 50) + 
  ggsn::scalebar(x.min = -125, x.max = -60, y.min = 35, y.max = 50, dist = 500, 
                 dd2km = TRUE, model = "WGS84", location = "bottomright", st.size = 1.5, 
                 st.dist = 0.05, st.bottom = FALSE)

# Save the figure
save_file <- file.path(fig_dir, "trial_location_map.jpg")
ggsave(filename = save_file, plot = g_map1, width = 3.5, height = 2, dpi = 1000)











