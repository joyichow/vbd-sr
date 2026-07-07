## Panel A: geographical distributions of studies
## Panel B: stacked bar chart, timeline of study's publication year
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(rnaturalearth)
library(sf)
library(cowplot)

raw <- read_excel("vbd_S1 Table.xlsx", sheet = 1, skip = 1, col_names = TRUE)

names(raw)[1] <- "study_label"
country_col <- grep("Country|WHO", names(raw), value = TRUE)[1]
if (is.na(country_col)) country_col <- names(raw)[4]

# the table has intervention class as header rows (no country value)
current_class <- NA
studies <- list()

for (i in seq_len(nrow(raw))) {
  first   <- raw[[i, "study_label"]]
  country <- raw[[i, country_col]]
  if (is.na(first) || first == "") next
  if (is.na(country) || country == "") {
    current_class <- first
    next
  }
  studies[[length(studies)+1]] <- data.frame(
    study = first, country_raw = country,
    intervention_class = current_class,
    stringsAsFactors = FALSE)
}
df <- bind_rows(studies)

# parse year
df$year <- as.integer(str_extract(df$study, "\\d{4}$"))

# parse country
df$country_clean <- str_remove_all(df$country_raw, "\\s*\\([^)]*\\)")
df <- df %>%
  separate_rows(country_clean, sep = "[,;]\\s*") %>%
  mutate(country_clean = str_trim(country_clean)) %>%
  filter(country_clean != "")

suppression_studies <- c("de Castro Poncio, 2023", "Lim, 2023", "Lim, 2024", "Lim, 2025", "Lim, 2026")

df <- df %>%
  mutate(intervention_class = case_when(
    intervention_class == "Mosquito modification" & study %in% suppression_studies ~ "Mosquito modification (suppression)",
    intervention_class == "Mosquito modification"                                  ~ "Mosquito modification (replacement)",
    TRUE ~ intervention_class
  ))

# Panel A: Geographical distribution of studies
country_counts <- df %>%
  distinct(study, country_clean) %>%
  count(country_clean, name = "n_studies")

world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  select(name, geometry)

# match
world <- world %>%
  left_join(country_counts, by = c("name" = "country_clean"))

pA <- ggplot(world) +
  geom_sf(aes(fill = n_studies), colour = "grey70", linewidth = 0.15) +
  scale_fill_gradient(low = "#c6dbef", high = "#08306b",
                      na.value = "grey90", name = "Studies per country",
                      breaks = seq(2, 8, 2)) +
  theme_void() +
  theme(legend.position = c(0.15, 0.35),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        plot.margin = margin(5, 5, 0, 5))

# Panel B: stacked barchart of study publication timeline
ic_levels <- c("Mosquito modification (replacement)", "Mosquito modification (suppression)",
               "Insecticide-treated materials (ITMs)", "Integrated approaches", "Autodissemination",
               "Environmental management", "Larvivorous predators", "Mass trapping",
               "Mosquito repellents", "Direct larviciding", "Indoor/targeted residual spraying")

ic_colours <- c(
  "Mosquito modification (replacement)"   = "#3182bd", 
  "Mosquito modification (suppression)"   = "#9ecae1",  
  "Insecticide-treated materials (ITMs)"  = "#2ca02c",
  "Integrated approaches"                 = "#bcbe90",
  "Autodissemination"                     = "#ff7f0e",
  "Environmental management"              = "#ffbb78",
  "Larvivorous predators"                 = "#8c564b",
  "Mass trapping"                         = "#9467bd",
  "Mosquito repellents"                   = "#e377c2",
  "Direct larviciding"                    = "#c7b8ea",
  "Indoor/targeted residual spraying"     = "#7f7f7f")

bar_df <- df %>%
  distinct(study, year, intervention_class) %>%
  mutate(intervention_class = factor(intervention_class, levels = ic_levels)) %>%
  count(year, intervention_class, .drop = FALSE) %>%
  filter(n > 0)

all_years <- seq(min(bar_df$year), max(bar_df$year))

pB <- ggplot(bar_df, aes(x = year, y = n, fill = intervention_class)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_x_continuous(breaks = unique(bar_df$year)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = ic_colours, name = NULL, drop = FALSE) +
  labs(x = "Publication year", y = "Number of studies") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.35, "cm"),
        legend.position = c(0.02, 0.95),
        legend.justification = c(0, 1), 
        legend.background = element_rect(fill = "white", colour = "grey50", linewidth = 0.2),
        plot.margin = margin(5, 5, 5, 5)) +
  guides(fill = guide_legend(ncol = 2))

# combine
fig2 <- plot_grid(pA, pB, ncol = 1, rel_heights = c(1, 1),
                  labels = c("A", "B"), label_size = 14)

ggsave("fig2_map_bar.pdf", fig2, width = 10, height = 10, dpi = 300)
ggsave("fig2_map_bar.tif", fig2, width = 10, height = 10, dpi = 300)
