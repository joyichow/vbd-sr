# Nov 01, 2026

# Panel A: geographical distributions of studies
# Panel B: stacked bar chart, timeline of study's publication year
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(rnaturalearth)
library(sf)
library(cowplot)

raw <- read_excel("28092026_vbd_S1 Table.xlsx", sheet = 1, skip = 1)

names(raw)[1] <- "study_label"
country_col <- names(raw)[4]

df <- raw %>%
  select(study = study_label, country_raw = all_of(country_col)) %>%
  filter(!is.na(study), study != "") %>%
  mutate(is_header = is.na(country_raw) | country_raw == "",
         intervention_class = ifelse(is_header, study, NA)) %>%
  fill(intervention_class) %>%
  filter(!is_header) %>%
  select(-is_header)

df$year <- as.integer(str_extract(df$study, "\\d{4}$"))

df$country_clean <- str_remove_all(df$country_raw, "\\s*\\([^)]*\\)")
df <- df %>%
  separate_rows(country_clean, sep = "[,;]\\s*") %>%
  mutate(country_clean = str_trim(country_clean)) %>%
  filter(country_clean != "")

suppression_studies <- c("de Castro Poncio, 2023", "Lim, 2023", "Lim, 2024", "Lim, 2025", "Lim, 2026")

mm <- which(df$intervention_class == "Mosquito modification")
df$intervention_class[mm] <- ifelse(df$study[mm] %in% suppression_studies,
                                    "Mosquito modification (suppression)",
                                    "Mosquito modification (replacement)")

################################################################################
################################ Fig 2A: Map ###################################
################################################################################
country_counts <- df %>%
  distinct(study, country_clean) %>%
  count(country_clean, name = "n_studies")

world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  select(name, geometry) %>%
  left_join(country_counts, by = c("name" = "country_clean"))

pA <- ggplot(world) +
  geom_sf(aes(fill = n_studies), colour = "grey70", linewidth = 0.15) +
  scale_fill_gradient(low = "#c6dbef", high = "#08306b", na.value = "grey90",
                      name = "Studies per country", breaks = seq(2, 8, 2)) +
  theme_void() +
  theme(legend.position = c(0.15, 0.35),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        plot.margin = margin(5, 5, 0, 5))

################################################################################
############################# Fig 2B: Timeline #################################
################################################################################
ic_colours <- c(
  "Mosquito modification (replacement)" = "#3182bd",
  "Mosquito modification (suppression)" = "#9ecae1",
  "Insecticide-treated materials (ITMs)" = "#2ca02c",
  "Integrated approaches" = "#bcbe90",
  "Autodissemination" = "#ff7f0e",
  "Environmental management" = "#ffbb78",
  "Larvivorous predators" = "#8c564b",
  "Mass trapping" = "#9467bd",
  "Spatial repellents" = "#e377c2",
  "Direct larviciding" = "#c7b8ea",
  "Indoor/targeted residual spraying" = "#7f7f7f")

bar_df <- df %>%
  distinct(study, year, intervention_class) %>%
  mutate(intervention_class = factor(intervention_class, levels = names(ic_colours))) %>%
  count(year, intervention_class)

pB <- ggplot(bar_df, aes(x = year, y = n, fill = intervention_class)) +
  geom_col(width = 0.7) +
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

fig2 <- plot_grid(pA, pB, ncol = 1, labels = c("A", "B"), label_size = 14)

ggsave("fig2_map_bar.pdf", fig2, width = 10, height = 10)
ggsave("fig2_map_bar.tif", fig2, width = 10, height = 10, dpi = 300)
