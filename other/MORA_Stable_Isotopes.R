# Develop LMWL

# Read in Precipitation stable isotopes
precip_isotopes <- read_csv("data/park/MORA/stable_isotopes/MORA_modeled_precip.csv")

# 1. Clean the column names and pivot into a long format
isotopes_long <- precip_isotopes %>%
  # Fix the missing underscore in Jan2H so regex pattern matches perfectly
  rename_with(~ gsub("Jan2H", "Jan_2H", .x)) %>%
  # Pivot only the isotope data columns
  pivot_longer(
    cols = contains("18O") | contains("2H"),
    names_to = c("Time_Period", ".value"),
    names_pattern = "(.*)_(18O|2H)"
  ) %>%
  # Standardize column names for downstream math/plotting
  rename(d18O = `18O`, d2H = `2H`)

# 2. Calculate slopes, intercepts, and R² for EVERY time period
isotope_metrics <- isotopes_long %>%
  group_by(Time_Period) %>%
  dplyr::summarize(
    slope     = round(coef(lm(d2H ~ d18O))[2], 2),
    intercept = round(coef(lm(d2H ~ d18O))[1], 2),
    r_squared = round(summary(lm(d2H ~ d18O))$r.squared, 3),
    .groups = 'drop'
  )

# Print the resulting equations table to your console
cat("Meteoric Water Line metrics for all time periods:\n")
print(isotope_metrics, n = Inf)

# 3. Plot all categories on a single Delta-O18 vs Deuterium plot
ggplot(isotopes_long, aes(x = d18O, y = d2H, color = Time_Period)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.8) +
  labs(
    title = "Meteoric Water Lines by Time Period",
    x = expression(paste(delta, " "^18, "O (‰, VSMOW)")),
    y = expression(paste(delta, " "^2, "H (‰, VSMOW)")),
    color = "Time Period"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.text = element_text(size = 9)
  )

# Read in isotope sample data
# Select relevant sites
# Organize labels for plotting / interpreting
isotopes <- read_csv("data/park/MORA/stable_isotopes/MORA_StableIsotopes.csv") %>%
  st_as_sf(coords = c("Long_dd", "Lat_dd"), crs = 4326) %>%
  dplyr::filter(Unique_SiteID %in%
                  c("MORA_035", # EC HW Spr
                    "MORA_036", # EC HW Spr
                    "MORA_037", # EC HW Spr
                    "MORA_047", # P HW Spr
                    "MORA_048", # P HW Spr
                    "MORA_049", # P HW Spr
                    "MORA_026", # EC
                    "MORA_046", # EC PV
                    "MORA_045", # EC PV
                    "MORA_044", # EC PV
                    "MORA_041",
                    "MORA_043",
                    "MORA_010" # Muir Peak Snow
                    )) %>%
  dplyr::mutate(
    label = case_when(
      Unique_SampleID == "MORA_010_001" ~ "Ice Field Nisqually HW",
      Unique_SiteID %in% c("MORA_035", "MORA_036", "MORA_037") ~ "Edith HW Spring",
      Unique_SampleID %in% c("MORA_047_002", "MORA_048_001") ~ "Paradise HW Spring",
      Unique_SiteID == "MORA_026" ~ "Edith @ Rd.",
      Unique_SampleID %in% c("MORA_047_001", "MORA_046_001") ~ "Edith above Paradise",
      Unique_SampleID == c("MORA_045_001") ~ "Paradise below Edith",
      Unique_SiteID %in% c("MORA_044", "MORA_041") ~ "Paradise Valley Springs",
      Unique_SiteID == "MORA_043" ~ "Paradise Lower Valley"
      )
    ) 

mapview(isotopes)

# Select period for LMWL
# Time_Period   slope intercept r_squared
# AnnualAverage  7.49      0.32     0.996
# Apr            7.59      2.36     0.996
# Aug            6.81     -5.02     0.939
# Dec            7.29      1.14     0.997
# Feb            8.13      6.4      0.997
# Jan            7.3      -0.45     0.993
# Jul            4.79    -26.7      0.969
# Jun            5.12    -21.4      0.823
# Mar            7.86      4.12     0.998
# May            7.28     -3.86     0.985
# Nov            7.74      4.02     0.998
# Oct            7.47      2.72     0.99 
# Sep            6.3      -5.4      0.876
# Summer         5.49    -16.6      0.923
# Winter         7.62      2.81     0.997


sel_lmwl <- isotope_metrics %>%
  dplyr::filter(Time_Period == "Dec")

# Some Meteoric Water Lines to Consider
# 1. GMWL: Global meteoric water line: slope = 8, intercept = 10
# 2. LMWL from modeled data above.
# 3. (Not plotted) LMWL for North Cascades in Oregon: slope = 8.3, intercept = 17 (Jefferson et al., 2006 https://doi.org/10.1029/2005WR004812)

# Plot selected isotopes along with MWLs
ggplot(isotopes, aes(x = delta_18O, y = delta_2H)) +
  # Add h and vline showing max in precip for contextualizing actual variability
  # in sample records
  #geom_hline(yintercept = max(isotopes_long$d2H), color = "gray80") +
  #geom_vline(xintercept = max(isotopes_long$d18O), color = "gray80") +
  #annotate("text", x = max(isotopes_long$d18O)-.25, y = -100, 
  #         label = "Max in Precip.", angle = 90, color = "gray40") +
  #annotate("text", y = max(isotopes_long$d2H)-2, x = -13, 
  #         label = "Max in Precip.", angle = 0, color = "gray40") +
  # Add meteoric water lines
  geom_abline(intercept = 10, slope = 8,
                            color = "grey50", linetype = "dashed", linewidth = 0.7) +
  geom_abline(intercept = sel_lmwl$intercept, slope = sel_lmwl$slope,
              color = "grey20", linetype = "solid", linewidth = 0.7) +
  annotate("text", x = -15, y = -111, label = "GMWL",
           color = "grey50", fontface = "italic", hjust = 0, angle = 30) +
  annotate("text", x = -15, y = -106, label = "LMWL",
           color = "grey20", fontface = "italic", hjust = 0, angle = 30) +
  # Add sample data
  geom_point(aes(color = label, shape = Sample_Type),
                          size = 3.5, alpha = 0.85) +
  
  scale_color_brewer(palette = "Accent") +
  labs(
   x = expression(delta^{18}*O ~ "(\u2030 vs. VSMOW)"),
    y = expression(delta^{2}*H ~ "(\u2030 vs. VSMOW)"),
    color = "Location",
    shape = "Sample Type",
    caption = "Observations: Edith and Paradise headwater springs have similar "
    ) +
  theme_minimal(base_size = 12) +
  theme(
   panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
   axis.line    = element_blank(),
   panel.grid.minor = element_blank(),
    legend.position  = "right"
  )

