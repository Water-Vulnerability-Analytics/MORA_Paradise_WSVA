# Shared setup for the MORA Paradise WSVA R Markdown files
# (MORA_Deep_Dive_workflow.Rmd and Explore_MORA_Climate_Grid.Rmd).
# Each Rmd sources this from its setup chunk:
#   source(here::here("Reports/report_setup.R"))
# Paths use here::here() so this works whether it runs from the project root
# or from Reports/ (where knitr starts before root.dir takes effect).

# Load functions after park files in case functions were updated since last save.
source(here::here("setup.R"), chdir = TRUE)
library(rnoaa)

# Set mapview basemap options
mapviewOptions(fgb = FALSE, 
               georaster = FALSE, 
               basemaps = c("Esri.WorldTopoMap",
                            "Esri.WorldImagery"))

# --------------------------------------------------------------------------
# Color system (validated for WCAG contrast + colorblind-safe separation)
#
# The core categorical hues below were checked with an OKLCH-based
# validator against: lightness band, chroma floor, CVD separation
# (protanopia/deuteranopia simulation), a normal-vision separation floor,
# and contrast vs. a white page (relevant since this report renders as a
# static Word document - figures are flat images, not interactive).
#
# Two documented, deliberate exceptions to the validator's defaults:
#   1. The neutral "Historical" gray (#595959) reads slightly below the
#      validator's chroma floor - expected and accepted for a pure neutral
#      used in a structural/reference role (a baseline period, not a
#      scenario), same treatment the validator itself allows for
#      diverging-palette midpoints. Its old value (gray80) was replaced
#      because gray80 fails contrast against a white page outright (WCAG
#      violation, and nearly invisible in a printed/grayscale copy).
#   2. "Projected" orange (#a1631f) was validated against the other three
#      core colors in "adjacent" mode (the check that matches how these
#      colors are actually used here - as neighbors in a line/bar/ribbon
#      legend), not the stricter "any two colors simultaneously visible"
#      mode meant for scatter/choropleth work. Orange and the CF red never
#      appear together in the same chart in this report, so the stricter
#      check isn't the relevant one for this palette's actual use.
# --------------------------------------------------------------------------

# Non CF colors (still under review - not yet reconciled to the validated
# system below; pal$orange and pal$darkblue ARE now part of the validated
# system, see note beside plot_colors)
pal <- list(
  black = "black", 
  gold = "#dca237", 
  lightblue = "#6eb2e4", 
  teal = "#459a74", 
  yellow = "#eee461", 
  darkblue = "#1b6ca8",   # validated - now shared with nps_pal$blue / "Snowy Wet"
  orange = "#a1631f",     # validated "Projected" - replaces the old #c66426
  pink = "#c07da5",
  lightpink = "#ECD8E4"
)

# Named Colors for general 2 tone cf
# grey/lightgrey below are not yet reconciled to the validated system
# (pending the full palette sweep); blue/red are validated core colors.
nps_pal <-list(blue = "#1b6ca8",   # validated
               red = "#b02a24",    # validated
               grey = "#a9a9a9", 
               lightgrey = "#dcdcdc")


# Colors for CF plots - validated core categorical palette.
# "Observed" stays pure black (max contrast, deliberately distinct from
# "Historical" by lightness alone since both are neutral/structural roles
# rather than scenario hues).
plot_colors = c("Less Snowy Dry" = nps_pal$red, 
                "Snowy Wet" = nps_pal$blue,  
                "Observed" = "black", 
                "Historical" = "#595959")

# --------------------------------------------------------------------------
# Shared plot theme
#
# One theme function used across figures for a consistent, professional
# look: legible text sizing, a restrained gridline, legend on the bottom
# (keeps the plot area wide for the Word-doc page), and clean panel
# borders instead of a boxed-in look.
# --------------------------------------------------------------------------
theme_wsva <- function(base_size = 11, legend = "bottom") {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.1), margin = margin(b = 4)),
      plot.subtitle = element_text(color = "gray30", size = rel(0.95), margin = margin(b = 8)),
      plot.caption = element_text(color = "gray40", size = rel(0.75), hjust = 0),
      axis.title = element_text(face = "bold", size = rel(0.95)),
      axis.text = element_text(color = "gray20"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
      strip.text = element_text(face = "bold", size = rel(0.95)),
      strip.background = element_rect(fill = "gray95", color = NA),
      legend.position = legend,
      legend.title = element_text(face = "bold", size = rel(0.9)),
      legend.text = element_text(size = rel(0.9))
    )
}

# Create output directories
dir.create(here::here("data/park/MORA/Figs"), showWarnings = FALSE, recursive = TRUE)
dir.create(here::here("data/park/MORA/Tables"), showWarnings = FALSE, recursive = TRUE)


# Define conversions
gal_to_cf <- 0.133681
cf_to_gal <- 1/gal_to_cf

# Define time periods
# BOTH ends of each window are INCLUSIVE: historical baseline = 1980-2012
# (1979 is the WBM spin-up year and is excluded); future = 2035-2064 (ends just
# before 2065-01-01). Always filter with year >= x[1] & year <= x[2].
hist_yrs <- c(1980, 2012)
fut_yrs <- c(2035, 2064)

# Pick CRS for project
crs_use <- "EPSG:5070"
