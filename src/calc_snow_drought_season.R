#' Calculate Snow Drought Status Across a Full Season
#'
#' Classifies snow drought day-by-day across the cool season, for one or
#' more groups (e.g. SNOTEL stations, or GCM/climate futures), following
#' the Hatchett, Rhoades & McEvoy (2022) framework (note, NIDIS uses
#' 20% as defining snow drought (similar to other drought vars)). Builds a
#' smoothed day-of-year percentile climatology (a +/- `window`-day moving
#' window) per group, then classifies every day in the record against it -
#' not just a single date or a handful of sampled dates.
#'
#' Works for observed station data (e.g. SNOTEL - group by station name,
#' each station's own full record is its climatology) or modeled climate
#' projections (e.g. WBM output - group by GCM, with `climatology_filter`
#' identifying the historical run so future scenarios are scored against
#' historical conditions rather than building their own climatology).
#'
#' @param df A data frame with at least date, SWE, and precip columns.
#' @param date_col,swe_col,precip_col Column names (as strings) for date,
#'   snow water equivalent, and water-year-to-date cumulative precip.
#'   Defaults match this project's existing naming (`date`, `swe_in`,
#'   `precip_in`).
#' @param group_col Column name (as a string) to compute separate
#'   climatologies for (e.g. `"name"` for stations, `"GCM"` for climate
#'   projections). Default `NULL` treats `df` as a single series.
#' @param climatology_filter Optional logical vector, same length as `df`,
#'   TRUE for rows that should be used to build the climatology (e.g.
#'   `df$CF == "Historical"`). All rows (including FALSE ones) are still
#'   classified against that climatology. Default `NULL` uses every row
#'   within its own group as its own climatology - appropriate for
#'   observed station records.
#' @param cool_season_months Months to classify. Default Nov-May, when
#'   snow drought is meaningful.
#' @param min_years Minimum years of climatology history desired per group
#'   before treating percentiles as reliable. NIDIS convention wants 20+.
#' @param window Number of days on either side of a given day-of-year to
#'   pool together when building the percentile climatology. Default 3
#'   (i.e. a 7-day window).
#'
#' @return A tibble with one row per date x group across the full cool
#'   season (not a single date or sampled weeks), with SWE/precip
#'   percentiles and `drought_type` classification.
calc_snow_drought_season <- function(df,
                                     date_col = "date",
                                     swe_col = "swe_in",
                                     precip_col = "precip_in",
                                     group_col = NULL,
                                     climatology_filter = NULL,
                                     cool_season_months = c(11, 12, 1, 2, 3, 4, 5),
                                     min_years = 20,
                                     window = 3) {

  df <- df %>%
    dplyr::rename(date = dplyr::all_of(date_col),
                  swe_in = dplyr::all_of(swe_col),
                  precip_in = dplyr::all_of(precip_col))

  df <- if (is.null(group_col)) {
    df %>% dplyr::mutate(.group = "all")
  } else {
    df %>% dplyr::mutate(.group = .data[[group_col]])
  }

  if (is.null(climatology_filter)) climatology_filter <- rep(TRUE, nrow(df))

  df <- df %>%
    dplyr::mutate(.in_climatology = climatology_filter,
                 water_year = calc_water_year(date),
                 month_day = format(date, "%m-%d")) %>%
    dplyr::filter(lubridate::month(date) %in% cool_season_months,
                 !is.na(swe_in), !is.na(precip_in))

  # Check years of climatology history per group
  year_counts <- df %>%
    dplyr::filter(.in_climatology) %>%
    dplyr::group_by(.group) %>%
    dplyr::summarise(n_years = dplyr::n_distinct(water_year), .groups = "drop")

  thin_groups <- year_counts %>% dplyr::filter(n_years < min_years)
  if (nrow(thin_groups) > 0) {
    warning(paste0("Fewer than ", min_years, " years of climatology history for: ",
                   paste(thin_groups$.group, collapse = ", ")))
  }

  # +/- window-day pool of month_days for a given focal day, anchored to a
  # fixed reference year so wraparound (e.g. Dec -> Jan) works
  get_window_mds <- function(focal_md, all_mds) {
    # turns the "MM-DD" string into a real Date object by 
    # anchoring it to an arbitrary placeholder year (2000, chosen 
    # specifically because it's a leap year, so "02-29" doesn't break). 
    focal_date <- as.Date(paste0("2000-", focal_md))
    window_dates <- seq(focal_date - window, focal_date + window, by = "day")
    window_mds <- format(window_dates, "%m-%d")
    window_mds[window_mds %in% all_mds]
  }

  # Percentile breakpoints for one variable, per group, from climatology rows only
  build_climatology <- function(clim_df, variable) {
    all_md <- clim_df %>% dplyr::distinct(month_day) %>% dplyr::arrange(month_day) %>% dplyr::pull(month_day)

    clim_df %>%
      dplyr::group_split(.group) %>%
      purrr::map(function(grp_df) {
        grp_name <- grp_df$.group[1]
        purrr::map(all_md, function(md) {
          window_mds <- get_window_mds(md, all_md)
          pool <- grp_df %>% dplyr::filter(month_day %in% window_mds) %>% dplyr::pull({{ variable }})
          pool <- pool[!is.na(pool)]
          if (length(pool) < 5) return(NULL)
          tibble::tibble(
            .group = grp_name,
            month_day = md,
            p10 = quantile(pool, 0.10), p20 = quantile(pool, 0.20),
            p30 = quantile(pool, 0.30), p50 = quantile(pool, 0.50),
            p70 = quantile(pool, 0.70), p80 = quantile(pool, 0.80),
            p90 = quantile(pool, 0.90)
          )
        }) %>%
          purrr::list_rbind()
      }) %>%
      purrr::list_rbind()
  }

  climatology_rows <- df %>% dplyr::filter(.in_climatology)

  clim_swe <- build_climatology(climatology_rows, swe_in)
  clim_precip <- build_climatology(climatology_rows, precip_in)

  interpolate_pctile <- function(x, p10, p20, p30, p50, p70, p80, p90) {
    breaks <- c(p10, p20, p30, p50, p70, p80, p90)
    probs  <- c(0.10, 0.20, 0.30, 0.50, 0.70, 0.80, 0.90)

    if (is.na(x)) return(NA_real_)
    if (x <= breaks[1]) return(0)
    if (x >= breaks[7]) return(100)

    idx <- findInterval(x, breaks)
    lo <- breaks[idx]; hi <- breaks[idx + 1]
    plo <- probs[idx]; phi <- probs[idx + 1]

    round((plo + (x - lo) / (hi - lo) * (phi - plo)) * 100, 1)
  }

  # Climatology is built from climatology_filter rows only, but applied to
  # every row in the group - so future/non-climatology rows still get scored
  result <- df %>%
    dplyr::left_join(clim_swe, by = c(".group", "month_day"), suffix = c("", "_swe")) %>%
    dplyr::left_join(clim_precip, by = c(".group", "month_day"), suffix = c("", "_prc")) %>%
    dplyr::mutate(
      swe_pctile = purrr::pmap_dbl(list(swe_in, p10, p20, p30, p50, p70, p80, p90), interpolate_pctile),
      precip_pctile = purrr::pmap_dbl(list(precip_in, p10_prc, p20_prc, p30_prc, p50_prc, p70_prc, p80_prc, p90_prc), interpolate_pctile),
      drought_type = dplyr::case_when(
        is.na(swe_pctile) | is.na(precip_pctile) ~ NA_character_,
        swe_pctile >= 30 ~ "No drought",
        swe_pctile < 30 & precip_pctile < 30 ~ "Dry",
        swe_pctile < 30 & precip_pctile >= 50 ~ "Warm",
        swe_pctile < 30 & precip_pctile >= 30 & precip_pctile < 50 ~ "Warm-Dry",
        TRUE ~ NA_character_
      ),
      drought_type = factor(drought_type, levels = c("No drought", "Dry", "Warm-Dry", "Warm"))
    )

  if (!is.null(group_col)) {
    # group_col's original column was never touched (.group was just a copy
    # of it), so select it directly instead of renaming .group back - avoids
    # a duplicate-name collision when group_col already exists in df
    result <- result %>%
      dplyr::select(date, dplyr::all_of(group_col), water_year, month_day,
                    swe_in, precip_in, swe_pctile, precip_pctile, drought_type)
  } else {
    result <- result %>%
      dplyr::select(date, water_year, month_day,
                    swe_in, precip_in, swe_pctile, precip_pctile, drought_type)
  }

  result
}
