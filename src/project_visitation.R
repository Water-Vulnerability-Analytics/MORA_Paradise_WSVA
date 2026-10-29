# project_visitation() is a thin dispatcher over named visitation-projection
# METHODS (see .visitation_methods at the bottom of this file). Five methods
# are implemented, ordered here from simplest to most complex: "linear_trend"
# (simple/robust linear extrapolation, optionally climate-adjusted),
# "annual_growth_rate" (a simple fixed compounding annual growth rate, no
# statistical fitting or climate dependence -- included because it mirrors
# how parks have historically hand-estimated visitation growth), "gam_logistic"
# (logistic population cap + annual GAM on climate predictors), "arima_climate"
# (ARIMA forecast, optionally with climate variables as external regressors),
# and "fisichelli" (Fisichelli et al. 2015 monthly temperature-response
# method with capacity ceilings). All five
# share the same (df, ...) input contract and the same output shape (a list
# with at least $annual / $monthly / $result, each with visitors_low/mid/high
# + source columns), so they're interchangeable from a caller's perspective --
# pick whichever fits the park's data.
#
# All five also share several internal helper functions (prefixed with `.`)
# for input validation, annual aggregation, monthly-proportion calculation,
# and monthly disaggregation, since that scaffolding is identical regardless
# of which model produces the underlying annual projection.
#
# NONE of this has been run against live data -- no R runtime was available
# while writing it. Test each method against a known case before trusting it
# for a report, and watch for model-fitting warnings (nls() convergence,
# auto.arima() order-selection issues, etc.) -- those are the most likely
# things to need per-park tuning.

# ============================================================================
# ---- Shared internal helpers ----------------------------------------------
# ============================================================================

#' @keywords internal
.prep_vis_input <- function(df, date_col, visitor_col, climate_vars, scenario_col,
                            group_cols, historical_label, future_end_year) {
  
  if (is.null(climate_vars)) climate_vars <- character(0)
  
  required_cols <- unique(c(date_col, visitor_col, climate_vars, scenario_col, group_cols))
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("df is missing required column(s): ", paste(missing_cols, collapse = ", "))
  }
  
  # ungroup() defensively -- if the caller's df still carries a group_by()
  # from however vis_clim_df/wbm_monthly was assembled (e.g. a grouped
  # mutate() that was never piped through ungroup()), every dplyr verb
  # downstream (including plain mutate() calls in the method bodies, not
  # just group_by()-using helpers) evaluates per-group instead of over the
  # whole table. That's silent almost everywhere in this file EXCEPT when a
  # mutate() column comes from predict()-ing on a newdata frame built from
  # the ungrouped row count -- there the per-group recycling fails loudly
  # ("`pred_raw` must be size 1, not 1800"), which is what surfaces the
  # problem. Stripping any incoming grouping here fixes it at the source for
  # all four methods, not just the one where it happens to throw.
  df <- df %>%
    dplyr::ungroup() %>%
    dplyr::rename(.ym = dplyr::all_of(date_col),
                  .visitors = dplyr::all_of(visitor_col)) %>%
    dplyr::mutate(.visitors = ifelse(.visitors == 0, NA, .visitors),
                  year  = lubridate::year(.ym),
                  month = lubridate::month(.ym))
  
  # Soft check that df is actually monthly (a hard requirement of the shared
  # contract) -- warns rather than errors, since a few legitimate gaps
  # shouldn't necessarily block a run.
  hist_dates <- sort(unique(df$.ym[df[[scenario_col]] == historical_label]))
  if (length(hist_dates) > 1 && any(diff(as.numeric(hist_dates)) > 45)) {
    warning("Historical rows have gaps > 45 days between observations -- ",
            "confirm `df` is a monthly time series as expected.")
  }
  
  list(
    df = df,
    climate_vars = climate_vars,
    hist_monthly = df %>% dplyr::filter(.data[[scenario_col]] == historical_label),
    fut_monthly  = df %>% dplyr::filter(.data[[scenario_col]] != historical_label,
                                        year <= future_end_year)
  )
}

#' @keywords internal
.annual_visitors <- function(hist_monthly, visitor_agg_fun) {
  hist_monthly %>%
    dplyr::group_by(year) %>%
    dplyr::summarise(.visitors = visitor_agg_fun(.visitors, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(.visitors = ifelse(.visitors == 0, NA, .visitors))
}

#' @keywords internal
.monthly_proportions <- function(hist_monthly, exclude_years) {
  hist_monthly %>%
    dplyr::filter(!year %in% exclude_years) %>%
    dplyr::group_by(month) %>%
    dplyr::summarize(mean_vis = mean(.visitors, na.rm = TRUE),
                     max_vis  = max(.visitors,  na.rm = TRUE),
                     min_vis  = min(.visitors,  na.rm = TRUE),
                     sd_vis   = sd(.visitors,   na.rm = TRUE),
                     cv_vis   = sd_vis / mean_vis,
                     .groups = "drop") %>%
    dplyr::mutate(month_date = lubridate::as_date(paste0("2000-", month, "-01")),
                  fraction = mean_vis / sum(mean_vis))
}

#' @keywords internal
# Monthly proportions computed separately for one or more named year ranges
# (e.g. a historical baseline period vs. a future projection window),
# rather than collapsed into a single fraction across the whole series the
# way .monthly_proportions() is. `dat` must have `year`, `month`, and a
# numeric `vis` column -- see .proportions_source() below, which builds that
# table from a method's historical + projected monthly output.
# `summary_years` is a named list of numeric year vectors, e.g.
# list(historical = 1980:2012, future = 2035:2064) -- the names become the
# `period` column in the returned (stacked, not nested) tibble.
.monthly_proportions_by_period <- function(dat, summary_years) {

  if (is.null(summary_years)) return(NULL)
  if (!is.list(summary_years) || is.null(names(summary_years)) ||
      any(names(summary_years) == "")) {
    stop("`summary_years` must be a named list of year vectors, e.g. ",
         "list(historical = 1980:2012, future = 2035:2064).")
  }

  purrr::imap_dfr(summary_years, function(yrs, period_name) {
    period_dat <- dat %>% dplyr::filter(year %in% yrs)

    if (nrow(period_dat) == 0 || all(is.na(period_dat$vis))) {
      warning("summary_years period '", period_name, "' (years ", min(yrs), "-",
              max(yrs), ") has no visitation data in that range -- its monthly ",
              "proportions will be NA. Check that this method's historical/",
              "projected years actually cover that period.")
    }

    period_dat %>%
      dplyr::group_by(month) %>%
      dplyr::summarize(mean_vis = mean(vis, na.rm = TRUE),
                       max_vis  = max(vis,  na.rm = TRUE),
                       min_vis  = min(vis,  na.rm = TRUE),
                       sd_vis   = sd(vis,   na.rm = TRUE),
                       cv_vis   = sd_vis / mean_vis,
                       .groups = "drop") %>%
      tidyr::complete(month = 1:12) %>%
      dplyr::mutate(period = period_name,
                    month_date = lubridate::as_date(paste0("2000-", month, "-01")),
                    fraction = mean_vis / sum(mean_vis, na.rm = TRUE)) %>%
      dplyr::relocate(period, month, month_date)
  })
}

#' @keywords internal
# Build the table .monthly_proportions_by_period() summarizes: true observed
# monthly visitation for historical months (year/month/.visitors, straight
# from .prep_vis_input()'s hist_monthly), and the projected monthly envelope
# (visitors_mid) for months at/after obs_cutoff_year. Kept separate from
# .disaggregate_monthly()'s output because that function re-derives even
# OBSERVED months as fraction * annual total (so the whole series is built
# the same way) -- for a seasonality summary we want the real observed
# monthly values where they exist, not a fraction-smoothed reconstruction of
# them. `exclude_years` (the same argument every method already takes) is
# applied only to the historical portion, matching how it's used everywhere
# else in this file -- the projected portion has nothing to exclude by
# construction.
.proportions_source <- function(hist_monthly, monthly_result, exclude_years = NULL) {
  hist_part <- hist_monthly
  if (!is.null(exclude_years)) {
    hist_part <- hist_part %>% dplyr::filter(!year %in% exclude_years)
  }

  dplyr::bind_rows(
    hist_part %>% dplyr::transmute(year, month, vis = .visitors),
    monthly_result %>%
      dplyr::filter(source == "projected") %>%
      dplyr::transmute(year = lubridate::year(ym), month = lubridate::month(ym),
                       vis = visitors_mid)
  )
}

#' @keywords internal
# Aggregate monthly climate_vars to annual using climate_agg_fun (applied
# uniformly to every variable), grouped by whatever grouping columns are
# passed in. With climate_vars = character(0), returns just the grouping
# columns (used by the pure trend/ARIMA methods when run without climate).
.annual_climate_agg <- function(data, climate_vars, climate_agg_fun, group_vars) {
  data %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::summarise(dplyr::across(dplyr::all_of(climate_vars),
                                   ~ climate_agg_fun(.x, na.rm = TRUE)),
                     .groups = "drop")
}

#' @keywords internal
# Combine observed annual visitation (before obs_cutoff_year) with the
# projected low/mid/high envelope (at/after obs_cutoff_year) into one table.
.assemble_annual_result <- function(annual_hist_visitors, annual_envelope, obs_cutoff_year) {
  dplyr::bind_rows(
    annual_hist_visitors %>%
      dplyr::filter(year < obs_cutoff_year) %>%
      dplyr::transmute(year, visitors_low = .visitors, visitors_mid = .visitors,
                       visitors_high = .visitors, source = "observed"),
    annual_envelope %>%
      dplyr::filter(year >= obs_cutoff_year) %>%
      dplyr::transmute(year, visitors_low, visitors_mid, visitors_high, source = "projected")
  ) %>%
    dplyr::arrange(year)
}

#' @keywords internal
# Disaggregate annual observed + projected visitation to monthly using
# historical monthly proportions. Observed monthly values are used before
# obs_cutoff_year; the annual envelope (scaled by monthly fraction) is used
# after it.
.disaggregate_monthly <- function(min_date, annual_hist_visitors, annual_envelope,
                                  monthly_proportions, obs_cutoff_year, future_end_year) {
  
  monthly_scaffold <- data.frame(
    .ym = seq.Date(min_date, lubridate::as_date(paste0(future_end_year, "-12-01")), by = "month")
  ) %>%
    dplyr::mutate(year = lubridate::year(.ym), month = lubridate::month(.ym))
  
  # NOTE: annual_hist_visitors contributes .visitors (observed) and
  # annual_envelope contributes visitors_low/mid/high (projected) -- no name
  # collision between the two joins, so no suffix is needed. The mutate()
  # below is self-referencing (RHS reads the just-joined annual value, LHS
  # overwrites it with the blended monthly value) -- standard dplyr, not a bug.
  monthly_scaffold %>%
    dplyr::left_join(annual_hist_visitors, by = "year") %>%
    dplyr::left_join(annual_envelope %>% dplyr::select(year, visitors_low, visitors_mid, visitors_high),
                     by = "year") %>%
    dplyr::left_join(monthly_proportions %>% dplyr::select(month, fraction), by = "month") %>%
    dplyr::mutate(
      visitors_mid  = ifelse(year < obs_cutoff_year, fraction * .visitors, fraction * visitors_mid),
      visitors_high = ifelse(year < obs_cutoff_year, fraction * .visitors, fraction * visitors_high),
      visitors_low  = ifelse(year < obs_cutoff_year, fraction * .visitors, fraction * visitors_low),
      source = ifelse(year < obs_cutoff_year, "observed", "projected")
    ) %>%
    dplyr::rename(ym = .ym) %>%
    dplyr::select(ym, visitors_low, visitors_mid, visitors_high, source)
}

# ============================================================================
# ---- Method 1: linear_trend -------------------------------------------------
# ============================================================================

#' Project future park visitation via (optionally climate-adjusted) linear trend
#'
#' The simplest method in this family: fits `lm(visitors ~ year)` (or
#' `lm(visitors ~ year + climate_vars)` if `climate_vars` is supplied) on
#' annual historical visitation, then extrapolates using prediction
#' intervals for the low/high bounds. No population-growth cap -- a pure
#' linear extrapolation will eventually produce implausible values for long
#' horizons, so this method is best suited to shorter `future_end_year`
#' horizons or as a simple baseline/sanity-check against `gam_logistic`.
#'
#' @inheritParams gam_logistic
#' @param climate_vars Optional character vector of monthly climate columns
#'   to include as additional linear predictors. Default NULL (pure
#'   `visitors ~ year` trend).
#' @param exclude_years_lm Integer vector of years to exclude from the `lm()`
#'   fit specifically (in addition to `exclude_years`, which still governs
#'   monthly proportions). Default NULL.
#' @param prediction_level Confidence level for the prediction interval used
#'   as the low/high bounds. Default 0.95.
#'
#' @return Same shape as `gam_logistic()`, with `$lm_model` in place of
#'   `$gam_model`/`$logistic_models`/`$gam_diagnostics` (no growth cap or GAM
#'   diagnostics to report for this method). `$monthly_proportions_summary`
#'   is included whenever `summary_years` is supplied (see `gam_logistic()`).
#'
#' @details Untested against live data -- see the note at the top of this
#'   file.
#'
#' @examples
#' \dontrun{
#' vis_proj <- linear_trend(
#'   df          = zion_visitation_climate,
#'   visitor_col = "TotalVisitors"
#' )
#' vis_proj_climate <- linear_trend(
#'   df           = zion_visitation_climate,
#'   visitor_col  = "TotalVisitors",
#'   climate_vars = c("tmean_F")
#' )
#' }
linear_trend <- function(df,
                         visitor_col,
                         climate_vars = NULL,
                         date_col = "ym",
                         scenario_col = "CF",
                         historical_label = "Historical",
                         group_cols = c("GCM", "CF"),
                         climate_agg_fun = mean,
                         visitor_agg_fun = sum,
                         exclude_years = 2020,
                         exclude_years_lm = NULL,
                         prediction_level = 0.95,
                         obs_cutoff_year = as.numeric(format(Sys.Date(), "%Y")),
                         future_end_year = 2100,
                         output_resolution = c("monthly", "annual"),
                         summary_years = list(historical = 1980:2012, future = 2035:2064)) {
  
  output_resolution <- match.arg(output_resolution)
  prep <- .prep_vis_input(df, date_col, visitor_col, climate_vars, scenario_col,
                          group_cols, historical_label, future_end_year)
  climate_vars <- prep$climate_vars
  hist_monthly <- prep$hist_monthly
  fut_monthly  <- prep$fut_monthly
  has_climate  <- length(climate_vars) > 0
  
  annual_hist_visitors <- .annual_visitors(hist_monthly, visitor_agg_fun)
  monthly_proportions  <- .monthly_proportions(hist_monthly, exclude_years)
  
  ## ---- Fit lm() on historical annual data ------------------------------------
  if (has_climate) {
    annual_hist_climate <- .annual_climate_agg(hist_monthly, climate_vars, climate_agg_fun, group_vars = "year")
    lm_dat <- annual_hist_visitors %>% dplyr::left_join(annual_hist_climate, by = "year")
  } else {
    lm_dat <- annual_hist_visitors
  }
  
  if (!is.null(exclude_years_lm)) {
    lm_dat <- lm_dat %>% dplyr::filter(!year %in% exclude_years_lm)
  }
  
  lm_formula <- reformulate(c("year", climate_vars), response = ".visitors")
  lm_model   <- lm(lm_formula, data = lm_dat)
  
  ## ---- Project forward --------------------------------------------------------
  if (has_climate) {
    newdata <- .annual_climate_agg(fut_monthly, climate_vars, climate_agg_fun,
                                   group_vars = unique(c("year", group_cols))) %>%
      dplyr::filter(year <= future_end_year)
  } else {
    newdata <- data.frame(year = (max(annual_hist_visitors$year, na.rm = TRUE) + 1):future_end_year)
  }
  
  preds <- predict(lm_model, newdata = newdata, interval = "prediction", level = prediction_level)
  
  annual_proj <- newdata %>%
    dplyr::mutate(vis_pred      = pmax(0, preds[, "fit"]),
                  vis_pred_low  = pmax(0, preds[, "lwr"]),
                  vis_pred_high = preds[, "upr"])
  
  annual_envelope <- annual_proj %>%
    dplyr::group_by(year) %>%
    dplyr::summarize(visitors_mid  = mean(vis_pred,      na.rm = TRUE),
                     visitors_high = max(vis_pred_high,  na.rm = TRUE),
                     visitors_low  = min(vis_pred_low,   na.rm = TRUE),
                     .groups = "drop")
  
  annual_result <- .assemble_annual_result(annual_hist_visitors, annual_envelope, obs_cutoff_year)
  
  out <- list(
    annual = annual_result,
    monthly_proportions = monthly_proportions,
    lm_model = lm_model
  )
  
  need_monthly <- output_resolution == "monthly" || !is.null(summary_years)

  if (need_monthly) {
    monthly_result <- .disaggregate_monthly(min(prep$df$.ym), annual_hist_visitors, annual_envelope,
                                            monthly_proportions, obs_cutoff_year, future_end_year)
  }

  if (output_resolution == "monthly") {
    out$monthly <- monthly_result
    out$result  <- monthly_result
  } else {
    out$result <- annual_result
  }

  if (!is.null(summary_years)) {
    out$monthly_proportions_summary <- .monthly_proportions_by_period(
      .proportions_source(hist_monthly, monthly_result, exclude_years),
      summary_years
    )
  }
  
  out
}

# ============================================================================
# ---- Method 2: annual_growth_rate -------------------------------------------
# ============================================================================

#' Project future park visitation via a simple fixed annual growth rate
#'
#' The simplest method in this family, included because it mirrors how many
#' parks have historically estimated future visitation by hand: annual
#' visitation is assumed to grow (or shrink) by a constant percentage each
#' year, compounded off of the most recent observed annual total. There is
#' no statistical fitting, no climate dependence, and no population ceiling
#' -- it is a fixed compounding rate, nothing more:
#' `visitors(year) = base_visitors * (1 + rate) ^ (year - base_year)`.
#'
#' Because the growth rate does not depend on climate, the projection is
#' IDENTICAL across every climate scenario/GCM in `df` -- this method is
#' climate-blind by design. It's most useful as a transparent, easy-to-
#' explain baseline against which to compare the climate-responsive methods
#' (`gam_logistic`, `arima_climate`, `fisichelli`) or the pure statistical
#' trend (`linear_trend`) -- e.g. "how does the park's traditional 5%/year
#' assumption compare to what a climate-informed model projects?"
#'
#' @inheritParams linear_trend
#' @param growth_rate_mid Numeric. Fixed fractional annual growth rate (e.g.
#'   `0.05` for 5%/year growth) applied to the central ("mid") projection.
#'   Required -- there is no sensible default for a park-specific
#'   assumption like this.
#' @param growth_rate_low,growth_rate_high Optional numeric fractional
#'   annual growth rates for the low/high bounds. Default NULL, which uses
#'   `growth_rate_mid` for both (i.e. no band -- low = mid = high, matching
#'   how this method is typically used as a single point assumption). Set
#'   these to show sensitivity to the assumed rate (e.g. 0.03/0.07 around a
#'   0.05 mid-case).
#' @param base_year Integer or NULL. The historical year whose observed
#'   annual visitation is used as the compounding base. Default NULL, which
#'   uses the most recent historical year with valid (non-NA, non-excluded)
#'   annual visitation.
#'
#' @return Same shape as `linear_trend()`; there is no model object to
#'   return since nothing is fit. `$growth_rates` reports the low/mid/high
#'   rates used, and `$base_year`/`$base_visitors` report the compounding
#'   base actually selected. `$monthly_proportions_summary` is included
#'   whenever `summary_years` is supplied (see `gam_logistic()`).
#'
#' @details Untested against live data -- see the note at the top of this
#'   file. This method ignores `climate_vars`/`temp_var` entirely (it has
#'   none) -- it accepts the same `(df, ...)` input contract purely so it
#'   can be swapped in via `project_visitation(method =
#'   "annual_growth_rate")` alongside the other methods, not because
#'   climate data is actually used.
#'
#' @examples
#' \dontrun{
#' vis_proj <- annual_growth_rate(
#'   df               = zion_visitation_climate,
#'   visitor_col      = "TotalVisitors",
#'   growth_rate_mid  = 0.05,
#'   growth_rate_low  = 0.03,
#'   growth_rate_high = 0.07
#' )
#' vis_proj$growth_rates
#' vis_proj$result
#' }
annual_growth_rate <- function(df,
                               visitor_col,
                               date_col = "ym",
                               scenario_col = "CF",
                               historical_label = "Historical",
                               group_cols = c("GCM", "CF"),
                               growth_rate_mid,
                               growth_rate_low = NULL,
                               growth_rate_high = NULL,
                               base_year = NULL,
                               visitor_agg_fun = sum,
                               exclude_years = NULL,
                               obs_cutoff_year = as.numeric(format(Sys.Date(), "%Y")),
                               future_end_year = 2100,
                               output_resolution = c("monthly", "annual"),
                               summary_years = list(historical = 1980:2012, future = 2035:2064)) {

  if (missing(growth_rate_mid) || is.null(growth_rate_mid)) {
    stop("annual_growth_rate(): `growth_rate_mid` is required (e.g. 0.05 for 5%/year growth).")
  }
  if (is.null(growth_rate_low))  growth_rate_low  <- growth_rate_mid
  if (is.null(growth_rate_high)) growth_rate_high <- growth_rate_mid

  output_resolution <- match.arg(output_resolution)
  # climate_vars = NULL -- this method has no climate inputs, but still
  # routes through the same .prep_vis_input() contract as every other
  # method so it can be swapped in via project_visitation(method = ...)
  # without the caller needing to change how `df` is built.
  prep <- .prep_vis_input(df, date_col, visitor_col, climate_vars = NULL, scenario_col,
                          group_cols, historical_label, future_end_year)
  hist_monthly <- prep$hist_monthly

  annual_hist_visitors <- .annual_visitors(hist_monthly, visitor_agg_fun)
  monthly_proportions  <- .monthly_proportions(hist_monthly, exclude_years)

  fit_dat <- annual_hist_visitors %>%
    dplyr::filter(!year %in% exclude_years) %>%
    tidyr::drop_na(.visitors, year)

  if (nrow(fit_dat) == 0) {
    stop("annual_growth_rate(): no valid (non-NA, non-excluded) annual historical ",
         "observations available to use as a compounding base -- check `df` for years ",
         "with no visitor data.")
  }

  if (is.null(base_year)) {
    base_row <- fit_dat %>% dplyr::filter(year == max(year))
  } else {
    base_row <- fit_dat %>% dplyr::filter(year == base_year)
    if (nrow(base_row) == 0) {
      stop("annual_growth_rate(): `base_year` = ", base_year, " has no valid annual ",
           "historical visitation in `df` (after removing NA/excluded years).")
    }
  }
  base_year_used <- base_row$year[1]
  base_visitors  <- base_row$.visitors[1]

  if (base_year_used >= future_end_year) {
    stop("annual_growth_rate(): the compounding base year (", base_year_used, ") is not ",
         "before `future_end_year` (", future_end_year, ") -- `", base_year_used + 1, ":",
         future_end_year, "` would produce a descending/invalid sequence. Check `base_year` ",
         "and `future_end_year`.")
  }

  future_years <- (base_year_used + 1):future_end_year

  # Climate-blind by design -- the same compounding curve applies to every
  # future year regardless of climate scenario/GCM, so (unlike the other
  # four methods) there is no cross-scenario min/max envelope here; the
  # low/high bounds come purely from growth_rate_low/growth_rate_high.
  annual_envelope <- data.frame(year = future_years) %>%
    dplyr::mutate(
      raw_mid  = base_visitors * (1 + growth_rate_mid)  ^ (year - base_year_used),
      raw_low  = base_visitors * (1 + growth_rate_low)  ^ (year - base_year_used),
      raw_high = base_visitors * (1 + growth_rate_high) ^ (year - base_year_used)
    ) %>%
    dplyr::transmute(
      year,
      visitors_mid  = raw_mid,
      visitors_low  = pmin(raw_low, raw_high),
      visitors_high = pmax(raw_low, raw_high)
    )

  annual_result <- .assemble_annual_result(annual_hist_visitors, annual_envelope, obs_cutoff_year)

  out <- list(
    annual = annual_result,
    monthly_proportions = monthly_proportions,
    growth_rates = c(low = growth_rate_low, mid = growth_rate_mid, high = growth_rate_high),
    base_year = base_year_used,
    base_visitors = base_visitors
  )

  need_monthly <- output_resolution == "monthly" || !is.null(summary_years)

  if (need_monthly) {
    monthly_result <- .disaggregate_monthly(min(prep$df$.ym), annual_hist_visitors, annual_envelope,
                                            monthly_proportions, obs_cutoff_year, future_end_year)
  }

  if (output_resolution == "monthly") {
    out$monthly <- monthly_result
    out$result  <- monthly_result
  } else {
    out$result <- annual_result
  }

  if (!is.null(summary_years)) {
    out$monthly_proportions_summary <- .monthly_proportions_by_period(
      .proportions_source(hist_monthly, monthly_result, exclude_years),
      summary_years
    )
  }
  
  out
}

# ============================================================================
# ---- Method 3: gam_logistic -------------------------------------------------
# ============================================================================

#' Project future park visitation via logistic growth + climate GAM
#'
#' Fits low/mid/high logistic population-growth curves (capped at a
#' user-specified % growth above the historical maximum), fits an annual GAM
#' relating visitation to the logistic population term and any number of
#' climate predictors, projects that GAM across future climate
#' scenarios/GCMs, smooths the resulting envelope with logistic curves, and
#' (if requested) disaggregates to monthly using historical monthly
#' proportions.
#'
#' The population ceiling (`growth_pct`/`low_pct`/`high_pct`) is a genuine
#' hard cap on the final output, not just an input to the GAM: the
#' climate-driven GAM projection is smoothed by re-fitting a logistic curve
#' with K FIXED at K_low/K_mid/K_high, and that functional form is
#' mathematically bounded by K for every year -- so however much the climate
#' terms push the raw GAM prediction, `$result` can never exceed the
#' population ceiling you set. A warning is issued if the raw GAM projection
#' exceeds K in any year, since that means the final smoothed curve is
#' saturating rather than tracking the climate signal.
#'
#' @param df A single MONTHLY time series data frame containing both
#'   visitation and climate data already joined at the same (monthly) time
#'   step. One row per month for the historical period; for the future
#'   period, one row per month per climate scenario/GCM (identified by
#'   `group_cols`). The visitor column may be NA for future rows. Historical
#'   rows are identified via `df[[scenario_col]] == historical_label`; every
#'   other value of `scenario_col` is treated as a future climate scenario.
#' @param visitor_col String. Name of the monthly visitor-count column in
#'   `df` (e.g. "TotalVisitors", "LoganVisitors").
#' @param climate_vars Character vector of monthly climate column names in
#'   `df` to use as GAM predictors (e.g. c("tmean_F", "snowpack_in_wbm")).
#'   Aggregated to annual (see `climate_agg_fun`) before fitting. Default
#'   NULL (population-only model: `.visitors ~ pop`).
#' @param date_col,scenario_col,historical_label,group_cols See details above
#'   / the package-level note at the top of this file. Defaults "ym", "CF",
#'   "Historical", c("GCM", "CF").
#' @param growth_pct Numeric. Fractional growth above the historical maximum
#'   used to set the *mid* logistic carrying capacity, i.e.
#'   `K_mid = max(observed annual visitors) * (1 + growth_pct)`. Default 0.5.
#' @param low_pct,high_pct Numeric multipliers applied to `K_mid` to get the
#'   low/high carrying-capacity bounds. Defaults 0.8 / 1.2.
#' @param climate_agg_fun Function used to aggregate each monthly climate
#'   variable to annual (applied uniformly to all `climate_vars`). Default
#'   `mean`.
#' @param visitor_agg_fun Function used to aggregate monthly visitors to
#'   annual. Default `sum`.
#' @param exclude_years Integer vector of calendar years to exclude when
#'   fitting the logistic models and computing monthly proportions. Default
#'   2020.
#' @param exclude_years_gam Integer vector of years to exclude when fitting
#'   the GAM specifically. Default NULL.
#' @param smooth_terms `TRUE`, `FALSE` (default), or `"auto"`. If `TRUE`,
#'   continuous climate predictors are wrapped in `mgcv::s(..., k =
#'   smooth_k)`; if `FALSE` they enter as linear terms. If `"auto"`, both a
#'   linear and a smoothed model are fit (via REML) and whichever has the
#'   lower AIC is used -- letting the historical data decide rather than
#'   requiring a manual choice, at the cost of fitting two models instead of
#'   one. Has no effect if `climate_vars` is empty (nothing to smooth).
#' @param smooth_k Integer. Basis dimension (`k`) for each `s()` term when
#'   smoothing is used. Default 4 (deliberately conservative). Annual
#'   visitation series are typically short (a few dozen years at most), and
#'   mgcv's own default (`k = 10`) can consume more effective degrees of
#'   freedom per climate predictor than a small annual series can support,
#'   especially with more than one `climate_vars` entry -- each additional
#'   unconstrained smooth term increases the risk of overfitting the
#'   historical noise rather than capturing a real climate response. Raise
#'   this only if diagnostics (e.g. `mgcv::gam.check()` on `$gam_model`)
#'   indicate the fit is too constrained for a genuinely nonlinear
#'   relationship, and always check against the sample size actually
#'   available (`nrow(gam_dat)`, visible via `summary($gam_model)`).
#' @param hist_nls_start,proj_nls_start Named lists of `nls()` starting
#'   values (`r`, `t0`). Defaults `t0 = 2010` / `t0 = 2050`; may need
#'   per-park tuning if `nls()` fails to converge.
#' @param proj_period Length-2 numeric vector, future summary period for
#'   `$proj_period_summary`. Default c(2035, 2064).
#' @param obs_cutoff_year Integer. Years before this use observed
#'   visitation; years at/after it use the projection. Default current year.
#' @param future_end_year Integer. Last year to project through. Default
#'   2100.
#' @param output_resolution "monthly" or "annual" -- which resolution to
#'   compute and return as `$result`. Input `df` must be monthly either way.
#'   Default "monthly".
#' @param park_pal Optional named list/vector of colors for the diagnostic
#'   plot (expects e.g. `$pink`); falls back to "tomato".
#' @param summary_years Optional named list of numeric year vectors giving
#'   one or more periods to summarize monthly seasonality for -- e.g. the
#'   default `list(historical = 1980:2012, future = 2035:2064)` compares
#'   historical monthly proportions against those implied by the projection
#'   for a future window. Each period's proportions are computed from real
#'   observed monthly visitation where the period falls in the historical
#'   record, and from the projected monthly envelope (`visitors_mid`) where
#'   it falls in the future -- see `$monthly_proportions_summary` below. Set
#'   to `NULL` to skip this (no extra computation, output unchanged).
#'
#' @return A list: `$result` (table matching `output_resolution`), `$annual`,
#'   `$monthly` (if requested), `$monthly_proportions`,
#'   `$monthly_proportions_summary` (if `summary_years` is supplied -- one
#'   stacked table with a `period` column, same columns as
#'   `$monthly_proportions` but computed separately per named period instead
#'   of across the whole series), `$logistic_models`
#'   (low/mid/high nls objects), `$K_low`/`$K_mid`/`$K_high` (the population
#'   ceilings actually used, for reporting/QA), `$gam_model`,
#'   `$gam_used_smooth` (logical, whether the fitted model ended up using
#'   smooth climate terms -- always matches `smooth_terms` unless it was
#'   `"auto"`), `$gam_concurvity` (from `mgcv::concurvity()`, NULL if there
#'   are no `climate_vars`; see Details), `$gam_diagnostics` (ggarrange
#'   panel), `$proj_period_summary`.
#'
#' @details Untested against live data -- see the note at the top of this
#'   file.
#'
#'   `$gam_concurvity` is worth checking whenever `climate_vars` is
#'   non-empty: concurvity is the GAM analog of multicollinearity -- it
#'   measures how well one smooth term could be approximated by the others.
#'   Because `pop` (the logistic population term) and a climate variable
#'   like mean temperature can both trend upward over the same historical
#'   period, the model may not always be able to cleanly separate "this is
#'   population growth" from "this is a climate effect." High concurvity
#'   (values approaching 1 in the "worst" row) is a sign the split between
#'   those two terms in the fitted coefficients is less trustworthy than the
#'   model's R2 alone would suggest, even if the overall fit looks good.
#'
#' @examples
#' \dontrun{
#' vis_proj <- gam_logistic(
#'   df           = zion_visitation_climate,   # ym, TotalVisitors, tmean_F, GCM, CF
#'   visitor_col  = "TotalVisitors",
#'   climate_vars = c("tmean_F", "snowpack_in_wbm"),
#'   growth_pct   = 0.3
#' )
#' vis_proj$gam_diagnostics
#' vis_proj$result
#' }
gam_logistic <- function(df,
                         visitor_col,
                         climate_vars = NULL,
                         date_col = "ym",
                         scenario_col = "CF",
                         historical_label = "Historical",
                         group_cols = c("GCM", "CF"),
                         growth_pct = 0.5,
                         low_pct = 0.8,
                         high_pct = 1.2,
                         climate_agg_fun = mean,
                         visitor_agg_fun = sum,
                         exclude_years = 2020,
                         exclude_years_gam = NULL,
                         smooth_terms = FALSE,
                         smooth_k = 4,
                         hist_nls_start = list(r = 0.05, t0 = 2010),
                         proj_nls_start = list(r = 0.05, t0 = 2050),
                         proj_period = c(2035, 2064),
                         obs_cutoff_year = as.numeric(format(Sys.Date(), "%Y")),
                         future_end_year = 2100,
                         output_resolution = c("monthly", "annual"),
                         park_pal = NULL,
                         summary_years = list(historical = 1980:2012, future = 2035:2064)) {

  output_resolution <- match.arg(output_resolution)
  if (!(isTRUE(smooth_terms) || isFALSE(smooth_terms) || identical(smooth_terms, "auto"))) {
    stop("gam_logistic(): `smooth_terms` must be TRUE, FALSE, or \"auto\" (got: ",
         paste(smooth_terms, collapse = ", "), ").")
  }
  prep <- .prep_vis_input(df, date_col, visitor_col, climate_vars, scenario_col,
                          group_cols, historical_label, future_end_year)
  climate_vars <- prep$climate_vars
  hist_monthly <- prep$hist_monthly
  fut_monthly  <- prep$fut_monthly

  accent_col <- if (!is.null(park_pal) && !is.null(park_pal$pink)) park_pal$pink else "tomato"
  
  annual_hist_visitors <- .annual_visitors(hist_monthly, visitor_agg_fun)
  annual_hist_climate  <- .annual_climate_agg(hist_monthly, climate_vars, climate_agg_fun, group_vars = "year")
  monthly_proportions  <- .monthly_proportions(hist_monthly, exclude_years)
  
  ## ---- Logistic growth model (low/mid/high carrying capacity) -------------
  logistic_mod <- function(K, r, t0, year) K / (1 + exp(-r * (year - t0)))
  
  K_mid  <- max(annual_hist_visitors$.visitors, na.rm = TRUE) * (1 + growth_pct)
  K_low  <- K_mid * low_pct
  K_high <- K_mid * high_pct
  
  # algorithm = "port" + explicit lower/upper bounds keeps the fit from
  # wandering into a region where exp() overflows and the Gauss-Newton
  # numericDeriv() step blows up ("Missing value or an infinity produced") --
  # the unconstrained default algorithm is fragile for logistic curves unless
  # the starting values happen to already be close to the solution.
  fit_hist_logistic <- function(K, dat) {
    nls(.visitors ~ logistic_mod(K, r, t0, year), data = dat, start = hist_nls_start,
        algorithm = "port",
        lower = c(r = 1e-4, t0 = min(dat$year) - 50),
        upper = c(r = 5,    t0 = max(dat$year) + 100),
        control = nls.control(maxiter = 200, warnOnly = TRUE))
  }
  
  vis_an_fit <- annual_hist_visitors %>%
    dplyr::filter(!year %in% exclude_years) %>%
    tidyr::drop_na(.visitors, year)
  
  if (nrow(vis_an_fit) < 3) {
    stop("gam_logistic(): fewer than 3 valid annual observations available for the ",
         "logistic fit after removing NA and excluded years. Check `df` for years with ",
         "no visitor data (a fully-NA/zero year collapses to NA after annual aggregation).")
  }
  
  logistic_low  <- fit_hist_logistic(K_low,  vis_an_fit)
  logistic_mid  <- fit_hist_logistic(K_mid,  vis_an_fit)
  logistic_high <- fit_hist_logistic(K_high, vis_an_fit)
  
  ## ---- Annual GAM: visitors ~ climate + population -------------------------
  gam_dat <- annual_hist_visitors %>%
    dplyr::left_join(annual_hist_climate, by = "year") %>%
    dplyr::arrange(year) %>%
    dplyr::mutate(pop = logistic_mod(K_mid, coef(logistic_mid)["r"], coef(logistic_mid)["t0"], year))
  
  if (!is.null(exclude_years_gam)) {
    gam_dat <- gam_dat %>% dplyr::filter(!year %in% exclude_years_gam)
  }
  
  # k is capped (smooth_k, default 4) rather than left at mgcv's default
  # (k = 10) -- annual visitation series are typically only a few dozen
  # points, and an unconstrained s() per climate_vars entry can consume more
  # effective degrees of freedom than a short series can support, especially
  # with more than one climate predictor. method = "REML" (rather than
  # mgcv::gam()'s older GCV.Cp default) gives more stable smoothness
  # selection and is what makes the AIC comparison below valid when
  # smooth_terms = "auto" -- both candidates need to be fit the same way for
  # their AICs to be comparable.
  build_gam_formula <- function(use_smooth) {
    terms <- if (use_smooth && length(climate_vars) > 0) {
      paste0("s(", climate_vars, ", k = ", smooth_k, ")")
    } else {
      climate_vars
    }
    reformulate(c(terms, "pop"), response = ".visitors")
  }

  if (identical(smooth_terms, "auto")) {
    if (length(climate_vars) == 0) {
      # Nothing to smooth -- "auto" degenerates to the linear (pop-only) fit.
      use_smooth <- FALSE
    } else {
      fit_linear <- mgcv::gam(build_gam_formula(FALSE), data = gam_dat, method = "REML")
      fit_smooth <- mgcv::gam(build_gam_formula(TRUE),  data = gam_dat, method = "REML")
      use_smooth <- AIC(fit_smooth) < AIC(fit_linear)
    }
  } else {
    use_smooth <- isTRUE(smooth_terms)
  }

  gam_formula <- build_gam_formula(use_smooth)
  gam_model   <- mgcv::gam(gam_formula, data = gam_dat, method = "REML")

  gam_concurvity <- if (length(climate_vars) > 0) {
    tryCatch(mgcv::concurvity(gam_model, full = TRUE), error = function(e) NULL)
  } else {
    NULL
  }

  gam_dat <- gam_dat %>%
    dplyr::mutate(pred_vis = pmax(0, predict(gam_model, newdata = gam_dat)),
                  pred_vis_residual = .visitors - pred_vis)
  
  r2_val <- round(summary(gam_model)$r.sq, 2)
  p_val  <- format.pval(mean(summary(gam_model)$p.pv), digits = 2, eps = 0.001)
  
  gam_diagnostics <- ggpubr::ggarrange(
    ggplot(gam_dat) +
      geom_line(aes(x = year, y = .visitors, color = "Observed"), linewidth = 1) +
      geom_line(aes(x = year, y = pred_vis,  color = "Modeled"),  linewidth = 1) +
      theme_bw() +
      theme(legend.position = c(0.85, 0.25),
            legend.background = element_rect(fill = scales::alpha("white", 0.7), color = NA),
            legend.key = element_blank()) +
      annotate("text",
               label = as.character(as.expression(bquote(R^2 == .(r2_val) ~~~ p < .(p_val)))),
               parse = TRUE, x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
               size = 4, fontface = "italic") +
      labs(x = "Year", y = "Annual Visitors", title = "Annual Model Performance") +
      scale_color_manual("", values = c("Modeled" = accent_col, "Observed" = "black")),
    ggplot(gam_dat, aes(x = year, y = pred_vis_residual)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      geom_line() + geom_point() +
      labs(y = "Residual (# visitors)", x = "Year", title = "Annual Residuals") +
      theme_bw(),
    ggplot(gam_dat, aes(x = .visitors, y = pred_vis)) +
      geom_abline(intercept = 0, slope = 1, linetype = "dotted", linewidth = 1) +
      geom_point(alpha = 0.7, size = 2) +
      theme_bw() +
      labs(x = "Observed Annual Visitors", y = "Modeled Annual Visitors", title = "1:1 Accuracy"),
    ncol = 1, align = "hv"
  )
  
  ## ---- Project across future climate scenarios ------------------------------
  # predict.gam() is row-wise and doesn't depend on grouping, so we can call
  # predict() once on the whole future annual table rather than looping per
  # GCM/CF group.
  group_vars_annual  <- unique(c("year", group_cols))
  annual_fut_climate <- .annual_climate_agg(fut_monthly, climate_vars, climate_agg_fun,
                                            group_vars = group_vars_annual)
  
  annual_proj <- annual_fut_climate %>%
    dplyr::mutate(pop = logistic_mod(K_mid, coef(logistic_mid)["r"], coef(logistic_mid)["t0"], year))
  
  preds <- predict(gam_model, newdata = annual_proj, se.fit = TRUE)
  
  annual_proj <- annual_proj %>%
    dplyr::mutate(vis_pred      = as.numeric(preds$fit),
                  vis_pred_se   = as.numeric(preds$se.fit),
                  vis_pred_low  = pmax(0, vis_pred - vis_pred_se * 1.96),
                  vis_pred_high = vis_pred + vis_pred_se * 1.96) %>%
    dplyr::group_by(year) %>%
    dplyr::mutate(vis_yearly_max  = max(vis_pred_high, na.rm = TRUE),
                  vis_yearly_min  = min(vis_pred_low,  na.rm = TRUE),
                  vis_yearly_mean = mean(c(vis_yearly_max, vis_yearly_min), na.rm = TRUE)) %>%
    dplyr::ungroup()
  
  ## ---- Smooth logistic envelopes over the noisy GAM projections -------------
  # K is passed in explicitly (K_low/K_mid/K_high, computed above from
  # growth_pct/low_pct/high_pct) rather than re-derived from the GAM output
  # (as an earlier version did via K_start <- max(y_var)). Fixing K at fit
  # time is what actually makes it a ceiling: the logistic form
  # K/(1+exp(-r*(year-t0))) is mathematically bounded by K for every year, so
  # the smoothed curve can never exceed the population cap you set. Deriving
  # K from the GAM's own predictions instead would just describe whatever
  # the GAM happened to produce -- growth_pct/low_pct/high_pct would then
  # have no real effect on the final output, only on the "pop" covariate fed
  # into the GAM (which the GAM is free to override via the climate terms).
  fit_proj_logistic <- function(y_var, dat, K) {
    if (max(y_var, na.rm = TRUE) > K) {
      warning("gam_logistic(): the GAM's climate-driven projection exceeds the population ",
              "ceiling (K = ", round(K), ") for at least one year -- the smoothed logistic ",
              "curve will saturate at K rather than following the raw GAM values. Consider ",
              "raising growth_pct/high_pct if this isn't intended, or check the climate ",
              "predictors' influence on the GAM.")
    }
    nls(y_var ~ logistic_mod(K, r, t0, year), data = dat, start = proj_nls_start,
        algorithm = "port",
        lower = c(r = 1e-4, t0 = min(dat$year) - 50),
        upper = c(r = 5,    t0 = max(dat$year) + 100),
        control = nls.control(maxiter = 200, warnOnly = TRUE))
  }
  
  annual_envelope <- annual_proj %>%
    dplyr::distinct(year, vis_yearly_mean, vis_yearly_max, vis_yearly_min)
  
  if (nrow(annual_envelope) == 0) {
    stop("gam_logistic(): no future climate-scenario rows found in `df` (rows where ",
         "df[[scenario_col]] != historical_label). This commonly happens when `df` is ",
         "built by left_join()-ing a climate table ONTO a visitation table -- if the ",
         "visitation table only has historical months, future rows get dropped by the ",
         "join. Build `df` the other way around: start from the full historical+future ",
         "climate time series and left_join() visitation onto IT, so future months are ",
         "kept (with NA in `visitor_col`).")
  }
  
  log_mid  <- fit_proj_logistic(annual_envelope$vis_yearly_mean, annual_envelope, K = K_mid)
  log_high <- fit_proj_logistic(annual_envelope$vis_yearly_max,  annual_envelope, K = K_high)
  log_low  <- fit_proj_logistic(annual_envelope$vis_yearly_min,  annual_envelope, K = K_low)
  
  annual_envelope <- annual_envelope %>%
    dplyr::mutate(visitors_mid  = as.numeric(predict(log_mid,  newdata = .)),
                  visitors_high = as.numeric(predict(log_high, newdata = .)),
                  visitors_low  = as.numeric(predict(log_low,  newdata = .)))
  
  proj_period_summary <- annual_envelope %>%
    dplyr::filter(year >= proj_period[1], year <= proj_period[2]) %>%
    dplyr::summarize(mean_vis = mean(visitors_mid,  na.rm = TRUE),
                     max_vis  = max(visitors_high, na.rm = TRUE),
                     min_vis  = min(visitors_low,  na.rm = TRUE))
  
  annual_result <- .assemble_annual_result(annual_hist_visitors, annual_envelope, obs_cutoff_year)
  
  out <- list(
    annual = annual_result,
    monthly_proportions = monthly_proportions,
    logistic_models = list(low = logistic_low, mid = logistic_mid, high = logistic_high),
    K_low = K_low, K_mid = K_mid, K_high = K_high,
    gam_model = gam_model,
    gam_used_smooth = use_smooth,
    gam_concurvity = gam_concurvity,
    gam_diagnostics = gam_diagnostics,
    proj_period_summary = proj_period_summary
  )
  
  need_monthly <- output_resolution == "monthly" || !is.null(summary_years)

  if (need_monthly) {
    monthly_result <- .disaggregate_monthly(min(prep$df$.ym), annual_hist_visitors, annual_envelope,
                                            monthly_proportions, obs_cutoff_year, future_end_year)
  }

  if (output_resolution == "monthly") {
    out$monthly <- monthly_result
    out$result  <- monthly_result
  } else {
    out$result <- annual_result
  }

  if (!is.null(summary_years)) {
    out$monthly_proportions_summary <- .monthly_proportions_by_period(
      .proportions_source(hist_monthly, monthly_result, exclude_years),
      summary_years
    )
  }
  
  out
}

# ============================================================================
# ---- Method 4: arima_climate -------------------------------------------------
# ============================================================================

#' Project future park visitation via ARIMA (optionally with climate as xreg)
#'
#' Fits an automatically-selected ARIMA model (`forecast::auto.arima()`) on
#' annual historical visitation, optionally using annual climate variables as
#' external regressors (xreg). Requires the {forecast} package.
#'
#' Unlike `gam_logistic()`/`linear_trend()`, `exclude_years` here is used
#' ONLY for the monthly-proportion calculation, not for the ARIMA fit itself
#' -- `auto.arima()`'s automatic order selection is not reliably robust to
#' internal missing values, so introducing an artificial gap for an excluded
#' year is more likely to break the fit than help it. Pre-filter `df` before
#' calling if you need a year fully removed from the fit.
#'
#' @inheritParams linear_trend
#' @param forecast_level Confidence level (as a percentage, per the
#'   {forecast} package convention, e.g. 95 not 0.95) used for the
#'   prediction interval that becomes the low/high bounds. Default 95.
#' @param arima_args Named list of additional arguments passed to
#'   `forecast::auto.arima()` (e.g. `list(seasonal = FALSE)`). Default
#'   `list()`.
#'
#' @return Same shape as `linear_trend()`, with `$arima_model` (the fitted
#'   `Arima` object) in place of `$lm_model`. `$monthly_proportions_summary`
#'   is included whenever `summary_years` is supplied (see `gam_logistic()`).
#'
#' @details Untested against live data -- see the note at the top of this
#'   file. If `climate_vars` is supplied, one `forecast()` call is made per
#'   future climate scenario group (`group_cols`), since each scenario has
#'   its own xreg trajectory; if `climate_vars` is NULL, a single forecast is
#'   made and applied identically to all scenarios (there's nothing to
#'   differentiate them).
#'
#' @examples
#' \dontrun{
#' vis_proj <- arima_climate(
#'   df           = zion_visitation_climate,
#'   visitor_col  = "TotalVisitors",
#'   climate_vars = c("tmean_F")
#' )
#' vis_proj$arima_model
#' vis_proj$result
#' }
arima_climate <- function(df,
                          visitor_col,
                          climate_vars = NULL,
                          date_col = "ym",
                          scenario_col = "CF",
                          historical_label = "Historical",
                          group_cols = c("GCM", "CF"),
                          climate_agg_fun = mean,
                          visitor_agg_fun = sum,
                          exclude_years = 2020,
                          forecast_level = 95,
                          arima_args = list(),
                          obs_cutoff_year = as.numeric(format(Sys.Date(), "%Y")),
                          future_end_year = 2100,
                          output_resolution = c("monthly", "annual"),
                          summary_years = list(historical = 1980:2012, future = 2035:2064)) {
  
  if (!requireNamespace("forecast", quietly = TRUE)) {
    stop("arima_climate() requires the {forecast} package. Install it with ",
         "install.packages('forecast').")
  }
  
  output_resolution <- match.arg(output_resolution)
  prep <- .prep_vis_input(df, date_col, visitor_col, climate_vars, scenario_col,
                          group_cols, historical_label, future_end_year)
  climate_vars <- prep$climate_vars
  hist_monthly <- prep$hist_monthly
  fut_monthly  <- prep$fut_monthly
  has_climate  <- length(climate_vars) > 0
  
  annual_hist_visitors <- .annual_visitors(hist_monthly, visitor_agg_fun)
  monthly_proportions  <- .monthly_proportions(hist_monthly, exclude_years)
  
  ## ---- Build a gap-free annual ts() over the full historical year range -----
  year_range <- min(annual_hist_visitors$year):max(annual_hist_visitors$year)
  
  vis_ts_df <- data.frame(year = year_range) %>%
    dplyr::left_join(annual_hist_visitors, by = "year")
  
  vis_ts <- stats::ts(vis_ts_df$.visitors, start = min(year_range), frequency = 1)
  
  if (has_climate) {
    annual_hist_climate <- .annual_climate_agg(hist_monthly, climate_vars, climate_agg_fun, group_vars = "year")
    xreg_hist <- data.frame(year = year_range) %>%
      dplyr::left_join(annual_hist_climate, by = "year") %>%
      dplyr::select(dplyr::all_of(climate_vars)) %>%
      as.matrix()
    
    arima_fit <- do.call(forecast::auto.arima, c(list(y = vis_ts, xreg = xreg_hist), arima_args))
  } else {
    arima_fit <- do.call(forecast::auto.arima, c(list(y = vis_ts), arima_args))
  }
  
  ## ---- Project forward ---------------------------------------------------------
  if (has_climate) {
    annual_fut_climate <- .annual_climate_agg(fut_monthly, climate_vars, climate_agg_fun,
                                              group_vars = unique(c("year", group_cols))) %>%
      dplyr::filter(year <= future_end_year)
    
    annual_proj <- annual_fut_climate %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
      dplyr::group_split() %>%
      purrr::map_dfr(function(this_scenario) {
        this_scenario <- this_scenario %>% dplyr::arrange(year)
        xreg_fut <- as.matrix(this_scenario[, climate_vars, drop = FALSE])
        fc <- forecast::forecast(arima_fit, xreg = xreg_fut, level = forecast_level)
        this_scenario %>%
          dplyr::mutate(vis_pred      = as.numeric(fc$mean),
                        vis_pred_low  = pmax(0, as.numeric(fc$lower[, 1])),
                        vis_pred_high = as.numeric(fc$upper[, 1]))
      })
  } else {
    future_years <- (max(year_range) + 1):future_end_year
    fc <- forecast::forecast(arima_fit, h = length(future_years), level = forecast_level)
    annual_proj <- data.frame(year = future_years,
                              vis_pred      = as.numeric(fc$mean),
                              vis_pred_low  = pmax(0, as.numeric(fc$lower[, 1])),
                              vis_pred_high = as.numeric(fc$upper[, 1]))
  }
  
  annual_envelope <- annual_proj %>%
    dplyr::group_by(year) %>%
    dplyr::summarize(visitors_mid  = mean(vis_pred,      na.rm = TRUE),
                     visitors_high = max(vis_pred_high, na.rm = TRUE),
                     visitors_low  = min(vis_pred_low,  na.rm = TRUE),
                     .groups = "drop")
  
  annual_result <- .assemble_annual_result(annual_hist_visitors, annual_envelope, obs_cutoff_year)
  
  out <- list(
    annual = annual_result,
    monthly_proportions = monthly_proportions,
    arima_model = arima_fit
  )
  
  need_monthly <- output_resolution == "monthly" || !is.null(summary_years)

  if (need_monthly) {
    monthly_result <- .disaggregate_monthly(min(prep$df$.ym), annual_hist_visitors, annual_envelope,
                                            monthly_proportions, obs_cutoff_year, future_end_year)
  }

  if (output_resolution == "monthly") {
    out$monthly <- monthly_result
    out$result  <- monthly_result
  } else {
    out$result <- annual_result
  }

  if (!is.null(summary_years)) {
    out$monthly_proportions_summary <- .monthly_proportions_by_period(
      .proportions_source(hist_monthly, monthly_result, exclude_years),
      summary_years
    )
  }
  
  out
}

# ============================================================================
# ---- Method 5: fisichelli ---------------------------------------------------
# ============================================================================

#' Project future park visitation via the Fisichelli et al. (2015) method
#'
#' Follows Fisichelli, Schuurman, Monahan & Ziesler (2015, PLoS ONE), "Protected
#' Area Tourism in a Changing Climate: Will Visitation at US National Parks Warm
#' Up or Overheat?": fits a monthly temperature-visitation curve using
#' long-term calendar-month climatology (one point per month -- Jan average,
#' Feb average, etc. -- rather than every individual year-month), selects
#' among null/linear/quadratic/cubic candidate models by BIC, and only uses
#' the temperature relationship for projection if it passes a screening rule
#' (temperature actually in the selected model, and adjusted R2 >= a
#' threshold). Parks/periods that fail the screen carry forward their
#' historical monthly climatology unchanged (no inferred climate effect, but
#' the seasonal pattern is preserved).
#'
#' This is fundamentally different from `gam_logistic()` /
#' `linear_trend()` / `arima_climate()` in two ways, both intentional and
#' worth understanding before comparing results:
#'
#' 1. There is no year-over-year growth term (no logistic cap, no `year`
#'    predictor, no ARIMA trend) -- the model relates monthly temperature to
#'    monthly visitation, full stop. Any long-term growth in your data should
#'    be handled separately if you need it (e.g. blend this method's climate
#'    signal with one of the other methods' population trend).
#' 2. `visitors_low`/`visitors_mid`/`visitors_high` are NOT statistical
#'    confidence/uncertainty bounds. Per Fisichelli et al., they're three
#'    different GROWTH CEILINGS applied to the same underlying
#'    temperature-driven prediction, based on the park's historically
#'    busiest month: a conservative cap (`cap_low_pct`, default 5% above
#'    that historical busiest-month average), a middle cap (`cap_mid_pct`,
#'    default the midpoint between `cap_low_pct` and `cap_high_pct`), and a
#'    generous cap (`cap_high_pct`, default 25% above). The idea is that
#'    visitation in already-at-capacity months can't grow indefinitely just
#'    because temperature is projected to keep changing -- the ceilings
#'    represent alternate capacity-constraint assumptions, not statistical
#'    uncertainty. Because all three are the same raw prediction passed
#'    through three caps of increasing size, `visitors_low <= visitors_mid
#'    <= visitors_high` is guaranteed for every year -- unlike an earlier
#'    version of this function, which left `visitors_mid` uncapped and could
#'    show it exceeding `visitors_high`.
#'
#' @param df A single MONTHLY time series data frame -- see `gam_logistic()`
#'   for the shared input contract (historical rows flagged via
#'   `scenario_col`/`historical_label`; future rows one per month per
#'   scenario, identified by `group_cols`).
#' @param visitor_col String. Name of the monthly visitor-count column.
#' @param temp_var String. Name of the single monthly temperature column to
#'   use as the predictor (Fisichelli et al. is a univariate
#'   temperature-only method -- unlike the other three methods here, this
#'   isn't a vector of arbitrary climate variables).
#' @param date_col,scenario_col,historical_label,group_cols See
#'   `gam_logistic()`. Same defaults.
#' @param degree_max Integer. Highest polynomial degree to try (in addition
#'   to the null/intercept-only model). Default 3 (null, linear, quadratic,
#'   cubic), matching the original method.
#' @param r2_threshold Adjusted-R2 screening threshold below which the
#'   temperature model is rejected in favor of the historical monthly
#'   climatology. Default 0.5, matching Fisichelli et al.
#' @param cap_low_pct,cap_high_pct Fractional growth above the historical
#'   busiest-month average used for the conservative/generous visitation
#'   ceilings (see Details). Defaults 0.05 / 0.25, matching the original.
#' @param cap_mid_pct Fractional growth above the historical busiest-month
#'   average used for the middle visitation ceiling that produces
#'   `visitors_mid` (see Details). Default NULL, which uses the midpoint of
#'   `cap_low_pct` and `cap_high_pct`.
#' @param exclude_years Optional integer vector of years to exclude from the
#'   monthly climatology and model fit. Default NULL (the original method
#'   does not exclude any particular year).
#' @param visitor_agg_fun Function used to aggregate monthly visitors to
#'   annual for the `$annual`/observed portion of the output. Default `sum`.
#' @param obs_cutoff_year,future_end_year,output_resolution See
#'   `gam_logistic()`. Same defaults.
#' @param summary_years See `gam_logistic()`. Worth noting for this method
#'   specifically: unlike the other four methods (which disaggregate to
#'   monthly using ONE fixed historical fraction per calendar month, so
#'   their future-period proportions can only differ from the historical
#'   ones if the requested periods pull from different years), `fisichelli()`
#'   predicts each future month directly from its own projected
#'   temperature -- so a real shift in seasonality (e.g. earlier/later
#'   shoulder-season visitation under warming) can actually show up as a
#'   difference between the `historical` and `future` rows here.
#'
#' @return Same overall shape as the other methods (`$result`, `$annual`,
#'   `$monthly` if requested), plus `$temp_model` (the selected lm object),
#'   `$model_name` ("null"/"linear"/"quadratic"/"cubic"), `$temp_sensitive`
#'   (logical, whether the screening rule passed), `$adj_r2`,
#'   `$monthly_climatology` (the 12-point, or fewer, training table), and
#'   `$monthly_proportions_summary` (if `summary_years` is supplied).
#'
#' @details Untested against live data -- see the note at the top of this
#'   file. The polynomial models are fit on only 12 (or fewer) climatology
#'   points, so cubic/quadratic fits have very few residual degrees of
#'   freedom; candidate degrees are skipped automatically if there isn't
#'   enough data to support them (`degree + 2` climatology points minimum),
#'   but a `summary(fit$temp_model)` sanity check is worth doing regardless.
#'
#' @examples
#' \dontrun{
#' vis_proj <- fisichelli(
#'   df          = zion_visitation_climate,   # ym, TotalVisitors, tmean_F, GCM, CF
#'   visitor_col = "TotalVisitors",
#'   temp_var    = "tmean_F"
#' )
#' vis_proj$temp_sensitive
#' vis_proj$result
#' }
fisichelli <- function(df,
                       visitor_col,
                       temp_var,
                       date_col = "ym",
                       scenario_col = "CF",
                       historical_label = "Historical",
                       group_cols = c("GCM", "CF"),
                       degree_max = 3,
                       r2_threshold = 0.5,
                       cap_low_pct = 0.05,
                       cap_mid_pct = NULL,
                       cap_high_pct = 0.25,
                       exclude_years = NULL,
                       visitor_agg_fun = sum,
                       obs_cutoff_year = as.numeric(format(Sys.Date(), "%Y")),
                       future_end_year = 2100,
                       output_resolution = c("monthly", "annual"),
                       summary_years = list(historical = 1980:2012, future = 2035:2064)) {
  
  if (length(temp_var) != 1) {
    stop("fisichelli(): `temp_var` must name a single temperature column, not a vector.")
  }
  
  # Runs `expr`, and on error re-throws with a stage label prefixed, so a
  # failure says WHICH part of the method broke instead of a bare message
  # bubbling up from deep inside lm()/predict()/dplyr.
  .stage <- function(label, expr) {
    tryCatch(expr, error = function(e) {
      stop("fisichelli(): failed during '", label, "' -- ", conditionMessage(e), call. = FALSE)
    })
  }
  
  output_resolution <- match.arg(output_resolution)
  prep <- .stage("input prep", .prep_vis_input(df, date_col, visitor_col, temp_var, scenario_col,
                                               group_cols, historical_label, future_end_year))
  hist_monthly <- prep$hist_monthly
  fut_monthly  <- prep$fut_monthly
  
  if (nrow(fut_monthly) == 0) {
    stop("fisichelli(): no future climate-scenario rows found in `df` (rows where ",
         "df[[scenario_col]] != historical_label). See the join-direction note from ",
         "gam_logistic() troubleshooting -- `df` likely needs to be built starting from ",
         "the full historical+future climate series, not from the visitation table.")
  }
  
  ## ---- Long-term monthly climatology (one point per calendar month) ---------
  hist_for_climatology <- if (!is.null(exclude_years)) {
    hist_monthly %>% dplyr::filter(!year %in% exclude_years)
  } else {
    hist_monthly
  }
  
  monthly_climatology <- .stage("climatology construction", {
    mc <- hist_for_climatology %>%
      dplyr::group_by(month) %>%
      dplyr::summarize(mean_temp     = mean(.data[[temp_var]], na.rm = TRUE),
                       mean_visitors = mean(.visitors, na.rm = TRUE),
                       n_years       = dplyr::n_distinct(year),
                       .groups = "drop")
    
    # Fill any missing calendar months with the across-month mean, so the
    # climatology (and downstream fallback projection) always has 12 rows.
    data.frame(month = 1:12) %>%
      dplyr::left_join(mc, by = "month") %>%
      dplyr::mutate(mean_visitors = ifelse(is.na(mean_visitors),
                                           mean(mean_visitors, na.rm = TRUE),
                                           mean_visitors))
  })
  
  if (all(is.na(monthly_climatology$mean_temp))) {
    stop("fisichelli(): `temp_var` ('", temp_var, "') is NA for every historical month -- ",
         "check that column in `df` before calling.")
  }
  if (all(is.na(monthly_climatology$mean_visitors))) {
    stop("fisichelli(): `visitor_col` ('", visitor_col, "') is NA for every historical month -- ",
         "check that column in `df` before calling.")
  }
  
  ## ---- Fit null/linear/quadratic/cubic candidates, select by BIC ------------
  fit_dat <- monthly_climatology %>% tidyr::drop_na(mean_temp, mean_visitors)
  
  fit_and_select_temp_model <- function(dat, degree_max) {
    n_obs <- nrow(dat)
    degree_names <- c("linear", "quadratic", "cubic")
    
    candidates <- list(null = lm(mean_visitors ~ 1, data = dat))
    for (d in seq_len(degree_max)) {
      if (n_obs >= d + 2) {
        nm <- if (d <= length(degree_names)) degree_names[d] else paste0("degree", d)
        # Bake the literal degree into the formula text (via reformulate())
        # rather than referencing the loop variable `d` inside poly(). If `d`
        # is left as a variable, poly() stores the SYMBOL in the fitted
        # formula, not its value at fit time -- by the time predict() is
        # called later, `d` has already advanced to its final loop value,
        # so e.g. a quadratic model ends up trying to reconstruct a cubic
        # basis ("was fitted with type nmatrix.2 but type nmatrix.3 was
        # supplied"). Literal text in the formula sidesteps that entirely.
        poly_formula <- reformulate(paste0("poly(mean_temp, ", d, ")"), response = "mean_visitors")
        candidates[[nm]] <- tryCatch(lm(poly_formula, data = dat),
                                     error = function(e) NULL)
      }
    }
    candidates <- candidates[!vapply(candidates, is.null, logical(1))]
    bics <- vapply(candidates, BIC, numeric(1))
    best_name  <- names(which.min(bics))
    best_model <- candidates[[best_name]]
    adj_r2 <- if (best_name == "null") NA_real_ else summary(best_model)$adj.r.squared
    
    list(model_name = best_name, model = best_model,
         bic = unname(bics[best_name]), adj_r2 = adj_r2)
  }
  
  fit_result <- .stage("temperature model fit/selection", fit_and_select_temp_model(fit_dat, degree_max))
  
  temp_in_model  <- fit_result$model_name != "null"
  temp_sensitive <- temp_in_model && !is.na(fit_result$adj_r2) && fit_result$adj_r2 >= r2_threshold
  
  ## ---- Visitation growth ceilings, from the historically busiest month -----
  # arrange()+slice(1) instead of filter(x == max(x)) -- avoids any float
  # equality edge case and always returns exactly one row.
  busiest <- .stage("busiest-month cap calculation", {
    monthly_climatology %>%
      dplyr::arrange(dplyr::desc(mean_visitors)) %>%
      dplyr::slice(1)
  })
  
  if (is.null(cap_mid_pct)) cap_mid_pct <- mean(c(cap_low_pct, cap_high_pct))

  cap_upper <- busiest$mean_visitors * (1 + cap_high_pct)
  cap_mid   <- busiest$mean_visitors * (1 + cap_mid_pct)
  cap_lower <- busiest$mean_visitors * (1 + cap_low_pct)

  ## ---- Project forward: apply the model (or climatology fallback) ----------
  fut_monthly <- .stage("future projection", {
    if (temp_sensitive) {
      newdata_temp <- data.frame(mean_temp = fut_monthly[[temp_var]])
      fut_monthly %>%
        dplyr::mutate(pred_raw = pmax(0, predict(fit_result$model, newdata = newdata_temp)))
    } else {
      fut_monthly %>%
        dplyr::left_join(monthly_climatology %>% dplyr::select(month, mean_visitors), by = "month") %>%
        dplyr::mutate(pred_raw = pmax(0, mean_visitors)) %>%
        dplyr::select(-mean_visitors)
    }
  })

  # All three of pred_low/pred_mid/pred_high are the SAME raw prediction run
  # through three caps of increasing size (cap_lower <= cap_mid <= cap_upper),
  # and pmin() is monotonic in the cap -- so pred_low <= pred_mid <= pred_high
  # holds pointwise for every row. That ordering survives sum() (monthly ->
  # annual) and, since every per-scenario value is bounded between the
  # cross-scenario min(...low) and max(...high), it survives mean()-vs-min/max
  # across scenarios too. This guarantees visitors_low <= visitors_mid <=
  # visitors_high in the final output -- unlike an earlier version, which
  # used the uncapped pred_raw for the mid column and could show it exceeding
  # the capped high column.
  fut_monthly <- fut_monthly %>%
    dplyr::mutate(pred_high = pmin(pred_raw, cap_upper),
                  pred_mid  = pmin(pred_raw, cap_mid),
                  pred_low  = pmin(pred_raw, cap_lower))

  ## ---- Monthly envelope across future climate scenarios ---------------------
  monthly_envelope <- .stage("monthly envelope assembly", {
    fut_monthly %>%
      dplyr::filter(year >= obs_cutoff_year) %>%
      dplyr::group_by(year, month) %>%
      dplyr::summarize(visitors_mid  = mean(pred_mid,  na.rm = TRUE),
                       visitors_high = max(pred_high, na.rm = TRUE),
                       visitors_low  = min(pred_low,  na.rm = TRUE),
                       .groups = "drop") %>%
      dplyr::mutate(ym = lubridate::ym(paste0(year, "-", month)),
                    source = "projected") %>%
      dplyr::select(ym, visitors_low, visitors_mid, visitors_high, source)
  })
  
  monthly_hist <- hist_monthly %>%
    dplyr::filter(year < obs_cutoff_year) %>%
    dplyr::transmute(ym = .ym, visitors_low = .visitors, visitors_mid = .visitors,
                     visitors_high = .visitors, source = "observed")
  
  monthly_result <- dplyr::bind_rows(monthly_hist, monthly_envelope) %>%
    dplyr::arrange(ym)
  
  ## ---- Annual: sum monthly within each scenario trajectory first, then -----
  ## ---- take the envelope across scenarios (matches the original's         --
  ## ---- per-park/gcm/ssp/run annual sum, then cross-scenario summary).      --
  annual_by_scenario <- .stage("annual aggregation by scenario", {
    fut_monthly %>%
      dplyr::filter(year >= obs_cutoff_year) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(group_cols)), year) %>%
      dplyr::summarize(annual_mid  = sum(pred_mid,  na.rm = TRUE),
                       annual_high = sum(pred_high, na.rm = TRUE),
                       annual_low  = sum(pred_low,  na.rm = TRUE),
                       n_months    = dplyr::n(),
                       .groups = "drop") %>%
      dplyr::filter(n_months == 12)  # drop partial years (e.g. a truncated final year)
  })

  if (nrow(annual_by_scenario) == 0) {
    stop("fisichelli(): every future year had fewer than 12 months of data per scenario ",
         "after aggregation, so none survived the complete-year filter -- check `df` for ",
         "gaps in the future climate rows (e.g. a truncated final year, or group_cols not ",
         "uniquely identifying each scenario trajectory).")
  }

  annual_envelope <- annual_by_scenario %>%
    dplyr::group_by(year) %>%
    dplyr::summarize(visitors_mid  = mean(annual_mid,  na.rm = TRUE),
                     visitors_high = max(annual_high, na.rm = TRUE),
                     visitors_low  = min(annual_low,  na.rm = TRUE),
                     .groups = "drop")
  
  annual_hist_visitors <- .annual_visitors(hist_monthly, visitor_agg_fun)
  annual_result <- .assemble_annual_result(annual_hist_visitors, annual_envelope, obs_cutoff_year)
  
  out <- list(
    annual = annual_result,
    monthly_climatology = monthly_climatology,
    temp_model = fit_result$model,
    model_name = fit_result$model_name,
    temp_sensitive = temp_sensitive,
    adj_r2 = fit_result$adj_r2,
    visit_cap_upper = cap_upper,
    visit_cap_mid   = cap_mid,
    visit_cap_lower = cap_lower
  )
  
  if (output_resolution == "monthly") {
    out$monthly <- monthly_result
    out$result  <- monthly_result
  } else {
    out$result <- annual_result
  }

  if (!is.null(summary_years)) {
    out$monthly_proportions_summary <- .monthly_proportions_by_period(
      .proportions_source(hist_monthly, monthly_result, exclude_years),
      summary_years
    )
  }
  
  out
}

# ============================================================================
# ---- Method registry + dispatcher ------------------------------------------
# ============================================================================

#' @keywords internal
.visitation_methods <- list(
  linear_trend       = linear_trend,
  annual_growth_rate = annual_growth_rate,
  gam_logistic       = gam_logistic,
  arima_climate      = arima_climate,
  fisichelli         = fisichelli
  # Add future methods here, e.g.:
  # prophet_climate = prophet_climate
)

#' Project future park visitation (dispatcher)
#'
#' Thin wrapper that calls a named visitation-projection method from
#' `.visitation_methods`: `"linear_trend"`, `"annual_growth_rate"`,
#' `"gam_logistic"`, `"arima_climate"`, or `"fisichelli"` (see each
#' function's own documentation for its specific arguments -- note
#' `fisichelli()` takes `temp_var` (a single column) instead of
#' `climate_vars`, and `annual_growth_rate()` takes `growth_rate_mid` and
#' has no climate argument at all). Add new methods by writing a function
#' with the same `(df, ...)` input contract and `$result`/`$annual`/
#' `$monthly` output shape, then registering it in `.visitation_methods`.
#'
#' @param df Single monthly time series data frame -- see the chosen
#'   method's documentation for the expected columns.
#' @param method String naming the method to use. Default "gam_logistic".
#' @param ... Passed through to the chosen method's function.
#'
#' @examples
#' \dontrun{
#' project_visitation(df, visitor_col = "TotalVisitors",
#'                     climate_vars = c("tmean_F"), method = "arima_climate")
#' }
project_visitation <- function(df, ..., method = "gam_logistic") {
  if (!method %in% names(.visitation_methods)) {
    stop("Unknown visitation projection method: '", method, "'. Available: ",
         paste(names(.visitation_methods), collapse = ", "))
  }
  .visitation_methods[[method]](df, ...)
}

