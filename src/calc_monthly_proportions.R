#' Calculate monthly proportions / climatology for a time series
#'
#' @description
#' Collapses a (daily, sub-monthly, or already-monthly) time series to a
#' monthly total/mean/max, then summarizes across years to a calendar-month
#' climatology -- mean/min/max/sd plus each month's fraction of the annual
#' total -- optionally split by one or more grouping columns already in the
#' data (e.g. `CF`, `GCM`, `system`, a pre- vs. post-event `period`) and/or
#' by one or more named year windows (e.g. a historical vs. mid-century
#' comparison window).
#'
#' @details
#' This generalizes the "monthly proportions" pattern used repeatedly in
#' this workflow (e.g. `monthly_visitor_proportions`, `monthly_leak_summary`,
#' `supply_monthly`, and `.monthly_proportions_by_period()` in
#' `project_visitation.R`) into one reusable function, so the same
#' month-of-year climatology + fraction-of-annual-total logic doesn't need
#' to be hand-rolled again for every new variable.
#'
#' Two steps happen:
#' 1. **Collapse to monthly.** Rows are grouped by year-month (and any
#'    `group_vars`) and reduced to one value per month via
#'    `monthly_agg_fun` ("sum" for additive fluxes like runoff or water-use
#'    gallons, "mean" for intensities like temperature, "max" for a state
#'    variable like snowpack). If `df` is already one row per year-month per
#'    group, this step is a no-op.
#' 2. **Summarize across years.** Those monthly values are summarized by
#'    calendar month (1-12), separately within each `group_vars`
#'    combination and each named window in `year_windows` (if supplied),
#'    giving `mean_val`/`min_val`/`max_val`/`sd_val`/`cv_val`/`n`, plus
#'    `fraction` = that month's share of the annual total within its
#'    group/window. `tidyr::complete()` guarantees 12 rows per group/window
#'    even where some months have no data in that window.
#'
#' @param df A data frame/tibble with at least a date column and a numeric
#'   value column.
#' @param date_col Character; name of the date column in `df`. Default
#'   `"date"`.
#' @param value_col Character; name of the numeric column to summarize.
#' @param group_vars Character vector of column name(s) in `df` to group by
#'   in addition to month (e.g. `c("CF", "system")`, or a pre-existing
#'   `"period"` column). Default `NULL` (no additional grouping).
#' @param monthly_agg_fun One of `"sum"`, `"mean"`, `"max"`; how to collapse
#'   rows within a year-month (and `group_vars`) to a single monthly value
#'   before summarizing across years. No default -- pick deliberately, since
#'   the wrong choice silently changes units/meaning (e.g. averaging daily
#'   flow when a monthly total was intended).
#' @param year_windows Optional named list of length-2 numeric vectors, e.g.
#'   `list(Historical = c(1980, 2012), Future = c(2035, 2065))`. Each window
#'   is inclusive of its start year and exclusive of its end year
#'   (`year >= start & year < end`), matching the `hist_yrs`/`fut_yrs`
#'   convention used elsewhere in this workflow. When supplied, the
#'   climatology is computed separately within each window and results are
#'   stacked with a `period` column (the list's names). When `NULL`
#'   (default), the climatology uses every year present in `df`.
#' @param min_n Minimum number of monthly observations required for a
#'   month/group/period cell to be summarized rather than returned as NA.
#'   Default `1` (no filtering).
#'
#' @return A tibble with `period` (if `year_windows` supplied), `group_vars`
#'   (if any), `month` (1-12), `month_date` (a dummy `2000-MM-01` Date for
#'   x-axis plotting), `mean_val`, `min_val`, `max_val`, `sd_val`, `cv_val`,
#'   `n`, and `fraction`.
#'
#' @examples
#' \dontrun{
#' # Historical vs. future monthly runoff climatology, per CF -- equivalent
#' # to the supply_monthly block in the D:S rewrite:
#' calc_monthly_proportions(
#'   df = wbm, date_col = "date", value_col = "runoff_gal_edith_wbm",
#'   group_vars = "CF", monthly_agg_fun = "sum",
#'   year_windows = list(Historical = c(hist_yrs[1], hist_yrs[2]),
#'                        Future     = c(fut_yrs[1],  fut_yrs[2]))
#' )
#'
#' # Monthly visitation climatology, full record, no extra grouping --
#' # equivalent to the original monthly_visitor_proportions block:
#' calc_monthly_proportions(
#'   df = visitors, date_col = "ym", value_col = "TotalVisitors",
#'   monthly_agg_fun = "mean"
#' )
#'
#' # Pre- vs. post-leak monthly water use -- equivalent to
#' # monthly_leak_summary, using a pre-existing "period" column instead of
#' # year_windows:
#' wateruse_monthly %>%
#'   dplyr::mutate(period = dplyr::if_else(ym < as.Date("2019-06-01"),
#'                                        "Pre-Leak", "Post-Leak")) %>%
#'   calc_monthly_proportions(date_col = "ym", value_col = "total_use",
#'                            group_vars = "period", monthly_agg_fun = "sum")
#' }
#' @export
calc_monthly_proportions <- function(df,
                                      date_col = "date",
                                      value_col,
                                      group_vars = NULL,
                                      monthly_agg_fun,
                                      year_windows = NULL,
                                      min_n = 1) {

  monthly_agg_fun <- match.arg(monthly_agg_fun, c("sum", "mean", "max"))
  agg_fn <- switch(monthly_agg_fun,
                   sum  = function(x) sum(x, na.rm = TRUE),
                   mean = function(x) mean(x, na.rm = TRUE),
                   max  = function(x) max(x, na.rm = TRUE))

  if (is.null(group_vars)) group_vars <- character(0)

  if (!date_col %in% names(df)) {
    stop("`date_col` ('", date_col, "') not found in `df`.")
  }
  if (!value_col %in% names(df)) {
    stop("`value_col` ('", value_col, "') not found in `df`.")
  }
  missing_groups <- setdiff(group_vars, names(df))
  if (length(missing_groups) > 0) {
    stop("`group_vars` not found in `df`: ", paste(missing_groups, collapse = ", "))
  }

  group_syms <- rlang::syms(group_vars)

  dat <- df %>%
    dplyr::rename(.date = dplyr::all_of(date_col),
                  .value = dplyr::all_of(value_col)) %>%
    dplyr::mutate(.ym = lubridate::floor_date(as.Date(.date), "month"))

  # --- Step 1: collapse to one value per year-month (x group_vars) ----------
  monthly_dat <- dat %>%
    dplyr::group_by(.ym, !!!group_syms) %>%
    dplyr::summarize(.value = agg_fn(.value), .groups = "drop") %>%
    dplyr::mutate(year = lubridate::year(.ym),
                  month = lubridate::month(.ym))

  # --- Step 2: summarize across years, by month (x group_vars x period) -----
  summarize_period <- function(dat_sub, period_name = NULL) {

    out <- dat_sub %>%
      dplyr::group_by(!!!group_syms, month) %>%
      dplyr::summarize(
        mean_val = mean(.value, na.rm = TRUE),
        min_val  = min(.value,  na.rm = TRUE),
        max_val  = max(.value,  na.rm = TRUE),
        sd_val   = stats::sd(.value, na.rm = TRUE),
        n        = sum(!is.na(.value)),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        cv_val = sd_val / mean_val,
        dplyr::across(c(mean_val, min_val, max_val, sd_val, cv_val),
                      ~ ifelse(n < min_n, NA_real_, .x))
      ) %>%
      tidyr::complete(month = 1:12, !!!group_syms) %>%
      dplyr::mutate(month_date = lubridate::as_date(paste0("2000-", month, "-01")))

    if (!is.null(period_name)) {
      out <- out %>%
        dplyr::mutate(period = period_name) %>%
        dplyr::relocate(period)
    }

    # fraction = each month's share of its group's (and period's) annual mean
    out %>%
      dplyr::group_by(!!!group_syms, dplyr::across(dplyr::any_of("period"))) %>%
      dplyr::mutate(fraction = mean_val / sum(mean_val, na.rm = TRUE)) %>%
      dplyr::ungroup()
  }

  if (is.null(year_windows)) {

    result <- summarize_period(monthly_dat)

  } else {

    if (!is.list(year_windows) || is.null(names(year_windows)) ||
        any(names(year_windows) == "")) {
      stop("`year_windows` must be a named list of length-2 numeric vectors, ",
           "e.g. list(Historical = c(1980, 2012), Future = c(2035, 2065)).")
    }

    result <- purrr::imap_dfr(year_windows, function(w, nm) {
      if (length(w) != 2) {
        stop("year_windows[['", nm, "']] must have exactly 2 values (start, end).")
      }
      sub <- monthly_dat %>% dplyr::filter(year >= w[1], year < w[2])
      if (nrow(sub) == 0) {
        warning("year_windows['", nm, "'] (", w[1], "-", w[2],
                ") matched no rows in `df` -- its output will be all NA.")
      }
      summarize_period(sub, period_name = nm)
    })
  }

  result %>%
    dplyr::relocate(dplyr::any_of("period"), dplyr::all_of(group_vars),
                    month, month_date)
}
