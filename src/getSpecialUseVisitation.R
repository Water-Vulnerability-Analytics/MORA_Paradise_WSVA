#' Retrieve Special Use Visitation Data for a National Park
#'
#' Scrapes monthly special use visitation data from the NPS IRMA Statistics
#' portal for a given park unit code, optionally saving results to CSV. This function was
#' made by ChatGPT Codex agentic AI.
#'
#' @param park Character. Four-letter NPS park unit code (e.g., \code{"ZION"}).
#' @param output_dir Character. Directory path for saving CSV output.
#'   Defaults to \code{"data/all/special_use_visitation/"}.
#' @param save_csv Logical. Whether to save results as a CSV file. Default \code{TRUE}.
#' @param retries Integer. Number of retry attempts per HTTP request. Default \code{6}.
#' @param quiet Logical. Suppress progress messages if \code{TRUE}. Default \code{FALSE}.
#'
#' @return A \code{tibble} with one row per month-metric combination, containing
#'   columns for park name, park code, month labels, metric name, this-month count,
#'   year-to-date count, and source URL. Attributes include \code{park},
#'   \code{failed_months}, and (if saved) \code{output_csv}.
#'   
getSpecialUseVisitation <- function(park,
                                    output_dir = "data/all/special_use_visitation/",
                                    save_csv = TRUE,
                                    retries = 6,
                                    quiet = FALSE) {
  
  park <- toupper(trimws(park))
  if (!nzchar(park)) stop("park is required.")
  
  landing_url <- sprintf(
    "https://irma.nps.gov/Stats/SSRSReports/Park%%20Specific%%20Reports/Monthly%%20Public%%20Use?Park=%s",
    utils::URLencode(park, reserved = TRUE)
  )
  
  cookie_file <- tempfile(fileext = ".cookies")
  h <- curl::new_handle()
  curl::handle_setheaders(h, "User-Agent" = "Mozilla/5.0-R-script-curl")
  curl::handle_setopt(h, followlocation = 1L, cookiefile = cookie_file, cookiejar = cookie_file)
  
  fetch_text <- function(url) {
    for (i in seq_len(retries)) {
      res <- try(curl::curl_fetch_memory(url, handle = h), silent = TRUE)
      if (!inherits(res, "try-error") && identical(res$status_code, 200L)) {
        txt <- rawToChar(res$content)
        bad <- grepl("Service unavailable|MvcReportViewer Error", txt, ignore.case = TRUE)
        if (!bad) return(txt)
      }
      Sys.sleep(min(2 * i, 10))
    }
    stop("Failed to fetch text: ", url)
  }
  
  fetch_raw <- function(url) {
    for (i in seq_len(retries)) {
      res <- try(curl::curl_fetch_memory(url, handle = h), silent = TRUE)
      if (!inherits(res, "try-error") && identical(res$status_code, 200L) && length(res$content) > 0) {
        return(res$content)
      }
      Sys.sleep(min(2 * i, 10))
    }
    stop("Failed to fetch raw: ", url)
  }
  
  set_rpt_date <- function(url, month_label) {
    enc <- utils::URLencode(month_label, reserved = TRUE)
    if (grepl("([?&])RptDate=", url)) {
      sub("([?&])RptDate=[^&]*", paste0("\\1RptDate=", enc), url)
    } else {
      paste0(url, "&RptDate=", enc)
    }
  }
  
  extract_export_xml_url <- function(viewer_html) {
    rel <- stringr::str_match(viewer_html, "\"ExportUrlBase\":\"([^\"]+Format=)\"")[, 2]
    if (is.na(rel)) stop("Could not find ExportUrlBase.")
    rel <- stringr::str_replace_all(rel, "\\\\u0026", "&")
    rel <- stringr::str_replace_all(rel, "\\\\/", "/")
    paste0("https://irma.nps.gov", rel, "XML")
  }
  
  parse_special_use_xml <- function(raw_xml, rpt_date_label, source_xml_url) {
    doc <- xml2::read_xml(raw_xml)
    xml2::xml_ns_strip(doc)
    
    hdr <- xml2::xml_find_first(doc, ".//Details")
    rows <- xml2::xml_find_all(doc, ".//Details1")
    if (length(rows) == 0) return(tibble::tibble())
    
    tibble::tibble(
      park_name = xml2::xml_attr(hdr, "ParkName"),
      park_code = xml2::xml_attr(hdr, "ParkCode"),
      month_label = xml2::xml_attr(hdr, "Textbox5"),
      month_mm_yyyy = xml2::xml_attr(hdr, "Textbox33"),
      month_start = as.Date(paste0("01 ", xml2::xml_attr(hdr, "Textbox5")), format = "%d %B %Y"),
      rpt_date_param = rpt_date_label,
      metric = xml2::xml_attr(rows, "Field38"),
      this_month = as.numeric(xml2::xml_attr(rows, "Field39")),
      year_to_date = as.numeric(xml2::xml_attr(rows, "Field40")),
      source_xml_url = source_xml_url
    )
  }
  
  # 1) Landing page -> iframe viewer URL
  landing_html <- fetch_text(landing_url)
  iframe_src <- landing_html |>
    xml2::read_html() |>
    rvest::html_element("iframe") |>
    rvest::html_attr("src")
  
  if (is.na(iframe_src)) stop("Could not find iframe src.")
  viewer_url <- if (startsWith(iframe_src, "http")) iframe_src else paste0("https://irma.nps.gov", iframe_src)
  viewer_url <- gsub("&amp;", "&", viewer_url, fixed = TRUE)
  if (!grepl("([?&])Park=", viewer_url)) viewer_url <- paste0(viewer_url, "&Park=", park)
  
  # 2) Get available month labels
  first_html <- fetch_text(viewer_url)
  rpt_dates <- first_html |>
    xml2::read_html() |>
    rvest::html_elements("select#ReportViewer_ctl04_ctl05_ddValue option") |>
    rvest::html_text2() |>
    stringr::str_squish()
  
  rpt_dates <- rpt_dates[!grepl("^<Select", rpt_dates)]
  if (length(rpt_dates) == 0) stop("No RptDate options found.")
  
  # 3) Loop months
  failed_months <- character(0)
  
  out <- purrr::map_dfr(rpt_dates, function(m) {
    if (!quiet) message("Fetching: ", m)
    tryCatch({
      month_url <- set_rpt_date(viewer_url, m)
      month_html <- fetch_text(month_url)
      xml_url <- extract_export_xml_url(month_html)
      xml_raw <- fetch_raw(xml_url)
      parse_special_use_xml(xml_raw, m, xml_url)
    }, error = function(e) {
      failed_months <<- c(failed_months, m)
      if (!quiet) message("  Failed ", m, ": ", conditionMessage(e))
      tibble::tibble()
    })
  })
  
  out <- dplyr::arrange(out, .data$month_start, .data$metric) %>%
    mutate(park_unit = park)
  
  if (save_csv) {
    if (is.null(output_dir)) output_dir <- paste0(tolower(park), "_special_use_downloads")
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    out_csv <- file.path(output_dir, paste0(park, ".csv"))
    readr::write_csv(out, out_csv)
    attr(out, "output_csv") <- out_csv
  }
  
  attr(out, "park") <- park
  attr(out, "failed_months") <- failed_months
  out
}

# parks <- c("ZION", "GLAC", "BRCA", "CEBR") %>%
#   map_dfr(~getSpecialUseVisitation(park = .x))
