# Regenerate man/figures/downloads.png from CRAN download logs.
# Rerun periodically (e.g. before a release) to refresh the README figure.

daily <- jsonlite::fromJSON(
  "https://cranlogs.r-pkg.org/downloads/daily/2000-01-01:2100-01-01/tReeTraits"
)$downloads[[1]]
stopifnot(is.data.frame(daily), nrow(daily) > 0)  # fail the run, don't commit a broken figure
daily$day <- as.Date(daily$day)

archive_date <- as.Date("2026-06-09")

week_start <- function(x) as.Date(cut(x, "week", start.on.monday = FALSE))

weekly <- daily |>
  dplyr::mutate(week = week_start(day)) |>
  dplyr::group_by(week) |>
  dplyr::summarize(downloads = sum(downloads), .groups = "drop")

# cranlogs omits days with no downloads, so fill in zero-download weeks
# through the current one; otherwise the chart stops at the last download.
weekly <- data.frame(
  week = seq(min(weekly$week), week_start(max(Sys.Date(), max(daily$day))), by = "7 days")
) |>
  dplyr::left_join(weekly, by = "week") |>
  dplyr::mutate(downloads = dplyr::coalesce(downloads, 0L))

p <- ggplot2::ggplot(weekly, ggplot2::aes(week, downloads)) +
  ggplot2::geom_col(fill = "#2a78d6", width = 6) +
  ggplot2::geom_vline(xintercept = archive_date,
                       linetype = "dashed", color = "#eb6834") +
  ggplot2::annotate("text", x = archive_date, y = max(weekly$downloads) * 1.05,
                     label = "Archived", color = "#eb6834", size = 3, hjust = 1) +
  ggplot2::labs(x = NULL, y = "Weekly downloads",
                title = "tReeTraits CRAN downloads") +
  ggplot2::theme_bw() +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                 panel.grid.major.x = ggplot2::element_blank())

ggplot2::ggsave("man/figures/downloads.png", p, width = 7, height = 3, dpi = 200)
