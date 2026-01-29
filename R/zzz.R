.onAttach <- function(libname, pkgname) {
  banner <- paste(
    "  ======================================================================",
    "      ██╗     ███████╗ ██████╗            ██╗   ██╗██╗  ██╗██████╗",
    "      ██║     ██╔════╝██╔═══██╗           ██║   ██║██║ ██╔╝██╔══██╗",
    "      ██║     █████╗  ██║   ██║    ██     ██║   ██║█████╔╝ ██████╔╝",
    "      ██║     ██╔══╝  ██║   ██║           ██║   ██║██╔═██╗ ██╔══██╗",
    "      ███████╗███████╗╚██████╔╝           ╚██████╔╝██║  ██╗██████╔╝",
    "      ╚══════╝╚══════╝ ╚═════╝            ╚═════╝ ╚═╝  ╚═╝╚═════╝",
    "  ======================================================================",
    sep = "\n"
  )

  version <- as.character(utils::packageVersion(pkgname))
  banner_width <- max(nchar(strsplit(banner, "\n", fixed = TRUE)[[1]]), 0)
  welcome_line <- "Welcome to leo.ukb!"
  info_line <- paste0("leo.ukb v", version, " | UK Biobank toolkit")
  time_line <- paste0("System time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  meta_lines <- vapply(
    c(welcome_line, info_line, time_line),
    .ukb_center_line,
    width = banner_width,
    FUN.VALUE = character(1)
  )

  startup_text <- paste(c(banner, meta_lines), collapse = "\n")
  .ukb_startup_message(startup_text)
}

.ukb_startup_message <- function(text) {
  if (.ukb_in_rstudio()) {
    leo.basic::leo_message(text, color = "94")
    return(invisible(NULL))
  }
  if (.ukb_supports_truecolor()) {
    packageStartupMessage(.ukb_gradient(text))
  } else {
    packageStartupMessage(text)
  }
}

.ukb_in_rstudio <- function() {
  nzchar(Sys.getenv("RSTUDIO"))
}

.ukb_supports_truecolor <- function() {
  if (nzchar(Sys.getenv("NO_COLOR"))) {
    return(FALSE)
  }
  term <- Sys.getenv("TERM")
  if (identical(term, "dumb")) {
    return(FALSE)
  }
  colorterm <- Sys.getenv("COLORTERM")
  if (nzchar(colorterm) && grepl("truecolor|24bit", colorterm, ignore.case = TRUE)) {
    return(TRUE)
  }
  FALSE
}

.ukb_center_line <- function(text, width) {
  if (width <= 0) {
    return(text)
  }
  pad <- floor((width - nchar(text)) / 2)
  if (pad < 0) {
    pad <- 0
  }
  paste0(strrep(" ", pad), text)
}

.ukb_gradient <- function(text,
                          start = c(80, 200, 255),
                          end = c(0, 80, 255)) {
  lines <- strsplit(text, "\n", fixed = TRUE)[[1]]
  colored <- vapply(lines, function(line) {
    if (!nzchar(line)) {
      return("")
    }

    chars <- strsplit(line, "", fixed = TRUE)[[1]]
    n <- length(chars)
    if (n == 1L) {
      t <- 0
    } else {
      t <- (seq_len(n) - 1) / (n - 1)
    }

    r <- round(start[1] + (end[1] - start[1]) * t)
    g <- round(start[2] + (end[2] - start[2]) * t)
    b <- round(start[3] + (end[3] - start[3]) * t)

    paste0(
      paste0(sprintf("\033[38;2;%d;%d;%dm%s", r, g, b, chars), collapse = ""),
      "\033[0m"
    )
  }, character(1))

  paste(colored, collapse = "\n")
}
