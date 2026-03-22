site_dir <- if (identical(Sys.getenv("GITHUB_ACTIONS"), "true")) {
  "docs-site"
} else {
  file.path(tempdir(), "leo-ukb-pkgdown-site")
}

pkgdown::build_site(
  pkg = ".",
  install = TRUE,
  new_process = FALSE,
  override = list(destination = site_dir)
)

zh_articles <- list.files(file.path(site_dir, "articles"), pattern = "-zh\\.html$", full.names = TRUE)

for (file in zh_articles) {
  html <- readLines(file, warn = FALSE, encoding = "UTF-8")
  html <- sub("<html lang=\"en\">", "<html lang=\"zh-CN\">", html, fixed = TRUE)
  writeLines(html, file, useBytes = TRUE)
}

article_checks <- stats::setNames(
  c("Survival Analysis with leo.ukb", "leo.ukb 生存分析教程"),
  c(
    file.path(site_dir, "articles", "survival-analysis.html"),
    file.path(site_dir, "articles", "survival-analysis-zh.html")
  )
)

for (file in names(article_checks)) {
  if (!file.exists(file)) {
    stop("pkgdown article was not created: ", file, call. = FALSE)
  }

  html <- readLines(file, warn = FALSE, encoding = "UTF-8")
  if (!grepl(article_checks[[file]], paste(html, collapse = "\n"), fixed = TRUE)) {
    stop(
      "pkgdown article build check failed for ", file,
      ". Expected to find article title/content marker: ", article_checks[[file]],
      call. = FALSE
    )
  }
}

message("pkgdown site written to: ", site_dir)
