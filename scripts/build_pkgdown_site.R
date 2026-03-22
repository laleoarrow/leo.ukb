devtools::build_vignettes(pkg = ".")

pkgdown::build_site(
  pkg = ".",
  install = TRUE,
  new_process = FALSE,
  override = list(destination = "docs-site")
)

zh_articles <- list.files("docs-site/articles", pattern = "-zh\\.html$", full.names = TRUE)

for (file in zh_articles) {
  html <- readLines(file, warn = FALSE, encoding = "UTF-8")
  html <- sub("<html lang=\"en\">", "<html lang=\"zh-CN\">", html, fixed = TRUE)
  writeLines(html, file, useBytes = TRUE)
}

article_checks <- c(
  "docs-site/articles/clinical-regression.html" = "Clinical Regression Workflows",
  "docs-site/articles/clinical-regression-zh.html" = "临床回归分析工作流"
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
