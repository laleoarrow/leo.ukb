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
  c(
    "Survival Analysis Step 1: Main Cox Analysis",
    "Survival Analysis Step 2: Subgroup and Interaction",
    "生存分析 Step 1：主 Cox 分析",
    "生存分析 Step 2：亚组分析与交互效应"
  ),
  c(
    file.path(site_dir, "articles", "survival-analysis-step1.html"),
    file.path(site_dir, "articles", "survival-analysis-step2.html"),
    file.path(site_dir, "articles", "survival-analysis-step1-zh.html"),
    file.path(site_dir, "articles", "survival-analysis-step2-zh.html")
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
