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
