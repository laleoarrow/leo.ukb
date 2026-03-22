pkgdown::build_site_github_pages(dest_dir = "docs-site", new_process = F, install = TRUE)

zh_articles <- list.files("docs-site/articles", pattern = "-zh\\.html$", full.names = TRUE)

for (file in zh_articles) {
  html <- readLines(file, warn = F, encoding = "UTF-8")
  html <- sub("<html lang=\"en\">", "<html lang=\"zh-CN\">", html, fixed = TRUE)
  writeLines(html, file, useBytes = TRUE)
}
