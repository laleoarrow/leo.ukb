# helper function to determine default node roles based on DAG structure
dag_node_roles <- function(dag_obj, exposure_name, outcome_name) {
  if (!inherits(dag_obj, "dagitty")) stop("`dag_obj` must inherit from 'dagitty'.", call. = F)
  leo.basic::leo_log("Automated delocate node roles based on their relationships to exposure and outcome.")
  leo.basic::leo_log("You can mannally adapt ...")
  ancestors_exposure <- setdiff(dagitty::ancestors(dag_obj, exposure_name), exposure_name)
  ancestors_outcome <- setdiff(dagitty::ancestors(dag_obj, outcome_name), outcome_name)
  descendants_exposure <- setdiff(dagitty::descendants(dag_obj, exposure_name), exposure_name)

  dplyr::mutate(
    dplyr::tibble(name = names(dag_obj)),
    node_type = dplyr::case_when(
      name == exposure_name ~ "Exposure",
      name == outcome_name ~ "Outcome",
      name %in% descendants_exposure & name %in% ancestors_outcome ~ "Mediator",
      name %in% ancestors_exposure & name %in% ancestors_outcome ~ "Confounder candidate",
      name %in% ancestors_exposure & !(name %in% ancestors_outcome) ~ "Exposure ancestor",
      !(name %in% ancestors_exposure) & name %in% ancestors_outcome & !(name %in% descendants_exposure) ~ "Outcome-only predictor",
      name %in% descendants_exposure & !(name %in% ancestors_outcome) ~ "Exposure descendant",
      TRUE ~ "Other node"
    )
  )
}

#' Create a `leo_dag` object
#'
#' @param dag_obj A `dagitty::dagitty()` object representing the DAG.
#' @param exposure_name Exposure node name.
#' @param outcome_name Outcome node name.
#' @param dag_name DAG name used in output metadata.
#' @param labels Optional named character vector for display labels.
#' @param layout Layout passed to `ggdag::tidy_dagitty()`.
#' @param use_existing_coords If `TRUE`, reuse coordinates stored in `dag_obj` when available.
#' @param seed Layout seed.
#' @return A `leo_dag` object.
#' @export
#' @examples 
#' # # Example usage
# # Example DAG kept for local development
# dag_t2d_total <- dagitty::dagitty(
#   "dag {
#     Age -> Diabetes
#     Age -> Iridocyclitis
#     Sex -> Diabetes
#     Sex -> Iridocyclitis
#     Ethnicity -> Diabetes
#     Ethnicity -> Iridocyclitis
#     TDI -> Lifestyle
#     TDI -> Diabetes
#     TDI -> Iridocyclitis
#     Education -> Lifestyle
#     Education -> Diabetes
#     Education -> Iridocyclitis
#     Lifestyle -> BMI
#     Lifestyle -> Diabetes
#     Lifestyle -> Iridocyclitis
#     BMI -> Diabetes
#     B27 -> Iridocyclitis
#     B27 -> AutoimmuneComorbidity
#     AutoimmuneComorbidity -> Iridocyclitis
#     Diabetes -> Iridocyclitis
#   }"
# )
# exposure_name <- "Diabetes"
# outcome_name <- "Iridocyclitis"
# x <- leo_dag(dag_t2d_total, exposure_name, outcome_name, dag_name = "t2d_total")
# print(x)
# plot(x)
# leo_dag_write(x, figure_dir = "figure", output_dir = "output")
leo_dag <- function(dag_obj, exposure_name, outcome_name, dag_name = NULL, node_roles_list = NULL,
                    labels = NULL, layout = "kk", use_existing_coords = F, seed = 2026) {
  required_pkgs <- c("dagitty", "ggdag", "ggplot2", "dplyr", "tibble")
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = T)]
  if (length(missing_pkgs) > 0) stop("Missing package(s): ", paste(missing_pkgs, collapse = ", "), call. = F)
  if (!inherits(dag_obj, "dagitty")) stop("`dag_obj` must inherit from 'dagitty'.", call. = F)
  if (!is.null(labels) && is.null(names(labels))) stop("`labels` must be a named character vector.", call. = F)

  node_roles <- dag_node_roles(dag_obj, exposure_name, outcome_name)
  dag_name <- if (is.null(dag_name) || dag_name == "") paste0(exposure_name, "_", outcome_name, "_dag") else dag_name
  base_filename <- gsub("[^A-Za-z0-9_]", "_", dag_name)
  dag_tbl <- dplyr::mutate(
    dplyr::left_join(ggdag::pull_dag_data(ggdag::tidy_dagitty(dag_obj, seed = seed, layout = layout)), node_roles, by = "name"),
    node_type = dplyr::if_else(is.na(node_type), "Other node", node_type)
  )
  if (use_existing_coords) {
    coords <- dagitty::coordinates(dag_obj)
    if (length(coords$x) > 0 && length(coords$y) > 0 && !anyNA(coords$x) && !anyNA(coords$y)) {
      dag_tbl$x <- unname(coords$x[dag_tbl$name])
      dag_tbl$y <- unname(coords$y[dag_tbl$name])
      has_to <- !is.na(dag_tbl$to)
      dag_tbl$xend[has_to] <- unname(coords$x[dag_tbl$to[has_to]])
      dag_tbl$yend[has_to] <- unname(coords$y[dag_tbl$to[has_to]])
    }
  }
  dag_tbl$label <- if (is.null(labels)) dag_tbl$name else unname(labels[dag_tbl$name])
  dag_tbl$label[is.na(dag_tbl$label) | dag_tbl$label == ""] <- dag_tbl$name[is.na(dag_tbl$label) | dag_tbl$label == ""]
  plot_obj <- ggplot2::ggplot(dag_tbl, ggplot2::aes(x = x, y = y, xend = xend, yend = yend)) +
    ggdag::geom_dag_edges_link(edge_colour = "black", edge_width = 0.7) +
    ggdag::geom_dag_point(ggplot2::aes(fill = node_type), shape = 21, color = "black", size = 13) +
    ggdag::geom_dag_label_repel(
      ggplot2::aes(label = label),
      color = "black",
      fill = grDevices::adjustcolor("white", alpha.f = 0.7),
      size = 4,
      fontface = "bold",
      box.padding = grid::unit(0.7, "lines"),
      label.padding = grid::unit(0.2, "lines"),
      point.padding = grid::unit(1.1, "lines"),
      label.size = 0.25,
      segment.color = "gray45",
      segment.size = 0.4,
      force = 8,
      max.iter = 5000,
      nudge_y = 0.05,
      show.legend = F
    ) +
    ggplot2::scale_fill_manual(values = c(
      "Exposure" = "#f2614e", "Outcome" = "#14a9ca", "Confounder candidate" = "#34ef5c6a", "Mediator" = "#6232f3",
      "Exposure ancestor" = "#B2ABD2", "Outcome-only predictor" = "#91CF60", "Exposure descendant" = "#F46D43", "Other node" = "#D9D9D9"
    )) +
    ggplot2::labs(title = paste("DAG for", dag_name), fill = "Node role") +
    ggdag::theme_dag() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 15, face = "bold", hjust = 0.5))
  adj_sets <- dagitty::adjustmentSets(dag_obj, exposure = exposure_name, outcome = outcome_name)
  minimal_adjustment_sets <- if (length(adj_sets) == 0) list() else lapply(adj_sets, function(x) as.character(unlist(x)))
  suggested_adjustment_set <- if (length(minimal_adjustment_sets) == 0) NULL else minimal_adjustment_sets[[1]]

  structure(
    list(
      dag_obj = dag_obj,
      exposure_name = exposure_name,
      outcome_name = outcome_name,
      dag_name = dag_name,
      base_filename = base_filename,
      labels = labels,
      layout = layout,
      use_existing_coords = use_existing_coords,
      seed = seed,
      node_roles = node_roles,
      dag_tbl = dag_tbl,
      plot_obj = plot_obj,
      minimal_adjustment_sets = minimal_adjustment_sets,
      suggested_adjustment_set = suggested_adjustment_set
    ),
    class = "leo_dag"
  )
}

#' @export
print.leo_dag <- function(x, ...) {
  cat("<leo_dag>", x$dag_name, "\n")
  cat("Exposure:", x$exposure_name, "\n")
  cat("Outcome :", x$outcome_name, "\n")
  if (length(x$minimal_adjustment_sets) == 0) {
    cat("Minimal adjustment sets: none\n")
  } else if (length(x$minimal_adjustment_sets) == 1) {
    cat("Minimal adjustment set :", if (length(x$minimal_adjustment_sets[[1]]) == 0) "{}" else paste(x$minimal_adjustment_sets[[1]], collapse = ", "), "\n")
  } else {
    cat("Minimal adjustment sets:\n")
    for (i in seq_along(x$minimal_adjustment_sets)) cat(" ", i, ":", if (length(x$minimal_adjustment_sets[[i]]) == 0) "{}" else paste(x$minimal_adjustment_sets[[i]], collapse = ", "), "\n")
  }
  invisible(x)
}

#' @export
summary.leo_dag <- function(object, ...) {
  structure(
    list(
      dag_name = object$dag_name,
      exposure_name = object$exposure_name,
      outcome_name = object$outcome_name,
      n_nodes = nrow(object$node_roles),
      n_adjustment_sets = length(object$minimal_adjustment_sets),
      minimal_adjustment_sets = object$minimal_adjustment_sets,
      suggested_adjustment_set = object$suggested_adjustment_set
    ),
    class = "summary.leo_dag"
  )
}

#' @export
plot.leo_dag <- function(x, ...) {
  x$plot_obj
}

#' Write outputs for a `leo_dag` object
#'
#' @param x A `leo_dag` object.
#' @param ... Reserved for S3 methods.
#' @export
leo_dag_write <- function(x, ...) {
  UseMethod("leo_dag_write")
}

#' Write a `leo_dag` object to disk
#'
#' @inheritParams leo_dag_write
#' @param figure_dir Figure output directory.
#' @param output_dir Table output directory.
#' @param width Figure width.
#' @param height Figure height.
#' @return Invisibly returns output file paths.
#' @export
leo_dag_write.leo_dag <- function(x, figure_dir = "figure/dag", output_dir = "output/dag", width = 8, height = 6, ...) {
  if (!requireNamespace("readr", quietly = T)) stop("Missing package(s): readr", call. = F)
  dir.create(figure_dir, recursive = T, showWarnings = F)
  dir.create(output_dir, recursive = T, showWarnings = F)

  pdf_file <- file.path(figure_dir, paste0(x$base_filename, ".pdf"))
  adj_csv <- file.path(output_dir, paste0(x$base_filename, "_adjustment_sets.csv"))
  roles_csv <- file.path(output_dir, paste0(x$base_filename, "_roles.csv"))
  adj_df <- tibble::tibble(
    adjustment_set = vapply(
      x$minimal_adjustment_sets,
      function(x) if (length(x) == 0) "{}" else paste(x, collapse = ", "),
      FUN.VALUE = character(1)
    )
  )
  if (length(x$minimal_adjustment_sets) == 0) adj_df <- tibble::tibble(adjustment_set = character())

  ggplot2::ggsave(pdf_file, x$plot_obj, width = width, height = height)
  readr::write_csv(adj_df, adj_csv)
  readr::write_csv(x$node_roles, roles_csv)

  invisible(list(pdf = pdf_file, adjustment_sets = adj_csv, roles = roles_csv))
}
