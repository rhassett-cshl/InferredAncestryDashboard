library(pool)
library(dplyr)

loadTableData <- function() {
  db <- dbPool(
    RSQLite::SQLite(),
    dbname = "C:/Users/hassett/Documents/AncestryInferenceApp/ancestryinference.db"
  )

  on.exit(poolClose(db))

  base_name <- "ancestryCall"
  base_tbl <- db %>% tbl(base_name)

  fk <- DBI::dbGetQuery(
    db,
    sprintf("PRAGMA foreign_key_list('%s')", base_name)
  )

  join_parent <- function(x, parent) {
    if (!nrow(fk)) return(x)
    hit <- fk[fk$table == parent, , drop = FALSE]
    if (!nrow(hit)) return(x)
    by_map <- stats::setNames(hit$to, hit$from)
    y <- db %>% tbl(parent)
    dplyr::left_join(
      x,
      y,
      by = by_map,
      suffix = c("", paste0("_", parent))
    )
  }

  base_tbl <- join_parent(base_tbl, "molecularProfile")
  base_tbl <- join_parent(base_tbl, "inferenceMethodProperties")
  base_tbl <- join_parent(base_tbl, "populationDefinition")

  selected <- base_tbl %>%
    dplyr::select(dplyr::any_of(c(
      "populationDefinitionId",
      "experiment",
      "bioSample",
      "bioProject",
      "libraryStrategy",
      "name_populationDefinition",
      "accuracy",
      "accuracyQuantifier",
      "name_inferenceMethodProperties"
    ))) %>%
    collect()

  coords <- DBI::dbGetQuery(
    db,
    "SELECT populationDefinitionId, latitude, longitude
     FROM populationDefinition"
  )
  selected[["populationDefinitionId"]] <- as.integer(selected[["populationDefinitionId"]])
  coords[["populationDefinitionId"]] <- as.integer(coords[["populationDefinitionId"]])

  selected <- dplyr::left_join(
    selected,
    coords,
    by = "populationDefinitionId"
  )
  selected <- dplyr::select(selected, -dplyr::any_of("populationDefinitionId"))

  selected <- selected %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::any_of(c(
          "bioProject",
          "libraryStrategy",
          "inferenceMethodProperties_name",
          "populationDefinition_name"
        )),
        ~ as.character(.)
      ),
      latitude = suppressWarnings(as.numeric(as.character(.data$latitude))),
      longitude = suppressWarnings(as.numeric(as.character(.data$longitude)))
    )

  return(as.data.frame(selected))
}
