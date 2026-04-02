library(pool)
library(dplyr)

getPopSpecificAccuracyForAncestryCall <- function(ancestry_call_id) {
  db <- dbPool(
    RSQLite::SQLite(),
    dbname = "C:/Users/hassett/Documents/AncestryInferenceApp/ancestryinference_v0_8.db"
  )
  on.exit(poolClose(db))

  ancestry_call_id <- suppressWarnings(as.integer(ancestry_call_id)[1])
  if (is.na(ancestry_call_id)) {
    return(data.frame())
  }

  pop_cols <- DBI::dbListFields(db, "popSpecificAccuracy")
  pop_id_col <- if ("populationDefinitionId" %in% pop_cols) {
    "populationDefinitionId"
  } else if ("populationDefinitonId" %in% pop_cols) {
    "populationDefinitonId"
  } else {
    return(data.frame())
  }

  has_ci_lower <- "CILowerBound" %in% pop_cols
  has_ci_upper <- "CIUpperBound" %in% pop_cols

  ci_lower_sql <- if (has_ci_lower) "p.CILowerBound" else "NULL"
  ci_upper_sql <- if (has_ci_upper) "p.CIUpperBound" else "NULL"

  qry <- sprintf(
    paste(
      "SELECT pd.name AS populationDefinition,",
      "p.accuracy AS accuracy,",
      "%s AS CILowerBound,",
      "%s AS CIUpperBound",
      "FROM popSpecificAccuracy p",
      "INNER JOIN (",
      "  SELECT %s AS pop_id,",
      "  MIN(popSpecificAccuracyId) AS first_id",
      "  FROM popSpecificAccuracy",
      "  WHERE ancestryCallId = ?",
      "  GROUP BY %s",
      ") first_per_pop",
      "  ON p.popSpecificAccuracyId = first_per_pop.first_id",
      "LEFT JOIN populationDefinition pd",
      "  ON pd.populationDefinitionId = first_per_pop.pop_id",
      "ORDER BY pd.name"
    ),
    ci_lower_sql,
    ci_upper_sql,
    pop_id_col,
    pop_id_col
  )

  out <- DBI::dbGetQuery(db, qry, params = list(ancestry_call_id))
  if (!nrow(out)) return(data.frame())
  as.data.frame(out)
}

getAncestryCallAccuracyForPlot <- function(ancestry_call_id) {
  db <- dbPool(
    RSQLite::SQLite(),
    dbname = "C:/Users/hassett/Documents/AncestryInferenceApp/ancestryinference_v0_8.db"
  )
  on.exit(poolClose(db))

  ancestry_call_id <- suppressWarnings(as.integer(ancestry_call_id)[1])
  if (is.na(ancestry_call_id)) {
    return(data.frame())
  }

  ac_cols <- DBI::dbListFields(db, "ancestryCall")
  ci_lower_col <- if ("CILowerBound" %in% ac_cols) "CILowerBound" else NULL
  ci_upper_col <- if ("CIUpperBound" %in% ac_cols) "CIUpperBound" else NULL

  bio_col <- NULL
  if ("bioSample" %in% ac_cols) {
    bio_col <- "bioSample"
  }

  tables <- DBI::dbListTables(db)
  mp_join_sql <- ""
  if (is.null(bio_col) &&
      "molecularProfileId" %in% ac_cols &&
      "molecularProfile" %in% tables) {
    mp_cols <- DBI::dbListFields(db, "molecularProfile")
    if ("bioSample" %in% mp_cols) {
      bio_col <- "mp.bioSample"
      mp_join_sql <- paste(
        "LEFT JOIN molecularProfile mp",
        "ON a.molecularProfileId = mp.molecularProfileId"
      )
    }
  }

  bio_sql <- if (is.null(bio_col)) "NULL" else bio_col
  ci_lower_sql <- if (is.null(ci_lower_col)) "NULL" else
    paste0("a.", ci_lower_col)
  ci_upper_sql <- if (is.null(ci_upper_col)) "NULL" else
    paste0("a.", ci_upper_col)

  qry <- sprintf(
    paste(
      "SELECT",
      "a.ancestryCallId AS ancestryCallId,",
      "a.accuracy AS accuracy,",
      "%s AS CILowerBound,",
      "%s AS CIUpperBound,",
      "%s AS bioSample",
      "FROM ancestryCall a",
      "%s",
      "WHERE a.ancestryCallId = ?",
      "LIMIT 1"
    ),
    ci_lower_sql,
    ci_upper_sql,
    bio_sql,
    mp_join_sql
  )

  out <- DBI::dbGetQuery(db, qry, params = list(ancestry_call_id))
  if (!nrow(out)) return(data.frame())
  as.data.frame(out)
}

loadTableData <- function() {
  db <- dbPool(
    RSQLite::SQLite(),
    dbname = "C:/Users/hassett/Documents/AncestryInferenceApp/ancestryinference_v0_8.db"
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
    # Print join columns: child_col = parent_col
    message("join_parent(", parent, "): ", paste(
      names(by_map), "=", by_map, collapse = ", "
    ))
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
  base_tbl <- join_parent(base_tbl, "source")
  base_tbl <- join_parent(base_tbl, "popSpecificAccuracy")

  src <- db %>% tbl("source") %>%
    dplyr::select(dplyr::all_of(c("sourceId", "database")))
  base_tbl <- base_tbl %>%
    dplyr::left_join(src, by = "sourceId")
 
  # populationResolution links via populationDefinition.populationResolutionId
  base_tbl <- base_tbl %>%
    dplyr::left_join(
      db %>% tbl("populationResolution"),
      by = "populationResolutionId",
      suffix = c("", "_populationResolution")
    )

  print(base_tbl)

  selected <- base_tbl %>%
    dplyr::select(dplyr::any_of(c(
      "ancestryCallId",
      "populationDefinitionId",
      "experiment",
      "bioSample",
      "bioProject",
      "libraryStrategy",
      "tissueType",
      "name_populationDefinition",
      "name_populationResolution",
      "accuracy",
      "accuracyQuantifier",
      "name_inferenceMethodProperties",
      "database"
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
          "populationDefinition_name",
          "name_populationResolution"
        )),
        ~ as.character(.)
      ),
      latitude = suppressWarnings(as.numeric(as.character(.data$latitude))),
      longitude = suppressWarnings(as.numeric(as.character(.data$longitude)))
    )

  return(as.data.frame(selected))
}
