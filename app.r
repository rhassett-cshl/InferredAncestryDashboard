library(magrittr)
library(htmltools)
source("database.R")

ui <- fluidPage(
  # Header
  headerPanel(
    title="Inferred Ancestry Dashboard"
  ),

  #,

  # Input widgets
  tabsetPanel(
    id="main_tabs",
    tabPanel(
      "Map",
      column(
        12,
        shinycssloaders::withSpinner(
          leaflet::leafletOutput(
            "map",
            height = "calc(100vh - 120px)"
          ),
          size=2,
          color="#0080b7"
        )
      )
    ),
    tabPanel(
      "Table",
      column(
        12,
        h4("Click a site"),
        downloadButton("download_table_csv", "Export table to CSV"),
        div(
          DT::dataTableOutput("table_input"),
          style="font-size:70%"
        )
      )
    ),
    tabPanel(
      "Plots"
    ),
    tabPanel(
      "User guide",
      fluidRow(
        column(
          8,
          includeMarkdown('./user_guide/user_guide.rmd')
        )
      )
    )
  )
)

server <- function(input, output, session){

  # Loading modal to keep user out of trouble while ancestryCall loads
  showModal(modalDialog(
    title = "LOADING - PLEASE WAIT...",
    "Please wait for ancestry data to load before proceeding.",
    size = "l",
    footer = NULL
  ))

  # Remove modal when ancestryCall (loadTableData) is ready
  ancestry_ready <- reactiveValues(loaded = FALSE)
  observe({
    req(ancestry_ready$loaded)
    removeModal()
  })

  # Extract ancestry calls with population definition coordinates
  ancestryCall <- loadTableData()
  ancestry_ready$loaded <- TRUE
  # Human-friendly labels for ancestryCall columns used elsewhere
  labels_map <- c(
    experiment = "Experiment",
    bioSample = "BioSample",
    bioProject = "BioProject",
    libraryStrategy = "Library Strategy",
    populationDefinition_name = "Population Definition",
    accuracy = "Accuracy",
    accuracyQuantifier = "Accuracy Quantifier",
    inferenceMethodProperties_name = "Inference Method",
    latitude = "latitude",
    longitude = "longitude"
  )
  cols_to_rename <- intersect(names(labels_map), names(ancestryCall))
  names(ancestryCall)[match(cols_to_rename, names(ancestryCall))] <-
    labels_map[cols_to_rename]

  # Empty reactive values object
  reactive_objects=reactiveValues()

  # Resources for returning site info on click:
  ## https://stackoverflow.com/questions/28938642/marker-mouse-click-event-in-r-leaflet-for-shiny
  ## https://stackoverflow.com/questions/42613984/how-to-implement-inputmap-marker-click-correctly?noredirect=1&lq=1


  # Map: render once ancestryCall is loaded
  output$map <- leaflet::renderLeaflet({
    m <- build_ancestry_markers(ancestryCall)
    lng_r <- range(m$lng, na.rm = TRUE)
    lat_r <- range(m$lat, na.rm = TRUE)
    if (length(m$layerId) == 0) {
      return(leaflet::leaflet() %>%
        leaflet::addTiles(options = leaflet::tileOptions(noWrap = TRUE)))
    }
    if (diff(lng_r) == 0) lng_r <- lng_r + c(-0.5, 0.5)
    if (diff(lat_r) == 0) lat_r <- lat_r + c(-0.5, 0.5)
    # Slight margin so fitBounds has room; restrict panning to avoid gray
    margin <- 0.15 * max(diff(lng_r), diff(lat_r), 1)
    bounds <- list(
      c(lat_r[1] - margin, lng_r[1] - margin),
      c(lat_r[2] + margin, lng_r[2] + margin)
    )
    m_base <- leaflet::leaflet(options =
      leaflet::leafletOptions(
        maxBounds = bounds,
        maxBoundsViscosity = 1
      )
    ) %>%
      leaflet::addTiles(options = leaflet::tileOptions(noWrap = TRUE))
    lng_range <- lng_r
    lat_range <- lat_r
    m_base %>%
      leaflet::addCircleMarkers(
        lng = m$lng,
        lat = m$lat,
        layerId = m$layerId,
        radius = 5,
        fillOpacity = 0.8,
        stroke = FALSE,
        color = "#2C7BB6",
        label = m$label
      ) %>%
      leaflet::fitBounds(
        lng1 = lng_range[1],
        lat1 = lat_range[1],
        lng2 = lng_range[2],
        lat2 = lat_range[2],
        options = list(padding = c(40, 40))
      )
  })
  
  # ancestryCall already loaded and labeled above

  # Format Accuracy as percentage (two decimals)
  if ("Accuracy" %in% names(ancestryCall)){
    acc <- suppressWarnings(
      as.numeric(ancestryCall$Accuracy)
    )
    is_num <- !is.na(acc) & is.finite(acc)
    if (any(is_num)){
      max_v <- max(acc[is_num], na.rm = TRUE)
      scale_factor <- if (max_v <= 1) 100 else 1
      fmt <- rep(NA_character_, length(acc))
      fmt[is_num] <- sprintf(
        "%.2f%%",
        acc[is_num] * scale_factor
      )
      if (any(!is_num)){
        fmt[!is_num] <- as.character(
          ancestryCall$Accuracy[!is_num]
        )
      }
      ancestryCall$Accuracy <- fmt
    }
  }

  # Turn Experiment values into SRA links while displaying same text
  if ("Experiment" %in% names(ancestryCall)){
    ancestryCall$Experiment <- vapply(
      as.character(ancestryCall$Experiment),
      function(val){
        if (is.na(val) || nzchar(val) == FALSE) return("")
        url <- paste0(
          "https://www.ncbi.nlm.nih.gov/sra/?term=",
          utils::URLencode(val, reserved = TRUE)
        )
        paste0(
          '<a href="', url,
          '" target="_blank" rel="noopener">', val, "</a>"
        )
      },
      FUN.VALUE = character(1)
    )
  }

  # Table interface
  output$table_input=DT::renderDataTable({
    # Do not escape the Experiment column so links render; escape others
    escape_cols <- seq_len(ncol(ancestryCall))
    exp_idx <- which(names(ancestryCall) == "Experiment")
    if (length(exp_idx) == 1){
      escape_cols <- setdiff(escape_cols, exp_idx)
    }

    DT::datatable(
      ancestryCall,
      escape = escape_cols,
      selection = "multiple",
      rownames = FALSE,
      filter="top",
      options = list(
        dom = "ltipr",
        scrollX = TRUE,
        autoWidth = TRUE,
        paging = FALSE,
        searching = TRUE,
        orderCellsTop = TRUE,
        columnDefs = list(
          list(targets = "_all", className = "dt-left"),
          list(targets = "_all", searchable = TRUE)
        )
      )
    )
  }, server = FALSE)

  # Ancestry rows to show on map: all with coords, or only selected table rows
  map_ancestry_display <- reactive({
    ac <- ancestryCall
    sel <- input$table_input_rows_selected
    if (length(sel) > 0) ac <- ac[sel, , drop = FALSE]
    lat_col <- "latitude"
    lng_col <- "longitude"
    if (!lat_col %in% names(ac) || !lng_col %in% names(ac))
      return(ac[0, , drop = FALSE])
    lng_raw <- vapply(ac[[lng_col]], function(x)
      as.numeric(unlist(x))[1], numeric(1))
    lat_raw <- vapply(ac[[lat_col]], function(x)
      as.numeric(unlist(x))[1], numeric(1))
    keep <- is.finite(lng_raw) & is.finite(lat_raw)
    ac[keep, , drop = FALSE]
  })

  # Build lng, lat, layerId, label for map markers; label = HTML summary
  build_ancestry_markers <- function(ac) {
    lat_col <- "latitude"
    lng_col <- "longitude"
    if (!lat_col %in% names(ac) || !lng_col %in% names(ac) || nrow(ac) == 0)
      return(list(lng = numeric(0), lat = numeric(0),
        layerId = character(0), label = list()))
    lng_raw <- vapply(ac[[lng_col]], function(x)
      as.numeric(unlist(x))[1], numeric(1))
    lat_raw <- vapply(ac[[lat_col]], function(x)
      as.numeric(unlist(x))[1], numeric(1))
    keep <- is.finite(lng_raw) & is.finite(lat_raw)
    ac <- ac[keep, , drop = FALSE]
    lng_vals <- lng_raw[keep]
    lat_vals <- lat_raw[keep]
    label_col <- if ("Experiment" %in% names(ac)) "Experiment" else
      if ("experiment" %in% names(ac)) "experiment" else NULL
    ids <- if (is.null(label_col)) as.character(seq_len(nrow(ac))) else
      as.character(ac[[label_col]])
    summary_cols <- c(
      "Experiment", "Population Definition", "Accuracy",
      "Accuracy Quantifier", "BioSample", "BioProject"
    )
    summary_cols <- intersect(summary_cols, names(ac))
    labels <- lapply(seq_len(nrow(ac)), function(i) {
      row <- ac[i, ]
      parts <- vapply(summary_cols, function(c) {
        paste0("<b>", c, ":</b> ", as.character(row[[c]])[1])
      }, character(1))
      htmltools::HTML(paste(parts, collapse = "<br/>"))
    })
    list(lng = lng_vals, lat = lat_vals, layerId = ids, label = labels)
  }

  # Export Table filtered/sorted data to CSV
  output$download_table_csv <- downloadHandler(
    filename = function() {
      paste0("ancestry_table_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      dat <- ancestryCall
      state <- input$table_input_state
      if (!is.null(state)) {
        # Apply global search
        global <- state$search$search
        if (is.character(global) && nzchar(trimws(global))) {
          pat <- gsub(" ", ".*", trimws(global))
          keep <- apply(dat, 1, function(r) any(grepl(pat, as.character(r),
            ignore.case = TRUE)))
          dat <- dat[keep, , drop = FALSE]
        }
        # Apply column-level search
        cols <- state$columns
        if (is.list(cols)) {
          for (j in seq_along(cols)) {
            if (j > ncol(dat)) next
            s <- cols[[j]]$search$search
            if (is.character(s) && nzchar(trimws(s))) {
              pat <- trimws(s)
              keep <- grepl(pat, as.character(dat[[j]]), ignore.case = TRUE)
              dat <- dat[keep, , drop = FALSE]
            }
          }
        }
        # Apply sort order
        ord <- state$order
        if (is.list(ord) && length(ord) > 0 && nrow(dat) > 0) {
          ord_cols <- c()
          ord_dec <- c()
          for (o in ord) {
            col_idx <- as.integer(o[[1]]) + 1L
            if (col_idx >= 1 && col_idx <= ncol(dat)) {
              ord_cols <- c(ord_cols, col_idx)
              ord_dec <- c(ord_dec, identical(o[[2]], "desc"))
            }
          }
          if (length(ord_cols) > 0) {
            ord_args <- lapply(seq_along(ord_cols), function(i) {
              x <- dat[[ord_cols[i]]]
              if (ord_dec[i]) -xtfrm(x) else xtfrm(x)
            })
            dat <- dat[do.call(order, ord_args), , drop = FALSE]
          }
        }
      }
      utils::write.csv(dat, file, row.names = FALSE)
    }
  )

  # Fit map to marker bounds when Map tab is shown (map has correct size then)
  map_bounds <- reactive({
    m <- build_ancestry_markers(map_ancestry_display())
    if (length(m$layerId) == 0) return(NULL)
    lng_r <- range(m$lng, na.rm = TRUE)
    lat_r <- range(m$lat, na.rm = TRUE)
    if (diff(lng_r) == 0) lng_r <- lng_r + c(-0.01, 0.01)
    if (diff(lat_r) == 0) lat_r <- lat_r + c(-0.01, 0.01)
    list(lng1 = lng_r[1], lat1 = lat_r[1], lng2 = lng_r[2], lat2 = lat_r[2])
  })
  observeEvent(input$main_tabs, {
    if (input$main_tabs != "Map") return()
    b <- map_bounds()
    if (is.null(b)) return()
    leaflet::leafletProxy("map") %>%
      leaflet::fitBounds(
        lng1 = b$lng1,
        lat1 = b$lat1,
        lng2 = b$lng2,
        lat2 = b$lat2,
        options = list(padding = c(20, 20))
      )
  }, ignoreInit = FALSE)

  # Table selection: show only selected rows on map (or all if none selected)
  observeEvent(input$table_input_rows_selected, {
    md <- map_ancestry_display()
    m <- build_ancestry_markers(md)
    proxy <- leaflet::leafletProxy("map")
    proxy %>% leaflet::clearMarkers()
    if (length(m$layerId) > 0) {
      proxy %>%
        leaflet::addCircleMarkers(
          lng = m$lng,
          lat = m$lat,
          layerId = m$layerId,
          radius = 5,
          fillOpacity = 0.8,
          stroke = FALSE,
          color = "#2C7BB6",
          label = m$label
        )
    }
  }, ignoreNULL = FALSE)

  # Filter Table to match clicked sites from map
  input_table_proxy = DT::dataTableProxy('table_input')
  observeEvent(input$map_marker_click,{
    id <- input$map_marker_click$id
    if (is.null(id)) return()
    # If the click corresponds to an ancestryCall Experiment, filter Table by it
    exp_col <- if ("Experiment" %in% names(ancestryCall)) "Experiment" else if ("experiment" %in% names(ancestryCall)) "experiment" else NULL
    if (!is.null(exp_col) && id %in% ancestryCall[[exp_col]]){
      input_table_proxy %>% DT::clearSearch() %>% DT::updateSearch(
        keywords = list(
          global = paste(id)
        )
      )
    } else {
      # Otherwise, just clear global search
      input_table_proxy %>% DT::clearSearch() %>% DT::updateSearch(
        keywords = list(global = "")
      )
    }
  })


  ## Map polygon click
  #observe({
  #  req(profiles_long)
  #  au_click <- input$map_shape_click
  #  #if (is.null(au_click)){return()}
  #  auid=au_click$id
  #  print(au_click$id)
  #})

}

## run app
shinyApp(ui = ui, server = server)
