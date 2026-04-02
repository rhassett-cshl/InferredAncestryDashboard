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
      "Plots",
      column(
        12,
        h4("Population-specific accuracy for selected ancestry call"),
        p("Select one or more rows in the Table tab."),
        plotOutput("pop_specific_accuracy_plot", height = "66vh"),
        fluidRow(
          column(
            12,
            div(
              style = paste(
                "display:flex;",
                "gap:10px;",
                "align-items:center;",
                "justify-content:center;",
                "margin-top:8px;"
              ),
              actionButton("plot_prev", "\u2190 Previous"),
              strong(textOutput("plot_nav_status", inline = TRUE)),
              actionButton("plot_next", "Next \u2192")
            )
          )
        )
      )
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
    tissueType = "Tissue Type",
    name_populationDefinition = "Population Definition",
    name_populationResolution = "Population Resolution",
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
    m_base <- leaflet::leaflet() %>%
      leaflet::addTiles(options = leaflet::tileOptions(noWrap = TRUE))
    lng_range <- lng_r
    lat_range <- lat_r
    out <- m_base %>%
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
    if (length(m$clusters$lng) > 0) {
      out <- out %>%
        leaflet::addCircleMarkers(
          lng = m$clusters$lng,
          lat = m$clusters$lat,
          layerId = paste0("cluster_", seq_along(m$clusters$lng)),
          radius = 10,
          fillOpacity = 0.9,
          stroke = TRUE,
          weight = 2,
          color = "#333333",
          fillColor = "#E67E22",
          popup = m$clusters$popup
        )
    }
    out %>%
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

  # Show Experiment as SRA link only when source database contains "SRA"
  if ("Experiment" %in% names(ancestryCall)) {
    db_col <- if ("database" %in% names(ancestryCall)) "database" else NULL
    ancestryCall$Experiment <- vapply(
      seq_len(nrow(ancestryCall)),
      function(i) {
        val <- as.character(ancestryCall$Experiment[i])
        if (is.na(val) || !nzchar(trimws(val))) return("")
        is_sra <- FALSE
        if (!is.null(db_col)) {
          db_val <- as.character(ancestryCall[[db_col]][i])
          is_sra <- !is.na(db_val) && grepl("SRA", db_val, ignore.case = TRUE)
        }
        if (!is_sra) return(val)
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

  # Table interface (exclude lat/long and database from display)
  table_cols <- setdiff(
    names(ancestryCall),
    c("ancestryCallId", "latitude", "longitude", "database")
  )
  output$table_input=DT::renderDataTable({
    tbl <- ancestryCall[, table_cols, drop = FALSE]
    # Force all columns to character so DT uses text filter (searchable) in each
    for (col in names(tbl)) {
      if (col != "Experiment")
        tbl[[col]] <- as.character(tbl[[col]])
    }
    # Do not escape the Experiment column so links render; escape others
    escape_cols <- seq_len(ncol(tbl))
    exp_idx <- which(names(tbl) == "Experiment")
    if (length(exp_idx) == 1){
      escape_cols <- setdiff(escape_cols, exp_idx)
    }

    DT::datatable(
      tbl,
      escape = escape_cols,
      selection = "multiple",
      rownames = FALSE,
      filter = list(position = "top", clear = FALSE),
      options = list(
        dom = "ltipr",
        scrollX = TRUE,
        autoWidth = TRUE,
        paging = FALSE,
        searching = TRUE,
        orderCellsTop = TRUE,
        search = list(regex = FALSE, caseInsensitive = TRUE),
        columnDefs = list(
          list(targets = "_all", className = "dt-left"),
          list(targets = "_all", searchable = TRUE)
        ),
        initComplete = htmlwidgets::JS(
          "function(settings, json) {",
          "  var api = this.api();",
          "  var n = api.columns().count();",
          "  var $container = $(api.table().container());",
          "  var $filterRow = $container.find('thead tr').filter(function() {",
          "    return $(this).find('td').length > 0;",
          "  }).first();",
          "  if ($filterRow.length === 0) $filterRow = $container.find('thead tr').eq(1);",
          "  $filterRow.find('td').each(function(i) {",
          "    if (i >= n) return;",
          "    var $cell = $(this);",
          "    var oldVal = $cell.find('input').val() || $cell.find('select').val() || '';",
          "    $cell.empty();",
          "    var $input = $('<input type=\"text\" placeholder=\"Search\">');",
          "    $input.val(oldVal);",
          "    $input.css('width', '100%');",
          "    $cell.append($input);",
          "    $input.on('keyup change', function() {",
          "      api.column(i).search(this.value).draw();",
          "    });",
          "  });",
          "}"
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

  # Build lng, lat, layerId, label for map markers; label = HTML summary.
  # Clusters (overlapping points) get a center marker with popup (Previous/Next).
  build_ancestry_markers <- function(ac) {
    lat_col <- "latitude"
    lng_col <- "longitude"
    if (!lat_col %in% names(ac) || !lng_col %in% names(ac) || nrow(ac) == 0)
      return(list(lng = numeric(0), lat = numeric(0),
        layerId = character(0), label = list(),
        clusters = list(lng = numeric(0), lat = numeric(0), popup = list())))
    lng_raw <- vapply(ac[[lng_col]], function(x)
      as.numeric(unlist(x))[1], numeric(1))
    lat_raw <- vapply(ac[[lat_col]], function(x)
      as.numeric(unlist(x))[1], numeric(1))
    keep <- is.finite(lng_raw) & is.finite(lat_raw)
    ac <- ac[keep, , drop = FALSE]
    lng_vals <- lng_raw[keep]
    lat_vals <- lat_raw[keep]
    key <- paste(round(lng_vals, 6), round(lat_vals, 6))
    grps <- split(seq_along(key), key)
    # Cluster centers (original coords) for groups with n > 1
    cluster_centers <- list()
    for (idx in grps) {
      n <- length(idx)
      if (n <= 1) next
      cluster_centers[[length(cluster_centers) + 1]] <- list(
        lng = lng_vals[idx[1]],
        lat = lat_vals[idx[1]],
        idx = idx
      )
    }
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
    # Popup HTML for each cluster (Previous/Next to cycle through labels)
    cluster_lng <- numeric(0)
    cluster_lat <- numeric(0)
    cluster_popup <- list()
    prev_js <- paste0(
      "var p=this.closest('.cluster-popup');var n=parseInt(p.dataset.n,10);",
      "var cur=parseInt(p.dataset.cur||0,10);cur=(cur+n-1)%n;p.dataset.cur=cur;",
      "var items=p.querySelectorAll('.cluster-item');",
      "for(var i=0;i<items.length;i++)items[i].style.display=(i===cur)?'block':'none';",
      "p.querySelector('.cluster-num').textContent=cur+1;"
    )
    next_js <- paste0(
      "var p=this.closest('.cluster-popup');var n=parseInt(p.dataset.n,10);",
      "var cur=parseInt(p.dataset.cur||0,10);cur=(cur+1)%n;p.dataset.cur=cur;",
      "var items=p.querySelectorAll('.cluster-item');",
      "for(var i=0;i<items.length;i++)items[i].style.display=(i===cur)?'block':'none';",
      "p.querySelector('.cluster-num').textContent=cur+1;"
    )
    for (cc in cluster_centers) {
      idx <- cc$idx
      n <- length(idx)
      items_html <- character(n)
      for (j in seq_len(n)) {
        item_content <- as.character(labels[[idx[j]]])
        items_html[j] <- paste0(
          '<div class="cluster-item" style="display:',
          if (j == 1) "block" else "none", ';">', item_content, "</div>"
        )
      }
      popup_html <- paste0(
        '<div class="cluster-popup" data-n="', n, '" data-cur="0">',
        '<p>Sample <span class="cluster-num">1</span> of ',
        '<span class="cluster-total">', n, '</span></p>',
        '<div class="cluster-items">', paste(items_html, collapse = ""),
        '</div><br/><button type="button" onclick="', prev_js,
        '">Previous</button> <button type="button" onclick="', next_js,
        '">Next</button></div>'
      )
      cluster_lng <- c(cluster_lng, cc$lng)
      cluster_lat <- c(cluster_lat, cc$lat)
      cluster_popup <- c(cluster_popup, list(htmltools::HTML(popup_html)))
    }
    list(
      lng = lng_vals,
      lat = lat_vals,
      layerId = ids,
      label = labels,
      clusters = list(lng = cluster_lng, lat = cluster_lat, popup = cluster_popup)
    )
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

  selected_ancestry_call_ids <- reactive({
    sel <- input$table_input_rows_selected
    if (length(sel) < 1) return(integer(0))
    if (!"ancestryCallId" %in% names(ancestryCall)) return(integer(0))
    ids <- suppressWarnings(as.integer(ancestryCall$ancestryCallId[sel]))
    ids <- ids[is.finite(ids)]
    unique(ids)
  })

  plot_idx <- reactiveVal(1L)

  observeEvent(input$table_input_rows_selected, {
    ids <- selected_ancestry_call_ids()
    if (!length(ids)) {
      plot_idx(1L)
      return()
    }
    cur <- isolate(plot_idx())
    if (!is.finite(cur) || cur < 1L || cur > length(ids)) {
      plot_idx(1L)
    }
  }, ignoreNULL = FALSE)

  observeEvent(input$plot_prev, {
    ids <- selected_ancestry_call_ids()
    if (!length(ids)) return()
    cur <- isolate(plot_idx())
    if (!is.finite(cur)) cur <- 1L
    plot_idx(max(1L, cur - 1L))
  })

  observeEvent(input$plot_next, {
    ids <- selected_ancestry_call_ids()
    if (!length(ids)) return()
    cur <- isolate(plot_idx())
    if (!is.finite(cur)) cur <- 1L
    plot_idx(min(length(ids), cur + 1L))
  })

  current_ancestry_call_id <- reactive({
    ids <- selected_ancestry_call_ids()
    if (!length(ids)) return(NULL)
    idx <- plot_idx()
    if (!is.finite(idx)) idx <- 1L
    idx <- max(1L, min(as.integer(idx), length(ids)))
    ids[idx]
  })

  output$plot_nav_status <- renderText({
    ids <- selected_ancestry_call_ids()
    n <- length(ids)
    if (!n) return("Sample 0 of 0")
    idx <- plot_idx()
    if (!is.finite(idx)) idx <- 1L
    idx <- max(1L, min(as.integer(idx), n))
    paste("Sample", idx, "of", n)
  })

  output$pop_specific_accuracy_plot <- renderPlot({
    ac_id <- current_ancestry_call_id()
    if (is.null(ac_id)) {
      return(
        ggplot2::ggplot() +
          ggplot2::theme_void() +
          ggplot2::annotate(
            "text",
            x = 0.5,
            y = 0.5,
            label = "No ancestry call selected in the Table tab."
          )
      )
    }

    pop_acc <- getPopSpecificAccuracyForAncestryCall(ac_id)
    if (!nrow(pop_acc)) {
      return(
        ggplot2::ggplot() +
          ggplot2::theme_void() +
          ggplot2::annotate(
            "text",
            x = 0.5,
            y = 0.5,
            label = "No population-specific entries found."
          )
      )
    }

    req("populationDefinition" %in% names(pop_acc))
    req("accuracy" %in% names(pop_acc))

    pop_acc <- pop_acc %>%
      dplyr::mutate(
        series = "Population definition",
        populationDefinition = as.character(.data$populationDefinition),
        populationDefinition = dplyr::if_else(
          is.na(.data$populationDefinition) |
            !nzchar(trimws(.data$populationDefinition)),
          "Unknown population",
          .data$populationDefinition
        ),
        accuracy = suppressWarnings(as.numeric(as.character(.data$accuracy))),
        CILowerBound = suppressWarnings(
          as.numeric(as.character(.data$CILowerBound))
        ),
        CIUpperBound = suppressWarnings(
          as.numeric(as.character(.data$CIUpperBound))
        )
      )

    ac_overall <- getAncestryCallAccuracyForPlot(ac_id)
    if (nrow(ac_overall) > 0) {
      bio_sample <- as.character(ac_overall$bioSample[1])
      if (is.na(bio_sample) || !nzchar(trimws(bio_sample))) {
        bio_sample <- paste("ancestryCall", ac_id)
      }
      ac_row <- data.frame(
        series = "Selected BioSample",
        populationDefinition = paste("BioSample:", bio_sample),
        accuracy = suppressWarnings(
          as.numeric(as.character(ac_overall$accuracy[1]))
        ),
        CILowerBound = suppressWarnings(
          as.numeric(as.character(ac_overall$CILowerBound[1]))
        ),
        CIUpperBound = suppressWarnings(
          as.numeric(as.character(ac_overall$CIUpperBound[1]))
        ),
        stringsAsFactors = FALSE
      )
      pop_acc <- dplyr::bind_rows(pop_acc, ac_row)
    }

    pop_acc <- pop_acc %>%
      dplyr::mutate(
        CILowerBound = dplyr::coalesce(.data$CILowerBound, .data$accuracy),
        CIUpperBound = dplyr::coalesce(.data$CIUpperBound, .data$accuracy)
      )

    keep <- is.finite(pop_acc$accuracy)
    pop_acc <- pop_acc[keep, , drop = FALSE]

    if (!nrow(pop_acc)) {
      return(
        ggplot2::ggplot() +
          ggplot2::theme_void() +
          ggplot2::annotate(
            "text",
            x = 0.5,
            y = 0.5,
            label = "No numeric accuracy values available."
          )
      )
    }

    y_label <- "Accuracy (%)"
    pop_acc <- pop_acc %>%
      dplyr::mutate(
        accuracy = .data$accuracy * 100,
        CILowerBound = .data$CILowerBound * 100,
        CIUpperBound = .data$CIUpperBound * 100
      )

    pop_acc <- pop_acc %>%
      dplyr::arrange(.data$accuracy) %>%
      dplyr::mutate(
        populationDefinition = factor(
          .data$populationDefinition,
          levels = unique(.data$populationDefinition)
        )
      )

    ggplot2::ggplot(
      pop_acc,
      ggplot2::aes(
        x = .data$populationDefinition,
        y = .data$accuracy
      )
    ) +
      ggplot2::geom_errorbar(
        ggplot2::aes(
          ymin = .data$CILowerBound,
          ymax = .data$CIUpperBound,
          color = .data$series
        ),
        width = 0.22,
        linewidth = 1.0,
        alpha = 0.85,
        na.rm = TRUE
      ) +
      ggplot2::geom_point(
        ggplot2::aes(fill = .data$series),
        size = 3.2,
        shape = 21,
        color = "white",
        stroke = 0.7,
        na.rm = TRUE
      ) +
      ggplot2::scale_fill_manual(
        values = c(
          "Population definition" = "#2D7FB8",
          "Selected BioSample" = "#E87D2A"
        ),
        drop = FALSE
      ) +
      ggplot2::scale_color_manual(
        values = c(
          "Population definition" = "#1B1E23",
          "Selected BioSample" = "#9A4D12"
        ),
        drop = FALSE
      ) +
      ggplot2::labs(
        title = "Accuracy with confidence intervals",
        subtitle = paste("ancestryCallId:", ac_id),
        x = "Population definition",
        y = y_label,
        fill = NULL,
        color = NULL
      ) +
      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold"),
        axis.line.x = ggplot2::element_line(
          color = "#2D2D2D",
          linewidth = 0.6
        ),
        axis.line.y = ggplot2::element_line(
          color = "#2D2D2D",
          linewidth = 0.6
        ),
        axis.ticks = ggplot2::element_line(color = "#2D2D2D"),
        axis.text.x = ggplot2::element_text(
          size = 10,
          angle = 35,
          hjust = 1
        ),
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_blank()
      )
  })

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
    proxy <- leaflet::leafletProxy("map") %>% leaflet::clearMarkers()
    if (length(m$layerId) > 0) {
      proxy <- proxy %>%
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
    if (length(m$clusters$lng) > 0) {
      proxy <- proxy %>%
        leaflet::addCircleMarkers(
          lng = m$clusters$lng,
          lat = m$clusters$lat,
          layerId = paste0("cluster_", seq_along(m$clusters$lng)),
          radius = 10,
          fillOpacity = 0.9,
          stroke = TRUE,
          weight = 2,
          color = "#333333",
          fillColor = "#E67E22",
          popup = m$clusters$popup
        )
    }
  }, ignoreNULL = FALSE)

  # Filter Table to match clicked sites from map
  input_table_proxy = DT::dataTableProxy('table_input')
  observeEvent(input$map_marker_click, {
    id <- input$map_marker_click$id
    if (is.null(id) || grepl("^cluster_", id)) return()
    # If the click corresponds to an ancestryCall Experiment, filter Table by it
    exp_col <- if ("Experiment" %in% names(ancestryCall)) "Experiment" else
      if ("experiment" %in% names(ancestryCall)) "experiment" else NULL
    if (!is.null(exp_col) && id %in% ancestryCall[[exp_col]]) {
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
