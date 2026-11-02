#' Launch stGradient Interactive Polygon Selector
#'
#' Opens a Shiny app for interactively drawing polygons on spatial data
#'
#' @return Launches Shiny application
#' @export
#' @examples
#' \dontrun{
#' launch_stGradient()
#' }
launch_stGradient <- function() {
  
  ui <- shiny::fluidPage(
    shiny::titlePanel("stGradient: Tumor Nest Selector"),
    
    shiny::fluidRow(
      shiny::column(12, 
                    shiny::textInput("seurat_name", 
                                     "Seurat object name:", 
                                     value = "seurat_obj"),
                    shiny::textInput("fov_name", 
                                     "FOV name:", 
                                     value = "fov"),
                    shiny::actionButton("load_data", 
                                        "Load Data", 
                                        style = "margin-bottom: 20px;")
      )
    ),
    
    plotly::plotlyOutput("plot", height = "600px"),
    
    shiny::fluidRow(
      shiny::column(6, 
                    shiny::textInput("polygon_name", 
                                     "Polygon object name:", 
                                     value = "drawn_polygons")),
      shiny::column(6, 
                    shiny::actionButton("extract", 
                                        "Extract Polygons", 
                                        style = "margin-top: 25px;"))
    ),
    
    shiny::verbatimTextOutput("shapes")
  )
  
  server <- function(input, output, session) {
    
    rv <- shiny::reactiveValues(
      shapes = NULL,
      coords = NULL,
      data_loaded = FALSE
    )
    
    shiny::observeEvent(input$load_data, {
      tryCatch({
        seurat_obj <- get(input$seurat_name, envir = .GlobalEnv)
        coords <- get_cell_coords(seurat_obj, input$fov_name)
        coords$is_cancer <- seurat_obj@meta.data[coords$cell_id, "is_cancer"]
        
        rv$coords <- coords
        rv$shapes <- NULL
        rv$data_loaded <- TRUE
        
        output$shapes <- shiny::renderPrint({
          cat(paste0("✓ Loaded data from '", input$seurat_name, 
                     "' with FOV '", input$fov_name, "'\n"))
          cat(paste0("  ", nrow(coords), " cells loaded\n"))
        })
        
      }, error = function(e) {
        output$shapes <- shiny::renderPrint({
          cat("ERROR loading data:\n")
          cat(paste0("  ", e$message, "\n"))
          cat("\nMake sure the Seurat object name and FOV name are correct.\n")
        })
        rv$data_loaded <- FALSE
      })
    })
    
    output$plot <- plotly::renderPlotly({
      shiny::req(rv$data_loaded)
      
      coords <- rv$coords
      
      plotly::plot_ly(data = coords, 
                      x = ~x, 
                      y = ~y,
                      type = 'scatter', 
                      mode = 'markers',
                      color = ~is_cancer,
                      colors = c("TRUE" = "red", "FALSE" = "lightgray"),
                      marker = list(size = 3, opacity = 0.5),
                      source = "plot") %>%
        plotly::layout(
          dragmode = 'drawclosedpath',
          newshape = list(
            line = list(color = 'blue', width = 2),
            fillcolor = 'rgba(0,0,255,0.2)'
          ),
          shapes = list(),
          title = "Draw tumor regions"
        ) %>%
        plotly::config(
          modeBarButtonsToAdd = c('drawclosedpath', 'drawrect', 'drawcircle', 'eraseshape'),
          displayModeBar = TRUE,
          editable = TRUE
        )
    })
    
    shiny::observe({
      event_data <- plotly::event_data("plotly_relayout", source = "plot")
      if (!is.null(event_data)) {
        rv$shapes <- apply_relayout_shapes(shiny::isolate(rv$shapes), event_data)
      }
    })
    
    shiny::observeEvent(input$extract, {
      
      if (!rv$data_loaded) {
        output$shapes <- shiny::renderPrint({
          cat("ERROR: Please load data first by clicking 'Load Data'.\n")
        })
        return()
      }
      
      shapes_data <- rv$shapes
      obj_name <- input$polygon_name
      
      if (obj_name == "" || !grepl("^[a-zA-Z][a-zA-Z0-9._]*$", obj_name)) {
        output$shapes <- shiny::renderPrint({
          cat("ERROR: Invalid object name.\n")
          cat("Object name must start with a letter and contain only letters, numbers, dots, or underscores.\n")
        })
        return()
      }
      
      output$shapes <- shiny::renderPrint({
        cat("Attempting to extract shapes...\n\n")
        
        if (is.null(shapes_data)) {
          cat("No shapes data captured yet.\n")
          return()
        }
        
        if (length(shapes_data) > 0) {
          cat(paste("Found", length(shapes_data), "shape(s)\n\n"))
          
          polygons_list <- list()
          
          for (i in seq_along(shapes_data)) {
            cat(paste("Processing shape", i, "\n"))
            
            coords_parsed <- shape_to_coords(shapes_data[[i]])
            
            if (!is.null(coords_parsed)) {
              polygon_sf <- sf::st_polygon(list(as.matrix(coords_parsed))) %>%
                sf::st_sfc() %>%
                sf::st_sf(region_id = i)
              
              polygons_list[[i]] <- polygon_sf
              cat(paste("  ✓ Shape", i, "parsed successfully\n"))
            } else {
              cat(paste("  ✗ Failed to parse shape", i, "\n"))
            }
          }
          
          if (length(polygons_list) > 0) {
            all_polygons <- do.call(rbind, polygons_list)
            
            assign(obj_name, all_polygons, envir = .GlobalEnv)
            
            cat(paste0("\n✓ SUCCESS! Polygons created as '", obj_name, 
                       "' in your environment\n\n"))
            print(all_polygons)
          } else {
            cat("\nFailed to parse any polygons from the shape data.\n")
          }
        } else {
          cat("No shapes found. Draw at least one shape before extracting.\n")
        }
      })
    })
  }
  
  shiny::shinyApp(ui, server)
}

# Merge a plotly_relayout event into the current list of shapes.
# Drawing/erasing sends the full `shapes` array; editing sends keys like
# "shapes[0].path" or "shapes[1].x0"; zoom/pan sends only axis keys.
apply_relayout_shapes <- function(current, event_data) {
  if (is.null(current)) current <- list()
  keys <- names(event_data)
  
  if ("shapes" %in% keys) {
    shapes <- event_data$shapes
    if (is.data.frame(shapes)) {
      current <- lapply(seq_len(nrow(shapes)), function(i) {
        as.list(shapes[i, intersect(c("type", "path", "x0", "y0", "x1", "y1"),
                                    names(shapes)), drop = FALSE])
      })
    } else if (length(shapes) == 0) {
      current <- list()
    } else {
      current <- lapply(shapes, as.list)
    }
  }
  
  edit_pattern <- "^shapes\\[([0-9]+)\\]\\.(path|x0|y0|x1|y1)$"
  for (key in grep(edit_pattern, keys, value = TRUE)) {
    idx <- as.integer(sub(edit_pattern, "\\1", key)) + 1
    field <- sub(edit_pattern, "\\2", key)
    if (idx <= length(current)) {
      current[[idx]][[field]] <- event_data[[key]]
    }
  }
  
  current
}

# Convert a plotly shape (path, rect or circle) to a closed coordinate matrix.
shape_to_coords <- function(shape, n_circle_points = 64) {
  type <- if (is.null(shape$type)) "path" else shape$type
  
  if (type == "path") {
    if (is.null(shape$path) || is.na(shape$path)) return(NULL)
    return(parse_svg_path(shape$path))
  }
  
  bounds <- suppressWarnings(as.numeric(c(shape$x0, shape$y0, shape$x1, shape$y1)))
  if (length(bounds) != 4 || any(is.na(bounds))) return(NULL)
  x0 <- bounds[1]; y0 <- bounds[2]; x1 <- bounds[3]; y1 <- bounds[4]
  
  if (type == "rect") {
    coords <- cbind(x = c(x0, x1, x1, x0, x0), y = c(y0, y0, y1, y1, y0))
  } else if (type == "circle") {
    theta <- seq(0, 2 * pi, length.out = n_circle_points + 1)
    coords <- cbind(x = (x0 + x1) / 2 + abs(x1 - x0) / 2 * cos(theta),
                    y = (y0 + y1) / 2 + abs(y1 - y0) / 2 * sin(theta))
    coords[nrow(coords), ] <- coords[1, ]
  } else {
    return(NULL)
  }
  
  coords
}
