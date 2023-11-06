# Mesh generation ----
# ||||||||||||||||||||

## Utilities ----
## ||||||||||||||

# Returns the coordinates of the vertices of a SpatialPolygon

pol_coords <- function(pol){
  coords <- pol@polygons[[1]]@Polygons[[1]]@coords
  colnames(coords) <- c("x", "y")
  SpatialPoints(coords[-nrow(coords),])
  
  return(SpatialPoints(coords[-nrow(coords),]))
}

# Returns the distance between a point p and all the points contained in points
dist_point_from_points <- function(p, points){
  points <- data.frame(points)
  points$x <- points$x - p$x
  points$y <- points$y - p$y
  dist <- sqrt(points$x^2 + points$y^2)
  
  return(dist)
}

# Cleans a SpatialPolygons object keeping only the "main polygon" with
# larger area. It returns a SpatialPolygons object organized like that:
# - Main Polygon
# - Polygon hole 1
# - Polygon hole 2 
# - Polygon hole ...

pol_clean <- function(pol, remove_holes = FALSE, minimum_area_hole = NULL){
  
  max_area <- 0
  holes <- list()
  domain <- NULL
  
  for(i in 1:length(pol@polygons[[1]]@Polygons)){
    if(pol@polygons[[1]]@Polygons[[i]]@hole){
      if(!is.null(minimum_area_hole)){
        if(pol@polygons[[1]]@Polygons[[i]]@area > minimum_area_hole){
          holes <- c(holes, pol@polygons[[1]]@Polygons[[i]])
        } else {
          if(VERBOSE){
            cat("\nAn hole has been removed!\n")
          }
        }
      } else {
        holes <- c(holes, pol@polygons[[1]]@Polygons[[i]])
      }
    } else{
      if(pol@polygons[[1]]@Polygons[[i]]@area > max_area){
        max_area <- pol@polygons[[1]]@Polygons[[i]]@area
        domain <- pol@polygons[[1]]@Polygons[[i]]
      } else{
        if(VERBOSE){
          cat("\nSome outliers have been removed!\n")
        }
      }
    }
  }
  
  if(!remove_holes){
    domain_new <- SpatialPolygons(list(Polygons(c(domain, holes), ID = "domain")))
  } else {
    domain_new <- SpatialPolygons(list(Polygons(list(domain), ID = "domain")))
  }
  
  return(domain_new)
}

# Returns the number of holes in a SpatialPolygon
# !! It assumes that there exist only one main polygon and that all the
# !! others polygons are holes

pol_holes_number <- function(pol){
  length(pol@polygons[[1]]@Polygons)-1
}

# Returns the coordinates of the id^th hole
# !! It assumes that there exist only one main polygon and that all the
# !! others polygons are holes

pol_holes_coords <- function(pol, id){
  
  if(length(pol@polygons[[1]]@Polygons)-1 < id){
    cat(paste("Error: the number of holes is less than", id))
  }
  
  coords <- pol@polygons[[1]]@Polygons[[1+id]]@coords
  colnames(coords) <- c("x", "y")
  SpatialPoints(coords[-nrow(coords),])
  
  return(SpatialPoints(coords[-nrow(coords),]))
}


## Hexagons ----
## |||||||||||||

# Generates a hexagon given the central point and the radius

hex <- function(p, h) {
  
  alpha_vect <- seq(pi/6, 2*pi-pi/6, by = pi/3)
  
  hex <- Polygon(data.frame(x = p$x + h*cos(alpha_vect),
                            y = p$y + h*sin(alpha_vect)))
  hex <- SpatialPolygons(list(Polygons(list(hex), ID = "hexagon")))
  
  return(hex)  
}

# Generates a grid (either hexagonal or square) of step h on the rectangular 
# box identified by bbox

generate_grid <- function(bbox, h, seed_point, type = "hexagonal") {
  
  xmin <- bbox[1,1] - 2*h
  xmax <- bbox[1,2] + 2*h
  ymin <- bbox[2,1] - 2*h
  ymax <- bbox[2,2] + 2*h
  
  if(is.null(seed_point)){
    seed_point <- SpatialPoints(data.frame(x =  bbox[1,1], y = bbox[2,2]))
  }
  
  if(type == "square"){
    type = "regular"
  }
  
  # Create a polygon representing the rectangle
  bbox <- Polygon(cbind(c(xmin, xmax, xmax, xmin, xmin), c(ymin, ymin, ymax, ymax, ymin)))
  bbox_sp <- SpatialPolygons(list(Polygons(list(bbox), ID = "rectangle")))
  
  # Create a grid within the rectangle
  grid <- spsample(bbox_sp, type = type, cellsize = h)@coords
  colnames(grid) <- c("x", "y")
  grid <- SpatialPoints(grid)
  
  # Find the index of the point with the minimum distance from the seed
  distances <- dist_point_from_points(seed_point, grid)
  closest_point_index <- which.min(distances)
  closest_point <- grid[closest_point_index, ]
  translation <- data.frame(seed_point@coords - closest_point@coords)
  
  # Translation to match the seed
  grid <- data.frame(grid)
  grid$x <- grid$x + translation$x
  grid$y <- grid$y + translation$y
  grid <- SpatialPoints(grid)
  
  return(grid)
}

# Generates a lattice (either hexagonal or square) on the smallest rectangular
# box containing the locations points cloud.
# The polygons that do not contain any point are then discarded and the domain 
# is generated as the union of the remaining polygons.

generate_lattice <- function(locations, h, bbox, seed_point, type = "hexagonal") {
  
  if(is.null(bbox)){
    bbox <- locations@bbox
  }
  
  grid <- generate_grid(bbox, h, seed_point, type)
  
  check <- c()
  polygons <- list()
  for(i in 1:nrow(grid@coords)){
    point <- grid[i,]
    if(type == "hexagonal") polygon <- hex(point, h/sqrt(3) + 1e-9) else
      if(type == "square") polygon <- square(point, h/sqrt(2) + 1e-9)
      check[i] <- any(!is.na(over(locations, polygon)))
      if(check[i]){
        polygons <- c(polygons, polygon@polygons[[1]]@Polygons[[1]])
      }
  }
  
  # Domain
  lattice <- SpatialPolygons(list(Polygons(polygons, ID = "lattice")))
  domain <- as(st_union(st_as_sf(lattice), by_feature = TRUE, is_coverage = TRUE), "Spatial")
  
  return(list(grid = grid,
              domain = domain))
}


## Squares ----
## |||||||||||

# Generates a square grid of step h on the rectangular box 
# identified by bbox

square_grid <- function(bbox, h, seed_point = NULL) {
  
  xmin <- bbox[1,1]
  xmax <- bbox[1,2]
  ymin <- bbox[2,1]
  ymax <- bbox[2,2]
  
  x <- seq(xmin, xmax, by = h)
  y <- seq(ymin, ymax, by = h)
  grid <- expand.grid(x = x, y = y)
  
  return(grid)
}

# Generates a square given the central point and the length of the side

square <- function(p, h) {
  
  alpha_vect <- seq(pi/4, 2*pi-pi/4, by = pi/2)
  
  square <- Polygon(data.frame(x = p$x + h*cos(alpha_vect),
                               y = p$y + h*sin(alpha_vect)))
  square <- SpatialPolygons(list(Polygons(list(square), ID = "square")))
  
  return(square)  
}


## Domain ----
## |||||||||||

simplify_domain <- function(lattice, keep,
                            remove_holes = FALSE, minimum_area_hole = NULL,
                            keep_hole = 1) {
  
  domain <- lattice$domain
  grid <- lattice$grid
  
  # Remove outliers
  domain <- pol_clean(domain, remove_holes, minimum_area_hole)
  
  # Simplify outer boundary
  boundary <- SpatialPolygons(list(Polygons(list(domain@polygons[[1]]@Polygons[[1]]), ID = "boundary")))
  new_boundary <- ms_simplify(boundary, keep = keep)
  domain@polygons[[1]]@Polygons[[1]] <- new_boundary@polygons[[1]]@Polygons[[1]]
  
  # Boundary points
  points.b <- pol_coords(domain)
  
  # Holes
  points.h <- list()
  if(pol_holes_number(domain) >= 1){
    for(i in 1:pol_holes_number(domain)){
      
      # Simplify hole boundary
      boundary <- SpatialPolygons(list(Polygons(list(domain@polygons[[1]]@Polygons[[1+i]]), ID = "hole")))
      new_boundary <- ms_simplify(boundary, keep = keep_hole)
      domain@polygons[[1]]@Polygons[[1+i]] <- new_boundary@polygons[[1]]@Polygons[[1]]
      domain@polygons[[1]]@Polygons[[1+i]]@hole = TRUE
      
      # Holes boundary points
      points.h[[i]] <- pol_holes_coords(domain, i)
    }
  }
  
  return(list(domain = domain,
              points.b = points.b,
              points.h = points.h))
}


## Mesh ----
## |||||||||

# Generates a delaunay triangulation of the domain contained in lattice
# whose elements must satisfy the maximum area condition provided.
# It also handle the presence of holes inside the domain.

generate_mesh <- function(lattice, maximum_area) {
  
  domain <- lattice$domain        # Discretized domain \Omega_h
  # points.a <- lattice$points      # All the points
  points.b <- lattice$points.b    # Outer boundary points
  points.h <- lattice$points.h    # Holes boundary points
  points.g <- lattice$grid        # Grid points
  
  # Extract the coordinates of the vertices
  nodes.boundary <- points.b@coords
  segments <- cbind(1:nrow(nodes.boundary),
                    c(2:nrow(nodes.boundary), 1))
  holes <- data.frame(x = c(), y = c())
  if(length(points.h) >= 1){
    for(i in 1:length(points.h)){
      first <- nrow(nodes.boundary) + 1
      
      # Update nodes
      nodes.boundary <- rbind(nodes.boundary,
                              points.h[[i]]@coords)
      
      # Update segments
      segments <- rbind(segments,
                        cbind(first:nrow(nodes.boundary),
                              c((first+1):nrow(nodes.boundary), first)))
      
      # Update holes
      pol_hole <- Polygon(points.h[[i]]@coords)
      pol_hole <- SpatialPolygons(list(Polygons(list(pol_hole), ID = "hole")))
      inner_point <- as.numeric(spsample(pol_hole, type = "random", n = 1)@coords)
      holes <- rbind(holes, inner_point)
    }
  }
  
  # Final nodes
  nodes <- data.frame(rbind(nodes.boundary))#, points.g@coords, points.a@coords))
  rownames(nodes) <- 1:nrow(nodes)
  
  # Remove duplicates
  nodes <- as(unique(st_as_sf(SpatialPoints(nodes))), "Spatial")@coords
  
  # Mesh
  mesh <- create.mesh.2D(nodes, segments = segments, holes = holes)
  mesh <- refine.mesh.2D(mesh, 30, maximum_area)
  
  return(mesh)
}
