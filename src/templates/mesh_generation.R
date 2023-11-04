# Mesh generation ----
# ||||||||||||||||||||

set.seed(0)
locations <- SpatialPoints(locations)

## Domain ----
## |||||||||||

# Lattice
# h <- 10
# bbox <- NULL
# seed_point <- NULL
# type <- "hexagonal"
lattice <- generate_lattice(locations, h, bbox, seed_point, type)

# Domain simplification
# simplification <- 0.09
lattice_simplified <- simplify_domain(lattice, simplification)


## Discarded locations ----
## ||||||||||||||||||||||||

indexes.discarded_locations <- is.na(over(locations, lattice_simplified$domain))
names.locations <- names.locations[!is.na(over(locations, lattice_simplified$domain))]
locations.final <- locations.initial[names.locations,]


## Mesh ----
## |||||||||

# maximum_area <- 1000
set.seed(0)
mesh <- generate_mesh(lattice_simplified, maximum_area)


## Save mesh ----
## ||||||||||||||

lattice <- lattice_simplified

save(locations.final,
     mesh,
     lattice,
     # Saving options
     file = paste(directory.results, name.dataset, "_mesh", ".RData", sep = ""))


## Clean ----
## ||||||||||

rm(lattice, lattice_simplified, mesh,
   locations, locations.final,
   indexes.discarded_locations, names.locations,
   plot)



