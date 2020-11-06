Steps_PVA <-  function(patchMap,
                       threshold = 0,
                       trans_mat,
                       st_dev_mat,
                       stages = NULL,
                       Rmax_pop = NULL,
                       max_ind = 3,
                       k_scale = 0.1,
                       k_location = 0.6,
                       max_cells = NULL,
                       repeats = 1,
                       cores = 1,
                       directory,
                       init_pop,
                       scenario, classes) {
  
  filepath <- file.path(directory, scenario)
  if (!dir.exists(filepath)) {
    if (!dir.exists(directory)) dir.create(directory)
    dir.create(filepath)
  }
  mask <- patchMap
  mask[!is.na(patchMap)] <- 1
  
  if (is.null(stages)) stages <- paste0("Stage_", 1:ncol(trans_mat))
  
  colnames(trans_mat) <- rownames(trans_mat) <- stages
  colnames(st_dev_mat) <- rownames(st_dev_mat) <- stages
  
  # Determine max growth rate from stage matrix
  if (is.null(Rmax_pop)) Rmax_pop <- abs(eigen(trans_mat)$values[1])
  # Calculate the stable states of the population
  stable_states <- abs(eigen(trans_mat)$vectors[,1] / sum(eigen(trans_mat)$vectors[,1])) 
  # Get maximum number of individuals per cell
  popN <- stack(replicate(ncol(trans_mat), ceiling(mask[[1]]*max_ind*init_pop)))
  # Adjust life-stage populations based on stable state population
  popN <- popN * stable_states
  # Which cells are non-NA and below a threshold habitat suitability
  idx <- which(!is.na(patchMap[[1]][]) & patchMap[[1]][] < threshold)
  pop_str <- stack(
    foreach(i = 1:nlayers(popN)) %do% {
      # set.seed(123)
      m <- ceiling(cellStats(popN[[i]], max, na.rm = TRUE))
      pop <- popN[[i]]
      pop[] <- sapply(pop[], function(x) if(is.na(x)) x else rbinom(prob = (x/m),
                                                                    size = m,
                                                                    n = 1))
      pop[idx] <- 0 # Thin the initial population to cells with habitat suitability above threshold
      pop
    })
  names(pop_str) <- colnames(trans_mat)
  
  # DIAGNOSTIC: initial populations check
  TotpopN <- sum(cellStats(pop_str, 'sum', na.rm = TRUE)) # Get total population size to check if sensible
  pop_dist <- cellStats(pop_str, sum)/TotpopN # What proportions of individuals in each life-stage?
  
  # Define function to relate carrying capacity to habitat suitability
  k_function <- function (landscape, timestep) {
    
    suit <- landscape$suitability
    
    if (raster::nlayers(suit) > 1) {
      suit <- suit[[timestep]]
    }
    
    idx <- which(!is.na(suit[]))
    
    k <- suit
    k[idx] <- plogis(suit[idx], scale = k_scale, location = k_location) * max_ind
    
    k
  }
  
  # DIAGNOSTIC: Plot relationship of carrying capacity to habitat suitability
  #x <- seq(cellStats(patchMap[[1]], min, na.rm = TRUE), cellStats(patchMap[[1]], max, na.rm = TRUE), 0.1)
  #plot(x, plogis(x, scale = k_scale, location = k_location) * k_mult, type = 'l', xlab = "Habitat Quality", ylab = "Max Individuals")
  
  # Set up steps objects
  landscape <- landscape(population = pop_str,
                         suitability = patchMap,
                         carrying_capacity = k_function)
  
  if (is.null(max_cells)) max_cells <- rep(10, length(stages))
  
  pop_dynamics <- population_dynamics(
    change = growth(transition_matrix = trans_mat,
                    global_stochasticity = st_dev_mat,
                    transition_function = list(competition_density(stages = classes, # this regulates the population based on approach to CC for certain stages
                                                                   R_max = Rmax_pop)
                    )
    ),
    dispersal = #kernel_dispersal(max_distance = 30, # depends on resolution see below 
      # dispersal_kernel = exponential_dispersal_kernel(distance_decay = 30),#Distance where the mean proportional of individuals have dispersed
      # dispersal_proportion = all_dispersing(c(0,1,0)) # Liklihood of age class dispersing
      cellular_automata_dispersal(max_cells = max_cells, # Low dispersers this value should be lower e.g. SQG = 30, PO = 200
                                  dispersal_proportion = density_dependence_dispersing()
      ), # max_cells also depends on the resolution of the map. E.g. 100 m resolution with max_cells = 500 is a max beeline distance of 50 km
    modification = NULL,
    density_dependence = NULL
  )
  
  future::plan(multiprocess, workers = cores)
  
  #  system.time(
  all_pop_data <- foreach (i = 1:repeats, .combine = rbind) %do% {
    
    #profvis(      
    results <- simulation(landscape = landscape,
                          population_dynamics = pop_dynamics,
                          habitat_dynamics = NULL,
                          timesteps = 100,
                          replicates = cores,
                          # demo_stochasticity = "none", 
                          future.globals = list(trans_mat = trans_mat,
                                                st_dev_mat = st_dev_mat,
                                                Rmax_pop = Rmax_pop,
                                                k_location = k_location,
                                                k_scale = k_scale,
                                                max_ind = max_ind,
                                                max_cells = max_cells,
                                                cores = cores))
    #)
    
    pop_data <- return_pop_data(results)
    rm(results)
    gc()
    pop_data
  }
  # )
  
  save(all_pop_data, file = paste0(filepath, "/Pop_data_", scenario, ".RData"))
  
}


return_pop_data <- function(results_object){
  n_sims <- length(results_object)
  n_timesteps <- length(results_object[[1]])
  idx <- which(!is.na(raster::getValues(results_object[[1]][[1]]$population[[1]])))
  pop_mat <- matrix(NA, nrow = n_sims, ncol = n_timesteps)
  for (i in seq_len(n_sims)){
    for (j in seq_len(n_timesteps)) {
      pop_mat[i, j] <- sum(raster::extract(results_object[[i]][[j]]$population, idx))
    }
  }
  pop_mat
}

# Plot population trend

plot_pop <- function(pop_data, upper_limit, title = "") {
  max_pop <- max(pop_data)
  min_pop <- min(pop_data)
  mean_trend <- colSums(pop_data) / nrow(pop_data)
  matplot(t(pop_data),
          type = 'l',
          lwd = 0.5,
          lty = 1,
          ylim = c(0, upper_limit),
          col = 'grey',
          xlab = 'Timestep',
          ylab = 'Total Landscape Population')
  lines(1:100, mean_trend, lwd = 1.5)
  text(1, upper_limit * .98, labels = title, pos = 4, cex = 0.9)
}

# Plot expected minimum populations

plot_ema <- function (pop_data_list, interval = 95, all_points = FALSE) {
  
  # extract number of simulations
  n_species <- length(pop_data_list)
  species_names <- names(pop_data_list)
  n_sims <- length(pop_data_list[[1]])
  sim_names <- names(pop_data_list[[1]])
  
  # calculate interval range
  interval_range <- c((100 - interval) / 2, 100 - (100 - interval) / 2) / 100
  
  # populate table with ema mean and error values
  df <- foreach (i = seq_len(n_species), .combine = rbind) %do% {
    foreach (j = seq_len(n_sims), .combine = rbind) %do% {
      pops <- pop_data_list[[i]][[j]]
      min_total_pops <- apply(pops, 1, function(x) min(x))
      ema_mean <- mean(min_total_pops)
      ema_lower <- stats::quantile(min_total_pops, interval_range)[1]
      ema_upper <- stats::quantile(min_total_pops, interval_range)[2]
      data.frame("name" = as.factor(paste0(species_names[[i]])),
                 "sim" = as.factor(paste0(sim_names[[j]])),
                 "ema_mean" = ema_mean,
                 "ema_lower" = ema_lower,
                 "ema_upper" = ema_upper)
    }
  }
  
  max_limit <- ceiling(max(df$ema_upper))
  
  # get the raw values also if all_points is set to TRUE
  if (all_points) {
    df_raw <- foreach (i = seq_len(n_species), .combine = rbind) %do% {
      foreach (j = seq_len(n_sims), .combine = rbind) %do% {
        pops <- pop_data_list[[i]][[j]]
        min_total_pops <- apply(pops, 1, function(x) min(x))
        data.frame("name" = as.factor(paste0(species_names[[i]])),
                   "sim" = as.factor(paste0(sim_names[[j]])),
                   "min_pops" = min_total_pops)
      }
    }
    max_limit <- ceiling(max(df_raw[, -c(1:2)]))
  }
  
  ggplot(df, aes(x = sim, y = ema_mean, ymin = ema_lower, ymax = ema_upper)) +
    {if (all_points)geom_jitter(data = df_raw, mapping = aes(y = min_pops), position = position_jitter(0.2), col = "lightgrey", size = 0.1)} +
    geom_errorbar(width = .2, position = position_dodge(0.05)) +
    geom_point(mapping = aes(shape = sim, fill = sim), color = "black", size = 2.0) +
    geom_hline(yintercept = 0, size = 1.2) +
    facet_grid(.~name,
               switch = "x", # Moves the labels from the top to the bottom
               #labeller = label_both # Adds the labels to the year and X variables
    ) + 
    xlab("") +
    ylab(paste0("Expected Minimum Abundance (", interval, "% Interval)")) +
    scale_y_continuous(limits = c(0, max_limit), expand = c(0, 0)) +
    scale_fill_manual(values = c('black', 'white', 'black', 'white')) +
    scale_shape_manual(values = c(21, 21, 24, 24)) +
    theme(
      strip.background = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white"),
      panel.border = element_blank(),
      axis.line.y = element_line(color="black"),
      legend.title = element_blank(),
      legend.key = element_blank()
    )
}
