#' Match trees based on NN distance
#'
#' Function to match detected tree tops with stem mapped trees. The 2D or 3D distance between each stem-mapped tree and detected tree is calculated and use to match the trees that are the closest to each others
#'
#' @param mapped_trees data.frame of mapped tree tops with at least the 5 following columns (from col 1 to 5 and in order): plot ID, unique tree ID, X, Y, Z. Other columns can be added after the 5th column.
#' @param detected_trees data.frame of detected tree tops with at least the 5 following columns (from col 1 to 5 and in order): plot ID, unique tree ID, X, Y, Z. Other columns can be added after the 5th column.
#' @param  distance_type Either "2D" or "3D".
#'
#' @export

match_nn_trees <- function(mapped_trees,
                           detected_trees,
                           distance_type = "3D") {

  # Check arguments
  if (length(distance_type) != 1) {
    stop("distance_type must be either \"2D\" or \"3D\"")
  }

  if (!distance_type %in% c("2D", "3D")) {
    stop("distance_type must be either \"2D\" or \"3D\"")
  }

  if(!is.data.frame(mapped_trees)) {
    stop("mapped_trees is not a data.frame")
  }

  if(ncol(mapped_trees) < 5) {
    stop("mapped_trees must have at least 5 columns in the following order: plot ID, tree ID, X, Y, Z")
  }

  if(!is.data.frame(detected_trees)) {
    stop("detected_trees is not a data.frame")
  }

  if(ncol(detected_trees) < 5) {
    stop("detected_trees must have at least 5 columns in the following order: plot ID, tree ID, X, Y, Z")
  }

  # Standardize column names
  colnames(mapped_trees)[1:5] <- c("plot_id", "tree_id", "X", "Y", "Z")
  colnames(detected_trees)[1:5] <- c("plot_id", "tree_id", "X", "Y", "Z")

  # Convert plot_id to character
  mapped_trees <- mapped_trees %>%
    mutate(across(c(plot_id, tree_id), as.character))

  detected_trees <- detected_trees %>%
    mutate(across(c(plot_id, tree_id), as.character))


  # Check for duplicated plot ID and tree ID
  find_duplicates <- duplicated(mapped_trees[, c("plot_id", "tree_id")])

  if(any(find_duplicates)) {
    stop(sprintf("Duplicated combination of plot_id and tree_id found in mapped_trees for plots %s and tree id %s", paste(mapped_trees$plot_id[find_duplicates], collapse = ","), paste(mapped_trees$tree_id[find_duplicates], collapse = ",")))
  }

  find_duplicates <- duplicated(detected_trees[, c("plot_id", "tree_id")])

  if(any(find_duplicates)) {
    stop(sprintf("Duplicated combination of plot_id and tree_id found in detected_trees for plots %s and tree id %s", paste(detected_trees$plot_id[find_duplicates], collapse = ","), paste(detected_trees$tree_id[find_duplicates], collapse = ",")))
  }

  # Keep only plots that are in both mapped and detected trees

  plot_id_mapped <- unique(mapped_trees$plot_id)
  plot_id_detected <- unique(detected_trees$plot_id)

  plot_id_inter <- intersect(plot_id_mapped, plot_id_detected)

  if(length(plot_id_inter) == 0) {
    stop("mapped_trees and detected trees do not share any plot ID")
  }

  mapped_trees <- mapped_trees %>%
    dplyr::filter(plot_id %in% plot_id_inter)

  detected_trees <- detected_trees %>%
    dplyr::filter(plot_id %in% plot_id_inter)

  #Initiate list that will store matching for all plots
  trees_matched <- list()

  # Loop through each plot
  for (n in 1:length(plot_id_inter)) {

    # Filter for current plot
    mapped_trees_plot <- dplyr::filter(mapped_trees,
                                       plot_id == plot_id_inter[n])

    detected_trees_plot <- dplyr::filter(detected_trees,
                                         plot_id == plot_id_inter[n])

    # Make ppp objects for spatstats functions
    mapped_trees_pp3 <- spatstat.geom::pp3(x = mapped_trees_plot$X,
                                           y = mapped_trees_plot$Y,
                                           z = mapped_trees_plot$Z,
                                           c(range(mapped_trees_plot$X),
                                             range(mapped_trees_plot$Y),
                                             range(mapped_trees_plot$Z)))

    mapped_trees_ppp <- spatstat.geom::ppp(x = mapped_trees_plot$X,
                                           y = mapped_trees_plot$Y,
                                           window = spatstat.geom::owin(range(mapped_trees_plot$X),
                                                         range(mapped_trees_plot$Y)))

    detected_trees_pp3 <- spatstat.geom::pp3(x = detected_trees_plot$X,
                                             y = detected_trees_plot$Y,
                                             z = detected_trees_plot$Z,
                                             c(range(detected_trees_plot$X),
                                               range(detected_trees_plot$Y),
                                               range(detected_trees_plot$Z)))

    detected_trees_ppp <- spatstat.geom::ppp(x = detected_trees_plot$X,
                                             y = detected_trees_plot$Y,
                                             window = spatstat.geom::owin(range(detected_trees_plot$X),
                                                           range(detected_trees_plot$Y)))

    # Cross distance between each detected tree and each stem mapped tree: 2D, 3D and vertical Zdist distance matrix

    dist_mtx <- list(dist_mtx_2D = spatstat.geom::crossdist(detected_trees_ppp, mapped_trees_ppp),
                     dist_mtx_3D = spatstat.geom::crossdist(detected_trees_pp3, mapped_trees_pp3),
                     Zdist_mtx = abs(outer(detected_trees_pp3$data$z, mapped_trees_pp3$data$z, "-")))

    # Add a dummy matrix (used internaly)
    dist_mtx$mtx <- matrix(1,
                           nrow = nrow(dist_mtx$dist_mtx_2D),
                           ncol = ncol(dist_mtx$dist_mtx_2D))

    # Set col and row names to tree IDs
    dist_mtx <- lapply(dist_mtx, function(x){
      colnames(x) <- mapped_trees_plot$tree_id
      rownames(x) <- detected_trees_plot$tree_id
      return(x)
    })

    # Iteratively look for NN

    trees_matched_list <- list() #Initate empty list
    counter_idx <- 1

    # As long that there is a detected tree to match and a stem-mapped tree available
    while(nrow(dist_mtx$mtx) > 0 & ncol(dist_mtx$mtx) > 0) {

      # Find smallest distance
      if(distance_type == "3D") {
        min_dist_idx <- which(dist_mtx$dist_mtx_3D == min(dist_mtx$dist_mtx_3D), arr.ind = TRUE)
      }else if(distance_type == "2D") {
        min_dist_idx <- which(dist_mtx$dist_mtx_2D == min(dist_mtx$dist_mtx_2D), arr.ind = TRUE)
      }

      # Add tree to tree_matched_list
      trees_matched_list[[counter_idx]] <- data.frame(tree_id_mapped = colnames(dist_mtx$mtx)[min_dist_idx[1,2]],
                                                      tree_id_detected = rownames(dist_mtx$mtx)[min_dist_idx[1,1]],
                                                      dist_2D = dist_mtx$dist_mtx_2D[min_dist_idx],
                                                      dist_3D = dist_mtx$dist_mtx_3D[min_dist_idx],
                                                      dist_Z = dist_mtx$Zdist_mtx[min_dist_idx])

      # Remove matched trees from dist matrices
      saved_dimnames <- dimnames(dist_mtx$mtx)

      if (nrow(dist_mtx$mtx) == 2) {
        dist_mtx <- lapply(dist_mtx, function(x) {
          out <- matrix(x[-min_dist_idx[1, 1], -min_dist_idx[1, 2]], nrow = 1)
          rownames(out) <- saved_dimnames[[1]][-min_dist_idx[1, 1]]
          colnames(out) <- saved_dimnames[[2]][-min_dist_idx[1, 2]]
          return(out)
        })
      }else if (ncol(dist_mtx$mtx) == 2) {
        dist_mtx <- lapply(dist_mtx, function(x) {
          out <- matrix(matrix(x[-min_dist_idx[1, 1], -min_dist_idx[1, 2]], ncol = 1))
          rownames(out) <- saved_dimnames[[1]][-min_dist_idx[1, 1]]
          colnames(out) <- saved_dimnames[[2]][-min_dist_idx[1, 2]]
          return(out)
        })
      }else{
        dist_mtx <- lapply(dist_mtx, function(x) {
          as.matrix(x[-min_dist_idx[1, 1], -min_dist_idx[1, 2]], rownames.force = TRUE)
        })
      }

      # Add 1 to counter for next iteration
      counter_idx <- counter_idx + 1

      #END OF WHILE LOOP

    }

    # Rbind all matched trees
    trees_matched_plot <- dplyr::bind_rows(trees_matched_list)

    # Check if some trees were not matched
    unmatched_mapped <- mapped_trees_plot$tree_id[which(!mapped_trees_plot$tree_id %in% trees_matched_plot$tree_id_mapped)]

    if (length(unmatched_mapped) > 0) {
      trees_matched_plot <- rbind(trees_matched_plot,
                                  data.frame(tree_id_mapped = unmatched_mapped,
                                             tree_id_detected = NA,
                                             dist_2D = NA,
                                             dist_3D = NA,
                                             dist_Z = NA))
    }

    unmatched_detected <- detected_trees_plot$tree_id[which(!detected_trees_plot$tree_id %in% trees_matched_plot$tree_id_detected)]

    if (length(unmatched_detected) > 0) {
      trees_matched_plot <- rbind(trees_matched_plot,
                                  data.frame(tree_id_mapped = NA,
                                             tree_id_detected = unmatched_detected,
                                             dist_2D = NA,
                                             dist_3D = NA,
                                             dist_Z = NA))
    }

    trees_matched[[n]] <- trees_matched_plot

    # END OF FOR LOOP (per plot)
  }

  # Return binded matched trees

  names(trees_matched) <- plot_id_inter

  out <- dplyr::bind_rows(trees_matched, .id = "plot_id")

  # Retrieve extra columns (if any)
  out <- dplyr::left_join(out, dplyr::select(mapped_trees, !c("X", "Y", "Z")), by = c("plot_id", "tree_id_mapped" = "tree_id"))

  out <- dplyr::left_join(out, dplyr::select(detected_trees, !c("X", "Y", "Z")), by = c("plot_id", "tree_id_detected" = "tree_id"))

  # Try to convert plot_id to integer. Otherwise keep as integer

  out <- tryCatch(
    {
      out_int <- out %>%
        mutate(plot_id = as.integer(plot_id))
      return(out_int)
    },
    warning = function(cond) {
      return(out)
    }
  )

  # Try to convert tree ids to integer. Otherwise keep as integer

  out <- tryCatch(
    {
      out_int <- out %>%
        mutate(across(contains(tree_id), as.integer))
      return(out_int)
    },
    warning = function(cond) {
      return(out)
    }
  )

  return(dplyr::as_tibble(out))
}
