#' Buffer function for rasters
#'
#' Buffr calculates a buffer around any user-defined cell value for RasterLayer objects.
#'
#' @name buffr
#'
#' @usage buffr(r, distance, units = "geographic", target_value = 1,
#'   mask = TRUE, mask_value = NA, max_rows_at_once = 1e+08,
#'   verbose = TRUE)
#'
#' @param r RasterLayer.
#' @param distance numeric>0. The buffer distance. The units are defined in the units argument.
#' @param units character. This argument must equal either "cell" or "geographic". If set to "cell", then the buffer distance is interpreted as the number of cells. If set to "geographic", then the distance is interpreted as being in the same units as the resolution of the raster file.
#' @param target_value numeric. Cell value to perform the buffer operation on. Buffered cells will take on this value.
#' @param mask logical. If TRUE then all cells that do not contain the target value are masked with the value defined in the mask_value argument (NA by default). If FALSE then all cells that do not contain the target value are returned with their original value.
#' @param mask_value numeric. Value to use for masking when the mask argument is set to TRUE.
#' @param max_rows_at_once numeric>0. Sets the maximum number of rows for R to process at once, which can cause larger jobs to be broken into smaller loops. Recommended to leave at the default value of 100000000.
#' @param verbose logical. Should progress indicators be supressed?
#'
#' @return RasterLayer
#' @author Ariel Fridman
#'
#' @details
#' Buffr assigns a cell value based on whether the center of the cell falls within the buffer radius.
#'
#' Buffr offers several advantages over the buffer function in the raster package v2.6-7:
#' \enumerate{
#'   \item It is substantially faster, especially for large rasters and buffers. This is partially because the computationally intensive operations are implemented in C++ using the Rcpp package.
#'   \item The user can define which cell value to buffer.
#'   \item It includes progress indicators.
#' }
#'
#' @examples
#' library(raster)
#' r <- raster(ncol = 36, nrow = 18)
#' r[] <- NA
#' r[500] <- 1
#' buff <- buffr(r = r, distance = 20, units = "geographic", target_value = 1)
#' plot(buff)
#'
#' r[100] <- 2
#' buff <- buffr(r = r, distance = 20, units = "geographic", target_value = 2, mask = FALSE)
#' plot(buff)
#'
#' @useDynLib buffr
#' @importFrom Rcpp sourceCpp
NULL
#' @export

buffr <- function (r, distance, units = "geographic", target_value = 1, mask = TRUE, mask_value = NA, max_rows_at_once = 100000000, verbose = TRUE) {

  if (units=="cell") {
    xres <- 1
    yres <- 1
  } else if (units=="geographic") {
    xres <- raster::xres(r)
    yres <- raster::yres(r)
  } else {
    if (verbose==TRUE) {warning("Units not correctly specified. Must be either \"cell\" or \"geographic\".\n")}
    return(r)
  }

  width_cell_x <- ceiling(distance/xres)
  width_cell_y <- ceiling(distance/yres)

  a <- r[]
  a <- matrix(data = a, nrow = dim(r)[1], ncol = dim(r)[2], byrow = TRUE)

  #Calculate edge
  edge <- which(as.array(a==target_value), arr.ind = TRUE)

  if (mask==TRUE) {
    a <- matrix(data = mask_value, nrow = dim(r)[1], ncol = dim(r)[2])
  }

  if (nrow(edge)==0) {
    if (verbose==TRUE) {warning("No cells with buffer value found.\n")}
    return(r)
  } else if (nrow(edge)==1) {
    edge_buffer_all_points <- edge
    edge_buffer_edge_points <- data.frame()
    rm(edge)
  } else if (nrow(edge)>1) {
    row_col <- paste_int(row = edge[,1], col = edge[,2])
    edge <- cbind(edge, edge = ifelse(pinp_num(paste_int(row = edge[,1]+1, col = edge[,2]+1), row_col) &
                                        pinp_num(paste_int(row = edge[,1], col = edge[,2]+1), row_col) &
                                        pinp_num(paste_int(row = edge[,1]-1, col = edge[,2]+1), row_col) &
                                        pinp_num(paste_int(row = edge[,1]+1, col = edge[,2]), row_col) &
                                        pinp_num(paste_int(row = edge[,1]-1, col = edge[,2]), row_col) &
                                        pinp_num(paste_int(row = edge[,1]+1, col = edge[,2]-1), row_col) &
                                        pinp_num(paste_int(row = edge[,1], col = edge[,2]-1), row_col) &
                                        pinp_num(paste_int(row = edge[,1]-1, col = edge[,2]-1), row_col), 0, 1))
    edge <- matrix(edge[edge[,3]==1,1:2], ncol = 2, dimnames = list(NULL, c("row", "col")))

    #Calculate edge of edge points
    row_col <- paste_int(row = edge[,1], col = edge[,2])

    buffer_all_points_col <- !(pinp_num(row_col, paste_int(row = edge[,1], col = edge[,2]+1)))
    row_col_rev <- row_col[buffer_all_points_col]
    edge_buffer_all_points <- matrix(edge[buffer_all_points_col,], ncol = 2, dimnames = list(NULL, c("row", "col")))

    buffer_all_points_row <- !(pinp_num(row_col_rev, paste_int(row = edge_buffer_all_points[,1]+1, col = edge_buffer_all_points[,2])))
    row_col_rev <- row_col_rev[buffer_all_points_row]
    edge_buffer_all_points <- matrix(edge_buffer_all_points[buffer_all_points_row,], ncol = 2, dimnames = list(NULL, c("row", "col")))

    edge_buffer_edge_points <- matrix(edge[!(pinp_num(row_col, row_col_rev)),], ncol = 2, dimnames = list(NULL, c("row", "col")))

    rm(row_col, row_col_rev, edge, buffer_all_points_col, buffer_all_points_row)
  }

  #buffer distance

  #all buffer points
  buffer_around <- expand.grid.alt(-width_cell_x:width_cell_x, -width_cell_y:width_cell_y)
  buffer_around <- cbind(buffer_around, distance = sqrt((buffer_around[,1]*xres)^2+(buffer_around[,2]*yres)^2))
  buffer_around <- buffer_around[buffer_around[,3]<=distance, 1:2]

  #edge buffer points
  buffer_around_shifted_coord <- cbind(buffer_around[,1]+abs(min(buffer_around[,1]))+1, buffer_around[,2]+abs(min(buffer_around[,2]))+1)
  row_col <- paste_int(row = buffer_around_shifted_coord[,1], col = buffer_around_shifted_coord[,2])
  buffer_around_edge <- cbind(buffer_around_shifted_coord, edge = ifelse(pinp_num(paste_int(row = buffer_around_shifted_coord[,1]+1, col = buffer_around_shifted_coord[,2]+1), row_col) &
                                                                           pinp_num(paste_int(row = buffer_around_shifted_coord[,1], col = buffer_around_shifted_coord[,2]+1), row_col) &
                                                                           pinp_num(paste_int(row = buffer_around_shifted_coord[,1]-1, col = buffer_around_shifted_coord[,2]+1), row_col) &
                                                                           pinp_num(paste_int(row = buffer_around_shifted_coord[,1]+1, col = buffer_around_shifted_coord[,2]), row_col) &
                                                                           pinp_num(paste_int(row = buffer_around_shifted_coord[,1]-1, col = buffer_around_shifted_coord[,2]), row_col) &
                                                                           pinp_num(paste_int(row = buffer_around_shifted_coord[,1]+1, col = buffer_around_shifted_coord[,2]-1), row_col) &
                                                                           pinp_num(paste_int(row = buffer_around_shifted_coord[,1], col = buffer_around_shifted_coord[,2]-1), row_col) &
                                                                           pinp_num(paste_int(row = buffer_around_shifted_coord[,1]-1, col = buffer_around_shifted_coord[,2]-1), row_col), 0, 1))
  buffer_around_edge <- buffer_around_edge[buffer_around_edge[,3]==1,1:2]
  buffer_around_edge <- cbind(Var1 = buffer_around_edge[,1]-(abs(min(buffer_around[,1]))+1), Var2 = buffer_around_edge[,2]-(abs(min(buffer_around[,2]))+1))

  #RLE for all buffer points
  buffer_around_shifted_mat <- matrix(0, nrow = max(buffer_around_shifted_coord[,1]), ncol = max(buffer_around_shifted_coord[,2]))
  buffer_around_shifted_mat[buffer_around_shifted_coord] <- 1
  run_encoded <- rle(as.integer(buffer_around_shifted_mat))
  w <- which(run_encoded$values==1)
  run_encoded$lengths[w+1] <- run_encoded$lengths[w+1]+run_encoded$lengths[w]-1
  run_encoded$values[w] <- run_encoded$lengths[w]
  run_encoded$lengths[w] <- 1
  buffer_around_shifted_mat <- matrix(inverse.rle(run_encoded), nrow = max(buffer_around_shifted_coord[,1]), ncol = max(buffer_around_shifted_coord[,2]))
  rle_coord <- which(as.array(buffer_around_shifted_mat>0), arr.ind = TRUE)
  buffer_around_rle <- cbind(rle_coord, values = buffer_around_shifted_mat[rle_coord])
  buffer_around_rle[,1] <- buffer_around_rle[,1]-(abs(min(buffer_around[,1]))+1)
  buffer_around_rle[,2] <- buffer_around_rle[,2]-(abs(min(buffer_around[,2]))+1)

  rm(buffer_around_shifted_coord, buffer_around_shifted_mat, run_encoded, w, rle_coord, row_col)

  if (verbose==TRUE) {cat("Starting Part 1 of 3\n")}

  for_rle <- matrix(0, nrow = nrow(a)+width_cell_x, ncol = ncol(a))

  max_rows_at_once <- ifelse(nrow(buffer_around_edge)>max_rows_at_once, nrow(buffer_around_edge), max_rows_at_once)
  edge_size <- floor(max_rows_at_once/nrow(buffer_around_edge))

  for (i in 1:ceiling(nrow(edge_buffer_all_points)/edge_size)) {

    if (nrow(edge_buffer_all_points)==1) {
      edge_temp <- edge_buffer_all_points
    } else if (nrow(edge_buffer_all_points)>1) {
      edge_temp <- edge_buffer_all_points[(i*edge_size-edge_size+1):min(nrow(edge_buffer_all_points), edge_size*i),]
    }

    gc()

    added_buffer <- expand_gridC_2(vec1_1 = edge_temp[,1], vec1_2 = edge_temp[,2], vec2_1 = buffer_around_rle[,1], vec2_2 = buffer_around_rle[,2])

    added_buffer <- cbind(rowSumsC(added_buffer),
                          value = rep.int(buffer_around_rle[,3], rep.int(nrow(edge_temp),nrow(buffer_around_rle))))

    added_buffer <- added_buffer[f2(mat = added_buffer, max_col = ncol(a)),]
    added_buffer[,3] <- f3(mat = added_buffer)
    added_buffer[,1] <- f4(vec = added_buffer[,1])

    added_buffer <- added_buffer[rev(order(added_buffer[,3])),]
    added_buffer <- added_buffer[!(dupe(paste_int(added_buffer[,1], added_buffer[,2]))),]

    for_rle[added_buffer[,1:2]] <- added_buffer[,3]

    rm(edge_temp, added_buffer)

    if (verbose==TRUE) {cat(i, "of", ceiling(nrow(edge_buffer_all_points)/edge_size), "iterations complete\n")}
  }

  if (verbose==TRUE) {cat("Starting Part 2 of 3\n")}

  run_encoded <- rle(as.integer(for_rle))
  run_encoded_save <- run_encoded

  first_one <- ifelse(run_encoded$values[1]==1,1,2)
  running_count <- 0
  running_count_original_i <- 0

  for (i in 1:length(run_encoded$lengths)) {
    val <- run_encoded$values[i]
    len <- run_encoded$lengths[i]
    if (val>0) { #value is positive
      if (val+(len-1)<running_count & len<running_count) { #value is positive - less than or equal to running count
        running_count <- running_count-len
        run_encoded$lengths[i] <- 0
      } else if (val+(len-1)<running_count & len>=running_count) {
        run_encoded$lengths[i] <- len-running_count
        running_count <- 0
      } else { #value is positive - greater than running count
        running_count <- val+len-1
        run_encoded$lengths[i] <- running_count
        run_encoded$values[i] <- 1
        run_encoded$lengths[running_count_original_i] <- ifelse(sum(run_encoded_save$lengths[running_count_original_i:(i-1)])>run_encoded$lengths[running_count_original_i],run_encoded$lengths[running_count_original_i], sum(run_encoded_save$lengths[(running_count_original_i):(i-1)]))
        running_count_original_i <- i
        running_count <- running_count-len
      }
    } else { #value is 0
      if (len<=running_count) { #value is 0 - less than or equal to running count
        running_count <- running_count-len
        run_encoded$lengths[i] <- 0
      } else { #value is 0 - greater than running count
        run_encoded$lengths[i] <- len-running_count
        running_count <- 0
      }
    }
  }

  for_rle <- matrix(inverse.rle(run_encoded), nrow = nrow(a)+width_cell_x, ncol = ncol(a))
  for_rle <- for_rle[1:nrow(a),]

  a[which(as.array(for_rle==1), arr.ind = TRUE)] <- target_value

  rm(for_rle, run_encoded, run_encoded_save, i)

  if (verbose==TRUE) {cat("Starting Part 3 of 3\n")}

  if (nrow(edge_buffer_edge_points)>0) {
    for (i in 1:ceiling(nrow(edge_buffer_edge_points)/edge_size)) {

      edge_temp <- matrix(edge_buffer_edge_points[(i*edge_size-edge_size+1):min(nrow(edge_buffer_edge_points), edge_size*i),], ncol = 2, dimnames = list(NULL, c("row", "col")))

      gc()

      added_buffer <- expand_gridC_2(vec1_1 = edge_temp[,1], vec1_2 = edge_temp[,2], vec2_1 = buffer_around_edge[,1], vec2_2 = buffer_around_edge[,2])

      added_buffer <- rowSumsC(added_buffer)
      added_buffer <- added_buffer[f1(mat = added_buffer, max_row = nrow(a)),]
      a[added_buffer] <- target_value

      rm(edge_temp, added_buffer)

      if (verbose==TRUE) {cat(i, "of", ceiling(nrow(edge_buffer_edge_points)/edge_size), "iterations complete\n")}
    }
  }

  b <- as.integer(t(a))
  r[] <- b
  return(r)
}
