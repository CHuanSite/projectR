#' Probability Distance Matrix
#'
#' Compute the probability distance matrix 
#'
#' @param x A matrix of dimension n * p, which is the embedded feature or a selection of genes
#' @param sigma_est Take value of NULL or a manually set number
#'
#' @importFrom pdist dist
#'
#' @return A computed probability distance matrix
#'
#' @keywords probability distance matrix
#'
#' @examples 
#' x = matrix(rnorm(100 * 3), nrow = 100) 
#' dist_mat = probDistance(x)
#' 
#' @export

probDistance <- function(x, sigma_est = NULL){
    x = x - apply(x, 2, mean)
    x = x / apply(x, 2, sd)
    if(is.null(sigma_est)){
        sigma = mean(apply(x, 2, sd))
    }else{
        sigma = sigma_est
    }
    dist_matrix = dist(x, diag = TRUE, upper = TRUE) %>% as.matrix
    exp_dist_matrix = exp(-1 * (dist_matrix * dist_matrix) / (sigma * sigma))
    diag(exp_dist_matrix) = 0
    exp_dist_sum = apply(exp_dist_matrix, 1, sum)
    exp_dist = exp_dist_matrix / exp_dist_sum
    return(exp_dist)
}
