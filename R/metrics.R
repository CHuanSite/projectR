#' @title Probability Distance Matrix
#'
#' @description Compute the probability distance matrix 
#'
#' @param x A matrix of dimension n * p, which is the embedded feature or a selection of genes
#' @param sigma_est Take value of NULL or a manually set number
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
    # x = x / apply(x, 2, sd)
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

#' @title KL Divergence between Two Probability Distance Matrix
#'
#' @description Compute the KL Divergence between two embedded probability distance matrix
#'
#' @param a The first embeded probability distance matrix
#' @param b The second embeded probability distance matrix
#'
#' @return KL Divergence between two probability distance matrix
#'
#' @keywords KL Divergence, Difference
#' 
#' @examples 
#' dist_a = probDistance(matrix(rnorm(100 * 3), nrow = 100) )
#' dist_b = probDistance(matrix(rnorm(100 * 5), nrow = 100) )
#' out = KL(dist_a, dist_b)
#' 
#' @export


KL <- function(dist_a, dist_b){
    x = dist_a * log(dist_a / dist_b)
    x[which(is.na(x))] = 0
    return(sum(x))
}


#' @title Measurement of Similarity
#'
#' @description Compute the similarity of two embedded latent spaces
#'
#' @param latent1 First embedded latent space
#' @param latent2 Second embedded latent space
#' @param sigma The sigma used in the maximum mean discrepancy (MMD) between two latent spaces
#'
#' @importFrom pdist pdist
#'
#' @return Computed MMD between two latent spaces
#'
#' @keywords Similarity, MMD 
#' 
#' @examples 
#' latent1 = matrix(rnorm(100 * 3), nrow = 100)
#' latent2 = matrix(rnorm(500 * 3), nrow = 500)
#' out = similarityMeasure(latent1, latent2, 1)
#'
#' @export

similarityMeasure <- function(latent1, latent2, sigma = 1){
    latent1 = (latent1 - apply(latent1, 2, mean)) / apply(latent1, 2, sd)
    latent2 = (latent2 - apply(latent2, 2, mean)) / apply(latent2, 2, sd)
    
    n1 = nrow(latent1)
    n2 = nrow(latent2)
    
    diff_1 = dist(latent1, diag = TRUE, upper = TRUE) %>% as.matrix
    diff_2 = dist(latent2, diag = TRUE, upper = TRUE) %>% as.matrix
    diff_12 = pdist(latent1, latent2) %>% as.matrix
    
    return(sum(exp(-1 * (diff_1 * diff_1) / sigma)) / (n1^2) + sum(exp(-1 * diff_2 * diff_2 / sigma)) / (n2^2) - 2 * sum(exp(-1 * diff_12 * diff_12 / sigma)) / (n1 * n2))
    
}
