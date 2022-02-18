#' @title DEET_feature_extract
#'
#' @description Identify which genes are associated with pieces of
#' metadata that a researcher queries.
#'
#' @param mat A gene-by-study matrix populated by the coefficients of that study.
#' By default, the coefficient is the log2Fold-change of genes as long as they are
#' differentially expressed (cutoff = padj < 0.05).
#' @param response A vector (binomial, categorical, or continuous) that is used to
#' associated the DEGs within the studies.
#' @param datatype indication of whether the response variable is binomial,
#' categorical, or continuous.
#'
#'
#' @return Named list given the elastic net coefficients and the eleastic net
#' regression between the response variable and the DEGs within DEET.
#' It also outputs the correlation, ANOVA, and wilcoxon test of every gene
#' against the response variable based on if it's continuous, categorical, or
#' binomial in nature.
#'
#' \itemize{
#'   \item elastic_net_coefficients - Association that a gene has with the
#'   response variable based on the elastic net regression.
#'   \item elastic_net - Output of the elastic net regression
#'   \item - basic_features gives the output of the
#'   correlation, ANOVA, and wilcoxon test of every gene against the
#'   response variable.
#' }
#'
#' @author Dustin Sokolowski, Jedid Ahn
#'
#' @examples
#'
#' data(DEET_feature_extract_example_matrix)
#' data(DEET_feature_extract_example_response)
#' single1 <- DEET_feature_extract(DEET_feature_extract_example_matrix,
#' DEET_feature_extract_example_response,"categorical")
#'
#' @references
#'
#' @export
#' @importFrom glmnet glmnet coef.glmnet
#' @importFrom stats cor.test p.adjust aov complete.cases wilcox.test
#' @importFrom utils data
#'
DEET_feature_extract <- function(mat, response, datatype) {

# Make sure that the type of variable where we'll extract features from is a type that we can work with
if(!(datatype %in% c("continuous","categorical","binomial"))) stop("Response variable must be 'continuous', 'categorical', or 'binomial'")



if(datatype == "continuous") {
  # Run an elastic net of the continuous variable
  if (is.factor(response)) response <- (as.numeric(levels(response))[response])

  tst1 <- glmnet::glmnet(x=t(mat), y=response, family = "gaussian", alpha=0.5)

  # run correlation between continuous response and mat variable.
  cor_spears <- t(apply(mat,1,function(x) {
    spear <- stats::cor.test(response, x)
    return(c(spear$estimate, spear$p.value))
  }
  ))
  colnames(cor_spears) <- c("Rho","P")
  cor_spears <- as.data.frame(cor_spears)
  cor_spears$FDR <- stats::p.adjust(cor_spears[,"P"])
  cor_spears <- cor_spears[cor_spears$FDR < 0.05,]
  basic <- list(res = cor_spears, test = "Spearman")


}
if(datatype == "categorical") {



  if(any(table(response) < 2)) {
    warning("Some of the response variables have fewer than two comparisons, removed")
  }
  response_names <- names(table(response))[table(response) > 2]
  mat <- mat[,(response %in% response_names)]
  response <- response[response %in% response_names]
  #mat <- mat[,(response %in% response_names)]
    l <- unique(response)

  if(length(l) == 2) stop("There are two conditions, datatype should be 'binomial'.")


  tst1 <- glmnet::glmnet(x=t(mat), y=response, family = "multinomial", alpha=0.5)


  aov_comp <- t(apply(mat,1,function(x) {
    sum1 <- summary(stats::aov(x ~ as.factor(response)))
    return(c(sum1[[1]]$`F value`[1], sum1[[1]]$`Pr(>F)`[1]))
  }))

  colnames(aov_comp) <- c("F", "P")
  aov_comp <- as.data.frame(aov_comp)
  aov_comp$FDR <- stats::p.adjust(aov_comp[,"P"])
  aov_comp <- aov_comp[order(aov_comp$FDR, decreasing = F),]
  aov_comp <- aov_comp[aov_comp$FDR < 0.05,]
  aov_comp <- aov_comp[stats::complete.cases(aov_comp),]
  basic <- list(res = aov_comp, test = "ANOVA")
}


if(datatype == "binomial") {
  tst1 <- glmnet::glmnet(x=t(mat), y=response, family = "binomial", alpha=0.5)
  # Could also add Wilcoxon's test
  l <- unique(y)
  if(length(l) != 2) stop("There are more than two conditions, datatype should be 'binomial' or 'continuous'.")

  wilcoxon_comp <- t(apply(mat,1,function(x) {
    in_one <- x[response == l[1]]
    other_one <- x[response == l[2]]
    wilcox <- stats::wilcox.test(in_one, other_one)$p.value
    if(mean(in_one) > mean(other_one)) {
      bias <- l[1]
    } else {
      bias <- l[2]
    }
    mean1 <- mean(in_one)
    mean2 <- mean(other_one)
    return(c(mean1,mean2,wilcox,bias))
  }))
  wilcoxon_comp <- as.data.frame(wilcoxon_comp, stringsAsFactors=FALSE)
  for(v in 1:3) wilcoxon_comp[,v] <- as.numeric(wilcoxon_comp[,v])
  colnames(wilcoxon_comp) <- c(paste0("mean_",l[1]), paste0("mean_",l[2]),"P","Bias" )
  wilcoxon_comp$FDR <- stats::p.adjust(wilcoxon_comp$P)
  wilcoxon_comp <- wilcoxon_comp[order(wilcoxon_comp$FDR), ]
  wilcoxon_comp <- wilcoxon_comp[wilcoxon_comp$FDR < 0.05,]
  wilcoxon_comp <- wilcoxon_comp[stats::complete.cases(wilcoxon_comp),]
  basic <- list(res = wilcoxon_comp, test = "Wilcoxon")
}


tst1_coef <- glmnet::coef.glmnet(tst1,s=0.1)

if(class( tst1_coef)== "list") { # This will be for binomial or gaussian methods
  getCoefs <- list()
  for(i in unique(response)) {
    RM <- apply(tst1_coef[[i]],1,mean)
    RM <- RM[abs(RM) > 0]
    getCoefs[[i]] <- RM
  }
} else { # this will be for continuous
  RM <- apply(tst1_coef,1,mean)
  RM <- RM[abs(RM) > 0]
  getCoefs <- RM
  # Genes associated with continuous variable in study

}

##
###
  outlist <- list(elastic_net_coefficients = getCoefs,
                  elastic_net = tst1,
                  basic_features = basic)
###
##

return(outlist)
}
