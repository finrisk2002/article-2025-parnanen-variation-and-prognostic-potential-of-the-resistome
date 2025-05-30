myplot <- function (A, feat, tse, base.size=20) {
  tse$Level <- A[feat, ]

  # Swap PC2
  # reducedDim(tse, "PCoA_ARG")[,2] <- -reducedDim(tse, "PCoA_ARG")[,2]

  scater::plotReducedDim(tse, "PCoA_ARG",
                            point_size=0.8, point_alpha=0.7,
                            colour_by = "Level") +			    
  labs(x="PC1", y="PC2", title=feat) +
  scale_color_gradient(low="white", high="black", name="Weight") +
  coord_equal() +
  theme(plot.margin  = unit(c(1, 1, 1, 1), "pt"),
	plot.title = element_text(size = base.size),  
        legend.position=c(0.8, 0.8),
        legend.text = element_text(size=0.7*base.size),
	legend.title = element_text(size=base.size),
                                  axis.text.x = element_text(size = base.size),
                                  axis.text.y = element_text(size = base.size),
                                  axis.title.x = element_text(size = base.size),
                                  axis.title.y = element_text(size = base.size)	  
				  ) 
}


#' @title Radial Theta Function
#' @description Adapted from \pkg{NeatMap} and \pkg{phyloseq} packages
#'   but not exported and hence not available via phyloseq. Completely
#'   rewritten to avoid license conflicts. Vectorized to gain efficiency;
#'   only calculates theta and omits r.
#' @param x position parameter
#' @return theta
#' @keywords internal
radial_theta <- function(x) {
    
    x <- as(x, "matrix")
    theta <- atan2((x[, 2] - mean(x[, 2])), (x[, 1] - mean(x[, 1])))
    names(theta) <- rownames(x)
    
    theta
    
}


get.dfplot2 <- function (fit, inner.quantile=0.05, outer.quantile=0.025) {

  library(tidybayes)
    # gather_draws(b_ARGload_log10,b_Family,b_Age,b_MEN,b_CURR_SMOKE,b_BL_USE_RX_L,b_BL_USE_RX_J01,b_PREVAL_DIAB,b_BMI,b_BP_TREAT,b_SysBP,b_TULOT,b_KY100_14) %>%
    # summarise_draws() %>%

  dfplot <- t(apply(as_draws_array(fit, variable = "^b_", regex = TRUE), 3, function (x) keyvars2(x, iq=inner.quantile, oq=outer.quantile))) %>%
    as.data.frame() %>%
    arrange(median) %>%
    select(-mean) %>%    
    # Convert to the original domain 
    mutate(median=exp(median)) %>%
    mutate(outer.quantile.low=exp(outer.quantile.low)) %>%
    mutate(outer.quantile.high=exp(outer.quantile.high)) %>%
    mutate(inner.quantile.low=exp(inner.quantile.low)) %>%
    mutate(inner.quantile.high=exp(inner.quantile.high))     

  library(reshape2)
  dfplot$variable <- rownames(dfplot)
  dfplot$variable <- str_remove(dfplot$variable, "b_")
  
  dfplot <- dfplot %>% mutate(variable_print = stringr::str_replace_all(variable, substitutions))
  dfplot <- dfplot %>% mutate(variable_print=factor(variable_print, levels=dfplot$variable_print)) 
  dfplot$variable <- NULL

  dfplot
  
}


keyvars2 <- function (x, iq=0.1, oq=0.04) {

  # Extract mean, inner and outer quantiles 
  c(
    outer.quantile.low=quantile(x, oq/2, names=FALSE),
    inner.quantile.low=quantile(x, iq/2, names=FALSE), 
    mean=mean(x, na.rm=TRUE, names=FALSE),
    median=median(x, na.rm=TRUE, names=FALSE),    
    # mean=log(mean(exp(x), names=FALSE)), # Mean in the original domain
    inner.quantile.high=quantile(x, 1-iq/2, names=FALSE),
    outer.quantile.high=quantile(x, 1-oq/2, names=FALSE)
  )    
}

get.dfplot <- function (fit) {

  library(tidybayes)
  
  # get_variables(fit)
  # paste( get_variables(fit)[2:14], collapse=",")

  dfplot <- fit %>% gather_draws(b_ARGload_log10,b_Enterobacteriaceae,b_Age,b_MEN,b_CURR_SMOKE,b_BL_USE_RX_L,b_BL_USE_RX_J01,b_PREVAL_DIAB,b_BMI,b_BP_TREAT,b_SysBP,b_TULOT,b_KY100_14) %>%
    summarise_draws() %>%
    arrange(median) %>%
    select(-mean) %>%
    # Convert to the original domain 
    mutate(median=exp(median)) %>%
    mutate(q5=exp(q5)) %>% # this is lower 95% i.e. quantile 0.025
    mutate(q95=exp(q95))   # this is upper 95% i.e. quantile 0.975

  # Sort by median
  dfplot$variable <- str_remove(dfplot$.variable, "b_")
  dfplot <- dfplot %>% mutate(variable_print = stringr::str_replace_all(variable, substitutions))
  dfplot <- dfplot %>% mutate(variable_print=factor(variable_print, levels=dfplot$variable_print)) 
  dfplot$.variable <- NULL

  dfplot
  
}

get.dfplot3 <- function (fit) {

  library(tidybayes)
  
  # get_variables(fit)
  # paste( get_variables(fit)[2:14], collapse=",")

  dfplot <- fit %>% gather_draws(b_ARGload_log10,b_Enterobacteriaceae,b_Age,b_MEN,b_CURR_SMOKE,b_BL_USE_RX_L,b_BL_USE_RX_J01,b_PREVAL_DIAB,b_BMI,b_BP_TREAT,b_SysBP,b_TULOT,b_food) %>%
    summarise_draws() %>%
    arrange(median) %>%
    select(-mean) %>%
    # Convert to the original domain 
    mutate(median=exp(median)) %>%
    mutate(q5=exp(q5)) %>%
    mutate(q95=exp(q95))

  # Sort by median
  dfplot$variable <- str_remove(dfplot$.variable, "b_")
  dfplot <- dfplot %>% mutate(variable_print = stringr::str_replace_all(variable, substitutions))
  dfplot <- dfplot %>% mutate(variable_print=factor(variable_print, levels=dfplot$variable_print)) 
  dfplot$.variable <- NULL

  dfplot
  
}





get.fitted.model <- function (d, prob=0.95) {

  # Cox model with brms
  library(brms)
  set.seed(6322)
  fit <- NULL
  fit <- brm(TIMEDIFF | cens(1 - EVENT) ~ 1 + .,
           data = d,
	   family = brmsfamily("cox"))
  # summary(fit) # cox_model_sepsisx
  fit
}



get.covariates <- function () {

  # Q57 Exercise at free-time
  # ((1 [I don't move much]) (2 [walk, cycle, etc. At least 4 h a week]) (3 [excercise at least 3 h a week]) (4 [regularly excercise  for competitive sports]))

  mortality.variables <- c("Enterobacteriaceae", # Enterobact.
                           "Age",            # age
			   "MEN",            # gender
			   "CURR_SMOKE",     # smoking
			   "BL_USE_RX_L",    # use of antineoplastic and immunomodulating agents
			   "BL_USE_RX_J01",  # antibiotics use in the past 6 months
			   "PREVAL_DIAB",    # prevalent diabetes
			   "BMI",            # bmi
			   "BP_TREAT", 	     # self-reported antihypertensive medication		   
			   "SysBP")          # blood pressure
  lifestyle.variables <- c("TULOT",          # income
                           "KY100_14"        #,  # salad (KY100_14) or healthy food score (HFC_score)
			   # "Q57"           # exercise, weak measure
			   )            

   

   varnames <- c(mortality.variables, lifestyle.variables)
   varnames

}


prepare.data.for.survival.analysis <- function (TSE, family="Enterobacteriaceae") {

  # TODO define some of this in Carpentry_TreeSE.R directly

  # Prepare data frame for the survival analysis
  DF <- cbind(as.data.frame(colData(TSE)), as.data.frame(t(assay(altExp(TSE, "Family"), "relabundance"))))

  # NOTE:
  # Add pseudocount (min non zero value / 2; take log10; then scale; not scaling relabundances directly)
  DF$Family <- log10(DF[, family] + min(DF[, family][DF[, family]>0])/2)

  # It is more reasonable to represent for interpretable effect sizes:
  DF$Age <- DF$BL_AGE/10  # age per 10 years
  DF$SysBP <- DF$SYSTM/10 # systolic BP per 10 hhmg
  DF$V2 <- DF$V2/1e6 # Library size - transform to Million reads

  # Variables that are more appropriate in log10 domain
  DF$ARGload_log10 <- log10(DF$SUM_norm)
  DF$vaesto  <- log10(DF$vaesto)
  DF[[family]] <- DF[["Family"]]; DF[["Family"]] <- NULL
  
  DF

}




age.labels <- function () {

  age.labs <- c("24-\n29",
             "30-\n39",
             "40-\n49",
             "50-\n59",
             "60-\n69",
             "70-\n75")

  age.labs

}

# Plot
brplot <- function (post, base.size=25) {
  dw <- 0.2
  ggplot(post, aes(x=varname, color=sex)) +
    geom_point(aes(y=mean), size=3,
      position=position_dodge(width=dw)
    ) +
    geom_line(aes(x=varname, y=mean, group=sex)) +     
    geom_linerange(aes(ymin=outer.quantile.low, ymax=outer.quantile.high), size=1,
      position=position_dodge(width=dw)) +
    scale_y_continuous(limits=range(unlist(lapply(posteriors, function (x) {x[, 2:6]})))) +    
    scale_color_manual(
      values = c("chocolate4", "cornflowerblue"),
      name  = "",
      labels = c("Women", "Men")
    ) +  
    labs(y="ARG load (RPKM)", x="Level") +
    theme_set(theme_classic(base.size)) +
      theme(
        plot.title = element_text(size = base.size),
                                  axis.text.x = element_text(size = 0.9*base.size),
                                  axis.text.y = element_text(size = base.size),
                                  axis.title.x = element_text(size = base.size),
                                  axis.title.y = element_text(size = base.size),
	                          legend.text = element_text(size=0.9*base.size)
				  ) 			 
}

# Plot
brplot2 <- function (post, base.size=25) {
  dw <- 0.2
  ggplot(post, aes(x=varname, color=Sex)) +
    geom_point(aes(x=varname, y=mean, color=Sex), size=3,
      position=position_dodge(width=dw)
    ) +
    geom_line(aes(x=varname, y=mean, color=Sex, group=Sex)) +     
    geom_linerange(aes(ymin=outer.quantile.low, ymax=outer.quantile.high, color=Sex), size=1,
      position=position_dodge(width=dw)) +
    #scale_y_continuous(limits=range(unlist(lapply(posteriors, function (x) {x[, 2:6]})))) +
    scale_y_continuous(limits=range(unlist(post[, 2:6]))) +        
    scale_color_manual(
      values = rev(c("chocolate4", "cornflowerblue")),
      name = ""#,
      # labels = c("Women", "Men")
    ) +
    facet_grid(. ~ Region) +
    # facet_wrap(Region ~ .) +     
    labs(y="ARG load (RPKM)", x="Level") +
    theme_set(theme_classic(base.size)) +
      theme(
        plot.title = element_text(size = base.size),
                                  axis.text.x = element_text(size = 0.9*base.size),
                                  axis.text.y = element_text(size = base.size),
                                  axis.title.x = element_text(size = base.size),
                                  axis.title.y = element_text(size = base.size),
	                          legend.text = element_text(size=0.9*base.size)
				  ) 			 
}



brplot3 <- function (post, base.size=25) {
  dw <- 0.2
  ggplot(post, aes(x=varname, color=Sex)) +
    geom_point(aes(y=mean), size=3,
      position=position_dodge(width=dw)
    ) +
    geom_linerange(aes(ymin=outer.quantile.low, ymax=outer.quantile.high), size=1,
      position=position_dodge(width=dw)) +
    scale_y_continuous(limits=range(unlist(lapply(posteriors, function (x) {x[!x$varname=="b_Intercept", 2:6]})))) +    
    scale_color_manual(
      values = c("chocolate4", "cornflowerblue"),
      name  = "",
      labels = c("Women", "Men")
    ) +  
    labs(y="ARG load (RPKM)", x="Level") +
    theme_set(theme_classic(base.size)) +
      theme(
        plot.title = element_text(size = base.size),
                                  axis.text.x = element_text(size = 0.9*base.size),
                                  axis.text.y = element_text(size = base.size),
                                  axis.title.x = element_text(size = base.size),
                                  axis.title.y = element_text(size = base.size),
	                          legend.text = element_text(size=0.9*base.size)
				  ) +
    geom_hline(yintercept=0, linetype=2) +     				  
    coord_flip() +
    facet_grid(Region ~ .)

				  
}

get.post <- function (df2, inner.quantile=0.1, outer.quantile=0.025) {

  # Control same covariates than in the linear model
  m0 <- lapply(unique(df2$group), function (g) {brm(log(SUM_norm) ~ varname + PREVAL_RX_J01_NEVT + PREVAL_RX_J01A_NEVT + PREVAL_RX_J01A + PREVAL_RX_J01F + BL_USE_RX_J01,
                                                         family = gaussian(), data = (df2 %>% filter(group==g)))})

  names(m0) <- unique(df2$group)
    
  library(tidyverse)
  # Extract posterior draws and summaries for all variables
  post <- do.call("rbind", lapply(names(m0), function (nam) {
    #exp(t(apply(as_draws_array(m0[[nam]], variable = "^b_", regex = TRUE), 3, function (x) keyvars(x, iq=inner.quantile, oq=outer.quantile)))) %>% # lognormal
    t(apply(as_draws_array(m0[[nam]], variable = "^b_", regex = TRUE), 3, function (x) keyvars(x, iq=inner.quantile, oq=outer.quantile))) %>%       # gaussian
    as.data.frame() %>%
    tibble::rownames_to_column(var="varname") %>%
    mutate(varname=str_replace(varname, "b_factorvarname", "")) %>%
    mutate(group=rep(nam, n()))}
    )) %>%
    separate("group", c("Sex", "Region")) %>%
    mutate(Sex=factor(Sex, levels=c("Women","Men"))) %>%
    mutate(Region=factor(Region))      
    # mutate(Region=factor(Region, levels=c("Helsinki","Turku","Oulu","Kuopio","Karelia","Lapland")))  

  post$variable <- rep(myvar, nrow(post))

  res <- list(m0=m0, post=post)

}


brplot4 <- function (df, base.size=20) {
  dw <- 0.2
  df[, 2:6] <- exp(df[, 2:6])
  ggplot(df, aes(x=Region, color=Sex)) +
    geom_point(aes(y=mean), size=3,
      position=position_dodge(width=dw)
    ) +
    geom_linerange(aes(ymin=outer.quantile.low, ymax=outer.quantile.high), size=1,
      position=position_dodge(width=dw)) +
    scale_y_log10(limits=range(df[, 2:6]),
                  #breaks=seq(min(round(range(df[, 2:6]), 2)), max(round(range(df[, 2:6]), 2)), 0.03)) +
                  #breaks=seq(0.94, 1.16, 0.04)) +
                  breaks=seq(0.95, 1.15, 0.05)) +    		  		  
    scale_color_manual(
      values = c("chocolate4", "cornflowerblue"),
      name  = ""
    ) +  
    labs(y="Effect size", x="") +    
    theme_set(theme_classic(base.size)) +
      theme(
        plot.title = element_text(size = base.size),
                                  # title.text = element_text(size = base.size),
                                  axis.text.x = element_text(size = base.size),				  
                                  axis.text.y = element_text(size = base.size),
                                  axis.title.x = element_text(size = base.size),
                                  axis.title.y = element_text(size = base.size),
	                          legend.text = element_text(size=base.size)
				  ) +
    geom_hline(yintercept=1, linetype=2) +     				  
    coord_flip() +
    # facet_grid(Region ~ variable)
    #facet_wrap(. ~ variable, nrow=2) +
    facet_wrap(variable ~ ., ncol=1) +    
    theme(legend.position = "top")

}

# -----------------



dic.varnames <- function () {
  c(KY100_14="Fresh vegetables",
  KY100_21="Poultry consumption",
  TULOT="Household income",
  vaesto_quantile="Population density",
  vaesto="Population density",  
  age_class="Age",
  BL_AGE="Age"  
  )
}

get.arg.palette <- function () {
  require(RColorBrewer)
  arg.palette <- brewer.pal(name="Set1", n=9)[-c(3,6)] # remove green and yellow
  names(arg.palette) <-  c(
         "Tetracycline",
         "Beta-lactam",
	 "Macrolide, Aminoglycoside, Tetracycline, Quinolone, Amphenicol, Rifamycin",
	 "Macrolide, Lincosamide, Streptogramin B",	 
	 "Amphenicol",	 
	 "Aminoglycoside",
	 "Other"
	 )
  arg.palette["Other"] <- "gray"

  arg.palette
}


get.arg.gene.palette <- function () {
  # remotes::install_github("KarstensLab/microshades")
  #pal <- RColorBrewer::brewer.pal(6, "Spectral") # Feel free to change :-)
  arg.palette <- get.arg.palette()    
  pal <- c("tet(O)"="darkred",
         "tet(Q)"="red",
	 "tet(W)"="pink",
	 "cfxA6"=unname(arg.palette["Beta-lactam"]),
	 "erm(B)"=unname(arg.palette["Macrolide, Lincosamide, Streptogramin B"]),
	 "Other"=unname(arg.palette["Other"]))
  pal
}





calculateEnrichment.brms <- function (tse, varname, assay.type, features=NULL, nbreaks=5) {

  A0 <- assay(tse, assay.type=assay.type)
  
  if (is.null(features)) {features <- rownames(A0)}  
  g <- NULL
  gs <- levels(colData(tse)[, varname])
 
  dat <- NULL
  pics <- list()
  for (feat in features) {
    vals <- NULL
    A <- A0[feat,]
    # Split in bins with min value (0 or very small usually) and then the rest on equally sized 4 quartiles
    br <- unique(c(min(A), quantile(A[A>min(A)], probs=seq(0, 1, length=nbreaks))))
    cuts <- cut(A, breaks=br, include.lowest=TRUE, right=FALSE)
    for (g in gs){
      # Observed / Expected gives the enrichment
      ptab <- table(cuts, colData(tse)[, varname])[, tolower(g)] / (mean(colData(tse)[, varname]==tolower(g)) * table(cuts))
      vals <- rbind(vals, unname(ptab))
    }

    # Bayes estimate for high ARG mode enrichment
    ddd <- data.frame(Bin=cuts, Mode=(colData(tse)[, varname] == "high")+0)
    set.seed(5452)
    library(brms)
    fit <- brm(formula = Mode ~ Bin - 1, data = ddd, family = bernoulli(link="logit")) 
    ddd2 <- as.data.frame(exp(fixef(fit)))[, c(1,3,4)] 
    ddd2$Bin <- factor(paste0("B", 0:4))
    # Proportional enrichments could be added .. 
    #ddd2$Enrichment <- ddd2$Estimate/0.1 - 1
    #ddd2$Enrichment <- ddd2$Estimate/0.1 - 1     
    pics[[feat]] <- ggplot(ddd2, aes(x=Bin, y=Estimate)) + 
      geom_errorbar(aes(ymin=Q2.5, ymax=Q97.5), width=.1) +
      geom_line() +
      geom_point() +
      geom_line(aes(x=Bin, y=Estimate, group=1)) +
      geom_hline(aes(yintercept=0.1), linetype=2, color="gray") +
      labs(x="", y="High ARG load prevalence (%)", title=feat) +
      scale_y_continuous(label=scales::percent)

  }

  pics

}

calculateEnrichment <- function (tse, varname, assay.type, features=NULL, nbreaks=5) {

  A0 <- assay(tse, assay.type=assay.type)
  
  if (is.null(features)) {features <- rownames(A0)}  
  g <- NULL
  gs <- levels(colData(tse)[, varname])
 
  dat <- NULL 
  for (feat in features) {
    vals <- NULL
    A <- A0[feat,]
    # Split in bins with min value (0 or very small usually) and then the rest on equally sized 4 quartiles
    br <- unique(c(min(A), quantile(A[A>min(A)], probs=seq(0, 1, length=nbreaks))))
    cuts <- cut(A, breaks=br, include.lowest=TRUE, right=FALSE)
    for (g in gs){
      # Observed / Expected gives the enrichment
      ptab <- table(cuts, colData(tse)[, varname])[, tolower(g)] / (mean(colData(tse)[, varname]==tolower(g)) * table(cuts))
      vals <- rbind(vals, unname(ptab))
    }
    colnames(vals) <- paste0("B", 0:4)
    vals <- as.data.frame(vals)
    vals$Level <- rep(feat, length(gs))
    vals$Mode <- gs
    dat <- rbind(dat, vals)

  }

  df <- reshape2::melt(dat %>% filter(Mode=="high") %>% dplyr::select(-Mode), id.vars=c("Level"))
  df$value <- df$value-1
  df$Level <- factor(df$Level, levels=features)
  df$Bin <- factor(df$variable); df$variable <- NULL

  df

}




calculateEnrichment2 <- function (tse, varname, assay.type, features=NULL, nbreaks=5) {

  A0 <- assay(tse, assay.type=assay.type)
  
  if (is.null(features)) {features <- rownames(A0)}  
  g <- NULL
  gs <- levels(colData(tse)[, varname])
 
  dat <- NULL 
  for (feat in features) {
    vals <- NULL
    A <- A0[feat,]
    # Split in bins with min value (0 or very small usually) and then the rest on equally sized 4 quartiles
    br <- unique(c(min(A), quantile(A[A>min(A)], probs=seq(0, 1, length=nbreaks))))
    cuts <- cut(A, breaks=br, include.lowest=TRUE, right=FALSE)
    dat <- cbind(dat, cuts)
  }

  colnames(dat) <- features
  rownames(dat) <- colnames(tse)
  
  dat

}


calculateEnrichment3 <- function (tse, varname, assay.type, features=NULL, nbreaks=5) {

  A0 <- assay(tse, assay.type=assay.type)
  
  if (is.null(features)) {features <- rownames(A0)}  
  g <- NULL
  gs <- levels(colData(tse)[, varname])
 
  dat <- NULL
  libsizes <- NULL # median lib sizes
  for (feat in features) {
    vals <- NULL
    A <- A0[feat,]
    # Split in bins with min value (0 or very small usually) and then the rest on equally sized 4 quartiles
    br <- unique(c(min(A), quantile(A[A>min(A)], probs=seq(0, 1, length=nbreaks))))
    cuts <- cut(A, breaks=br, include.lowest=TRUE, right=FALSE)
    for (g in gs){
      # Observed / Expected gives the enrichment
      ptab <- table(cuts, colData(tse)[, varname])[, tolower(g)] / (mean(colData(tse)[, varname]==tolower(g)) * table(cuts))
      vals <- rbind(vals, unname(ptab))
    }
    colnames(vals) <- paste0("B", 0:4)
    vals <- as.data.frame(vals)
    vals$Level <- rep(feat, length(gs))
    vals$Mode <- gs
    dat <- rbind(dat, vals)

    # Median lib sizes for the given cuts
    libsizes <- rbind(libsizes, sapply(split(tse$V2,cuts), median))

  }

  rownames(libsizes) <- features
  colnames(libsizes) <- colnames(vals)[1:5]

  df <- reshape2::melt(dat %>% filter(Mode=="high") %>% select(-Mode), id.vars=c("Level"))
  df$value <- df$value-1
  df$Level <- factor(df$Level, levels=features)
  df$Bin <- factor(df$variable); df$variable <- NULL


  #df
  libsizes %>% as.data.frame

}




#### Connecting regions and municipality geometries ####

#### Connecting regions and municipality geometries ####
# Make a function to avoid repetition
# NOTE: OBJECT mun MUST BE IN GLOBAL ENVIRONMENT BEFORE RUNNING THESE FUNCTIONS
clean_and_join_regions <- function(region_name = "Lapland", index_number = 1, color="HighARGRatio", color2="HighARGPrevalence") {

  inds <- which(mun$maakunta_name_en %in% region_name)
  
  # note: deep assignment arrow
  # print(index_number)
  mun$freq.relative[inds] <<- dfr[[color]][index_number]

  # note: deep assignment arrow
  mun$freq[inds] <<- dfr[[color2]][index_number]
  # mun$mean.rpkm[inds] <- dfr$mean.rpkm[2]
  vaesto_prob <- mun$vaesto[which(mun$maakunta_name_en %in% region_name)] / sum(mun$vaesto[which(mun$maakunta_name_en %in% region_name)])
  n_cases <- sample(x = inds, size = dfr$N[index_number], replace = TRUE, prob = vaesto_prob)
  xyz <- table(n_cases)
  xyz_names <- names(xyz)
  for (j in seq_along(xyz)) {
    index <- as.numeric(xyz_names[j])
    # note: deep assignment arrow
    mun$share_of_cases[index] <<- xyz[j]
  }
}

clean_and_join_municipalities <- function(municipality_name = "Helsinki", index_number = 1, color="HighARGRatio", color2="HighARGPrevalence") {
  # if (!exists(object)) stop("Give an existing object")
  # if (!(class(object) %in% "sf")) stop("Give an sf object")
  inds <- which(mun$municipality_name_fi %in% municipality_name)
  # note: deep assignment arrow
  mun$freq.relative[inds] <<- dfr[[color]][index_number]
  # note: deep assignment arrow
  mun$freq[inds] <<- dfr[[color2]][index_number]
  vaesto_prob <- mun$vaesto[which(mun$municipality_name_fi %in% municipality_name)] / sum(mun$vaesto[which(mun$municipality_name_fi %in% municipality_name)])
  n_cases <- sample(x = inds, size = dfr$N[index_number], replace = TRUE, prob = vaesto_prob)
  xyz <- table(n_cases)
  xyz_names <- names(xyz)
  for (j in seq_along(xyz)) {
    index <- as.numeric(xyz_names[j])
    # note: deep assignment arrow
    mun$share_of_cases[index] <<- xyz[j]
  }
}


keyvars <- function (x, iq=0.1, oq=0.04) {

  # Extract mean, inner and outer quantiles 
  c(
    outer.quantile.low=quantile(x, oq/2, names=FALSE),
    inner.quantile.low=quantile(x, iq/2, names=FALSE), 
    mean=mean(x, na.rm=TRUE, names=FALSE),
    # mean=log(mean(exp(x), names=FALSE)), # Mean in the original domain
    inner.quantile.high=quantile(x, 1-iq/2, names=FALSE),
    outer.quantile.high=quantile(x, 1-oq/2, names=FALSE)
  )    
}

# Inverse logit
toprob <- function (logodds) {

  # Odds
  o <- exp(logodds)

  # Probs
  p <- o/(1+o)

  return(p)
  
}




pal.family <- function () { 

  # Colorblind friendly
  palo <- palette.colors(palette = "Okabe-Ito")

  #black        orange       skyblue   bluishgreen        yellow 
  #    "#000000"     "#E69F00"     "#56B4E9"     "#009E73"     "#F0E442" 
  #         blue    vermillion reddishpurple          gray 
  #    "#0072B2"     "#D55E00"     "#CC79A7"     "#999999" 

  c(
    #"Bacteroidaceae"="#663333",
    #"Bifidobacteriaceae"="#6699CC",
    #"Enterobacteriaceae"="khaki1",    
    #"Lachnospiraceae"="#330000",
    #"Prevotellaceae"="#CC9999",
    #"Ruminococcaceae"="lightcyan4",
    #"Coriobacteriaceae"="#888888"
    #"Other"="#CC9966",
    #"Eubacteriaceae"="#CCCFFF",
    #"Rikenellaceae"="khaki4",
    
    "Bacteroidaceae"=unname(palo["reddishpurple"]),    
    "Bifidobacteriaceae"=unname(palo["skyblue"]),    
    "Enterobacteriaceae"=unname(palo["black"]),    
    "Eubacteriaceae"=unname(palo["vermillion"]),
    "Lachnospiraceae"=unname(palo["orange"]),    
    "Prevotellaceae"=unname(palo["bluishgreen"]),        
    #"Rikenellaceae"=unname(palo["yellow"]),
    "Ruminococcaceae"=unname(palo["blue"]),    
    "Rikenellaceae"=unname(palo["gray"]),
    # "Coriobacteriaceae"=unname(palo["gray"]),    
    "Other"="lightgray"

  )
}


pal.es <- function () {

  palo <- pal.family()
  
  c(
         "ES-Bact"=unname(palo[["Bacteroidaceae"]]),
         "ES-Firm"=unname(palo[["Lachnospiraceae"]]),
	 "ES-Prev"=unname(palo[["Prevotellaceae"]]),
	 "ES-Bifi"=unname(palo[["Bifidobacteriaceae"]]),
	 "ES-Esch"=unname(palo[["Enterobacteriaceae"]])
  )
}


plotfunes <- function (es, nbreaks=5, esname) {
  A <- colData(TSE)[, paste0("nmf",es)]
  br <- unique(c(min(A), quantile(A[A>min(A)], probs=seq(0, 1, length=nbreaks))))
  ptab <- prop.table(table(cut(A, breaks=br, include.lowest=TRUE, right=FALSE), colData(TSE)$ARG_burd), 1)
  df <- reshape2::melt(ptab)
  pal <- rev(c("red", "pink", "lightblue", "blue"))
  ggplot(df, aes(x=Var1, y=value, group=Var2, color=Var2)) +
    geom_line() +
    geom_point() +  
    scale_color_manual(values=pal) +
    scale_y_continuous(label=scales::percent, limits=c(0,0.7)) +    
    labs(title=esname, x="NMF quantile", y="Share of the ARG group (%)") 
}

plotfun1 <- function (tax, nbreaks=5) {
  A <- assay(altExp(TSE, "Genus"), "relabundance")[tax,]
  br <- c(0, quantile(A[A>0], probs=seq(0, 1, length=nbreaks)))
  ptab <- prop.table(table(cut(A, breaks=br, include.lowest=TRUE, right=FALSE), colData(TSE)$ARG_burd), 1)
  df <- reshape2::melt(ptab)
  pal <- rev(c("red", "pink", "lightblue", "blue"))

  ggplot(df, aes(x=Var1, y=value, group=Var2, color=Var2)) +
    geom_line() +
    geom_point() +  
    scale_color_manual(values=pal) +
    scale_y_continuous(label=scales::percent) +    
    labs(title=tax, x="Abundance quantile", y="Share of the ARG group (%)") 
}

plotfun1f <- function (tax, nbreaks=5) {
  A <- assay(altExp(TSE, "Family"), "relabundance")[tax,]
  br <- unique(c(min(A), quantile(A[A>min(A)], probs=seq(0, 1, length=nbreaks))))
  ptab <- prop.table(table(cut(A, breaks=br, include.lowest=TRUE, right=FALSE), colData(TSE)$ARG_burd), 1)
  df <- reshape2::melt(ptab)
  pal <- rev(c("red", "pink", "lightblue", "blue"))

  p <- ggplot(df, aes(x=Var1, y=value, group=Var2, color=Var2)) +
    geom_line() +
    geom_point() +  
    scale_color_manual(values=pal) +
    scale_y_continuous(label=scales::percent, limits=c(0,1)) +    
    labs(title=tax, x="Abundance quantile", y="Share of the ARG group (%)")

  return(p)
}

plotfun2 <- function (tax, nbreaks=5) {
  A <- assay(altExp(TSE, "Genus"), "relabundance")[tax,]
  br <- c(0, quantile(A[A>0], probs=seq(0, 1, length=nbreaks)))
  ptab <- prop.table(table(cut(A, breaks=br, include.lowest=TRUE, right=FALSE), colData(TSE)$ARG_burd), 1)
  df <- reshape2::melt(ptab)
  pal <- rev(c("red", "pink", "lightblue", "blue"))

  # Barplot version
  ggplot(df, aes(x=Var1, y=value, fill=Var2)) +
    geom_bar(stat="identity", position="stack", color="black") +
    scale_fill_manual(values=pal) +
    labs(title=tax, x="Abundance quantile", y="Share of the ARG group (%)")     
}


plotfun3 <- function (tax, nbreaks=5) {
  A <- assay(altExp(TSE, "Genus"), "relabundance")[tax,]
  df <- as.data.frame(colData(TSE))
  br <- c(0, quantile(A[A>0], probs=seq(0, 1, length=nbreaks)))
  df$Quantile <- cut(A, breaks=br, include.lowest=TRUE, right=FALSE)

  # Barplot version
  ggplot(df, aes(x=Quantile, y=SUM_norm)) +
    geom_boxplot(outlier.shape=NA) +
    geom_jitter(alpha=0.2, width=0.2) +
    scale_y_log10() + 
    labs(title=tax, x="Abundance quantile", y="ARG (RPMK)")   
}

plotfun4 <- function (tax, nbreaks=5) {
  A <- assay(altExp(TSE, "Genus"), "relabundance")[tax,]
  df <- as.data.frame(colData(TSE))
  br <- c(0, quantile(A[A>0], probs=seq(0, 1, length=nbreaks)))
  df$Abundance <- A

  # Barplot version
  ggplot(df, aes(x=Abundance, y=SUM_norm)) +
    geom_point(alpha=0.2, width=0.2) +
    scale_x_log10(label=scales::percent) +
    scale_y_log10() +     
    labs(title=tax, x="Relative Abundance (%)", y="ARG (RPMK)")   
}

plotfun5 <- function (tax) {
  df <- as.data.frame(colData(TSE))
  df$counts <- assay(altExp(TSE, "Genus"), "counts")[tax,]
  df$clr <- assay(altExp(TSE, "Genus"), "clr")[tax,]
  df$relabundance <- assay(altExp(TSE, "Genus"), "relabundance")[tax,]
  df$Abundance <- df$clr
  # df <- df %>% filter(counts>0)
  ggplot(df, aes(x = ARG_burd, y = Abundance)) +
       geom_boxplot() +       
       geom_jitter(width=0.2, alpha=0.2) +
       # scale_y_log10(labels=scales::percent) +               
       labs(title=tax)    
}

plotfun <- function (i, DF) {
  DF$abundance <- DF[,i]

  # For plotting, replace zeroes with small constant
  # but so that it doesn't affect line fit
  DF$abundance2 <- DF[,i] <- DF$abundance
  DF$abundance2[DF$abundance2==0] <- min(DF$abundance2[DF$abundance2>0])
  
  ggplot(DF, aes(x = abundance, y = ARG_load)) +
    geom_point(aes(x=abundance2), alpha=0.1, size=1, color="gray") +
    geom_smooth(color = "blue", fill = "lightblue", method = "lm") +  
    geom_hline(yintercept=median(DF$ARG_load), linetype="dashed") +
    scale_y_continuous(trans="log", breaks=c(10, 100, 1000),
      limits=c(
      floor(min(DF$ARG_load)),
      ceiling(max(DF$ARG_load)))) +
    scale_x_continuous(trans="log10", breaks=c(0.1, 1, 10, 100)/100, labels=scales::percent) +     
    labs(x=paste(i, "(%)"), y="ARG load (RPKM)") +  
    theme_classic() +
    # xlim(4,7) +
    theme(legend.position="none")
}

find.top.taxa <- function(x,taxa){
 top.taxa <- tax_glom(x, taxa)
 otu <- otu_table(t(top.taxa)) %>% microbiome::transform("clr") # remove the transformation if using a merge_sample object
 tax <- tax_table(top.taxa)
 j <- apply(otu,1,which.max)
 k <- j[!duplicated(j)]
 l <- (tax[k,])
 m <- data.frame(otu[,k])
 s <- as.name(taxa)
 colnames(m)  = l[,taxa]
 n <- colnames(m)[apply(m,1,which.max)]
 m[,taxa] <- n
 return(m)
 gc()
}


## Cox ********************************************* ####

## Generate formula, fit Cox PH and wrangle results.
# data  = data.frame with all needed variables as columns
# predictors  = vector of variables whose effect we are interested in
# covariates  = vector of controlled covariates
# status  = 0/1 corresponding to no event/event; alive/dead etc.
# time_to_event  = follow-up time until event or censoring
# splines  = TRUE/FALSE; should the effect of predictor be modelled non-linearly?
# If TRUE, both linear and non-linear models are fit and for each predictor the fit with lower p-value is chosen.
# normalize  = TRUE/FALSE; standardize (x-mean(x)/sd(x)) predictor before fitting model?
# This makes comparison of hazard ratios easier
# test_ph_assumption  = TRUE/FALSE; if TRUE the function also returns a data frame with variable pairs that potentially
# do not meet the assumption of proportional hazards.


cox_regression <- function(data,
   predictors,
   covariates,
   status  = "DEATH",
   time_to_event  = "DEATH_AGEDIFF",
   splines,
   normalize) {


 if(normalize) {

 if(class(data[, predictors])  ==  "numeric") {

  x <- data[, predictors]
  data[, predictors] <- (x - mean(x, na.rm  = T))/sd(x, na.rm  = T)
 } else {
  data[, predictors] <- apply(data[, predictors], 2, FUN  = function(x) {(x - mean(x, na.rm  = T))/sd(x, na.rm  = T) })
 }

 }

 ## Formulas ***************************
 linear_formulas <- lapply(predictors, function(x) {
 formula_data <- deparse(substitute(data))

 formula <- paste0("Surv(", formula_data, "$", time_to_event, ", ",formula_data,"$", status, ")  ~  ",paste(covariates, collapse  = " + "), "  +  ",x)

 return(formula)

 }) %>%
 set_names(predictors)

 if(splines) {
 spline_formulas <- lapply(predictors, function(x) {
  formula_data <- deparse(substitute(data))

  formula <- paste0("Surv(", formula_data, "$", time_to_event, ", ",formula_data,"$", status, ")  ~  ",paste(covariates, collapse  = " + "), "  +  pspline(",x, ")")

  return(formula)

 }) %>%
  set_names(predictors)
 }

 ## Cox regression *********************
 print("Cox")
 linear_cox_fit <- lapply(linear_formulas, function(x) {

 coxph(as.formula(x), data = data, x = TRUE)

 })

 if(splines) {
 spline_cox_fit <- lapply(spline_formulas, function(x) {

  coxph(as.formula(x), data = data, x = TRUE)

 })

 return(list(linear_cox_fit  = linear_cox_fit, spline_cox_fit  = spline_cox_fit))
 } else {


 return(linear_cox_fit  = linear_cox_fit)
 }



}


cox_results <- function(data,
 predictors,
 covariates,
 status  = "DEATH",
 time_to_event  = "DEATH_AGEDIFF",
 alpha_level,
 splines,
 fit_list) {


 linear_cox_fit <- fit_list[["linear_cox_fit"]]
 if(splines) {spline_cox_fit <- fit_list[["spline_cox_fit"]]}


 ## Results *****************************
 print("Results")
 results <- lapply(predictors, function(x) {

 df <- summary(linear_cox_fit[[x]])$coefficients %>% as.data.frame()
 df <- df[nrow(df), ] %>%
  select(coef, "se(coef)", "z", "Pr(>|z|)") %>%
  set_colnames(c("coef", "se_coef", "test_statistic_value", "p")) %>%
  mutate(test_statistic  = "Wald")



 if(splines) {
  spline_df <- summary(spline_cox_fit[[x]])$coefficients %>%
 as.data.frame() %>%
 select(coef, "se(coef)", "Chisq", p) %>%
 set_colnames(c("coef", "se_coef", "test_statistic_value", "p")) %>%
 mutate(test_statistic  = "Chisq") %>%
 slice(nrow(.))

  # filter less significant association
  df <- rbind(df, spline_df) %>%
 mutate(association  = c("linear", "non-linear"))
 }



 df <- df %>%
  mutate(predictor  = x)




 }) %>%
 do.call(rbind, .)

 # Multiple testing correction
 results <- results %>%
 group_by(association) %>%
 mutate(P_adjusted  = p.adjust(p, "BH")) %>%
 ungroup() %>%
 group_by(predictor)


 # Results in neat form for presentation
 neat_results <- results %>%
 # filter(p  ==  min(p)) %>%
 ungroup() %>%
 mutate(HR  = round(exp(coef),3)) %>%
 mutate(HR_lower_95  = round(exp(coef - 1.96*se_coef), 3),
 HR_upper_95  = round(exp(coef  +  1.96*se_coef), 3),
 P  = round(p, 5),
 Coefficient  = round(coef, 3),
 "Coefficient SE"  = round(se_coef, 3)) %>%
 mutate(HR  = paste0(HR, " (95% CI, ", HR_lower_95, "-", HR_upper_95, ")")) %>%
 select(Predictor  = predictor, Coefficient, "Coefficient SE", HR, "P_adjusted", test_statistic_value, test_statistic) %>%
 mutate(HR  = ifelse(is.na(Coefficient), NA, HR), Association  = ifelse(is.na(Coefficient), "non-linear", "linear"))  %>%
 filter(P_adjusted < alpha_level) %>%
 arrange(Association, P_adjusted) %>%
 set_colnames(c("Predictor", "Coefficient", "Coefficient SE", "HR", "P (adjusted)", "Test Statistic Value", "Test Statistic", "Association"))


 # Results in a form more convenient for further manipulations
 results <- results %>%
 ungroup %>%
 mutate(PH  = exp(coef)) %>%
 mutate(p_adj  = P_adjusted) %>%
 mutate(linearity  = ifelse(is.na(coef), "non-linear", "linear"),
 posneg  = ifelse(is.na(coef), "nonlinear",
   ifelse(coef > 0, "bad", "good")))

 linear_results <- results %>% filter(linearity  ==  "linear")
 non_linear_results <- results %>% filter(linearity  ==  "non-linear")

 if(nrow(neat_results)  ==  0) {
 return(list(results  = results,
  linear_results  = linear_results,
  non_linear_results  = non_linear_results))

 }

 return(list(neat_results  = neat_results,
 results  = results,
 linear_results  = linear_results,
 non_linear_results  = non_linear_results))


}

cox_ph_assumptions <- function(data,
   predictors,
   covariates,
   status  = "DEATH",
   time_to_event  = "DEATH_AGEDIFF",
   splines,
   fit_list) {


 linear_cox_fit <- fit_list[["linear_cox_fit"]]
 if(splines) {spline_cox_fit <- fit_list[["spline_cox_fit"]]}


 ph_assumption <-  lapply(predictors, function(m) {
 test <- cox.zph(linear_cox_fit[[m]])

 p_values <- test$table[, "p"]

 # significant cases
 x <- which(p_values < 1)

 if(length(x)  ==  0) {
  return(NULL)
 }

 df <- data.frame(feature  = m, variable_not_ph  = names(x), p_value  = p_values[x])

 }) %>%
 do.call(rbind, .) %>%
 mutate(p_adj  = p.adjust(p_value, "BH")) %>%
 filter(p_value < alpha_level)


 return(ph_assumption)

}


# Wrap the above three functions into one
cox_wrapper <- function(data,
 predictors,
 covariates,
 status  = "DEATH",
 time_to_event  = "DEATH_AGEDIFF",
 alpha_level  = alpha_level,
 splines,
 normalize,
 test_ph_assumption  = FALSE) {


 if(normalize) {

 if(class(data[, predictors])  ==  "numeric") {

  x <- data[, predictors]
  data[, predictors] <- (x - mean(x, na.rm  = T))/sd(x, na.rm  = T)
 } else {
  data[, predictors] <- apply(data[, predictors], 2, FUN  = function(x) {(x - mean(x, na.rm  = T))/sd(x, na.rm  = T) })
 }

 }

 ## Formulas ***************************
 linear_formulas <- lapply(predictors, function(x) {
 formula_data <- deparse(substitute(data))

 formula <- paste0("Surv(", formula_data, "$", time_to_event, ", ",formula_data,"$", status, ")  ~  ",paste(covariates, collapse  = " + "), "  +  ",x)

 return(formula)

 }) %>%
 set_names(predictors)

 if(splines) {
 spline_formulas <- lapply(predictors, function(x) {
  formula_data <- deparse(substitute(data))

  formula <- paste0("Surv(", formula_data, "$", time_to_event, ", ",formula_data,"$", status, ")  ~  ",paste(covariates, collapse  = " + "), "  +  pspline(",x, ")")

  return(formula)

 }) %>%
  set_names(predictors)
 }

 ## Cox regression *********************
 print("Cox")
 linear_cox_fit <- lapply(linear_formulas, function(x) {

 coxph(as.formula(x), data = data, x = TRUE)

 })

 if(splines) {
 spline_cox_fit <- lapply(spline_formulas, function(x) {

  coxph(as.formula(x), data = data, x = TRUE)

 })
 }


 ## Check PH assumptions ****************
 if(test_ph_assumption) {
 print("PH assumptions")
 ph_assumption <-  lapply(predictors, function(m) {
  test <- cox.zph(linear_cox_fit[[m]])

  p_values <- test$table[, "p"]

  # significant cases
  x <- which(p_values < 1)

  if(length(x)  ==  0) {
 return(NULL)
  }

  df <- data.frame(feature  = m, variable_not_ph  = names(x), p_value  = p_values[x])

 }) %>%
  do.call(rbind, .) %>%
  mutate(p_adj  = p.adjust(p_value, "BH")) %>%
  filter(p_value < alpha_level)

 }

 ## Results *****************************
 print("Results")
 results <- lapply(predictors, function(x) {

 df <- summary(linear_cox_fit[[x]])$coefficients %>% as.data.frame()
 df <- df[nrow(df), ] %>%
  select(coef, "se(coef)", "z", "Pr(>|z|)") %>%
  set_colnames(c("coef", "se_coef", "test_stat_value", "p")) %>%
  mutate(test_stat  = "Wald")



 if(splines) {
  spline_df <- summary(spline_cox_fit[[x]])$coefficients %>%
 as.data.frame() %>%
 select(coef, "se(coef)", Chisq, p) %>%
 set_colnames(c("coef", "se_coef", "test_stat_value", "p")) %>%
 mutate(test_stat  = "Chisq") %>%
 slice(nrow(.))

  df <- rbind(df, spline_df) %>%
 mutate(association  = c("linear", "non-linear"))
 }



 df <- df %>%
  mutate(predictor  = x)




 }) %>%
 do.call(rbind, .)

 # Multiple testing correction
 results <- results %>%
 group_by(association) %>%
 mutate(P_adjusted  = p.adjust(p, "BH")) %>%
 ungroup() %>%
 group_by(predictor)








 # Results in neat form for presentation
 neat_results <- results %>%
 # filter(p  ==  min(p)) %>%
 ungroup() %>%
 mutate(HR  = round(exp(coef),3)) %>%
 mutate(HR_lower_95  = round(exp(coef - 1.96*se_coef), 3),
 HR_upper_95  = round(exp(coef  +  1.96*se_coef), 3),
 P  = round(p, 5),
 Coefficient  = round(coef, 3),
 "Coefficient SE"  = round(se_coef, 3)) %>%
 mutate(HR  = paste0(HR, " (95% CI, ", HR_lower_95, "-", HR_upper_95, ")")) %>%
 select(Predictor  = predictor, Coefficient, "Coefficient SE", HR, "P_adjusted", "test_stat_value", "test_stat") %>%
 mutate(HR  = ifelse(is.na(Coefficient), NA, HR), Association  = ifelse(is.na(Coefficient), "non-linear", "linear"))  %>%
 filter(P_adjusted < alpha_level) %>%
 arrange(Association, P_adjusted) %>%
 set_colnames(c("Predictor", "Coefficient", "Coefficient SE", "HR", "P (adjusted)", "Test Statistic Value", "Test Statistic", "Association"))


 # Results in a form more convenient for further manipulations
 results <- results %>%
 ungroup %>%
 mutate(PH  = exp(coef)) %>%
 mutate(p_adj  = P_adjusted) %>%
 mutate(linearity  = ifelse(is.na(coef), "non-linear", "linear"),
 posneg  = ifelse(is.na(coef), "nonlinear",
   ifelse(coef < 0, "bad", "good")))






 if(nrow(neat_results)  ==  0) {
 return(list(results  = results))
 }

 if(test_ph_assumption) {

 if(nrow(neat_results)  ==  0) {
  return(list(results  = results,
 ph_assumption  = ph_assumption))
 }

 return(list(neat_results  = neat_results,
  results  = results,
  ph_assumption  = ph_assumption))
 }

 return(list(neat_results  = neat_results,
 results  = results))

}


# wrapper2: returns covariate results as well
cox_wrapper_with_covariates <- function(data,
   predictors ,
   covariates,
   status  = "DEATH",
   time_to_event  = "DEATH_AGEDIFF",
   alpha_level  = alpha_level,
   splines,
   normalize,
   p_adjust  = TRUE) {


 if(normalize) {

 if(class(data[, predictors])  ==  "numeric") {

  x <- data[, predictors]
  data[, predictors] <- (x - mean(x, na.rm  = T))/sd(x, na.rm  = T)
 } else {
  data[, predictors] <- apply(data[, predictors],
   2,
   FUN  = function(x) {
   (x - mean(x, na.rm  = T))/sd(x, na.rm  = T)
   })
 }


 for(co in covariates) {

  if(class(data[, co])  ==  "numeric") {

 x <- data[, co]
 data[, co] <- (x - mean(x, na.rm  = T))/sd(x, na.rm  = T)
  }

 }

 }

 ## Formulas ***************************
 linear_formulas <- lapply(predictors, function(x) {
 formula_data <- deparse(substitute(data))

 formula <- paste0("Surv(", formula_data, "$", time_to_event, ", ",formula_data,"$", status, ")  ~  ",paste(covariates, collapse  = " + "), "  +  ",x)

 return(formula)

 }) %>%
 set_names(predictors)

 if(splines) {
 spline_formulas <- lapply(predictors, function(x) {
  formula_data <- deparse(substitute(data))

  formula <- paste0("Surv(", formula_data, "$", time_to_event, ", ",formula_data,"$", status, ")  ~  ",paste(covariates, collapse  = " + "), "  +  pspline(",x, ")")

  return(formula)

 }) %>%
  set_names(predictors)
 }

 ## Cox regression *********************
 print("Cox")
 linear_cox_fit <- lapply(linear_formulas, function(x) {

 coxph(as.formula(x), data = data, x = TRUE)

 })

 if(splines) {
 spline_cox_fit <- lapply(spline_formulas, function(x) {

  coxph(as.formula(x), data = data, x = TRUE)

 })
 }


 ## Results *****************************
 print("Results")
 results <- lapply(predictors, function(x) {

 df <- summary(linear_cox_fit[[x]])$coefficients %>% as.data.frame()
 df <- df %>%
  rownames_to_column(var  = "variable") %>%
  select(variable, coef, "se(coef)", "z", "Pr(>|z|)") %>%
  set_colnames(c("variable", "coef", "se_coef", "test_stat_value", "p")) %>%
  mutate(test_stat  = "Wald")



 # if(splines) {
 #  spline_df <- summary(spline_cox_fit[[x]])$coefficients %>%
 # as.data.frame() %>%
 # rownames_to_column(var  = "variable") %>%
 # select(variable, coef, "se(coef)", Chisq, p) %>%
 # set_colnames(c("variable", "coef", "se_coef", "test_stat_value", "p")) %>%
 # mutate(test_stat  = "Chisq") %>%
 # slice(-(nrow(.) - 1))
 #
 #  df <- rbind(df, spline_df) %>%
 # mutate(association  = c("linear", "non-linear"))
 # }



 return(df)

 }) %>%
 do.call(rbind, .)

 # Multiple testing correction
 # results <- results %>%
 # mutate(P_adjusted  = p.adjust(p, "BH"))











 # Results in neat form for presentation
 neat_results <- results %>%
 # filter(p  ==  min(p)) %>%
 ungroup() %>%
 mutate(HR  = round(exp(coef),3)) %>%
 mutate(HR_lower_95  = round(exp(coef - 1.96*se_coef), 3),
 HR_upper_95  = round(exp(coef  +  1.96*se_coef), 3),
 # P  = round(p, 5),
 P  = p,
 Coefficient  = round(coef, 3),
 "Coefficient SE"  = round(se_coef, 3),
 Variable  = variable) %>%
 mutate(HR  = paste0(HR, " (95% CI, ", HR_lower_95, "-", HR_upper_95, ")")) %>%
 select(Variable, Coefficient, "Coefficient SE", HR, "P", "test_stat_value", "test_stat") %>%
 mutate(HR  = ifelse(is.na(Coefficient), NA, HR), Association  = ifelse(is.na(Coefficient), "non-linear", "linear"))  %>%
 # filter(P_adjusted < alpha_level) %>%
 arrange(Association, P) %>%
 set_colnames(c("Variable", "Coefficient", "Coefficient SE", "HR", "P", "Test Statistic Value", "Test Statistic", "Association"))


 # Results in a form more convenient for further manipulations
 results <- results %>%
 ungroup %>%
 mutate(PH  = exp(coef)) %>%
 mutate(linearity  = ifelse(is.na(coef), "non-linear", "linear"))

 if(nrow(neat_results)  ==  0) {
 return(list(results  = results))
 }

 if(test_ph_assumption) {

 if(nrow(neat_results)  ==  0) {
  return(list(results  = results,
 ph_assumption  = ph_assumption))
 }

 return(list(neat_results  = neat_results,
  results  = results,
  ph_assumption  = ph_assumption))
 }

 return(list(neat_results  = neat_results,
 results  = results))

}

# Wrapper for categorical predictor. Compares last level to the first
cox_wrapper_categorical <- function(data,
   predictors ,
   covariates,
   status  = "DEATH",
   time_to_event  = "DEATH_AGEDIFF",
   alpha_level  = alpha_level,
   normalize) {


 if(normalize) {

 if(class(data[, predictors])  ==  "numeric") {

  x <- data[, predictors]
  data[, predictors] <- (x - mean(x, na.rm  = T))/sd(x, na.rm  = T)
 }

 for(co in covariates) {

  if(class(data[, co])  ==  "numeric") {

 x <- data[, co]
 data[, co] <- (x - mean(x, na.rm  = T))/sd(x, na.rm  = T)
  }

 }

 }

 ## Formulas ***************************
 print("formulas")
 linear_formulas <- lapply(predictors, function(x) {
 formula_data <- deparse(substitute(data))

 formula <- paste0("Surv(", formula_data, "$", time_to_event, ", ",formula_data,"$", status, ")  ~  ",paste(covariates, collapse  = " + "), "  +  ",x)

 return(formula)

 }) %>%
 set_names(predictors)


 ## Cox regression *********************
 print("Cox")
 linear_cox_fit <- lapply(linear_formulas, function(x) {

 coxph(as.formula(x), data = data, x = TRUE)

 })


 ## Results *****************************
 print("Results")
 results <- lapply(predictors, function(x) {

 df <- summary(linear_cox_fit[[x]])$coefficients %>% as.data.frame()
 df <- df %>%
  rownames_to_column(var  = "variable") %>%
  select(variable, coef, "se(coef)", "z", "Pr(>|z|)") %>%
  set_colnames(c("variable", "coef", "se_coef", "test_stat_value", "p")) %>%
  mutate(test_stat  = "Wald")



 # if(splines) {
 #  spline_df <- summary(spline_cox_fit[[x]])$coefficients %>%
 # as.data.frame() %>%
 # rownames_to_column(var  = "variable") %>%
 # select(variable, coef, "se(coef)", Chisq, p) %>%
 # set_colnames(c("variable", "coef", "se_coef", "test_stat_value", "p")) %>%
 # mutate(test_stat  = "Chisq") %>%
 # slice(-(nrow(.) - 1))
 #
 #  df <- rbind(df, spline_df) %>%
 # mutate(association  = c("linear", "non-linear"))
 # }



 return(df)

 }) %>%
 do.call(rbind, .)

 # Results in neat form for presentation
 print("neat_resiults")
 neat_results <- results %>%
 # filter(p  ==  min(p)) %>%
 ungroup() %>%
 mutate(HR  = round(exp(coef),3)) %>%
 mutate(HR_lower_95  = round(exp(coef - 1.96*se_coef), 3),
 HR_upper_95  = round(exp(coef  +  1.96*se_coef), 3),
 # P  = round(p, 5),
 P  = p,
 Coefficient  = round(coef, 3),
 "Coefficient SE"  = round(se_coef, 3),
 Variable  = variable) %>%
 mutate(HR  = paste0(HR, " (95% CI, ", HR_lower_95, "-", HR_upper_95, ")")) %>%
 select(Variable, Coefficient, "Coefficient SE", HR, "P", "test_stat_value", "test_stat") %>%
 mutate(HR  = ifelse(is.na(Coefficient), NA, HR), Association  = ifelse(is.na(Coefficient), "non-linear", "linear"))  %>%
 # filter(P_adjusted < alpha_level) %>%
 arrange(Association, P) %>%
 set_colnames(c("Variable", "Coefficient", "Coefficient SE", "HR", "P", "Test Statistic Value", "Test Statistic", "Association"))


 # Results in a form more convenient for further manipulations
 print("results")
 results <- results %>%
 ungroup %>%
 mutate(PH  = exp(coef)) %>%
 mutate(linearity  = ifelse(is.na(coef), "non-linear", "linear"))




 print("retursn")

 if(nrow(neat_results)  ==  0) {
 return(list(results  = results))
 }


 return(list(neat_results  = neat_results,
 results  = results))

}


# Hazard ratio plotter
# Note that the fit needs to be explicitely written and not from a loop.
plot_HR <- function(fit,
 pred,
 data = microbiome::meta(pseq),
 ref_value = NULL,
 ribbon_color  = "#5C88DA",
 title  = TRUE,
 return_data  = FALSE) {


 # clean predicator name
 pred_name <- clean_genus_names(pred)

 # rename data
 mydata <- data

 # not sure if this is redundant
 hr <- smoothHR(data = (mydata)[, c(covariates, pred, "DEATH_AGEDIFF", "DEATH")], coxfit  = fit)

 # use predictor median as reference value, unless defined in function call
 if(is.null(ref_value)) {
 ref_value <- median(mydata[, pred])
 }

 # range of predictor for graph
 r <- range(hr$dataset[, pred])



 # simulate values and edit data frame
 fit_sim <- predict(hr,
 predictor  = pred,
 pred.value  = ref_value,
 prediction.values  = seq(from = r[1], to = r[2], length.out = 100), data = mydata)

 fit_sim_df <- fit_sim %>% as.data.frame()
 colnames(fit_sim_df) <- c("pred", "ln_hr", "ln_lower95", "ln_upper95")



 # natural logarithm
 fit_sim_df <- fit_sim_df %>%
 mutate(hr  = exp(ln_hr), lower95 = exp(ln_lower95), upper95 = exp(ln_upper95)) %>%
 mutate(log2_hr = log2(hr), log2_lower95 = log2(lower95), log2_upper95 = log2(upper95))

 if(return_data) {
 return(fit_sim_df)
 }

 # plot
 if(!isTRUE(title)) {
 p <- ggplot(data.frame(fit_sim_df), aes(x = pred, y = log2_hr))  +
  # geom_line(size  = 1)  +
  geom_line()  +
  geom_hline(yintercept  = 0, linetype = "dashed")  +
  geom_ribbon(aes(ymin  = fit_sim_df$log2_lower95, ymax  = fit_sim_df$log2_upper95),
 linetype  = 2,
 alpha  = .5,
 fill  = ribbon_color)  +
  scale_y_continuous(breaks = -2:6, labels  = 2^(-2:6))  +
  labs(y = "HR for Death")  +
  coord_cartesian(xlim  = c(r[1], r[2]))



 } else {

 p <- ggplot(data.frame(fit_sim_df), aes(x = pred, y = log2_hr))  +
  # geom_line(size  = 1)  +
  geom_line()  +
  geom_hline(yintercept  = 0, linetype = "dashed")  +
  geom_ribbon(aes(ymin  = fit_sim_df$log2_lower95, ymax  = fit_sim_df$log2_upper95),
 linetype  = 2,
 alpha  = .5,
 fill  = ribbon_color)  +
  scale_y_continuous(breaks = -2:6, labels  = 2^(-2:6))  +
  labs(y = "HR for Death", title = pred_name)  +
  coord_cartesian(xlim  = c(r[1], r[2]))
 }



 return(p)

}







# Generate Cox model formulas
cox_formula <- function(data  = meta, predictor, covariates, splines  = TRUE) {

 pred <- ifelse(class(data[, predictor]) == "numeric" & splines  ==  TRUE, paste0(" +  pspline(", predictor,")"), paste(" + ",predictor))

 formula_data <- deparse(substitute(data))

 formula <- paste0("Surv(",formula_data,"$DEATH_AGEDIFF, ",formula_data,"$Death)  ~ ",paste(covariates, collapse  = " + "),pred)

 return(formula)

}


functional_cox_formula <- function(data, predictor, covariates) {

 pred <- ifelse(class(data[, predictor]) == "numeric", paste0(" +  pspline(", predictor,")"), paste(" + ",predictor))

 formula <- paste("Surv(meta$DEATH_AGEDIFF, meta$Death)  ~ ",paste(covariates, collapse  = " + "),pred)

 return(formula)

}


# Get results
get_cox_results <- function(var, model_list, data = NULL) {

 if(!is.null(data)) {
 meta <- data
 }

 # separate cases for numeric and factor variables
 if(class(meta[, var]) == "numeric") {

 # number of rows in the summary.
 # Necessary to extract rows form summary by index as
 # long predictor names can be cut


 df <- summary(model_list[[var]])$coefficients %>%
  as.data.frame %>%
  rownames_to_column(var  = "predictor")

 df %>%
  filter(predictor  ==  var) %>%
  select(predictor, coef, "se(coef)", "Pr(>|z|)") %>%
  set_colnames(c("predictor", "coef", "se_coef", "p"))

 } else if(class(meta[, var]) == "factor") {

 summary(model_list[[var]])$coefficients %>%
  as.data.frame %>%
  rownames_to_column(var  = "predictor") %>%
  filter(grepl(var, predictor)) %>%
  select("Pr(>|z|)", coef, "se(coef)") %>%
  set_colnames(c("p", "coef", "se_coef")) %>%
  mutate(predictor = var, linearity = NA, predictor = var, pred_class = "factor")

 } else warning("Check variable class")



}



# Cox regression cross validatation. Computes the ROC and AUC.
cox_cv <- function(data,
  predictor,
  response_var,
  response_time  = "DEATH_AGEDIFF",
  covariates  = covariates,
  n_folds  = 5,
  seed  = 11235) {

 # set seed for reproducibility
 set.seed(seed)

 # Data **************************************************************************

 # mydata is data with used variables filtered
 mydata <- data[, c(response_var,
 response_time,
 predictor,
 covariates)] %>%
 drop_na()


 # modify variable names so upcoming function calls will work
 mydata$response_time  = mydata[, response_time]
 mydata[, response_time] <- NULL
 mydata$response_var  = mydata[, response_var]
 mydata[, response_var] <- NULL



 # Folds *************************************************************************
 # Randomly shuffle the data
 mydata <- mydata[sample(nrow(mydata)),]

 # Create n_folds equally sized folds
 folds <- cut(seq(1, nrow(mydata)), breaks  = n_folds, labels  = FALSE)


 # Cross validate ***************************************************************
 print("CV")
 cox_cv <- lapply(1:n_folds, function(fold) {
 print(fold)

 # training samples
 train <- folds != fold

 # fit cox

 formula <- paste0("Surv(response_time, response_var)  ~  BL_AGE  +  BMI  +  MEN  +  pspline(",predictor,")") %>%
  as.formula()

 cox_train <- coxph(formula  = formula,
  data  = mydata[train, ],
  x  = TRUE)

 # Predict in test set
 cox_pred <- predict(cox_train, newdata  = mydata[!train, ])

 return(cox_pred)

 })


 # ROC & AUC ********************************************************************
 print("ROC & AUC")
 ROC_df <- lapply(1:n_folds, function(fold) {
 print(fold)

 train <- folds != fold

 # Combine
 dat <- cbind(mydata[!train, ], survival_pred  = cox_cv[[fold]])

 roc_res <- with(dat,
 survivalROC(Stime  = response_time,
  status   = response_var,
  marker   = survival_pred,
  predict.time  = max(mydata[, "response_time"]),
  method   = "KM"))  # KM method without smoothing

 df <- data.frame(predictor  = predictor,
 FPR  = roc_res[["FP"]],
 TPR  = roc_res[["TP"]],
 fold  = rep(fold, length(roc_res[["TP"]])),
 AUC  = rep(roc_res[["AUC"]],length(roc_res[["TP"]])))
 return(df)

 }) %>% do.call(rbind, .)

 return(ROC_df)
}

## Survival Random Forest ************************** ####

# Compute the harrell's c in n_folds cross-validation.
random_survival_forest_cv <- function(data,
  predictor_list,
  response_var,
  response_time  = "DEATH_AGEDIFF",
  n_folds  = 5,
  seed  = 11235,
  save_file  = "") {




 # set seed for reproducibility
 set.seed(seed)

 # Data **************************************************************************
 # save results to
 survival_RF <- list()

 # mydata is data with used variables filtered
 mydata <- data[, c(response_var,
 response_time,
 predictor_list %>% unlist %>% unique())] %>%
 drop_na()

 # drop subjects with no response time
 mydata <- mydata %>% drop_na(response_time)

 # modify variable names so upcoming function calls will work
 mydata$response_time  = mydata[, response_time]
 mydata[, response_time] <- NULL
 mydata$response_var  = mydata[, response_var]
 mydata[, response_var] <- NULL


 # Folds *************************************************************************
 # Randomly shuffle the data
 mydata <- mydata[sample(nrow(mydata)),]

 # Create n_folds equally sized folds
 folds <- cut(seq(1, nrow(mydata)), breaks  = n_folds, labels  = FALSE)


 # Loop over list of predictors. Use for() to save intermediate steps and ease debugging
 for(subs in names(predictor_list)) {
 print(subs)

 # subset data
 rf_data <- mydata[, c("response_time", "response_var", predictor_list[[subs]])]


 # Cross validate ***************************************************************
 print("CV")
 rf_cv <- lapply(1:n_folds, function(fold) {
  print(fold)

  # training samples
  train <- folds != fold

  # grow random forest
  forest <- rfsrc(Surv(response_time, response_var)  ~  .,
  data  = rf_data[train, ],
  importance  = TRUE)

  # predict in test set
  pred <- predict(forest,
  newdata  = rf_data[!train, ])

  return(list(forest_importance  = forest$importance,
 prediction_values  = pred$predicted,
 train_samples  = train,
 harrell_c  = 1 - forest$err.rate[length(forest$err.rate)]))

 })

 # ROC & AUC ********************************************************************
 print("ROC & AUC")
 ROC_df <- lapply(1:n_folds, function(fold) {
  print(fold)
  train <- folds != fold

  roc_data <- cbind(rf_data[!train, ], survival_pred  = rf_cv[[fold]]$prediction_values)
  roc_res <- with(roc_data,
  survivalROC(Stime  = response_time,
   status  = response_var,
   marker  = survival_pred,
   predict.time  = max(mydata[, "response_time"]),
   method  = "KM"))  # KM method without smoothing

  df <- data.frame(predictor  = subs,
  FPR  = roc_res[["FP"]],
  TPR  = roc_res[["TP"]],
  fold  = rep(fold, length(roc_res[["TP"]])),
  AUC  = rep(roc_res[["AUC"]],length(roc_res[["TP"]])))
  return(df)

 }) %>% do.call(rbind, .)


 # Results to list
 survival_RF[[subs]] <- list(cv  = rf_cv, roc  = ROC_df)

 }

 return(survival_RF)

}







## MISC ******************************************** ####


get_causes_of_death <- function(x) {

 if(class(x)  ==  "phyloseq") {
 causes <- meta(x)[, "K_TPKS"]
 } else {
 causes <- x
 }

 # Remove number part
 causes <- gsub("[0-9] + ", "", causes)

 for(i in 1:length(causes)) {
 if(is.na(causes[i])) {
  causes[i] <- NA
 } else if(causes[i]  ==  "C") {
  causes[i] <- "Cancer"
 } else if(causes[i]  ==  "G") {
  causes[i] <- "Neurological"
 } else if(causes[i]  ==  "I") {
  causes[i] <- "Cardiovascular"
 } else if(causes[i]  ==  "K") {
  causes[i] <- "Gastrointestinal"
 } else if(causes[i]  ==  "J") {
  causes[i] <- "Respiratory"
 } else if(!is.na(causes[i])) {
  causes[i] <- "Other"
 }
 }

 if(class(x)  ==  "phyloseq") {

 cause_of_death <- causes

 sample_data(x) <- cbind(meta(x), cause_of_death)
 return(x)

 } else {
 return(causes)
 }


}


# Phyloseq editor
# Takes a phyloseq object and adds pca components, diversity measures, edit taxa names

add_many_things_to_meta_data <- function(pseq, species_pseq) {

 my_pseq <- pseq


 ## Diversities
 my_pseq <- add_diversities_to_meta_data(my_pseq)


 ## Taxa names: Bacteroides (Bacteria) --> Bacteroides_Bacteria
 taxa_names(my_pseq) <- dirty_genus_names(taxa_names(my_pseq))


 ## CLR transformed abundances
 my_pseq <- add_clr_abundances_to_meta_data(my_pseq)


 ## Scaled PCA components, continuous
 my_pseq <- add_pca_components_to_meta_data(my_pseq, species_pseq)


 return(my_pseq)

}

add_pca_components_to_meta_data <- function(pseq, species_pseq, n_axes  = 3, check_existence  = TRUE) {

 set.seed(12345)

 my_pseq <- pseq

 # Check if PCA output file exists. If not, then compute
 if(file.exists("output/raw_output/species_pca.rds") & check_existence) {
 pca_clr <- readRDS(file  = "output/raw_output/species_pca.rds")
 } else {

 print("Doing PCA. This may take a while!")

 pca_clr <- prcomp(microbiome::transform(
  microbiome::transform(species_pseq,
  "compositional"),
  "clr"
 ) %>%
  abundances %>%
  t)

 # Flip PC1 to get a positive association with mortality
 pca_clr$rotation[, "PC1"] <- -1*pca_clr$rotation[, "PC1"]
 pca_clr$x[, "PC1"] <- -1*pca_clr$x[, "PC1"]

 }




 sample_data(my_pseq) <- cbind(sample_data(my_pseq), pca_clr$x %>%
   as_tibble %>%
   select(paste0("PC", 1:n_axes)))


 return(my_pseq)

}

add_clr_abundances_to_meta_data <- function(pseq) {

 my_pseq <- pseq

 otus_clr <- abundances(my_pseq) %>%
 microbiome::transform("compositional") %>%
 microbiome::transform("clr") %>%
 t

 colnames(otus_clr) <- colnames(otus_clr) %>%
 gsub("\\(", "_", .) %>%
 gsub("\\)", "", .) %>%
 gsub(" ", "", .)

 sample_data(my_pseq) <- cbind(meta(my_pseq), otus_clr)

 return(my_pseq)


}




add_diversities_to_meta_data <- function(pseq, species_pseq) {

 my_pseq <- pseq

 ## Diversities
 # species_pseq <- readRDS(species_level_phyloseq_path)
 species_diversities <- estimate_richness(species_pseq, measures  = c("Observed", "Shannon"))
 sample_data(my_pseq) <- cbind(meta(my_pseq), species_diversities[, c("Observed", "Shannon")])

 return(my_pseq)
}

get_pca <- function(species_pseq, check_existence  = TRUE) {

 set.seed(12345)


 # Check if PCA output file exists. If not, then compute
 if(file.exists("output/raw_output/species_pca.Rdata") & check_existence) {
 pca_clr <- readRDS(file  = "output/raw_output/species_pca.Rdata")
 } else {

 print("Doing PCA. This may take a while!")

 pca_clr <- prcomp(microbiome::transform(
  microbiome::transform(species_pseq,
  "compositional"),
  "clr"
 ) %>%
  abundances %>%
  t)

 # Flip PC1 to get a positive association with mortality
 pca_clr$rotation[, "PC1"] <- -1*pca_clr$rotation[, "PC1"]
 pca_clr$x[, "PC1"] <- -1*pca_clr$x[, "PC1"]

 }

 return(pca_clr)

}


# get subnet abundances, with CLR transformation clr  = TRUE
subnet_abundances <- function(pseq, subnet_components, clr  = "TRUE") {

 my_pseq <- pseq


 ## Add subnet components
 if(!exists("subnet_components")) subnet_components <- readRDS(network_identities_path)

 subnet_abundances <- abundances(my_pseq) %>% t %>% as.data.frame()

 # Aggregate subnet reads and remove the constituents
 for(i in unique(subnet_components$componentID)) {

 # Get subnet taxa
 subnet_taxa <- subnet_components %>%
  filter(componentID  ==  i) %>%
  pull(OTUunderscored)

 # Aggregate reads
 subnet_reads <- subnet_abundances[, subnet_taxa] %>%
  rowSums()

 # Remove taxa from otu table
 subnet_abundances <- subnet_abundances[, !(colnames(subnet_abundances) %in% subnet_taxa)]

 # Add aggregated reads
 subnet_abundances[, paste0("subnet_", i)] <- subnet_reads

 }


 if(isTRUE(clr)) {

 # CLR transform subnet abundances
 abundances <- subnet_abundances %>%
  t %>%
  microbiome::transform("compositional") %>%
  microbiome::transform("clr") %>%
  t

 } else {
 abundances <- subnet_abundances
 }



 return(abundances)
}


# Z-normalization
z_normalize <- function(x, m  = TRUE) {

 if(m  ==  TRUE) {

 return((x-mean(x, na.rm  = T))/sd(x, na.rm  = T))

 } else {
 return(x/sd(x, na.rm  = T))

 }


}

# Get percent of explained variation from prcomp x for axis n
PC_variation_explained <- function(x, n, round) {
 eigs <- x$sdev^2
 var <- eigs[n] / sum(eigs)
 return(100*round(var, round))
}

# wrapper: matrix to tibble with rownames to column
m_neat <- function(x, colnames  = NULL) {
 x <- x %>%
 as.data.frame() %>%
 rownames_to_column()

 if(!is.null(colnames)) {
 x <- x %>% set_colnames(colnames)
 }

 return(x)
}

# change variable names according to a df with variable explanations
change_colnames <- function(df, legend_df, col = 3) {
 cols <- colnames(df)
 new_cols <- legend_df[cols, col]
 colnames(df) <- new_cols
 return(df)
}
change_varname <- function(var, legend_df, col = 3) {
 row <- which(rownames(legend_df)  ==  var)
 sober_varname <- legend_df[row, col]
 return(sober_varname)
}

# remove empty rows and columns
trim_empty_rows <- function(df, value = 0) {

 df_reduced <- df[, colSums(df != value) > 0]
 df_reduced <- df_reduced[rowSums(df_reduced != value) > 0, ]

}

# clean taxa names: Bacteroides_Bacteria --> Bacteroides (Bacteria)
clean_genus_names <- function(genera, kingdom  = TRUE) {

 if(isTRUE(kingdom)){

 title <- gsub("_", " ", genera)
 title <- gsub("Bacteria", "(Bacteria)", title)
 title <- gsub("Viruses", "(Viruses)", title)
 title <- gsub("Archaea", "(Archaea)", title)

 for(i in 1:length(title)) {

  if(grepl("Plasmid", title[i])) {
 title[i] <- gsub(")Plasmid", "Plasmid)", title[i])
  }
 }
 title <- gsub("Plasmid", " Plasmid", title)
 } else {

 title <- gsub("_", "", genera)
 title <- gsub("Bacteria", "", title)
 title <- gsub("Viruses", " (Viruses)", title)
 title <- gsub("Archaea", " (Archaea)", title)
 title <- gsub("Plasmid", " (Plasmid)", title)
 title <- gsub("\\) \\(", " ", title)


 }



 return(title)
}
# dirty taxa names
dirty_genus_names <- function(x) {

 if(class(x)  ==  "phyloseq") {
 genera <- taxa_names(x)

 } else {
 x <- genera
 }

 genera_edited <- genera %>%
 gsub(" ", "_", .) %>%
 gsub("\\(", "", .) %>%
 sub("\\)", "", .) %>%
 gsub("-", "_", .) %>%
 gsub(",", "_", .)

 if(class(x)  ==  "phyloseq") {

 taxa_names(x) <- genera_edited
 return(x)

 } else {
 return(genera_edited)
 }


}

# Clean coefficient table
cox_clean <- function(cox_model) {
cox_DEATH_coefficients <- summary(cox_model)$coefficients
df_p <- summary(cox_model)$conf.int %>% data.frame()
df_p <- df_p[,1:4] %>% round(digits = 2)
df_p$Variables <- rownames(df_p)
tmp <- as.data.frame(cox_DEATH_coefficients) %>% dplyr::select("Pr(>|z|)")
df_p$p <- tmp$"Pr(>|z|)"
return(df_p)
}


# Clean coefficient table
cox_clean2 <- function(cox_model) {
  cox_DEATH_coefficients <- summary(cox_model)$coefficients
  df_p <- summary(cox_model)$conf.int %>% data.frame()
  df_p <- df_p[,1:4] %>% round(digits = 2)
  df_p$Variables <- rownames(df_p)
  tmp <- as.data.frame(cox_DEATH_coefficients) %>% dplyr::select("Pr(>|z|)")
  df_p$p <- tmp$"Pr(>|z|)"
  return(df_p)
}


ord_plot <- function(
 ord,
 group,
 ellipse  = TRUE,
 spyder  = FALSE,
 centroids  = FALSE,
 CI  = 0.75,
 palette  = NULL
 ){

 if(!is.factor(group)) group <- factor(group)

 # Define color palette: ----
 if(is.null(palette)){
 cols <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")
 palette <- colorRampPalette(cols)(nlevels(group))
 }

 # Define plotting frame: ----
 if(is.matrix(ord)){
 labels <- paste('Axis',1:2,sep = '-')
 Z  = ord
 } else if('dbrda' %in% class(ord)) {
 Z <- summary(ord)$sites[,1:2]
 labels <- paste('RDA',1:2,sep = '-')
 } else if(grepl('MDS',class(ord))){
 Z <- ord$points
 labels <- paste('MDS',1:2,sep = '-')
 }

 if(is.list(ord)){
 if('dbrda' %in% class(ord)){
  labels <- paste('RDA',1:2,sep = '-')
  s <- summary(ord)
  x <- s$cont$importance['Proportion Explained',1:2]
 } else if(class(ord) == 'list'){
  labels <- paste('Axis',1:2,sep = '-')
  x <- ord$eig[1:2] / sum(ord$eig[ord$eig>0])
  Z <- ord$points
 } else if (class(ord) == 'pcoa') {
  labels <- paste('Axis', 1:2, sep = '-')
  x <- ord$values$Eigenvalues[1:2] / sum(ord$values$Eigenvalues[ord$values$Eigenvalues>0])
  Z <- ord$vectors
 }
 labels <- paste(labels," (", as.character(round(100*x,1)),"%)", sep = '')
 }
 df <- data.frame(Z, group)
 names(df)[1:2] <- paste0('X',1:2)

 # Add centroids:
 group_by(df,group) %>%
 mutate(centroid_X  = mean(X1),
 centroid_Y  = mean(X2)) %>%
 arrange(group) -> df

 # Make visualization: ----
 ggplot(df,aes(X1,X2,fill = group))  +
 # geom_point(shape = 21,size = 3)  +
 coord_equal()  +
 xlab(labels[1])  +
 ylab(labels[2]) -> p

 if(ellipse){
 p <- p  +  stat_ellipse(level  = CI,
  geom = 'polygon',alpha = 0.1)
 }

 if(centroids){
 p <- p  +  geom_point(aes(x = centroid_X,
  y = centroid_Y,
  color = group), size = 3)
 }
 if(spyder){
 p <- p  +  geom_segment(data = df,
  aes(x = centroid_X, xend = X1, y = centroid_Y,
   yend = X2, color = group),
  show.legend  = FALSE, size = .1)
 }

 p <- p  +  theme_classic()  +
 scale_fill_manual(values = palette,name = '')  +
 scale_color_manual(values = palette,name = '')  +
 theme(legend.position  = 'top')

 return(p)
}

## function pairwise.adonis()
pairwise.adonis <- function(x,factors, sim.function  = 'vegdist', sim.method  = 'horn', p.adjust.m  = 'fdr')
{
library(vegan)


co  = combn(unique(as.character(factors)),2)
pairs  = c()
F.Model  = c()
R2  = c()
p.value  = c()
Df.factor = c()
Df.residual = c()


for(elem in 1:ncol(co)){
if(sim.function  ==  'daisy'){
library(cluster); x1  = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric = sim.method)
} else{x1  = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method = sim.method)}

ad  = adonis(x1  ~  factors[factors %in% c(co[1,elem],co[2,elem])] );
pairs  = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
F.Model  = c(F.Model,ad$aov.tab[1,4]);
R2  = c(R2,ad$aov.tab[1,5]);
p.value  = c(p.value,ad$aov.tab[1,6])
Df.factor    = c(Df.factor, ad$aov.tab[1,1])
Df.residual  = c(Df.residual, ad$aov.tab[2,1])

}
p.adjusted  = p.adjust(p.value,method = p.adjust.m)
sig  = c(rep('',length(p.adjusted)))
sig[p.adjusted <= 0.05] <- '.'
sig[p.adjusted <= 0.01] <- '*'
sig[p.adjusted <= 0.001] <- '**'
sig[p.adjusted <= 0.0001] <- '***'

pairw.res  = data.frame(pairs, Df.factor, Df.residual, F.Model, R2, p.value, p.adjusted, sig)
print("Signif. codes: 0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1")
return(pairw.res)

}

## 

#-----------------------
# Super dictionary
source("Dictionary.R")