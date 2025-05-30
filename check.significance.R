    row.names(DF) <- DF$Row.names

    # Linear regression for ARG load and metadata
    result <-
      do.call(rbind, lapply(vars, function(x) {
      
        #fit <-
        #lm(ARGload_log10 ~ DF[, x]+
        #   PREVAL_RX_J01_NEVT + PREVAL_RX_J01A_NEVT +
        #   PREVAL_RX_J01A + PREVAL_RX_J01F + BL_USE_RX_J01, data = DF)

        # Ignore covariates for the region/gender specific models due to low power;
	# control these in full data models
        fit <- lm(ARGload_log10 ~ DF[, x], data = DF)

      lm_with_cov.result <- fit %>% summary()
      CI <- fit %>% confint()
      CI_low <- 10^CI[2, 1]
      CI_up  <- 10^CI[2, 2]
      cov.estimate <- lm_with_cov.result$coefficients[2, 1]
      cov.pvalue  <- lm_with_cov.result$coefficients[2, 4]
      return(
        data.frame(
          variable = x,
          cov.pvalue = cov.pvalue,
          cov.estimate = 10^cov.estimate,
          CI2.5 = CI_low,
          CI97.5 = CI_up
        )
      )
    }))
    rownames(result) <- vars
    result0 <- result

    result <-
      result %>% t() %>%
      data.frame() %>%
      dplyr::select(-contains(c("sum", "SUM_norm", "scaled", "ARG", "_KAIKKI", "GAMMA", "RX_J"))) %>%
      t() %>%
      data.frame()
      
    result$padj <-
      p.adjust(result$cov.pvalue, method = "fdr")

    result$VARIABLE <- rownames(result)
    result <- result %>% mutate_if(is.numeric, signif, digits=3)

    df <-
      data.frame(VARIABLE = FR02names$VARIABLE, LONGNAME = FR02names$LONGNAME)

    result <- dplyr::right_join(df, result, by = join_by(VARIABLE))

    sigs <- result# [result$padj<0.05,] 
    df <- sigs %>% replace_variable_with_dictionary("VARIABLE", substitutions) %>% 
      mutate(VARIABLE=factor(VARIABLE, levels=VARIABLE)) 

    df$gender <- rep(men, nrow(df))
