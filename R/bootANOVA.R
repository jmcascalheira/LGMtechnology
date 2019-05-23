-------------------------------------------------------------------------
  # Project: Bootstrap AN(C)OVA with Posthoc Pairwise Comparisons
  # Coded by: Sam Mancuso
  # Mail: sammancuso.wordpress.com
  # Date: 27 May 2015
  # Version: 2.0
  # Revision Date: 2 June 2015
  # -------------------------------------------------------------------------


# Bootstrapped AN(C)OVA Options -------------------------------------------
#
# formula:      Formula for AN(C)OVA model
#
# conf.int:     Confidence intervals (default = 0.95)
#
# dec:          Number of decimal places (default = 3)
#
# reps:         Number of bootstrap replications (default = 1000)
#
# pw.comp:      List of variables for pairwise comparisons (default = NULL)
#
# seed:         Seed to replicate results (default = 1234)

# Bootstrap statistic -----------------------------------------------------
anboot <- function(formula, data, pw.comp, i){

  dataResamp <- data[i,] # Resample rows

  lmFit <- lm(formula, dataResamp)

  FCoef <- Anova(lmFit, type = 3)[["F value"]] # Extract F-values
  FCoef <- FCoef[1:length(FCoef) - 1] # Removes NA entry

  if(is.null(pw.comp) == TRUE) {

    return <- c(FCoef)

    # If pairwise comparisons have been specified
  } else {

    postHocCoeff <- NULL # To prevent errors

    for(i in 1:length(pw.comp)) {

      # Parses expression to be passed to lsmeans()
      exprString <- paste("pairwise ~", pw.comp[i])
      exprEval <- eval(parse(text = exprString))

      suppressMessages(lsmFit <- lsmeans(lmFit, exprEval, data = dataResamp))

      suppressMessages(postHocTemp<- summary(lsmFit))

      # Get comparison coeffecient(s)
      postHocCoeff <- c(postHocCoeff, as.numeric(postHocTemp$contrasts[["estimate"]]))

    }

    return <- c(FCoef, postHocCoeff)
  }
}

# Function to run bootstrap and provide output ----------------------------
boot.anova <- function(data, formula, conf.int = 0.95, dec = 2, reps = 1000,
                       pw.comp = NULL, seed = 1234) {
  # Set seed (for replication)
  set.seed(seed)

  # Fit non-bootstrapped ANOVA to obtain Sum Sq and df
  lmFit <- lm(formula, data = data)
  # Get DV name
  dvName <- colnames(model.frame(lmFit))[1]
  # Get IV names
  ivNames <- attr(lmFit$terms, "term.labels")

  # Fit ANOVA model with Type III Sum of Squares
  anovaFit <- Anova(lmFit, type = 3)

  # Extract F
  Fssq <- anovaFit[["Sum Sq"]]           # Sum of Squares (Type III)
  FValues <- anovaFit[["F value"]]       # F-values
  Fdf <- anovaFit[["Df"]]                # df

  # Additional information for use later
  FLength <- length(FValues) - 1         # Number of F-statistics - Residual

  # Pairwise comparisons (if requested)
  if(is.null(pw.comp) == FALSE) {
    # Initialise variables to prevent errors
    postHocNames <- NULL
    postHocCoeff <- NULL
    postHocComparisons <- NULL
    postHocDF <- NULL
    CompName <- NULL

    for(i in 1:length(pw.comp)) {
      # Parses expression to be passed to lsmeans()
      exprString <- paste("pairwise ~", pw.comp[i])
      exprEval <- eval(parse(text = exprString))
      suppressMessages(postHocRes <- summary(lsmeans(lmFit, exprEval, data = data)))

      # If interaction term specified
      if(grepl("|", pw.comp[i], fixed = TRUE) == FALSE) {
        CompName <- as.character(postHocRes$contrasts[, 1])
      } else {
        # Comparison name: "A: B"
        CompName <- paste(as.character(postHocRes$contrasts[, 1]),
                          as.character(postHocRes$contrasts[, 2]),
                          sep = ": ")
      }

      # Comparison(s) name(s)
      postHocNames <- c(postHocNames, CompName)

      # Get comparison coeffecient(s)
      postHocCoeff <- c(postHocCoeff, as.numeric(postHocRes$contrasts[["estimate"]]))

      # Get number of comparisons (to adjust confidence interval for
      # multiple comparisons)
      postHocComparisons <- c(postHocComparisons,
                              length(as.numeric(postHocRes$contrasts[["estimate"]])))

      # Get df
      postHocDF <- c(postHocDF, as.numeric(postHocRes$contrasts[["df"]]))
    }
  }

  # Call anboot function
  anbootRes <- boot(statistic = anboot, formula = formula,
                    data = data, R = reps, pw.comp = pw.comp)

  # Get F coefficients (Fvalue from above)
  FSE <- summary(anbootRes)$bootSE[1:FLength]     # Bootstrap SE
  FZval <- FValues[1:FLength]/FSE                 # z-values
  FZp <- 2*pnorm(-abs(FZval))                     # Pr(>|z|)

  # Critical F-Value (excludes the 'Residuals' term)
  Fcrit <- NULL # To prevent errors
  ResidDf <- Fdf[length(Fdf)] # Get Residual df (last df)

  for(i in 1:FLength) {
    Fcrit <- c(Fcrit, qf(conf.int, Fdf[i], ResidDf))
  }

  # Get Bootstrapped confidence intervals
  FCILower <- NULL
  FCIUpper <- NULL

  for(i in 1:FLength){
    FbootBca <- boot.ci(anbootRes, index = i, type = "bca", conf = conf.int)$bca
    FCILower <- c(FCILower, FbootBca[4])
    FCIUpper <- c(FCIUpper, FbootBca[5])
  }

  # p-value stars
  FZstar <- NULL

  for(i in 1:length(FZp)) {
    if(FZp[i] < .001) {
      FZstar[i] <- "***"
    } else if(FZp[i] < .01) {
      FZstar[i] <- "**"
    } else if(FZp[i] < .05) {
      FZstar[i] <- "*"
    } else if(FZp[i] < .10){
      FZstar[i] <- "."
    } else {
      FZstar[i] <- ""
    }
  }


  # Create a Data Table for Anova Table

  # Column names
  cnames <- c("Sum Sq", "Df", "F value", "SE", "LB", "UB", "Crit F", "z", "Pr(>|z|)", "")

  # Row names
  rnames <- c("(Intercept)", ivNames, "Residuals")

  # Turn warnings off as the variables are not equal length when creating output.df
  options(warn = -1)

  # Create data frame to use as output
  outputAnova <- as.data.frame(cbind(round(Fssq, dec),
                                     Fdf,
                                     format(round(FValues, dec), dec),
                                     format(round(FSE, dec), dec),
                                     format(round(FCILower, dec), dec),
                                     format(round(FCIUpper, dec), dec),
                                     format(round(Fcrit, dec), dec),
                                     format(round(FZval, dec), dec),
                                     signif(FZp, 3),
                                     FZstar))

  # Remove redundant values for the Residuals row
  outputAnova[nrow(outputAnova), 3:10] <- NA

  # Set column and row names
  colnames(outputAnova) <- cnames
  rownames(outputAnova) <- rnames

  # Print results
  cat("BCa Bootstrap Anova Table (Type III tests)",
      "\n",
      "\n", "Response: ", dvName,
      "\n",
      sep ="")

  print(outputAnova, justify = "right", na.print = "" , quote = FALSE )

  cat("---",
      "\n", "Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1",
      sep = "")

  # Turn warnings on
  options(warn = 0)

  # Multiple Comparisons if requested
  if(is.null(pw.comp) == FALSE) {
    # Get Multiple comparison information
    BootSE <- summary(anbootRes)$bootSE

    MDiffStart <- FLength + 1
    MDiffLength <- length(BootSE)

    MDiffSE <- BootSE[MDiffStart:MDiffLength]

    # Calculate t value
    MDifft <- abs(postHocCoeff/MDiffSE)
    MDifftp <- 2*pt(-abs(MDifft),df = postHocDF)

    # Create a vector of confidence intervals
    confIntAlpha <- NULL

    for(i in 1:length(postHocComparisons)){
      confIntAlpha <- c(confIntAlpha,
                        rep((1 - conf.int)/postHocComparisons[i], length = postHocComparisons[i]))
    }

    confIntList <- 1 - confIntAlpha

    # Contrasts star
    MDiffStar <- NULL

    for(i in 1:length(MDifftp)) {
      if(MDifftp[i] < confIntAlpha[i]) {
        MDiffStar[i] <- "*"
      } else {
        MDiffStar[i] = ""
      }
    }

    MDiffCILower <- NULL
    MDiffCIUpper <- NULL

    # If single comparison
    if(MDiffStart == MDiffLength) {
      MDiffbootBca <- boot.ci(anbootRes, index = i, type = "bca", conf = confIntList[1])$bca
      MDiffCILower <- MDiffbootBca[4]
      MDiffCIUpper <- MDiffbootBca[5]
    } else {
      j <- 1

      for(i in MDiffStart:MDiffLength) {
        MDiffbootBca <- boot.ci(anbootRes, index = i, type = "bca", conf = confIntList[j])$bca
        MDiffCILower <- c(MDiffCILower, MDiffbootBca[4])
        MDiffCIUpper <- c(MDiffCIUpper, MDiffbootBca[5])
        j <- j + 1
      }
    }

    # Create data table for pairwise comparisons
    # Turn warnings off as the variables are not equal length when creating output.df
    options(warn = -1)

    # Create data frame to use as output
    cat("\n",
        "\n",
        "Post hoc Pairwise Contrasts with Bonferroni Adjustment",
        "\n",
        sep ="")

    # Output each specified contrast
    startVal <- 0
    endVal <- 0

    for(i in 1:length(pw.comp)) {
      cat("\n", "Contrast: ", pw.comp[i],
          "\n", sep = "")

      if(i == 1) {
        startVal <- 1
      } else {
        startVal <- endVal + 1
      }

      endVal <- endVal + postHocComparisons[i]
      outputPH <- as.data.frame(cbind(round(postHocCoeff[startVal:endVal], dec),
                                      format(round(MDiffSE[startVal:endVal], dec), dec),
                                      format(round(MDiffCILower[startVal:endVal], dec), dec),
                                      format(round(MDiffCIUpper[startVal:endVal], dec), dec),
                                      format(round(confIntList[startVal:endVal]*100, dec), dec),
                                      format(round(MDifft[startVal:endVal], dec), dec),
                                      postHocDF[startVal:endVal],
                                      signif(MDifftp[startVal:endVal], 3),
                                      signif(confIntAlpha[startVal:endVal], 3),
                                      MDiffStar[startVal:endVal]))
      # Column Names
      cnames <- c("Estimate", "SE", "LB", "UB", "% CI", "t value", "Df", "Pr(>|t|)", "Adj. Alpha", "")

      # Row names
      rnames <- postHocNames[startVal:endVal]

      colnames(outputPH) <- cnames
      rownames(outputPH) <- rnames

      # Turn warnings on
      options(warn = 0)

      print(outputPH, justify = "left", na.print = "" , quote = FALSE )

      cat("---",
          "\n", "Signif. codes:  ‘*’ Adj. p ‘ ’ 1",
          "\n",
          sep = "")
    }
  }
}
