parse.propensity.functions <- function(propensity.funs, x, parms) {
  names(propensity.funs) <- NULL

  x.names <- names(x)
  parms.names <- names(parms)

  strings <- sapply(seq_along(propensity.funs), function(i) {
    string <- propensity.funs[[i]]

    aag <- strsplit(gsub("([A-Za-z0-9_\\.][A-Za-z0-9_\\.]*)", " \\1 ", string), " ")[[1]]

    # check for x
    aag.match <- match(aag, x.names)
    ix <- !is.na(aag.match)
    aag[ix] <- paste0("x[[", aag.match[ix], "]]")

    # check for parms
    aag.match <- match(aag, parms.names)
    ix <- !is.na(aag.match)
    aag[ix] <- paste0("parms[[", aag.match[ix], "]]")

    # create line of code
    paste0(
      "  a[[", i, "]] <- ",
      paste(aag, collapse = ""),
      "\n"
    )
  })

  funtext <- paste0(
    "function(x, parms) {\n",
    "  a <- rep(NA, ", length(propensity.funs), ")\n",
    paste(strings, collapse = ""),
    "  a\n",
    "}\n"
  )

  eval(parse(text = funtext))
}
