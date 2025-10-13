gc_count <- function(formula, data, group, effect="ATE", model, param.tune=NULL, cv=10, boot.type="bcv",
                          boot.number=500, boot.tune=FALSE, boot.mi = FALSE, progress=TRUE, seed=NULL, m=5, ...) {
  
  if (!is.logical(boot.mi)) {stop("The argument \"boot.mi\" needs to be a logical TRUE or FALSE")}
  if (boot.mi == FALSE) {
    return(.gc_count(formula=formula, data=data, group=group, effect=effect, model=model, param.tune=param.tune, cv=cv,
                          boot.type=boot.type,boot.number=boot.number, boot.tune=boot.tune, progress=progress ,seed=seed))
  }
  if (boot.mi == TRUE) {
    return(.miboot_gc_count(formula=formula, data=data, group=group, effect=effect, model=model, param.tune=param.tune, cv=cv,
                                 boot.type=boot.type,boot.number=boot.number, boot.tune=boot.tune, progress=progress, seed=seed,
                                 m=m, ...))
  }
}
