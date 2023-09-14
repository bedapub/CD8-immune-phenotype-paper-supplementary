# load Rdata object without polluting the global environment
my_load <- function(fname) {
  e <- new.env()
  load(fname, e)
  as.list(e)
}

# apply make.names to eset sample (column) names
make_names_eset <- function(eset) {
  colnames(eset) <- make.names(colnames(eset))
  eset
}

pheno_col <- c(
  inflamed = "#d95f02",
  excluded = "#1b9e77",
  desert = "#7570b3"
)

atezo_study_names <- c(
  "GO29293" = "IMvigor 210",
  "GO28915" = "OAK",
  "GO29294" = "IMvigor 211"
)
