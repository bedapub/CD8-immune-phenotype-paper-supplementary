# Tumor-agnostic transcriptome-based classifier identifies spatial infiltration patterns of CD8+ T cells in the tumor microenvironment and predicts clinical # outcome in early- and late-phase clinical trials

This is a repository for accompanying analysis code.

The trained model can predict CD8 immune phenotype from RNA-Seq counts data
can be found [here](https://github.com/bedapub/cd8ippred).

To rerun the analysis, you should generally follow `.Rmd` files in the
`analysis/` directory in ascending order.

It is expected that data is located in an S3-compatible storage; replace
`"UUID"` with actual UUIDs in the YAML headers of the Rmarkdowns.
Alternatively, you could replace [`aws.s3`](https://github.com/cloudyr/aws.s3)
calls with respective file reading code. E.g., the following code:

``` r
eset <- aws.s3::s3readRDS(
  "eset_batch_corrected.Rds",
  bucket = params$output_collection
)

pheno_data <- aws.s3::s3read_using(
  read_tsv,
  col_types = cols(),
  object = "samples_anno.tsv",
  bucket = params$data_collection
)
```

will become:

``` r
eset <- readRDS(here::here("output/eset_batch_corrected.Rds"))

pheno_data <- read_tsv(
  here::here("data/samples_anno.tsv"),
  col_types = cols()
)
```

