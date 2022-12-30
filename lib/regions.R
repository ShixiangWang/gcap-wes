# library(modules)
import(magrittr)
import(cli)

get_intersect_size <- function(x.start, x.end, y.start, y.end, use_furrr = FALSE) {
  if (use_furrr) {
    future::plan("multisession", workers = 20) # future::availableCores() - 1L)
    r <- furrr::future_pmap(list(x.start = x.start, x.end = x.end, y.start = y.start, y.end = y.end), function(x.start, x.end,
                                                                                                               y.start, y.end) {
      IRanges::width(IRanges::intersect(IRanges::IRanges(x.start, x.end), IRanges::IRanges(y.start, y.end)))
    }, .progress = TRUE)
  } else {
    r <- sigminer:::get_intersect_size(x.start, x.end, y.start, y.end)
  }
  unlist(r)
}

# x = c('raw/CircleMap/tcga/sample1_realign.bed', 'raw/CircleMap/tcga/sample1_repeats.bed')
collapse_to_genes <- function(x, ref_file = "data/hg38_genes.bed", cores = 4, use_furrr = FALSE, save = FALSE) {
  if (!is.data.frame(x)) {
    if (!all(file.exists(x))) {
      cli_abort("Can't find file{?s}: {.file {x[!file.exists(x)]}}")
    }
  }
  if (!file.exists(ref_file)) {
    cli_abort("Can't find file {.file {ref_file}}!")
  }

  if (!is.data.frame(x) && length(x) > 1) {
    cli_alert_info("Multiple input files detected.")
    out <- parallel::mclapply(x, collapse_to_genes, ref_file = ref_file, save = save, mc.cores = cores)
    names(out) <- x
  } else {
    if (is.data.frame(x)) {
      cli_alert_info("A data.frame detected. Please make sure the first 3 columns are for chr, start, end.")
      x <- data.table::as.data.table(x)
    } else {
      cli_alert_info("Processing {x} ...")
      x2 <- x
      x <- data.table::fread(x, header = FALSE)
    }

    cli_alert("reading reference file")
    y <- data.table::fread(ref_file, header = FALSE)
    colnames(x)[1:3] <- colnames(y)[1:3] <- c("chr", "start", "end")
    colnames(y)[4] <- "gene_id"

    cli_alert("finding overlaps")
    # Determine the intersect size
    data.table::setkey(y, chr, start, end)
    out <- data.table::foverlaps(x, y)[!is.na(gene_id)]

    cli_alert("calculating intersect size")
    out[, `:=`(intersect_size, get_intersect_size(i.start, i.end, start, end, use_furrr = use_furrr))]
    # Calculate the region cov ratio
    out[, `:=`(intersect_ratio, intersect_size / (abs(end - start) + 1))]

    if (save) {
      saveRDS(out, file = paste0(basename(x2), ".rds"))
    }
  }
  cli_alert("done")
  if (save) {
    return(invisible())
  } else {
    out
  }
}

# x <- c('raw/CircleMap/tcga/sample1_realign.bed', 'raw/CircleMap/tcga/sample1_repeats.bed') out <- collapse_to_genes(x)

overlaps <- function(x, y, use_furrr = FALSE) {
  ## Overlaps genome regions from x and y @param x region data to explore.  @param y region data to be overlapped.  @param
  ## use_furrr if `TRUE`, use `furrr` package to parallel calculation of intersection size.  Note, set `use_furrr = FALSE` is
  ## much faster now.  @return a `data.table`

  cli_alert_info("x and y must have chr, start, end at the first 3 columns")

  if (!is.data.frame(x)) {
    cli_abort("x must be a {.code data.frame}")
  }
  if (!is.data.frame(y)) {
    cli_abort("y must be a {.code data.frame}")
  }

  x <- data.table::as.data.table(x)
  y <- data.table::as.data.table(y)

  colnames(x)[1:3] <- colnames(y)[1:3] <- c("chr", "start", "end")

  # Determine the intersect size
  data.table::setkey(y, chr, start, end)
  out <- data.table::foverlaps(x, y)[!is.na(start)]

  cli_alert("calculating intersect size and ratio")
  out[, `:=`(intersect_size, get_intersect_size(i.start, i.end, start, end, use_furrr = use_furrr))]
  # Calculate the region cov ratio
  out[, `:=`(intersect_ratio, intersect_size / (abs(end - start) + 1))]

  cli_alert("done")
  out
}
