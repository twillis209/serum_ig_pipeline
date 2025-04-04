library(data.table)
setDTthreads(snakemake@threads)

flank <- snakemake@params$window/2

iso <- snakemake@wildcards$isotype
non_ig <- snakemake@wildcards$non_ig

ig <- fread(snakemake@input$ig)
ig[, is_candidate := FALSE]
ig[, non_ig_lead_rsid := '']
ig[, non_ig_min_p := 1.0]

merged <- fread(snakemake@input$merged)
merged[, names(.SD) := lapply(.SD, as.numeric), .SDcols = patterns('^p_value')]

for(i in seq_len(nrow(ig))) {
  non_ig_neighbours <- merged[
    chromosome == ig[i, chromosome] &
    base_pair_location %between% c(max(ig[i, base_pair_location] - flank, 0),
                                   ig[i, base_pair_location] + flank)
  ]

  if(non_ig_neighbours[, .N] > 0) {
    if(non_ig_neighbours[, min(p.non_ig, na.rm = TRUE),
                      env = list(p.non_ig = sprintf("p_value.%s", non_ig))] < 5e-8) {
    non_ig_lead_rsid <- non_ig_neighbours[which.min(p.non_ig), rsid,
                                          env = list(p.non_ig = sprintf("p_value.%s", non_ig))]
    non_ig_min_p <- non_ig_neighbours[, min(p.non_ig, na.rm = T),
                      env = list(p.non_ig = sprintf("p_value.%s", non_ig))]
    set(ig, i, 'is_candidate', TRUE)
    set(ig, i, 'non_ig_lead_rsid', non_ig_lead_rsid)
    set(ig, i, 'non_ig_min_p', non_ig_min_p)
    }
  }
}


fwrite(ig[is_candidate == TRUE], file = snakemake@output[[1]], sep = '\t')
