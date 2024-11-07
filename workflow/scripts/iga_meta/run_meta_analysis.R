library(data.table)
setDTthreads(snakemake@threads)

chr_col <- snakemake@params[['chr_col']]
bp_col <- snakemake@params[['bp_col']]
ref_col <- snakemake@params[['ref_col']]
alt_col <- snakemake@params[['alt_col']]
beta_col <- snakemake@params[['beta_col']]
se_col <- snakemake@params[['se_col']]
p_col <- snakemake@params[['p_col']]

cols <- c(chr_col, bp_col, ref_col, alt_col, beta_col, se_col, p_col)

lyons <- fread(snakemake@input[['lyons']], sep = '\t', header = T, select = cols)
setnames(lyons, c('BETA', 'SE', 'P'), c('BETA.lyons', 'SE.lyons', 'P.lyons'))
lyons <- unique(lyons, by = c(chr_col, bp_col, ref_col, alt_col))

liu <- fread(snakemake@input[['liu']], sep = '\t', header = T, select = cols)
setnames(liu, c('BETA', 'SE', 'P'), c('BETA.liu', 'SE.liu', 'P.liu'))
liu <- unique(liu, by = c(chr_col, bp_col, ref_col, alt_col))

liu_decode <- fread(snakemake@input[['liu_decode']], sep = '\t', header = T, select = cols)
setnames(liu_decode, c('BETA', 'SE', 'P'), c('BETA.liu_decode', 'SE.liu_decode', 'P.liu_decode'))
liu_decode <- unique(liu_decode, by = c(chr_col, bp_col, ref_col, alt_col))

dennis <- fread(snakemake@input[['dennis']], sep = '\t', header = T, select = cols)
setnames(dennis, c('BETA', 'SE', 'P'), c('BETA.dennis', 'SE.dennis', 'P.dennis'))
dennis <- unique(dennis, by = c(chr_col, bp_col, ref_col, alt_col))

merged <- merge(liu, liu_decode, by = c(chr_col, bp_col, ref_col, alt_col), all = T)
merged <- merge(merged, lyons, by = c(chr_col, bp_col, ref_col, alt_col), all.x = T)
merged <- merge(merged, dennis, by = c(chr_col, bp_col, ref_col, alt_col), all.x = T)
merged <- unique(merged, by = c(chr_col, bp_col, ref_col, alt_col))

merged[, SNPID := paste(get(chr_col), get(bp_col), get(ref_col), get(alt_col), sep = ':')]

merged[, `:=` (wt.lyons = SE.lyons^-2, wt.liu = SE.liu^-2, wt.liu_decode = SE.liu_decode^-2, wt.dennis = SE.dennis^-2)]

merged[!is.na(BETA.lyons) & !is.na(BETA.liu), `:=` (BETA.liu_lyons = (BETA.liu * wt.liu + BETA.lyons * wt.lyons)/(wt.liu + wt.lyons), SE.liu_lyons = (wt.liu + wt.lyons)^-0.5)]
merged[is.na(BETA.lyons) & !is.na(BETA.liu), `:=` (BETA.liu_lyons = BETA.liu , SE.liu_lyons = SE.liu)]
merged[!is.na(BETA.lyons) & is.na(BETA.liu), `:=` (BETA.liu_lyons = BETA.lyons, SE.liu_lyons = SE.lyons)]
merged[!is.na(BETA.liu_lyons), Z2.liu_lyons := (BETA.liu_lyons/SE.liu_lyons)^2]
merged[!is.na(BETA.liu_lyons), P.liu_lyons := pchisq(Z2.liu_lyons, df = 1, lower.tail = F)]

fwrite(merged[!is.na(BETA.liu_lyons), .(CHR38 = get(chr_col), BP38 = get(bp_col), REF = get(ref_col), ALT = get(alt_col), SNPID, BETA = BETA.liu_lyons, SE = SE.liu_lyons, Z2 = Z2.liu_lyons, P = P.liu_lyons)], file = snakemake@output[['liu_lyons']], sep = '\t')

merged[!is.na(BETA.lyons) & !is.na(BETA.liu_decode), `:=` (BETA.liu_decode_lyons = (BETA.liu_decode * wt.liu_decode + BETA.lyons * wt.lyons)/(wt.liu_decode + wt.lyons), SE.liu_decode_lyons = (wt.liu_decode + wt.lyons)^-0.5)]
merged[is.na(BETA.lyons) & !is.na(BETA.liu_decode), `:=` (BETA.liu_decode_lyons = BETA.liu_decode , SE.liu_decode_lyons = SE.liu_decode)]
merged[!is.na(BETA.lyons) & is.na(BETA.liu_decode), `:=` (BETA.liu_decode_lyons = BETA.lyons , SE.liu_decode_lyons = SE.lyons)]
merged[!is.na(BETA.liu_decode_lyons), Z2.liu_decode_lyons := (BETA.liu_decode_lyons/SE.liu_decode_lyons)^2]
merged[!is.na(BETA.liu_decode_lyons), P.liu_decode_lyons := pchisq(Z2.liu_decode_lyons, df = 1, lower.tail = F)]

fwrite(merged[!is.na(BETA.liu_decode_lyons), .(CHR38 = get(chr_col), BP38 = get(bp_col), REF = get(ref_col), ALT = get(alt_col), SNPID, BETA = BETA.liu_decode_lyons, SE = SE.liu_decode_lyons, Z2 = Z2.liu_decode_lyons, P = P.liu_decode_lyons)], file = snakemake@output[['liu_decode_lyons']], sep = '\t')

merged[, wt.liu_decode_lyons := SE.liu_decode_lyons^-2]

merged[!is.na(BETA.liu_decode_lyons) & !is.na(BETA.dennis), `:=` (BETA.liu_decode_lyons_dennis = (BETA.liu_decode_lyons * wt.liu_decode_lyons + BETA.dennis * wt.dennis)/(wt.liu_decode_lyons + wt.dennis), SE.liu_decode_lyons_dennis = (wt.liu_decode_lyons + wt.dennis)^-0.5)]
merged[!is.na(BETA.liu_decode_lyons) & is.na(BETA.dennis), `:=` (BETA.liu_decode_lyons_dennis = BETA.liu_decode_lyons, SE.liu_decode_lyons_dennis = SE.liu_decode_lyons)]
merged[!is.na(BETA.liu_decode_lyons_dennis), Z2.liu_decode_lyons_dennis := (BETA.liu_decode_lyons_dennis/SE.liu_decode_lyons_dennis)^2]
merged[!is.na(BETA.liu_decode_lyons_dennis), P.liu_decode_lyons_dennis := pchisq(Z2.liu_decode_lyons_dennis, df = 1, lower.tail = F)]

fwrite(merged[!is.na(BETA.liu_decode_lyons_dennis), .(CHR38 = get(chr_col), BP38 = get(bp_col), REF = get(ref_col), ALT = get(alt_col), SNPID, BETA = BETA.liu_decode_lyons_dennis, SE = SE.liu_decode_lyons_dennis, Z2 = Z2.liu_decode_lyons_dennis, P = P.liu_decode_lyons_dennis)], file = snakemake@output[['liu_decode_lyons_dennis']], sep = '\t')

fwrite(merged[, .(CHR38 = get(chr_col),
                  BP38 = get(bp_col),
                  REF = get(ref_col),
                  ALT = get(alt_col),
                  SNPID,
                  BETA.meta = BETA.liu_decode_lyons_dennis, SE.meta = SE.liu_decode_lyons_dennis, P.meta = P.liu_decode_lyons_dennis,
                  BETA.liu_decode = BETA.liu_decode, SE.liu_decode = SE.liu_decode, P.liu_decode = P.liu_decode,
                  BETA.lyons = BETA.lyons, SE.lyons = SE.lyons, P.lyons = P.lyons,
                  BETA.dennis = BETA.dennis, SE.dennis = SE.dennis, P.dennis = P.dennis
                  )], file = snakemake@output[['merged']], sep = '\t')
