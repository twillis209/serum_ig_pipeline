library(data.table)

lead <- fread(snakemake@input$lead)

novel <- fread(snakemake@input$novel)

lead[, Novel := FALSE]

lead[novel, on = .(chromosome, base_pair_location), Novel := TRUE]

fwrite(lead, file = snakemake@output[[1]], sep = '\t')
