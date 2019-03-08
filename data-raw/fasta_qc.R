library(ggplot2)

# get the header lines
fasta_file <- readLines("inst/extdata/NewYork_A_H3N2_20190307.fasta")
header_lines <- which(grepl("^>", fasta_file, fixed = F))
headers <- fasta_file[header_lines]
fields <- purrr::map(headers, ~as.data.frame(matrix(strsplit(gsub("_|>", " ", .x), "|", fixed = T)[[1]], nrow = 1)))
fields <- dplyr::bind_rows(fields)
colnames(fields) <- c("isolate_name", "collection_date", "last_modified",
                      "type", "lineage", "originating lab", "submitting lab")
fields$isolate_name <- gsub(" ", "", fields$isolate_name)
fields$collection_date <- gsub(" ", "", fields$collection_date)
fields$last_modified <- gsub(" ", "", fields$last_modified)
fields$type <- gsub(" ", "", fields$type)
fields$lineage <- gsub(" ", "", fields$lineage)
readr::write_csv(fields, path = "NewYork_A_H3N2_20190307.csv")

# remove samples with missing date information
line_ranges <- data.frame(id = seq_along(header_lines), start = header_lines)
line_ranges$end <- c(line_ranges$start[-1] - 1, length(fasta_file))
remove_isolates <- which(grepl("unknown", fields$collection_date)) # 134 are removed
remove_ranges <- line_ranges[remove_isolates,]
remove_ranges$remove_idx <- with(remove_ranges, purrr::map2(start, end, ~ seq(.x, .y))
fasta_filtered <- fasta_file[-unlist(remove_ranges$remove_idx)]
writeLines(fasta_filtered, "inst/extdata/NewYork_A_H3N2_20190307_filtered.fasta")
fields_filtered <- fields[-remove_isolates,]
readr::write_tsv(fields_filtered, path = "inst/extdata/NewYork_A_H3N2_20190307_filtered.tsv")

# plot histogram information
fields_filtered$collection_date <- as.POSIXct(strptime(fields_filtered$collection_date,
                                            format = "%Y-%m-%d"))
p <- ggplot(data = fields_filtered) +
  geom_histogram(aes(x = collection_date), bins = 100)
ggsave(p, file = "~/Desktop/H3N2_NewYork.png", width = 5, height = 4)
