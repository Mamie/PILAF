library(ggplot2)
library(dplyr)
library(gridBase)
library(grid)

preprocess_fasta <- function(path) {
  fasta_file <- readLines(path)
  header_lines <- which(grepl("^>", fasta_file, fixed = F))
  headers <- fasta_file[header_lines]
  fields <- purrr::map(headers, ~as.data.frame(matrix(strsplit(gsub("_|>", " ", .x), "|", fixed = T)[[1]], nrow = 1)))
  fields <- dplyr::bind_rows(fields)
  colnames(fields) <- c("isolate_name", "collection_date", "last_modified",
                        "type", "lineage", "originating lab", "submitting lab")
  fields$isolate_name <- gsub(" ", "", fields$isolate_name)
  fields$collection_date <- as.POSIXct(strptime(gsub(" ", "", fields$collection_date),
                                                format = "%Y-%m-%d"))
  fields$last_modified <- as.POSIXct(strptime(gsub(" ", "", fields$last_modified),
                                              format = "%Y-%m-%d"))
  fields$type <- gsub(" ", "", fields$type)
  fields$lineage <- gsub(" ", "", fields$lineage)
  fields$id <- seq(nrow(fields))
  fields <- fields %>%
    arrange(isolate_name, desc(last_modified)) %>%
    mutate(duplicate = duplicated(isolate_name)) %>%
    arrange(id)
  return(list(fasta_file = fasta_file, fields = fields))
}

filter_date <- function(fasta, start_date, field_outpath, data_output) {
  fasta_file <- fasta$fasta_file
  header_lines <- which(grepl("^>", fasta_file, fixed = F))
  fields <- fasta$fields
  line_ranges <- data.frame(id = seq_along(header_lines), start = header_lines)
  line_ranges$end <- c(line_ranges$start[-1] - 1, length(fasta_file))
  remove_isolates <- which(is.na(fields$collection_date) | fields$duplicate | fields$collection_date < as.POSIXct(start_date))
  remove_ranges <- line_ranges[remove_isolates,]
  remove_ranges$remove_idx <- with(remove_ranges, purrr::map2(start, end, ~ seq(.x, .y)))
  fasta_filtered <- fasta_file[-unlist(remove_ranges$remove_idx)]
  fields_filtered <- fields[-remove_isolates,]
  readr::write_tsv(fields_filtered, path = field_outpath)
  filtered_headers_idx <- which(grepl(">", fasta_filtered))
  filtered_headers <- fasta_filtered[filtered_headers_idx]
  fasta_filtered[filtered_headers_idx] <- paste0(">", fields_filtered$isolate_name, "|", fields_filtered$collection_date)
  writeLines(fasta_filtered, data_output)
  fields_filtered$collection_date <- as.POSIXct(strptime(fields_filtered$collection_date,
                                                         format = "%Y-%m-%d"))
  return(list(fasta_file = fasta_filtered, fields = fields_filtered))
}

# H1N1
H1N1 <- preprocess_fasta("~/Downloads/20180318_NewYork_H1N1.fasta")
H1N1_filtered <- filter_date(H1N1, "2014-01-01", "inst/extdata/H1N1/NewYork_A_H1N1_20190318_filtered.tsv", "inst/extdata/H1N1/NewYork_A_H1N1_20190318_filtered.fasta") # 265
H1N1_last_time <- 2019.0931506849315

# B
Yamagata <- preprocess_fasta("~/Downloads/20190315_NewYork_B_Yamagata.fasta")
Yamagata_filtered <- filter_date(Yamagata, "2014-01-01", "inst/extdata/Yamagata/NewYork_B_Yamagata_20190318_filtered.tsv", "inst/extdata/Yamagata/NewYork_B_Yamagata_20190318_filtered.fasta") # 69
Yamagata_last_time <- 2019.0739726027398

# B
Victoria <- preprocess_fasta("~/Downloads/20190315_NewYork_B_Victoria.fasta")
Victoria_filtered <- filter_date(Victoria, "2014-01-01", "inst/extdata/Victoria/NewYork_B_Victoria_20190318_filtered.tsv", "inst/extdata/Victoria/NewYork_B_Victoria_20190318_filtered.fasta") # 52
Victoria_last_time <- 2018.317808219178

# get the header lines
fasta_file <- readLines("inst/extdata/NewYork_A_H3N2_20190307.fasta")
header_lines <- which(grepl("^>", fasta_file, fixed = F))
headers <- fasta_file[header_lines]
fields <- purrr::map(headers, ~as.data.frame(matrix(strsplit(gsub("_|>", " ", .x), "|", fixed = T)[[1]], nrow = 1)))
fields <- dplyr::bind_rows(fields)
colnames(fields) <- c("isolate_name", "collection_date", "last_modified",
                      "type", "lineage", "originating lab", "submitting lab")
fields$isolate_name <- gsub(" ", "", fields$isolate_name)
fields$collection_date <- as.POSIXct(strptime(gsub(" ", "", fields$collection_date),
                                              format = "%Y-%m-%d"))
fields$last_modified <- as.POSIXct(strptime(gsub(" ", "", fields$last_modified),
                                            format = "%Y-%m-%d"))
fields$type <- gsub(" ", "", fields$type)
fields$lineage <- gsub(" ", "", fields$lineage)
fields$id <- seq(nrow(fields))
# readr::write_csv(fields, path = "NewYork_A_H3N2_20190307.csv")
fields <- fields %>%
  arrange(isolate_name, desc(last_modified)) %>%
  mutate(duplicate = duplicated(isolate_name)) %>%
  arrange(id)

# remove samples with missing date information
line_ranges <- data.frame(id = seq_along(header_lines), start = header_lines)
line_ranges$end <- c(line_ranges$start[-1] - 1, length(fasta_file))
remove_isolates <- which(is.na(fields$collection_date) | fields$duplicate | fields$collection_date < as.POSIXct("2014-01-01"))
remove_ranges <- line_ranges[remove_isolates,]
remove_ranges$remove_idx <- with(remove_ranges, purrr::map2(start, end, ~ seq(.x, .y)))
fasta_filtered <- fasta_file[-unlist(remove_ranges$remove_idx)]
fields_filtered <- fields[-remove_isolates,]
readr::write_tsv(fields_filtered, path = "inst/extdata/NewYork_A_H3N2_20190307_filtered.tsv")
filtered_headers_idx <- which(grepl(">", fasta_filtered))
filtered_headers <- fasta_filtered[filtered_headers_idx]
fasta_filtered[filtered_headers_idx] <- paste0(">", fields_filtered$isolate_name, "|", fields_filtered$collection_date)
writeLines(fasta_filtered, "inst/extdata/NewYork_A_H3N2_20190307_filtered.fasta")

# plot histogram information
fields_filtered$collection_date <- as.POSIXct(strptime(fields_filtered$collection_date,
                                            format = "%Y-%m-%d"))
p <- ggplot(data = fields_filtered) +
  geom_histogram(aes(x = collection_date), bins = 100)
ggsave(p, file = "~/Desktop/H3N2_NewYork.png", width = 5, height = 4)

# read the tree after running the BEAST
H3N2_MCC_tree <- ape::read.nexus("data-raw/20190308_H3N2_SkyGrid/MCC.tree")
H3N2_MCC_phylo <- phylodyn::summarize_phylo(H3N2_MCC_tree)
last_time <- 2019.0739726027398
root_time = last_time - max(H3N2_MCC_phylo$coal_times)
estimated_root_height <- last_time - 6.687

# read output from skygrid reconstruction

skygrid <- readr::read_tsv("data-raw/20190308_H3N2_SkyGrid/skygrid_reconstruction.tsv", skip = 1) %>%
  mutate(Time = last_time - Time)

pdf(file = "~/Desktop/H3N2.pdf", width = 5, height = 6)
plot_MCC_Ne(H3N2_MCC_tree, root_time, last_time, skygrid, breaks = seq(2013, 2019))
dev.off()

# ILI index
ny_data <- cdcfluview::ilinet(region = "state") %>%
  filter(region == 'New York')
cdc_data <- select(ny_data, year, week, ilitotal)
cdc_data$Time <-
