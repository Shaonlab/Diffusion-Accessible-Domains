#!/usr/bin/Rscript

# Use "Rscript script.r -h" to see the help message and the options associated with the script 

suppressPackageStartupMessages(library(dplyr, quietly = TRUE))
suppressPackageStartupMessages(library(stats, quietly = TRUE))
suppressPackageStartupMessages(library(bedr, quietly = TRUE))
suppressPackageStartupMessages(library(tidyr, quietly = TRUE))
suppressPackageStartupMessages(library(data.table, quietly = TRUE))
suppressPackageStartupMessages(library(ggplot2, quietly = TRUE))
suppressPackageStartupMessages(library(mclust, quietly = TRUE))
suppressPackageStartupMessages(library(argparser, quietly = TRUE))

options(dplyr.summarise.inform = FALSE)

# Create a parser
p <- arg_parser("*** Calculating difference in epigenetic signal across the genome using ChIP-Seq, Replication-Seq and HiC-Seq data ***")

# Add command line arguments
p <- add_argument(p, "-m", "--modification", help = "Name of the epigenetic modification")
p <- add_argument(p, "-e", "--early-replication-file", help = "Name of Early-replication timing file (.bed)")
p <- add_argument(p, "-l", "--late-replication-file", help = "Name of Late-replication timing file (.bed)")
p <- add_argument(p, "-c", "--correlation-matrices", help = "Complete path to directory of .txt files for HiC data (inter & intra chromosomal)")
p <- add_argument(p, "-a", "--early-median", help = "Median value for early fragments")
p <- add_argument(p, "-b", "--late-median", help = "Median value for late fragments")
p <- add_argument(p, "-o", "--out", help = "Complete path to directory to save generated plots")

# Parse the command line arguments
argv <- parse_args(p)

# Before running this script, please check whether the current wd has the files
mod <- argv$m
hist_HiC_df <- bedr(engine = "bedtools", params = c("-wo"), input = list(a = "HiC_all_chr_loci.bed", b = paste0("GM12878_", mod, "_hg19.bed")), method = "intersect")
hist_HiC_df$V11 <- as.numeric(hist_HiC_df$V11)
hist_HiC_df$V17 <- as.numeric(hist_HiC_df$V17)
Signal <- hist_HiC_df$V11 * hist_HiC_df$V17 / 1000000
hist_HiC_df <- mutate(hist_HiC_df, Signal)
ab <- group_by(hist_HiC_df, chr, start, end)
signal_sum <- summarize(ab, sum_Signal = sum(Signal))
signal_sum_df <- as.data.frame(signal_sum)
chr_start <- as.numeric(signal_sum_df$start)
chr_start <- format(chr_start + 1, scientific = FALSE)
signal_sum_df <- format(mutate(signal_sum_df, chr_start), scientific = FALSE)
drops <- c("start")
signal_sum_df <- signal_sum_df[, !(names(signal_sum_df) %in% drops)]
signal_sum_df <- signal_sum_df %>%
    rename(
        start = chr_start
    )
signal_sum_df <- signal_sum_df[c("chr", "start", "end", "sum_Signal")]

rep_early <- read.csv(argv$e, header = FALSE, sep = "\t")
ba <- group_by(rep_early, V1, V2, V3)
frag_sum_early <- summarize(ba, sum_frag = sum(V10))
rep_early_df <- as.data.frame(frag_sum_early)
chr_start <- as.numeric(rep_early_df$V2)
chr_start <- format(chr_start + 1, scientific = FALSE)
rep_early_df <- format(mutate(rep_early_df, chr_start), scientific = FALSE)
drops <- c("V2")
rep_early_df <- rep_early_df[, !(names(rep_early_df) %in% drops)]
rep_early_df <- rep_early_df %>%
    rename(
        start = chr_start
    )
rep_early_df <- rep_early_df[c("V1", "start", "V3", "sum_frag")]
rep_early_df$combine <- as.character(interaction(rep_early_df$V1, rep_early_df$start, rep_early_df$V3))


# Late replicating fragments
rep_late <- read.csv(argv$l, header = FALSE, sep = "\t")
dc <- group_by(rep_late, V1, V2, V3)
frag_sum_late <- summarize(dc, sum_frag = sum(V10))
rep_late_df <- as.data.frame(frag_sum_late)
chr_start <- as.numeric(rep_late_df$V2)
chr_start <- format(chr_start + 1, scientific = FALSE)
rep_late_df <- format(mutate(rep_late_df, chr_start), scientific = FALSE)
drops <- c("V2")
rep_late_df <- rep_late_df[, !(names(rep_late_df) %in% drops)]
rep_late_df <- rep_late_df %>%
    rename(
        start = chr_start
    )
rep_late_df <- rep_late_df[c("V1", "start", "V3", "sum_frag")]
rep_late_df$combine <- as.character(interaction(rep_late_df$V1, rep_late_df$start, rep_late_df$V3))

# new method to find unique fragments (mixture-model)
# to check overlap between the 2 replication seq datasets
early_frag <- semi_join(rep_early_df, rep_late_df, by = "combine")
late_frag <- semi_join(rep_late_df, rep_early_df, by = "combine")
# to find unqiue fragments in the 2 replication seq datasets
unique_early_frag <- anti_join(rep_early_df, rep_late_df, by = "combine")
unique_early_frag <- unique_early_frag[c(1:4)]
unique_late_frag <- anti_join(rep_late_df, rep_early_df, by = "combine")
unique_late_frag <- unique_late_frag[c(1:4)]

df.frag <- cbind(early_frag, late_frag)
df.frag$difference <- as.numeric(early_frag$sum_frag) - as.numeric(late_frag$sum_frag)
df.frag <- subset(df.frag, df.frag$difference < -1000 | df.frag$difference > 1000)

rep_early_df <- subset(df.frag, df.frag$difference > 0)
rep_early_df <- rep_early_df[c(1:4)]
rep_early_df <- rbind(rep_early_df, unique_early_frag)
rep_late_df <- subset(df.frag, df.frag$difference < 0)
rep_late_df <- rep_late_df[c(1:4)]
rep_late_df <- rbind(rep_late_df, unique_late_frag)

# Calculating mean signal for Epigenetic signal normalization
colnames(rep_early_df) <- c("chr", "start", "end", "sum_frag")
colnames(rep_late_df) <- c("chr", "start", "end", "sum_frag")

early_fragments <- semi_join(signal_sum_df, rep_early_df, by = c("chr", "start", "end"))
mean_early <- mean(as.numeric(early_fragments$sum_Signal))
late_fragments <- semi_join(signal_sum_df, rep_late_df, by = c("chr", "start", "end"))
mean_late <- mean(as.numeric(late_fragments$sum_Signal))

colnames(rep_early_df) <- c("V1", "start", "V3", "sum_frag")
colnames(rep_late_df) <- c("V1", "start", "V3", "sum_frag")

all_Dij_Early <- c()
all_HiC_Early <- c()
all_Dij_Late <- c()
all_HiC_Late <- c()
setwd(argv$c)
cat("\n")
cat("Processing data... \n")
for (p in 1:22) {
  for (q in (p):22) {

    # Chromosome A- histone modification dataset
    chrNo <- paste("chr", p, sep = "")
    # get subset of chrA histone modification data
    hist_mod_chromA <- subset(signal_sum_df, signal_sum_df$chr == chrNo)

    # Replication timing dataset for Chromosome A
    # Early Replication
    early_rep_chromA <- subset(rep_early_df, rep_early_df$V1 == chrNo)
    early_rep_hist_chromA <- inner_join(early_rep_chromA, hist_mod_chromA, by = "start")
    drops <- c("chr", "end")
    early_rep_hist_chromA <- early_rep_hist_chromA[, !(names(early_rep_hist_chromA) %in% drops)]
    # convert sum_signal in hist_mod_chromA into a numeric form for use in calculating D(i.j) matrix
    # later
    early_rep_hist_chromA$sum_Signal <- as.numeric(early_rep_hist_chromA$sum_Signal)
    early_rep_hist_chromA$sum_Signal <- early_rep_hist_chromA$sum_Signal / mean_early
    chr_A_early <- paste0(early_rep_hist_chromA$V1, ":", early_rep_hist_chromA$start, "-", early_rep_hist_chromA$V3)
    chr_A_early <- gsub(" ", "", chr_A_early)

    # Late Replication
    late_rep_chromA <- subset(rep_late_df, rep_late_df$V1 == chrNo)
    late_rep_hist_chromA <- inner_join(late_rep_chromA, hist_mod_chromA, by = "start")
    drops <- c("chr", "end")
    late_rep_hist_chromA <- late_rep_hist_chromA[, !(names(late_rep_hist_chromA) %in% drops)]
    late_rep_hist_chromA$sum_Signal <- as.numeric(late_rep_hist_chromA$sum_Signal)
    late_rep_hist_chromA$sum_Signal <- late_rep_hist_chromA$sum_Signal / mean_late
    chr_A_late <- paste0(late_rep_hist_chromA$V1, ":", late_rep_hist_chromA$start, "-", late_rep_hist_chromA$V3)
    chr_A_late <- gsub(" ", "", chr_A_late)

    # Chromosome B- histone modification dataset
    chrNo <- paste("chr", q, sep = "")
    # get subset of chrB histone modification data
    hist_mod_chromB <- subset(signal_sum_df, signal_sum_df$chr == chrNo)
    # convert sum_signal in hist_mod_chromB into a numeric form for use in calculating D(i.j) matrix
    # later
    # hist_mod_chromB$sum_Signal <- as.numeric(hist_mod_chromB$sum_Signal)
    chr_B <- paste0(hist_mod_chromB$chr, ".", hist_mod_chromB$start, ".", hist_mod_chromB$end)
    chr_B <- gsub(" ", "", chr_B)

    # Replication timing dataset for Chromosome B
    # Early Replication
    early_rep_chromB <- subset(rep_early_df, rep_early_df$V1 == chrNo)
    early_rep_hist_chromB <- inner_join(early_rep_chromB, hist_mod_chromB, by = "start")
    drops <- c("chr", "end")
    early_rep_hist_chromB <- early_rep_hist_chromB[, !(names(early_rep_hist_chromB) %in% drops)]
    early_rep_hist_chromB$sum_Signal <- as.numeric(early_rep_hist_chromB$sum_Signal)
    early_rep_hist_chromB$sum_Signal <- early_rep_hist_chromB$sum_Signal / mean_early
    chr_B_early <- paste0(early_rep_hist_chromB$V1, ".", early_rep_hist_chromB$start, ".", early_rep_hist_chromB$V3)
    chr_B_early <- gsub(" ", "", chr_B_early)

    # Late Replication
    late_rep_chromB <- subset(rep_late_df, rep_late_df$V1 == chrNo)
    late_rep_hist_chromB <- inner_join(late_rep_chromB, hist_mod_chromB, by = "start")
    drops <- c("chr", "end")
    late_rep_hist_chromB <- late_rep_hist_chromB[, !(names(late_rep_hist_chromB) %in% drops)]
    late_rep_hist_chromB$sum_Signal <- as.numeric(late_rep_hist_chromB$sum_Signal)
    late_rep_hist_chromB$sum_Signal <- late_rep_hist_chromB$sum_Signal / mean_late
    chr_B_late <- paste0(late_rep_hist_chromB$V1, ".", late_rep_hist_chromB$start, ".", late_rep_hist_chromB$V3)
    chr_B_late <- gsub(" ", "", chr_B_late)

    # HiC dataset
    file <- paste0("HIC_gm06690_chr", p, "_", "chr", q, "_", "1000000_", "pearson", ".txt")
    whole_matrix <- readLines(file)
    ignore <- whole_matrix[-c(1:1)]

    # For Early Replicating fragments
    chrAB_HiC_E <- read.csv(textConnection(ignore), sep = "\t", stringsAsFactors = FALSE, row.names = 1)
    drops <- c("X.1")
    chrAB_HiC_E <- chrAB_HiC_E[, !(names(chrAB_HiC_E) %in% drops)]
    chromoA <- rownames(chrAB_HiC_E)
    chromoB <- names(chrAB_HiC_E)


    chrAB_HiC_E <- chrAB_HiC_E[grep(paste(chr_A_early, collapse = "|"), chromoA, value = TRUE),]
    chrAB_HiC_E <- t(chrAB_HiC_E)
    chrAB_HiC_E <- chrAB_HiC_E[grep(paste(chr_B_early, collapse = "|"), chromoB, value = TRUE),]
    chrAB_HiC_E <- t(chrAB_HiC_E)
    chrAB_HiC_E <- as.matrix(chrAB_HiC_E)
    chrAB_HiC_E <- matrix(chrAB_HiC_E, ncol = ncol(chrAB_HiC_E), dimnames = NULL)
    early_hic_vec <- c()
    for (x in 1:nrow(chrAB_HiC_E)) {
      if (p == q) {
        if (x == nrow(chrAB_HiC_E)) {
          next
        }
        for (y in (x + 1):ncol(chrAB_HiC_E)) {
          early_hic_vec <- c(early_hic_vec, chrAB_HiC_E[x, y])
        }
      } else {
        for (y in 1:ncol(chrAB_HiC_E)) {
          early_hic_vec <- c(early_hic_vec, chrAB_HiC_E[x, y])
        }
      }
    }

    early_dij_vec <- c()
    for (x in 1:nrow(early_rep_hist_chromA)) {
      if (p == q) {
        if (x == nrow(early_rep_hist_chromA)) {
          next
        }
        for (y in (x + 1):nrow(early_rep_hist_chromB)) {
          early_dij_vec <- c(early_dij_vec, abs(log10(early_rep_hist_chromA[x, "sum_Signal"] / early_rep_hist_chromB[y, "sum_Signal"])))
        }
      } else {
        for (y in 1:nrow(early_rep_hist_chromB)) {
          early_dij_vec <- c(early_dij_vec, abs(log10(early_rep_hist_chromA[x, "sum_Signal"] / early_rep_hist_chromB[y, "sum_Signal"])))
        }
      }
    }
    all_Dij_Early <- append(all_Dij_Early, early_dij_vec)
    all_HiC_Early <- append(all_HiC_Early, early_hic_vec)


    # For Late Replicating fragments
    chrAB_HiC_L <- read.csv(textConnection(ignore), sep = "\t", stringsAsFactors = FALSE, row.names = 1)
    drops <- c("X.1")
    chrAB_HiC_L <- chrAB_HiC_L[, !(names(chrAB_HiC_L) %in% drops)]
    chromoA <- rownames(chrAB_HiC_L)
    chromoB <- names(chrAB_HiC_L)


    chrAB_HiC_L <- chrAB_HiC_L[grep(paste(chr_A_late, collapse = "|"), chromoA, value = TRUE),]
    chrAB_HiC_L <- t(chrAB_HiC_L)
    # chrAB_HiC <- chrAB_HiC[,which(chromoB %in% chr_B)]
    chrAB_HiC_L <- chrAB_HiC_L[grep(paste(chr_B_late, collapse = "|"), chromoB, value = TRUE),]
    chrAB_HiC_L <- t(chrAB_HiC_L)
    chrAB_HiC_L <- as.matrix(chrAB_HiC_L)
    chrAB_HiC_L <- matrix(chrAB_HiC_L, ncol = ncol(chrAB_HiC_L), dimnames = NULL)
    late_hic_vec <- c()
    for (x in 1:nrow(chrAB_HiC_L)) {
      if (p == q) {
        if (x == nrow(chrAB_HiC_L)) {
          next
        }
        for (y in (x + 1):ncol(chrAB_HiC_L)) {
          late_hic_vec <- c(late_hic_vec, chrAB_HiC_L[x, y])
        }
      } else {
        for (y in 1:ncol(chrAB_HiC_L)) {
          late_hic_vec <- c(late_hic_vec, chrAB_HiC_L[x, y])
        }
      }
    }

    late_dij_vec <- c()
    for (x in 1:nrow(late_rep_hist_chromA)) {
      if (p == q) {
        if (x == nrow(late_rep_hist_chromA)) {
          next
        }
        for (y in (x + 1):nrow(late_rep_hist_chromB)) {
          late_dij_vec <- c(late_dij_vec, abs(log10(late_rep_hist_chromA[x, "sum_Signal"] / late_rep_hist_chromB[y, "sum_Signal"])))
        }
      } else {
        for (y in 1:nrow(late_rep_hist_chromB)) {
          late_dij_vec <- c(late_dij_vec, abs(log10(late_rep_hist_chromA[x, "sum_Signal"] / late_rep_hist_chromB[y, "sum_Signal"])))
        }
      }
    }

    all_Dij_Late <- append(all_Dij_Late, late_dij_vec)
    all_HiC_Late <- append(all_HiC_Late, late_hic_vec)

    # print(paste0("Chromosome A: ", p))
    # print(paste0("Chromosome B: ", q))
  }
}
cat("\n")
cat("Generating plot... \n")
setwd(argv$o)
minimum <- -0.32
maximum <- 0.57

# for early replicating DNA
HiC_D <- data.frame(HiC_distance = all_HiC_Early, D = all_Dij_Early)
HiC_D <- subset(HiC_D, HiC_distance > minimum)
HiC_D <- subset(HiC_D, HiC_distance < maximum)
intervals <- 18.0
intervalWidth <- (maximum - minimum) / intervals
HiC_D <- mutate(HiC_D, Grouped_HiC_distance = ceiling((HiC_distance - minimum) / intervalWidth))
intervalNames <- seq((minimum + intervalWidth / 2.0), (maximum - intervalWidth / 2.0), intervalWidth)
Rep_time <- rep(c("Early-Early"), nrow(HiC_D))
HiC_D <- mutate(HiC_D, Rep_time)
HiC_D_E <- HiC_D

# for late replicating DNA
HiC_D <- data.frame(HiC_distance = all_HiC_Late, D = all_Dij_Late)
HiC_D <- subset(HiC_D, HiC_distance > minimum)
HiC_D <- subset(HiC_D, HiC_distance < maximum)
intervals <- 18.0
intervalWidth <- (maximum - minimum) / intervals
HiC_D <- mutate(HiC_D, Grouped_HiC_distance = ceiling((HiC_distance - minimum) / intervalWidth))
intervalNames <- seq((minimum + intervalWidth / 2.0), (maximum - intervalWidth / 2.0), intervalWidth)
Rep_time <- rep(c("Late-Late"), nrow(HiC_D))
HiC_D <- mutate(HiC_D, Rep_time)
HiC_D_L <- HiC_D

# combining datasets
comb_data <- rbind(HiC_D_E, HiC_D_L)

median_early <- as.numeric(argv$a)
median_late <- as.numeric(argv$b)
# ggplot combined boxplot
myplot <- ggplot(comb_data, aes(x = HiC_distance, y = D, fill = Replication - Time)) +
    xlab("Spatial Proximity") +
    ylab("Q(i,j)") +
    ggtitle(paste0(mod)) +
    theme_bw() +
    theme(axis.text = element_text(size = 15), legend.position = c(0.70, 0.90), legend.text = element_text(size = 20), legend.title = element_blank(), axis.title = element_text(size = 30), legend.key = element_blank(), legend.background = element_blank()) +
    geom_boxplot(aes(fill = factor(Rep_time, levels = c("Early-Early", "Late-Late")), cut_width(HiC_distance, 0.2)), outlier.shape = NA, outlier.alpha = 0.1, coef = 0) + coord_cartesian(ylim = c(0.1, 1.4)) + geom_hline(aes(yintercept = median_early), color = "#F8766D", linetype = "dashed", size=1.5) + geom_hline(aes(yintercept = median_late), color = "#00BFC4", linetype = "dashed", size=1.5)
pdf(paste0("ggplot_combined_", mod, "_1000_rep-time_cutoff.pdf"))
suppressWarnings(print(myplot))
invisible(dev.off())
cat(paste0("\nAll Plots Saved in:", argv$o))
cat("\n")