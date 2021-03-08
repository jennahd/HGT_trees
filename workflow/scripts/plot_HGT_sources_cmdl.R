#!/usr/bin/env Rscript
#Load libraries
suppressPackageStartupMessages(library("ggforce"))
suppressPackageStartupMessages(library("gridExtra"))
suppressPackageStartupMessages(library("ggpubr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("argparse"))

#Command-line Implemetation
parser <- ArgumentParser(description="Plot taxonomy (domain, superphylum, & phylum) of sister and nested clades to taxa of interest.")

parser$add_argument("-in", "--input_file", help = "Path to the input tsv file containing cluster ID, cutoffs, sister taxonomy, and nested taxonomy.")
parser$add_argument("-out", "--output_file", help = "Path to the pdf file that will be output with alluvial plots of taxonomy at different cutoffs.")
args <- parser$parse_args()
print(args)

#Superphylum colours
super = c("SS de-novo"="#434343",
          "SS cellular_organisms"="#ff7f00",
          "SS Bacteria"="#33a02c",
          "SS Archaea"="#1f78b4",
          "SS DPANN_group"="#1f78b4",
          "SS Euryarchaeota"="#1f78b4",
          "SS Asgard_group"="#1f78b4",
          "SS TACK_group"="#1f78b4",
          "SS Eukaryota"="#e31a1c",
          "SS Amoebozoa" = "#e31a1c",
          "SS CRuMs" = "#e31a1c",
          "SS Discoba" = "#e31a1c",
          "SS Glaucocystophyceae" = "#e31a1c",
          "SS Haptista" = "#e31a1c",
          "SS Rhodophyta"="#e31a1c",
          "SS Sar"="#e31a1c",
          "SS Opisthokonta"="#e31a1c",
          "SS Viridiplantae" = "#e31a1c",
          "SS PVC_group"="#8dd3c7",
          "SS FCB_group"="#ffffb3",
          "SS Terrabacteria_group"="#bebada",
          "SS Nitrospirae"="#fb8072",
          "SS Acidobacteria" = "#fdb462",
          "SS Proteobacteria"="#80b1d3",
          "SS Aquificae" = "#b3de69",
          "SS Thermotogae" = "#fccde5",
          "SS Spirochaetes" = "#bc80bd",
          "SS Elusimicrobia" = "#ccebc5",
          "SS Synergistetes" = "#ffed6f",
          "SS Caldiserica_Cryosericota_group" = "#d9d9d9",
          "SS Calditrichaeota" = "#d9d9d9",
          "SS Chrysiogenetes" = "#d9d9d9",
          "SS Coprothermobacterota" = "#d9d9d9",
          "SS Deferribacteres" = "#d9d9d9",
          "SS Dictyoglomi" = "#d9d9d9",
          "SS Krumholzibacteriota" = "#d9d9d9",
          "SS Thermodesulfobacteria" = "#d9d9d9",
          "SS Fusobacteria" = "#d9d9d9",
          "SS unclassified_Bacteria"="#d9d9d9",
          "SS Bacteria_incertae_sedis"="#d9d9d9")

HGT<-read.table(args$input_file, header=F, sep="\t")
names(HGT) <- c("ID", "cluster", "cutoff",
                "sister_domain","sister_superphylum","sister_phylum",
                "nested_domain","nested_superphylum","nested_phylum")
HGT <- as_tibble(HGT)

cutoffs <- c("0.9", "0.75", "0.5")
pltList <- list()

for (i in cutoffs) {
  #get cutoff data
  strict_HGT <- filter(HGT, cutoff == i)
  strict_HGT <- strict_HGT[, c(4,5,7,8)]
  #Add identifier for taxonomic level
  strict_HGT$sister_domain = paste0('SD ', strict_HGT$sister_domain)
  strict_HGT$sister_superphylum = paste0('SS ', strict_HGT$sister_superphylum)
  strict_HGT$nested_domain = paste0('ND ', strict_HGT$nested_domain)
  strict_HGT$nested_superphylum = paste0('NS ', strict_HGT$nested_superphylum)
  #get frequency and re-order factors
  freq <- strict_HGT %>% group_by(sister_domain, sister_superphylum,
                                  nested_superphylum, nested_domain) %>% tally()
  freq2 <- gather_set_data(freq, 1:4)
  freq2$y <- factor(freq2$y, levels=unique(freq2$y))
  freq2$x <- factor(freq2$x, levels=c("sister_domain","sister_superphylum",
                                      "nested_superphylum","nested_domain"))
  #make plot
  sister_strict <- ggplot(freq2, aes(x, id = id, split = y, value = n)) +
    geom_parallel_sets(aes(fill = sister_superphylum), alpha = 0.5, axis.width = 0.6) +
    geom_parallel_sets_axes(axis.width = 0.6) +
    geom_parallel_sets_labels(colour = 'grey', angle = 0, size = 3, vjust=-0.3) +
    scale_fill_manual(values = super) +
    theme_transparent() +
    theme(
      axis.text = element_text(size = 0),
      axis.text.x = element_text(size = 10, colour = "#343434", margin=margin(5,5,5,5,"pt")),
      axis.text.y = element_text(size = 8, colour = "#343434"),
      axis.line.y = element_line(colour = "#343434"),
      axis.ticks.y = element_line(colour = "#343434"),
      axis.ticks.x = element_line(colour = "#343434"),
      axis.ticks.length.x = unit(10, "pt"),
      legend.position = "none",
      plot.title=element_text(size=16,face="bold", hjust=1)) +
    ggtitle(i)
  pltList[[i]]<- print(sister_strict)
}

pdf(args$output_file,width=20,height=25)

grid.arrange(grobs = pltList[1:3], nrow=3)

dev.off()
