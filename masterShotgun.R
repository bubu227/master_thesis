
# libraries
library(tidyverse)
library(ggplot2)
library(reshape2)
library(readxl)
library(data.table)
library(ggvenn)
library(stringr)
library(rtracklayer)

# plotting functions
wd = getwd()
source("/home/niklas/mount/netscratch/dep_psl/grp_rgo/kieln/plotting_functions/plotting_functions.R")
source("/home/niklas/mount/netscratch/dep_psl/grp_rgo/kieln/reanalysis_assay_7_shotgun/analysis/analysis_functions.R")

# map the snps to mutations
gff = loadGff()
gffCols = c("ID", "product", "KO", "description")
gffSubset = gff[, gffCols]

#load KO terms
koTerms = load_ko()

# load mapping data
assay_dir = "/home/niklas/mount/netscratch/dep_psl/grp_rgo/kieln/reanalysis_assay_7_shotgun"
mapping_path = paste0(assay_dir, "/mapping/design_file.xlsx")
design = read_excel(mapping_path, sheet = "design")

design$condition = NA
idx = design$syncom == "at" & design$plant == "col0" | design$syncom == "lj" & design$plant == "gifu"
design$condition[idx] = "native"

idx = design$syncom == "at" & design$plant == "gifu" | design$syncom == "lj" & design$plant == "col0"
design$condition[idx] = "non-native"

#extract sample names
samples = design$library_number

# figure dir
fig_dir = paste0(assay_dir, "/figures_master")
dir.create(fig_dir, showWarnings = FALSE)

# load the parsed vcf files
variant_dir = paste0(assay_dir, "/variants_gatk_bowtie")
paths = list.files(variant_dir, pattern= glob2rx("*_mutect_mapped.tsv"), full.names=T, recursive=T)
df_list = sapply(paths, read.delim, simplify=FALSE, sep = "\t", header = TRUE, comment.char = '#')
snp_df = rbindlist(df_list, idcol="folder", fill = TRUE) %>% mutate(folder = basename(folder))
snp_df$sample = gsub("*_mutect_mapped.tsv", "", snp_df$folder)

# reformat EFF column
snp_df$EFF = gsub("\\(", "|", snp_df$EFF)
snp_df$EFF = gsub("\\)", "", snp_df$EFF)
effects = str_split_fixed(snp_df$EFF, "\\|", 11)
effects = data.frame(effects)
effect_cols = c("type", "impact", "effect", "nt", "aa", "gene", "gene_type", "transcript_type", "transcript", "exon", "warning")
colnames(effects) = effect_cols

# add effects to dataframe
snp_df = snp_df %>% select(-EFF)
snp_df = bind_cols(snp_df, effects) 

# add meta data
snp_df$strain = gsub("_.*","", snp_df$CHROM) 
snp_df = snp_df %>% 
            left_join(taxo_table, by = c("strain" = "Sample_ID")) %>%
            left_join(color_family, by = c("family" = "ID")) %>% 
            left_join(design, by = c("sample" = "library_number"))

# filter data according to mapping quality
colnames(snp_df)
dim(snp_df)

# filter low quality calls
idx = sapply(str_split(snp_df$MBQ, ","), FUN = function(x) min(as.numeric(x))) > 30
snp_df = snp_df[idx, ]

idx = sapply(str_split(snp_df$MMQ, ","), FUN = function(x) min(as.numeric(x))) > 30
snp_df = snp_df[idx, ]


# type of different mutations
type_counts = snp_df %>% count(type) %>% arrange(-n)
type_counts$type = factor(type_counts$type, levels = unique(type_counts$type))

p = ggplot(type_counts, aes(x = type, y = n, fill = type)) + 
    geom_bar(stat="identity", alpha = 0.4, color = "black") + 
    main_theme + 
    theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
    theme(legend.position = "none",
            plot.margin = margin(10, 10, 10, 100))

fig_name = "MasterMutCounts"
save_fig(fig_name, p)


# remove synonymous coding
rm_types = c("INTERGENIC", "TRANSCRIPT", "SYNONYMOUS_STOP", "SYNONYMOUS_CODING")
non_syn = snp_df %>% 
            filter(!type %in% rm_types)

non_syn_counts = non_syn %>% count(type) %>% arrange(-n)
non_syn_counts = change_type(non_syn_counts)
non_syn_counts$type = factor(non_syn_counts$type, levels = unique(non_syn_counts$type))
#non_syn_counts$n = log2(non_syn_counts$n + 1)

p = ggplot(non_syn_counts, aes(x = type, y = n, fill = type)) + 
    geom_bar(stat="identity", alpha = 0.4) + 
    main_theme + 
    theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
    theme(legend.position = "none",
            plot.margin = margin(10, 10, 10, 20)) + 
    ggtitle(paste0("non-synonymous SNPs: ", sum(non_syn_counts$n))) + 
    ylab("Number of mutations") + 
    geom_text(aes(label = n)) + 
    xlab("")

fig_name = "MasterMutCountsNonSym"
save_fig(fig_name, p, width = 4, height = 6)



# look at number of mutations per strain
counts = non_syn %>%
    count(strain) %>% 
    left_join(taxo_table, by = c("strain" = "Sample_ID")) %>%
    arrange(-n) %>%
    left_join(color_family, by = c("family" = "ID"))

counts$strain = factor(counts$strain, levels = unique(counts$strain))
colors = counts$Color
names(colors) = counts$family

p = ggplot(counts, aes(x = strain, y = n, fill = family)) + 
    geom_bar(stat = "identity") + 
    main_theme + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
    theme(legend.position = "top",
          legend.direction = "vertical") + 
    scale_fill_manual("Taxonomic family", values=colors) + 
    ylab("Mutations per strain")

fig_name = "mutationsPerStrain"
save_fig(fig_name, p, width = 5, height = 8)


# splitting by host
# enrichment per host venn diagram
non_syn = non_syn %>% 
            left_join(gffSubset, by = c("gene" = "ID")) %>% 
            left_join(koTerms, by = c("KO" = "ko")) %>% 
            filter(product != "hypothetical protein")

vennCol0 = non_syn %>% filter(plant == "col0")
vennGifu = non_syn %>% filter(plant == "gifu")

# VENN using KO terms
library(VennDiagram)

vennList = list(
  Col0 = vennCol0$KO,
  Gifu = vennGifu$KO
)

venn.diagram(
  x = vennList,
  category.names = c("Col-0" , "Gifu"),
  filename = paste0(fig_dir, "/vennPlantKO.png"),
  output = TRUE,
  disable.logging = TRUE,
  imagetype = "png",

  # Circles
  lwd = 1,
  fill = colors_host$color[1:2],

  # Numbers
  cex = 1.5,
  fontface = "plain",
  fontfamily = "sans",
  print.mode = c("raw", "percent"),

  # Set names
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.default.pos = "text",
  cat.fontfamily = "sans",
  cat.just = list(
                  c(0.5, -12),
                  c(0.5, -12))
)


# extract the overlapping terms
intersection = Reduce(intersect, vennList)
unique(koTerms$C[koTerms$ko %in% intersection[1]])


# KO enrichment for hosts
koCounts = non_syn %>% 
            filter(!KO %in% intersection) %>%
            count(C, plant) %>% 
            arrange(-n) %>% 
            drop_na(C) %>%
            group_by(plant) %>%
            slice_max(order_by = n, n = 5)

koCounts$C[koCounts$C == "Biofilm formation - Pseudomonas aeruginosa"] = "Biofilm formation"
koCounts$plant[koCounts$plant == "col0"] = "Col-0"
koCounts$plant[koCounts$plant == "gifu"] = "Gifu"

p =  ggplot(koCounts, aes(x = reorder_within(C, -n, plant), y = n, fill = C)) + 
  geom_bar(stat = "identity", alpha = 0.3, color = "black") + 
  facet_wrap(~plant, scales = "free_x") + 
  main_theme  + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  scale_x_reordered() + 
  theme(legend.position = "none") + 
  xlab("") + 
  ylab("Counts of KO-terms")

fig_name = "koEnrichmentHost"
save_fig(fig_name, p, width = 5, height = 5)


# exclude KO terms that are present in both hosts
koCounts = non_syn %>% 
            filter(!KO %in% intersection) %>%
            count(C, plant) %>% 
            arrange(-n) %>% 
            drop_na(C) %>%
            group_by(plant) %>%
            slice_max(order_by = n, n = 20)

koCountsCol0 = koCounts %>% filter(plant == "col0")
koCountsGifu = koCounts %>% filter(plant == "gifu")

intersectionPlant = Reduce(intersect, list(koCountsCol0$C, koCountsGifu$C))
koCountsUnique = koCounts %>% 
                    filter(!C %in% intersectionPlant) %>% 
                    arrange(-n) %>% 
                    group_by(plant) %>%
                    slice_max(order_by = n, n = 3)

koCountsUnique$C[koCountsUnique$C == "Biofilm formation - Pseudomonas aeruginosa"] = "Biofilm formation"
koCountsUnique$plant[koCountsUnique$plant == "col0"] = "Col-0"
koCountsUnique$plant[koCountsUnique$plant == "gifu"] = "Gifu"

p =  ggplot(koCountsUnique, aes(x = reorder_within(C, -n, plant), y = n, fill = C)) + 
  geom_bar(stat = "identity", alpha = 0.3, color = "black") + 
  facet_wrap(~plant, scales = "free_x") + 
  main_theme  + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  scale_x_reordered() + 
  theme(legend.position = "none") + 
  xlab("") + 
  ylab("Counts of KO-terms")

fig_name = "koEnrichmentHostUnique"
save_fig(fig_name, p, width = 6, height = 6)


# Gene enrichments

# remove genes present in both hosts
vennList = list(
  Col0 = vennCol0$gene.y,
  Gifu = vennGifu$gene.y
)

intersection = Reduce(intersect, vennList)

non_syn$product = sapply(non_syn$product, FUN = function(x) paste0(x, collapse = "/"))
non_syn$product = gsub(" ", "", non_syn$product)

GeneCounts = non_syn %>% 
            filter(!gene.y %in% intersection) %>%
            count(gene.y, plant) %>% 
            arrange(-n) %>% 
            drop_na(gene.y) %>%
            group_by(plant) %>%
            slice_max(order_by = n, n = 10)

GeneCounts$plant[GeneCounts$plant == "col0"] = "Col-0"
GeneCounts$plant[GeneCounts$plant == "gifu"] = "Gifu"
 
# reformat product names


# add the KO pathways
koSubset = non_syn %>% 
            filter(!gene.y %in% intersection) %>%
            select(gene.y, B)

toRemove = c("Aging", "Endocrine system")

getLast = function(x) return(x[1])

geneKOC = GeneCounts %>% left_join(koSubset) %>% distinct()
geneKOC = geneKOC %>% 
            filter(!str_detect(B, "Protein families:")) %>%
            filter(!str_detect(B, "Brite Hierarchies")) %>%
            filter(!B %in% toRemove) %>%
            group_by(gene.y, plant, n) %>% 
            summarise("KO-pathway" = getLast(B))

unique(geneKOC$`KO-pathway`)
geneKOC$gene.y = gsub(" \\[.*", "", geneKOC$gene.y)

p =  ggplot(geneKOC, aes(x = reorder_within(gene.y, n, plant), y = n, fill = `KO-pathway`)) + 
  geom_bar(stat = "identity", alpha = 0.6, color = "black") + 
  facet_wrap(~plant, scales = "free_y", ncol = 1) + 
  main_theme  +
  theme(strip.text.x = element_text(size = 24),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 22),
        legend.direction="vertical",
        legend.title = element_text(size=22),
        legend.text = element_text(size=20)) +
  scale_x_reordered() + 
  theme(legend.position = "bottom") + 
  xlab("") + 
  ylab("Number of variants per gene") + 
  coord_flip() + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_fill_manual("KO-terms", values = c25)

fig_name = "GeneEnrichmentHost"
save_fig(fig_name, p, width = 14, height = 18)



