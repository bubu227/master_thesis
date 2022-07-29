# set wd bc r on macos is stupid
#setwd("/Users/nkiel/MPIPZ/biodata/dep_psl/grp_rgo/kieln/inv_experimental_evolution/study_isolation_protocol/assay_6_miseq/analysis/")

# load packages
library(ggplot2)
library(tidyverse)
library(reshape)
library(psych)

# load data and functions
wd = getwd()
source("load_data.R")


#########################################################################################################################
#does inoculation of two week old plants also result in a host preference phenotype???
idx = str_detect(design$syncom, "mixed")
design_sub = design[idx, ]
otu_list = by_condition(idx)
otu_tab = otu_table_norm[, idx]
otu_tab_abs = otu_table[, idx]


at_strains_ra <- colSums(otu_tab[rownames(otu_tab) %in% at_strains, ])
lj_strains_ra <- colSums(otu_tab[rownames(otu_tab) %in% lj_strains, ])

host_pref_ra = rbind(otu_tab, at_strains_ra, lj_strains_ra)


df = melt(host_pref_ra)
colnames(df) <- c("strain", "SampleID", "RA")
df = df %>% left_join(design, by = c("SampleID" = "#SampleID"))

idx <- match(df$SampleID, design$`#SampleID`)
df$treatment <- design$plant[idx]

df$family <- taxonomy$family[match(df$strain, taxonomy$strain)]
df$family <- as.character(df$family)
df$family[is.na(df$family)] <- "All strains"

df$strain <- as.character(df$strain)
df$strain[df$strain=="at_strains_ra"] <- "At-SPHERE strains"
df$strain[df$strain=="lj_strains_ra"] <- "Lj-SPHERE strains"

df$family_strain <- paste(df$family, " (", df$strain, ")", sep="")

df$RA_log <- log10(df$RA * 1000 + 1)


fmt_dcimals <- function(decimals=0) {
  function(x) format(x, nsmall=decimals, scientific=F)
}

family <- "All strains"

df_family <- df[which(df$family==family), ]

if (family!="All strains") df_family$strain <- factor(df_family$strain, levels=rev(sort(unique(df_family$strain))))

#colors <- colors[colors$group %in% df_family$treatment, ]
df_family$treatment <- factor(df_family$treatment, levels=colors$group)
df_family <- unique(df_family)

df_family_at <- df_family[df_family$strain %in% c("At-SPHERE strains", at_strains), ]
df_family_lj <- df_family[df_family$strain %in% c("Lj-SPHERE strains", lj_strains), ] 

if (dim(df_family_at)[1] > 0 & dim(df_family_lj)[1] > 0) {
  
  # Kruskal-Wallis test of group differences
  
  pval <- kruskal.test(RA ~ treatment, data=df_family_at)$p.value
  
  lim_padding <- 0.15
  
  idx <- df_family_at$treatment %in% c("col0", "gifu")
  lim <- mean(df_family_at$RA[idx])
  h_lim <- min(1, lim+lim_padding)
  l_lim <- max(0, lim-lim_padding)
  
  my_y_title <- expression(paste(italic("At"), "-SPHERE strains"))
  
  p1 <- ggplot(df_family_at, aes(x=treatment, y=RA, color=treatment)) +
    geom_boxplot(alpha=1, outlier.size=0, size=1, width=.2*length(unique(df_family$treatment)), fill = NA) +
    geom_jitter(position=position_jitterdodge(jitter.width=.1*length(unique(df_family$treatment)),
                                              dodge.width=.2*length(unique(df_family$treatment))), size=1, alpha=0.5) +
    scale_y_continuous(position="right", labels=percent, limits=c(l_lim, h_lim), breaks=seq(0, 1, by=.1)) +
    scale_colour_manual(values=as.character(colors$color)) +
    labs(x="", y="Relative Abundance") +
    theme(axis.title=element_text(size=9)) +
    theme(plot.margin=unit(c(5, 0, 2, 0), "mm")) +
    theme(plot.title = element_text(hjust = 0.5)) + 
    main_theme +
    theme(legend.position="none") + 
    ggtitle(my_y_title) + 
    theme(panel.background = element_rect(fill = "mistyrose",
                                          colour = "mistyrose"))+ 
    scale_x_discrete(labels= c("Gifu", "Col-0"))
  
  
  idx <- df_family_lj$treatment %in% c("col0", "gifu")
  lim <- mean(df_family_lj$RA[idx])
  h_lim <- min(1, lim+lim_padding)
  l_lim <- max(0, lim-lim_padding)
  my_y_title <- expression(paste(italic("Lj"), "-SPHERE strains"))
  
  p2 <- ggplot(df_family_lj, aes(x=treatment, y=RA, color=treatment)) +
    geom_boxplot(alpha=1, outlier.size=0, size=1, width=.2*length(unique(df_family$treatment)),fill = NA) +
    geom_jitter(position=position_jitterdodge(jitter.width=.1*length(unique(df_family$treatment)),
                                              dodge.width=.2*length(unique(df_family$treatment))), size=1, alpha=0.5) +
    scale_y_continuous(position="left", labels=percent, limits=c(l_lim, h_lim), breaks=seq(0, 1, by=.1)) +
    scale_colour_manual(values=as.character(colors$color)) +
    labs(x="", y="Relative Abundance") +
    theme(axis.title=element_text(size=9)) +
    theme(plot.margin=unit(c(5, 0, 2, 5), "mm")) +
    theme(plot.title = element_text(hjust = 0.5)) + 
    main_theme +
    theme(legend.position="none") + 
    ggtitle(my_y_title) + 
    theme(panel.background = element_rect(fill = "azure",
                                          colour = "azure")) + 
    scale_x_discrete(labels= c("Col-0", "Gifu"))
  
  
  gA <- ggplotGrob(p2)
  gB <- ggplotGrob(p1)
  maxHeight = grid::unit.pmax(gA$heights, gB$heights)
  gA$heights <- as.list(maxHeight)
  gB$heights <- as.list(maxHeight)
  ttl = "Host preference at t=0\n"
  pg1 <- grid.arrange(gA, gB, ncol=2, top=textGrob(paste(ttl, "pval=", format(pval, digits=3), sep=""), gp=gpar(fontsize=20)))

  fname = "host_preference_gen_0"
  f = paste0(fig_dir, fname, ".pdf")
  ggsave(f, pg1, width = 10, height = 7)
  f = paste0(fig_dir, fname, ".png")
  ggsave(f, pg1, width = 10, height = 7)
}


######################################################################################################
# check for host preference after ONE generation
idx = str_detect(design$sample_name, "comp")
otu = otu_table_norm[, idx]


at_strains_ra <- colSums(otu[rownames(otu) %in% at_strains, ])
lj_strains_ra <- colSums(otu[rownames(otu) %in% lj_strains, ])

host_pref_ra = rbind(otu, at_strains_ra, lj_strains_ra)

df = melt(host_pref_ra)
colnames(df) <- c("strain", "SampleID", "RA")
df = df %>% left_join(design, by = c("SampleID" = "#SampleID"))

idx <- match(df$SampleID, design$`#SampleID`)
df$treatment <- design$plant[idx]

df$family <- taxonomy$family[match(df$strain, taxonomy$strain)]
df$family <- as.character(df$family)
df$family[is.na(df$family)] <- "All strains"

df$strain <- as.character(df$strain)
df$strain[df$strain=="at_strains_ra"] <- "At-SPHERE strains"
df$strain[df$strain=="lj_strains_ra"] <- "Lj-SPHERE strains"

df$family_strain <- paste(df$family, " (", df$strain, ")", sep="")

df$RA_log <- log10(df$RA * 1000 + 1)


fmt_dcimals <- function(decimals=0) {
  function(x) format(x, nsmall=decimals, scientific=F)
}

family <- "All strains"

df_family <- df[which(df$family==family), ]

if (family!="All strains") df_family$strain <- factor(df_family$strain, levels=rev(sort(unique(df_family$strain))))

#colors <- colors[colors$group %in% df_family$treatment, ]
df_family$treatment <- factor(df_family$treatment, levels=colors$group)
df_family <- unique(df_family)

df_family_at <- df_family[df_family$strain %in% c("At-SPHERE strains", at_strains), ]
df_family_lj <- df_family[df_family$strain %in% c("Lj-SPHERE strains", lj_strains), ] 

if (dim(df_family_at)[1] > 0 & dim(df_family_lj)[1] > 0) {
  
  # Kruskal-Wallis test of group differences
  
  pval <- kruskal.test(RA ~ treatment, data=df_family_at)$p.value
  
  lim_padding <- 0.3
  
  idx <- df_family_at$treatment %in% c("col0", "gifu")
  lim <- mean(df_family_at$RA[idx])
  h_lim <- min(1, lim+lim_padding)
  l_lim <- max(0, lim-lim_padding)
  
  my_y_title <- expression(paste(italic("At"), "-SPHERE strains"))
  
  p1 <- ggplot(df_family_at, aes(x=treatment, y=RA, color=treatment)) +
    geom_boxplot(alpha=1, outlier.size=0, size=1, width=.2*length(unique(df_family$treatment)), fill = NA) +
    geom_jitter(position=position_jitterdodge(jitter.width=.1*length(unique(df_family$treatment)),
                                              dodge.width=.2*length(unique(df_family$treatment))), size=1, alpha=0.5) +
    scale_y_continuous(position="right", labels=percent, limits=c(l_lim, h_lim), breaks=seq(0, 1, by=.1)) +
    scale_colour_manual(values=as.character(colors$color)) +
    labs(x="", y="Relative Abundance") +
    theme(axis.title=element_text(size=9)) +
    theme(plot.margin=unit(c(5, 0, 2, 0), "mm")) +
    theme(plot.title = element_text(hjust = 0.5)) + 
    main_theme +
    theme(legend.position="none") + 
    ggtitle(my_y_title) + 
    theme(panel.background = element_rect(fill = "mistyrose",
                                          colour = "mistyrose"))+ 
    scale_x_discrete(labels= c("Gifu", "Col-0"))
  
  idx <- df_family_lj$treatment %in% c("col0", "gifu")
  lim <- mean(df_family_lj$RA[idx])
  h_lim <- min(1, lim+lim_padding)
  l_lim <- max(0, lim-lim_padding)
  my_y_title <- expression(paste(italic("Lj"), "-SPHERE strains"))
  
  p2 <- ggplot(df_family_lj, aes(x=treatment, y=RA, color=treatment)) +
    geom_boxplot(alpha=1, outlier.size=0, size=1, width=.2*length(unique(df_family$treatment)),fill = NA) +
    geom_jitter(position=position_jitterdodge(jitter.width=.1*length(unique(df_family$treatment)),
                                              dodge.width=.2*length(unique(df_family$treatment))), size=1, alpha=0.5) +
    scale_y_continuous(position="left", labels=percent, limits=c(l_lim, h_lim), breaks=seq(0, 1, by=.1)) +
    scale_colour_manual(values=as.character(colors$color)) +
    labs(x="", y="Relative Abundance") +
    theme(axis.title=element_text(size=9)) +
    theme(plot.margin=unit(c(5, 0, 2, 5), "mm")) +
    theme(plot.title = element_text(hjust = 0.5)) + 
    main_theme +
    theme(legend.position="none") + 
    ggtitle(my_y_title) + 
    theme(panel.background = element_rect(fill = "azure",
                                          colour = "azure")) + 
    scale_x_discrete(labels= c("Col-0", "Gifu"))
  
  gA <- ggplotGrob(p2)
  gB <- ggplotGrob(p1)
  maxHeight = grid::unit.pmax(gA$heights, gB$heights)
  gA$heights <- as.list(maxHeight)
  gB$heights <- as.list(maxHeight)
  ttle = "Host preference at t=1\n"
  pg1 <- grid.arrange(gA, gB, ncol=2, top=textGrob(paste(ttle, "pval=", format(pval, digits=3), sep=""), gp=gpar(fontsize=20)))
  
  fname = "host_preference_gen_1"
  f = paste0(fig_dir, fname, ".pdf")
  ggsave(f, pg1, width = 10, height = 7)
  f = paste0(fig_dir, fname, ".png")
  ggsave(f, pg1, width = 10, height = 7)
}


# ggmuller plots
library(tidyr)
hp_subset = design[match(df$SampleID, design$`#SampleID`), ]
lines = unique(hp_subset$syncom)
numbers = extract_numeric(lines)
lines = as.numeric(unlist(strsplit(as.character(numbers), "")))

# muller plots if lines of interest
idx = (design$Plate == 4 & str_detect(design$sample_name, "line")) | (design$Plate == 3 & str_detect(design$sample_name, "root"))
design[idx & str_detect(design$sample_name, "7"), ]

otu_list_split = split_syncoms(idx, otu_table)
otu_list_split = add_geneneration(otu_list_split)
otu_list_split$generation = as.numeric(otu_list_split$generation)

# get inputs
at_input = design[design$sample_name == "input-at__" & design$Plate == 3, ]
lj_input = design[design$sample_name == "input-lj__" & design$Plate == 3, ]
at_nodule_input = design[design$sample_name == "input-at-nodule" & design$Plate == 3, ]


# plot selected generations
p1 = plot_profile_over_time(lines[1])
p2 = plot_profile_over_time(lines[2])
p3 = plot_profile_over_time(lines[3])
p4 = plot_profile_over_time(lines[4])


# save the plots
g = ggarrange(p2, p4, p1, p3, ncol=2, nrow=2, common.legend = TRUE, legend="right")
f_name = "master_community_ver_time"
save_fig(g, f_name, width=8, height=7)




### HPI
# load packages
library(ggplot2)
library(tidyverse)
library(reshape)
library(psych)
library(ggrepel)
library(rstatix)



#########################################################################################################################
# host preference index of strains after two generations
idx = str_detect(design$syncom, "mixed")
gen0_hpi = calc_hpi(idx, otu_table_norm)
gen0_hpi$gen = 0


idx = str_detect(design$sample_name, "comp")
gen1_hpi = calc_hpi(idx, otu_table_norm)
gen1_hpi$gen = 1

hpi_merged = rbind(gen0_hpi, gen1_hpi)

lim = 6
hpi_merged$hpi_ratio[hpi_merged$hpi_ratio <= 1] = 0
hpi_merged$hpi_ratio[hpi_merged$hpi_ratio >= lim] = lim

hpi_merged$gen[hpi_merged$gen == 0] = "Generation 0"
hpi_merged$gen[hpi_merged$gen == 1] = "Generation 1"

p1 = ggplot(hpi_merged, aes(x=sphere, y=family)) +
             geom_point(aes(size=log2(mean_ra_host), fill=hpi_ratio, color=sig), shape=21) +
             labs(x="", y="") +
             scale_fill_gradientn(colors=c("darkgrey", "yellow", "red"), values=rescale(c(0, 1, lim)), limits=c(0, lim)) +
             scale_color_manual(values=c("transparent", "black"), guide = "none") +
             theme(axis.text.x=element_text(face = "italic")) +
             theme(axis.title=element_text(size=9)) +
             theme(plot.margin=unit(c(5, 0, 2, 5), "mm")) +
             main_theme +
             theme(legend.position="right") + 
             facet_wrap(~gen) + 
             scale_x_discrete(labels = c("At-SPHERE", "Lj-SPHERE")) + 
             guides(size=guide_legend(title="RA")) + 
             labs(fill = "HPI")

g = ggplot_build(p1)
size_keys = g$plot$scales$scales[[4]]$get_labels()
size_labels = paste0(round(2**as.numeric(size_keys), 3)*100, "%")
p1 = p1 + scale_size(labels = size_labels)

f_name = "master_hpi_dot"
save_fig(p1, f_name, width=10, height=7)



# plot with connected lines
idx = str_detect(design$syncom, "mixed")
gen0_hpi = calc_hpi(idx, otu_table_norm)
gen0_hpi$gen = 0
gen0_hpi$syncom = gen0_hpi$sphere

idx = str_detect(design$sample_name, "comp")
gen1_hpi = calc_hpi(idx, otu_table_norm)
gen1_hpi$gen = 1
gen1_hpi$syncom = gen1_hpi$sphere

hpi_merged = rbind(gen0_hpi, gen1_hpi)
hpi_merged$sphere = paste0(hpi_merged$sphere, hpi_merged$gen)

#hpi_merged$hpi_ratio[hpi_merged$hpi_ratio <= 1] = 0
hpi_merged$hpi_ratio[hpi_merged$hpi_ratio >= lim] = 6

#set order
order = c("at0", "lj0", "at1", "lj1")
hpi_merged$sphere = factor(hpi_merged$sphere, levels = order)

# get df for plottig
hpi_merged$hpi_log2 = log2(hpi_merged$hpi_ratio)
hpi_label = hpi_merged %>% filter(gen == 1)
x_offset = max(hpi_label$gen) #+ 0.3

# renaming for plotting
hpi_merged$syncom[hpi_merged$syncom == "at"] = "At-SPHERE"
hpi_merged$syncom[hpi_merged$syncom == "lj"] = "Lj-SPHERE"


p1 = ggplot(hpi_merged, aes(x=gen, y=hpi_log2, color = family)) + 
    geom_line() + 
    geom_point(aes(x=gen, y=hpi_log2, color = family, size = RA)) + 
    facet_wrap(~ syncom, ncol = 1, scales = "free") + 
    main_theme +
    scale_colour_manual(values = color_df$Color) +
    scale_x_continuous(limits = c(0, x_offset), breaks = c(0:max(hpi_label$gen))) +
    theme(legend.position="right") + 
    xlab("Generation") + 
    ylab("Host preference index (log2-transformed)") + 
    theme(axis.line.x=element_line(color = "black")) +
    theme(strip.text = element_text(face = "italic")) + 
    guides(color=guide_legend(title="Taxonomic family")) + 
    geom_hline(yintercept=0, linetype="dashed", color = "red")

f_name = "master_hpi_line"
save_fig(p1, f_name, width=7, height=9)



# calculate relative change of hpi
taxonomy = taxonomy[, c("strain", "sphere")] %>% 
  left_join(taxo_table, by = c("strain" = "ID")) %>% 
  select(family, order, class, phylum) %>%
  distinct()

hpi_fam$hpi_ratio[hpi_fam$hpi_ratio < 0.01] = 0.01

hpi_0 = hpi_fam %>% filter(gen == 0)
hpi_1 = hpi_fam %>% filter(gen == 1)

cols = c("family", "syncom", "hpi_ratio")

# merge stuff together
hpi_01 = hpi_0 %>% 
  left_join(hpi_1[, cols], by = c("family", "syncom")) %>% 
  mutate(ratio_dt = hpi_ratio.y / hpi_ratio.x) %>% 
  left_join(taxonomy)

hpi_01$ratio_dt[hpi_01$ratio_dt > lim] = lim
hpi_01$ratio_dt = log2(hpi_01$ratio_dt)

# reordering and factoring
hpi_01 = hpi_01 %>% arrange(desc(class))
hpi_01$family = factor(hpi_01$family, unique(hpi_01$family))

# get min and max vals
plt_max = ceiling(max(hpi_01$ratio_dt))
plt_min = -ceiling(abs(min(hpi_01$ratio_dt)))

# plot change over generation
p1 = ggplot(hpi_01, aes(x=ratio_dt, y=family, color = class)) + 
    geom_segment(aes(x = ratio_dt, xend=0, yend=family)) +
    geom_vline(xintercept=0) + 
    geom_point(size =3) + 
    facet_wrap(~syncom) + 
    main_theme + 
    theme(legend.position = "right") + 
    xlab("log2 of fold change of host preference index") + 
    theme(panel.grid.major = element_line(colour="grey90", size=0.5)) + 
    guides(color=guide_legend(title="Taxonomic class")) + 
    ylab("Taxonomic family")+
    theme(strip.text = element_text(face = "italic")) + 
    scale_x_continuous(limits = c(plt_min, plt_max))

f_name = "master_hpi_fold_change"
save_fig(p1, f_name, width=8, height=7)



# Re-isolation of limiting dilutions
idx = design$compartment == "collection"
design_subset = design[idx, ]
otu_subset = t(otu_table_norm)[idx, ]

otu_melt = melt(otu_subset)
colnames(otu_melt) = c("sample_id", "strain", "RA")

# at sphere and lj specific plots
otu_melt = otu_melt %>% left_join(design, by = c("sample_id" = "#SampleID"))
otu_melt$syncom = gsub("At", "at", otu_melt$syncom)
otu_melt$syncom = gsub("Lj", "lj", otu_melt$syncom)


# at and lf
otu_melt_tax = otu_melt %>% 
  left_join(taxonomy) %>%
  filter(!(sphere == "at" & syncom == "lj")) %>% # only correctly sequenced results
  filter(!(sphere == "lj" & syncom == "at")) %>% # only correctly sequenced results
  filter(!(RA == 0)) %>% # remove everything with RA == 0
  filter(!str_detect(sample_name, "nodule")) %>% # no nodule samples because why?
  filter(Plate == 2) #only second plate --> more reliable


otu_melt_tax$strain = reorder(otu_melt_tax$strain, -otu_melt_tax$RA, median)
otu_melt_tax$family = reorder(otu_melt_tax$family, -otu_melt_tax$RA, median)

coloring = data.frame(strain = levels(otu_melt_tax$strain))
coloring = coloring %>% left_join(taxonomy) %>% left_join(fam_colors)

otu_melt_tax$syncom[otu_melt_tax$syncom == "lj"] = "Lj-SPHERE"
otu_melt_tax$syncom[otu_melt_tax$syncom == "at"] = "At-SPHERE"

lwd = 0.8
point_size = 0.4

# plot that stuff
p1 = ggplot(data = otu_melt_tax, aes(x= family, y = RA, color = syncom)) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.6, lwd = lwd) + 
  main_theme + 
  geom_point(position=position_jitterdodge(jitter.width = 0.1), size = point_size) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "right") + 
  scale_y_continuous(trans='log10') + 
  xlab("") + 
  ylab("Relative abundance") +
  scale_color_discrete(
          "SynCom",
          labels = c(expression(italic("At-SPHERE")), expression(italic("Lj-SPHERE"))))

f_name = "master_reisoation"
save_fig(p1, f_name, width=10, height=7)




# for at: 2000, 18000
plate_1 = c("4_6",
            "13_15",
            "22_24",
            "32_33",
            "40_42",
            "50_51")
plate_2 = c("1_3",
            "4_6",
            "11_12",
            "14_15",
            "20_21",
            "22_24",
            "29_30",
            "32_33",
            "38_39",
            "40_41",
            "47_48", # ??
            "50_51")
vec = unique(c(plate_1, plate_2))
vec = c(vec, stringi::stri_reverse(vec))
vec = paste0(vec, "*")
vec = paste0(vec, collapse = "|")

idx = str_detect(design$sample_name, vec) & design$Plate == 2  & !(str_detect(design$sample_name, "nodule") | str_detect(design$sample_name, "pre")) & (design$compartment == "root" | design$compartment == "collection")
otu_list = by_condition(idx, otu_table_norm)
otu_tab_s = t(otu_table_norm)[idx, ]
design_subset = design[idx, ]

# PCoA Bray-Curtis
bray_curtis <- vegdist(otu_tab_s, method="bray")

k <- 2
pcoa <- cmdscale(bray_curtis, k=k, eig=T)
points <- pcoa$points
eig <- pcoa$eig
points <- as.data.frame(points)
colnames(points) <- c("x", "y")

points <- cbind(points, design[match(rownames(points), design$`#SampleID`), ])

colors_pcoa <- colors[colors$group %in% points$plant, ]
points$plant <- factor(points$plant, levels=colors_pcoa$group)
shapes_pcoa <- shapes[shapes$group %in% points$compartment, ]
points$compartment <- factor(points$compartment, levels=shapes_pcoa$group)
points$Plate = as.character(points$Plate)

pcoa_alpha = 0.6
pcoa_size = 2

points$compartment = as.character(points$compartment)
points$compartment[points$compartment == "collection"] = "Limiting dilution"
points$compartment[points$compartment == "root"] = "Root"
shapes_pcoa$group[shapes_pcoa$group == "collection"] = "Limiting dilution"
shapes_pcoa$group[shapes_pcoa$group == "root"] = "Root"

p = ggplot(points, aes(x=x, y=y, color=syncom, shape=compartment)) +
  geom_point(alpha=pcoa_alpha, size=pcoa_size) +
  scale_colour_manual(values=as.character(colors_pcoa$color)) +
  scale_shape_manual(values=rev(shapes_pcoa$shape)) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
  main_theme +
  theme(legend.position="right") +
  scale_color_discrete(
          "SynCom",
          labels = c(expression(italic("At-SPHERE")), expression(italic("Lj-SPHERE"))))


f_name = "master_reisoation_pcoa"
save_fig(p, f_name, width=7, height=5)


# investigate points that are a little off

##################
#CPCOA
idx = str_detect(design$sample_name, vec) & design$Plate == 2  & !(str_detect(design$sample_name, "nodule") | str_detect(design$sample_name, "pre")) & (design$compartment == "root" | design$compartment == "collection")
otu_list = by_condition(idx, otu_table_norm)
otu_tab_s = t(otu_table_norm)[idx, ]
design_subset = design[idx, ]

bray_curtis <- vegdist(otu_tab_s, method="bray")
capscale.gen <- capscale(bray_curtis ~ syncom*compartment*plant, data=design_subset, add=F, sqrt.dist=T)

# ANOVA-like permutation analysis
perm_anova.gen <- anova.cca(capscale.gen, permutations = 9999)
print(perm_anova.gen)

# generate variability tables and calculate confidence intervals for the variance
options(scipen = 999)
var_tbl.gen <- variability_table(capscale.gen)

eig <- capscale.gen$CCA$eig

variance <- capscale.gen$CCA$tot.chi / capscale.gen$tot.chi

variance <- var_tbl.gen["constrained", "proportion"]
p.val <- perm_anova.gen[1, 4]

points_cpcoa <- capscale.gen$CCA$wa[, 1:2]
colnames(points_cpcoa) <- c("x", "y")
points_cpcoa <- cbind(points_cpcoa, design_subset[match(rownames(points_cpcoa), design_subset$`#SampleID`), ])

# plot CPCo 1 and 2
colors_cpcoa <- colors[colors$group %in% points_cpcoa$plant, ]
points_cpcoa$genotype <- factor(points_cpcoa$plant, levels=colors_cpcoa$group)
shapes_cpcoa <- shapes[shapes$group %in% points_cpcoa$compartment, ]
points_cpcoa$compartment <- factor(points_cpcoa$compartment, levels=shapes_cpcoa$group)
points_cpcoa$Plate = as.character(points_cpcoa$Plate)

points_cpcoa$compartment = as.character(points_cpcoa$compartment)
points_cpcoa$compartment[points_cpcoa$compartment == "collection"] = "Limiting dilution"
points_cpcoa$compartment[points_cpcoa$compartment == "root"] = "Root"
shapes_pcoa$group[shapes_pcoa$group == "collection"] = "Limiting dilution"
shapes_pcoa$group[shapes_pcoa$group == "root"] = "Root"

labels = c(expression(italic("At-SPHERE")), expression(italic("Lj-SPHERE")))

p <- ggplot(points_cpcoa, aes(x=x, y=y, color=syncom, shape=compartment)) +
  geom_point(alpha=pcoa_alpha, size=pcoa_size) +
  scale_colour_manual("SynCom", values=as.character(colors_cpcoa$color), labels = labels) +
  theme(axis.text.x=element_text(size=12)) +
  theme(axis.text.y=element_text(size=12)) +
  theme(axis.title=element_text(size=13)) +
  labs(x=paste("CPCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("CPCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) +
  ggtitle(paste(format(100 * variance, digits=3), " % of variance; p-value = ", format(p.val, digits=5), sep="")) +
  main_theme +
  theme(legend.position="right") +
  scale_shape_discrete("Sample")

f_name = "masterReisolationCPCoA"
save_fig(p, f_name, width=7, height=5)


