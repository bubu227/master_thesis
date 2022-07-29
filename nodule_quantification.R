
# clean-up of data
rm(list = ls())

# load packages
library(tidyverse)
library(readxl)
library(stringr)
library(ggplot2)
library(rstatix)
library(rcompanion)
library(dunn.test)
library(rcompanion)

source("/home/niklas/mount/netscratch/dep_psl/grp_rgo/kieln/plotting_functions/plotting_functions.R")

# define paths
fw_excel = "../tables/lotus_nodules.xlsx"
mapping_file = "../tables/mapping_#5.xlsx"
ppp_file = "../processed_data/plants_per_pot.xlsx"

#create figure folder if not existing
fig_folder = "../figures/"
dir.create(fig_folder, showWarnings = FALSE)

# load the raw data as excel files
fw_table = read_excel(fw_excel)
mapping = read_excel(mapping_file)
ppp_table = read_excel(ppp_file)

# add missing value
lj6 = data.frame("ID" = "Lj6", "plants_number" = 2)
ppp_table = bind_rows(ppp_table, lj6)

# preprocess the dataframe
idx = str_starts(colnames(fw_table), "nodule")
nodules = fw_table[, idx]

# create dataframe with nodule number per flowpot
total_nodules = rowSums(nodules,na.rm=T)
per_flowpot = data.frame("ID" = fw_table$ID, nodule_number = total_nodules)

# append data together
mapping$SynCom = gsub("At", "At-SPHERE", mapping$SynCom)
mapping$SynCom = gsub("Lj", "Lj-SPHERE", mapping$SynCom)
mapping$species[mapping$species == "Lj"] = "Gifu"

per_flowpot = per_flowpot %>% left_join(mapping, by = "ID") %>% left_join(ppp_table)

per_flowpot$nodules_per_plant = per_flowpot$nodule_number / per_flowpot$plants_number


# coloring
group.colors = c("-" = "#DDDDDD", 
                 "At-SPHERE" = "#EE6677",
                 "Lj-SPHERE" = "#66CCEE")

Tol_light <- c('#BBCC33', 
               '#AAAA00', 
               '#77AADD', 
               '#EE8866', 
               '#EEDD88', 
               '#FFAABB', 
               '#99DDFF', 
               '#44BB99', 
               '#DDDDDD')

Tol_bright <- c('#EE6677', 
                '#228833', 
                '#4477AA', 
                '#CCBB44', 
                '#66CCEE', 
                '#AA3377', 
                '#BBBBBB')


# join two dataframes together
cols = c("ID", "species", "SynCom", "Media", "nodule_number", "nodules_per_plant")
per_flowpot = per_flowpot[, cols]

weight_table = per_flowpot

# remove mixed syncom experiments
idx = weight_table$SynCom == "mixed"
weight_table = weight_table[!idx, ]

# different grouping
groups = weight_table %>% group_by(SynCom, Media) %>% group_indices()
weight_table$group = as.character(groups)
kruskal_ind = weight_table %>% group_by(SynCom) %>% kruskal_test(nodules_per_plant ~ Media)

#extract only significant groups
idx = kruskal_ind$p < 0.05
significant.groups = kruskal_ind[idx, c("SynCom")]

treatments = expand.grid(species = unique(weight_table$species),
                         syncom = unique(weight_table$SynCom))

# post hoc test
species_wise_dunn = function(counter){
  syncom = significant.groups$SynCom[counter]
  idx = weight_table$SynCom == syncom
  subset = weight_table[idx, ]
  dt = dunn.test(subset$value, subset$group, method="bh")
  l = cldList(P.adjusted ~ comparisons, data=dt, threshold=0.05)
  return(l)
}

if (nrow(significant.groups) > 1){
  dt = lapply(c(1:nrow(significant.groups)), species_wise_dunn)
  l = bind_rows(dt)

  colnames(l) = c("group", "letter", "mono_letter")
  l

  # add statistical data to plot
  weight_table = weight_table %>% left_join(l, by="group")
  weight_table = weight_table %>% left_join(kruskal_ind) %>% 
    mutate(title = paste0(SynCom, "; p-val = ", p))
} else {
  weight_table = weight_table %>% left_join(kruskal_ind) %>% 
    mutate(title = paste0(SynCom, "; p-val = ", p))
  weight_table$letter = ""
}

ymax = round(max(weight_table$nodules_per_plant)*1.1)

p = ggplot(data = weight_table, 
       aes(x = Media, y = nodules_per_plant, fill=SynCom)) + 
  facet_grid(c("species", "SynCom"), scales = "fixed") + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(size = 0.5, width = 0.2) + 
  scale_fill_manual(values=group.colors) + 
  main_theme +
  scale_y_continuous(breaks= seq(0,ymax, 2), limits=c(0, ymax)) +
  ylab("Mean nodules per plant") +
  theme(strip.background = element_blank(),
        strip.text.y = element_text(size = 15),
        strip.text.x = element_text(face = "italic")) + 
  geom_text(
    data = weight_table %>% group_by(group, letter, Media, SynCom, species) %>% 
      summarise(y = max(nodules_per_plant)* 1.05) %>% ungroup(),
    aes(x = Media, y = y, label = letter), 
    size = 4,  
    parse = TRUE) + 
  geom_text(data = kruskal_ind  %>% 
            mutate(y = max(weight_table$nodules_per_plant) * 1.1),
            aes(x = 2, y = y, label = paste0("p-val: ", p)))

fig_name = paste0(fig_folder, "./nodule_quant.pdf")
ggsave(file = fig_name, plot = p, device = "pdf",
       width=10, height = 8)
fig_name = paste0(fig_folder, "./nodule_quant.png")
ggsave(file = fig_name, plot = p, width=10, height = 8)


################# all syncoms in a bin
kruskal_ind = weight_table %>% kruskal_test(nodules_per_plant ~ Media)

ymax = round(max(weight_table$nodules_per_plant)*1.1)

letter_df = weight_table %>% group_by(Media, letter) %>% summarise(y = max(nodules_per_plant)* 1.05)

p = ggplot(data = weight_table, 
       aes(x = Media, y = nodules_per_plant)) + 
  facet_grid(c("species"), scales = "fixed") + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(size = 0.5, width = 0.2) + 
  scale_fill_manual(values=group.colors) + 
  theme_bw() +
  scale_y_continuous(breaks= seq(0,ymax, 2), limits=c(0, ymax)) +
  ylab("Mean nodules per plant") +
  theme(strip.text = element_text(size = 12),
        axis.text.x = element_text(size = 10)) + 
  geom_text(
    data = weight_table %>% group_by(group, letter, Media, SynCom, species) %>% 
      summarise(y = max(nodules_per_plant)* 1.05) %>% ungroup(),
    aes(x = Media, y = y, label = letter), 
    size = 4,  
    parse = TRUE) + 
  geom_text(data = kruskal_ind  %>% 
            mutate(y = max(weight_table$nodules_per_plant) * 1.1),
            aes(x = 2, y = y, label = paste0("p-val: ", p)))



fig_name = paste0(fig_folder, "./nodule_quant_collapsed.pdf")
ggsave(file = fig_name, plot = p, device = "pdf",
       width=10, height = 8)
