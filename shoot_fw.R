
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
library(writexl)

# plotting theme
source("/home/niklas/mount/netscratch/dep_psl/grp_rgo/kieln/plotting_functions/plotting_functions.R")

# define paths
fw_excel = "../tables/shoot_fw_5.xlsx"
mapping_file = "../tables/mapping_#5.xlsx"

#create figure folder if not existing
fig_folder = "../figures/"
processed_folder = "../processed_data/"
dir.create(fig_folder, showWarnings = FALSE)
dir.create(processed_folder, showWarnings = FALSE)

# load the raw data as excel files
fw_table = read_excel(fw_excel)
mapping = read_excel(mapping_file)

# preprocess the dataframe
idx = str_starts(colnames(fw_table), "weight")
weights = fw_table[, idx]

#remove parts of mapping that are not required
idx = c(1:58)
mapping = mapping[idx, ]

mapping$SynCom[mapping$SynCom == "At"] = "At-SPHERE"
mapping$SynCom[mapping$SynCom == "Lj"] = "Lj-SPHERE"
mapping$species[mapping$species == "At"] = "Col-0"
mapping$species[mapping$species == "Lj"] = "Gifu"

# substract first column from second and so on
weights <- weights[ , colSums(is.na(weights)) < nrow(weights)]

rolled_df = weights[, c(1:ncol(weights)-1)]
corrected_weights = weights[, 2:ncol(weights)] - rolled_df
corrected_weights$weight_1 = rolled_df$weight_1

#save the intermediate table to disk
plants_per_pot = rowSums(!is.na(corrected_weights))
intermediate_res = data.frame("ID" = fw_table$ID, "plants_number"= plants_per_pot)
file_name = paste0(processed_folder, "plants_per_pot.xlsx")
write_xlsx(intermediate_res, file_name)

# calculate the mean
means = lapply(c(1:nrow(corrected_weights)), function(x) mean(as.numeric(corrected_weights[x,]), na.rm = TRUE))
corrected_weights$ID = fw_table$ID
corrected_weights$mean_fw = unlist(means)


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


# now with single measurements
corrected_weights = corrected_weights %>% select(-mean_fw)
cols_to_pivot = colnames(corrected_weights)[str_detect(colnames(corrected_weights), "weight")]

#pivot longer
weight_table = pivot_longer(corrected_weights,
             cols = all_of(cols_to_pivot))


# remove na measurements
idx = is.na(weight_table$value)
weight_table = weight_table[!idx, ]

# join two dataframes together
reduced_mapping = mapping[, c("number", "ID", "species", "SynCom", "Media")]
weight_table = weight_table %>% 
  left_join(reduced_mapping, by = "ID")


# remove mixed syncom experiments
idx = weight_table$SynCom == "mixed"
weight_table = weight_table[!idx, ]

# different grouping
groups = weight_table %>% group_by(SynCom, species, Media) %>% group_indices()
weight_table$group = as.character(groups)
kruskal_ind = weight_table %>% group_by(SynCom, species) %>% kruskal_test(value ~ Media)

#extract only significant groups
idx = kruskal_ind$p < 0.05
significant.groups = kruskal_ind[idx, c("species", "SynCom")]

treatments = expand.grid(species = unique(weight_table$species),
                         syncom = unique(weight_table$SynCom))

# post hoc test
species_wise_dunn = function(counter){
  species = significant.groups$species[counter]
  syncom = significant.groups$SynCom[counter]
  idx = weight_table$species == species & weight_table$SynCom == syncom
  subset = weight_table[idx, ]
  dt = dunn.test(subset$value, subset$group, method="bh")
  l = cldList(P.adjusted ~ comparisons, data=dt, threshold=0.05)
  return(l)
}

dt = lapply(c(1:nrow(significant.groups)), species_wise_dunn)
l = bind_rows(dt)

colnames(l) = c("group", "letter", "mono_letter")
l

# add statistical data to plot
weight_table = weight_table %>% left_join(l, by="group")
weight_table = weight_table %>% left_join(kruskal_ind) %>% 
  mutate(title = paste0(SynCom, "; p-val = ", p))

# renaming of stuff
weight_table_p = weight_table


p = ggplot(data = weight_table_p, 
       aes(x = Media, y = value, fill=SynCom)) + 
  facet_grid(c("species", "SynCom"), scales = "free") + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(size = 0.5, width = 0.2) + 
  scale_fill_manual(values=group.colors) + 
  theme_bw() +
  ylab("Shoot fresh weight [mg]") +
  main_theme + 
  theme(strip.background = element_blank(),
        strip.text.y = element_text(size = 15),
        strip.text.x =element_text(face = "italic")) + 
  geom_text(
    data = weight_table_p %>% group_by(group, letter, Media, SynCom, species) %>% 
      summarise(y = max(value)* 1.1) %>% ungroup(),
    aes(x = Media, y = y, label = letter), 
    size = 4,  
    parse = TRUE) + 
  geom_text(data = weight_table_p %>% group_by(species) %>% 
              summarise(y = max(value) * 1.2) %>% ungroup() %>% 
              left_join(kruskal_ind),
            aes(x = 2, y = y, label = paste0("p-val: ", p)))
p

fig_name = paste0(fig_folder, "./shoot_fw.pdf")
ggsave(file = fig_name, plot = p, device = "pdf",
       width=10, height = 8)
fig_name = paste0(fig_folder, "./shoot_fw.png")
ggsave(file = fig_name, plot = p, width=10, height = 8)


