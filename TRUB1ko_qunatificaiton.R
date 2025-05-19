# packages and global parameters, function
library(tidyverse)
library(ggrepel)
local <- '/Users/shiyuanguo/Library/CloudStorage/GoogleDrive-sguo039@ucr.edu/My Drive/PhD_study/TRUB1 pseudoU crosstalk/TRUB1 manuscript'

# linear equation compatible with ggplot
lm_equation <- function(y, x, data, ...){
  # annotating the equation on ggplot. 
  formula <- as.formula(paste(y,paste(x, collapse = ' + '),sep = ' ~ '))
  m <- lm(formula, data, ...) 
  substituteList <- list(a = format(abs(unname(coef(m)[1])), digits = 3),
                         b = format(unname(coef(m)[2]), digits = 3),
                         r2 = format(summary(m)$r.squared, digits = 3))
  if (coef(m)[1] > 0){
    eq <- substitute(italic(y) == b %.% italic(x)+a*";"~~italic(r)^2~"="~r2, substituteList)
  } else {
    eq <- substitute(italic(y) == b %.% italic(x)-a*";"~~italic(r)^2~"="~r2, substituteList)
  }
  as.character(as.expression(eq))
}


# quantification of RWE proteomics -----
# obtain the injection info from both replicates
injectionfiles = paste0(file.path(local, 'data', 'methods_43txt_3', '3methods_3txts'), '_000', 1:3, '.csv')
injectionList <- lapply(injectionfiles, function(n){
  readr::read_csv(n) %>%
    dplyr::mutate(Comment = stringr::str_remove_all(Comment,' \\((light|heavy)\\)') %>% 
                    stringr::str_remove_all('\\[[:graph:][:digit:]{2}\\.[:digit:]{6}\\]')) %>%
    .[['Comment']] %>% 
    unique(.)
})

# 1) loading SILAC_peptidesRatio_TRUB1KO 
# skyline document grid > SILAC_peptideRatio > export as 'SILAC_peptideRatio.csv'
# TRUB1 KO 293T vs ctrl 
# F = H:TRUB1KO + L:293T; R = L:TRUB1KO + H:293T
ptbl <- read_csv(file.path(local, 'data', 'SILAC_peptideRatio_TRUB1KO.csv')) %>% 
  filter(Protein != 'sp|P02769|ALBU_BOVIN-standard') %>%
  ############ filter replicate base on injection number.
  mutate(prmlist = case_when(
    Peptide %in% injectionList[[1]] ~ 1,
    Peptide %in% injectionList[[2]] ~ 2,
    Peptide %in% injectionList[[3]] ~ 3,
    TRUE ~ 0)) %>%
  group_by(prmlist) %>%
  mutate(repli = str_extract(Replicate, '\\d$') %>% as.numeric(),
         fr = str_extract(Replicate, '[:alpha:]')) %>% 
  filter(repli == prmlist) %>% 
  ungroup() %>%
  ############
  mutate(`Total Area` = as.numeric(`Total Area`))%>% 
  select(Peptide, `Protein Name`, fr, `Total Area`, `Isotope Label Type`) %>% 
  pivot_wider(names_from = c(fr, `Isotope Label Type`), values_from = `Total Area`)%>%
  mutate(Fwd =   F_heavy / F_light, 
         Rvs =   R_light / R_heavy) %>% 
  mutate(`Protein Name` = if_else(grepl('ENSP00000409983', `Protein Name`), paste('IGF2BP3', `Protein Name`, sep = ' '), `Protein Name`)) %>% 
  mutate(`Protein Name` = if_else(grepl('ENSP00000384369', `Protein Name`), paste('METTL15', `Protein Name`, sep = ' '), `Protein Name`)) %>% 
  separate_wider_regex('Protein Name', c(proteinName = ".*?", " |/", ".*"))

# 2) protein level (w\o normalization) 
prottbl <- ptbl %>%
  group_by(proteinName) %>% # using `group_by` followed by `summarise` instead of `nest` to get two list-col for fwd and rvs respectively. 
  summarise(across(matches('Fwd|Rvs'), ~list(.x))) %>% 
  ungroup() %>%
  mutate(across(matches('Fwd|Rvs'), ~map_dbl(., ~mean(.x, na.rm = TRUE)))) %>% 
  rowwise() %>% 
  mutate(log2fc = log2(mean(c(Fwd, Rvs), na.rm = TRUE))) %>% 
  arrange(desc(log2fc))

# 3) output peptide table
th <- 1.5 # include YRDC
gt_ls <- ptbl %>% filter(Fwd > th & Rvs > th) %>% pull(proteinName)
gt_tbl <- ptbl %>% 
  group_by(proteinName) %>% 
  filter(proteinName %in% gt_ls) %>% 
  ungroup()
lt_ls <- prottbl %>% filter(Fwd < 1/th & Rvs < 1/th) %>% pull(proteinName)
lt_tbl <- ptbl %>% 
  group_by(proteinName) %>% 
  filter(proteinName %in% lt_ls) %>% 
  ungroup()
wb <- openxlsx::createWorkbook(); 
tn <- 'PRM_RWE_fulllist'; openxlsx::addWorksheet(wb, tn); openxlsx::writeData(wb, tn, ptbl);
tn <- 'gt'; openxlsx::addWorksheet(wb, tn); openxlsx::writeData(wb, tn, gt_tbl); 
tn <- 'lt'; openxlsx::addWorksheet(wb, tn); openxlsx::writeData(wb, tn, lt_tbl); 
openxlsx::openXL(wb)

# 4) output protein table arranged in desc order protein level fc
pplist <- ptbl %>%
  rowwise() %>% 
  mutate(mean_peptidelvl = mean(c(Fwd, Rvs), na.rm = TRUE),
         log2fc_peptidelvl = log2(mean_peptidelvl)) %>% 
  ungroup() %>% 
  group_by(proteinName) %>% # using `group_by` followed by `summarise` instead of `nest` to get two list-col for fwd and rvs respectively. 
  mutate(across(matches('Fwd|Rvs'), ~list(.x))) %>% 
  ungroup() %>%
  rowwise() %>% 
  mutate(mean_proteinlvl = mean(c(Fwd, Rvs), na.rm = TRUE),
         log2fc_proteinlvl = log2(mean_proteinlvl)) %>% 
  dplyr::arrange(desc(log2fc_proteinlvl)) %>% 
  ungroup() %>%
  group_by(proteinName) %>%
  nest() %>%
  ungroup() %>%
  mutate(rnum = row_number()) %>%
  mutate(shade = case_when(
    rnum %% 2 == 1 ~ 0,
    rnum %% 2 == 0 ~ 1,
  )) %>% 
  unnest(data) %>% 
  dplyr::select(-Fwd, -Rvs)
wb <- openxlsx::createWorkbook(); 
openxlsx::addWorksheet(wb, 'test'); 
lsty1 <-  openxlsx::createStyle(bgFill = "#d6d2d2")
openxlsx::conditionalFormatting(wb, 'test', rows = 1:(nrow(pplist)+1), cols = 1:ncol(pplist), rule = paste0("$", LETTERS[which(colnames(pplist) == "shade")], "1==", 1), style = lsty1)
openxlsx::writeData(wb, 'test', pplist); 
openxlsx::openXL(wb)

# peptide level plot
# cutoff <- 1.5
# lab_cutoff <- 1.5
# p <- ptbl %>% 
# ggplot(aes(x = Fwd, y = Rvs, group = proteinName))+
#   geom_point(color = 'grey', size = 1, alpha = 0.5)+
#   geom_point(data = . %>% filter((Fwd > cutoff & Rvs > cutoff) | (Fwd < 1/cutoff & Rvs < 1/cutoff)), color = 'red', size = 1)+
#   ggrepel::geom_label_repel(data = . %>% filter((proteinName %in% c('YRDC', 'METTL17')) | (Fwd > lab_cutoff & Rvs > lab_cutoff) | (Fwd < 1/lab_cutoff & Rvs < 1/lab_cutoff) ), aes(label = proteinName), box.padding = 0.5, max.overlaps = Inf)+
#   # geom_point(data = . %>% filter(nchar(zincFinger) !=0 & `num<1.5` <= 2 & mean > 1.5) %>% filter(gn %in% c('POLD1')), aes(color = `Arsenite replaceable zinc-finger domain`), size = 1)+
#   ggh4x::coord_axes_inside(labels_inside = TRUE, xlim = c(0, 2.5), ylim = c(0, 2.5), xintercept = 1, yintercept = 1, ratio = 1)+
#   # labs(x = "1/Fwd =  WT(L) / TRUB1KO(H)", y = "Rvs =  WT(H) / TRUB1KO(L)")+
#   labs(x = "Fwd = TRUB1KO(H) / WT(L)", y = "1/Rvs = TRUB1KO(L) / WT(H)")+
#   theme_classic()+ # increase x limit
#   theme(legend.position = "none")
# plotly::ggplotly(p)

# 5) protein ratio plot (log2 scale with normal ratios)
plot_shape = 20; plot_size = 0.5; plot_alpha = 0.4; plot_pos = '#F8766D'; global_linewidth <- 0.3 # ggploting parameters
# p <- prottbl %>% 
#   ggplot(aes(x = Fwd, y = Rvs, group = proteinName), shape = plot_shape)+
#   geom_point(color = 'grey', size = plot_size, alpha = plot_alpha)+
#   geom_point(data = . %>% filter(proteinName %in% c(gt_ls, lt_ls)), color = plot_pos, size = plot_size)+
#   ggrepel::geom_text_repel(data = . %>% filter(proteinName %in% c(gt_ls, tail(lt_ls, 6))), aes(label = proteinName), size= 3, force_pull = 10, force = 1, min.segment.length = 0.1, max.overlaps = Inf)+
#   scale_x_continuous(trans = scales::log2_trans(), breaks = c(.5,1 ,2))+
#   scale_y_continuous(trans = scales::log2_trans(),breaks = c(.5 ,2))+
#   ggh4x::coord_axes_inside(labels_inside = TRUE, xlim = c(0.4, 3), ylim = c(0.4, 3), xintercept = 0, yintercept = 0, ratio = 1)+
#   labs(x = "TRUB1KO(H) / WT(L)", y = "TRUB1KO(L) / WT(H)")+
#   theme_classic()+ # increase x limit
#   theme(legend.position = "none",
#         # panel.border = element_rect(colour = "black", fill=NA),
#         legend.box.background = element_rect(colour = "black"),
#         legend.background = element_blank())

# 5.1) normal scale with log2(ratios)
p <- prottbl %>% 
  mutate(logfwd = log2(Fwd), logrvs = log2(Rvs)) %>% 
  ggplot(aes(x = logfwd, y = logrvs, group = proteinName), shape = plot_shape)+
  geom_point(color = 'grey', size = plot_size, alpha = plot_alpha)+
  geom_point(data = . %>% filter(proteinName %in% c(gt_ls, lt_ls)), color = plot_pos, size = plot_size)+
  ggrepel::geom_text_repel(data = . %>% filter(proteinName %in% c(gt_ls, tail(lt_ls, 6))), aes(label = proteinName), size= 3, force_pull = 10, force = 1, min.segment.length = 0.1, max.overlaps = Inf)+
  scale_x_continuous(breaks = c(-1, 0, 1))+
  scale_y_continuous(breaks = c(-1, 1))+
  ggh4x::coord_axes_inside(labels_inside = TRUE, xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5), xintercept = 0, yintercept = 0, ratio = 1)+
  labs(x = "Log2(TRUB1KO(H) / WT(L))", y = "log2(TRUB1KO(L) / WT(H))")+
  theme_classic()+ # increase x limit
  theme(legend.position = "none",
        # panel.border = element_rect(colour = "black", fill=NA),
        legend.box.background = element_rect(colour = "black"),
        legend.background = element_blank())
# plotly::ggplotly(p)
ggsave(filename = file.path(local, 'figure', 'protein_dotplot.svg'), device = 'svg', plot = p, width = 3, height = 3, units = 'in')

# Gm calibration ----
# rG quantification from 250513; Gm quantification from 250422
Aratio_rG <- readxl::read_excel(file.path(local, 'data', 'calibration.xlsx'), sheet = '250513_rGcalibration_magicmagic', range = readxl::cell_rows(1:5)) %>% 
  mutate(labeled = ifelse(str_detect(analytes, '-'), 'labeled', 'unlab')) %>%
  mutate(analytes = case_when(
    str_detect(analytes, 'rG') ~'rG',
    str_detect(analytes, '(Gm|m22G)') ~ 'Gm'
  )) %>% # split name to `analytes` and `label`
  select(-IV193_5_20250514154638) %>% # weird peak shape 
  select(analytes, labeled, starts_with('IV193_')) %>%
  pivot_longer(where(is.double), names_to = 'sample', values_to = 'amt') %>%
  pivot_wider(names_from = 'labeled', values_from = 'amt') %>% 
  mutate(ratio = unlab / labeled) %>%
  filter(!is.na(ratio)) %>% 
  mutate(sample = str_extract(sample, '^([^_-]+_[^_-]+)?')) %>% # from gemini!! 
  group_by(analytes, sample) %>% 
  summarize(Aratio = mean(ratio, na.rm = TRUE), # Aratio for AUC ratio
            Aratiorsd = sd(ratio, na.rm = TRUE) / Aratio) %>% 
  ungroup() %>% 
  filter(analytes == 'rG')

Aratio_Gm <- readxl::read_excel(file.path(local, 'data', 'calibration.xlsx'), sheet = '250422_Gmcalibration_magicmagic', range = readxl::cell_rows(1:5)) %>% 
  mutate(labeled = ifelse(str_detect(analytes, '-'), 'labeled', 'unlab')) %>%
  mutate(analytes = case_when(
    str_detect(analytes, 'rG') ~'rG',
    str_detect(analytes, '(Gm|m22G)') ~ 'Gm'
  )) %>% # split name to `analytes` and `label`
  select(analytes, labeled, starts_with('IV193_')) %>%
  pivot_longer(where(is.double), names_to = 'sample', values_to = 'amt') %>%
  pivot_wider(names_from = 'labeled', values_from = 'amt') %>% 
  mutate(ratio = unlab / labeled) %>%
  filter(!is.na(ratio)) %>% 
  mutate(sample = str_extract(sample, '^([^_-]+_[^_-]+)?')) %>% # from gemini!! 
  group_by(analytes, sample) %>% 
  summarize(Aratio = mean(ratio, na.rm = TRUE), # Aratio for AUC ratio
            Aratiorsd = sd(ratio, na.rm = TRUE) / Aratio) %>% 
  ungroup() %>% 
  filter(analytes == 'Gm') 

range(Aratio_rG$Aratio)
range(Aratio_Gm$Aratio)

# Mratio for Amount ratio
rG_unlab = 0.785; rG_label = 0.165; # uM
Gm_unlab = 2.22; Gm_label = 47.8 # nM 
Mratio <- tibble(rG_unlab = c(rep(rG_unlab, 3), rep(rG_unlab * 10, 3), rep(rG_unlab * 100, 3)),
       Gm_unlab = c(rep(Gm_unlab, 3), rep(Gm_unlab * 10, 3), rep(Gm_unlab * 100, 3)),
       vol_unlab = rep(c(1,2,5), 3), 
       rG_label = rG_label,
       Gm_label = Gm_label,
       vol_label = 1 # 5 standard mix, add 5uL per sample 
       ) %>% 
  mutate(rG_unlab = rG_unlab * vol_unlab,
         Gm_unlab = Gm_unlab * vol_unlab,
         rG_label = rG_label * vol_label,
         Gm_label = Gm_label * vol_label) %>% 
  mutate(Mratio_rG = rG_unlab / rG_label,
         Mratio_Gm = Gm_unlab / Gm_label) %>% 
  pivot_longer(starts_with('Mratio'), names_to = 'analytes', names_pattern = 'Mratio_(.*)', values_to = 'Mratio') %>% 
  select(analytes, Mratio)

Gm_cal <- bind_cols(Aratio_Gm, Mratio %>% filter(analytes == 'Gm')) %>% 
  filter(sample != 'IV193_1')
Gmplot <- Gm_cal %>% 
  ggplot(aes(x = Mratio, y = Aratio)) +
  geom_point()+
  geom_smooth(method = 'lm', se = FALSE, formula = 'y~x', size = 0.5,alpha = 0.3) +
  annotate(geom="text", x = max(Gm_cal$Mratio)*0.05, y = max(Gm_cal$Aratio)*1.1, hjust = 0, label=lm_equation('Aratio', 'Mratio', data = Gm_cal), parse=TRUE)+
  ggh4x::coord_axes_inside(labels_inside = TRUE, xlim = c(0, max(Gm_cal$Mratio)*1.1), ylim = c(0, max(Gm_cal$Aratio)*1.1), xintercept = 0, yintercept = 0, ratio = max(Gm_cal$Mratio)/max(Gm_cal$Aratio) / 1.5)+ # ratio = max(Gm_cal$Mratio)/max(Gm_cal$Aratio)
  labs(x = "Amount (Gm / D6-m2,2G)", y = "Area (Gm / D6-m2,2G)")+
  theme_classic()+ # increase x limit
  theme(legend.position = "none",
        # panel.border = element_rect(colour = "black", fill=NA),
        legend.box.background = element_rect(colour = "black"),
        legend.background = element_blank())
# ggsave(filename = file.path(local, 'figure', 'GmCalibration.pdf'), device = 'pdf', plot = Gmplot, width = 5, height = 5 * 0.58, units = 'in')

rG_cal <- bind_cols(Aratio_rG, Mratio %>% filter(analytes == 'rG'))
rGplot <- rG_cal %>% 
  ggplot(aes(x = Mratio, y = Aratio)) +
  geom_point()+
  geom_smooth(method = 'lm', se = FALSE, formula = 'y~x', color = 'red', size = 0.5,alpha = 0.3) +
  annotate(geom="text", x = max(rG_cal$Mratio)*0.05, y = max(rG_cal$Aratio)*1.1, hjust = 0, label=lm_equation('Aratio', 'Mratio', data = rG_cal), parse=TRUE)+
  ggh4x::coord_axes_inside(labels_inside = TRUE, xlim = c(0, max(rG_cal$Mratio)*1.1), ylim = c(0, max(rG_cal$Aratio)*1.1), xintercept = 0, yintercept = 0, ratio = max(rG_cal$Mratio)/max(rG_cal$Aratio) / 1.5)+ # ratio = max(rG_cal$Mratio)/max(rG_cal$Aratio)
  labs(x = "Amount (rG / 15N5-rG)", y = "Area (rG / 15N5-rG)")+
  theme_classic()+ # increase x limit
  theme(legend.position = "none",
        # panel.border = element_rect(colour = "black", fill=NA),
        legend.box.background = element_rect(colour = "black"),
        legend.background = element_blank())
# ggsave(filename = file.path(local, 'figure', 'rGCalibration.pdf'), device = 'pdf', plot = rGplot, width = 5, height = 5 * 0.58, units = 'in')


# Gm quantification ----
Gmlm <- lm(Mratio ~ Aratio, data = Gm_cal)
rGlm <- lm(Mratio ~ Aratio, data = rG_cal) 
rG_label = 16.5; # uM
Gm_label = 4780 # nM 

ko <- readxl::read_excel(file.path(local, 'data', 'calibration.xlsx'), sheet = '250401_Gm_t6A', range = readxl::cell_rows(1:5)) %>% 
  mutate(labeled = ifelse(str_detect(analytes, '-'), 'labeled', 'unlab')) %>%
  mutate(analytes = case_when(
    str_detect(analytes, 'rG') ~'rG',
    str_detect(analytes, '(Gm|m22G)') ~ 'Gm'
  )) %>% # split name to `analytes` and `label`
  select(analytes, labeled, starts_with('IV189_')) %>%
  pivot_longer(where(is.double), names_to = 'sample', values_to = 'amt') %>%
  pivot_wider(names_from = 'labeled', values_from = 'amt') %>% 
  mutate(Aratio = unlab / labeled)
trub1ko <- ko %>% 
  mutate(Mratio = c(predict(rGlm, newdata = ko[ko$analytes == 'rG','Aratio']),
         predict(Gmlm, newdata = ko[ko$analytes == 'Gm','Aratio']))) %>% 
  mutate(unlab = case_when(
    analytes == 'Gm' ~ Mratio * Gm_label,
    analytes == 'rG' ~ Mratio * rG_label
    )) %>% 
  select(analytes,sample,unlab) %>% 
  pivot_wider(names_from = 'analytes',  values_from = 'unlab') %>% 
  mutate(`Gm/rG (ppm)` = Gm / (rG * 1000) * 1000000)

wb <- openxlsx::createWorkbook(); openxlsx::addWorksheet(wb, 'test'); openxlsx::writeData(wb, 'test', trub1ko); openxlsx::openXL(wb)

# p <- trub1ko %>%
#   mutate(trt = case_when(
#     str_detect(sample,'c') ~ 'HEK293T',
#     str_detect(sample,'k') ~ 'TRUB1 KO'
# )) %>% 
#   ggplot(aes(x = trt, y = `Gm/rG (ppm)`))+
#   stat_summary(geom = 'col', fun = function(x) mean(x))+
#   geom_point(aes(color = sample))
# plotly::ggplotly((p))


