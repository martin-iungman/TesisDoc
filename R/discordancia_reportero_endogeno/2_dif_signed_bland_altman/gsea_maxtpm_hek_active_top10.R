# Top 10 GSEA terms in each direction (NES), for the GSEA already
# computed in remap_tf_maxtpm_hek_active_gsea_ora.R (same folder) - the
# only one of the two new combinations with GSEA hits (30/30 significant).
# The earlier plot took top 20 by |NES|, which happened to be almost all
# negative (only ~6 positive terms exist among the 30 significant).

library(tidyverse)
library(readxl)

out_dir <- "figures/discordancia_reportero_endogeno/2_dif_signed_bland_altman"
gse <- read_excel(file.path(out_dir, "GSEA_remap_rank_dif_maxtpm_hek_active.xlsx"), sheet = "GSEA") %>%
  filter(p.adjust < 0.05)

top_pos <- gse %>% filter(NES > 0) %>% slice_max(order_by = NES, n = 10)
top_neg <- gse %>% filter(NES < 0) %>% slice_min(order_by = NES, n = 10)

message("Terminos positivos disponibles: ", sum(gse$NES > 0), " (de 10 pedidos)")
message("Terminos negativos disponibles: ", sum(gse$NES < 0), " (de 10 pedidos)")

p <- bind_rows(top_pos, top_neg) %>%
  mutate(
    dif = ifelse(NES > 0, "Reportero sobreestima", "Reportero subestima"),
    Description = fct_reorder(Description, NES)
  ) %>%
  ggplot(aes(x = NES, y = Description, fill = dif)) +
  geom_col() +
  ggpubr::theme_pubr(base_size = 12) +
  theme(legend.position = "top") +
  labs(fill = "Efecto", x = "NES", y = NULL,
       title = "Top 10 términos GO por NES en cada sentido (GSEA)",
       subtitle = "max_tpm + todo ChIP de ReMap, solo promotores activos en HEK293") +
  scale_fill_manual(values = c("Reportero sobreestima" = "#D6741F", "Reportero subestima" = "#7FB800"))

ggsave(file.path(out_dir, "GSEA_top10_cada_sentido_maxtpm_hek_active.jpg"), p, width = 10, height = 7, units = "in")
message("Saved to ", out_dir)
