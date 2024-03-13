library(cowplot)
library(gridExtra)

# FIGURE 2:
#################################################
# Plot A
set.seed(5)
plotA <- plotPathClusters(
  enrichment = go_enrich@result,
  sim = clusters$similarity,
  clusters = clusters$clusters,
  fontSize = 8,
  outerCutoff = 0.1,
  drawEllipses = TRUE,
  nodeSize = "Count",
  colorBy = "pvalue",
  colorType = 'pval',
  repelLabels = TRUE
)

# Plot B
set.seed(4)
plotB <- plotPathClusters(
  enrichment = go_enrich_no_it@result,
  sim = clusters_no_it$similarity,
  clusters = clusters_no_it$clusters,
  fontSize = 8,
  outerCutoff = 0.1,
  drawEllipses = TRUE,
  nodeSize = "Count",
  colorBy = "pvalue",
  colorType = 'pval',
  repelLabels = TRUE
)

# Plot C
plotC <- ggplot(df.m, aes(x = group, y = value, fill = group)) +
  geom_boxplot() +
  labs(x = "group", y = "Expression Level") +
  theme_bw() +
  theme(legend.position = "none") +
  facet_wrap(~ variable, scales = "free", ncol = 8) +
  geom_text(
    data = df.m,
    aes(x = 1, y = max(value), 
        label = paste("p-value =", p_value)),
    hjust = 0.2, vjust = 1, size = 5.5
  )


multi_panel_figure <- cowplot::plot_grid(
  cowplot::plot_grid(plotA, plotB, nrow =  1, ncol =  2, labels = c("A)", "B)"), label_size =  20),
  cowplot::plot_grid(NULL, plotC, NULL, nrow =  1, labels = c("C)"), label_size =  20, label_x =  0, rel_widths = c(0.5, 7, 0.5)),
  nrow =  2)

save_plot(
  "C://Users/Gerard/Desktop/Figure2_12.png",
  multi_panel_figure,
  base_width = 22,
  base_height = 17,
  bg = "white",
  dpi = 300
)
#################################################



# FIGURE 3:
#################################################
plotA <-
  ggplot(df, aes(x = "", y = count.Freq, fill = count.event)) +
  geom_bar(stat = "identity", color = "black", size = 0.8) +
  labs(title = "") +
  coord_polar("y", start = 80) +
  geom_text(aes(label = paste0(round(count.Freq, 2), '%')),
            position = position_stack(vjust = 0.5),
            size = 6.5) +
  scale_fill_brewer(palette = "Pastel1",
                    name = "",
                    label = new_labels) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    legend.position = "bottom",
    legend.box.spacing = unit(-3, "lines"),
    legend.text = element_text(size = 21)) +
  guides(fill=guide_legend(nrow=2, byrow=TRUE)) 


plotA


png("C://Users/Gerard/Desktop/Figure3.png", units = "in", width = 22, height = 22, res = 300)
multi <- grid.arrange(
  plotA, v.splicing, v.dia,
  layout_matrix = rbind(c(1, 2), c(3, NA)))

grid.text("ASAH1", x = 0.89, y = 0.77, gp = gpar(fontface = "italic", fontsize = 22, col = "black", fontfamily = "serif"))
grid.text("CYCS", x = 0.89, y = 0.75, gp = gpar(fontface = "italic", fontsize = 22, col = "black", fontfamily = "serif"))
grid.text("HMGB1", x = 0.89, y = 0.73, gp = gpar(fontface = "italic", fontsize = 22, col = "black", fontfamily = "serif"))
grid.text("SELENOP", x = 0.89, y = 0.71, gp = gpar(fontface = "italic", fontsize = 22, col = "black", fontfamily = "serif"))
dev.off()

multi <- cowplot::plot_grid(
  cowplot::plot_grid(plotA, v.splicing, nrow =  1, ncol =  2, labels = c("A)", "B)"), label_size =  24),
  cowplot::plot_grid(NULL, v.dia, NULL, nrow =  1, labels = c("C)"), label_size =  24, label_x =  1, rel_widths = c(0.5, 1, 0.5)),
  nrow =  2)


png("C://Users/Gerard/Desktop/Figure3.png", units = "in", width =  22, height =  22, res =  300)
multi
grid.text("ASAH1", x =  0.89, y =  0.77, gp = gpar(fontface = "italic", fontsize =  23, col = "black", fontfamily = "serif"))
grid.text("CYCS", x =  0.89, y =  0.75, gp = gpar(fontface = "italic", fontsize =  23, col = "black", fontfamily = "serif"))
grid.text("HMGB1", x =  0.89, y =  0.73, gp = gpar(fontface = "italic", fontsize =  23, col = "black", fontfamily = "serif"))
grid.text("SELENOP", x =  0.89, y =  0.71, gp = gpar(fontface = "italic", fontsize =  23, col = "black", fontfamily = "serif"))

grid.text("FHL1", x = 0.75, y = 0.8, gp = gpar(fontface = "italic", fontsize = 23, col = "black", fontfamily = "serif"))
grid.text("GNAS", x = 0.75, y = 0.78, gp = gpar(fontface = "italic", fontsize = 23, col = "black", fontfamily = "serif"))
grid.text("SPP1", x = 0.75, y = 0.76, gp = gpar(fontface = "italic", fontsize = 23, col = "black", fontfamily = "serif"))
grid.text("ARL1", x = 0.75, y = 0.74, gp = gpar(fontface = "italic", fontsize = 23, col = "black", fontfamily = "serif"))
grid.text("MORF4L2", x = 0.75, y = 0.72, gp = gpar(fontface = "italic", fontsize = 23, col = "black", fontfamily = "serif"))
grid.text("HMGN1", x = 0.75, y = 0.7, gp = gpar(fontface = "italic", fontsize = 23, col = "black", fontfamily = "serif"))
grid.text("RNASE4", x = 0.75, y = 0.68, gp = gpar(fontface = "italic", fontsize = 23, col = "black", fontfamily = "serif"))

grid.text("EXTL3", x = 0.5, y = 0.3, gp = gpar(fontface = "italic", fontsize = 23, col = "black", fontfamily = "serif"))
grid.text("ZFR", x = 0.5, y = 0.28, gp = gpar(fontface = "italic", fontsize = 23, col = "black", fontfamily = "serif"))
grid.text("VPS37C", x = 0.5, y = 0.26, gp = gpar(fontface = "italic", fontsize = 23, col = "black", fontfamily = "serif"))
grid.text("ZNF784", x = 0.5, y = 0.24, gp = gpar(fontface = "italic", fontsize = 23, col = "black", fontfamily = "serif"))
grid.text("DUSP8", x = 0.5, y = 0.22, gp = gpar(fontface = "italic", fontsize = 23, col = "black", fontfamily = "serif"))
grid.text("USP33", x = 0.5, y = 0.20, gp = gpar(fontface = "italic", fontsize = 23, col = "black", fontfamily = "serif"))
grid.text("RFX1", x = 0.5, y = 0.18, gp = gpar(fontface = "italic", fontsize = 23, col = "black", fontfamily = "serif"))
grid.text("DISP1", x = 0.5, y = 0.16, gp = gpar(fontface = "italic", fontsize = 23, col = "black", fontfamily = "serif"))
dev.off()
#################################################


aaa2 <- ifelse(aaa == 1, "AAA", "Control")

# FIGURE 4:
#################################################
plotA <- fviz_pca_ind(
  my_pca,
  habillage = aaa2,
  palette = c("#00AFBB", "#FC4E07"),
  repel = TRUE,
  addEllipses = TRUE,
  ellipse.level = 0.95,
  geom = "point"
) + labs(title = "") +
  theme(axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30),
        legend.text = element_text(size = 30),
        legend.title = element_blank(),
        aspect.ratio = 1)


set.seed(28)
plotB <- plotPathClusters(
  enrichment = go_enrich@result,
  sim = clusters$similarity,
  clusters = clusters$clusters,
  fontSize = 16,
  outerCutoff = 0.1,
  drawEllipses = TRUE,
  nodeSize = "Count",
  colorBy = "pvalue",
  colorType = 'pval',
  repelLabels = TRUE
) +
  theme(legend.text = element_text(size = 24),
        legend.title = element_text(size = 26))

multi_panel_figure <- cowplot::plot_grid(
  plotA, plotB, nrow =  1, ncol =  2, labels = c("A)", "B)"), label_size =  26) + 
  theme(plot.margin = margin(t =  0, r =  10, b =  10, l =  10, "mm"))


multi_panel_figure <- cowplot::plot_grid(
  plotA, nrow =  1, ncol =  1, labels = c("A)"), label_size =  26)

multi_panel_figure <- cowplot::plot_grid(
  plotA, label_size =  26)

multi_panel_figure <- cowplot::plot_grid(
  plotB, label_size =  26)

multi_panel_figure <- cowplot::plot_grid(
  plotA, plotB, labels = c("A)", "B)"), label_y =  0.8, label_size =  26)

save_plot(
  "C://Users/Gerard/Desktop/Figure4.png",
  multi_panel_figure,
  base_width = 30,
  base_height = 26,
  bg = "white",
  dpi = 250
)
#################################################





