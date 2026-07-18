## =========================
## 0. 加载 R 包
## =========================
library(sf)
library(dplyr)
library(ggplot2)
library(ggspatial)
library(patchwork)
library(viridis)
library(ggpubr)
library(readr)
library(lwgeom)
library(svglite)
library(ragg)

## =========================
## 1. 文件路径
## =========================
province_path <- "省界_Project.shp"
country_path  <- "国界_Project.shp"

result_rds <- "Huperzia_Spatial_Phylo_Master_Results.rds"
result_csv <- "Huperzia_Spatial_Phylo_Master_Results.csv"

out_prefix <- "Huperzia_Spatial_Phylo_Stable"

## =========================
## 2. 读取空间系统发育结果
## =========================
if (file.exists(result_rds)) {
  master_data <- readRDS(result_rds)
} else if (file.exists(result_csv)) {
  master_data <- read.csv(result_csv, stringsAsFactors = FALSE, check.names = FALSE)
} else {
  stop("未找到结果文件：Huperzia_Spatial_Phylo_Master_Results.rds 或 Huperzia_Spatial_Phylo_Master_Results.csv")
}

if (!all(c("Lon", "Lat") %in% names(master_data))) {
  if (all(c("lon", "lat") %in% names(master_data))) {
    master_data <- master_data %>% rename(Longitude = lon, Latitude = lat)
  } else if (all(c("x", "y") %in% names(master_data))) {
    master_data <- master_data %>% rename(Longitude = x, Latitude = y)
  } else {
    stop("master_data 中没有找到 Longitude / Latitude 列。请检查经纬度列名。")
  }
}

required_cols <- c("Lon", "Lat", "PE", "RPE", "Endemism_Type", "Mean_ED_Age")
missing_cols <- setdiff(required_cols, names(master_data))
if (length(missing_cols) > 0) {
  stop("master_data 缺少以下列：", paste(missing_cols, collapse = ", "))
}

endemism_levels <- c(
  "Neo-endemism",
  "Paleo-endemism",
  "Mixed-endemism",
  "Super-endemism",
  "Not Significant"
)

master_data <- master_data %>%
  mutate(
    Longitude = as.numeric(Lon),
    Latitude = as.numeric(Lat),
    PE = as.numeric(PE),
    RPE = as.numeric(RPE),
    Mean_ED_Age = as.numeric(Mean_ED_Age),
    Endemism_Type = factor(Endemism_Type, levels = endemism_levels)
  )

## =========================
## 3. 读取并修复 shp
## =========================
read_fix_sf <- function(path, encoding = "GBK", force_crs_if_missing = NA) {
  if (!file.exists(path)) {
    stop("文件不存在：", path)
  }

  x <- sf::st_read(path, quiet = TRUE, options = paste0("ENCODING=", encoding))

  if (is.na(sf::st_crs(x))) {
    if (is.na(force_crs_if_missing)) {
      stop(
        "底图缺少 CRS / prj 文件：", path,
        "\n请确认同目录下是否有 .prj 文件；如果没有，请在 read_fix_sf() 中设置 force_crs_if_missing。"
      )
    } else {
      sf::st_crs(x) <- force_crs_if_missing
    }
  }

  x <- x %>%
    sf::st_zm(drop = TRUE, what = "ZM") %>%
    sf::st_transform(4326) %>%
    sf::st_make_valid()

  ## 面数据用 buffer(0) 修复；线数据不要 buffer，否则会变成面
  geom_dim <- unique(as.integer(sf::st_dimension(x)))
  if (all(geom_dim == 2, na.rm = TRUE)) {
    x <- suppressWarnings(sf::st_buffer(x, 0))
    x <- sf::st_make_valid(x)
    x <- suppressWarnings(sf::st_collection_extract(x, "POLYGON"))
  }

  return(x)
}

china_province <- read_fix_sf(province_path, encoding = "GBK")
china_country_raw <- read_fix_sf(country_path, encoding = "GBK")

message("省界范围：")
print(sf::st_bbox(china_province))
message("国界范围：")
print(sf::st_bbox(china_country_raw))

## =========================
## 4. 投影与范围
## =========================
china_crs <- "+proj=aea +lat_1=25 +lat_2=47 +lon_0=105 +datum=WGS84 +units=m +no_defs"

main_xlim <- c(73, 135)
main_ylim <- c(17, 55)

south_china_sea_xlim <- c(105, 125)
south_china_sea_ylim <- c(3, 25)

grid_width  <- 0.2
grid_height <- 0.2

## =========================
## 5. 构建 AI 更稳定的边界线
## =========================
china_province_proj <- china_province %>%
  sf::st_transform(china_crs) %>%
  sf::st_make_valid()

china_province_fill <- china_province_proj

china_province_boundary <- sf::st_sf(
  geometry = sf::st_boundary(sf::st_geometry(china_province_proj)),
  crs = sf::st_crs(china_province_proj)
)

china_country_outline <- suppressWarnings(sf::st_union(sf::st_geometry(china_province_proj)))
china_country_outline <- sf::st_make_valid(china_country_outline)
china_country_outline <- sf::st_sf(
  geometry = sf::st_boundary(china_country_outline),
  crs = sf::st_crs(china_province_proj)
)

## 国界原始文件只用于南海插图，保留其中可能存在的南海诸岛/九段线等信息
china_country_raw_proj <- china_country_raw %>%
  sf::st_transform(china_crs) %>%
  sf::st_make_valid()

## =========================
## 6. 统一主题
## =========================
theme_map_china <- function() {
  theme_bw(base_family = "Times New Roman") +
    theme(
      text = element_text(family = "Times New Roman"),
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      plot.tag = element_text(face = "bold", size = 14, family = "Times New Roman"),
      plot.tag.position = c(0.02, 0.98),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8, color = "black"),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      legend.key.height = unit(0.45, "cm"),
      legend.key.width = unit(0.30, "cm"),
      plot.margin = margin(5, 5, 5, 5)
    )
}

theme_age_plot <- function() {
  theme_bw(base_family = "Times New Roman") +
    theme(
      text = element_text(family = "Times New Roman"),
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      plot.tag = element_text(face = "bold", size = 14, family = "Times New Roman"),
      plot.tag.position = c(0.02, 0.98),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      axis.title.y = element_text(size = 9),
      legend.position = "none",
      panel.grid.major = element_line(color = "grey88", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      plot.margin = margin(5, 5, 5, 5)
    )
}

## =========================
## 7. 标准中国地图函数
## =========================
plot_china_grid <- function(data,
                            value_col,
                            title_text,
                            legend_title,
                            fill_scale,
                            tag_letter,
                            show_tile_na = FALSE) {

  plot_data <- data %>%
    filter(!is.na(Longitude), !is.na(Latitude))

  if (!show_tile_na) {
    plot_data <- plot_data %>% filter(!is.na(.data[[value_col]]))
  }

  p_main <- ggplot() +
    geom_sf(
      data = china_province_fill,
      fill = "white",
      color = NA
    ) +
    geom_tile(
      data = plot_data,
      aes(x = Longitude, y = Latitude, fill = .data[[value_col]]),
      width = grid_width,
      height = grid_height,
      alpha = 0.92
    ) +
    geom_sf(
      data = china_province_boundary,
      fill = NA,
      color = "grey60",
      linewidth = 0.18,
      lineend = "round",
      linejoin = "round"
    ) +
    geom_sf(
      data = china_country_outline,
      fill = NA,
      color = "black",
      linewidth = 0.42,
      lineend = "round",
      linejoin = "round"
    ) +
    fill_scale +
    coord_sf(
      crs = china_crs,
      default_crs = sf::st_crs(4326),
      xlim = main_xlim,
      ylim = main_ylim,
      expand = FALSE
    ) +
    annotation_scale(
      location = "bl",
      width_hint = 0.20,
      text_cex = 0.65,
      line_width = 0.45,
      pad_x = unit(0.25, "cm"),
      pad_y = unit(0.25, "cm")
    ) +
    annotation_north_arrow(
      location = "tl",
      which_north = "grid",
      pad_x = unit(0.35, "cm"),
      pad_y = unit(0.35, "cm"),
      style = north_arrow_nautical(
        fill = c("black", "white"),
        line_col = "black",
        text_col = "black"
      ),
      height = unit(1.0, "cm"),
      width  = unit(1.0, "cm")
    ) +
    labs(
      title = title_text,
      tag = tag_letter,
      x = "Longitude",
      y = "Latitude",
      fill = legend_title
    ) +
    theme_map_china()

  p_inset <- ggplot() +
    geom_sf(
      data = china_province_fill,
      fill = "white",
      color = "grey70",
      linewidth = 0.12,
      lineend = "round",
      linejoin = "round"
    ) +
    geom_sf(
      data = china_country_raw_proj,
      fill = NA,
      color = "black",
      linewidth = 0.28,
      lineend = "round",
      linejoin = "round"
    ) +
    coord_sf(
      crs = china_crs,
      default_crs = sf::st_crs(4326),
      xlim = south_china_sea_xlim,
      ylim = south_china_sea_ylim,
      expand = FALSE
    ) +
    theme_void(base_family = "Times New Roman") +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.40),
      plot.margin = margin(1, 1, 1, 1)
    )

  p_final <- p_main +
    patchwork::inset_element(
      p_inset,
      left = 0.865,
      bottom = 0.035,
      right = 0.998,
      top = 0.265,
      align_to = "panel"
    )

  return(p_final)
}

## =========================
## 8. 颜色设置
## =========================
canape_colors <- c(
  "Neo-endemism" = "#E41A1C",
  "Paleo-endemism" = "#377EB8",
  "Mixed-endemism" = "#984EA3",
  "Super-endemism" = "#FF7F00",
  "Not Significant" = "grey85"
)

## =========================
## 9. A 图：PE
## =========================
p_pe <- plot_china_grid(
  data = master_data,
  value_col = "PE",
  title_text = "Phylogenetic Endemism (PE)",
  legend_title = "PE",
  fill_scale = scale_fill_viridis_c(option = "plasma", na.value = "transparent"),
  tag_letter = "A"
)

## =========================
## 10. B 图：RPE
## =========================
p_rpe <- plot_china_grid(
  data = master_data,
  value_col = "RPE",
  title_text = "Relative Phylogenetic Endemism (RPE)",
  legend_title = "RPE",
  fill_scale = scale_fill_viridis_c(option = "viridis", na.value = "transparent"),
  tag_letter = "B"
)

## =========================
## 11. C 图：CANAPE
## =========================
p_canape <- plot_china_grid(
  data = master_data,
  value_col = "Endemism_Type",
  title_text = "CANAPE Categorical Endemism",
  legend_title = "Type",
  fill_scale = scale_fill_manual(
    values = canape_colors,
    drop = FALSE,
    na.value = "transparent"
  ),
  tag_letter = "C"
)

## =========================
## 12. D 图：Mean ED 年龄统计图
## =========================
age_data <- master_data %>%
  filter(!is.na(Endemism_Type), !is.na(Mean_ED_Age)) %>%
  mutate(Endemism_Type = factor(Endemism_Type, levels = endemism_levels))

kw_test <- kruskal.test(Mean_ED_Age ~ Endemism_Type, data = age_data)
kw_label <- paste0(
  "Kruskal-Wallis, p = ",
  format.pval(kw_test$p.value, digits = 3, eps = 2.2e-16)
)

valid_groups <- age_data %>%
  count(Endemism_Type) %>%
  filter(n >= 2) %>%
  pull(Endemism_Type) %>%
  as.character()

candidate_comparisons <- list(
  c("Neo-endemism", "Paleo-endemism"),
  c("Paleo-endemism", "Mixed-endemism"),
  c("Paleo-endemism", "Not Significant")
)

my_comparisons <- Filter(
  function(x) all(x %in% valid_groups),
  candidate_comparisons
)

p_age <- ggplot(
  age_data,
  aes(x = Endemism_Type, y = Mean_ED_Age, fill = Endemism_Type)
) +
  geom_boxplot(
    width = 0.45,
    outlier.shape = NA,
    alpha = 0.75,
    linewidth = 0.35
  ) +
  geom_jitter(
    aes(color = Endemism_Type),
    width = 0.12,
    alpha = 0.55,
    size = 1.45,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = canape_colors, drop = FALSE) +
  scale_color_manual(values = canape_colors, drop = FALSE) +
  labs(
    title = "Mean Evolutionary Age across Endemism Types",
    tag = "D",
    x = NULL,
    y = "Mean Evolutionary Age (ED, Ma)"
  ) +
  annotate(
    "text",
    x = 1,
    y = max(age_data$Mean_ED_Age, na.rm = TRUE) * 0.98,
    label = kw_label,
    hjust = 0,
    size = 3.2,
    family = "Times New Roman"
  ) +
  theme_age_plot()

if (length(my_comparisons) > 0) {
  p_age <- p_age +
    ggpubr::stat_compare_means(
      comparisons = my_comparisons,
      method = "wilcox.test",
      method.args = list(exact = FALSE),
      p.adjust.method = "BH",
      label = "p.signif",
      size = 3,
      family = "Times New Roman"
    )
}

## =========================
## 13. 组合最终图
## =========================
final_panel <- (p_pe | p_rpe) / (p_canape | p_age)

## =========================
## 14. 导出：PDF / SVG / PNG
## =========================
## 单图 PDF
suppressWarnings(ggsave(paste0(out_prefix, "_A_PE.pdf"), p_pe, width = 8.5, height = 6.5, dpi = 300, device = cairo_pdf))
suppressWarnings(ggsave(paste0(out_prefix, "_B_RPE.pdf"), p_rpe, width = 8.5, height = 6.5, dpi = 300, device = cairo_pdf))
suppressWarnings(ggsave(paste0(out_prefix, "_C_CANAPE.pdf"), p_canape, width = 8.5, height = 6.5, dpi = 300, device = cairo_pdf))
suppressWarnings(ggsave(paste0(out_prefix, "_D_Mean_ED_Age.pdf"), p_age, width = 7.8, height = 6.5, dpi = 300, device = cairo_pdf))

## 最终 PDF
suppressWarnings(ggsave(
  paste0(out_prefix, "_Final_Panel.pdf"),
  final_panel,
  width = 17,
  height = 12,
  dpi = 300,
  device = cairo_pdf
))

## 最终 SVG
suppressWarnings(ggsave(
  paste0(out_prefix, "_Final_Panel.svg"),
  final_panel,
  width = 17,
  height = 12,
  device = svglite::svglite
))

## 最终 PNG
ragg::agg_png(
  filename = paste0(out_prefix, "_Final_Panel_600dpi.png"),
  width = 17,
  height = 12,
  units = "in",
  res = 600,
  scaling = 1
)
print(final_panel)
dev.off()

## =========================
## 15. 输出汇总表
## =========================
summary_table <- master_data %>%
  filter(!is.na(Endemism_Type)) %>%
  group_by(Endemism_Type) %>%
  summarise(
    N_grid = n(),
    Mean_PE = mean(PE, na.rm = TRUE),
    Median_PE = median(PE, na.rm = TRUE),
    Mean_RPE = mean(RPE, na.rm = TRUE),
    Median_RPE = median(RPE, na.rm = TRUE),
    Mean_ED_Age = mean(Mean_ED_Age, na.rm = TRUE),
    Median_ED_Age = median(Mean_ED_Age, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  summary_table,
  paste0(out_prefix, "_CANAPE_Group_Summary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

message("全部完成！输出文件前缀为：", out_prefix)