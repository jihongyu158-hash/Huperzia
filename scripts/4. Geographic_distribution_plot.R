############################################################
## 中国植被类型 + 物种分布点标准制图
## 使用：国界 + 省界 + 市界 + 植被图 + 物种分布点
############################################################

library(sf)
library(dplyr)
library(ggplot2)
library(ggspatial)
library(cowplot)

############################################################
## 1. 路径设置 ####
############################################################

veg_path <- "vegetation_china.shp"

species_path <- "Huperzia_Data.csv"

## 这里改成你自己的 shp 路径
country_boundary_path  <- "国界_Project.shp"
province_boundary_path <- "省界_Project.shp"
city_boundary_path     <- "市界_Project.shp"

## 如果你有九段线或南海诸岛图层，可以填；没有就设为 NA
nine_dash_line_path <- NA
## nine_dash_line_path <- "D:/a5bc0-main/china_SHP/九段线.shp"

out_dir <- "XXX"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

############################################################
## 2. 通用函数 ####
############################################################

read_shp_gbk <- function(path) {
  
  if (is.na(path) || is.null(path) || !file.exists(path)) {
    warning("文件不存在或路径为空：", path)
    return(NULL)
  }
  
  x <- st_read(path, options = "ENCODING=GBK", quiet = TRUE)
  return(x)
}

make_valid_quiet <- function(x) {
  
  if (is.null(x)) return(NULL)
  
  x <- suppressWarnings(st_make_valid(x))
  x <- x[!st_is_empty(x), ]
  
  return(x)
}

transform_to_crs <- function(x, crs_target) {
  
  if (is.null(x)) return(NULL)
  
  if (is.na(st_crs(x))) {
    stop("错误：某个 shp 文件没有 CRS，请先在 GIS 软件中确认坐标系统。")
  }
  
  st_transform(x, crs_target)
}

############################################################
## 3. 设置投影 ####
############################################################

crs_wgs84 <- st_crs(4326)

## 中国常用 Albers 等积投影
crs_china_albers <- st_crs(
  "+proj=aea +lat_1=25 +lat_2=47 +lat_0=0 +lon_0=105 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
)

############################################################
## 4. 读取植被图层 ####
############################################################

veg <- read_shp_gbk(veg_path)
veg <- make_valid_quiet(veg)

cat("\n================ 植被图层字段 ================\n")
print(names(veg))

code_field <- "植被型编号"

if (!code_field %in% names(veg)) {
  stop("错误：植被图层中没有字段：", code_field, "。请先运行 names(veg) 查看真实字段名。")
}

cat("\n================ 植被编码 ================\n")
print(sort(unique(as.character(veg[[code_field]]))))

############################################################
## 5. 植被类型编码映射 ####
############################################################

code_map <- c(
  `6`  = "温带针叶、落叶阔叶混交林",
  `15` = "热带雨林",
  `12` = "亚热带常绿阔叶林",
  `5`  = "亚热带和热带山地针叶林"
 )

veg <- veg %>%
  mutate(
    veg_code = as.character(.data[[code_field]]),
    veg_name = dplyr::recode(
      veg_code,
      !!!code_map,
      .default = NA_character_
    )
  )

target_types <- c(
  "温带针叶、落叶阔叶混交林",
  "亚热带和热带山地针叶林",
  "热带雨林",
  "亚热带常绿阔叶林"
)

veg_sel <- veg %>%
  filter(veg_name %in% target_types)

if (nrow(veg_sel) == 0) {
  stop("错误：veg_sel 为空，请检查 code_map 中的编码是否与植被 shp 中的编码一致。")
}

cat("\n================ 筛选后的植被类型 ================\n")
print(table(veg_sel$veg_name))

############################################################
## 6. 读取国界、省界、市界 ####
############################################################

china_country <- read_shp_gbk(country_boundary_path)
china_province <- read_shp_gbk(province_boundary_path)
china_city <- read_shp_gbk(city_boundary_path)
nine_dash_line <- read_shp_gbk(nine_dash_line_path)

china_country <- make_valid_quiet(china_country)
china_province <- make_valid_quiet(china_province)
china_city <- make_valid_quiet(china_city)
nine_dash_line <- make_valid_quiet(nine_dash_line)

if (is.null(china_country)) {
  stop("错误：没有成功读取国界 shp。请检查 country_boundary_path。")
}

if (is.null(china_province)) {
  warning("没有读取到省界 shp，后续将不绘制省界。")
}

if (is.null(china_city)) {
  warning("没有读取到市界 shp，后续将不绘制市界。")
}

############################################################
## 7. 统一投影 ####
############################################################

veg_sel <- transform_to_crs(veg_sel, crs_china_albers)
china_country <- transform_to_crs(china_country, crs_china_albers)
china_province <- transform_to_crs(china_province, crs_china_albers)
china_city <- transform_to_crs(china_city, crs_china_albers)
nine_dash_line <- transform_to_crs(nine_dash_line, crs_china_albers)

############################################################
## 8. 读取物种分布点 ####
############################################################

species_data <- read.csv(
  species_path,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("\n================ 物种数据字段 ================\n")
print(names(species_data))
print(head(species_data))

if (!all(c("Lon", "Lat") %in% names(species_data))) {
  stop("错误：物种分布数据中没有 Lon 和 Lat 两列。")
}

if (!"species" %in% names(species_data)) {
  stop("错误：物种分布数据中没有 species 列。")
}

species_data <- species_data %>%
  mutate(
    Lon = as.numeric(Lon),
    Lat = as.numeric(Lat),
    species = as.character(species)
  ) %>%
  filter(
    !is.na(Lon),
    !is.na(Lat),
    Lon >= 70,
    Lon <= 140,
    Lat >= 0,
    Lat <= 60
  )

## 物种分组
## 注意：这里的名字必须和 species_data$species 完全一致
species_to_group <- c(
  "HuperziaA" = "Group1",
  "HuperziaB" = "Group2",
  "HuperziaC" = "Group3"
)

species_data <- species_data %>%
  mutate(
    group = dplyr::recode(
      species,
      !!!species_to_group,
      .default = "未分组"
    )
  )

cat("\n================ 物种分组统计 ================\n")
print(table(species_data$group, useNA = "ifany"))

species_sf <- st_as_sf(
  species_data,
  coords = c("Lon", "Lat"),
  crs = crs_wgs84,
  remove = FALSE
)

species_sf <- st_transform(species_sf, crs_china_albers)

############################################################
## 9. 点位叠加植被类型 ####
############################################################

species_vegetation <- st_join(
  species_sf,
  veg_sel[, c("veg_name")],
  join = st_intersects,
  left = TRUE
)

write.csv(
  st_drop_geometry(species_vegetation),
  file = file.path(out_dir, "物种点位_叠加植被类型结果.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

cat("\n================ 点位对应植被类型统计 ================\n")
print(table(species_vegetation$veg_name, useNA = "ifany"))

############################################################
## 10. 配色 ####
############################################################

veg_colors <- c(
  "温带针叶、落叶阔叶混交林" = "#61D04F",
  "热带雨林" = "#923262",
  "亚热带常绿阔叶林" = "#DF536B",
  "亚热带和热带山地针叶林" = "#2297E6",
)

veg_colors <- veg_colors[names(veg_colors) %in% unique(veg_sel$veg_name)]

group_colors <- c(
  "Group1" = "#155C00",
  "Group2" = "#0300A4",
  "Group3" = "#FFA810",
  "未分组" = "black"
)

############################################################
## 11. 设置主图和南海小图显示范围：修正版 ####
############################################################

## 关键修改：
## 不再把经纬度 bbox 手动转换成 Albers 后再裁剪。
## coord_sf 里使用 default_crs = crs_wgs84，让 xlim / ylim 按经纬度解释。
## 这样海南岛、黑龙江北部都更不容易被裁掉。

main_xlim <- c(72, 137)
main_ylim <- c(15, 56)

## 海南岛约在 108–111°E, 18–20°N
## 所以 main_ylim 最低不能高于 18，建议保留到 15 或更低。
## 如果你想连南海小图也更自然，可以设为 c(3, 56)，但主图会变得很长。
south_xlim <- c(105, 125)
south_ylim <- c(3, 25)


############################################################
## 12. 主图：完整显示大陆、海南岛和东北 ####
############################################################

main_map <- ggplot() +
  
  ## 目标植被类型
  geom_sf(
    data = veg_sel,
    aes(fill = veg_name),
    color = NA,
    alpha = 0.95
  ) +
  
  ## 市界，颜色浅一些
  {
    if (!is.null(china_city)) {
      geom_sf(
        data = china_city,
        fill = NA,
        color = "grey82",
        linewidth = 0.12
      )
    }
  } +
  
  ## 省界，稍深一些
  {
    if (!is.null(china_province)) {
      geom_sf(
        data = china_province,
        fill = NA,
        color = "grey45",
        linewidth = 0.35
      )
    }
  } +
  
  ## 国界，黑色加粗
  geom_sf(
    data = china_country,
    fill = NA,
    color = "black",
    linewidth = 0.75
  ) +
  
  ## 九段线，如果有
  {
    if (!is.null(nine_dash_line)) {
      geom_sf(
        data = nine_dash_line,
        color = "black",
        linewidth = 0.45,
        linetype = "longdash"
      )
    }
  } +
  
  ## 物种分布点
  geom_sf(
    data = species_sf,
    aes(color = group),
    size = 1.4,
    shape = 17,
    alpha = 0.9
  ) +
  
  scale_fill_manual(
    values = veg_colors,
    name = "植被类型",
    drop = FALSE
  ) +
  
  scale_color_manual(
    values = group_colors,
    name = "物种类别",
    drop = FALSE
  ) +
  
  annotation_scale(
    location = "bl",
    width_hint = 0.25,
    text_cex = 0.8,
    line_width = 0.5
  ) +
  
  annotation_north_arrow(
    location = "tl",
    which_north = "true",
    pad_x = unit(0.35, "cm"),
    pad_y = unit(0.35, "cm"),
    style = north_arrow_fancy_orienteering
  ) +
  
  ## 这里是关键：xlim / ylim 用经纬度解释
  coord_sf(
    crs = crs_china_albers,
    default_crs = crs_wgs84,
    xlim = main_xlim,
    ylim = main_ylim,
    expand = FALSE,
    clip = "on"
  ) +
  
  labs(
    title = "中国主要植被类型与石杉类物种分布",
    subtitle = "Vegetation types and species occurrence records in China",
    x = NULL,
    y = NULL,
    caption = "注：正式发表或公开展示时，请使用标准地图底图并按要求标注审图号。"
  ) +
  
  theme_bw(base_size = 13) +
  theme(
    panel.grid.major = element_line(color = "grey88", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    axis.text = element_text(size = 8),
    axis.ticks = element_line(linewidth = 0.25),
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 9),
    legend.key.height = unit(0.45, "cm"),
    legend.key.width = unit(0.55, "cm"),
    plot.caption = element_text(size = 8, hjust = 0)
  )


############################################################
## 13. 南海诸岛小图：修正版 ####
############################################################

south_china_sea_map <- ggplot() +
  
  geom_sf(
    data = veg_sel,
    aes(fill = veg_name),
    color = NA,
    alpha = 0.95,
    show.legend = FALSE
  ) +
  
  {
    if (!is.null(china_city)) {
      geom_sf(
        data = china_city,
        fill = NA,
        color = "grey82",
        linewidth = 0.1,
        show.legend = FALSE
      )
    }
  } +
  
  {
    if (!is.null(china_province)) {
      geom_sf(
        data = china_province,
        fill = NA,
        color = "grey55",
        linewidth = 0.2,
        show.legend = FALSE
      )
    }
  } +
  
  geom_sf(
    data = china_country,
    fill = NA,
    color = "black",
    linewidth = 0.5,
    show.legend = FALSE
  ) +
  
  {
    if (!is.null(nine_dash_line)) {
      geom_sf(
        data = nine_dash_line,
        color = "black",
        linewidth = 0.4,
        linetype = "longdash",
        show.legend = FALSE
      )
    }
  } +
  
  scale_fill_manual(values = veg_colors) +
  
  coord_sf(
    crs = crs_china_albers,
    default_crs = crs_wgs84,
    xlim = south_xlim,
    ylim = south_ylim,
    expand = FALSE
  ) +
  
  labs(title = "南海诸岛") +
  
  theme_bw(base_size = 8) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 8),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(1, 1, 1, 1, "mm")
  )


############################################################
## 14. 合成主图和南海小图：避免遮挡海南岛 ####
############################################################

## 之前小图位置可能偏左、偏大，会遮挡华南—海南附近区域。
## 这里把南海小图缩小，并放到更靠右下的位置。
final_map <- ggdraw() +
  draw_plot(main_map, x = 0, y = 0, width = 1, height = 1) +
  draw_plot(
    south_china_sea_map,
    x = 0.73,
    y = 0.06,
    width = 0.18,
    height = 0.25
  )

############################################################
## 15. 导出 ####
############################################################

ggsave(
  filename = file.path(out_dir, "中国植被与物种分布图_完整国界版.pdf"),
  plot = final_map,
  device = cairo_pdf,
  width = 12,
  height = 8,
  dpi = 300,
  limitsize = FALSE
)

ggsave(
  filename = file.path(out_dir, "中国植被与物种分布图_完整国界版.png"),
  plot = final_map,
  width = 12,
  height = 8,
  dpi = 600,
  limitsize = FALSE
)

cat("\n============================================================\n")
cat("绘图完成！\n")
cat("输出文件：\n")
cat(file.path(out_dir, "中国植被与物种分布图_完整国界版.pdf"), "\n")
cat(file.path(out_dir, "中国植被与物种分布图_完整国界版.png"), "\n")
cat("============================================================\n")