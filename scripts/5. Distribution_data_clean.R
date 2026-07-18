# ==============================================================
# 步骤 0：安装并加载必备包
# ==============================================================
library(ape)
library(CoordinateCleaner)
library(dplyr)
library(tidyr)

# ==============================================================
# 步骤 1：读取系统树与分布数据
# ==============================================================
phy <- read.tree("23species_cds_beast_time_tree_onlyHuperzia.tree")
df_raw <- read.csv("Huperzia_Data.csv")


# ==============================================================
# 步骤 2：数据匹配与空间大清洗 (Data Cleaning)
# ==============================================================
# 2.1 强制过滤：只保留在系统树(phy)上存在的物种分布点
df_matched <- df_raw %>% 
  filter(species %in% phy$tip.label)

cat("与系统树匹配后，剩余记录数:", nrow(df_matched), "\n")

# 2.2 框定中国及周边区域 (经度 73~136，纬度 18~54)
df_china <- df_matched %>%
  filter(Lon >= 73 & Lon <= 136) %>%
  filter(Lat >= 18 & Lat <= 54) %>%
  filter(!is.na(Lon) & !is.na(Lat))

#2.3 CoordinateCleaner 严苛排雷 (去除首都、海洋、零点等异常坐标)
flags <- clean_coordinates(x = df_china, 
                           lon = "Lon", lat = "Lat", 
                           tests = c("capitals", "centroids", "equal", "gbif", "institutions", "seas", "zeros"))
df_clean <- df_china[flags$.summary, ]

cat("经过坐标合规性清洗后，剩余有效分布点数:", nrow(df_clean), "\n")

# ==============================================================
# 步骤 3：栅格化与生成 CANAPE 群落矩阵 (Community Matrix)
# ==============================================================
grid_size <- 0.2 # 设定 0.2度 × 0.2度 网格

df_grid <- df_clean %>%
  mutate(
    grid_Lon = round(Lon / grid_size) * grid_size,
    grid_Lat = round(Lat / grid_size) * grid_size,
    Grid_ID = paste("Grid", grid_Lon, grid_Lat, sep = "_")
  )

# 空间去重：同一个网格内，同一个物种只保留1条记录
df_unique <- df_grid %>%
  distinct(Grid_ID, species, .keep_all = TRUE)

# 交叉表生成群落矩阵：行是Grid，列是物种
comm_matrix <- table(df_unique$Grid_ID, df_unique$species)
comm_matrix <- as.data.frame.matrix(comm_matrix)

# 转为二元矩阵：大于0的均设为1
comm_matrix[comm_matrix > 0] <- 1

# ==============================================================
# 步骤 4：CANAPE 终极安全检查 (极其重要！)
# ==============================================================
missing_in_tree <- setdiff(colnames(comm_matrix), phy$tip.label)
missing_in_matrix <- setdiff(phy$tip.label, colnames(comm_matrix))

if(length(missing_in_tree) > 0) {
  cat("\n⚠️ 警告：以下物种在矩阵中，但不在树上（请检查拼写）：\n", missing_in_tree)
} else {
  cat("\n✅ 完美：矩阵中所有物种都在系统树中！")
}

if(length(missing_in_matrix) > 0) {
  cat("\n⚠️ 致命警告：以下系统树上的物种，在清洗后的分布数据中没有留下任何点！\n", missing_in_matrix, "\n")
  cat("👉 解决办法：CANAPE 不允许树上有的物种但在地图上毫无分布。你必须在树中把这些没分布点的 Tip 剪掉 (使用 drop.tip 函数)。")
} else {
  cat("\n✅ 完美：系统树上的所有 24 个物种，都在地图上拥有至少一个有效分布网格！可以安全运行 CANAPE！\n")
}

# 导出矩阵
write.csv(comm_matrix, "Huperzia_Community_Matrix0.2.csv", row.names = TRUE)

