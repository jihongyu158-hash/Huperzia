# ============================================================================
# 空间系统发育 / CANAPE 计算脚本
# 输入：Huperzia_Community_Matrix0.2.csv
# 输出：Huperzia_Spatial_Phylo_Master_Results.csv / .rds / 统计汇总表
# ============================================================================

# ------------------------------
# 0. 环境准备
# ------------------------------
required_pkgs <- c("canaper", "ape", "dplyr", "tidyr", "tibble", "picante")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("请先安装缺失的 R 包：", paste(missing_pkgs, collapse = ", "), call. = FALSE)
}

library(canaper)
library(ape)
library(dplyr)
library(tidyr)
library(tibble)
library(picante)

set.seed(123)
getwd()
# ------------------------------
# 1. 参数区：按需修改
# ------------------------------
comm_file <- "Huperzia_Community_Matrix0.2.csv"
n_reps <- 999
null_model <- "curveball"

# 24 物种系统树
# 注意：树的 tip.label 必须与群落矩阵列名一致
# 如果你的列名中有空格、横线等，请先统一成与树 tip.label 完全一致的格式
# ------------------------------
tree_text <- "((((((((Huperzia_liangshanica:1.3457692936435421,Huperzia_crassifolia:1.3457692936435421):8.434076113852882,(Huperzia_emeiensis:0.2619059635104577,Huperzia_nanchuanensis:0.2619059635104577):9.517939443985966):5.807332679813925,((Huperzia_nanlingensis:2.5345343240414593,Huperzia_sutchueniana:2.5345343240414593):6.177774723224815,Huperzia_javanica:8.712309047266274):6.874869040044075):2.2482391459068936,Huperzia_crispata:17.835417233217242):3.467424375982432,(Huperzia_herteriana:2.1712315891857656,Huperzia_quasipolytrichoides:2.1712315891857656):19.13161002001391):8.501452008721174,((((Huperzia_kunmingensis:1.791129648400812,Huperzia_somae:1.791129648400812):1.4106207301329619,Huperzia_rubicaulis:3.2017503785337738):14.0776059876504,Huperzia_delavayi:17.279356366184174):6.678505484553284,(Huperzia_chinensis:1.6072114905779529,Huperzia_tibetica:1.6072114905779529):22.350650360159506):5.84643176718339):5.116831982274107,(Huperzia_lucidula:9.876098801097669,Huperzia_serrata:9.876098801097669):25.045026799097286):38.548436995926636,((((Huperzia_appressa:0.048575143630387174,Huperzia_arctica:0.048575143630387174):2.0525355073167475,Huperzia_asiatica:2.1011106509471347):3.7804668649243354,Huperzia_selago:5.88157751587147):6.055809756218764,Huperzia_miyoshiana:11.937387272090234):61.53217532403136);"
phy <- read.tree(text = tree_text)

# ------------------------------
# 2. 读取并检查群落矩阵
# ------------------------------
if (!file.exists(comm_file)) {
  stop("找不到群落矩阵文件：", comm_file, call. = FALSE)
}

comm_matrix <- read.csv(comm_file, row.names = 1, check.names = FALSE)
comm_matrix <- as.matrix(comm_matrix)
storage.mode(comm_matrix) <- "numeric"

# 转为 0/1 矩阵，避免 abundance 被误当作 presence/absence
comm_matrix[is.na(comm_matrix)] <- 0
comm_matrix <- ifelse(comm_matrix > 0, 1, 0)

# 删除空网格和空物种
comm_matrix <- comm_matrix[rowSums(comm_matrix) > 0, colSums(comm_matrix) > 0, drop = FALSE]

# 检查物种名是否匹配
sp_in_comm <- colnames(comm_matrix)
sp_in_tree <- phy$tip.label
missing_in_tree <- setdiff(sp_in_comm, sp_in_tree)
missing_in_comm <- setdiff(sp_in_tree, sp_in_comm)

if (length(missing_in_tree) > 0) {
  stop(
    "群落矩阵中有物种不在系统树中：\n",
    paste(missing_in_tree, collapse = ", "),
    "\n请统一群落矩阵列名和树 tip.label。",
    call. = FALSE
  )
}

if (length(missing_in_comm) > 0) {
  message("系统树中有物种不在群落矩阵中，将从树上剪掉：", paste(missing_in_comm, collapse = ", "))
  phy <- drop.tip(phy, missing_in_comm)
}

# 让群落矩阵列顺序与系统树 tip.label 一致
comm_matrix <- comm_matrix[, phy$tip.label, drop = FALSE]

cat("成功加载数据：共有 ", nrow(comm_matrix), " 个有效网格，", ncol(comm_matrix), " 个物种。\n", sep = "")

# ------------------------------
# 3. 独立指标：物种丰富度与平均 ED age
# ------------------------------
richness_df <- data.frame(
  Grid_ID = rownames(comm_matrix),
  Richness = rowSums(comm_matrix > 0),
  stringsAsFactors = FALSE
)

# picante::evol.distinct 不同版本返回列名可能略有差别，这里做兼容处理
ed_values <- evol.distinct(phy, type = "fair.proportion")
ed_values <- as.data.frame(ed_values)

if (!"Species" %in% names(ed_values)) {
  ed_values$Species <- rownames(ed_values)
}

ed_col <- dplyr::case_when(
  "w" %in% names(ed_values) ~ "w",
  "ED" %in% names(ed_values) ~ "ED",
  "ed" %in% names(ed_values) ~ "ed",
  TRUE ~ NA_character_
)

if (is.na(ed_col)) {
  stop("无法在 evol.distinct() 结果中识别 ED 列。请检查 ed_values 的列名：", paste(names(ed_values), collapse = ", "), call. = FALSE)
}

ed_tbl <- ed_values %>%
  transmute(Species = .data$Species, ED_Age = as.numeric(.data[[ed_col]]))

age_df <- data.frame(
  Grid_ID = rownames(comm_matrix),
  Mean_ED_Age = apply(comm_matrix, 1, function(x) {
    sp_in_grid <- names(x)[x > 0]
    mean(ed_tbl$ED_Age[ed_tbl$Species %in% sp_in_grid], na.rm = TRUE)
  }),
  stringsAsFactors = FALSE
)

# ------------------------------
# 4. CANAPE 核心计算
# ------------------------------
cat("正在运行 CANAPE 零模型随机化检验：null_model = ", null_model,
    ", n_reps = ", n_reps, "。\n", sep = "")

rand_res <- cpr_rand_test(
  comm = comm_matrix,
  phy = phy,
  null_model = null_model,
  n_reps = n_reps
)

endem_res <- cpr_classify_endem(rand_res)

# ------------------------------
# 5. 合并 master_data
# ------------------------------
master_data <- endem_res %>%
  as.data.frame() %>%
  rownames_to_column("Grid_ID") %>%
  mutate(
    PE = as.numeric(pe_obs),
    RPE = if_else(!is.na(pe_alt_obs) & pe_alt_obs > 0, pe_obs / pe_alt_obs, NA_real_),
    Endemism_Type = factor(
      endem_type,
      levels = c("neo", "paleo", "mixed", "super", "not significant"),
      labels = c(
        "Neo-endemism",
        "Paleo-endemism",
        "Mixed-endemism",
        "Super-endemism",
        "Not Significant"
      )
    )
  ) %>%
  left_join(richness_df, by = "Grid_ID") %>%
  left_join(age_df, by = "Grid_ID") %>%
  tidyr::extract(
    Grid_ID,
    into = c("Lon", "Lat"),
    regex = "^.*?(-?\\d+(?:\\.\\d+)?)_(-?\\d+(?:\\.\\d+)?)$",
    remove = FALSE,
    convert = TRUE
  ) %>%
  mutate(
    Lon = as.numeric(Lon),
    Lat = as.numeric(Lat)
  )

# 如果 Grid_ID 不是 Grid_经度_纬度 这种格式，这里会提示你修改解析规则
if (any(is.na(master_data$Lon) | is.na(master_data$Lat))) {
  warning(
    "部分 Grid_ID 未能解析出 Lon/Lat。请检查 Grid_ID 是否类似 Grid_120.1_30.1 或 120.1_30.1。",
    call. = FALSE
  )
}

# ------------------------------
# 6. 导出结果和统计表
# ------------------------------
write.csv(master_data, "Huperzia_Spatial_Phylo_Master_Results.csv", row.names = FALSE)
saveRDS(master_data, "Huperzia_Spatial_Phylo_Master_Results.rds")

cat("\n核心计算完成：已导出 Huperzia_Spatial_Phylo_Master_Results.csv 和 .rds。\n")

cat("\n=======================================================\n")
cat("物种丰富度频数：\n")
print(table(master_data$Richness, useNA = "ifany"))

cat("\n=======================================================\n")
cat("各 CANAPE 类型的丰富度与演化年龄摘要：\n")
final_stat <- master_data %>%
  group_by(Endemism_Type) %>%
  summarise(
    Grid_Count = n(),
    Mean_Richness = round(mean(Richness, na.rm = TRUE), 2),
    Median_Richness = round(median(Richness, na.rm = TRUE), 2),
    Grids_with_1_Sp = sum(Richness == 1, na.rm = TRUE),
    Mean_PE = round(mean(PE, na.rm = TRUE), 4),
    Mean_RPE = round(mean(RPE, na.rm = TRUE), 4),
    Mean_Evolutionary_Age = round(mean(Mean_ED_Age, na.rm = TRUE), 2),
    Median_Evolutionary_Age = round(median(Mean_ED_Age, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  arrange(desc(Grid_Count))

print(final_stat)
write.csv(final_stat, "Huperzia_CANAPE_Group_Summary.csv", row.names = FALSE)
cat("=======================================================\n")

# ------------------------------
# 7. 统计检验
# ------------------------------
stats_data <- master_data %>%
  filter(!is.na(Endemism_Type), !is.na(Mean_ED_Age)) %>%
  droplevels()

if (nlevels(stats_data$Endemism_Type) >= 2) {
  cat("\nKruskal-Wallis test for Mean_ED_Age ~ Endemism_Type:\n")
  print(kruskal.test(Mean_ED_Age ~ Endemism_Type, data = stats_data))

  cat("\nPairwise Wilcoxon test, BH adjusted; exact = FALSE 避免 ties 警告：\n")
  pairwise_age <- pairwise.wilcox.test(
    x = stats_data$Mean_ED_Age,
    g = stats_data$Endemism_Type,
    p.adjust.method = "BH",
    exact = FALSE
  )
  print(pairwise_age)
} else {
  message("有效 Endemism_Type 少于 2 组，跳过组间统计检验。")
}

