library(ape)
library(phytools)
library(geiger)

phy <- read.tree("23species_cds_beast_time_tree_onlyHuperzia.tree")

## 检查树是否有枝长
if (is.null(phy$edge.length)) {
  stop("错误：系统发育树没有 branch lengths，无法进行后续模型拟合和祖先状态重建。")
}

## 检查是否有重复 tip
if (any(duplicated(phy$tip.label))) {
  dup_tips <- phy$tip.label[duplicated(phy$tip.label)]
  print(dup_tips)
  stop("错误：树中存在重复 tip.label，请先处理。")
}

## 读取生态位数据 ####
read_trait_table <- function(file) {
  
  dat1 <- try(
    read.csv(
      file,
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    silent = TRUE
  )
  
  if (!inherits(dat1, "try-error") && ncol(dat1) >= 2) {
    return(dat1)
  }
  
  dat2 <- try(
    read.delim(
      file,
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    silent = TRUE
  )
  
  if (!inherits(dat2, "try-error") && ncol(dat2) >= 2) {
    return(dat2)
  }
  
  dat3 <- read.table(
    file,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  return(dat3)
}

traits <- read_trait_table("biome_data1.csv")

## 去掉列名首尾空格
colnames(traits) <- trimws(colnames(traits))

cat("\n================ 读取到的生态位数据 ================\n")
cat("数据维度：\n")
print(dim(traits))

cat("列名：\n")
print(colnames(traits))

cat("前几行：\n")
print(head(traits))


##检查并提取 species 和 biome 两列 ####

## 自动识别 species 和 biome 列，大小写不敏感
species_col <- grep("^species$", colnames(traits), ignore.case = TRUE, value = TRUE)
biome_col   <- grep("^biome$",   colnames(traits), ignore.case = TRUE, value = TRUE)

if (length(species_col) != 1) {
  stop("错误：没有找到唯一的 species 列。请确认表头是否为 species。")
}

if (length(biome_col) != 1) {
  stop("错误：没有找到唯一的 biome 列。请确认表头是否为 biome。")
}

## 只保留 species 和 biome 两列
traits <- traits[, c(species_col, biome_col)]
colnames(traits) <- c("species", "biome")

## 去掉首尾空格
traits$species <- trimws(as.character(traits$species))
traits$biome   <- trimws(as.character(traits$biome))

## 删除缺失值和空值
traits <- traits[
  !is.na(traits$species) & traits$species != "" &
    !is.na(traits$biome) & traits$biome != "",
]

## 检查是否有重复物种
if (any(duplicated(traits$species))) {
  dup_species <- traits$species[duplicated(traits$species)]
  cat("\n以下物种在 biome_data1.csv 中重复：\n")
  print(dup_species)
  stop("错误：性状表中存在重复 species，请先检查并去重。")
}

cat("\n================ 清理后的生态位数据 ================\n")
cat("物种数：", nrow(traits), "\n")
cat("生态位状态统计：\n")
print(table(traits$biome))

## 匹配系统发育树和生态位数据 ####

tree_species  <- phy$tip.label
trait_species <- traits$species

common_species <- intersect(tree_species, trait_species)

cat("\n================ 物种匹配情况 ================\n")
cat("树中的物种数：", length(tree_species), "\n")
cat("性状表中的物种数：", length(trait_species), "\n")
cat("成功匹配的物种数：", length(common_species), "\n")

tree_only <- setdiff(tree_species, trait_species)
trait_only <- setdiff(trait_species, tree_species)

cat("只在树中、不在性状表中的物种数：", length(tree_only), "\n")
cat("只在性状表中、不在树中的物种数：", length(trait_only), "\n")

## 输出未匹配物种，方便排查
write.csv(
  data.frame(tree_only = tree_only),
  file = "Species_only_in_tree_not_in_biome_data.csv",
  row.names = FALSE
)

write.csv(
  data.frame(trait_only = trait_only),
  file = "Species_only_in_biome_data_not_in_tree.csv",
  row.names = FALSE
)

if (length(common_species) < 2) {
  stop("错误：树和性状表匹配到的物种少于 2 个，无法分析。请检查物种名是否一致。")
}

## 剪切系统发育树，只保留有生态位数据的物种
phy <- drop.tip(phy, setdiff(phy$tip.label, common_species))

## 按照剪切后的树 tip.label 顺序排列性状
traits_matched <- traits[match(phy$tip.label, traits$species), ]

## 再次检查顺序
if (!all(traits_matched$species == phy$tip.label)) {
  stop("错误：性状数据顺序与树 tip.label 顺序不一致。")
}

## 生成命名向量 Hab
Hab <- traits_matched$biome
names(Hab) <- traits_matched$species

## 如果 biome 都是数字，例如 1,2,3,4,5，则按数字顺序设定 factor level
## 如果 biome 是文字，则按字母顺序设定 level
if (all(grepl("^[0-9]+$", Hab))) {
  Hab_levels <- as.character(sort(as.numeric(unique(Hab))))
} else {
  Hab_levels <- sort(unique(Hab))
}

Hab <- factor(Hab, levels = Hab_levels)
Hab <- droplevels(Hab)

## 检查最终数据
if (length(Hab) != Ntip(phy)) {
  stop("错误：Hab 长度与树的 tip 数量不一致。")
}

if (any(is.na(Hab))) {
  print(Hab[is.na(Hab)])
  stop("错误：Hab 中存在 NA。")
}

if (nlevels(Hab) < 2) {
  stop("错误：生态位状态少于 2 类，无法进行离散性状祖先状态重建。")
}

cat("\n================ 最终用于分析的数据 ================\n")
cat("最终树的物种数：", Ntip(phy), "\n")
cat("最终生态位状态数量：", nlevels(Hab), "\n")
cat("各生态位状态的物种数量：\n")
print(table(Hab))

##设置颜色 ####
base_cols <- c(
  "#ff3b30", "#ff9500", "#ffcc00",
  "#4cd964", "#5ac8fa", "#5856d6",
  "#007aff", "#af52de", "#8e8e93", "#34c759"
)

if (nlevels(Hab) <= length(base_cols)) {
  colors <- setNames(base_cols[1:nlevels(Hab)], levels(Hab))
} else {
  colors <- setNames(
    colorRampPalette(base_cols)(nlevels(Hab)),
    levels(Hab)
  )
}

cat("\n================ 生态位状态颜色 ================\n")
print(colors)

##绘制 tip 生态位状态 ####

pdf("Habitat_tip_states.pdf", width = 8, height = 10)

par(mar = c(5.5, 5.5, 5.5, 5.5))
plotTree(phy, fsize = 0.4)

tiplabels(
  pie = to.matrix(Hab, levels(Hab)),
  piecol = colors,
  cex = 0.5
)

usr <- par("usr")
x_safe <- usr[1] + 0.05 * (usr[2] - usr[1])
y_safe <- usr[3] + 0.01 * (usr[4] - usr[3])

add.simmap.legend(
  colors = colors,
  prompt = FALSE,
  x = x_safe,
  y = y_safe,
  fsize = 0.8,
  vertical = FALSE
)

dev.off()



##拟合 ER / SYM / ARD 模型 ####

models_to_test <- c("ER", "SYM", "ARD")

safe_fitDiscrete <- function(model_name) {
  
  cat("\n正在拟合模型：", model_name, "\n")
  
  fit <- tryCatch(
    {
      fitDiscrete(
        phy = phy,
        dat = Hab,
        model = model_name
      )
    },
    error = function(e) {
      message("模型 ", model_name, " 拟合失败：", e$message)
      return(NULL)
    }
  )
  
  return(fit)
}

fits <- setNames(
  lapply(models_to_test, safe_fitDiscrete),
  models_to_test
)


## 提取 log-likelihood、AIC、AICc、free parameters ####

get_free_transition_rates <- function(model_name, n_states) {
  
  if (model_name == "ER") {
    return(1)
  }
  
  if (model_name == "SYM") {
    return(n_states * (n_states - 1) / 2)
  }
  
  if (model_name == "ARD") {
    return(n_states * (n_states - 1))
  }
  
  return(NA)
}

extract_fit_stats <- function(fit, model_name, n_states, n_tips) {
  
  if (is.null(fit)) {
    return(data.frame(
      Model = model_name,
      log_likelihood = NA,
      AIC = NA,
      AICc = NA,
      free_parameters = NA,
      free_transition_rates = get_free_transition_rates(model_name, n_states),
      stringsAsFactors = FALSE
    ))
  }
  
  opt <- fit$opt
  
  logLik_value <- if ("lnL" %in% names(opt)) {
    opt$lnL
  } else if ("logLik" %in% names(opt)) {
    opt$logLik
  } else {
    NA
  }
  
  ## geiger::fitDiscrete 中的 k 是模型实际用于 AIC/AICc 的参数数量
  k_value <- if ("k" %in% names(opt)) {
    opt$k
  } else if ("npars" %in% names(opt)) {
    opt$npars
  } else {
    NA
  }
  
  AIC_value <- if ("aic" %in% names(opt)) {
    opt$aic
  } else if (!is.na(logLik_value) && !is.na(k_value)) {
    -2 * logLik_value + 2 * k_value
  } else {
    NA
  }
  
  AICc_value <- if ("aicc" %in% names(opt)) {
    opt$aicc
  } else if (!is.na(AIC_value) && !is.na(k_value)) {
    
    if ((n_tips - k_value - 1) > 0) {
      AIC_value + (2 * k_value * (k_value + 1)) / (n_tips - k_value - 1)
    } else {
      NA
    }
    
  } else {
    NA
  }
  
  data.frame(
    Model = model_name,
    log_likelihood = logLik_value,
    AIC = AIC_value,
    AICc = AICc_value,
    free_parameters = k_value,
    free_transition_rates = get_free_transition_rates(model_name, n_states),
    stringsAsFactors = FALSE
  )
}

model_stats <- do.call(
  rbind,
  lapply(
    models_to_test,
    function(m) {
      extract_fit_stats(
        fit = fits[[m]],
        model_name = m,
        n_states = nlevels(Hab),
        n_tips = Ntip(phy)
      )
    }
  )
)


##计算 delta AIC / AIC weight / delta AICc / AICc weight ####

add_delta_weight <- function(tab, criterion) {
  
  values <- tab[[criterion]]
  valid <- is.finite(values)
  
  delta_name <- paste0("delta_", criterion)
  weight_name <- paste0(criterion, "_weight")
  
  tab[[delta_name]] <- NA
  tab[[weight_name]] <- NA
  
  if (any(valid)) {
    
    delta <- rep(NA, length(values))
    delta[valid] <- values[valid] - min(values[valid], na.rm = TRUE)
    
    weight <- rep(NA, length(values))
    weight[valid] <- exp(-0.5 * delta[valid]) /
      sum(exp(-0.5 * delta[valid]), na.rm = TRUE)
    
    tab[[delta_name]] <- delta
    tab[[weight_name]] <- weight
  }
  
  return(tab)
}

model_stats <- add_delta_weight(model_stats, "AIC")
model_stats <- add_delta_weight(model_stats, "AICc")

model_stats_AICc <- model_stats[order(model_stats$AICc), ]
model_stats_AIC  <- model_stats[order(model_stats$AIC), ]

cat("\n================ 模型比较结果：按 AICc 排序 ================\n")
print(model_stats_AICc)

cat("\n================ 模型比较结果：按 AIC 排序 ================\n")
print(model_stats_AIC)


## 输出模型比较结果
write.csv(
  model_stats_AICc,
  file = "Habitat_model_comparison_logLik_AIC_AICc_free_parameters.csv",
  row.names = FALSE
)

capture.output(
  {
    cat("============================================================\n")
    cat("Discrete character model comparison for habitat / biome\n")
    cat("============================================================\n\n")
    
    cat("Tree file:\n")
    cat("23cdsjiujiu290new.tree\n\n")
    
    cat("Trait file:\n")
    cat("biome_data1.csv\n\n")
    
    cat("Number of tips used in analysis:\n")
    print(Ntip(phy))
    
    cat("\nTrait states:\n")
    print(table(Hab))
    
    cat("\n\n================ ER model ================\n")
    print(fits$ER)
    
    cat("\n\n================ SYM model ================\n")
    print(fits$SYM)
    
    cat("\n\n================ ARD model ================\n")
    print(fits$ARD)
    
    cat("\n\n================ Model comparison table, sorted by AICc ================\n")
    print(model_stats_AICc)
    
    cat("\n\n================ Model comparison table, sorted by AIC ================\n")
    print(model_stats_AIC)
  },
  file = "Habitat_model_comparison_full_output.txt"
)


##自动选择最佳模型

selection_criterion <- "AICc"

## 如果 AICc 全部不可用，则自动改用 AIC
if (!any(is.finite(model_stats[[selection_criterion]]))) {
  warning("AICc 全部不可用，自动改用 AIC 选择最佳模型。")
  selection_criterion <- "AIC"
}

valid_models <- model_stats[is.finite(model_stats[[selection_criterion]]), ]

if (nrow(valid_models) == 0) {
  stop("错误：所有模型的 AIC/AICc 都不可用，无法选择最佳模型。")
}

best_row <- valid_models[which.min(valid_models[[selection_criterion]]), ]
best_model <- best_row$Model
best_fit <- fits[[best_model]]

cat("\n============================================================\n")
cat("最终选择模型：", best_model, "\n")
cat("选择依据：最小 ", selection_criterion, "\n")
cat("log-likelihood：", best_row$log_likelihood, "\n")
cat("AIC：", best_row$AIC, "\n")
cat("AICc：", best_row$AICc, "\n")
cat("free parameters：", best_row$free_parameters, "\n")
cat("free transition rates：", best_row$free_transition_rates, "\n")
cat("============================================================\n")

writeLines(
  c(
    "============================================================",
    "Best model selected for habitat ancestral state reconstruction",
    "============================================================",
    paste("Selection criterion:", selection_criterion),
    paste("Best model:", best_model),
    paste("log-likelihood:", best_row$log_likelihood),
    paste("AIC:", best_row$AIC),
    paste("AICc:", best_row$AICc),
    paste("free parameters:", best_row$free_parameters),
    paste("free transition rates:", best_row$free_transition_rates)
  ),
  con = "Habitat_best_model_selected.txt"
)



##绘制 ER / SYM / ARD 模型拟合图 ####

pdf("Habitat_fitted_models_ER_SYM_ARD.pdf", width = 12, height = 4)

par(mfrow = c(1, 3), mar = c(3, 3, 4, 1))

for (m in models_to_test) {
  
  if (!is.null(fits[[m]])) {
    
    tryCatch(
      {
        plot(fits[[m]], signif = 5)
        title(main = paste0("Fitted ", m, " model"), line = 1)
      },
      error = function(e) {
        plot.new()
        title(main = paste0(m, " plot failed"), line = 1)
        message("模型 ", m, " 绘图失败：", e$message)
      }
    )
    
  } else {
    
    plot.new()
    title(main = paste0(m, " model failed"), line = 1)
  }
}

dev.off()


##使用最佳模型进行 ace 祖先状态重建 ####

fitACE <- ace(
  x = Hab,
  phy = phy,
  model = best_model,
  type = "discrete"
)

cat("\n================ ace 祖先状态重建结果 ================\n")
print(fitACE)

write.csv(
  fitACE$lik.anc,
  file = paste0("Habitat_ancestral_state_likelihood_ace_", best_model, ".csv")
)

capture.output(
  {
    cat("============================================================\n")
    cat("Ancestral state reconstruction using ace\n")
    cat("============================================================\n\n")
    
    cat("Best model used:\n")
    print(best_model)
    
    cat("\nFull ace output:\n")
    print(fitACE)
    
    cat("\nAncestral state likelihoods:\n")
    print(fitACE$lik.anc)
  },
  file = paste0("Habitat_ancestral_state_reconstruction_ace_", best_model, "_output.txt")
)


## 使用最佳模型进行 stochastic character mapping

set.seed(123)

mtrees <- make.simmap(
  tree = phy,
  x = Hab,
  model = best_model,
  nsim = 100,
  message = FALSE
)

pd <- summary(mtrees)

capture.output(
  {
    cat("============================================================\n")
    cat("Summary of stochastic character mapping\n")
    cat("============================================================\n\n")
    
    cat("Best model used:\n")
    print(best_model)
    
    cat("\nSummary:\n")
    print(pd)
  },
  file = paste0("Habitat_simmap_summary_", best_model, "_output.txt")
)

write.csv(
  pd$ace,
  file = paste0("Habitat_simmap_node_posterior_probabilities_", best_model, ".csv")
)

if (!is.null(pd$times)) {
  write.csv(
    pd$times,
    file = paste0("Habitat_simmap_state_times_", best_model, ".csv")
  )
}


## 绘制 100 次 stochastic mapping 结果 ####

pdf(
  paste0("Habitat_simmap_100_replicates_", best_model, ".pdf"),
  width = 12,
  height = 12
)

par(mfrow = c(10, 10), mar = c(0, 0, 0, 0))

invisible(
  sapply(
    mtrees,
    plot,
    colors = colors,
    lwd = 0.4,
    ftype = "off"
  )
)

dev.off()


##绘制 stochastic mapping summary 结果 ####

pdf(
  paste0("Habitat_simmap_summary_", best_model, ".pdf"),
  width = 8,
  height = 10
)

plot(
  pd,
  colors = colors,
  pie = pd$ace,
  cex = 0.5,
  ftype = "off"
)

usr <- par("usr")
x_safe <- usr[1] + 0.05 * (usr[2] - usr[1])
y_safe <- usr[3] + 0.01 * (usr[4] - usr[3])

add.simmap.legend(
  colors = colors,
  prompt = FALSE,
  x = x_safe,
  y = y_safe,
  fsize = 0.8,
  vertical = FALSE
)

dev.off()


## densityMap：仅在二元性状时绘制 ####

if (nlevels(Hab) == 2) {
  
  pdf(
    paste0("Habitat_densityMap_", best_model, ".pdf"),
    width = 8,
    height = 10
  )
  
  densityMap(
    mtrees,
    colors = colors,
    lwd = 1,
    fsize = 0.3,
    outline = FALSE
  )
  
  dev.off()
  
} else {
  
  message(
    "注意：当前 biome 有 ", nlevels(Hab),
    " 个状态，不是二元性状，因此跳过 densityMap。多状态性状建议使用 simmap summary pie 图。"
  )
}


## 保存##

saveRDS(phy, file = "Habitat_matched_phylogeny.rds")
saveRDS(Hab, file = "Habitat_matched_trait_vector.rds")
saveRDS(fits, file = "Habitat_fitDiscrete_ER_SYM_ARD.rds")
saveRDS(fitACE, file = paste0("Habitat_ace_", best_model, ".rds"))
saveRDS(mtrees, file = paste0("Habitat_simmap_", best_model, ".rds"))


cat("\n============================================================\n")
cat("分析完成！\n")
cat("最终使用模型：", best_model, "\n")
cat("选择依据：", selection_criterion, "\n")
cat("模型比较表已输出：Habitat_model_comparison_logLik_AIC_AICc_free_parameters.csv\n")
cat("============================================================\n")