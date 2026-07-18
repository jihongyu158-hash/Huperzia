library(ape)
library(BAMMtools)
library(coda)

tree <- read.tree("23species_cds_beast_time_tree_onlyHuperzia.tree")
plot(tree, cex = 0.15)
nodelabels(frame = "none", cex = 0.15)
save.image("tree.rda")

#BAMM计算先验值，更改control文件参数
priors <- setBAMMpriors(phy = tree, total.taxa = 23, outfile = NULL)

#读取结果文件
mcmcout <- read.csv("mcmc_out.txt", header = TRUE)
ed <- getEventData(tree, "event_data.txt", burnin = 0.10)

#验证收敛性
plot(mcmcout$logLik ~ mcmcout$generation)
burnstart <- floor(0.1 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]

effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)

#检查
check_sample_adequacy <- function(event_data_file) 
  ed <- getEventData(tree, eventdata = event_data_file, burnin = 0.1)
nsamples <- length(ed$numberEvents)

print(paste("可用样本点数:", nsamples))

if(nsamples < 1000) {
  print("警告: 样本点可能不足，建议降低eventDataWriteFreq")
} else if(nsamples > 10000) {
  print("样本点充足，但文件可能过大，可考虑提高eventDataWriteFreq")
} else {
  print("样本点数量适中")
}

#node检查
clade1 <- extract.clade(tree, node=36)
clade2 <- extract.clade(tree, node=28)

print("Clade 1 tips:")
print(clade1$tip.label)
print("Clade 2 tips:")
print(clade2$tip.label)

# 方法1：分别绘制两个节点的速率，使用不同的颜色
par(mfrow = c(1, 1)) # 重置图形参数

# 绘制整个树的背景（浅灰色）
plotRateThroughTime(ed, ratetype = "speciation", 
                    intervalCol = "lightgray", avgCol = "lightgray",
                    ylim = c(0, 0.5)) # 根据您的数据调整y轴范围

# 叠加支系A的速率（红色）
plotRateThroughTime(ed, ratetype = "speciation", 
                    node = 36, nodetype = "include",
                    intervalCol = NA, avgCol = "red", 
                    add = TRUE)

# 叠加支系B的速率（蓝色）
plotRateThroughTime(ed, ratetype = "speciation", 
                    node = 28, nodetype = "exclude",
                    intervalCol = NA, avgCol = "blue", 
                    add = TRUE)

# 添加图例
legend("topright", legend = c("支系A", "支系B"), 
       col = c("red", "blue"), lwd = 2)

# 使用 getCladeRates 进行定量比较
rates_A <- getCladeRates(ed, node = 36)
rates_B <- getCladeRates(ed, node = 28)

#计算所有物种的形成速率
rates <- getCladeRates(ed)
cat("平均物种形成速率:", mean(rates$lambda), "\n")
cat("平均物种形成速率:", mean(rates$mu), "\n")
netdiv <- rates$lambda - rates$mu
cat("平均净多样化速率:", mean(netdiv), "\n")

# 比较物种形成速率
cat("支系A平均物种形成速率:", mean(rates_A$lambda), "\n")
cat("支系B平均物种形成速率:", mean(rates_B$lambda), "\n")

# 比较灭绝速率
cat("支系A平均灭绝速率:", mean(rates_A$mu), "\n")
cat("支系B平均灭绝速率:", mean(rates_B$mu), "\n")

# 比较净多样化速率
netdiv_A <- rates_A$lambda - rates_A$mu
netdiv_B <- rates_B$lambda - rates_B$mu

cat("支系A平均净多样化速率:", mean(netdiv_A), "\n")
cat("支系B平均净多样化速率:", mean(netdiv_B), "\n")

# 统计检验：两个支系速率是否显著不同？
t_test_lambda <- t.test(rates_A$lambda, rates_B$lambda)
t_test_netdiv <- t.test(netdiv_A, netdiv_B)
t_test_mu <- t.test(rates_A$mu, rates_B$mu)
cat("物种形成速率t检验p值:", t_test_lambda$p.value, "\n")
cat("净多样化速率t检验p值:", t_test_netdiv$p.value, "\n")
cat("物种灭绝速率t检验p值:", t_test_mu$p.value, "\n")

#绘制BAMM物种多样化速率图
bamm.tree <- plot(ed,lwd=2,labels = T,cex = 0.35)
addBAMMshifts(ed,cex = 2)
addBAMMlegend(bamm.tree,location = c(-30,-10,-15,25),nTicks = 6,side = 4,las=1)

#绘制所有物种的物种多样化速率图
par(new=TRUE)
plotRateThroughTime(ed, ratetype="speciation",intervalCol="#E72D00", avgCol="#E72D00")
plotRateThroughTime(ed, ratetype="extinction", intervalCol="#1F00A2", avgCol="#1F00A2",add = TRUE)
plotRateThroughTime(ed, ratetype="netdiv",intervalCol="#29A200", avgCol="#29A200",add = TRUE)

#绘制各区域的物种多样化速率图
plot.new()
par(mfrow=c(3,3))
st <- max(branching.times(tree))
#speciation
plotRateThroughTime(ed, ratetype="speciation", intervalCol="#E72D00", avgCol="#E72D00", start.time=st, ylim=c(0,0.5), cex.axis=2)
text(x=30, y= 0.4, label="All", font=4, cex=2.0, pos=4)
plotRateThroughTime(ed, ratetype="speciation", intervalCol="#E72D00", avgCol="#E72D00", start.time=st, node=28, ylim=c(0,0.5),cex.axis=1.5)
text(x=30, y= 0.4, label="EBLFs", font=4, cex=2.0, pos=4)
plotRateThroughTime(ed, ratetype="speciation", intervalCol="#E72D00", avgCol="#E72D00", start.time=st, node=36, ylim=c(0,0.5),cex.axis=1.5)
text(x=30, y= 0.4, label="MNLFs", font=4, cex=2.0, pos=4)

#extinction
plotRateThroughTime(ed, ratetype="extinction", intervalCol="#1F00A2", avgCol="#1F00A2", start.time=st, ylim=c(0,0.5), cex.axis=2)
text(x=30, y= 0.4, label="All", font=4, cex=2.0, pos=4)
plotRateThroughTime(ed, ratetype="extinction", intervalCol="#1F00A2", avgCol="#1F00A2", start.time=st, node=28, ylim=c(0,0.5),cex.axis=1.5)
text(x=30, y= 0.4, label="EBLFs", font=4, cex=2.0, pos=4)
plotRateThroughTime(ed, ratetype="extinction", intervalCol="#1F00A2", avgCol="#1F00A2", start.time=st, node=36, ylim=c(0,0.5),cex.axis=1.5)
text(x=30, y= 0.4, label="MNLFs", font=4, cex=2.0, pos=4)

#netdiv
plotRateThroughTime(ed, ratetype="netdiv", intervalCol="#29A200", avgCol="#29A200", start.time=st, ylim=c(0,0.5), cex.axis=2)
text(x=30, y= 0.4, label="All", font=4, cex=2.0, pos=4)
plotRateThroughTime(ed, ratetype="netdiv", intervalCol="#29A200", avgCol="#29A200", start.time=st, node=28, ylim=c(0,0.5),cex.axis=1.5)
text(x=30, y= 0.4, label="EBLFs", font=4, cex=2.0, pos=4)
plotRateThroughTime(ed, ratetype="netdiv", intervalCol="#29A200", avgCol="#29A200", start.time=st, node=36, ylim=c(0,0.5),cex.axis=1.5)
text(x=30, y= 0.4, label="MNLFs", font=4, cex=2.0, pos=4)