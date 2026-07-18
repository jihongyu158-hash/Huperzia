import os
import re
import subprocess
import shutil
from scipy.stats import chi2
from statsmodels.stats.multitest import multipletests

# =========================
# 基础设置
# =========================
genes = [f for f in os.listdir('.') if f.endswith('.phy')]
tree_file = "tree"


# =========================
# 工具函数
# =========================
def run_codeml(ctl_name, ctl_text):
    """
    指定专属的 .ctl 文件运行 codeml，避免多基因串行或并发时死锁/覆盖
    """
    with open(ctl_name, "w") as f:
        f.write(ctl_text)

    # PAML 支持直接将 ctl 文件作为参数传入
    res = subprocess.run(
        ["codeml", ctl_name],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    return res.returncode == 0


def extract_lnL(text):
    """
    更具鲁棒性地提取 lnL 值的正则
    """
    vals = re.findall(r"lnL.*?(-?\d+\.\d+)", text)
    return float(vals[-1]) if vals else None


def classify_beb(text):
    """
    解析 Bayes Empirical Bayes (BEB) 结果
    """
    strong, moderate = [], []

    if "Bayes Empirical Bayes" not in text:
        return strong, moderate

    sec = text.split("Bayes Empirical Bayes")[1]
    # 兼容有无星号、以及包含可能存在的 mean w 列的排版
    sites = re.findall(r"(\d+)\s+([A-Za-z])\s+(\d+\.\d+)(\*+)?", sec)

    for pos, aa, pp, star in sites:
        pp = float(pp)
        # PAML 习惯：** 代表 PP >= 0.99, * 代表 0.95 <= PP < 0.99
        if pp >= 0.99:
            strong.append(f"{pos}{aa}(PP={pp:.3f})")
        elif pp >= 0.95:
            moderate.append(f"{pos}{aa}(PP={pp:.3f})")

    return strong, moderate


def extract_omega2(text):
    """
    精准抓取 Branch-site Model A 中前景枝的 w2 值
    PAML 典型输出：foreground w     0.00010  1.00000  4.32100  4.32100
    """
    
    match = re.search(r"foreground\s+w?\s+[\d\.]+\s+[\d\.]+\s+([\d\.]+)", text)
    if match:
        return match.group(1)
    
    
    match_alt = re.search(r"w2\s*=\s*([\d\.]+)", text)
    if match_alt:
        return match_alt.group(1)
        
    return "NA"


def clean_paml_junk(gene_name):
    """
    清理 PAML 运行过程中在当前目录下生成的垃圾中间文件，保持工作区干净
    """
    junk_files = ["2NG.D", "2NG.DN", "2NG.EDS", "rst", "rub", "lnf", "rst1"]
    for jf in junk_files:
        if os.path.exists(jf):
            try:
                os.remove(jf)
            except OSError:
                pass
    # 可选：删除专属的 ctl 文件
    for ext in [".alt.ctl", ".null.ctl"]:
        target = f"{gene_name}{ext}"
        if os.path.exists(target):
            os.remove(target)


# =========================
# 主循环
# =========================
for gene in genes:
    name = gene.replace(".phy", "")
    print(f"正在分析基因: {name} ...")

    alt_ctl_fn = f"{name}.alt.ctl"
    null_ctl_fn = f"{name}.null.ctl"
    alt_out = f"{name}.alt.out"
    null_out = f"{name}.null.out"

    # =====================
    # ALT model (branch-site A)
    # =====================
    alt_ctl = f"""
seqfile = {gene}
treefile = {tree_file}
outfile = {alt_out}

noisy = 0
verbose = 0
runmode = 0

seqtype = 1
CodonFreq = 2
clock = 0

model = 2
NSsites = 2

icode = 0
fix_kappa = 0
kappa = 2

fix_omega = 0
omega = 2

cleandata = 1
"""

    ok1 = run_codeml(alt_ctl_fn, alt_ctl)

    # =====================
    # NULL model
    # =====================
    null_ctl = f"""
seqfile = {gene}
treefile = {tree_file}
outfile = {null_out}

noisy = 0
verbose = 0
runmode = 0

seqtype = 1
CodonFreq = 2
clock = 0

model = 2
NSsites = 2

icode = 0
fix_kappa = 0
kappa = 2

fix_omega = 1
omega = 1

cleandata = 1
"""

    ok2 = run_codeml(null_ctl_fn, null_ctl)

    if not (ok1 and ok2):
        print(f"❌ codeml 运行失败，跳过基因: {name}")
        clean_paml_junk(name)
        continue

    # =====================
    # 读取与解析结果
    # =====================
    lnL_alt = lnL_null = None
    strong_beb = []
    moderate_beb = []
    omega2 = "NA"

    if os.path.exists(alt_out):
        with open(alt_out) as f:
            alt_text = f.read()

        lnL_alt = extract_lnL(alt_text)
        strong_beb, moderate_beb = classify_beb(alt_text)
        omega2 = extract_omega2(alt_text)

    if os.path.exists(null_out):
        with open(null_out) as f:
            null_text = f.read()

        lnL_null = extract_lnL(null_text)

    # =====================
    # LRT 似然比检验计算
    # =====================
    if lnL_alt is not None and lnL_null is not None:
        delta = max(0, lnL_alt - lnL_null)
        LRT = 2 * delta
        
        # 使用卡方分布（自由度 df=1）转换为原始 p 值
        p_value = chi2.sf(LRT, df=1) 
        
        summary.append([
            name,          # 0
            lnL_alt,       # 1
            lnL_null,      # 2
            LRT,           # 3
            p_value,       # 4
            None,          # 5 (预留给 q-value)
            None,          # 6 (预留给 FDR 判断)
            omega2,        # 7
            strong_beb,    # 8
            moderate_beb   # 9
        ])

        print(f"✔ 成功解析 -> LRT: {LRT:.3f}, 原始 p-value: {p_value:.5f}, Foreground w2: {omega2}")
    else:
        print(f"⚠️ 结果文件解析异常 (lnL 未抓取成功): {name}")

    # 清理每次运行生成的临时垃圾文件
    clean_paml_junk(name)


# =========================
# FDR 多重假设检验矫正
# =========================
if summary:
    # 提取所有基因的原始 p-value
    pvals = [r[4] for r in summary]

    # 使用 Benjamini-Hochberg (BH) 方法矫正
    reject, qvals, _, _ = multipletests(pvals, method="fdr_bh")

    for i in range(len(summary)):
        summary[i][5] = qvals[i]
        summary[i][6] = "Yes" if reject[i] else "No"

    # =========================
    # 导出统计结果 CSV
    # =========================
    output_fn = "branch_site_summary_FDR.csv"
    with open(output_fn, "w") as f:
        f.write(
            "Gene,lnL_Alt,lnL_Null,LRT,p_value,q_value,FDR_Significant,"
            "Foreground_w2,Strong_BEB,Moderate_BEB\n"
        )

        for r in summary:
            strong = ";".join(r[8]) if r[8] else "None"
            moderate = ";".join(r[9]) if r[9] else "None"

            f.write(
                f"{r[0]},{r[1]},{r[2]},{r[3]:.4f},"
                f"{r[4]:.6e},{r[5]:.6e},{r[6]},"
                f"{r[7]},\"{strong}\",\"{moderate}\"\n"
            )

    print(f"\n🎉 恭喜！数据全部处理完毕。综合报告已保存至: {output_fn}")
else:
    print("\n❌ 错误：未成功获得任何有效基因的分析数据，请检查输入序列或 PAML 路径。")