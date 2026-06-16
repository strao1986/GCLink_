import scanpy as sc
import numpy as np
import pandas as pd
from scTenifold import scTenifoldKnk

### ACO2 gene in lung tissue ###
#adata = sc.read_h5ad("/public/jiangjw/02.anxiety_ADs_EUR/10.MAGMA_Celltyping/sc_dataset/lung/lung.cellxgene.h5ad")

#print(np.max(adata.X)) #如果最大值 > 50 或 >100说明是 raw counts

#adata_counts = adata.raw.to_adata().copy()

## 筛选细胞类型
#celltypes = ["T_CD8_CytT", "Plasma_cells", "Mast_cells"]

#adata_counts = adata_counts[adata_counts.obs['Celltypes_updated_July_2020'] == "T_CD8_CytT"].copy()

## 高变基因筛选
#adata_tmp = adata_counts.copy()
#sc.pp.normalize_total(adata_tmp, target_sum=1e4)
#sc.pp.log1p(adata_tmp)
#sc.pp.highly_variable_genes(
#     adata_tmp,
#     n_top_genes=3000,
#     flavor="seurat"
# )
#hvg_genes = adata_tmp.var[adata_tmp.var.highly_variable].index
#hvg_genes = list(hvg_genes)
#hvg_genes.append("ACO2")
#print("HVG number:", len(hvg_genes))
#adata_hvg = adata_counts[:, hvg_genes].copy()
#print(adata_hvg)

## 转置为行是细胞，列是基因
#X = adata_hvg.X.toarray() if not isinstance(adata_hvg.X, np.ndarray) else adata_hvg.X
#df = pd.DataFrame(
#     X.T,
#     index=adata_hvg.var_names,   # gene
#     columns=adata_hvg.obs_names  # cell
# )
#print(df.shape)

## 虚拟敲除
#sc_knk = scTenifoldKnk(
#     data=df,
#     ko_method="default",
#     ko_genes=["ACO2"],
#     qc_kws={
#         "min_lib_size": 10,
#         "min_percent": 0.0005
#     }
# )
#result = sc_knk.build()
#result.to_csv("/public/jiangjw/02.anxiety_ADs_EUR/14.Virtual_Knockout/ACO2_KO_T_CD8_results.csv")


## 循环处理 ##
# 读取数据
adata = sc.read_h5ad("/public/jiangjw/02.anxiety_ADs_EUR/10.MAGMA_Celltyping/sc_dataset/lung/lung.cellxgene.h5ad")

# 使用 raw counts
adata_counts = adata.raw.to_adata().copy()

# 需要分析的细胞类型
#celltypes = ["T_CD8_CytT", "Plasma_cells", "Mast_cells"]
celltypes = ["Mast_cells"]

for ct in celltypes:
    
    print(f"\n===== Processing {ct} =====")
    
    # 1️⃣ 筛选细胞类型
    adata_ct = adata_counts[
        adata_counts.obs['Celltypes_updated_July_2020'] == ct
    ].copy()
    #adata_ct = adata_counts.copy()
    
    print("Cell number:", adata_ct.n_obs)
    
    # 如果细胞太少，跳过
    if adata_ct.n_obs < 200:
        print(f"{ct} cell number too small, skip.")
        continue
    
    # 2️⃣ 计算高变基因（只用于筛选，不影响原始count）
    adata_tmp = adata_ct.copy()
    sc.pp.normalize_total(adata_tmp, target_sum=1e4)
    sc.pp.log1p(adata_tmp)
    sc.pp.highly_variable_genes(
        adata_tmp,
        n_top_genes=5000,
        flavor="seurat"
    )
    
    hvg_genes = list(
        adata_tmp.var[adata_tmp.var.highly_variable].index
    )
    
    # 强制加入 ACO2（如果存在）
    if "ACO2" in adata_ct.var_names:
        if "ACO2" not in hvg_genes:
            hvg_genes.append("ACO2")
    else:
        print("ACO2 not in dataset, skip.")
        continue
    
    print("HVG number:", len(hvg_genes))
    
    # 3️⃣ 取原始 count
    adata_hvg = adata_ct[:, hvg_genes].copy()
    
    X = adata_hvg.X.toarray() if not isinstance(
        adata_hvg.X, np.ndarray
    ) else adata_hvg.X
    
    df = pd.DataFrame(
        X.T,
        index=adata_hvg.var_names,
        columns=adata_hvg.obs_names
    )
    
    print("Matrix shape:", df.shape)
    
    # 4️⃣ 虚拟敲除
    sc_knk = scTenifoldKnk(
        data=df,
        ko_method="default",
        ko_genes=["ACO2"],
        qc_kws={
            "min_lib_size": 10,
            "min_percent": 0.0005
        }
    )
    
    result = sc_knk.build()
    
    # 5️⃣保存结果
    output_path = f"/public/jiangjw/02.anxiety_ADs_EUR/14.Virtual_Knockout/ACO2_KO_{ct}_HVG5000_results.csv"
    result.to_csv(output_path)
    
    print(f"Saved: {output_path}")
