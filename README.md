# SUDE单细胞降维工具

## 概述

SUDE (Scalable Uniform Distributed Embedding) 是一个专为大规模高维数据设计的降维算法，特别适合单细胞RNA测序数据分析。本工具提供了与scanpy生态系统兼容的SUDE实现。

## 算法特点

### 核心技术
- **共享最近邻 (SNN) 距离修正**：通过分析细胞间的共享邻居信息，优化距离度量
- **Landmark采样策略**：使用PPS (Plum Pudding Sampling) 智能选择代表性细胞点
- **概率分布优化**：通过最小化KL散度学习最优的低维表示
- **分块计算**：支持大数据集的内存高效处理

### 优势
- 🚀 **高可扩展性**：通过landmark采样处理大规模数据集
- 🎯 **保持结构**：同时保持局部和全局数据结构
- 🔧 **灵活性**：多种初始化方法和参数可调
- 🏗️ **集成性**：完美兼容scanpy工作流程

## 安装要求

```bash
# 基本依赖
pip install scanpy pandas numpy scikit-learn scipy matplotlib

# 可选依赖 (用于深度学习集成)
pip install scvi-tools
```

## 快速开始

### 基本用法

```python
import scanpy as sc
from sude_single_cell import SUDE

# 加载数据
adata = sc.datasets.pbmc3k_processed()

# 确保有邻居图 (通常在预处理中计算)
sc.pp.neighbors(adata, use_rep='X_pca')

# 应用SUDE降维
SUDE(adata, use_rep='X_pca', n_comps=2)

# 可视化结果
sc.pl.embedding(adata, basis='sude', color='louvain')
```

### 完整的单细胞分析流程

```python
import scanpy as sc
from sude_single_cell import SUDE

# 1. 数据加载和预处理
adata = sc.read_h5ad('your_data.h5ad')

# 2. 质量控制
sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# 3. 标准化和特征选择
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata.raw = adata
adata = adata[:, adata.var.highly_variable]

# 4. 降维
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, use_rep='X_pca')

# 5. SUDE降维
SUDE(adata, use_rep='X_pca', n_comps=2, k_neighbors=20)

# 6. 聚类和可视化
sc.tl.leiden(adata)
sc.pl.embedding(adata, basis='sude', color='leiden')
```

## 参数详解

### 主要参数

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `n_comps` | int | 2 | 降维后的维度数 |
| `use_rep` | str | None | 使用的数据表示，如'X_pca'、'X_scvi' |
| `k_neighbors` | int | 20 | PPS采样的邻居数，影响landmark选择 |
| `large_data` | bool | False | 是否为大数据集启用分块计算 |
| `initialization` | str | 'le' | 初始化方法：'le'、'pca'、'mds' |
| `aggregation_coef` | float | 1.2 | 聚合系数，控制SNN权重 |
| `max_epochs` | int | 50 | 最大训练轮数 |

### 参数调优建议

#### 数据大小相关
- **小数据集** (<5000 细胞): 使用默认参数
- **中等数据集** (5000-20000 细胞): `k_neighbors=30`
- **大数据集** (>20000 细胞): `large_data=True`, `k_neighbors=50`

#### 初始化方法选择
- **'le' (拉普拉斯特征映射)**: 最好的局部结构保持，推荐用于大多数情况
- **'pca' (主成分分析)**: 计算最快，适合初步探索
- **'mds' (多维标度)**: 最好的距离保持，计算较慢

#### 聚合系数调节
- **0.8-1.0**: 更重视原始距离
- **1.2-1.5**: 平衡原始距离和共享邻居信息（推荐）
- **1.5-2.0**: 更重视共享邻居信息

## 高级用法

### 1. 大数据集优化

```python
# 对于超大数据集 (>50000 细胞)
SUDE(
    adata, 
    use_rep='X_pca', 
    n_comps=2,
    k_neighbors=50,           # 增加邻居数以提高代表性
    large_data=True,          # 启用分块计算
    max_epochs=30,            # 减少训练轮数以加快速度
    aggregation_coef=1.5      # 增强邻居信息权重
)
```

### 2. 与深度学习方法集成

```python
import scvi

# 使用scVI进行降维
scvi.model.SCVI.setup_anndata(adata)
model = scvi.model.SCVI(adata, n_latent=30)
model.train()
adata.obsm["X_scvi"] = model.get_latent_representation()

# 在scVI表示基础上应用SUDE
SUDE(adata, use_rep='X_scvi', n_comps=2)
```

### 3. 批次效应处理

```python
# 对于有批次效应的数据
# 先使用批次矫正方法 (如Harmony, scVI等)
import scanpy.external as sce
sce.pp.harmony_integrate(adata, key='batch')

# 然后应用SUDE
SUDE(adata, use_rep='X_pca_harmony', n_comps=2)
```

### 4. 参数网格搜索

```python
# 系统性地测试不同参数组合
import itertools

k_neighbors_list = [20, 30, 50]
agg_coef_list = [1.0, 1.2, 1.5]
init_methods = ['le', 'pca']

best_score = -np.inf
best_params = None

for k_neighbors, agg_coef, init_method in itertools.product(
    k_neighbors_list, agg_coef_list, init_methods
):
    print(f"测试参数: k_neighbors={k_neighbors}, agg_coef={agg_coef}, init={init_method}")
    
    # 应用SUDE
    adata_test = adata.copy()
    SUDE(
        adata_test, 
        use_rep='X_pca',
        k_neighbors=k_neighbors,
        aggregation_coef=agg_coef,
        initialization=init_method,
        key_added=f'X_sude_test'
    )
    
    # 评估质量 (例如使用聚类质量指标)
    # score = evaluate_embedding_quality(adata_test)
    # if score > best_score:
    #     best_score = score
    #     best_params = (k_neighbors, agg_coef, init_method)
```

## 结果解释和评估

### 1. 可视化比较

```python
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# 原始方法比较
sc.pl.pca(adata, color='leiden', ax=axes[0], show=False)
axes[0].set_title('PCA')

sc.pl.umap(adata, color='leiden', ax=axes[1], show=False)  
axes[1].set_title('UMAP')

sc.pl.embedding(adata, basis='sude', color='leiden', ax=axes[2], show=False)
axes[2].set_title('SUDE')

plt.tight_layout()
plt.show()
```

### 2. 定量评估

```python
from sklearn.metrics import adjusted_rand_score, silhouette_score

# 聚类一致性评估
ari_score = adjusted_rand_score(adata.obs['leiden'], adata.obs['true_labels'])
print(f"调整后兰德指数: {ari_score:.3f}")

# 轮廓系数评估
sil_score = silhouette_score(adata.obsm['X_sude'], adata.obs['leiden'])
print(f"轮廓系数: {sil_score:.3f}")
```

## 故障排除

### 常见问题

1. **内存错误**
   ```python
   # 解决方案：启用大数据模式
   SUDE(adata, large_data=True)
   ```

2. **计算时间过长**
   ```python
   # 解决方案：减少训练轮数和增加邻居数
   SUDE(adata, max_epochs=20, k_neighbors=50)
   ```

3. **结果质量不佳**
   ```python
   # 解决方案：调整聚合系数和初始化方法
   SUDE(adata, aggregation_coef=1.5, initialization='le')
   ```

### 性能优化建议

- 对于大数据集，建议先用PCA降维到50-100维再应用SUDE
- 使用适当的`k_neighbors`值：太小会丢失全局结构，太大会增加计算成本
- 对于探索性分析，可以先用较少的训练轮数快速获得结果

## 算法原理

SUDE算法包含以下关键步骤：

1. **距离修正**: 使用共享最近邻信息修正原始欧几里得距离
2. **Landmark选择**: 通过PPS算法选择具有高邻居密度的代表性点
3. **概率矩阵构建**: 构建高维空间的概率分布矩阵P
4. **梯度优化**: 通过最小化KL散度优化低维表示
5. **非Landmark映射**: 使用CLLE将剩余点映射到低维空间

## 引用

如果您在研究中使用了SUDE算法，请考虑引用相关论文。

## 贡献和反馈

欢迎提交Issue和Pull Request来改进这个工具！

---

**更多信息和高级用法请参考示例文件 `sude_example.py`** 