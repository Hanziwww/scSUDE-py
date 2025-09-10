import numpy as np
import pandas as pd
from sude_py.sude import sude
import scanpy as sc
from typing import Optional, Union, Literal
import warnings


def SUDE(
    adata,
    n_comps: int = 2,
    use_rep: Optional[str] = None,
    key_added: str = 'X_sude',
    k_neighbors: int = 20,
    large_data: bool = False,
    initialization: Literal['le', 'pca', 'mds'] = 'le',
    aggregation_coef: float = 1.2,
    max_epochs: int = 50,
    normalize: bool = True,
    copy: bool = False,
    random_state: Optional[int] = None
):
    """
    应用SUDE (Scalable Uniform Distributed Embedding) 降维算法到单细胞数据
    
    SUDE是一个基于流形学习的降维算法，特别适合处理大规模高维单细胞数据。
    该算法通过共享最近邻信息修正距离度量，并使用landmark采样策略提高可扩展性。
    
    Parameters
    ----------
    adata : AnnData
        带注释的数据矩阵，单细胞表达数据
    n_comps : int, default: 2
        降维后的维度数
    use_rep : str, optional
        使用的表示，如 'X_pca', 'X_scvi'等。如果为None，使用adata.X
    key_added : str, default: 'X_sude'
        结果存储在adata.obsm中的键名
    k_neighbors : int, default: 20
        PPS采样中使用的近邻数量，影响landmark选择
    large_data : bool, default: False
        是否为大数据集，如果为True，使用分块计算避免内存溢出
    initialization : {'le', 'pca', 'mds'}, default: 'le'
        初始化方法:
        - 'le': 拉普拉斯特征映射 (推荐用于保持局部结构)
        - 'pca': 主成分分析 (计算快速)
        - 'mds': 多维标度 (保持距离关系)
    aggregation_coef : float, default: 1.2
        聚合系数，控制共享邻居信息的权重
    max_epochs : int, default: 50
        最大训练轮数
    normalize : bool, default: True
        是否对数据进行min-max标准化
    copy : bool, default: False
        是否返回数据副本
    random_state : int, optional
        随机种子，用于结果可重现性
        
    Returns
    -------
    adata : AnnData
        如果copy=True，返回数据副本；否则就地修改adata
        降维结果存储在adata.obsm[key_added]中
        
    Examples
    --------
    基本用法:
    >>> import scanpy as sc
    >>> from sude_single_cell import SUDE
    >>> 
    >>> # 加载数据并预处理
    >>> adata = sc.datasets.pbmc3k_processed()
    >>> 
    >>> # 使用PCA表示进行SUDE降维
    >>> sc.pp.neighbors(adata, use_rep='X_pca')
    >>> SUDE(adata, use_rep='X_pca', n_comps=2)
    >>> 
    >>> # 可视化结果
    >>> sc.pl.embedding(adata, basis='sude')
    
    大数据集用法:
    >>> # 对于大于10000个细胞的数据集
    >>> SUDE(adata, use_rep='X_pca', large_data=True, k_neighbors=30)
    
    使用scVI表示:
    >>> # 如果已经有scVI等深度学习降维表示
    >>> SUDE(adata, use_rep='X_scvi', n_comps=2)
    
    Notes
    -----
    - 建议在运行SUDE之前先进行PCA降维或使用其他降维表示(如scVI)
    - 对于大数据集(>10000细胞)，建议设置large_data=True
    - k_neighbors参数会影响landmark选择，较大的值会选择更少的landmark
    - aggregation_coef控制共享邻居信息的重要性，通常在1.0-2.0之间
    """
    
    if copy:
        adata = adata.copy()
    
    # 设置随机种子
    if random_state is not None:
        np.random.seed(random_state)
    
    # 获取输入数据
    if use_rep is None:
        X = adata.X
        if hasattr(X, 'toarray'):  # 如果是稀疏矩阵
            X = X.toarray()
        rep_used = 'X'
    else:
        if use_rep not in adata.obsm:
            raise ValueError(f"指定的表示 '{use_rep}' 不存在于 adata.obsm 中")
        X = adata.obsm[use_rep]
        rep_used = use_rep
    
    # 检查数据格式
    if not isinstance(X, np.ndarray):
        X = np.array(X)
    
    # 数据验证
    if X.shape[0] != adata.n_obs:
        raise ValueError(f"输入数据的样本数 ({X.shape[0]}) 与adata的观测数 ({adata.n_obs}) 不匹配")
    
    if n_comps >= X.shape[1]:
        raise ValueError(f"降维目标维度 ({n_comps}) 不能大于等于输入数据维度 ({X.shape[1]})")
    
    # 数据大小提示
    n_cells, n_features = X.shape
    print(f"对 {n_cells} 个细胞、{n_features} 个特征进行SUDE降维...")
    print(f"使用表示: {rep_used}")
    print(f"目标维度: {n_comps}")
    
    # 对于大数据集的建议
    if n_cells > 10000 and not large_data:
        warnings.warn(
            f"检测到大数据集 ({n_cells} 个细胞)。建议设置 large_data=True 以优化内存使用。",
            UserWarning
        )
    
    # 自适应调整k_neighbors
    if k_neighbors >= n_cells:
        k_neighbors = min(n_cells - 1, 50)
        warnings.warn(
            f"k_neighbors 过大，已调整为 {k_neighbors}",
            UserWarning
        )
    
    # 尝试复用预计算的邻居图（来自 sc.pp.neighbors）
    precomputed_knn = None
    try:
        if 'neighbors' in adata.uns and ('distances' in adata.obsp or 'connectivities' in adata.obsp):
            use_dist = 'distances' in adata.obsp
            mat = adata.obsp['distances'] if use_dist else adata.obsp['connectivities']
            if hasattr(mat, 'tocsr'):
                mat = mat.tocsr()
            precomputed_knn = np.zeros((n_cells, k_neighbors + 1), dtype=int)
            for i in range(n_cells):
                row_start = mat.indptr[i]
                row_end = mat.indptr[i + 1]
                cols = mat.indices[row_start:row_end]
                data = mat.data[row_start:row_end]
                # ensure self included
                if i not in cols:
                    cols = np.append(cols, i)
                    data = np.append(data, 0.0 if use_dist else (np.max(data) if data.size > 0 else 1.0))
                # select top k by rule: distances ascending, connectivities descending
                order = np.argsort(data) if use_dist else np.argsort(-data)
                top_idx = cols[order][:k_neighbors + 1] if order.size > 0 else np.array([i])
                # pad if not enough neighbors
                if top_idx.size < k_neighbors + 1:
                    pad = np.full(k_neighbors + 1 - top_idx.size, i, dtype=int)
                    top_idx = np.concatenate([top_idx, pad])
                precomputed_knn[i] = top_idx
    except Exception as e:
        warnings.warn(f"复用预计算邻居失败，转为内部计算。原因: {e}")
        precomputed_knn = None
    
    try:
        # 执行SUDE降维
        print("开始SUDE降维计算...")
        Y = sude(
            X=X,
            no_dims=n_comps,
            k1=k_neighbors,
            normalize=normalize,
            large=large_data,
            initialize=initialization,
            agg_coef=aggregation_coef,
            T_epoch=max_epochs,
            precomputed_knn=precomputed_knn
        )
        
        # 存储结果
        adata.obsm[key_added] = Y
        
        # 记录参数信息
        adata.uns[f'{key_added}_params'] = {
            'n_comps': n_comps,
            'use_rep': use_rep,
            'k_neighbors': k_neighbors,
            'large_data': large_data,
            'initialization': initialization,
            'aggregation_coef': aggregation_coef,
            'max_epochs': max_epochs,
            'normalize': normalize,
            'random_state': random_state,
            'rep_used': rep_used
        }
        
        # 注册embedding以便scanpy绘图
        basis_name = key_added.replace('X_', '') if key_added.startswith('X_') else key_added
        adata.obsm[f'X_{basis_name}'] = Y
        adata.obsm_keys()
        
        if copy:
            return adata
        return None
    except Exception as e:
        raise RuntimeError(f"SUDE降维过程中发生错误: {e}") from e


def plot_sude_embedding(
    adata,
    color=None,
    basis='sude',
    **kwargs
):
    """
    可视化SUDE降维结果的便捷函数
    
    Parameters
    ----------
    adata : AnnData
        包含SUDE结果的数据
    color : str or list, optional
        用于着色的变量
    basis : str, default: 'sude'
        embedding的basis名称
    **kwargs
        传递给scanpy.pl.embedding的其他参数
    """
    if f'X_{basis}' not in adata.obsm:
        available_bases = [key for key in adata.obsm.keys() if key.startswith('X_')]
        raise ValueError(f"没有找到 'X_{basis}' embedding。可用的embedding: {available_bases}")
    
    sc.pl.embedding(adata, basis=basis, color=color, **kwargs)


# 别名函数，与用户期望的用法一致
def sude_embedding(adata, **kwargs):
    """
    SUDE降维的别名函数，与用户期望的用法一致
    """
    return SUDE(adata, **kwargs) 