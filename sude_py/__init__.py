"""
SUDE (Scalable Uniform Distributed Embedding) 降维算法包

这个包实现了SUDE算法，包含以下主要组件：
- sude: 主要的SUDE降维函数
- pps: Plum Pudding Sampling landmark选择
- learning_s: 小数据集学习算法
- learning_l: 大数据集学习算法
- clle: 约束局部线性embedding
- init_pca, pca, mds: 各种初始化方法
- opt_scale: 尺度优化
"""

from .sude import sude

__version__ = "1.0.0"
__author__ = "SUDE Contributors"

__all__ = ['sude'] 