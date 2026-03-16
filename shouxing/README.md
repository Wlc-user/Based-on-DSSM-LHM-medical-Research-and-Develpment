# 化学数据处理系统

## 项目概述
本项目实现了一个三层架构的化学数据处理系统：
1. **数据层**：爬取PubChem和ChEMBL的公开数据
2. **图谱层**：使用手性特征进行节点关联，构建分子相似性图谱
3. **应用层**：通过图谱关系快速筛选相似结构的分子

## 安装依赖
```bash
pip install -r requirements.txt
```

## 使用方法

### 运行完整流程
```bash
python main.py
```

### 单独使用各模块

#### 数据爬取
```python
from src.data_crawler import DataCrawler

crawler = DataCrawler()
crawler.crawl_pubchem('aspirin', limit=50)
crawler.crawl_chembl('caffeine', limit=50)
crawler.combine_data()
```

#### 图谱构建
```python
from src.graph_builder import GraphBuilder

builder = GraphBuilder()
builder.build_graph()
builder.save_graph()
```

#### 相似性筛选
```python
from src.similarity_filter import SimilarityFilter

filter = SimilarityFilter()
similar_molecules = filter.find_similar_molecules('CC(=O)OC1=CC=CC=C1C(=O)O', top_k=10)
for node, sim, data in similar_molecules:
    print(f"{data['name']}: {sim:.3f}")
```

## 项目结构
- `data/` : 存储爬取的数据和图谱文件
- `src/` : 源代码
  - `data_crawler.py` : 数据爬取模块
  - `graph_builder.py` : 图谱构建模块
  - `similarity_filter.py` : 相似性筛选模块
- `main.py` : 主入口文件
- `requirements.txt` : 依赖列表

## 技术栈
- **数据爬取**: PubChemPy, ChEMBL Web Resource Client
- **化学计算**: RDKit
- **图谱处理**: NetworkX
- **数据处理**: Pandas, NumPy

## 注意事项
- PubChem和ChEMBL的数据量很大，建议设置合理的limit参数
- 图谱构建对计算资源要求较高，大数据集需要优化
- 相似性计算基于Morgan指纹和手性中心数量的组合