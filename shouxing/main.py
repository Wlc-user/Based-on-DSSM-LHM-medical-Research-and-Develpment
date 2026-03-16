from src.data_crawler import DataCrawler
from src.graph_builder import GraphBuilder
from src.similarity_filter import SimilarityFilter
import os

def main():
    print("化学数据处理系统启动")
    
    # 1. 数据层：爬取数据
    crawler = DataCrawler()
    
    # 示例：爬取一些常见化合物
    queries = ['aspirin', 'caffeine', 'glucose']
    
    for query in queries:
        print(f"\n爬取 {query} 数据...")
        crawler.crawl_pubchem(query, limit=20)
        crawler.crawl_chembl(query, limit=20)
    
    # 合并数据
    combined_df = crawler.combine_data()
    print(f"总共收集了 {len(combined_df)} 个分子")
    
    # 2. 图谱层：构建图谱
    builder = GraphBuilder()
    builder.build_graph()
    builder.save_graph()
    
    # 3. 应用层：筛选相似分子
    filter = SimilarityFilter()
    
    # 示例：查找与阿司匹林相似的分子
    target_smiles = 'CC(=O)OC1=CC=CC=C1C(=O)O'  # 阿司匹林的SMILES
    similar = filter.find_similar_molecules(target_smiles, top_k=5)
    
    print(f"\n与目标分子最相似的5个分子：")
    for node, sim, data in similar:
        print(f"- {data['name']} (相似度: {sim:.3f}, 手性中心: {data['chiral_centers']})")
    
    # 筛选具有手性特征的分子
    chiral_molecules = filter.filter_by_chiral_features(min_chiral_centers=1)
    print(f"\n具有至少1个手性中心的分子数量: {len(chiral_molecules)}")

if __name__ == "__main__":
    main()