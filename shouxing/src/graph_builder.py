import networkx as nx
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
import pickle

class GraphBuilder:
    def __init__(self, data_file='data/combined_data.csv'):
        self.data_file = data_file
        self.graph = nx.Graph()

    def load_data(self):
        """加载数据"""
        try:
            self.df = pd.read_csv(self.data_file)
            print(f"加载了 {len(self.df)} 个分子")
            return True
        except FileNotFoundError:
            print("数据文件不存在，请先运行数据爬取")
            return False

    def calculate_chiral_similarity(self, mol1, mol2):
        """计算两个分子的手性相似性"""
        if not mol1 or not mol2:
            return 0
        
        # 获取手性中心
        chiral1 = Chem.FindMolChiralCenters(mol1, includeUnassigned=True)
        chiral2 = Chem.FindMolChiralCenters(mol2, includeUnassigned=True)
        
        if not chiral1 and not chiral2:
            return 1  # 都没有手性中心，相似
        if not chiral1 or not chiral2:
            return 0  # 一个有，一个没有，不相似
        
        # 计算手性中心数量相似性
        count_sim = 1 - abs(len(chiral1) - len(chiral2)) / max(len(chiral1), len(chiral2))
        
        # 计算Morgan指纹相似性
        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=1024)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=1024)
        fp_sim = Chem.DataStructs.TanimotoSimilarity(fp1, fp2)
        
        return (count_sim + fp_sim) / 2

    def build_graph(self):
        """构建图谱"""
        if not hasattr(self, 'df'):
            if not self.load_data():
                return
        
        # 添加节点
        for idx, row in self.df.iterrows():
            smiles = row['smiles']
            if pd.isna(smiles):
                continue
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                node_id = f"{row['source']}_{row.get('cid', row.get('chembl_id', idx))}"
                self.graph.add_node(node_id, 
                                  smiles=smiles,
                                  name=row.get('name', ''),
                                  chiral_centers=row['chiral_centers'],
                                  source=row['source'])
        
        # 添加边（基于手性相似性）
        nodes = list(self.graph.nodes())
        for i in range(len(nodes)):
            for j in range(i+1, len(nodes)):
                node1, node2 = nodes[i], nodes[j]
                smiles1 = self.graph.nodes[node1]['smiles']
                smiles2 = self.graph.nodes[node2]['smiles']
                
                mol1 = Chem.MolFromSmiles(smiles1)
                mol2 = Chem.MolFromSmiles(smiles2)
                
                similarity = self.calculate_chiral_similarity(mol1, mol2)
                
                if similarity > 0.5:  # 相似性阈值
                    self.graph.add_edge(node1, node2, weight=similarity)
        
        print(f"构建图谱完成：{len(self.graph.nodes())} 个节点，{len(self.graph.edges())} 条边")

    def save_graph(self, filename='data/molecule_graph.pkl'):
        """保存图谱"""
        with open(filename, 'wb') as f:
            pickle.dump(self.graph, f)
        print(f"图谱已保存到 {filename}")

    def load_graph(self, filename='data/molecule_graph.pkl'):
        """加载图谱"""
        try:
            with open(filename, 'rb') as f:
                self.graph = pickle.load(f)
            print(f"图谱已加载：{len(self.graph.nodes())} 个节点")
            return True
        except FileNotFoundError:
            print("图谱文件不存在")
            return False