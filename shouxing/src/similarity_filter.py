import networkx as nx
from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
import pickle

class SimilarityFilter:
    def __init__(self, graph_file='data/molecule_graph.pkl'):
        self.graph_file = graph_file
        self.graph = None

    def load_graph(self):
        """加载图谱"""
        try:
            with open(self.graph_file, 'rb') as f:
                self.graph = pickle.load(f)
            print(f"图谱已加载：{len(self.graph.nodes())} 个节点")
            return True
        except FileNotFoundError:
            print("图谱文件不存在，请先构建图谱")
            return False

    def find_similar_molecules(self, target_smiles, top_k=10):
        """通过图谱关系找到相似分子"""
        if not self.graph:
            if not self.load_graph():
                return []
        
        target_mol = Chem.MolFromSmiles(target_smiles)
        if not target_mol:
            print("无效的SMILES字符串")
            return []
        
        # 计算目标分子的手性中心
        target_chiral = Chem.FindMolChiralCenters(target_mol, includeUnassigned=True)
        
        # 找到图中与目标最相似的节点
        similarities = []
        for node, data in self.graph.nodes(data=True):
            node_smiles = data['smiles']
            node_mol = Chem.MolFromSmiles(node_smiles)
            if node_mol:
                # 计算相似性
                fp1 = AllChem.GetMorganFingerprintAsBitVect(target_mol, 2, nBits=1024)
                fp2 = AllChem.GetMorganFingerprintAsBitVect(node_mol, 2, nBits=1024)
                tanimoto_sim = Chem.DataStructs.TanimotoSimilarity(fp1, fp2)
                
                # 手性相似性
                node_chiral = Chem.FindMolChiralCenters(node_mol, includeUnassigned=True)
                chiral_sim = 1 - abs(len(target_chiral) - len(node_chiral)) / max(len(target_chiral), len(node_chiral)) if target_chiral or node_chiral else 1
                
                combined_sim = (tanimoto_sim + chiral_sim) / 2
                similarities.append((node, combined_sim, data))
        
        # 排序并返回前top_k个
        similarities.sort(key=lambda x: x[1], reverse=True)
        return similarities[:top_k]

    def filter_by_chiral_features(self, min_chiral_centers=1):
        """筛选具有特定手性特征的分子"""
        if not self.graph:
            if not self.load_graph():
                return []
        
        filtered = []
        for node, data in self.graph.nodes(data=True):
            if data['chiral_centers'] >= min_chiral_centers:
                filtered.append((node, data))
        
        return filtered

    def get_subgraph_by_similarity(self, target_smiles, threshold=0.7):
        """获取相似性子图"""
        if not self.graph:
            if not self.load_graph():
                return nx.Graph()
        
        similar_nodes = []
        similarities = self.find_similar_molecules(target_smiles, top_k=50)
        
        for node, sim, data in similarities:
            if sim >= threshold:
                similar_nodes.append(node)
        
        if similar_nodes:
            return self.graph.subgraph(similar_nodes)
        return nx.Graph()