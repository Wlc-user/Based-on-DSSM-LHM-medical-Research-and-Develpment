import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
import pickle
import os
from sklearn.metrics.pairwise import cosine_similarity

class MoleculeFingerprintRetriever:
    def __init__(self, data_file='data/pubchem_bulk_data.csv', fingerprint_size=256):
        self.data_file = data_file
        self.fingerprint_size = fingerprint_size
        self.fingerprints = []
        self.smiles_list = []
        self.names_list = []
        self.df = None
        
    def load_data(self):
        """加载分子数据"""
        if not os.path.exists(self.data_file):
            print(f"数据文件 {self.data_file} 不存在，请先运行数据爬取")
            return False
            
        self.df = pd.read_csv(self.data_file)
        print(f"加载了 {len(self.df)} 个分子")
        return True
    
    def generate_fingerprints(self):
        """生成Morgan指纹"""
        if self.df is None:
            if not self.load_data():
                return False
                
        print("生成Morgan指纹...")
        self.fingerprints = []
        self.smiles_list = []
        self.names_list = []
        
        for idx, row in self.df.iterrows():
            smiles = row['smiles']
            name = row.get('name', f'Molecule_{idx}')
            
            if pd.isna(smiles):
                continue
                
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                continue
                
            # 生成Morgan指纹
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=self.fingerprint_size)
            self.fingerprints.append(fp)
            self.smiles_list.append(smiles)
            self.names_list.append(name)
            
            if (idx + 1) % 100 == 0:
                print(f"已处理 {idx + 1} 个分子")
        
        print(f"成功生成 {len(self.fingerprints)} 个指纹")
        return True
    
    def save_fingerprints(self, filename='data/molecule_fingerprints.pkl'):
        """保存指纹数据"""
        data = {
            'fingerprints': self.fingerprints,
            'smiles_list': self.smiles_list,
            'names_list': self.names_list,
            'fingerprint_size': self.fingerprint_size
        }
        
        os.makedirs(os.path.dirname(filename), exist_ok=True)
        with open(filename, 'wb') as f:
            pickle.dump(data, f)
        print(f"指纹数据已保存到 {filename}")
    
    def load_fingerprints(self, filename='data/molecule_fingerprints.pkl'):
        """加载指纹数据"""
        try:
            with open(filename, 'rb') as f:
                data = pickle.load(f)
            
            self.fingerprints = data['fingerprints']
            self.smiles_list = data['smiles_list']
            self.names_list = data['names_list']
            self.fingerprint_size = data['fingerprint_size']
            
            print(f"加载了 {len(self.fingerprints)} 个指纹")
            return True
        except FileNotFoundError:
            print(f"指纹文件 {filename} 不存在")
            return False
    
    def calculate_similarity_matrix(self):
        """计算相似度矩阵（可选，用于分析）"""
        if not self.fingerprints:
            print("没有指纹数据，请先生成或加载指纹")
            return None
            
        n = len(self.fingerprints)
        similarity_matrix = np.zeros((n, n))
        
        print("计算相似度矩阵...")
        for i in range(n):
            for j in range(i, n):
                sim = DataStructs.TanimotoSimilarity(self.fingerprints[i], self.fingerprints[j])
                similarity_matrix[i, j] = sim
                similarity_matrix[j, i] = sim
                
            if (i + 1) % 50 == 0:
                print(f"已计算 {i + 1}/{n} 个分子的相似度")
        
        return similarity_matrix
    
    def retrieve_similar_molecules(self, query_smiles, top_k=10):
        """检索相似分子（召回逻辑）"""
        if not self.fingerprints:
            print("没有指纹数据，请先生成或加载指纹")
            return []
        
        # 生成查询分子的指纹
        query_mol = Chem.MolFromSmiles(query_smiles)
        if query_mol is None:
            print("无效的查询SMILES")
            return []
            
        query_fp = AllChem.GetMorganFingerprintAsBitVect(query_mol, 2, nBits=self.fingerprint_size)
        
        # 计算与所有候选分子的相似度
        similarities = []
        for i, fp in enumerate(self.fingerprints):
            sim = DataStructs.TanimotoSimilarity(query_fp, fp)
            similarities.append((i, sim))
        
        # 按相似度排序，取Top K
        similarities.sort(key=lambda x: x[1], reverse=True)
        top_results = similarities[:top_k]
        
        # 返回结果
        results = []
        for idx, sim in top_results:
            results.append({
                'name': self.names_list[idx],
                'smiles': self.smiles_list[idx],
                'similarity': sim,
                'rank': len(results) + 1
            })
        
        return results
    
    def get_similarity_distribution(self, query_smiles):
        """获取相似度分布"""
        if not self.fingerprints:
            return None
            
        query_mol = Chem.MolFromSmiles(query_smiles)
        if query_mol is None:
            return None
            
        query_fp = AllChem.GetMorganFingerprintAsBitVect(query_mol, 2, nBits=self.fingerprint_size)
        
        similarities = []
        for fp in self.fingerprints:
            sim = DataStructs.TanimotoSimilarity(query_fp, fp)
            similarities.append(sim)
            
        return similarities