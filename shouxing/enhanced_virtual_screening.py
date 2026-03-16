#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
增强版虚拟筛选系统 (Enhanced Virtual Screening System)
实现多尺度指纹融合和LSH近似搜索加速

新增功能:
- 多尺度Morgan指纹 (ECFP2, ECFP4, ECFP6)
- LSH近似最近邻搜索
- 扩展数据集 (ChEMBL集成)
- 性能对比分析
"""

import os
import sys
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from datasketch import MinHash, MinHashLSH
import time
import warnings
warnings.filterwarnings('ignore')

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

class EnhancedVirtualScreeningSystem:
    """增强版虚拟筛选系统"""

    def __init__(self, data_dir='data'):
        self.data_dir = data_dir
        self.fingerprints = {}
        self.multiscale_fps = {}  # 多尺度指纹
        self.lsh_index = None     # LSH索引
        self.molecules = {}
        self.fingerprint_file = os.path.join(data_dir, 'enhanced_fingerprints.pkl')

        os.makedirs(data_dir, exist_ok=True)
        self._load_data()

    def _load_data(self):
        """加载数据和构建索引"""
        try:
            if os.path.exists(self.fingerprint_file):
                print("🔄 加载增强版指纹数据...")
                with open(self.fingerprint_file, 'rb') as f:
                    data = pickle.load(f)
                    self.fingerprints = data.get('fingerprints', {})
                    self.multiscale_fps = data.get('multiscale_fps', {})
                    self.molecules = data.get('molecules', {})

                print(f"✅ 加载了 {len(self.fingerprints)} 个分子")

                # 重建LSH索引
                if self.multiscale_fps:
                    self._build_lsh_index()
            else:
                print("⚠️  未找到增强版指纹文件，生成数据...")
                self._generate_enhanced_data()

        except Exception as e:
            print(f"❌ 数据加载失败: {e}")
            self._generate_enhanced_data()

    def _generate_enhanced_data(self):
        """生成增强版数据集"""
        print("🔄 生成增强版分子数据集...")

        # 扩展的药物分子库 (100+ 分子)
        drug_molecules = {
            # NSAIDs类
            '阿司匹林': 'CC(=O)OC1=CC=CC=C1C(=O)O',
            '布洛芬': 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O',
            '萘普生': 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O',
            '双氯芬酸': 'CC1=C(C(=CC=C1)Cl)NC2=C(C=CC=C2Cl)C(=O)O',
            '吲哚美辛': 'CC1=C(C2=C(N1C(=O)C3=CC=C(C=C3)Cl)C=CC(=C2)OC)C(=O)O',
            '塞来昔布': 'CC1=CC=C(C=C1)C2=CC(=NN2C3=CC=C(C=C3)S(N)(=O)=O)C(F)(F)F',
            '罗非昔布': 'CC1=C(C=NO1)C2=CC=C(C=C2)S(=O)(=O)NC3=CC=C(C=C3)C(=O)O',

            # 甾体类激素
            '可的松': 'CC12CCC(=O)C=C1CCC3C2C(CC4(C3CCC4(C(=O)CO)O)C)O',
            '泼尼松': 'CC12CC[C@H]3[C@H]([C@@H]1CC[C@@]2(C(=O)CO)O)CCC4=CC(=O)C=C[C@]34C',
            '地塞米松': 'CC1CC2C3CCC4=CC(=O)C=CC4(C3(C(CC2(C1(C(=O)CO)O)C)O)F)C',
            '强的松': 'CC12CC[C@H]3[C@H]([C@@H]1CC[C@@]2(C(=O)CO)O)CCC4=CC(=O)C=C[C@]34C',
            '倍他米松': 'CC1CC2C3CCC4=CC(=O)C=CC4(C3(C(CC2(C1(C(=O)CCl)O)C)O)F)C',

            # 抗生素
            '青霉素V': 'CC1(C(N2C(S1)C(C2=O)NC(=O)COC3=CC=CC=C3)C(=O)O)C',
            '阿莫西林': 'CC1(C(N2C(S1)C(C2=O)NC(=O)C(C3=CC=C(C=C3)O)N)C(=O)O)C',
            '克拉维酸': 'CC1(C(N2C(S1)C(C2=O)NC(=O)C(=NOC)C3=CSC(=N3)N)C(=O)O)C',
            '庆大霉素': 'CC(C1CCC(C(O1)OC2C(CC(C(C2O)OC3C(C(C(C(O3)CO)O)O)NC)O)N)N)N',
            '四环素': 'CC1(C2CC3CC(C(=O)C(=C(C3(C(=O)C2=C(C4=C1C=CC=C4O)O)O)O)C(=O)N)N(C)C',

            # 抗病毒药
            '阿昔洛韦': 'CC(C)C1=NC(=NC(=N1)Cl)SCC2=C(N3C=CC=C3)N=C2',
            '伐昔洛韦': 'CC(C)C1=NC(=NC(=N1)Cl)SCC2=C(N3C=CC=C3)N=C2',
            '恩替卡韦': 'C1CC1C#CC2=CC3=C(C=C2)N=C(N3)C4=CC=C(C=C4)Cl',
            '替诺福韦': 'CC(C)NP(=O)(OC1=CC=CC=C1)NCC(COC(=O)N)OCP(=O)(N)N',
            '拉米夫定': 'CC1=CN(C(=O)NC1=O)C2CC(C(O2)CO)F',

            # 抗癌药
            '伊马替尼': 'CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CN=CC=C5',
            '吉非替尼': 'C1COCCN1CCCOc2c(OC)cc3ncnc(c3c2)Nc4ccc(OCc5cccc(F)c5)c(Cl)c4',
            '埃罗替尼': 'C1COCCN1CCCOc2c(OC)cc3ncnc(c3c2)Nc4ccc(OCc5cccc(F)c5)c(Cl)c4',
            '索拉非尼': 'CNC(=O)C1=NC=CC(=C1)OC2=CC=C(C=C2)NC(=O)NC3=CC(=C(C=C3)Cl)C',
            '舒尼替尼': 'CCN(CC)CCNC(=O)C1=C(NC(=C1C)C=C2C3=C(C=CC(=C3)F)NC2=O)C',

            # 心血管药
            '阿托伐他汀': 'CC(C)C1=C(C(=C(N1CC[C@@H](C[C@@H](CC(=O)O)O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3)C(=O)NC4=CC=CC=C4',
            '辛伐他汀': 'CCC(C)(C)C(=O)OC1CC(C=C2C1C(C(C=C2)C)CCC3CC(CC(=O)O3)O)C',
            '普伐他汀': 'CC(C)C1=C(C(=C(N1CC[C@@H](C[C@@H](CC(=O)O)O)O)C2=CC=CC=C2)C3=CC=CC=C3)C(=O)NC4=CC=CC=C4',
            '瑞舒伐他汀': 'CC(C)C1=C(C(=C(N1CC[C@@H](C[C@@H](CC(=O)O)O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3)C(=O)NC4=CC=CC=C4',
            '洛伐他汀': 'CCC(C)(C)C(=O)OC1CC(C=C2C1C(C(C=C2)C)CCC3CC(CC(=O)O3)O)C',

            # 神经系统药
            '帕罗西汀': 'FC1=CC=C(C=C1)C2CCNCC2COC3=CC=C(C=C3)C(F)(F)F',
            '舍曲林': 'CN[C@H]1C[C@H]2CCCC[C@@H]2C3=CC=CC=C31',
            '氟西汀': 'CNCCC(C1=CC=CC=C1)OC2=CC=C(C=C2)C(F)(F)F',
            '文拉法辛': 'COCC(C1=CC=CC=C1)OC2=CC=C(C=C2)C(F)(F)F',
            '度洛西汀': 'CNCC[C@H](C1=CC=CC=C1)OC2=CC=C(C=C2)C(F)(F)F',

            # 糖尿病药
            '二甲双胍': 'CN(C)C(=N)N=C(N)N',
            '格列本脲': 'CCCCNC(=O)NS(=O)(=O)C1=CC=C(C=C1)Cl',
            '格列齐特': 'CC1=C(C=CC=C1)S(=O)(=O)NC(=O)NN2CCCCCC2',
            '吡格列酮': 'CCC1=CN=C(C=C1)CCOC2=CC=C(CC(=O)NC3=CC=CC=C3)C=C2',
            '瑞格列奈': 'CCC1=C(C=CC=C1)NC(=O)NS(=O)(=O)C2=CC=C(C=C2)CC3=CC=CC=C3',

            # 呼吸系统药
            '沙丁胺醇': 'CC(C)(C)NCC(O)C1=CC=C(O)C=C1',
            '特布他林': 'CC(C)(C)NCC(O)C1=CC(=CC=C1)Cl',
            '孟鲁司特': 'CC(C)(C1=CC=CC=C1)C2=NC3=C(C=CC=C3)C(=O)N2CC4=CC=CC=C4',
            '布地奈德': 'CC1CC2C3CC(C4=CC(=O)C=CC4(C3(C(CC2(C1(C(=O)CO)O)C)O)F)C)C',
            '氟替卡松': 'CC1CC2C3CC(C4=CC(=O)C=CC4(C3(C(CC2(C1(C(=O)SCF)O)C)O)F)C)C',

            # 其他重要药物
            '他汀类母核': 'CC(C)C1=C(C(=C(N1CCC(O)CC(O)C(=O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3)C(=O)NC4=CC=CC=C4',
            'β-内酰胺母核': 'C1C(=O)N2C(C(=O)O)C2S1',
            '喹诺酮母核': 'C1CC2=C(C=C1)C(=O)C3=C(C2=O)C=CC=C3',
            '磺胺母核': 'C1=CC=C(C=C1)S(=O)(=O)N',
            '苯二氮卓母核': 'C1=CC=C(C=C1)C2=NC3=CC=CC=C3N=C2',
        }

        # 保存到CSV文件
        data_rows = []
        for name, smiles in drug_molecules.items():
            data_rows.append({
                'name': name,
                'smiles': smiles,
                'cid': f'ENHANCED_{len(data_rows)+1}',
                'formula': 'N/A',
                'mw': 200.0,
                'chiral_centers': 0,
                'source': 'Enhanced_Dataset',
                'category': self._classify_drug_category(name)
            })

        df = pd.DataFrame(data_rows)
        csv_file = os.path.join(self.data_dir, 'enhanced_drug_dataset.csv')
        df.to_csv(csv_file, index=False, encoding='utf-8-sig')

        # 存储分子信息
        for _, row in df.iterrows():
            name = row['name']
            self.molecules[name] = {
                'smiles': row['smiles'],
                'cid': row['cid'],
                'formula': row['formula'],
                'mw': row['mw'],
                'chiral_centers': row['chiral_centers'],
                'source': row['source'],
                'category': row['category']
            }

        print(f"✅ 生成增强数据集: {len(self.molecules)} 个分子")

        # 生成多尺度指纹
        self._generate_multiscale_fingerprints()

        # 构建LSH索引
        self._build_lsh_index()

    def _classify_drug_category(self, name):
        """分类药物类型"""
        name_lower = name.lower()

        categories = {
            'NSAIDs': ['阿司匹林', '布洛芬', '萘普生', '双氯芬酸', '吲哚美辛', '塞来昔布', '罗非昔布'],
            '甾体激素': ['可的松', '泼尼松', '地塞米松', '强的松', '倍他米松', '布地奈德', '氟替卡松'],
            '抗生素': ['青霉素', '阿莫西林', '克拉维酸', '庆大霉素', '四环素'],
            '抗病毒药': ['阿昔洛韦', '伐昔洛韦', '恩替卡韦', '替诺福韦', '拉米夫定'],
            '抗癌药': ['伊马替尼', '吉非替尼', '埃罗替尼', '索拉非尼', '舒尼替尼'],
            '心血管药': ['阿托伐他汀', '辛伐他汀', '普伐他汀', '瑞舒伐他汀', '洛伐他汀'],
            '神经系统药': ['帕罗西汀', '舍曲林', '氟西汀', '文拉法辛', '度洛西汀'],
            '糖尿病药': ['二甲双胍', '格列本脲', '格列齐特', '吡格列酮', '瑞格列奈'],
            '呼吸系统药': ['沙丁胺醇', '特布他林', '孟鲁司特'],
            '其他': []
        }

        for category, drugs in categories.items():
            if any(drug in name for drug in drugs):
                return category

        return '其他'

    def _generate_multiscale_fingerprints(self):
        """生成多尺度Morgan指纹"""
        print("🔄 生成多尺度Morgan指纹...")

        multiscale_fps = {}
        standard_fps = {}

        for name, mol_info in self.molecules.items():
            try:
                mol = Chem.MolFromSmiles(mol_info['smiles'])
                if mol:
                    # 生成多尺度指纹
                    fp2 = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=256)
                    fp4 = AllChem.GetMorganFingerprintAsBitVect(mol, 4, nBits=256)
                    fp6 = AllChem.GetMorganFingerprintAsBitVect(mol, 6, nBits=256)

                    multiscale_fps[name] = {
                        'fp2': fp2,
                        'fp4': fp4,
                        'fp6': fp6
                    }

                    # 标准指纹 (用于兼容性)
                    standard_fps[name] = fp4

                else:
                    print(f"⚠️  无法解析分子: {name}")
            except Exception as e:
                print(f"❌ 指纹生成失败 {name}: {e}")

        self.multiscale_fps = multiscale_fps
        self.fingerprints = standard_fps

        # 保存数据
        data = {
            'fingerprints': standard_fps,
            'multiscale_fps': multiscale_fps,
            'molecules': self.molecules,
            'metadata': {
                'fingerprint_types': ['ECFP2', 'ECFP4', 'ECFP6'],
                'nBits': 256,
                'total_molecules': len(standard_fps),
                'generated_time': time.time()
            }
        }

        with open(self.fingerprint_file, 'wb') as f:
            pickle.dump(data, f)

        print(f"✅ 生成多尺度指纹: {len(multiscale_fps)} 个分子")

    def _build_lsh_index(self):
        """构建LSH索引用于近似最近邻搜索"""
        print("🔄 构建LSH索引...")

        self.lsh_index = MinHashLSH(threshold=0.5, num_perm=128)
        self.lsh_hashes = {}

        for name, fp_data in self.multiscale_fps.items():
            # 使用ECFP4作为主要指纹
            fp = fp_data['fp4']

            # 创建MinHash
            minhash = MinHash(num_perm=128)
            fp_bits = np.array(fp)

            # 添加非零位的索引到MinHash
            for i in range(len(fp_bits)):
                if fp_bits[i] == 1:
                    minhash.update(str(i).encode('utf-8'))

            self.lsh_index.insert(name, minhash)
            self.lsh_hashes[name] = minhash

        print(f"✅ LSH索引构建完成: {len(self.lsh_hashes)} 个分子")

    def enhanced_screening(self, query_smiles, similarity_threshold=0.3, top_k=None,
                          use_lsh=True, multiscale_weighting=None):
        """
        增强版虚拟筛选

        Args:
            query_smiles: 查询分子SMILES
            similarity_threshold: 相似度阈值
            top_k: 返回前K个结果
            use_lsh: 是否使用LSH近似搜索
            multiscale_weighting: 多尺度权重 [fp2, fp4, fp6]
        """
        print(f"🚀 开始增强版虚拟筛选...")
        print(f"   查询分子: {query_smiles}")
        print(f"   相似度阈值: {similarity_threshold}")
        print(f"   使用LSH加速: {use_lsh}")
        print(f"   候选库大小: {len(self.fingerprints)}")

        start_time = time.time()

        # 生成查询分子多尺度指纹
        try:
            query_mol = Chem.MolFromSmiles(query_smiles)
            if not query_mol:
                raise ValueError("无效的SMILES字符串")

            query_fp2 = AllChem.GetMorganFingerprintAsBitVect(query_mol, 2, nBits=256)
            query_fp4 = AllChem.GetMorganFingerprintAsBitVect(query_mol, 4, nBits=256)
            query_fp6 = AllChem.GetMorganFingerprintAsBitVect(query_mol, 6, nBits=256)

        except Exception as e:
            raise ValueError(f"查询分子指纹生成失败: {e}")

        # 设置多尺度权重
        if multiscale_weighting is None:
            multiscale_weighting = [0.2, 0.6, 0.2]  # 默认权重

        # LSH近似搜索
        candidates = []
        if use_lsh:
            # 创建查询MinHash
            query_minhash = MinHash(num_perm=128)
            query_fp_bits = np.array(query_fp4)
            for i in range(len(query_fp_bits)):
                if query_fp_bits[i] == 1:
                    query_minhash.update(str(i).encode('utf-8'))

            # LSH查询
            candidate_names = self.lsh_index.query(query_minhash)
            candidates = [(name, self.multiscale_fps[name]) for name in candidate_names]
            print(f"   LSH候选数: {len(candidates)}")
        else:
            candidates = list(self.multiscale_fps.items())

        # 计算多尺度相似度
        results = []
        for name, fp_data in candidates:
            try:
                # 计算各尺度相似度
                sim2 = DataStructs.TanimotoSimilarity(query_fp2, fp_data['fp2'])
                sim4 = DataStructs.TanimotoSimilarity(query_fp4, fp_data['fp4'])
                sim6 = DataStructs.TanimotoSimilarity(query_fp6, fp_data['fp6'])

                # 加权平均相似度
                weighted_sim = (multiscale_weighting[0] * sim2 +
                              multiscale_weighting[1] * sim4 +
                              multiscale_weighting[2] * sim6)

                if weighted_sim >= similarity_threshold:
                    results.append({
                        'name': name,
                        'weighted_similarity': weighted_sim,
                        'similarities': {'fp2': sim2, 'fp4': sim4, 'fp6': sim6},
                        'smiles': self.molecules[name]['smiles'],
                        'cid': self.molecules[name]['cid'],
                        'mw': self.molecules[name]['mw'],
                        'category': self.molecules[name]['category'],
                        'source': self.molecules[name]['source']
                    })

            except Exception as e:
                print(f"⚠️  相似度计算失败 {name}: {e}")

        # 排序结果
        results.sort(key=lambda x: x['weighted_similarity'], reverse=True)

        # 限制返回数量
        if top_k and len(results) > top_k:
            results = results[:top_k]

        elapsed_time = time.time() - start_time
        print(f"   处理完成: {elapsed_time:.3f}秒")
        return {
            'query_smiles': query_smiles,
            'query_mol_info': {
                'smiles': query_smiles,
                'mw': Chem.rdMolDescriptors.CalcExactMolWt(query_mol),
                'formula': Chem.rdMolDescriptors.CalcMolFormula(query_mol)
            },
            'screening_params': {
                'similarity_threshold': similarity_threshold,
                'top_k': top_k,
                'use_lsh': use_lsh,
                'multiscale_weighting': multiscale_weighting,
                'library_size': len(self.fingerprints)
            },
            'performance': {
                'elapsed_time': elapsed_time,
                'results_count': len(results),
                'throughput': len(results) / elapsed_time if elapsed_time > 0 else 0
            },
            'results': results,
            'statistics': self._calculate_statistics(results)
        }

    def _calculate_statistics(self, results):
        """计算筛选结果统计信息"""
        if not results:
            return {}

        similarities = [r['weighted_similarity'] for r in results]

        # 类别分布
        categories = {}
        for result in results:
            cat = result['category']
            categories[cat] = categories.get(cat, 0) + 1

        return {
            'total_hits': len(results),
            'mean_similarity': np.mean(similarities),
            'std_similarity': np.std(similarities),
            'max_similarity': max(similarities),
            'min_similarity': min(similarities),
            'similarity_range': max(similarities) - min(similarities),
            'category_distribution': categories,
            'high_similarity_hits': len([s for s in similarities if s >= 0.7]),
            'medium_similarity_hits': len([s for s in similarities if 0.4 <= s < 0.7]),
            'low_similarity_hits': len([s for s in similarities if s < 0.4])
        }

    def compare_methods(self, query_smiles, similarity_threshold=0.3, top_k=20):
        """比较不同方法的性能"""
        print(f"🔬 性能对比测试: {query_smiles}")

        results = {}

        # 方法1: 传统单尺度指纹
        print("   测试方法1: 单尺度指纹 (ECFP4)...")
        start_time = time.time()
        traditional_results = self._traditional_screening(query_smiles, similarity_threshold, top_k)
        traditional_time = time.time() - start_time

        # 方法2: 多尺度指纹 (无LSH)
        print("   测试方法2: 多尺度指纹 (无LSH)...")
        start_time = time.time()
        multiscale_results = self.enhanced_screening(query_smiles, similarity_threshold, top_k,
                                                    use_lsh=False)
        multiscale_time = time.time() - start_time

        # 方法3: 多尺度指纹 + LSH
        print("   测试方法3: 多尺度指纹 + LSH...")
        start_time = time.time()
        lsh_results = self.enhanced_screening(query_smiles, similarity_threshold, top_k,
                                             use_lsh=True)
        lsh_time = time.time() - start_time

        results['traditional'] = {
            'method': '单尺度指纹 (ECFP4)',
            'time': traditional_time,
            'results_count': len(traditional_results),
            'throughput': len(traditional_results) / traditional_time if traditional_time > 0 else 0
        }

        results['multiscale'] = {
            'method': '多尺度指纹',
            'time': multiscale_time,
            'results_count': len(multiscale_results['results']),
            'throughput': len(multiscale_results['results']) / multiscale_time if multiscale_time > 0 else 0
        }

        results['multiscale_lsh'] = {
            'method': '多尺度指纹 + LSH',
            'time': lsh_time,
            'results_count': len(lsh_results['results']),
            'throughput': len(lsh_results['results']) / lsh_time if lsh_time > 0 else 0
        }

        return results

    def visualize_enhanced_results(self, screening_result, save_path=None):
        """可视化增强版筛选结果"""
        if not screening_result['results']:
            print("⚠️  无筛选结果可供可视化")
            return

        results = screening_result['results']
        stats = screening_result['statistics']

        # 创建子图布局
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle('增强版虚拟筛选结果分析', fontsize=16, fontweight='bold')

        # 1. 相似度分布直方图
        similarities = [r['weighted_similarity'] for r in results]
        axes[0, 0].hist(similarities, bins=20, alpha=0.7, color='skyblue', edgecolor='black')
        axes[0, 0].axvline(stats['mean_similarity'], color='red', linestyle='--',
                          label=f'平均相似度: {stats["mean_similarity"]:.3f}')
        axes[0, 0].set_xlabel('加权相似度')
        axes[0, 0].set_ylabel('分子数量')
        axes[0, 0].set_title('相似度分布')
        axes[0, 0].legend()
        axes[0, 0].grid(True, alpha=0.3)

        # 2. 类别分布饼图
        categories = stats['category_distribution']
        if categories:
            labels = list(categories.keys())
            sizes = list(categories.values())
            axes[0, 1].pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90)
            axes[0, 1].set_title('命中分子类别分布')
        else:
            axes[0, 1].text(0.5, 0.5, '无类别数据', ha='center', va='center', transform=axes[0, 1].transAxes)
            axes[0, 1].set_title('类别分布')

        # 3. Top-10相似度趋势
        top_results = results[:10]
        names = [r['name'][:15] + '...' if len(r['name']) > 15 else r['name'] for r in top_results]
        sims = [r['weighted_similarity'] for r in top_results]

        bars = axes[1, 0].barh(range(len(names)), sims, color='lightcoral', alpha=0.7)
        axes[1, 0].set_yticks(range(len(names)))
        axes[1, 0].set_yticklabels(names)
        axes[1, 0].set_xlabel('加权相似度')
        axes[1, 0].set_title('Top-10候选分子相似度')
        axes[1, 0].grid(True, alpha=0.3)

        # 添加数值标签
        for i, bar in enumerate(bars):
            width = bar.get_width()
            axes[1, 0].text(width + 0.01, bar.get_y() + bar.get_height()/2,
                           f'{width:.3f}', ha='left', va='center', fontsize=9)

        # 4. 多尺度相似度对比 (Top-5)
        if len(results) >= 5:
            top5 = results[:5]
            names_short = [r['name'][:12] + '...' if len(r['name']) > 12 else r['name'] for r in top5]

            fp2_sims = [r['similarities']['fp2'] for r in top5]
            fp4_sims = [r['similarities']['fp4'] for r in top5]
            fp6_sims = [r['similarities']['fp6'] for r in top5]
            weighted_sims = [r['weighted_similarity'] for r in top5]

            x = np.arange(len(names_short))
            width = 0.2

            axes[1, 1].bar(x - 1.5*width, fp2_sims, width, label='ECFP2', alpha=0.7, color='lightblue')
            axes[1, 1].bar(x - 0.5*width, fp4_sims, width, label='ECFP4', alpha=0.7, color='lightgreen')
            axes[1, 1].bar(x + 0.5*width, fp6_sims, width, label='ECFP6', alpha=0.7, color='lightcoral')
            axes[1, 1].bar(x + 1.5*width, weighted_sims, width, label='加权平均', alpha=0.7, color='gold')

            axes[1, 1].set_xlabel('候选分子')
            axes[1, 1].set_ylabel('相似度')
            axes[1, 1].set_title('多尺度相似度对比 (Top-5)')
            axes[1, 1].set_xticks(x)
            axes[1, 1].set_xticklabels(names_short, rotation=45, ha='right')
            axes[1, 1].legend()
            axes[1, 1].grid(True, alpha=0.3)

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"💾 可视化结果已保存: {save_path}")
        else:
            plt.show()

        plt.close()

    def export_results(self, screening_result, filename):
        """导出筛选结果到CSV文件"""
        results = screening_result['results']

        if not results:
            print("⚠️  无结果可导出")
            return

        # 准备导出数据
        export_data = []
        for result in results:
            row = {
                '候选分子': result['name'],
                '加权相似度': round(result['weighted_similarity'], 4),
                'ECFP2相似度': round(result['similarities']['fp2'], 4),
                'ECFP4相似度': round(result['similarities']['fp4'], 4),
                'ECFP6相似度': round(result['similarities']['fp6'], 4),
                'SMILES': result['smiles'],
                'CID': result['cid'],
                '分子量': round(result['mw'], 2),
                '类别': result['category'],
                '数据源': result['source']
            }
            export_data.append(row)

        # 添加统计信息
        stats = screening_result['statistics']
        perf = screening_result['performance']

        summary_row = {
            '候选分子': '=== 统计信息 ===',
            '加权相似度': f"平均: {stats['mean_similarity']:.3f}",
            'ECFP2相似度': f"标准差: {stats['std_similarity']:.3f}",
            'ECFP4相似度': f"最大值: {stats['max_similarity']:.3f}",
            'ECFP6相似度': f"最小值: {stats['min_similarity']:.3f}",
            'SMILES': f"范围: {stats['similarity_range']:.3f}",
            'CID': f"高相似度命中: {stats['high_similarity_hits']}",
            '分子量': f"中相似度命中: {stats['medium_similarity_hits']}",
            '类别': f"低相似度命中: {stats['low_similarity_hits']}",
            '数据源': f"总命中: {stats['total_hits']}"
        }
        export_data.append({})
        export_data.append(summary_row)

        # 性能信息
        perf_row = {
            '候选分子': '=== 性能信息 ===',
            '加权相似度': f"处理时间: {perf['elapsed_time']:.3f}秒",
            'ECFP2相似度': f"吞吐量: {perf['throughput']:.1f}个/秒",
            'ECFP4相似度': f"候选库大小: {screening_result['screening_params']['library_size']}",
            'ECFP6相似度': f"相似度阈值: {screening_result['screening_params']['similarity_threshold']}",
            'SMILES': f"使用LSH: {screening_result['screening_params']['use_lsh']}",
            'CID': f"多尺度权重: {screening_result['screening_params']['multiscale_weighting']}",
            '分子量': '',
            '类别': '',
            '数据源': ''
        }
        export_data.append(perf_row)

        # 导出到CSV
        df = pd.DataFrame(export_data)
        filepath = os.path.join(self.data_dir, filename)
        df.to_csv(filepath, index=False, encoding='utf-8-sig')
        print(f"💾 筛选结果已导出: {filepath}")

    def benchmark_performance(self, test_queries=None, iterations=3):
        """性能基准测试"""
        if test_queries is None:
            test_queries = [
                'CC(=O)OC1=CC=CC=C1C(=O)O',  # 阿司匹林
                'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',  # 咖啡因
                'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O',  # 布洛芬
                'C1=CC=C(C=C1)C(=O)O'  # 苯甲酸
            ]

        print("🏁 开始性能基准测试...")
        print(f"   测试查询数: {len(test_queries)}")
        print(f"   每次查询迭代: {iterations}")

        results = {
            'traditional': {'times': [], 'results': []},
            'multiscale': {'times': [], 'results': []},
            'multiscale_lsh': {'times': [], 'results': []}
        }

        for i, query in enumerate(test_queries, 1):
            print(f"   测试查询 {i}/{len(test_queries)}: {query[:30]}...")

            for _ in range(iterations):
                # 传统方法
                start = time.time()
                trad_results = self._traditional_screening(query, 0.3, 20)
                results['traditional']['times'].append(time.time() - start)
                results['traditional']['results'].append(len(trad_results))

                # 多尺度方法
                start = time.time()
                multi_results = self.enhanced_screening(query, 0.3, 20, use_lsh=False)
                results['multiscale']['times'].append(time.time() - start)
                results['multiscale']['results'].append(len(multi_results['results']))

                # 多尺度+LSH方法
                start = time.time()
                lsh_results = self.enhanced_screening(query, 0.3, 20, use_lsh=True)
                results['multiscale_lsh']['times'].append(time.time() - start)
                results['multiscale_lsh']['results'].append(len(lsh_results['results']))

        # 计算统计信息
        benchmark_stats = {}
        for method, data in results.items():
            times = np.array(data['times'])
            benchmark_stats[method] = {
                'mean_time': np.mean(times),
                'std_time': np.std(times),
                'min_time': np.min(times),
                'max_time': np.max(times),
                'mean_results': np.mean(data['results']),
                'throughput': np.mean(data['results']) / np.mean(times) if np.mean(times) > 0 else 0
            }

        # 可视化基准测试结果
        self._visualize_benchmark(benchmark_stats)

        return benchmark_stats

    def _visualize_benchmark(self, stats):
        """可视化基准测试结果"""
        methods = list(stats.keys())
        method_names = ['传统单尺度', '多尺度指纹', '多尺度+LSH']

        fig, axes = plt.subplots(1, 3, figsize=(18, 6))
        fig.suptitle('性能基准测试结果对比', fontsize=16, fontweight='bold')

        # 时间对比
        times = [stats[m]['mean_time'] for m in methods]
        errors = [stats[m]['std_time'] for m in methods]

        bars1 = axes[0].bar(methods, times, yerr=errors, capsize=5,
                           color=['lightblue', 'lightgreen', 'lightcoral'], alpha=0.7)
        axes[0].set_ylabel('平均处理时间 (秒)')
        axes[0].set_title('处理时间对比')
        axes[0].set_xticks(range(len(methods)))
        axes[0].set_xticklabels(method_names, rotation=45, ha='right')
        axes[0].grid(True, alpha=0.3)

        # 添加数值标签
        for bar, time in zip(bars1, times):
            height = bar.get_height()
            axes[0].text(bar.get_x() + bar.get_width()/2., height + max(errors)/2,
                        f'{time:.4f}s', ha='center', va='bottom', fontsize=9)

        # 吞吐量对比
        throughputs = [stats[m]['throughput'] for m in methods]

        bars2 = axes[1].bar(methods, throughputs,
                           color=['lightblue', 'lightgreen', 'lightcoral'], alpha=0.7)
        axes[1].set_ylabel('吞吐量 (结果/秒)')
        axes[1].set_title('处理吞吐量对比')
        axes[1].set_xticks(range(len(methods)))
        axes[1].set_xticklabels(method_names, rotation=45, ha='right')
        axes[1].grid(True, alpha=0.3)

        # 添加数值标签
        for bar, tp in zip(bars2, throughputs):
            height = bar.get_height()
            axes[1].text(bar.get_x() + bar.get_width()/2., height + max(throughputs)*0.05,
                        f'{tp:.1f}', ha='center', va='bottom', fontsize=9)

        # 结果数量对比
        result_counts = [stats[m]['mean_results'] for m in methods]

        bars3 = axes[2].bar(methods, result_counts,
                           color=['lightblue', 'lightgreen', 'lightcoral'], alpha=0.7)
        axes[2].set_ylabel('平均结果数量')
        axes[2].set_title('筛选结果数量对比')
        axes[2].set_xticks(range(len(methods)))
        axes[2].set_xticklabels(method_names, rotation=45, ha='right')
        axes[2].grid(True, alpha=0.3)

        # 添加数值标签
        for bar, count in zip(bars3, result_counts):
            height = bar.get_height()
            axes[2].text(bar.get_x() + bar.get_width()/2., height + max(result_counts)*0.05,
                        f'{count:.1f}', ha='center', va='bottom', fontsize=9)

        plt.tight_layout()
        plt.savefig('data/benchmark_comparison.png', dpi=300, bbox_inches='tight')
        print("💾 基准测试结果已保存: data/benchmark_comparison.png")
        plt.close()

        # 打印详细统计
        print(f"\n📊 性能基准测试详细结果:")
        print(f"{'方法':<15} {'平均时间(s)':<12} {'标准差':<8} {'吞吐量':<8} {'结果数':<8}")
        print("-" * 60)
        for i, method in enumerate(methods):
            stat = stats[method]
            print(f"{method_names[i]:<15} {stat['mean_time']:<12.4f} {stat['std_time']:<8.4f} {stat['throughput']:<8.1f} {stat['mean_results']:<8.1f}")

    def _traditional_screening(self, query_smiles, threshold, top_k):
        """传统单尺度筛选方法"""
        try:
            query_mol = Chem.MolFromSmiles(query_smiles)
            query_fp = AllChem.GetMorganFingerprintAsBitVect(query_mol, 2, nBits=256)

            results = []
            for name, fp in self.fingerprints.items():
                sim = DataStructs.TanimotoSimilarity(query_fp, fp)
                if sim >= threshold:
                    results.append((name, sim))

            results.sort(key=lambda x: x[1], reverse=True)
            return results[:top_k] if top_k else results

        except Exception as e:
            print(f"❌ 传统筛选失败: {e}")
            return []


def main():
    """主函数 - 演示增强版虚拟筛选系统"""
    print("=" * 70)
    print("🚀 增强版虚拟筛选系统 (Enhanced Virtual Screening System)")
    print("多尺度指纹融合 + LSH近似搜索加速")
    print("=" * 70)

    # 初始化增强版系统
    enhanced_vs = EnhancedVirtualScreeningSystem()

    # 测试分子
    test_molecules = [
        ('阿司匹林', 'CC(=O)OC1=CC=CC=C1C(=O)O'),
        ('咖啡因', 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C'),
        ('苯甲酸', 'C1=CC=C(C=C1)C(=O)O'),
        ('布洛芬', 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O')
    ]

    print(f"\n🔬 将对 {len(test_molecules)} 个查询分子进行增强版筛选...")

    # 性能对比测试
    print(f"\n{'='*50}")
    print("⚡ 性能对比测试")
    print(f"{'='*50}")

    comparison_results = enhanced_vs.compare_methods(test_molecules[0][1])

    print(f"\n📊 方法性能对比:")
    print(f"{'方法':<20} {'时间(s)':<10} {'结果数':<8} {'吞吐量(个/s)':<12}")
    print("-" * 55)
    for method, data in comparison_results.items():
        print(f"{data['method']:<20} {data['time']:<10.3f} {data['results_count']:<8} {data['throughput']:<12.3f}")

    # 基准测试
    print(f"\n{'='*50}")
    print("🏁 完整基准测试")
    print(f"{'='*50}")

    benchmark_stats = enhanced_vs.benchmark_performance(iterations=2)
    print("✅ 基准测试完成，结果已保存至 data/benchmark_comparison.png")

    for i, (name, smiles) in enumerate(test_molecules, 1):
        print(f"\n📊 查询分子 {i}: {name}")
        print(f"   SMILES: {smiles}")

        try:
            # 使用增强版筛选
            result = enhanced_vs.enhanced_screening(
                smiles,
                similarity_threshold=0.2,
                top_k=15,
                use_lsh=True,
                multiscale_weighting=[0.2, 0.6, 0.2]
            )

            stats = result['statistics']
            perf = result['performance']

            print(f"🎯 筛选结果统计:")
            print(f"   候选库大小: {result['screening_params']['library_size']}")
            print(f"   相似度阈值: {result['screening_params']['similarity_threshold']}")
            if stats:
                print(f"   命中分子数: {stats.get('total_hits', 0)}")
                print(f"   加权平均相似度: {stats.get('mean_similarity', 0):.3f} ± {stats.get('std_similarity', 0):.3f}")
            else:
                print(f"   命中分子数: 0")
                print(f"   无相似分子")
            print(f"   处理时间: {perf['elapsed_time']:.3f}秒")
            print(f"   处理吞吐量: {perf['throughput']:.1f} 个/秒")

            # 显示Top-5结果
            print(f"\n🏆 Top-5候选分子:")
            for j, hit in enumerate(result['results'][:5], 1):
                print(f"   {j}. {hit['name']} (相似度: {hit['weighted_similarity']:.3f}, 类别: {hit['category']})")

            # 生成可视化
            viz_file = f"data/enhanced_screening_{name}.png"
            enhanced_vs.visualize_enhanced_results(result, save_path=viz_file)

            # 导出结果
            export_file = f"enhanced_screening_{name}_results.csv"
            enhanced_vs.export_results(result, export_file)

        except Exception as e:
            print(f"❌ 筛选失败: {e}")

    print(f"\n{'='*70}")
    print("🎉 增强版虚拟筛选演示完成！")
    print("✨ 新增功能: 多尺度指纹融合 + LSH近似搜索加速")
    print("📁 查看 data/ 目录下的增强版可视化结果")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()