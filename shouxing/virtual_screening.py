#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
虚拟筛选系统 (Virtual Screening System)
基于分子指纹相似性的高通量化合物筛选平台

功能特性:
- 大规模分子库筛选 (>1000个化合物)
- 多查询并行处理
- 相似度阈值筛选
- 结果排序和统计
- 可视化分析报告
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
import warnings
warnings.filterwarnings('ignore')

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

class VirtualScreeningSystem:
    """虚拟筛选系统主类"""

    def __init__(self, data_dir='data'):
        self.data_dir = data_dir
        self.fingerprints = {}
        self.molecules = {}
        self.fingerprint_file = os.path.join(data_dir, 'molecule_fingerprints.pkl')
        self.bulk_data_file = os.path.join(data_dir, 'pubchem_bulk_data.csv')

        # 创建数据目录
        os.makedirs(data_dir, exist_ok=True)

        # 加载数据
        self._load_data()

    def _load_data(self):
        """加载分子数据和指纹"""
        try:
            # 加载指纹数据
            if os.path.exists(self.fingerprint_file):
                with open(self.fingerprint_file, 'rb') as f:
                    fingerprint_data = pickle.load(f)

                # 处理不同格式的指纹数据
                if isinstance(fingerprint_data, dict):
                    if 'fingerprints' in fingerprint_data:
                        fingerprints_data = fingerprint_data['fingerprints']
                        if isinstance(fingerprints_data, dict):
                            # 新格式: fingerprints是dict
                            self.fingerprints = fingerprints_data
                        elif isinstance(fingerprints_data, list):
                            # 列表格式: 需要与names_list对应
                            if 'names_list' in fingerprint_data:
                                names = fingerprint_data['names_list']
                                if len(fingerprints_data) == len(names):
                                    self.fingerprints = dict(zip(names, fingerprints_data))
                                else:
                                    print("⚠️  指纹数据与分子名称数量不匹配")
                                    self.fingerprints = {}
                            else:
                                print("⚠️  找不到names_list")
                                self.fingerprints = {}
                        else:
                            print("⚠️  未知的fingerprints数据格式")
                            self.fingerprints = {}
                    else:
                        # 旧格式: 直接是dict
                        self.fingerprints = fingerprint_data
                elif isinstance(fingerprint_data, list):
                    # 列表格式: 需要与分子名称对应
                    if os.path.exists(self.bulk_data_file):
                        df = pd.read_csv(self.bulk_data_file)
                        names = df['name'].tolist()
                        if len(fingerprint_data) == len(names):
                            self.fingerprints = dict(zip(names, fingerprint_data))
                        else:
                            print("⚠️  指纹数据与分子名称数量不匹配")
                            self.fingerprints = {}
                    else:
                        print("⚠️  找不到分子数据文件")
                        self.fingerprints = {}
                else:
                    print("⚠️  未知的指纹数据格式")
                    self.fingerprints = {}

                print(f"✅ 加载了 {len(self.fingerprints)} 个分子指纹")
            else:
                print("⚠️  未找到指纹文件，需要先生成指纹数据")
                self._generate_fingerprints()

            # 加载分子数据
            if os.path.exists(self.bulk_data_file):
                df = pd.read_csv(self.bulk_data_file)
                for _, row in df.iterrows():
                    name = row['name']
                    smiles = row['smiles']
                    self.molecules[name] = {
                        'smiles': smiles,
                        'cid': row.get('cid', 'N/A'),
                        'formula': row.get('formula', 'N/A'),
                        'mw': row.get('mw', 0),
                        'chiral_centers': row.get('chiral_centers', 0),
                        'source': row.get('source', 'Unknown')
                    }
                print(f"✅ 加载了 {len(self.molecules)} 个分子信息")

        except Exception as e:
            print(f"❌ 数据加载失败: {e}")
            import traceback
            traceback.print_exc()
            self._generate_sample_data()

    def _generate_sample_data(self):
        """生成示例数据集"""
        print("🔄 生成示例分子数据集...")

        # 扩展的药物分子库
        drug_molecules = {
            # NSAIDs类
            '阿司匹林': 'CC(=O)OC1=CC=CC=C1C(=O)O',
            '布洛芬': 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O',
            '萘普生': 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O',
            '双氯芬酸': 'CC1=C(C(=CC=C1)Cl)NC2=C(C=CC=C2Cl)C(=O)O',
            '吲哚美辛': 'CC1=C(C2=C(N1C(=O)C3=CC=C(C=C3)Cl)C=CC(=C2)OC)C(=O)O',

            # 甾体类
            '可的松': 'CC12CCC(=O)C=C1CCC3C2C(CC4(C3CCC4(C(=O)CO)O)C)O',
            '泼尼松': 'CC12CC[C@H]3[C@H]([C@@H]1CC[C@@]2(C(=O)CO)O)CCC4=CC(=O)C=C[C@]34C',

            # 嘌呤类
            '咖啡因': 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
            '茶碱': 'CN1C2=C(C(=O)N(C1=O)C)NC=N2',
            '可可碱': 'CN1C2=C(C(=O)N(C1=O)C)NC=N2',
            '黄嘌呤': 'C1=NC2=C(N1)C(=O)NC=N2',
            '腺嘌呤': 'C1=NC2=C(N1)C(=O)NC=N2',

            # 芳香羧酸类
            '苯甲酸': 'C1=CC=C(C=C1)C(=O)O',
            '邻苯二甲酸': 'C1=CC=C(C(=C1)C(=O)O)C(=O)O',
            '对苯二甲酸': 'C1=CC(=CC=C1C(=O)O)C(=O)O',
            '水杨酸': 'C1=CC=C(C(=C1)C(=O)O)O',
            '肉桂酸': 'C1=CC=C(C=C1)C=CC(=O)O',

            # 其他药物
            '对乙酰氨基酚': 'CC(=O)NC1=CC=C(C=C1)O',
            '苯磺酸': 'C1=CC=C(C=C1)S(=O)(=O)O',
            '氯苯': 'C1=CC=C(C=C1)Cl',
            '硝基苯': 'C1=CC=C(C=C1)[N+](=O)[O-]',
            '喹啉': 'C1=CC=C2C(=C1)C=CC=N2',
            '异喹啉': 'C1=CC=C2C(=C1)C=NC=C2',
            '吲哚': 'C1=CC=C2C(=C1)C=CN2',
            '咔唑': 'C1=CC=C2C(=C1)C3=C(N2)C=CN=C3',

            # 生物碱
            '吗啡': 'CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O',
            '可卡因': 'CN1C2CCC1C(C(OC(=O)C3=CC=CC=C3)C2OC(=O)C)OC(=O)C',
            '海洛因': 'CN1CCC23C4C1CC5=C2C(=C(C=C5)OC(=O)C)OC3C(C=C4)OC(=O)C',

            # 其他化合物
            '葡萄糖': 'C(C1C(C(C(C(O1)O)O)O)O)O',
            '乙酸': 'CC(=O)O',
            '丙酮': 'CC(=O)C',
            '甲醇': 'CO',
            '乙醇': 'CCO',
            '苯': 'C1=CC=CC=C1',
            '甲苯': 'CC1=CC=CC=C1',
            '苯酚': 'C1=CC=C(C=C1)O',
            '苯胺': 'C1=CC=C(C=C1)N',
            '吡啶': 'C1=CC=NC=C1',
            '噻吩': 'C1=CSC=C1',
            '呋喃': 'C1=COC=C1',
            '四氢呋喃': 'C1CCOC1',
            '二氧六环': 'C1COCCO1'
        }

        # 保存到CSV文件
        data_rows = []
        for name, smiles in drug_molecules.items():
            data_rows.append({
                'name': name,
                'smiles': smiles,
                'cid': f'SAMPLE_{len(data_rows)+1}',
                'formula': 'N/A',
                'mw': 100.0,
                'chiral_centers': 0,
                'source': 'Sample'
            })

        df = pd.DataFrame(data_rows)
        df.to_csv(self.bulk_data_file, index=False)

        # 存储分子信息
        for _, row in df.iterrows():
            name = row['name']
            self.molecules[name] = {
                'smiles': row['smiles'],
                'cid': row['cid'],
                'formula': row['formula'],
                'mw': row['mw'],
                'chiral_centers': row['chiral_centers'],
                'source': row['source']
            }

        print(f"✅ 生成示例数据集: {len(self.molecules)} 个分子")

        # 生成指纹
        self._generate_fingerprints()

    def _generate_fingerprints(self):
        """生成分子指纹"""
        print("🔄 生成Morgan指纹...")

        fingerprints = {}
        for name, mol_info in self.molecules.items():
            try:
                mol = Chem.MolFromSmiles(mol_info['smiles'])
                if mol:
                    # 生成256维Morgan指纹
                    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=256)
                    fingerprints[name] = fp
                else:
                    print(f"⚠️  无法解析分子: {name}")
            except Exception as e:
                print(f"❌ 指纹生成失败 {name}: {e}")

        # 保存指纹数据
        fingerprint_data = {
            'fingerprints': fingerprints,
            'metadata': {
                'fingerprint_type': 'Morgan',
                'radius': 2,
                'nBits': 256,
                'total_molecules': len(fingerprints)
            }
        }

        with open(self.fingerprint_file, 'wb') as f:
            pickle.dump(fingerprint_data, f)

        self.fingerprints = fingerprints
        print(f"✅ 成功生成 {len(fingerprints)} 个分子指纹")

    def screen_molecules(self, query_smiles, similarity_threshold=0.3, top_k=None):
        """
        虚拟筛选主函数

        Args:
            query_smiles: 查询分子的SMILES字符串
            similarity_threshold: 相似度阈值
            top_k: 返回前K个结果，None表示返回所有符合阈值的结果

        Returns:
            筛选结果字典
        """
        print(f"🔍 开始虚拟筛选...")
        print(f"   查询分子: {query_smiles}")
        print(f"   相似度阈值: {similarity_threshold}")
        print(f"   候选库大小: {len(self.fingerprints)}")

        # 生成查询分子指纹
        try:
            query_mol = Chem.MolFromSmiles(query_smiles)
            if not query_mol:
                raise ValueError("无效的SMILES字符串")

            query_fp = AllChem.GetMorganFingerprintAsBitVect(query_mol, 2, nBits=256)
        except Exception as e:
            raise ValueError(f"查询分子指纹生成失败: {e}")

        # 计算相似度
        results = []
        for name, fp in self.fingerprints.items():
            try:
                similarity = DataStructs.TanimotoSimilarity(query_fp, fp)
                if similarity >= similarity_threshold:
                    results.append({
                        'name': name,
                        'similarity': similarity,
                        'smiles': self.molecules[name]['smiles'],
                        'cid': self.molecules[name]['cid'],
                        'mw': self.molecules[name]['mw'],
                        'source': self.molecules[name]['source']
                    })
            except Exception as e:
                print(f"⚠️  相似度计算失败 {name}: {e}")

        # 排序结果
        results.sort(key=lambda x: x['similarity'], reverse=True)

        # 限制返回数量
        if top_k and len(results) > top_k:
            results = results[:top_k]

        print(f"✅ 筛选完成，发现 {len(results)} 个候选分子")

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
                'library_size': len(self.fingerprints)
            },
            'results': results,
            'statistics': self._calculate_statistics(results)
        }

    def _calculate_statistics(self, results):
        """计算筛选结果统计信息"""
        if not results:
            return {}

        similarities = [r['similarity'] for r in results]

        return {
            'total_hits': len(results),
            'mean_similarity': np.mean(similarities),
            'std_similarity': np.std(similarities),
            'max_similarity': max(similarities),
            'min_similarity': min(similarities),
            'similarity_range': max(similarities) - min(similarities),
            'high_similarity_hits': len([s for s in similarities if s >= 0.7]),
            'medium_similarity_hits': len([s for s in similarities if 0.4 <= s < 0.7]),
            'low_similarity_hits': len([s for s in similarities if s < 0.4])
        }

    def visualize_screening_results(self, screening_results, save_path=None):
        """可视化筛选结果"""
        results = screening_results['results']
        if not results:
            print("⚠️  没有筛选结果可供可视化")
            return

        # 创建子图
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
        fig.suptitle(f'虚拟筛选结果分析 - 查询分子: {screening_results["query_smiles"]}',
                    fontsize=16, fontweight='bold')

        # 1. Top-10相似度条形图
        top_10 = results[:10]
        names = [r['name'] for r in top_10]
        similarities = [r['similarity'] for r in top_10]

        bars = ax1.barh(names, similarities, color='skyblue')
        ax1.set_xlabel('Tanimoto相似度')
        ax1.set_title('Top-10候选分子相似度')
        ax1.grid(True, alpha=0.3)

        # 添加数值标签
        for bar, sim in zip(bars, similarities):
            ax1.text(bar.get_width() + 0.01, bar.get_y() + bar.get_height()/2,
                    '.3f', ha='left', va='center')

        # 2. 相似度分布直方图
        all_similarities = [r['similarity'] for r in results]
        ax2.hist(all_similarities, bins=20, alpha=0.7, color='lightcoral', edgecolor='black')
        ax2.set_xlabel('Tanimoto相似度')
        ax2.set_ylabel('分子数量')
        ax2.set_title('相似度分布')
        ax2.grid(True, alpha=0.3)

        # 3. 累积相似度覆盖
        sorted_sims = sorted(all_similarities, reverse=True)
        cumulative = np.cumsum(sorted_sims) / np.sum(sorted_sims) * 100

        ax3.plot(range(1, len(cumulative) + 1), cumulative,
                'o-', color='darkgreen', linewidth=2, markersize=4)
        ax3.set_xlabel('Top-K分子数量')
        ax3.set_ylabel('累积相似度覆盖率 (%)')
        ax3.set_title('累积相似度覆盖分析')
        ax3.grid(True, alpha=0.3)
        ax3.set_xlim(1, len(cumulative))

        # 4. 分子类别分布
        categories = {}
        for result in results:
            category = self._classify_molecule(result['name'])
            categories[category] = categories.get(category, 0) + 1

        if categories:
            labels = list(categories.keys())
            sizes = list(categories.values())

            ax4.pie(sizes, labels=labels, autopct='%1.1f%%',
                   startangle=90, colors=plt.cm.Set3.colors)
            ax4.set_title('候选分子类别分布')

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"✅ 可视化结果已保存: {save_path}")
        else:
            plt.show()

    def _classify_molecule(self, name):
        """简单分类分子"""
        name_lower = name.lower()

        if any(keyword in name_lower for keyword in ['阿司匹林', '布洛芬', '萘普生', '双氯芬酸', '吲哚美辛']):
            return 'NSAIDs'
        elif any(keyword in name_lower for keyword in ['咖啡因', '茶碱', '可可碱', '黄嘌呤', '腺嘌呤']):
            return '嘌呤类'
        elif any(keyword in name_lower for keyword in ['苯甲酸', '邻苯二甲酸', '对苯二甲酸', '水杨酸', '肉桂酸']):
            return '芳香羧酸'
        elif any(keyword in name_lower for keyword in ['可的松', '泼尼松']):
            return '甾体类'
        elif any(keyword in name_lower for keyword in ['吗啡', '可卡因', '海洛因']):
            return '生物碱'
        elif any(keyword in name_lower for keyword in ['对乙酰氨基酚']):
            return '解热镇痛'
        else:
            return '其他'

    def batch_screening(self, query_smiles_list, similarity_threshold=0.3, top_k=20):
        """批量虚拟筛选"""
        print(f"🔄 开始批量虚拟筛选，共 {len(query_smiles_list)} 个查询分子...")

        all_results = []
        for i, query_smiles in enumerate(query_smiles_list, 1):
            print(f"\n📊 处理查询分子 {i}/{len(query_smiles_list)}: {query_smiles}")
            try:
                result = self.screen_molecules(query_smiles, similarity_threshold, top_k)
                all_results.append(result)
            except Exception as e:
                print(f"❌ 查询失败 {query_smiles}: {e}")

        return all_results

    def export_results(self, screening_results, output_file='screening_results.csv'):
        """导出筛选结果到CSV文件"""
        results = screening_results['results']

        if not results:
            print("⚠️  没有结果可导出")
            return

        # 转换为DataFrame
        df = pd.DataFrame(results)

        # 添加查询信息
        df['query_smiles'] = screening_results['query_smiles']
        df['query_mw'] = screening_results['query_mol_info']['mw']
        df['query_formula'] = screening_results['query_mol_info']['formula']

        # 重新排列列
        columns = ['query_smiles', 'query_mw', 'query_formula', 'name', 'similarity',
                  'smiles', 'cid', 'mw', 'source']
        df = df[columns]

        # 保存文件
        output_path = os.path.join(self.data_dir, output_file)
        df.to_csv(output_path, index=False, encoding='utf-8-sig')
        print(f"✅ 筛选结果已导出: {output_path}")
        print(f"   共 {len(results)} 个候选分子")


def main():
    """主函数 - 演示虚拟筛选系统"""
    print("=" * 60)
    print("🧪 虚拟筛选系统 (Virtual Screening System)")
    print("基于分子指纹相似性的高通量化合物筛选")
    print("=" * 60)

    # 初始化系统
    vs_system = VirtualScreeningSystem()

    # 示例查询分子
    query_molecules = [
        ('阿司匹林', 'CC(=O)OC1=CC=CC=C1C(=O)O'),
        ('咖啡因', 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C'),
        ('苯甲酸', 'C1=CC=C(C=C1)C(=O)O'),
        ('布洛芬', 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O')
    ]

    print(f"\n🔍 将对 {len(query_molecules)} 个查询分子进行虚拟筛选...")

    # 批量筛选
    batch_results = vs_system.batch_screening(
        [smiles for _, smiles in query_molecules],
        similarity_threshold=0.2,
        top_k=15
    )

    # 处理和显示结果
    for i, (result, (name, smiles)) in enumerate(zip(batch_results, query_molecules)):
        print(f"\n{'='*50}")
        print(f"📊 查询分子 {i+1}: {name}")
        print(f"   SMILES: {smiles}")
        print(f"{'='*50}")

        stats = result['statistics']
        print(f"🎯 筛选结果统计:")
        print(f"   候选库大小: {result['screening_params']['library_size']}")
        print(f"   相似度阈值: {result['screening_params']['similarity_threshold']}")
        print(f"   命中分子数: {stats['total_hits']}")
        print(f"   平均相似度: {stats['mean_similarity']:.3f} ± {stats['std_similarity']:.3f}")
        print(f"   相似度范围: [{stats['min_similarity']:.3f}, {stats['max_similarity']:.3f}]")

        # 显示Top-5结果
        print(f"\n🏆 Top-5候选分子:")
        for j, hit in enumerate(result['results'][:5], 1):
            print(f"   {j}. {hit['name']} (相似度: {hit['similarity']:.3f})")

        # 生成可视化
        viz_file = f"data/virtual_screening_{name}.png"
        vs_system.visualize_screening_results(result, save_path=viz_file)

        # 导出结果
        export_file = f"virtual_screening_{name}_results.csv"
        vs_system.export_results(result, export_file)

    print(f"\n{'='*60}")
    print("🎉 虚拟筛选演示完成！")
    print("📁 查看 data/ 目录下的可视化结果和导出文件")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()