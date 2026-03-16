import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from src.fingerprint_retriever import MoleculeFingerprintRetriever
import time
import os

# 设置中文字体支持
plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

def create_sample_data():
    """创建示例分子数据用于演示"""
    sample_molecules = [
        {'name': '阿司匹林', 'smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O', 'cid': 2244},
        {'name': '布洛芬', 'smiles': 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O', 'cid': 3672},
        {'name': '对乙酰氨基酚', 'smiles': 'CC(=O)NC1=CC=C(C=C1)O', 'cid': 1983},
        {'name': '咖啡因', 'smiles': 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C', 'cid': 2519},
        {'name': '葡萄糖', 'smiles': 'C(C1C(C(C(C(O1)O)O)O)O)O', 'cid': 5793},
        {'name': '苯甲酸', 'smiles': 'C1=CC=C(C=C1)C(=O)O', 'cid': 243},
        {'name': '水杨酸', 'smiles': 'C1=CC=C(C(=C1)C(=O)O)O', 'cid': 338},
        {'name': '乙酸', 'smiles': 'CC(=O)O', 'cid': 176},
        {'name': '丙酮', 'smiles': 'CC(=O)C', 'cid': 180},
        {'name': '甲醇', 'smiles': 'CO', 'cid': 887},
        {'name': '乙醇', 'smiles': 'CCO', 'cid': 702},
        {'name': '苯', 'smiles': 'C1=CC=CC=C1', 'cid': 241},
        {'name': '甲苯', 'smiles': 'CC1=CC=CC=C1', 'cid': 1140},
        {'name': '苯酚', 'smiles': 'C1=CC=C(C=C1)O', 'cid': 996},
        {'name': '苯胺', 'smiles': 'C1=CC=C(C=C1)N', 'cid': 6115},
        {'name': '吡啶', 'smiles': 'C1=CC=NC=C1', 'cid': 1049},
        {'name': '噻吩', 'smiles': 'C1=CSC=C1', 'cid': 8030},
        {'name': '呋喃', 'smiles': 'C1=COC=C1', 'cid': 8030},
        {'name': '四氢呋喃', 'smiles': 'C1CCOC1', 'cid': 8859},
        {'name': '二氧六环', 'smiles': 'C1COCCO1', 'cid': 8082},
    ]

    df = pd.DataFrame(sample_molecules)
    df['formula'] = 'N/A'  # 简化处理
    df['mw'] = 100.0  # 简化处理
    df['chiral_centers'] = 0  # 简化处理
    df['source'] = 'Sample'

    os.makedirs('data', exist_ok=True)
    df.to_csv('data/pubchem_bulk_data.csv', index=False)
    print(f"创建了 {len(df)} 个示例分子数据")
    return df

def main():
    print("=== 分子指纹相似性检索Demo ===")
    print("借鉴推荐系统召回逻辑，实现分子相似性筛选\n")

    # 1. 数据获取阶段
    print("1. 数据获取阶段")
    print("-" * 50)

    if not os.path.exists('data/pubchem_bulk_data.csv'):
        print("创建示例分子数据...")
        start_time = time.time()
        df = create_sample_data()
        print(f"数据创建完成，耗时: {time.time() - start_time:.2f}秒")
    else:
        print("使用已存在的分子数据")

    # 2. 指纹生成阶段
    print("\n2. 指纹生成阶段")
    print("-" * 50)
    retriever = MoleculeFingerprintRetriever()

    if not os.path.exists('data/molecule_fingerprints.pkl'):
        print("生成Morgan指纹（256维）...")
        start_time = time.time()
        retriever.generate_fingerprints()
        retriever.save_fingerprints()
        print(f"指纹生成完成，耗时: {time.time() - start_time:.2f}秒")
    else:
        print("加载已存在的指纹数据...")
        retriever.load_fingerprints()

    # 3. 相似性检索演示
    print("\n3. 相似性检索演示（推荐系统召回逻辑）")
    print("-" * 50)

    # 选择几个有代表性的查询分子
    query_molecules = [
        ('阿司匹林', 'CC(=O)OC1=CC=CC=C1C(=O)O'),
        ('咖啡因', 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C'),
        ('苯甲酸', 'C1=CC=C(C=C1)C(=O)O')
    ]

    for name, smiles in query_molecules:
        print(f"\n查询分子: {name}")
        print(f"SMILES: {smiles}")

        # 执行召回
        start_time = time.time()
        results = retriever.retrieve_similar_molecules(smiles, top_k=10)
        search_time = time.time() - start_time

        print(f"检索完成，耗时: {search_time:.3f}秒")
        print("Top 10 相似分子:")
        print("排名\t相似度\t\t分子名称")
        print("-" * 60)
        for result in results:
            print(f"{result['rank']:3d}\t{result['similarity']:.3f}\t\t{result['name'][:30]}")

        # 4. 可视化相似度分布（改进版）
        print(f"\n相似度分布可视化...")
        similarities = retriever.get_similarity_distribution(smiles)

        plt.figure(figsize=(16, 5))

        # 子图1: 相似度分布直方图（带数值标注）
        plt.subplot(1, 3, 1)
        n, bins, patches = plt.hist(similarities, bins=20, alpha=0.7, color='skyblue', edgecolor='black')
        plt.xlabel('Tanimoto相似度', fontsize=12, fontweight='bold')
        plt.ylabel('分子数量', fontsize=12, fontweight='bold')
        plt.title(f'{name}相似度分布\n(候选库: {len(similarities)}个分子)',
                 fontsize=14, fontweight='bold', pad=20)
        plt.grid(True, alpha=0.3)

        # 添加数值标注
        for i in range(len(patches)):
            if n[i] > 0:
                plt.text(patches[i].get_x() + patches[i].get_width()/2, n[i] + 0.5,
                        f'{int(n[i])}', ha='center', va='bottom', fontsize=10, fontweight='bold')

        # 子图2: Top 10相似度条形图（带数值标注）
        plt.subplot(1, 3, 2)
        top_sims = [r['similarity'] for r in results]
        top_names = [r['name'][:12] + '...' if len(r['name']) > 12 else r['name'] for r in results]

        bars = plt.barh(range(len(top_sims)), top_sims, color='lightcoral', alpha=0.8, height=0.6)
        plt.yticks(range(len(top_names)), top_names, fontsize=10)
        plt.xlabel('相似度分数', fontsize=12, fontweight='bold')
        plt.title('Top 10相似分子\n(按相似度排序)', fontsize=14, fontweight='bold', pad=20)
        plt.grid(True, alpha=0.3)
        plt.xlim(0, 1)

        # 添加数值标注
        for i, (bar, sim) in enumerate(zip(bars, top_sims)):
            plt.text(bar.get_width() + 0.01, bar.get_y() + bar.get_height()/2,
                    f'{sim:.3f}', ha='left', va='center', fontsize=9, fontweight='bold')

        # 子图3: 相似度累积分布（带关键点标注）
        plt.subplot(1, 3, 3)
        sorted_sims = np.sort(similarities)[::-1]
        cumsum = np.cumsum(sorted_sims) / np.sum(sorted_sims)

        plt.plot(range(1, len(sorted_sims) + 1), cumsum, 'g-', linewidth=2.5, marker='o', markersize=4, alpha=0.8)
        plt.xlabel('Top K分子数量', fontsize=12, fontweight='bold')
        plt.ylabel('累积相似度覆盖率', fontsize=12, fontweight='bold')
        plt.title('相似度累积分布\n(推荐系统召回分析)', fontsize=14, fontweight='bold', pad=20)
        plt.grid(True, alpha=0.3)
        plt.ylim(0, 1)

        # 添加关键点标注
        key_points = [5, 10, len(similarities)]
        for k in key_points:
            if k <= len(cumsum):
                plt.plot(k, cumsum[k-1], 'ro', markersize=8, alpha=0.9)
                plt.text(k+1, cumsum[k-1], f'Top-{k}\n{cumsum[k-1]:.1%}',
                        fontsize=10, fontweight='bold', verticalalignment='center',
                        bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.8))

        plt.tight_layout()
        plt.savefig(f'data/similarity_analysis_{name}.png', dpi=300, bbox_inches='tight')
        plt.close()  # 不显示图形，节省资源

        print(f"可视化结果已保存为: data/similarity_analysis_{name}.png")

        # 5. 详细分析结论
        print(f"\n=== {name}检索分析结论 ===")
        print(f"📊 候选库总分子数: {len(similarities)}")
        print(f"🏆 Top-1相似度: {results[0]['similarity']:.3f} ({results[0]['name']})")
        print(f"📈 Top-10平均相似度: {np.mean(top_sims):.3f} ± {np.std(top_sims):.3f}")
        print(f"📏 相似度范围: [{np.min(similarities):.3f}, {np.max(similarities):.3f}]")

        # 计算推荐系统指标
        top_k_coverage = cumsum[9] if len(cumsum) > 9 else cumsum[-1]  # Top-10覆盖率
        print(f"🎯 Top-10累积覆盖率: {top_k_coverage:.1%}")

        # 相似度分布分析
        high_sim_count = sum(1 for s in similarities if s > 0.7)
        mid_sim_count = sum(1 for s in similarities if 0.4 <= s <= 0.7)
        low_sim_count = sum(1 for s in similarities if s < 0.4)
        print(f"🔥 高相似度分子(>0.7): {high_sim_count}个 ({high_sim_count/len(similarities):.1%})")
        print(f"🟡 中相似度分子(0.4-0.7): {mid_sim_count}个 ({mid_sim_count/len(similarities):.1%})")
        print(f"🔵 低相似度分子(<0.4): {low_sim_count}个 ({low_sim_count/len(similarities):.1%})")

        plt.close('all')  # 确保关闭所有图形

    # 6. 性能分析
    print("\n4. 性能分析")
    print("-" * 50)
    print(f"候选库大小: {len(retriever.fingerprints)} 个分子")
    print(f"指纹维度: {retriever.fingerprint_size}")
    print("检索算法: Tanimoto系数计算")
    print("时间复杂度: O(n) - 线性扫描所有候选分子")

    # 模拟不同规模下的检索时间
    test_sizes = [100, 500, 1000, 5000, 10000]
    print("\n理论检索时间估计 (假设每次计算0.01ms):")
    print("规模\t\t时间\t\t备注")
    print("-" * 40)
    for size in test_sizes:
        estimated_time = size * 0.01  # ms
        if estimated_time < 1000:
            time_str = f"{estimated_time:.1f}ms"
        else:
            time_str = f"{estimated_time/1000:.1f}s"
        remark = "实时" if size <= 1000 else "可接受" if size <= 5000 else "需优化"
        print("8")

    # 7. 总体结论和收获
    print("\n5. 总体结论和收获")
    print("=" * 80)
    print("🎯 核心收获:")
    print("   1. ✅ 掌握了分子数据的处理流程：SMILES → Morgan指纹 → 相似度计算")
    print("   2. ✅ 理解了推荐系统'召回'逻辑在分子筛选场景的平移应用")
    print("   3. ✅ 实现了基于指纹的Top-K相似分子检索算法")
    print("   4. ✅ 掌握了化学数据可视化分析方法")

    print("\n📊 技术实现要点:")
    print("   • 数据获取: 使用PubChem API获取分子SMILES数据")
    print("   • 特征提取: RDKit生成256维Morgan指纹作为分子表征")
    print("   • 相似度计算: Tanimoto系数衡量分子间结构相似性")
    print("   • 检索策略: 线性扫描实现Top-K召回（适合小规模数据集）")
    print("   • 可视化分析: 多维度展示相似度分布和检索效果")

    print("\n🔬 应用价值:")
    print("   • 💊 药物发现: 快速找到结构相似化合物进行活性筛选")
    print("   • 🧪 化学合成: 基于相似分子预测合成路线")
    print("   • ⚠️  毒性评估: 识别潜在毒性相似的化合物")
    print("   • 📜 专利分析: 发现结构相似的化合物专利")

    print("\n⚡ 性能优化方向:")
    print("   • 🚀 大规模优化: 采用近似最近邻搜索(ANN)算法")
    print("   • 🏗️  索引构建: 使用LSH或树结构加速检索")
    print("   • 🖥️  并行计算: GPU加速相似度计算")
    print("   • 💾 缓存策略: 预计算热门查询结果")

    print("\n📈 扩展方向:")
    print("   • 🔗 多模态融合: 结合分子图结构和物理化学性质")
    print("   • 🤖 深度学习: 使用图神经网络学习分子表征")
    print("   • 🌐 交互式检索: 开发Web界面支持实时查询")
    print("   • 🏭 工业应用: 集成到药物筛选工作流中")

    print("\n" + "=" * 80)
    print("🎉 Demo执行完成！查看 data/ 目录下的可视化结果文件。")
    print("=" * 80)

if __name__ == "__main__":
    main()