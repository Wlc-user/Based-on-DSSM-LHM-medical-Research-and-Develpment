#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
化学知识图谱演示系统
展示知识图谱的智能推理和关联推荐功能

演示内容:
1. 知识图谱构建过程
2. 分子上下文智能解释
3. 相关概念关联推荐
4. 知识图谱可视化
"""

import json
import matplotlib.pyplot as plt
from chemistry_knowledge_graph import ChemistryKnowledgeGraph
from chemistry_interpreter import ChemicalFormulaInterpreter

def demonstrate_knowledge_graph():
    """演示化学知识图谱功能"""
    print("🧠 化学知识图谱演示系统")
    print("=" * 50)

    # 1. 构建知识图谱
    print("🔄 步骤1: 构建知识图谱...")
    kg = ChemistryKnowledgeGraph()

    print(f"   📊 知识图谱统计:")
    print(f"      • 节点数量: {len(kg.graph.nodes)}")
    print(f"      • 边数量: {len(kg.graph.edges)}")

    # 统计不同类型节点
    node_types = {}
    for node, data in kg.graph.nodes(data=True):
        node_type = data.get('type', 'unknown')
        node_types[node_type] = node_types.get(node_type, 0) + 1

    print("      • 节点类型分布:")
    for node_type, count in node_types.items():
        print(f"        - {node_type}: {count}")

    # 2. 演示分子上下文解释
    print("\n🔍 步骤2: 分子上下文智能解释...")

    test_molecules = [
        ("阿司匹林", "CC(=O)OC1=CC=CC=C1C(=O)O", "经典的非甾体抗炎药"),
        ("咖啡因", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "中枢神经系统兴奋剂"),
        ("布洛芬", "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", "常用止痛药"),
        ("青霉素", "CC1(C(N2C(S1)C(C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C", "β-内酰胺类抗生素")
    ]

    results = []

    for name, smiles, description in test_molecules:
        print(f"\n   分析分子: {name} ({description})")

        context = kg.explain_molecule_context(smiles, name)
        if context and 'error' not in context:
            print(f"      化学式: {context['basic_info']['formula']}")
            print(f"      分子量: {context['basic_info']['molecular_weight']} g/mol")
            print(f"      官能团: {len(context['functional_groups'])} 个")

            # 显示官能团详情
            for fg in context['functional_groups'][:3]:  # 只显示前3个
                print(f"        • {fg['name']} ({fg['type']})")

            print(f"      智能推荐: {len(context['recommendations'])} 条")
            for rec in context['recommendations'][:2]:  # 只显示前2个推荐
                print(f"        • {rec['title']}: {rec['description']}")

            # 统计知识关联
            knowledge_links = len(context['knowledge_context'])
            print(f"      知识关联: {knowledge_links} 个概念")

            results.append({
                'name': name,
                'description': description,
                'context': context
            })
        else:
            print(f"      ❌ 分析失败: {context.get('error', '未知错误')}")

    # 3. 演示概念关联查询
    print("\n🔗 步骤3: 概念关联查询...")

    test_concepts = [
        ('element', 'C', '碳元素'),
        ('element', 'O', '氧元素'),
        ('fg', 'carboxylic_acid', '羧酸基'),
        ('fg', 'aromatic', '芳香环')
    ]

    for concept_type, concept_id, concept_name in test_concepts:
        print(f"\n   查询概念: {concept_name}")
        related = kg.get_related_concepts(concept_type, concept_id, max_depth=2)
        print(f"      相关概念: {len(related)} 个")

        for node, data in related[:5]:  # 只显示前5个
            node_type = data.get('type', 'unknown')
            node_name = data.get('name', node.split('_')[-1])
            print(f"        • {node_name} ({node_type})")

    # 4. 生成可视化
    print("\n📊 步骤4: 生成知识图谱可视化...")

    # 选择一个子集进行可视化（避免图太复杂）
    subset_nodes = []
    for node in kg.graph.nodes():
        if kg.graph.nodes[node].get('type') in ['element', 'functional_group', 'drug_category']:
            subset_nodes.append(node)
        if len(subset_nodes) >= 30:  # 限制节点数量
            break

    if subset_nodes:
        fig = kg.visualize_knowledge_graph(subset=subset_nodes, figsize=(14, 10))
        plt.savefig('chemistry_knowledge_graph_demo.png', dpi=150, bbox_inches='tight')
        plt.close()
        print("      ✅ 可视化图表已保存: chemistry_knowledge_graph_demo.png")

    # 5. 生成演示报告
    print("\n📝 步骤5: 生成演示报告...")

    report = generate_demo_report(results, kg)
    with open('KNOWLEDGE_GRAPH_DEMO_REPORT.md', 'w', encoding='utf-8') as f:
        f.write(report)

    print("      ✅ 演示报告已生成: KNOWLEDGE_GRAPH_DEMO_REPORT.md")

    return kg, results

def generate_demo_report(results, kg):
    """生成演示报告"""
    report = """# 化学知识图谱演示报告

## 概述

本报告展示了化学知识图谱系统的构建和应用，实现了智能的分子分析和概念关联功能。

## 知识图谱统计

- **节点数量**: {node_count}
- **边数量**: {edge_count}
- **节点类型分布**:
{node_types}

## 分子分析结果

以下是对测试分子的智能分析结果：

""".format(
        node_count=len(kg.graph.nodes),
        edge_count=len(kg.graph.edges),
        node_types="\n".join([f"  - {k}: {v}" for k, v in get_node_type_distribution(kg).items()])
    )

    for result in results:
        context = result['context']
        report += f"""
### {result['name']} - {result['description']}

- **SMILES**: `{context['basic_info']['smiles']}`
- **化学式**: {context['basic_info']['formula']}
- **分子量**: {context['basic_info']['molecular_weight']} g/mol

#### 化学式分析
{format_formula_analysis(context['formula_analysis'])}

#### 官能团识别
{format_functional_groups(context['functional_groups'])}

#### 智能推荐
{format_recommendations(context['recommendations'])}

#### 知识关联
发现了 {len(context['knowledge_context'])} 个相关概念的关联。

---
"""

    report += """
## 核心特性

### 1. 多层次知识表示
- **元素周期表知识**: 包含20种常见元素的性质和相互关系
- **分子分类体系**: 有机物/无机物、药物分类等
- **官能团关系网络**: 9种主要官能团及其反应关系
- **药物分类图谱**: 基于ATC标准的药物分类体系
- **性质关联网络**: 分子性质间的相关性和影响关系

### 2. 智能推理能力
- **上下文感知**: 根据分子结构提供相关的化学知识
- **关联推荐**: 发现分子间的潜在关系和反应可能性
- **性质预测**: 基于结构推断分子性质和行为

### 3. 可扩展架构
- **模块化设计**: 易于添加新的知识类型和关系
- **标准化接口**: 支持多种查询和推理操作
- **可视化支持**: 提供知识图谱的可视化展示

## 应用场景

1. **药物研发**: 快速识别分子结构特征和潜在作用机制
2. **化学教育**: 提供直观的分子结构解释和概念关联
3. **虚拟筛选**: 基于知识图谱的智能分子相似性搜索
4. **毒性预测**: 通过结构-性质关联进行安全性评估

## 技术实现

- **图数据库**: 使用NetworkX构建内存图数据库
- **知识表示**: 节点-边模型，支持多种关系类型
- **推理引擎**: 基于图遍历的关联查询和推荐算法
- **可视化**: Matplotlib和Seaborn进行图谱可视化

## 未来扩展方向

1. **知识库扩充**: 增加更多元素、官能团和分子数据
2. **关系推理**: 实现更复杂的逻辑推理和路径发现
3. **机器学习集成**: 结合ML模型进行性质预测和关联发现
4. **多语言支持**: 支持多种语言的化学术语和解释
5. **实时更新**: 支持知识图谱的动态更新和版本控制

---
*报告生成时间: 2026年3月16日*
"""

    return report

def get_node_type_distribution(kg):
    """获取节点类型分布"""
    node_types = {}
    for node, data in kg.graph.nodes(data=True):
        node_type = data.get('type', 'unknown')
        node_types[node_type] = node_types.get(node_type, 0) + 1
    return node_types

def format_formula_analysis(analysis):
    """格式化化学式分析结果"""
    if not analysis:
        return "暂无分析结果"

    elements = analysis.get('elements', {})
    element_str = ", ".join([f"{elem}: {data['count']}个" for elem, data in elements.items()])

    return f"""
- **元素组成**: {element_str}
- **总原子数**: {analysis.get('total_atoms', 0)}
- **解释**: {analysis.get('explanation', '暂无解释')}
"""

def format_functional_groups(groups):
    """格式化官能团信息"""
    if not groups:
        return "- **检测结果**: 未检测到明显官能团"

    fg_str = "\n".join([f"- **{fg['name']}** ({fg['type']}): {fg.get('description', '暂无描述')}" for fg in groups])
    return f"""
- **检测到 {len(groups)} 个官能团**:
{fg_str}
"""

def format_recommendations(recommendations):
    """格式化推荐信息"""
    if not recommendations:
        return "- **推荐**: 暂无特殊推荐"

    rec_str = "\n".join([f"- **{rec['title']}**: {rec['description']}" for rec in recommendations])
    return f"""
- **{len(recommendations)} 条智能推荐**:
{rec_str}
"""

if __name__ == "__main__":
    # 运行演示
    kg, results = demonstrate_knowledge_graph()

    print("\n🎉 演示完成!")
    print("📁 生成的文件:")
    print("   • chemistry_knowledge_graph_demo.png - 知识图谱可视化")
    print("   • KNOWLEDGE_GRAPH_DEMO_REPORT.md - 详细演示报告")
    print("   • chemistry_knowledge_graph.json - 知识图谱数据")

    print("\n💡 知识图谱现已集成到系统中，可以通过以下方式使用:")
    print("   • 分子搜索时会自动提供上下文解释")
    print("   • 相似性分析会考虑知识关联")
    print("   • Web界面支持知识图谱查询")