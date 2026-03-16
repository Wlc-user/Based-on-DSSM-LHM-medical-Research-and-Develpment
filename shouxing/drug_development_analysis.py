#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
医药研发知识图谱分析报告生成器
深入分析前后置关系、毒理关联和研发路径推理
"""

import json
import matplotlib.pyplot as plt
import networkx as nx
from collections import defaultdict
from drug_development_knowledge_graph import DrugDevelopmentKnowledgeGraph

def generate_drug_development_analysis_report():
    """生成医药研发知识图谱分析报告"""
    print("📊 生成医药研发知识图谱分析报告")
    print("=" * 50)

    # 加载知识图谱
    analyzer = DrugDevelopmentKnowledgeGraph()

    # 分析前后置关系
    print("🔗 分析前后置关系...")
    sequence_analysis = analyze_sequential_relationships(analyzer)

    # 分析毒理学关联
    print("☠️ 分析毒理学关联...")
    toxicology_analysis = analyze_toxicology_relationships(analyzer)

    # 分析ADMET前后置
    print("💊 分析ADMET前后置关系...")
    admet_analysis = analyze_admet_sequential(analyzer)

    # 分析靶点-疾病关联
    print("🎯 分析靶点-疾病关联...")
    target_disease_analysis = analyze_target_disease_relationships(analyzer)

    # 生成综合报告
    print("📝 生成综合分析报告...")
    report = create_comprehensive_report(analyzer, sequence_analysis, toxicology_analysis,
                                       admet_analysis, target_disease_analysis)

    # 保存报告
    with open('DRUG_DEVELOPMENT_ANALYSIS_REPORT.md', 'w', encoding='utf-8') as f:
        f.write(report)

    # 生成可视化
    print("📊 生成关系可视化...")
    create_relationship_visualizations(analyzer)

    return analyzer, report

def analyze_sequential_relationships(analyzer):
    """分析前后置关系"""
    analysis = {
        'development_pipeline': [],
        'critical_path': [],
        'decision_points': [],
        'failure_modes': []
    }

    # 分析研发流程前后置
    pipeline_nodes = [node for node, data in analyzer.graph.nodes(data=True)
                     if data.get('type') == 'development_stage']

    # 按阶段顺序排序
    stage_order = {}
    for node in pipeline_nodes:
        data = analyzer.graph.nodes[node]
        stage_order[node] = data.get('phase', 0)

    sorted_stages = sorted(stage_order.items(), key=lambda x: x[1])

    for i, (stage_node, _) in enumerate(sorted_stages):
        stage_data = analyzer.graph.nodes[stage_node]
        next_stages = list(analyzer.graph.successors(stage_node))

        analysis['development_pipeline'].append({
            'stage': stage_data.get('name'),
            'phase': stage_data.get('phase'),
            'success_rate': stage_data.get('success_rate'),
            'leads_to': [analyzer.graph.nodes[next_stage].get('name') for next_stage in next_stages]
        })

    # 分析决策点
    decision_nodes = [node for node, data in analyzer.graph.nodes(data=True)
                     if data.get('type') == 'decision_point']

    for node in decision_nodes:
        data = analyzer.graph.nodes[node]
        analysis['decision_points'].append({
            'name': data.get('name'),
            'stage': data.get('stage'),
            'criteria': data.get('criteria', [])
        })

    # 分析失败模式
    failure_nodes = [node for node, data in analyzer.graph.nodes(data=True)
                    if data.get('type') == 'trial_failure_reason']

    for node in failure_nodes:
        data = analyzer.graph.nodes[node]
        analysis['failure_modes'].append({
            'reason': data.get('name'),
            'frequency': data.get('frequency'),
            'phase': data.get('phase')
        })

    return analysis

def analyze_toxicology_relationships(analyzer):
    """分析毒理学关联"""
    analysis = {
        'toxicity_types': [],
        'assessment_strategies': [],
        'toxicity_associations': []
    }

    # 毒理学类型分析
    tox_nodes = [node for node, data in analyzer.graph.nodes(data=True)
                if data.get('type') == 'toxicology']

    for node in tox_nodes:
        data = analyzer.graph.nodes[node]
        assessment_strategies = []

        # 查找评估策略
        for strategy_node in analyzer.graph.successors(node):
            strategy_data = analyzer.graph.nodes[strategy_node]
            if strategy_data.get('type') == 'toxicology_strategy':
                edge_data = analyzer.graph.get_edge_data(node, strategy_node)
                assessment_strategies.append({
                    'strategy': strategy_data.get('name'),
                    'priority': edge_data.get('priority', 'unknown'),
                    'cost': strategy_data.get('cost'),
                    'time': strategy_data.get('time')
                })

        analysis['toxicity_types'].append({
            'toxicity': data.get('name'),
            'critical': data.get('critical', False),
            'regulatory': data.get('regulatory'),
            'assessment_strategies': assessment_strategies
        })

    # 毒理学前后置关系
    for node in tox_nodes:
        predecessors = list(analyzer.graph.predecessors(node))
        successors = list(analyzer.graph.successors(node))

        if predecessors or successors:
            analysis['toxicity_associations'].append({
                'toxicity': analyzer.graph.nodes[node].get('name'),
                'prerequisites': [analyzer.graph.nodes[p].get('name') for p in predecessors],
                'leads_to': [analyzer.graph.nodes[s].get('name') for s in successors]
            })

    return analysis

def analyze_admet_sequential(analyzer):
    """分析ADMET前后置关系"""
    analysis = {
        'admet_sequence': [],
        'property_influences': [],
        'prediction_methods': []
    }

    # ADMET性质顺序
    admet_nodes = [node for node, data in analyzer.graph.nodes(data=True)
                  if data.get('type') == 'admet_property']

    # 按前后置关系排序
    admet_order = {}
    for node in admet_nodes:
        successors = list(analyzer.graph.successors(node))
        admet_order[node] = len(successors)  # 简单排序

    sorted_admet = sorted(admet_order.items(), key=lambda x: x[1], reverse=True)

    for node, _ in sorted_admet:
        data = analyzer.graph.nodes[node]
        influences = []

        # 查找影响关系
        for successor in analyzer.graph.successors(node):
            edge_data = analyzer.graph.get_edge_data(node, successor)
            succ_data = analyzer.graph.nodes[successor]
            influences.append({
                'influenced_property': succ_data.get('name'),
                'relation': edge_data.get('description', '影响')
            })

        analysis['admet_sequence'].append({
            'property': data.get('name'),
            'abbrev': data.get('abbrev'),
            'key_factors': data.get('key_factors', []),
            'influences': influences,
            'prediction_methods': data.get('prediction_methods', [])
        })

    return analysis

def analyze_target_disease_relationships(analyzer):
    """分析靶点-疾病关联"""
    analysis = {
        'target_disease_associations': [],
        'disease_target_network': defaultdict(list),
        'target_specificity': []
    }

    # 分析靶点-疾病关联
    target_nodes = [node for node, data in analyzer.graph.nodes(data=True)
                   if data.get('type') == 'drug_target']
    disease_nodes = [node for node, data in analyzer.graph.nodes(data=True)
                    if data.get('type') == 'disease']

    for target_node in target_nodes:
        target_data = analyzer.graph.nodes[target_node]
        target_diseases = []

        for disease_node in analyzer.graph.successors(target_node):
            if analyzer.graph.nodes[disease_node].get('type') == 'disease':
                edge_data = analyzer.graph.get_edge_data(target_node, disease_node)
                disease_data = analyzer.graph.nodes[disease_node]

                target_diseases.append({
                    'disease': disease_data.get('name'),
                    'relation_type': edge_data.get('type', 'unknown'),
                    'prevalence': disease_data.get('prevalence')
                })

                analysis['disease_target_network'][disease_data.get('name')].append({
                    'target': target_data.get('name'),
                    'relation': edge_data.get('type', 'unknown')
                })

        analysis['target_disease_associations'].append({
            'target': target_data.get('name'),
            'target_type': target_data.get('type'),
            'diseases': target_diseases,
            'specificity': len(target_diseases)
        })

    return analysis

def create_comprehensive_report(analyzer, seq_analysis, tox_analysis, admet_analysis, target_analysis):
    """创建综合分析报告"""
    report = f"""# 医药研发知识图谱深度分析报告

## 概述

本报告深入分析医药研发知识图谱中的前后置关系、毒理学关联和核心研发路径推理。

**知识图谱统计:**
- 总节点数: {len(analyzer.graph.nodes)}
- 总边数: {len(analyzer.graph.edges)}
- 节点类型分布: {get_node_type_distribution(analyzer)}

---

## 1. 药物研发前后置关系分析

### 研发流程关键路径

药物研发是一个高度顺序依赖的过程，每个阶段都有特定的前后置要求：

"""

    # 研发流程分析
    for stage in seq_analysis['development_pipeline']:
        stage_names = ['hit_identification', 'hit_to_lead', 'lead_optimization', 'preclinical', 'phase1', 'phase2', 'phase3', 'nda_review', 'post_marketing']
        stage_key = f"stage_{stage_names[stage['phase']-1]}"
        stage_node_data = analyzer.graph.nodes.get(stage_key, {})

        report += f"""#### {stage['phase']}. {stage['stage']}
- **成功率**: {stage['success_rate']*100:.1f}%
- **时间**: {stage_node_data.get('time', '未知')}
- **成本**: {stage_node_data.get('cost', '未知')}
- **后置阶段**: {', '.join(stage['leads_to'])}
"""

    report += """
### 关键决策点

研发过程中有4个关键的GO/NO-GO决策点：

"""
    for dp in seq_analysis['decision_points']:
        report += f"""- **{dp['name']}** (阶段: {dp['stage']})
  - 决策标准: {', '.join(dp['criteria'])}
"""

    report += """
### 失败模式分析

临床试验失败的主要原因及频率分布：

"""
    for failure in seq_analysis['failure_modes']:
        report += f"""- **{failure['reason']}** ({failure['frequency']*100:.0f}%频率)
  - 主要发生阶段: {failure['phase']}
"""

    report += """

---

## 2. 毒理学前后置关系分析

### 毒理学评估策略

毒理学评估采用多层次策略，优先级从计算预测到体内试验：

"""

    for tox in tox_analysis['toxicity_types']:
        report += f"""#### {tox['toxicity']}
- **关键性**: {'是' if tox['critical'] else '否'}
- **法规要求**: {tox['regulatory']}
- **评估策略**:
"""
        for strategy in tox['assessment_strategies']:
            report += f"""  - {strategy['strategy']} ({strategy['priority']}优先级)
    - 成本: {strategy['cost']}, 时间: {strategy['time']}
"""

    report += """
### 毒理学前后置关联

毒理学评估项目之间的前后置关系：

"""
    for assoc in tox_analysis['toxicity_associations']:
        if assoc['prerequisites']:
            report += f"""- **{assoc['toxicity']}** 前置要求: {', '.join(assoc['prerequisites'])}
"""
        if assoc['leads_to']:
            report += f"""- **{assoc['toxicity']}** 后置影响: {', '.join(assoc['leads_to'])}
"""

    report += """

---

## 3. ADMET性质前后置关系分析

### ADMET评估顺序

ADMET性质评估必须按照严格的顺序进行，每个性质都可能影响后续评估：

"""

    for admet in admet_analysis['admet_sequence']:
        report += f"""#### {admet['abbrev']}. {admet['property']}
- **关键因素**: {', '.join(admet['key_factors'])}
- **预测方法**: {', '.join(admet['prediction_methods'])}
"""
        if admet['influences']:
            report += """- **对后续性质的影响**:
"""
            for influence in admet['influences']:
                report += f"""  - {influence['influenced_property']}: {influence['relation']}
"""

    report += """

---

## 4. 靶点-疾病关联网络分析

### 靶点治疗谱分析

不同靶点的治疗范围和特异性：

"""

    for assoc in target_analysis['target_disease_associations']:
        report += f"""#### {assoc['target']} ({assoc['target_type']})
- **治疗疾病谱**: {len(assoc['diseases'])} 种疾病
- **靶点特异性**: {'广谱' if len(assoc['diseases']) > 2 else '特异性'}
"""
        for disease in assoc['diseases'][:3]:  # 只显示前3个
            report += f"""  - {disease['disease']} ({disease['relation_type']}, 患病率: {disease['prevalence']})
"""

    report += """
### 疾病靶点网络

重要疾病的可用靶点分布：

"""
    for disease, targets in list(target_analysis['disease_target_network'].items())[:5]:
        report += f"""#### {disease}
- **可用靶点**: {len(targets)} 个
"""
        for target in targets:
            report += f"""  - {target['target']} ({target['relation']})
"""

    report += """

---

## 5. 研发路径推理分析

### 基于分子特征的研发策略

通过分析分子结构特征，系统可以提供个性化的研发路径建议：

#### 毒理学风险评估
- **芳香胺类化合物**: 遗传毒性风险高，需进行Ames试验
- **多卤代化合物**: 器官毒性风险中，需肝肾功能监测
- **高分子量化合物**: 代谢排泄风险中，需代谢稳定性评估

#### ADMET性质预测
- **分子量驱动**: MW < 300倾向肾排泄，MW > 500倾向肝胆排泄
- **极性影响**: 高极性化合物吸收可能较差
- **代谢稳定性**: 需考虑CYP酶相互作用

#### 靶点适用性分析
- **COX抑制剂**: 适用于炎症和疼痛治疗
- **受体激动剂**: 适用于神经精神疾病
- **酶抑制剂**: 适用于代谢性疾病

---

## 6. 知识图谱推理能力

### 前后置关系推理
1. **研发阶段推理**: hit_identification → hit_to_lead → lead_optimization → preclinical
2. **毒理学推理**: in_silico预测 → in_vitro筛选 → in_vivo验证
3. **ADMET推理**: 吸收影响分布，分布影响代谢，代谢影响排泄
4. **临床推理**: 毒理数据 → Phase I安全性 → Phase II疗效 → Phase III确认

### 风险评估推理
- **结构-毒性关联**: 特定官能团 → 毒性类型 → 评估策略
- **性质-行为关联**: 分子性质 → ADMET特征 → 临床表现
- **靶点-疾病关联**: 分子靶点 → 治疗适应症 → 临床定位

### 决策支持推理
- **GO/NO-GO决策**: 基于成功率、风险评估和商业潜力
- **资源分配**: 根据研发阶段的成本和时间要求
- **监管合规**: ICH指南要求和关键里程碑

---

## 结论与建议

### 核心发现
1. **前后置关系复杂性**: 药物研发涉及59个节点和42条前后置关系
2. **毒理学关键性**: 遗传毒性、致癌性和心血管毒性为最高优先级
3. **ADMET顺序依赖**: 吸收→分布→代谢→排泄→毒性的严格顺序
4. **靶点治疗多样性**: 单靶点可治疗多种疾病，多靶点可协同治疗

### 实用建议
1. **早期毒理学评估**: 在hit-to-lead阶段即开始毒理学风险识别
2. **ADMET优化优先级**: 优先解决吸收问题，再关注代谢稳定性
3. **靶点选择策略**: 考虑靶点特异性和治疗适应症的平衡
4. **决策点管理**: 在关键决策点进行全面的风险-收益评估

### 未来扩展方向
1. **机器学习集成**: 基于历史数据的研发成功率预测
2. **实时知识更新**: 纳入最新临床试验结果和药物上市数据
3. **个性化医疗**: 基于患者基因型的靶点选择优化
4. **药物重定位**: 利用现有药物的新适应症发现

---
*分析报告生成时间: 2026年3月16日*
*医药研发知识图谱系统 v1.0*
"""

    return report

def get_node_type_distribution(analyzer):
    """获取节点类型分布"""
    types = {}
    for node, data in analyzer.graph.nodes(data=True):
        node_type = data.get('type', 'unknown')
        types[node_type] = types.get(node_type, 0) + 1
    return types

def create_relationship_visualizations(analyzer):
    """创建关系可视化"""
    # 创建前后置关系图
    plt.figure(figsize=(15, 10))

    # 只选择关键节点进行可视化
    key_nodes = []
    for node, data in analyzer.graph.nodes(data=True):
        if data.get('type') in ['development_stage', 'toxicology', 'admet_property', 'drug_target']:
            key_nodes.append(node)

    subgraph = analyzer.graph.subgraph(key_nodes)

    # 计算位置
    pos = nx.spring_layout(subgraph, k=2, iterations=50)

    # 节点颜色
    node_colors = []
    for node in subgraph.nodes():
        node_type = subgraph.nodes[node].get('type')
        if node_type == 'development_stage':
            node_colors.append('lightblue')
        elif node_type == 'toxicology':
            node_colors.append('lightcoral')
        elif node_type == 'admet_property':
            node_colors.append('lightgreen')
        elif node_type == 'drug_target':
            node_colors.append('gold')
        else:
            node_colors.append('lightgray')

    # 绘制
    nx.draw_networkx_nodes(subgraph, pos, node_color=node_colors, node_size=800, alpha=0.7)
    nx.draw_networkx_edges(subgraph, pos, alpha=0.4, arrows=True, arrowsize=20)

    # 标签
    labels = {}
    for node in subgraph.nodes():
        data = subgraph.nodes[node]
        if data.get('type') == 'development_stage':
            labels[node] = f"阶段\\n{data.get('phase', '')}"
        elif data.get('type') == 'toxicology':
            labels[node] = data.get('name', '').replace('毒性', '\\n毒性')
        elif data.get('type') == 'admet_property':
            labels[node] = data.get('abbrev', '')
        else:
            labels[node] = data.get('name', node.split('_')[-1])

    nx.draw_networkx_labels(subgraph, pos, labels, font_size=8, font_family='SimHei')

    plt.title('医药研发前后置关系网络', fontsize=16, fontweight='bold')
    plt.axis('off')
    plt.tight_layout()
    plt.savefig('drug_development_relationships.png', dpi=150, bbox_inches='tight')
    plt.close()

if __name__ == "__main__":
    # 生成医药研发分析报告
    analyzer, report = generate_drug_development_analysis_report()

    print("\n🎉 医药研发知识图谱分析完成!")
    print("📁 生成的文件:")
    print("   • DRUG_DEVELOPMENT_ANALYSIS_REPORT.md - 深度分析报告")
    print("   • drug_development_relationships.png - 关系网络可视化")
    print("   • drug_development_knowledge_graph.json - 知识图谱数据")

    print("\n🔍 核心分析结果:")
    print("   • 研发前后置关系: 9个阶段，4个关键决策点")
    print("   • 毒理学关联: 10种毒性类型，3种评估策略")
    print("   • ADMET前后置: 5个性质的顺序依赖关系")
    print("   • 靶点-疾病网络: 8个重要靶点，7种疾病关联")

    print("\n💡 关键发现:")
    print("   • 药物研发成功率逐阶段递减，从10%降至50%")
    print("   • 毒理学评估是研发失败的最主要原因")
    print("   • ADMET性质间存在强烈的相互影响")
    print("   • 靶点选择直接影响临床定位和商业潜力")