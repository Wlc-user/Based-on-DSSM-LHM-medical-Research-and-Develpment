#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
化学知识图谱系统
构建结构化的化学知识网络，实现智能推理和关联推荐

知识图谱包含:
- 元素周期表知识图谱
- 分子分类体系图谱
- 官能团关系图谱
- 药物分类图谱
- 性质关联图谱
"""

import json
import networkx as nx
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
from rdkit import Chem
from chemistry_interpreter import ChemicalFormulaInterpreter

class ChemistryKnowledgeGraph:
    """化学知识图谱"""

    def __init__(self):
        self.graph = nx.DiGraph()
        self.interpreter = ChemicalFormulaInterpreter()
        self._build_knowledge_graph()

    def _build_knowledge_graph(self):
        """构建知识图谱"""
        print("🔄 构建化学知识图谱...")

        # 1. 元素周期表知识
        self._build_element_graph()

        # 2. 分子分类体系
        self._build_molecule_classification_graph()

        # 3. 官能团关系网络
        self._build_functional_group_graph()

        # 4. 药物分类体系
        self._build_drug_classification_graph()

        # 5. 性质关联网络
        self._build_property_association_graph()

        print(f"✅ 知识图谱构建完成: {len(self.graph.nodes)} 个节点, {len(self.graph.edges)} 条边")

    def _build_element_graph(self):
        """构建元素周期表知识图谱"""
        print("   📊 构建元素周期表图谱...")

        # 元素基本信息
        elements_data = {
            'H': {'name': '氢', 'atomic_number': 1, 'period': 1, 'group': 1, 'electronegativity': 2.2, 'common_oxidation_states': [1, -1]},
            'C': {'name': '碳', 'atomic_number': 6, 'period': 2, 'group': 14, 'electronegativity': 2.55, 'common_oxidation_states': [4, 2, -4]},
            'N': {'name': '氮', 'atomic_number': 7, 'period': 2, 'group': 15, 'electronegativity': 3.04, 'common_oxidation_states': [5, 3, 2, 1, -3]},
            'O': {'name': '氧', 'atomic_number': 8, 'period': 2, 'group': 16, 'electronegativity': 3.44, 'common_oxidation_states': [2, -2]},
            'F': {'name': '氟', 'atomic_number': 9, 'period': 2, 'group': 17, 'electronegativity': 3.98, 'common_oxidation_states': [-1]},
            'Cl': {'name': '氯', 'atomic_number': 17, 'period': 3, 'group': 17, 'electronegativity': 3.16, 'common_oxidation_states': [7, 5, 3, 1, -1]},
            'Br': {'name': '溴', 'atomic_number': 35, 'period': 4, 'group': 17, 'electronegativity': 2.96, 'common_oxidation_states': [7, 5, 3, 1, -1]},
            'I': {'name': '碘', 'atomic_number': 53, 'period': 5, 'group': 17, 'electronegativity': 2.66, 'common_oxidation_states': [7, 5, 3, 1, -1]},
            'S': {'name': '硫', 'atomic_number': 16, 'period': 3, 'group': 16, 'electronegativity': 2.58, 'common_oxidation_states': [6, 4, 2, -2]},
            'P': {'name': '磷', 'atomic_number': 15, 'period': 3, 'group': 15, 'electronegativity': 2.19, 'common_oxidation_states': [5, 3, -3]},
            'Na': {'name': '钠', 'atomic_number': 11, 'period': 3, 'group': 1, 'electronegativity': 0.93, 'common_oxidation_states': [1]},
            'K': {'name': '钾', 'atomic_number': 19, 'period': 4, 'group': 1, 'electronegativity': 0.82, 'common_oxidation_states': [1]},
            'Ca': {'name': '钙', 'atomic_number': 20, 'period': 4, 'group': 2, 'electronegativity': 1.00, 'common_oxidation_states': [2]},
            'Mg': {'name': '镁', 'atomic_number': 12, 'period': 3, 'group': 2, 'electronegativity': 1.31, 'common_oxidation_states': [2]},
            'Fe': {'name': '铁', 'atomic_number': 26, 'period': 4, 'group': 8, 'electronegativity': 1.83, 'common_oxidation_states': [3, 2]},
            'Cu': {'name': '铜', 'atomic_number': 29, 'period': 4, 'group': 11, 'electronegativity': 1.90, 'common_oxidation_states': [2, 1]},
            'Zn': {'name': '锌', 'atomic_number': 30, 'period': 4, 'group': 12, 'electronegativity': 1.65, 'common_oxidation_states': [2]},
            'Al': {'name': '铝', 'atomic_number': 13, 'period': 3, 'group': 13, 'electronegativity': 1.61, 'common_oxidation_states': [3]},
            'Si': {'name': '硅', 'atomic_number': 14, 'period': 3, 'group': 14, 'electronegativity': 1.90, 'common_oxidation_states': [4]},
            'B': {'name': '硼', 'atomic_number': 5, 'period': 2, 'group': 13, 'electronegativity': 2.04, 'common_oxidation_states': [3]}
        }

        # 添加元素节点
        for symbol, data in elements_data.items():
            self.graph.add_node(f"element_{symbol}",
                              type="element",
                              symbol=symbol,
                              **data)

        # 添加元素间关系
        for symbol1 in elements_data:
            for symbol2 in elements_data:
                if symbol1 != symbol2:
                    # 电负性关系
                    en1 = elements_data[symbol1]['electronegativity']
                    en2 = elements_data[symbol2]['electronegativity']
                    if abs(en1 - en2) < 0.5:
                        self.graph.add_edge(f"element_{symbol1}", f"element_{symbol2}",
                                          relation="similar_electronegativity",
                                          weight=1-abs(en1-en2))

                    # 周期表位置关系
                    if elements_data[symbol1]['period'] == elements_data[symbol2]['period']:
                        self.graph.add_edge(f"element_{symbol1}", f"element_{symbol2}",
                                          relation="same_period")
                    if elements_data[symbol1]['group'] == elements_data[symbol2]['group']:
                        self.graph.add_edge(f"element_{symbol1}", f"element_{symbol2}",
                                          relation="same_group")

    def _build_molecule_classification_graph(self):
        """构建分子分类体系图谱"""
        print("   🧬 构建分子分类体系图谱...")

        # 分子大类
        molecule_classes = {
            'organic': {'name': '有机化合物', 'description': '含碳化合物'},
            'inorganic': {'name': '无机化合物', 'description': '不含碳或含碳但不含氢的化合物'},
            'drug': {'name': '药物分子', 'description': '具有药理活性的化合物'},
            'natural_product': {'name': '天然产物', 'description': '从生物体中提取的化合物'},
            'synthetic': {'name': '合成化合物', 'description': '人工合成的化合物'}
        }

        # 药物小类
        drug_subclasses = {
            'nsaids': {'name': '非甾体抗炎药', 'parent': 'drug', 'description': '如阿司匹林、布洛芬'},
            'antibiotics': {'name': '抗生素', 'parent': 'drug', 'description': '如青霉素、红霉素'},
            'antiviral': {'name': '抗病毒药', 'parent': 'drug', 'description': '如阿昔洛韦'},
            'antineoplastic': {'name': '抗肿瘤药', 'parent': 'drug', 'description': '如顺铂、多柔比星'},
            'cardiovascular': {'name': '心血管药', 'parent': 'drug', 'description': '如硝酸甘油'},
            'psychoactive': {'name': '精神药物', 'parent': 'drug', 'description': '如咖啡因、吗啡'}
        }

        # 添加分类节点
        for class_id, data in molecule_classes.items():
            self.graph.add_node(f"class_{class_id}",
                              type="molecule_class",
                              class_id=class_id,
                              **data)

        for subclass_id, data in drug_subclasses.items():
            self.graph.add_node(f"class_{subclass_id}",
                              type="molecule_subclass",
                              class_id=subclass_id,
                              **data)
            # 连接到父类
            self.graph.add_edge(f"class_{data['parent']}", f"class_{subclass_id}",
                              relation="has_subclass")

    def _build_functional_group_graph(self):
        """构建官能团关系图谱"""
        print("   ⚗️ 构建官能团关系图谱...")

        functional_groups = {
            'carboxylic_acid': {
                'name': '羧酸基', 'formula': 'COOH', 'reactivity': '酸性',
                'common_reactions': ['酯化', '成盐', '脱羧'], 'polarity': 'high'
            },
            'ester': {
                'name': '酯基', 'formula': 'COOR', 'reactivity': '中等',
                'common_reactions': ['水解', '皂化'], 'polarity': 'medium'
            },
            'amide': {
                'name': '酰胺基', 'formula': 'CONH2', 'reactivity': '低',
                'common_reactions': ['水解'], 'polarity': 'high'
            },
            'alcohol': {
                'name': '羟基', 'formula': 'OH', 'reactivity': '中等',
                'common_reactions': ['酯化', '氧化', '脱水'], 'polarity': 'high'
            },
            'amine': {
                'name': '氨基', 'formula': 'NH2', 'reactivity': '碱性',
                'common_reactions': ['成盐', '酰化'], 'polarity': 'medium'
            },
            'aromatic': {
                'name': '芳香环', 'formula': 'Ar', 'reactivity': '亲电取代',
                'common_reactions': ['硝化', '卤化', '磺化'], 'polarity': 'low'
            },
            'alkene': {
                'name': '碳碳双键', 'formula': 'C=C', 'reactivity': '高',
                'common_reactions': ['加成', '氧化'], 'polarity': 'low'
            },
            'alkyne': {
                'name': '碳碳三键', 'formula': 'C≡C', 'reactivity': '高',
                'common_reactions': ['加成'], 'polarity': 'low'
            },
            'halogen': {
                'name': '卤素', 'formula': 'X', 'reactivity': '亲核取代',
                'common_reactions': ['取代', '消除'], 'polarity': 'medium'
            }
        }

        # 添加官能团节点
        for fg_id, data in functional_groups.items():
            self.graph.add_node(f"fg_{fg_id}",
                              type="functional_group",
                              fg_id=fg_id,
                              **data)

        # 添加官能团间关系
        fg_relations = [
            ('carboxylic_acid', 'ester', 'can_form', '酯化反应'),
            ('carboxylic_acid', 'amide', 'can_form', '酰胺化反应'),
            ('alcohol', 'ester', 'can_form', '酯化反应'),
            ('amine', 'amide', 'can_form', '酰胺化反应'),
            ('alcohol', 'carboxylic_acid', 'related', '羧酸还原'),
            ('aromatic', 'halogen', 'can_react', '卤化反应'),
            ('alkene', 'alcohol', 'can_form', '水合反应'),
            ('alkyne', 'alkene', 'can_reduce', '部分还原')
        ]

        for fg1, fg2, relation, description in fg_relations:
            self.graph.add_edge(f"fg_{fg1}", f"fg_{fg2}",
                              relation=relation,
                              description=description)

    def _build_drug_classification_graph(self):
        """构建药物分类图谱"""
        print("   💊 构建药物分类图谱...")

        # 药物分类体系 (ATC分类)
        drug_categories = {
            'A': {'name': '消化系统和代谢', 'description': '治疗消化和代谢疾病的药物'},
            'B': {'name': '血液和造血器官', 'description': '抗凝血药、造血药等'},
            'C': {'name': '心血管系统', 'description': '心血管疾病治疗药'},
            'D': {'name': '皮肤病用药', 'description': '皮肤病治疗药'},
            'G': {'name': '泌尿生殖系统和性激素', 'description': '泌尿生殖系统药'},
            'H': {'name': '全身激素制剂', 'description': '激素类药物'},
            'J': {'name': '抗感染药', 'description': '抗生素、抗病毒药等'},
            'L': {'name': '抗肿瘤药和免疫调节剂', 'description': '化疗药、免疫抑制剂'},
            'M': {'name': '肌肉骨骼系统', 'description': '风湿病治疗药'},
            'N': {'name': '神经系统', 'description': '中枢神经系统药'},
            'P': {'name': '抗寄生虫药', 'description': '驱虫药等'},
            'R': {'name': '呼吸系统', 'description': '呼吸道疾病治疗药'},
            'S': {'name': '感觉器官', 'description': '眼科、耳鼻喉科药'},
            'V': {'name': '其他', 'description': '其他类别药物'}
        }

        # 添加药物分类节点
        for cat_id, data in drug_categories.items():
            self.graph.add_node(f"drug_cat_{cat_id}",
                              type="drug_category",
                              category_id=cat_id,
                              **data)

        # 具体药物示例
        drug_examples = {
            'aspirin': {'name': '阿司匹林', 'category': 'B', 'class': 'nsaids', 'indication': '抗血小板、抗炎'},
            'ibuprofen': {'name': '布洛芬', 'category': 'M', 'class': 'nsaids', 'indication': '抗炎、镇痛'},
            'caffeine': {'name': '咖啡因', 'category': 'N', 'class': 'psychoactive', 'indication': '中枢兴奋'},
            'penicillin': {'name': '青霉素', 'category': 'J', 'class': 'antibiotics', 'indication': '抗细菌感染'}
        }

        for drug_id, data in drug_examples.items():
            self.graph.add_node(f"drug_{drug_id}",
                              type="drug",
                              drug_id=drug_id,
                              **data)
            # 连接到分类
            self.graph.add_edge(f"drug_cat_{data['category']}", f"drug_{drug_id}",
                              relation="contains")
            self.graph.add_edge(f"class_{data['class']}", f"drug_{drug_id}",
                              relation="instance_of")

    def _build_property_association_graph(self):
        """构建性质关联图谱"""
        print("   🔗 构建性质关联图谱...")

        # 分子性质节点
        properties = {
            'mw': {'name': '分子量', 'unit': 'g/mol', 'type': 'physical', 'description': '分子质量'},
            'logp': {'name': '脂水分配系数', 'unit': 'logP', 'type': 'physicochemical', 'description': '疏水性指标'},
            'tpsa': {'name': '拓扑极性表面积', 'unit': 'Å²', 'type': 'physicochemical', 'description': '分子极性'},
            'hbd': {'name': '氢键供体数', 'unit': '个', 'type': 'structural', 'description': '氢键供体'},
            'hba': {'name': '氢键受体数', 'unit': '个', 'type': 'structural', 'description': '氢键受体'},
            'rotatable_bonds': {'name': '可旋转键数', 'unit': '个', 'type': 'structural', 'description': '分子柔性'},
            'aromatic_rings': {'name': '芳香环数', 'unit': '个', 'type': 'structural', 'description': '芳香性'}
        }

        # 添加性质节点
        for prop_id, data in properties.items():
            self.graph.add_node(f"prop_{prop_id}",
                              prop_id=prop_id,
                              **data)

        # 性质间关联
        property_relations = [
            ('mw', 'logp', 'correlates', '分子量影响脂溶性'),
            ('tpsa', 'logp', 'anti_correlates', '极性与脂溶性负相关'),
            ('hbd', 'tpsa', 'correlates', '氢键供体增加极性'),
            ('hba', 'tpsa', 'correlates', '氢键受体增加极性'),
            ('aromatic_rings', 'logp', 'correlates', '芳香环增加疏水性'),
            ('rotatable_bonds', 'mw', 'correlates', '可旋转键反映分子大小')
        ]

        for prop1, prop2, relation, description in property_relations:
            self.graph.add_edge(f"prop_{prop1}", f"prop_{prop2}",
                              relation=relation,
                              description=description)

    def get_related_concepts(self, concept_type, concept_id, max_depth=2):
        """获取相关概念"""
        start_node = f"{concept_type}_{concept_id}"

        if start_node not in self.graph:
            return []

        # 使用BFS获取相关节点
        related = []
        visited = set()
        queue = [(start_node, 0)]

        while queue:
            node, depth = queue.pop(0)
            if node in visited or depth > max_depth:
                continue

            visited.add(node)
            node_data = self.graph.nodes[node]

            # 只收集不同类型的节点
            if node_data.get('type') != concept_type:
                related.append((node, node_data))

            # 添加邻居节点
            for neighbor in self.graph.neighbors(node):
                if neighbor not in visited:
                    queue.append((neighbor, depth + 1))

        return list(related)

    def explain_molecule_context(self, smiles, molecule_name=""):
        """基于知识图谱解释分子上下文"""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return None

            # 基本分子信息
            formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
            mw = Chem.Descriptors.MolWt(mol)

            # 化学式分析
            formula_analysis = self.interpreter.explain_formula(formula)

            # 官能团识别
            functional_groups = self.interpreter.identify_functional_groups(mol)

            # 基于知识图谱的智能解释
            context_explanation = {
                'basic_info': {
                    'name': molecule_name,
                    'smiles': smiles,
                    'formula': formula,
                    'molecular_weight': round(mw, 2)
                },
                'formula_analysis': formula_analysis,
                'functional_groups': functional_groups,
                'knowledge_context': {},
                'recommendations': []
            }

            # 元素相关知识
            elements_used = list(formula_analysis['elements'].keys())
            for element in elements_used:
                related = self.get_related_concepts('element', element, max_depth=1)
                if related:
                    context_explanation['knowledge_context'][f'element_{element}'] = related

            # 官能团相关知识
            for fg in functional_groups:
                fg_id = fg['type']
                related = self.get_related_concepts('fg', fg_id, max_depth=1)
                if related:
                    context_explanation['knowledge_context'][f'fg_{fg_id}'] = related

            # 生成智能推荐
            context_explanation['recommendations'] = self._generate_recommendations(context_explanation)

            return context_explanation

        except Exception as e:
            return {'error': str(e)}

    def _generate_recommendations(self, context):
        """基于上下文生成推荐"""
        recommendations = []

        # 基于官能团的反应推荐
        fg_types = [fg['type'] for fg in context['functional_groups']]
        if 'carboxylic_acid' in fg_types:
            recommendations.append({
                'type': 'reaction',
                'title': '酯化反应',
                'description': '羧酸基可以与醇反应生成酯',
                'related_groups': ['carboxylic_acid', 'alcohol']
            })

        if 'aromatic' in fg_types:
            recommendations.append({
                'type': 'reaction',
                'title': '亲电取代反应',
                'description': '芳香环可以进行硝化、卤化等反应',
                'related_groups': ['aromatic']
            })

        # 基于性质的建议
        mw = context['basic_info']['molecular_weight']
        if mw > 500:
            recommendations.append({
                'type': 'property',
                'title': '分子量较大',
                'description': '分子量超过500，可能影响生物利用度',
                'suggestion': '考虑分子简化或增加极性基团'
            })

        # 基于元素组成的建议
        elements = context['formula_analysis']['elements']
        if elements.get('F', {}).get('count', 0) > 5:
            recommendations.append({
                'type': 'toxicity',
                'title': '含氟较多',
                'description': '分子含多个氟原子，可能增加毒性',
                'suggestion': '注意毒性评估'
            })

        return recommendations

    def visualize_knowledge_graph(self, subset=None, figsize=(12, 8)):
        """可视化知识图谱"""
        plt.figure(figsize=figsize)

        if subset:
            subgraph = self.graph.subgraph(subset)
        else:
            subgraph = self.graph

        # 计算节点位置
        pos = nx.spring_layout(subgraph, k=1, iterations=50)

        # 不同类型节点的颜色
        node_colors = []
        for node in subgraph.nodes():
            node_type = subgraph.nodes[node].get('type', 'unknown')
            if node_type == 'element':
                node_colors.append('lightblue')
            elif node_type == 'functional_group':
                node_colors.append('lightgreen')
            elif node_type == 'drug':
                node_colors.append('lightcoral')
            elif node_type.startswith('class'):
                node_colors.append('gold')
            elif node_type.startswith('prop'):
                node_colors.append('violet')
            else:
                node_colors.append('lightgray')

        # 绘制节点
        nx.draw_networkx_nodes(subgraph, pos, node_color=node_colors,
                             node_size=500, alpha=0.7)

        # 绘制边
        nx.draw_networkx_edges(subgraph, pos, alpha=0.3, arrows=True, arrowsize=20)

        # 绘制标签
        labels = {}
        for node in subgraph.nodes():
            node_data = subgraph.nodes[node]
            if node_data.get('type') == 'element':
                labels[node] = node_data.get('symbol', node)
            elif node_data.get('type') == 'functional_group':
                labels[node] = node_data.get('name', node)
            elif node_data.get('type') == 'drug':
                labels[node] = node_data.get('name', node)
            else:
                labels[node] = node_data.get('name', node.split('_')[-1])

        nx.draw_networkx_labels(subgraph, pos, labels, font_size=8, font_family='SimHei')

        plt.title('化学知识图谱', fontsize=16, fontweight='bold')
        plt.axis('off')
        plt.tight_layout()

        return plt.gcf()

    def export_knowledge_graph(self, filename='chemistry_knowledge_graph.json'):
        """导出知识图谱"""
        # 转换为可序列化的格式
        graph_data = {
            'nodes': [],
            'edges': []
        }

        for node, data in self.graph.nodes(data=True):
            graph_data['nodes'].append({
                'id': node,
                'type': data.get('type'),
                'data': data
            })

        for source, target, data in self.graph.edges(data=True):
            graph_data['edges'].append({
                'source': source,
                'target': target,
                'relation': data.get('relation'),
                'data': data
            })

        with open(filename, 'w', encoding='utf-8') as f:
            json.dump(graph_data, f, ensure_ascii=False, indent=2)

        print(f"✅ 知识图谱已导出到: {filename}")

def create_chemistry_knowledge_base():
    """创建化学知识库"""
    kg = ChemistryKnowledgeGraph()

    # 导出知识图谱
    kg.export_knowledge_graph()

    # 可视化示例
    print("📊 生成知识图谱可视化...")
    fig = kg.visualize_knowledge_graph()
    fig.savefig('chemistry_knowledge_graph.png', dpi=150, bbox_inches='tight')
    plt.close()

    return kg

if __name__ == "__main__":
    # 创建化学知识图谱
    print("🧠 构建化学知识图谱系统")
    print("=" * 50)

    kg = create_chemistry_knowledge_base()

    # 测试知识图谱功能
    print("\n🔍 测试知识图谱功能")
    print("-" * 30)

    # 测试分子上下文解释
    test_molecules = [
        ("阿司匹林", "CC(=O)OC1=CC=CC=C1C(=O)O"),
        ("咖啡因", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
    ]

    for name, smiles in test_molecules:
        print(f"\n分析分子: {name}")
        context = kg.explain_molecule_context(smiles, name)
        if context and 'error' not in context:
            print(f"  化学式: {context['basic_info']['formula']}")
            print(f"  官能团: {len(context['functional_groups'])} 个")
            print(f"  推荐: {len(context['recommendations'])} 条")
            if context['recommendations']:
                print(f"    示例: {context['recommendations'][0]['title']}")
        else:
            print(f"  分析失败: {context.get('error', '未知错误')}")

    print("\n✅ 化学知识图谱系统构建完成!")
    print("📁 输出文件:")
    print("   • chemistry_knowledge_graph.json - 知识图谱数据")
    print("   • chemistry_knowledge_graph.png - 知识图谱可视化")