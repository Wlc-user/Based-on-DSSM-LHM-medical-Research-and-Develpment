#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
化学式解释器模块
帮助用户理解化学式和分子性质的知识体系

解决化学知识与用户常识之间的鸿沟：
- 化学式解析和可视化
- 分子组成元素解释
- 分子性质含义说明
- 中英文对照表
"""

import re
from collections import defaultdict
import json
import os
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

class ChemicalFormulaInterpreter:
    """化学式解释器"""

    def __init__(self):
        # 元素符号到中文名称的映射
        self.element_names = {
            'H': '氢',
            'C': '碳',
            'N': '氮',
            'O': '氧',
            'F': '氟',
            'Cl': '氯',
            'Br': '溴',
            'I': '碘',
            'S': '硫',
            'P': '磷',
            'Na': '钠',
            'K': '钾',
            'Ca': '钙',
            'Mg': '镁',
            'Fe': '铁',
            'Cu': '铜',
            'Zn': '锌',
            'Al': '铝',
            'Si': '硅',
            'B': '硼'
        }

        # 分子性质解释
        self.property_explanations = {
            'mw': {
                'name': '分子量',
                'unit': 'g/mol',
                'description': '分子中所有原子的质量总和',
                'importance': '影响药物的吸收、分布和代谢'
            },
            'logp': {
                'name': '脂水分配系数',
                'unit': 'logP',
                'description': '化合物在油水两相中的分配倾向',
                'importance': '影响药物的生物利用度和细胞膜通透性'
            },
            'tpsa': {
                'name': '拓扑极性表面积',
                'unit': 'Å²',
                'description': '分子中所有极性原子(氮、氧、硫等)的表面积',
                'importance': '影响药物的水溶性和生物活性'
            },
            'hbd': {
                'name': '氢键供体数',
                'unit': '个',
                'description': '能提供氢键的原子数(通常是-OH、-NH等)',
                'importance': '影响分子间的相互作用和药效'
            },
            'hba': {
                'name': '氢键受体数',
                'unit': '个',
                'description': '能接受氢键的原子数(通常是氮、氧等)',
                'importance': '影响分子间的相互作用'
            }
        }

        # 常见官能团识别
        self.functional_groups = {
            'carboxylic_acid': {
                'pattern': 'C(=O)O',
                'name': '羧酸基',
                'description': '具有酸性，可与碱反应'
            },
            'ester': {
                'pattern': 'C(=O)O',
                'name': '酯基',
                'description': '羧酸衍生物，具有酯的香味'
            },
            'amide': {
                'pattern': 'C(=O)N',
                'name': '酰胺基',
                'description': '蛋白质的基本结构单元'
            },
            'alcohol': {
                'pattern': 'O',
                'name': '羟基',
                'description': '具有-OH基团，可形成氢键'
            },
            'amine': {
                'pattern': 'N',
                'name': '氨基',
                'description': '含氮化合物，具有碱性'
            },
            'aromatic': {
                'pattern': 'c1ccccc1',
                'name': '芳香环',
                'description': '苯环结构，具有特殊稳定性'
            }
        }

    def parse_formula(self, formula):
        """解析化学式，返回元素组成"""
        if not formula:
            return {}

        # 使用正则表达式解析化学式
        # 匹配元素符号和数字
        pattern = r'([A-Z][a-z]?)(\d*)'
        matches = re.findall(pattern, formula)

        composition = defaultdict(int)
        for element, count in matches:
            count = int(count) if count else 1
            composition[element] += count

        return dict(composition)

    def explain_formula(self, formula):
        """解释化学式的组成"""
        composition = self.parse_formula(formula)

        explanation = {
            'formula': formula,
            'elements': {},
            'total_atoms': sum(composition.values()),
            'element_count': len(composition),
            'chinese_explanation': []
        }

        for element, count in composition.items():
            element_info = {
                'symbol': element,
                'count': count,
                'chinese_name': self.element_names.get(element, element),
                'percentage': round(count / explanation['total_atoms'] * 100, 1)
            }
            explanation['elements'][element] = element_info

            # 生成中文解释
            chinese_name = self.element_names.get(element, f'元素{element}')
            if count == 1:
                explanation['chinese_explanation'].append(f'1个{chinese_name}原子')
            else:
                explanation['chinese_explanation'].append(f'{count}个{chinese_name}原子')

        explanation['chinese_text'] = f"分子式{formula}表示该分子由{'、'.join(explanation['chinese_explanation'])}组成。"

        return explanation

    def explain_properties(self, properties):
        """解释分子性质"""
        explanations = {}

        for prop_key, value in properties.items():
            if prop_key in self.property_explanations:
                prop_info = self.property_explanations[prop_key].copy()
                prop_info['value'] = value
                explanations[prop_key] = prop_info

        return explanations

    def identify_functional_groups(self, mol):
        """识别分子中的官能团"""
        if not mol:
            return []

        groups_found = []

        # 检查各种官能团
        for group_key, group_info in self.functional_groups.items():
            try:
                # 使用SMARTS模式匹配
                if 'pattern' in group_info:
                    pattern = Chem.MolFromSmarts(group_info['pattern'])
                    if pattern and mol.HasSubstructMatch(pattern):
                        groups_found.append({
                            'type': group_key,
                            'name': group_info['name'],
                            'description': group_info['description'],
                            'count': len(mol.GetSubstructMatches(pattern))
                        })
            except:
                continue

        return groups_found

    def get_molecule_insights(self, smiles, properties=None):
        """获取分子的综合见解"""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return None

            insights = {
                'smiles': smiles,
                'basic_info': {
                    'num_atoms': mol.GetNumAtoms(),
                    'num_bonds': mol.GetNumBonds(),
                    'num_rings': mol.GetRingInfo().NumRings(),
                    'molecular_weight': Descriptors.MolWt(mol),
                    'formula': Chem.rdMolDescriptors.CalcMolFormula(mol)
                }
            }

            # 化学式解释
            formula_explanation = self.explain_formula(insights['basic_info']['formula'])
            insights['formula_explanation'] = formula_explanation

            # 性质解释
            if properties:
                insights['property_explanations'] = self.explain_properties(properties)

            # 官能团识别
            insights['functional_groups'] = self.identify_functional_groups(mol)

            # 生成用户友好的总结
            insights['user_summary'] = self._generate_user_summary(insights)

            return insights

        except Exception as e:
            return {'error': str(e)}

    def _generate_user_summary(self, insights):
        """生成用户友好的分子总结"""
        basic = insights['basic_info']
        formula_exp = insights['formula_explanation']

        summary = f"""
这个分子有{basic['num_atoms']}个原子，{basic['num_bonds']}个化学键。
{formula_exp['chinese_text']}
分子量为{basic['molecular_weight']:.1f}，表示这个分子相对较{'大' if basic['molecular_weight'] > 200 else '小'}。
"""

        if insights['functional_groups']:
            groups_text = [f"{g['name']}({g['count']}个)" for g in insights['functional_groups']]
            summary += f"分子中含有以下官能团：{', '.join(groups_text)}。\n"

        if 'property_explanations' in insights:
            summary += "分子性质分析：\n"
            for prop_key, prop_info in insights['property_explanations'].items():
                summary += f"- {prop_info['name']}: {prop_info['value']} {prop_info['unit']} - {prop_info['description']}\n"

        return summary.strip()

    def get_element_knowledge(self):
        """获取元素知识库"""
        return {
            'common_elements': self.element_names,
            'periodic_trends': {
                'electronegativity': '元素吸引电子的能力',
                'atomic_radius': '原子大小',
                'ionization_energy': '失去电子的难易程度'
            },
            'bio_relevance': {
                'C': '生命的基础元素，形成有机化合物',
                'H': '构成水和有机化合物',
                'O': '氧化反应，构成水和二氧化碳',
                'N': '构成蛋白质、DNA和RNA',
                'P': '构成DNA、RNA和ATP',
                'S': '构成蛋白质中的半胱氨酸和甲硫氨酸'
            }
        }

def create_chemistry_education_content():
    """创建化学教育内容"""
    interpreter = ChemicalFormulaInterpreter()

    # 示例分子分析
    examples = [
        {
            'name': '阿司匹林',
            'smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O',
            'properties': {'mw': 180.16, 'logp': 1.19, 'tpsa': 63.6}
        },
        {
            'name': '咖啡因',
            'smiles': 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
            'properties': {'mw': 194.19, 'logp': -0.07, 'tpsa': 58.2}
        },
        {
            'name': '葡萄糖',
            'smiles': 'OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O',
            'properties': {'mw': 180.16, 'logp': -3.24, 'tpsa': 110.4}
        }
    ]

    education_content = {
        'element_knowledge': interpreter.get_element_knowledge(),
        'examples': []
    }

    for example in examples:
        insights = interpreter.get_molecule_insights(
            example['smiles'],
            example['properties']
        )
        if insights:
            education_content['examples'].append({
                'name': example['name'],
                'insights': insights
            })

    return education_content

if __name__ == "__main__":
    # 测试化学式解释器
    interpreter = ChemicalFormulaInterpreter()

    # 测试化学式解析
    print("🧪 化学式解释器测试")
    print("=" * 50)

    test_formulas = ['C9H8O4', 'C8H10N4O2', 'C6H12O6']
    for formula in test_formulas:
        explanation = interpreter.explain_formula(formula)
        print(f"化学式: {formula}")
        print(f"中文解释: {explanation['chinese_text']}")
        print(f"元素组成: {explanation['elements']}")
        print()

    # 测试分子分析
    print("🔬 分子分析测试")
    print("=" * 50)

    test_molecules = [
        ('阿司匹林', 'CC(=O)OC1=CC=CC=C1C(=O)O'),
        ('咖啡因', 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C')
    ]

    for name, smiles in test_molecules:
        print(f"分子: {name}")
        insights = interpreter.get_molecule_insights(smiles)
        if insights and 'error' not in insights:
            print(f"基本信息: {insights['basic_info']}")
            print(f"用户总结: {insights['user_summary'][:200]}...")
        else:
            print(f"分析失败: {insights.get('error', '未知错误') if insights else '无结果'}")
        print()

    # 生成教育内容
    print("📚 生成化学教育内容...")
    education_content = create_chemistry_education_content()

    # 保存到文件
    output_file = 'chemistry_education_content.json'
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(education_content, f, ensure_ascii=False, indent=2)

    print(f"✅ 教育内容已保存到: {output_file}")