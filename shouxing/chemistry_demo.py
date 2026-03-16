#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
化学知识鸿沟演示脚本
展示如何弥合化学式与用户常识之间的知识差距
"""

import requests
import json
import webbrowser
from chemistry_interpreter import ChemicalFormulaInterpreter

def demonstrate_chemistry_gap():
    """演示化学知识鸿沟的弥合"""
    print("=" * 80)
    print("🧪 化学知识鸿沟弥合演示")
    print("=" * 80)

    interpreter = ChemicalFormulaInterpreter()

    # 示例1: 阿司匹林的化学式解释
    print("\n1️⃣ 阿司匹林的化学式解释")
    print("-" * 40)

    aspirin_formula = "C9H8O4"
    explanation = interpreter.explain_formula(aspirin_formula)

    print(f"化学式: {aspirin_formula}")
    print(f"中文解释: {explanation['chinese_text']}")
    print("元素组成:")
    for element, info in explanation['elements'].items():
        print(f"  {element} ({info['chinese_name']}): {info['count']}个原子 ({info['percentage']}%)")

    # 示例2: 分子结构分析
    print("\n2️⃣ 分子结构分析")
    print("-" * 40)

    molecules = [
        ("阿司匹林", "CC(=O)OC1=CC=CC=C1C(=O)O"),
        ("咖啡因", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"),
        ("葡萄糖", "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O")
    ]

    for name, smiles in molecules:
        print(f"\n分析分子: {name}")
        insights = interpreter.get_molecule_insights(smiles)
        if insights and 'error' not in insights:
            basic = insights['basic_info']
            print(f"  原子数: {basic['num_atoms']}, 化学键数: {basic['num_bonds']}")
            print(f"  分子量: {basic['molecular_weight']:.1f} g/mol")
            print(f"  化学式: {basic['formula']}")

            # 显示官能团
            if insights.get('functional_groups'):
                groups = [f"{g['name']}({g['count']}个)" for g in insights['functional_groups']]
                print(f"  官能团: {', '.join(groups)}")

            # 显示用户友好的总结
            summary = insights.get('user_summary', '').split('\n')[0]
            print(f"  总结: {summary}")
        else:
            print(f"  分析失败: {insights.get('error', '未知错误') if insights else '无结果'}")

    # 示例3: Web界面演示
    print("\n3️⃣ Web界面化学知识功能")
    print("-" * 40)

    base_url = "http://localhost:5000"

    # 测试化学式解释API
    try:
        response = requests.post(f"{base_url}/api/explain_formula",
                               json={"formula": "C8H10N4O2"}, timeout=5)
        if response.status_code == 200:
            result = response.json()
            print("✅ 化学式解释API: 正常")
            print(f"   示例: {result['chinese_text'][:60]}...")
        else:
            print("❌ 化学式解释API: 失败")
    except Exception as e:
        print(f"⚠️ 化学式解释API测试失败: {str(e)}")

    # 测试分子分析API
    try:
        response = requests.post(f"{base_url}/api/analyze_molecule",
                               json={"smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"}, timeout=5)
        if response.status_code == 200:
            result = response.json()
            print("✅ 分子分析API: 正常")
            print(f"   咖啡因分子式: {result['basic_info']['formula']}")
        else:
            print("❌ 分子分析API: 失败")
    except Exception as e:
        print(f"⚠️ 分子分析API测试失败: {str(e)}")

    # 测试元素知识库
    try:
        response = requests.get(f"{base_url}/api/element_knowledge", timeout=5)
        if response.status_code == 200:
            knowledge = response.json()
            print("✅ 元素知识库API: 正常")
            print(f"   常见元素数量: {len(knowledge['common_elements'])}")
        else:
            print("❌ 元素知识库API: 失败")
    except Exception as e:
        print(f"⚠️ 元素知识库API测试失败: {str(e)}")

def show_chemistry_education_content():
    """展示化学教育内容"""
    print("\n4️⃣ 化学教育内容示例")
    print("-" * 40)

    try:
        with open('chemistry_education_content.json', 'r', encoding='utf-8') as f:
            content = json.load(f)

        print("📚 元素知识库:")
        elements = content['element_knowledge']['common_elements']
        sample_elements = dict(list(elements.items())[:5])
        for symbol, name in sample_elements.items():
            print(f"   {symbol} -> {name}")

        print("\n🧬 示例分子分析:")
        if content['examples']:
            example = content['examples'][0]
            basic_info = example['insights']['basic_info']
            print(f"   分子: {example['name']}")
            print(f"   化学式: {basic_info['formula']}")
            print(f"   分子量: {basic_info['molecular_weight']:.1f}")

    except FileNotFoundError:
        print("⚠️ 化学教育内容文件不存在，请先运行 chemistry_interpreter.py")
    except Exception as e:
        print(f"⚠️ 读取教育内容失败: {str(e)}")

def open_chemistry_demo():
    """打开化学知识演示页面"""
    print("\n🌐 打开Web界面演示...")
    try:
        webbrowser.open("http://localhost:5000/search")
        print("✅ 已打开分子搜索页面")
        print("💡 提示: 在搜索结果中点击'解释'按钮查看化学式解释")
        print("🔬 提示: 点击'分析'按钮查看分子结构分析")
    except Exception as e:
        print(f"⚠️ 无法自动打开浏览器: {str(e)}")
        print("请手动访问: http://localhost:5000/search")

def demonstrate_knowledge_gap():
    """演示知识鸿沟的具体表现"""
    print("\n5️⃣ 知识鸿沟分析")
    print("-" * 40)

    print("🔍 鸿沟表现:")
    print("   • 用户看到: '阿司匹林' (熟悉的中文名称)")
    print("   • 系统显示: 'C9H8O4' (陌生的化学符号)")
    print("   • 用户疑惑: 这些字母和数字是什么意思?")
    print()
    print("🌉 弥合方案:")
    print("   • 化学式解释: 将'C9H8O4'转换为'9个碳原子、8个氢原子、4个氧原子'")
    print("   • 分子分析: 识别官能团、计算性质、生成结构图像")
    print("   • 教育内容: 提供元素知识和分子背景信息")
    print("   • 交互界面: 一键获取详细解释，无需专业知识")

    print("\n📊 效果评估:")
    print("   • 知识门槛: 从需要化学专业知识 → 普通用户可理解")
    print("   • 用户体验: 从困惑 → 清晰明了")
    print("   • 学习效果: 从被动接受 → 主动探索")

def main():
    """主函数"""
    print("🧬 分子相似性检索系统 - 化学知识鸿沟弥合演示")
    print(" bridging the gap between chemical formulas and common knowledge")

    # 演示知识鸿沟弥合
    demonstrate_chemistry_gap()

    # 展示教育内容
    show_chemistry_education_content()

    # 分析知识鸿沟
    demonstrate_knowledge_gap()

    # 打开Web演示
    open_chemistry_demo()

    print("\n" + "=" * 80)
    print("🎯 总结: 通过技术手段弥合化学知识鸿沟")
    print("   • 化学式可视化解释")
    print("   • 分子结构智能分析")
    print("   • 用户友好的交互界面")
    print("   • 渐进式的知识普及")
    print("=" * 80)

if __name__ == "__main__":
    main()