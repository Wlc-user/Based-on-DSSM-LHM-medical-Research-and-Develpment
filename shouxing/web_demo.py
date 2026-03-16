#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Web界面演示脚本
展示分子相似性检索系统的Web界面功能
"""

import requests
import json
import time
import webbrowser
from datetime import datetime

def test_web_interface():
    """测试Web界面功能"""
    print("=" * 70)
    print("🌐 Web界面功能演示")
    print("=" * 70)

    base_url = "http://localhost:5000"

    # 测试1: SMILES验证
    print("\n1️⃣ 测试SMILES验证功能")
    test_smiles = [
        "CC(=O)OC1=CC=CC=C1C(=O)O",  # 阿司匹林
        "CC1=CC=C(C=C1)C2=CC(=NN2C3=CC=C(C=C3)S(=O)(=O)N)C(F)(F)F",  # 洛索洛芬
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",  # 咖啡因
        "INVALID_SMILES"  # 无效SMILES
    ]

    for smiles in test_smiles:
        try:
            response = requests.post(f"{base_url}/api/validate_smiles",
                                   json={"smiles": smiles}, timeout=5)
            result = response.json()
            if result.get('valid'):
                props = result['properties']
                print(f"✅ {smiles[:30]}... | MW: {props['molecular_weight']:.1f} | 公式: {props['formula']}")
            else:
                print(f"❌ {smiles[:30]}... | 错误: {result.get('error', '未知错误')}")
        except Exception as e:
            print(f"⚠️  {smiles[:30]}... | 网络错误: {str(e)}")

    # 测试2: 分子库统计
    print("\n2️⃣ 测试分子库统计")
    try:
        response = requests.get(f"{base_url}/api/library_stats", timeout=5)
        stats = response.json()
        print(f"📊 分子库统计:")
        print(f"   总化合物数: {stats.get('total_compounds', 'N/A')}")
        print(f"   数据源: {list(stats.get('data_sources', {}).keys())}")
        print(f"   类别: {list(stats.get('categories', {}).keys())[:5]}...")  # 只显示前5个
    except Exception as e:
        print(f"⚠️ 库统计获取失败: {str(e)}")

    # 测试3: 分子相似性搜索
    print("\n3️⃣ 测试分子相似性搜索")
    search_tests = [
        {
            "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",  # 阿司匹林
            "name": "阿司匹林",
            "threshold": 0.3,
            "top_k": 3
        },
        {
            "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",  # 咖啡因
            "name": "咖啡因",
            "threshold": 0.2,
            "top_k": 3
        }
    ]

    for test in search_tests:
        try:
            print(f"\n🔍 搜索与 {test['name']} 相似的分子:")
            response = requests.post(f"{base_url}/api/search",
                                   json=test, timeout=10)
            result = response.json()

            if 'results' in result and result['results']:
                print(f"   找到 {len(result['results'])} 个相似分子")
                print(f"   处理时间: {result['performance']['elapsed_time']:.3f}秒")
                print(f"   吞吐量: {result['performance']['throughput']:.1f} 个/秒")

                for i, r in enumerate(result['results'][:3]):
                    print(f"   {i+1}. {r['name']} (相似度: {r['weighted_similarity']:.3f})")
            else:
                print("   未找到相似分子")

        except Exception as e:
            print(f"⚠️ 搜索失败: {str(e)}")

    # 测试4: 批量搜索
    print("\n4️⃣ 测试批量搜索功能")
    try:
        batch_data = {
            "smiles_list": [
                "CC(=O)OC1=CC=CC=C1C(=O)O",  # 阿司匹林
                "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",  # 咖啡因
            ],
            "threshold": 0.3,
            "top_k": 2
        }

        response = requests.post(f"{base_url}/api/batch_search",
                               json=batch_data, timeout=15)
        batch_result = response.json()

        print(f"📦 批量搜索完成: 处理了 {batch_result.get('batch_size', 0)} 个分子")
        for item in batch_result.get('results', []):
            if 'error' not in item:
                smiles = item['query_smiles'][:20] + "..."
                count = item.get('results_count', 0)
                print(f"   {smiles} -> 找到 {count} 个相似分子")
            else:
                print(f"   {item['query_smiles'][:20]}... -> 错误: {item['error']}")

    except Exception as e:
        print(f"⚠️ 批量搜索失败: {str(e)}")

    # 测试5: 分子结构图像生成
    print("\n5️⃣ 测试分子结构图像生成")
    try:
        test_smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"
        response = requests.get(f"{base_url}/api/molecule_image/{test_smiles}", timeout=5)
        if response.status_code == 200:
            print("✅ 分子结构图像生成成功")
        else:
            print("❌ 分子结构图像生成失败")
    except Exception as e:
        print(f"⚠️ 图像生成测试失败: {str(e)}")

    print("\n" + "=" * 70)
    print("🎉 Web界面功能测试完成!")
    print("📱 访问Web界面: http://localhost:5000")
    print("🔍 搜索页面: http://localhost:5000/search")
    print("=" * 70)

def open_web_interface():
    """打开Web界面"""
    print("\n🌐 正在打开Web界面...")
    try:
        webbrowser.open("http://localhost:5000")
        print("✅ Web界面已在浏览器中打开")
    except Exception as e:
        print(f"⚠️ 无法自动打开浏览器: {str(e)}")
        print("请手动访问: http://localhost:5000")

def main():
    """主函数"""
    print("🧬 分子相似性检索系统 - Web界面演示")
    print(f"⏰ 开始时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    # 等待Web服务器启动
    print("⏳ 等待Web服务器启动...")
    time.sleep(3)

    # 运行功能测试
    test_web_interface()

    # 打开Web界面
    open_web_interface()

    print(f"\n🏁 演示完成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

if __name__ == "__main__":
    main()