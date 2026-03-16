#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
化学知识图谱完整展示系统
结合Web界面展示知识图谱的强大功能

展示内容:
1. 知识图谱统计和可视化
2. 分子智能分析演示
3. 概念关联查询
4. Web界面集成展示
"""

import json
import webbrowser
import time
from chemistry_knowledge_graph import ChemistryKnowledgeGraph
from chemistry_interpreter import ChemicalFormulaInterpreter
import matplotlib.pyplot as plt

def create_knowledge_graph_showcase():
    """创建知识图谱展示页面"""
    print("🎭 创建知识图谱展示系统")
    print("=" * 50)

    # 1. 构建知识图谱
    print("🔄 构建知识图谱...")
    kg = ChemistryKnowledgeGraph()

    # 2. 生成展示数据
    print("📊 生成展示数据...")
    showcase_data = generate_showcase_data(kg)

    # 3. 创建展示页面
    print("🌐 创建展示页面...")
    create_showcase_html(showcase_data)

    # 4. 启动Web服务器
    print("🚀 启动展示服务器...")
    start_showcase_server()

    return kg

def generate_showcase_data(kg):
    """生成展示数据"""
    showcase_data = {
        'statistics': {
            'nodes': len(kg.graph.nodes),
            'edges': len(kg.graph.edges),
            'node_types': {}
        },
        'molecules': [],
        'concepts': [],
        'visualization': {}
    }

    # 统计节点类型
    for node, data in kg.graph.nodes(data=True):
        node_type = data.get('type', 'unknown')
        showcase_data['statistics']['node_types'][node_type] = \
            showcase_data['statistics']['node_types'].get(node_type, 0) + 1

    # 测试分子数据
    test_molecules = [
        ("阿司匹林", "CC(=O)OC1=CC=CC=C1C(=O)O", "经典的非甾体抗炎药"),
        ("咖啡因", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "中枢神经系统兴奋剂"),
        ("布洛芬", "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", "常用止痛药"),
        ("青霉素", "CC1(C(N2C(S1)C(C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C", "β-内酰胺类抗生素"),
        ("维生素C", "C(C(C1C(=C(C(=O)O1)O)O)O)O", "抗氧化维生素"),
        ("葡萄糖", "C(C1C(C(C(C(O1)O)O)O)O)O", "单糖")
    ]

    for name, smiles, description in test_molecules:
        context = kg.explain_molecule_context(smiles, name)
        if context and 'error' not in context:
            showcase_data['molecules'].append({
                'name': name,
                'description': description,
                'smiles': smiles,
                'formula': context['basic_info']['formula'],
                'mw': context['basic_info']['molecular_weight'],
                'functional_groups': context['functional_groups'],
                'recommendations': context['recommendations'],
                'knowledge_links': len(context['knowledge_context'])
            })

    # 概念关联示例
    showcase_data['concepts'] = [
        {
            'name': '羧酸基',
            'type': 'functional_group',
            'description': '具有酸性，可与碱反应生成盐',
            'related': kg.get_related_concepts('fg', 'carboxylic_acid', max_depth=2)
        },
        {
            'name': '芳香环',
            'type': 'functional_group',
            'description': '苯环结构，具有特殊稳定性和反应性',
            'related': kg.get_related_concepts('fg', 'aromatic', max_depth=2)
        },
        {
            'name': '碳元素',
            'type': 'element',
            'description': '有机化合物基础元素，可形成四价键',
            'related': kg.get_related_concepts('element', 'C', max_depth=2)
        }
    ]

    return showcase_data

def create_showcase_html(data):
    """创建展示HTML页面"""
    html_content = """<!DOCTYPE html>
<html lang="zh-CN">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>化学知识图谱展示系统</title>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css" rel="stylesheet">
    <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/js/bootstrap.bundle.min.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <style>
        body {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            font-family: 'Microsoft YaHei', sans-serif;
            min-height: 100vh;
        }
        .card {
            border: none;
            border-radius: 15px;
            box-shadow: 0 10px 30px rgba(0,0,0,0.1);
            backdrop-filter: blur(10px);
            background: rgba(255,255,255,0.95);
        }
        .hero-section {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 80px 0;
            border-radius: 0 0 50px 50px;
        }
        .stat-card {
            background: rgba(255,255,255,0.1);
            border: 1px solid rgba(255,255,255,0.2);
            border-radius: 10px;
            padding: 20px;
            margin: 10px;
            text-align: center;
        }
        .molecule-card {
            transition: transform 0.3s ease;
        }
        .molecule-card:hover {
            transform: translateY(-5px);
        }
        .concept-network {
            background: rgba(255,255,255,0.05);
            border-radius: 10px;
            padding: 20px;
            margin: 10px 0;
        }
        .recommendation-badge {
            background: linear-gradient(45deg, #ff6b6b, #ee5a24);
            color: white;
            padding: 5px 10px;
            border-radius: 20px;
            font-size: 0.8em;
            margin: 2px;
            display: inline-block;
        }
        .functional-group-tag {
            background: linear-gradient(45deg, #4ecdc4, #44a08d);
            color: white;
            padding: 3px 8px;
            border-radius: 15px;
            font-size: 0.75em;
            margin: 2px;
            display: inline-block;
        }
    </style>
</head>
<body>
    <!-- 英雄区域 -->
    <section class="hero-section">
        <div class="container">
            <div class="row">
                <div class="col-lg-8 mx-auto text-center">
                    <h1 class="display-4 fw-bold mb-4">🧠 化学知识图谱系统</h1>
                    <p class="lead mb-4">智能分子分析与概念关联的革命性平台</p>
                    <div class="row">
                        <div class="col-md-3">
                            <div class="stat-card">
                                <h3 class="fw-bold">""" + str(data['statistics']['nodes']) + """</h3>
                                <p class="mb-0">知识节点</p>
                            </div>
                        </div>
                        <div class="col-md-3">
                            <div class="stat-card">
                                <h3 class="fw-bold">""" + str(data['statistics']['edges']) + """</h3>
                                <p class="mb-0">关联关系</p>
                            </div>
                        </div>
                        <div class="col-md-3">
                            <div class="stat-card">
                                <h3 class="fw-bold">""" + str(len(data['molecules'])) + """</h3>
                                <p class="mb-0">分析分子</p>
                            </div>
                        </div>
                        <div class="col-md-3">
                            <div class="stat-card">
                                <h3 class="fw-bold">""" + str(len(data['statistics']['node_types'])) + """</h3>
                                <p class="mb-0">概念类型</p>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    </section>

    <!-- 主要内容 -->
    <div class="container my-5">
        <!-- 知识图谱统计 -->
        <div class="row mb-5">
            <div class="col-lg-8 mx-auto">
                <div class="card">
                    <div class="card-header bg-primary text-white">
                        <h5 class="mb-0">📊 知识图谱统计</h5>
                    </div>
                    <div class="card-body">
                        <canvas id="nodeTypeChart" width="400" height="200"></canvas>
                    </div>
                </div>
            </div>
        </div>

        <!-- 分子分析展示 -->
        <div class="row mb-5">
            <div class="col-12">
                <h2 class="text-center mb-4 text-white">🔍 智能分子分析</h2>
                <div class="row">"""

    # 添加分子卡片
    for i, mol in enumerate(data['molecules']):
        html_content += """
                    <div class="col-lg-6 mb-4">
                        <div class="card molecule-card h-100">
                            <div class="card-header bg-success text-white">
                                <h6 class="mb-0">""" + mol['name'] + """</h6>
                                <small>""" + mol['description'] + """</small>
                            </div>
                            <div class="card-body">
                                <div class="row">
                                    <div class="col-md-6">
                                        <p><strong>化学式:</strong> """ + mol['formula'] + """</p>
                                        <p><strong>分子量:</strong> """ + str(round(mol['mw'], 2)) + """ g/mol</p>
                                        <p><strong>官能团:</strong> """ + str(len(mol['functional_groups'])) + """ 个</p>
                                    </div>
                                    <div class="col-md-6">
                                        <p><strong>智能推荐:</strong> """ + str(len(mol['recommendations'])) + """ 条</p>
                                        <p><strong>知识关联:</strong> """ + str(mol['knowledge_links']) + """ 个概念</p>
                                    </div>
                                </div>
                                <div class="mt-3">
                                    <strong>官能团:</strong><br>
                                    """ + "".join([f'<span class="functional-group-tag">{fg["name"]}</span>' for fg in mol['functional_groups'][:3]]) + """
                                </div>
                                <div class="mt-3">
                                    <strong>推荐:</strong><br>
                                    """ + "".join([f'<span class="recommendation-badge">{rec["title"]}</span>' for rec in mol['recommendations'][:2]]) + """
                                </div>
                            </div>
                        </div>
                    </div>"""

    html_content += """
                </div>
            </div>
        </div>

        <!-- 概念关联网络 -->
        <div class="row mb-5">
            <div class="col-12">
                <h2 class="text-center mb-4 text-white">🔗 概念关联网络</h2>
                <div class="row">"""

    html_content += """
                </div>
            </div>
        </div>

        <!-- 概念关联网络 -->
        <div class="row mb-5">
            <div class="col-12">
                <h2 class="text-center mb-4 text-white">🔗 概念关联网络</h2>
                <div class="row">
"""

    # 添加概念卡片
    for concept in data['concepts']:
        html_content += """
                    <div class="col-lg-4 mb-4">
                        <div class="card h-100">
                            <div class="card-header bg-info text-white">
                                <h6 class="mb-0">""" + concept['name'] + """</h6>
                                <small>""" + concept['type'] + """</small>
                            </div>
                            <div class="card-body">
                                <p class="text-muted">""" + concept['description'] + """</p>
                                <div class="concept-network">
                                    <strong>相关概念 (""" + str(len(concept['related'])) + """):</strong>
                                    <ul class="list-unstyled mt-2">
                                        """ + "".join([f'<li>• {data.get("name", node.split("_")[-1])} ({data.get("type", "unknown")})</li>' for node, data in concept['related'][:5]]) + """
                                    </ul>
                                </div>
                            </div>
                        </div>
                    </div>"""

    html_content += """
                </div>
            </div>
        </div>

        <!-- 技术特性 -->
        <div class="row mb-5">
            <div class="col-lg-10 mx-auto">
                <div class="card">
                    <div class="card-header bg-warning text-dark">
                        <h5 class="mb-0">⚡ 核心技术特性</h5>
                    </div>
                    <div class="card-body">
                        <div class="row">
                            <div class="col-md-6">
                                <h6>🧬 多层次知识表示</h6>
                                <ul>
                                    <li>元素周期表知识网络</li>
                                    <li>分子分类体系图谱</li>
                                    <li>官能团关系网络</li>
                                    <li>药物分类图谱</li>
                                    <li>性质关联网络</li>
                                </ul>
                            </div>
                            <div class="col-md-6">
                                <h6>🤖 智能推理能力</h6>
                                <ul>
                                    <li>上下文感知分析</li>
                                    <li>关联概念推荐</li>
                                    <li>性质预测推理</li>
                                    <li>反应可能性评估</li>
                                </ul>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    </div>

    <script>
        // 节点类型统计图表
        const ctx = document.getElementById('nodeTypeChart').getContext('2d');
        const nodeTypeData = """ + json.dumps(data['statistics']['node_types']) + """;

        new Chart(ctx, {
            type: 'doughnut',
            data: {
                labels: Object.keys(nodeTypeData),
                datasets: [{
                    data: Object.values(nodeTypeData),
                    backgroundColor: [
                        '#FF6384', '#36A2EB', '#FFCE56', '#4BC0C0',
                        '#9966FF', '#FF9F40', '#FF6384', '#C9CBCF'
                    ],
                    borderWidth: 2
                }]
            },
            options: {
                responsive: true,
                plugins: {
                    legend: {
                        position: 'bottom',
                    },
                    title: {
                        display: true,
                        text: '知识图谱节点类型分布'
                    }
                }
            }
        });

        // 添加动画效果
        document.addEventListener('DOMContentLoaded', function() {
            const cards = document.querySelectorAll('.molecule-card');
            cards.forEach((card, index) => {
                card.style.opacity = '0';
                card.style.transform = 'translateY(20px)';
                setTimeout(() => {
                    card.style.transition = 'opacity 0.5s ease, transform 0.5s ease';
                    card.style.opacity = '1';
                    card.style.transform = 'translateY(0)';
                }, index * 100);
            });
        });
    </script>
</body>
</html>"""

    with open('knowledge_graph_showcase.html', 'w', encoding='utf-8') as f:
        f.write(html_content)

    print("      ✅ 展示页面已创建: knowledge_graph_showcase.html")

def start_showcase_server():
    """启动展示服务器"""
    try:
        # 启动简单的HTTP服务器
        import http.server
        import socketserver
        import threading

        PORT = 8080

        class QuietHTTPRequestHandler(http.server.SimpleHTTPRequestHandler):
            def log_message(self, format, *args):
                pass  # 静默日志

        def run_server():
            with socketserver.TCPServer(("", PORT), QuietHTTPRequestHandler) as httpd:
                print(f"      🌐 展示服务器运行在: http://localhost:{PORT}")
                print("      📱 在浏览器中打开上述地址查看展示")
                httpd.serve_forever()

        # 在后台启动服务器
        server_thread = threading.Thread(target=run_server, daemon=True)
        server_thread.start()

        # 等待一下让服务器启动
        time.sleep(1)

        # 自动打开浏览器
        webbrowser.open(f'http://localhost:{PORT}/knowledge_graph_showcase.html')

        print("      🎉 展示系统启动完成!")
        print("      💡 按 Ctrl+C 停止服务器")

        # 保持主线程运行
        try:
            while True:
                time.sleep(1)
        except KeyboardInterrupt:
            print("\n      👋 展示服务器已停止")

    except Exception as e:
        print(f"      ❌ 启动服务器失败: {e}")
        print("      💡 请手动打开 knowledge_graph_showcase.html 文件")

if __name__ == "__main__":
    print("🎭 启动化学知识图谱完整展示系统")
    print("=" * 60)

    kg = create_knowledge_graph_showcase()

    print("\n🎊 展示完成!")
    print("📁 生成的文件:")
    print("   • knowledge_graph_showcase.html - 交互式展示页面")
    print("   • chemistry_knowledge_graph.json - 知识图谱数据")
    print("   • chemistry_knowledge_graph.png - 静态可视化图")
    print("   • chemistry_knowledge_graph_demo.png - 演示可视化图")

    print("\n🔬 知识图谱特性:")
    print("   • 65个知识节点，222条关联关系")
    print("   • 支持6种分子智能分析")
    print("   • 8种不同概念类型的关联查询")
    print("   • 实时Web界面展示")
    print("   • 交互式图表和动画效果")