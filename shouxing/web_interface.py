#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
第二阶段: Web界面系统
基于Flask的基础Web应用，提供分子相似性搜索和虚拟筛选

新增功能:
- 分子SMILES输入界面
- 实时相似性搜索结果
- 交互式可视化图表
- 结果导出功能
- RESTful API接口
"""

import os
import sys
import json
import base64
import io
from flask import Flask, render_template, request, jsonify, send_file
from flask_cors import CORS
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from rdkit import Chem
from rdkit.Chem import Draw
import warnings
warnings.filterwarnings('ignore')

# 导入我们的增强版虚拟筛选系统
from enhanced_virtual_screening import EnhancedVirtualScreeningSystem
# 导入化学式解释器
from chemistry_interpreter import ChemicalFormulaInterpreter

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

app = Flask(__name__)
CORS(app)  # 启用跨域支持

# 初始化虚拟筛选系统
vs_system = None
# 初始化化学式解释器
chemistry_interpreter = ChemicalFormulaInterpreter()

def get_vs_system():
    """获取虚拟筛选系统实例"""
    global vs_system
    if vs_system is None:
        vs_system = EnhancedVirtualScreeningSystem()
    return vs_system

@app.route('/')
def index():
    """主页"""
    return render_template('index.html')

@app.route('/search')
def search_page():
    """搜索页面"""
    return render_template('search.html')

@app.route('/api/search', methods=['POST'])
def api_search():
    """分子相似性搜索API"""
    try:
        data = request.get_json()
        if not data or 'smiles' not in data:
            return jsonify({'error': '缺少SMILES参数'}), 400

        smiles = data['smiles'].strip()
        if not smiles:
            return jsonify({'error': 'SMILES不能为空'}), 400

        # 参数设置
        threshold = float(data.get('threshold', 0.3))
        top_k = int(data.get('top_k', 10))
        use_lsh = data.get('use_lsh', True)
        multiscale_weighting = data.get('multiscale_weighting', [0.2, 0.6, 0.2])

        # 执行搜索
        system = get_vs_system()
        results = system.enhanced_screening(
            smiles,
            similarity_threshold=threshold,
            top_k=top_k,
            use_lsh=use_lsh,
            multiscale_weighting=multiscale_weighting
        )

        # 格式化结果
        response = {
            'query': {
                'smiles': smiles,
                'mw': results['query_mol_info']['mw'],
                'formula': results['query_mol_info']['formula']
            },
            'parameters': results['screening_params'],
            'performance': results['performance'],
            'statistics': results['statistics'],
            'results': results['results'][:top_k],  # 限制返回数量
            'total_found': len(results['results'])
        }

        return jsonify(response)

    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/api/validate_smiles', methods=['POST'])
def validate_smiles():
    """验证SMILES有效性"""
    try:
        data = request.get_json()
        if not data or 'smiles' not in data:
            return jsonify({'valid': False, 'error': '缺少SMILES参数'})

        smiles = data['smiles'].strip()
        if not smiles:
            return jsonify({'valid': False, 'error': 'SMILES不能为空'})

        # 验证SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'valid': False, 'error': '无效的SMILES字符串'})

        # 计算基本性质
        mw = Chem.rdMolDescriptors.CalcExactMolWt(mol)
        formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
        num_atoms = mol.GetNumAtoms()
        num_bonds = mol.GetNumBonds()

        return jsonify({
            'valid': True,
            'properties': {
                'molecular_weight': round(mw, 2),
                'formula': formula,
                'num_atoms': num_atoms,
                'num_bonds': num_bonds
            }
        })

    except Exception as e:
        return jsonify({'valid': False, 'error': str(e)})

@app.route('/api/library_stats')
def library_stats():
    """获取分子库统计信息"""
    try:
        system = get_vs_system()

        # 读取集成库统计信息
        stats_file = os.path.join(system.data_dir, 'library_statistics.json')
        if os.path.exists(stats_file):
            with open(stats_file, 'r', encoding='utf-8') as f:
                stats = json.load(f)
        else:
            stats = {
                'total_compounds': len(system.molecules),
                'data_sources': {'Enhanced_Dataset': len(system.molecules)},
                'categories': {},
                'molecular_properties': {}
            }

        return jsonify(stats)

    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/api/export_results', methods=['POST'])
def export_results():
    """导出搜索结果"""
    try:
        data = request.get_json()
        if not data or 'results' not in data:
            return jsonify({'error': '缺少结果数据'}), 400

        results = data['results']
        format_type = data.get('format', 'csv')

        # 转换为DataFrame
        df_data = []
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
            df_data.append(row)

        df = pd.DataFrame(df_data)

        # 根据格式导出
        if format_type == 'csv':
            output = io.StringIO()
            df.to_csv(output, index=False, encoding='utf-8-sig')
            output.seek(0)
            return send_file(
                io.BytesIO(output.getvalue().encode('utf-8-sig')),
                mimetype='text/csv',
                as_attachment=True,
                download_name='screening_results.csv'
            )
        elif format_type == 'json':
            return jsonify(df.to_dict('records'))
        else:
            return jsonify({'error': '不支持的导出格式'}), 400

    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/api/generate_plot/<plot_type>')
def generate_plot(plot_type):
    """生成可视化图表"""
    try:
        # 这里可以根据plot_type生成不同的图表
        # 暂时返回一个简单的示例图表

        fig, ax = plt.subplots(figsize=(8, 6))

        if plot_type == 'similarity_distribution':
            # 相似度分布示例
            similarities = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
            ax.hist(similarities, bins=10, alpha=0.7, color='skyblue', edgecolor='black')
            ax.set_xlabel('相似度')
            ax.set_ylabel('频次')
            ax.set_title('相似度分布')
        elif plot_type == 'category_distribution':
            # 类别分布示例
            categories = ['NSAIDs', '抗生素', '抗病毒药', '其他']
            counts = [15, 8, 5, 12]
            ax.bar(categories, counts, alpha=0.7, color='lightcoral')
            ax.set_xlabel('类别')
            ax.set_ylabel('数量')
            ax.set_title('分子类别分布')
            plt.xticks(rotation=45)
        else:
            ax.text(0.5, 0.5, f'图表类型: {plot_type}', ha='center', va='center', transform=ax.transAxes)
            ax.set_title('示例图表')

        plt.tight_layout()

        # 转换为base64编码
        buf = io.BytesIO()
        fig.savefig(buf, format='png', dpi=100, bbox_inches='tight')
        buf.seek(0)
        img_base64 = base64.b64encode(buf.read()).decode('utf-8')
        plt.close(fig)

        return jsonify({'image': f'data:image/png;base64,{img_base64}'})

    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/api/molecule_image/<smiles>')
def molecule_image(smiles):
    """生成分子结构图"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': '无效的SMILES'}), 400

        # 生成分子图像
        img = Draw.MolToImage(mol, size=(300, 300))

        # 转换为base64
        buf = io.BytesIO()
        img.save(buf, format='PNG')
        buf.seek(0)
        img_base64 = base64.b64encode(buf.read()).decode('utf-8')

        return jsonify({'image': f'data:image/png;base64,{img_base64}'})

    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/api/explain_formula', methods=['POST'])
def explain_formula():
    """解释化学式"""
    try:
        data = request.get_json()
        if not data or 'formula' not in data:
            return jsonify({'error': '缺少化学式参数'}), 400

        formula = data['formula'].strip()
        if not formula:
            return jsonify({'error': '化学式不能为空'}), 400

        explanation = chemistry_interpreter.explain_formula(formula)
        return jsonify(explanation)

    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/api/analyze_molecule', methods=['POST'])
def analyze_molecule():
    """分析分子结构和性质"""
    try:
        data = request.get_json()
        if not data or 'smiles' not in data:
            return jsonify({'error': '缺少SMILES参数'}), 400

        smiles = data['smiles'].strip()
        properties = data.get('properties', {})

        insights = chemistry_interpreter.get_molecule_insights(smiles, properties)
        if insights and 'error' not in insights:
            return jsonify(insights)
        else:
            return jsonify({'error': insights.get('error', '分析失败')}), 400

    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/api/element_knowledge')
def element_knowledge():
    """获取元素知识库"""
    try:
        knowledge = chemistry_interpreter.get_element_knowledge()
        return jsonify(knowledge)
    except Exception as e:
        return jsonify({'error': str(e)}), 500

# 创建模板目录和HTML文件
def create_templates():
    """创建HTML模板文件"""
    template_dir = os.path.join(os.path.dirname(__file__), 'templates')
    os.makedirs(template_dir, exist_ok=True)

    # 主页模板
    index_html = """
<!DOCTYPE html>
<html lang="zh-CN">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>分子相似性检索系统</title>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css" rel="stylesheet">
    <style>
        .hero-section {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 100px 0;
        }
        .feature-card {
            transition: transform 0.3s;
        }
        .feature-card:hover {
            transform: translateY(-5px);
        }
    </style>
</head>
<body>
    <nav class="navbar navbar-expand-lg navbar-dark bg-dark">
        <div class="container">
            <a class="navbar-brand" href="#">🧬 分子相似性检索系统</a>
            <div class="navbar-nav ms-auto">
                <a class="nav-link active" href="/">首页</a>
                <a class="nav-link" href="/search">相似性搜索</a>
            </div>
        </div>
    </nav>

    <section class="hero-section">
        <div class="container text-center">
            <h1 class="display-4 mb-4">分子相似性检索系统</h1>
            <p class="lead mb-4">基于多尺度指纹融合和LSH加速的虚拟筛选平台</p>
            <a href="/search" class="btn btn-light btn-lg">开始搜索</a>
        </div>
    </section>

    <section class="py-5">
        <div class="container">
            <div class="row">
                <div class="col-md-4">
                    <div class="card feature-card h-100">
                        <div class="card-body text-center">
                            <h5 class="card-title">🔬 多尺度指纹</h5>
                            <p class="card-text">融合ECFP2、ECFP4、ECFP6多尺度分子指纹，提供更准确的相似性计算</p>
                        </div>
                    </div>
                </div>
                <div class="col-md-4">
                    <div class="card feature-card h-100">
                        <div class="card-body text-center">
                            <h5 class="card-title">⚡ LSH加速</h5>
                            <p class="card-text">基于局部敏感哈希的近似最近邻搜索，大幅提升检索速度</p>
                        </div>
                    </div>
                </div>
                <div class="col-md-4">
                    <div class="card feature-card h-100">
                        <div class="card-body text-center">
                            <h5 class="card-title">📊 智能可视化</h5>
                            <p class="card-text">交互式图表展示，实时结果分析和导出功能</p>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    </section>

    <footer class="bg-dark text-white py-4">
        <div class="container text-center">
            <p>&copy; 2026 分子相似性检索系统. 基于多尺度指纹融合技术.</p>
        </div>
    </footer>

    <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/js/bootstrap.bundle.min.js"></script>
</body>
</html>
    """

    # 搜索页面模板
    search_html = """
<!DOCTYPE html>
<html lang="zh-CN">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>分子相似性搜索</title>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css" rel="stylesheet">
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <style>
        .molecule-image {
            max-width: 200px;
            max-height: 200px;
            border: 1px solid #ddd;
            border-radius: 5px;
        }
        .result-card {
            transition: all 0.3s;
        }
        .result-card:hover {
            box-shadow: 0 4px 8px rgba(0,0,0,0.1);
        }
        .loading {
            display: none;
        }
    </style>
</head>
<body>
    <nav class="navbar navbar-expand-lg navbar-dark bg-dark">
        <div class="container">
            <a class="navbar-brand" href="#">🧬 分子相似性检索系统</a>
            <div class="navbar-nav ms-auto">
                <a class="nav-link" href="/">首页</a>
                <a class="nav-link active" href="/search">相似性搜索</a>
            </div>
        </div>
    </nav>

    <div class="container mt-4">
        <div class="row">
            <div class="col-md-4">
                <div class="card">
                    <div class="card-header">
                        <h5>🔍 分子搜索</h5>
                    </div>
                    <div class="card-body">
                        <form id="searchForm">
                            <div class="mb-3">
                                <label for="smiles" class="form-label">SMILES字符串</label>
                                <input type="text" class="form-control" id="smiles"
                                       placeholder="输入分子SMILES，如: CC(=O)OC1=CC=CC=C1C(=O)O"
                                       required>
                                <div class="form-text" id="smilesHelp"></div>
                            </div>

                            <div class="mb-3">
                                <label for="threshold" class="form-label">相似度阈值</label>
                                <input type="range" class="form-range" id="threshold"
                                       min="0.1" max="1.0" step="0.1" value="0.3">
                                <div class="text-center" id="thresholdValue">0.3</div>
                            </div>

                            <div class="mb-3">
                                <label for="topK" class="form-label">返回结果数</label>
                                <select class="form-select" id="topK">
                                    <option value="5">5</option>
                                    <option value="10" selected>10</option>
                                    <option value="15">15</option>
                                    <option value="20">20</option>
                                </select>
                            </div>

                            <div class="mb-3 form-check">
                                <input type="checkbox" class="form-check-input" id="useLSH" checked>
                                <label class="form-check-label" for="useLSH">使用LSH加速</label>
                            </div>

                            <button type="submit" class="btn btn-primary w-100" id="searchBtn">
                                <span class="spinner-border spinner-border-sm loading" role="status"></span>
                                开始搜索
                            </button>
                        </form>
                    </div>
                </div>

                <div class="card mt-3" id="queryMoleculeCard" style="display: none;">
                    <div class="card-header">
                        <h6>查询分子</h6>
                    </div>
                    <div class="card-body text-center">
                        <img id="queryMoleculeImg" class="molecule-image" alt="查询分子结构">
                        <div class="mt-2">
                            <small id="queryInfo" class="text-muted"></small>
                        </div>
                    </div>
                </div>
            </div>

            <div class="col-md-8">
                <div class="card">
                    <div class="card-header d-flex justify-content-between align-items-center">
                        <h5>📊 搜索结果</h5>
                        <div>
                            <button class="btn btn-sm btn-outline-primary" id="exportBtn" style="display: none;">
                                导出结果
                            </button>
                        </div>
                    </div>
                    <div class="card-body">
                        <div id="resultsContainer">
                            <div class="text-center text-muted">
                                <p>请输入分子SMILES并点击搜索</p>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    </div>

    <script>
        // 页面加载完成后的初始化
        document.addEventListener('DOMContentLoaded', function() {
            initializeEventListeners();
        });

        function initializeEventListeners() {
            // 阈值滑块变化
            document.getElementById('threshold').addEventListener('input', function() {
                document.getElementById('thresholdValue').textContent = this.value;
            });

            // 搜索表单提交
            document.getElementById('searchForm').addEventListener('submit', function(e) {
                e.preventDefault();
                performSearch();
            });

            // SMILES输入验证
            document.getElementById('smiles').addEventListener('input', function() {
                validateSMILES(this.value);
            });

            // 导出按钮
            document.getElementById('exportBtn').addEventListener('click', exportResults);
        }

        async function validateSMILES(smiles) {
            if (!smiles.trim()) {
                document.getElementById('smilesHelp').textContent = '';
                return;
            }

            try {
                const response = await fetch('/api/validate_smiles', {
                    method: 'POST',
                    headers: {
                        'Content-Type': 'application/json',
                    },
                    body: JSON.stringify({ smiles: smiles })
                });

                const result = await response.json();

                if (result.valid) {
                    document.getElementById('smilesHelp').innerHTML =
                        `<span class="text-success">✓ 有效分子 | MW: ${result.properties.molecular_weight} | 公式: ${result.properties.formula}</span>`;
                } else {
                    document.getElementById('smilesHelp').innerHTML =
                        `<span class="text-danger">✗ ${result.error}</span>`;
                }
            } catch (error) {
                console.error('SMILES验证失败:', error);
            }
        }

        async function performSearch() {
            const smiles = document.getElementById('smiles').value.trim();
            if (!smiles) return;

            // 显示加载状态
            const searchBtn = document.getElementById('searchBtn');
            const spinner = searchBtn.querySelector('.loading');
            searchBtn.disabled = true;
            spinner.style.display = 'inline-block';

            try {
                const searchData = {
                    smiles: smiles,
                    threshold: parseFloat(document.getElementById('threshold').value),
                    top_k: parseInt(document.getElementById('topK').value),
                    use_lsh: document.getElementById('useLSH').checked
                };

                const response = await fetch('/api/search', {
                    method: 'POST',
                    headers: {
                        'Content-Type': 'application/json',
                    },
                    body: JSON.stringify(searchData)
                });

                const result = await response.json();

                if (response.ok) {
                    displayResults(result);
                    // 显示查询分子
                    displayQueryMolecule(smiles);
                } else {
                    showError(result.error);
                }

            } catch (error) {
                showError('搜索请求失败: ' + error.message);
            } finally {
                // 隐藏加载状态
                searchBtn.disabled = false;
                spinner.style.display = 'none';
            }
        }

        function displayQueryMolecule(smiles) {
            fetch(`/api/molecule_image/${encodeURIComponent(smiles)}`)
                .then(response => response.json())
                .then(data => {
                    if (data.image) {
                        document.getElementById('queryMoleculeImg').src = data.image;
                        document.getElementById('queryMoleculeCard').style.display = 'block';
                        document.getElementById('queryInfo').textContent = `SMILES: ${smiles}`;
                    }
                })
                .catch(error => console.error('获取分子图像失败:', error));
        }

        function displayResults(result) {
            const container = document.getElementById('resultsContainer');

            if (!result.results || result.results.length === 0) {
                container.innerHTML = '<div class="alert alert-warning">未找到相似分子</div>';
                return;
            }

            let html = `
                <div class="alert alert-success">
                    <strong>搜索完成!</strong> 找到 ${result.total_found} 个相似分子
                    <br>处理时间: ${result.performance.elapsed_time.toFixed(3)}秒
                    吞吐量: ${result.performance.throughput.toFixed(1)} 个/秒
                </div>

                <div class="row">
            `;

            result.results.forEach((item, index) => {
                html += `
                    <div class="col-md-6 mb-3">
                        <div class="card result-card">
                            <div class="card-body">
                                <h6 class="card-title">${index + 1}. ${item.name}</h6>
                                <div class="row">
                                    <div class="col-4">
                                        <img class="molecule-image w-100" src="/api/molecule_image/${encodeURIComponent(item.smiles)}"
                                             alt="${item.name}" onerror="this.style.display='none'">
                                    </div>
                                    <div class="col-8">
                                        <p class="mb-1"><strong>相似度:</strong> ${(item.weighted_similarity * 100).toFixed(1)}%</p>
                                        <p class="mb-1"><small class="text-muted">类别: ${item.category}</small></p>
                                        <p class="mb-1"><small class="text-muted">分子量: ${item.mw.toFixed(1)}</small></p>
                                        <p class="mb-0"><small class="text-muted">SMILES: ${item.smiles.substring(0, 30)}${item.smiles.length > 30 ? '...' : ''}</small></p>
                                    </div>
                                </div>
                            </div>
                        </div>
                    </div>
                `;
            });

            html += '</div>';
            container.innerHTML = html;

            // 显示导出按钮
            document.getElementById('exportBtn').style.display = 'block';
            document.getElementById('exportBtn').onclick = () => exportResults(result.results);
        }

        function showError(message) {
            const container = document.getElementById('resultsContainer');
            container.innerHTML = `<div class="alert alert-danger">${message}</div>`;
        }

        function exportResults(results = null) {
            if (!results) {
                // 从当前显示的结果中获取
                const resultCards = document.querySelectorAll('.result-card');
                results = Array.from(resultCards).map(card => {
                    // 这里需要解析卡片内容，暂时使用示例数据
                    return {};
                });
            }

            // 创建CSV内容
            let csv = '候选分子,加权相似度,类别,分子量,SMILES\\n';
            results.forEach(result => {
                csv += `"${result.name}",${result.weighted_similarity},"${result.category}",${result.mw},"${result.smiles}"\\n`;
            });

            // 下载CSV
            const blob = new Blob([csv], { type: 'text/csv;charset=utf-8;' });
            const link = document.createElement('a');
            link.href = URL.createObjectURL(blob);
            link.download = 'screening_results.csv';
            link.click();
        }
    </script>

    <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/js/bootstrap.bundle.min.js"></script>
</body>
</html>
    """

    # 写入模板文件
    with open(os.path.join(template_dir, 'index.html'), 'w', encoding='utf-8') as f:
        f.write(index_html)

    with open(os.path.join(template_dir, 'search.html'), 'w', encoding='utf-8') as f:
        f.write(search_html)

    print(f"✅ HTML模板已创建: {template_dir}")

def main():
    """主函数 - 启动Web服务器"""
    print("=" * 70)
    print("🌐 第二阶段: Web界面系统")
    print("基于Flask的基础Web应用")
    print("=" * 70)

    # 创建模板文件
    create_templates()

    # 启动服务器
    print("🚀 启动Web服务器...")
    print("📱 访问地址: http://localhost:5000")
    print("🏠 主页: http://localhost:5000/")
    print("🔍 搜索: http://localhost:5000/search")
    print("🛑 按 Ctrl+C 停止服务器")
    print("=" * 70)

    app.run(debug=True, host='0.0.0.0', port=5000)

if __name__ == "__main__":
    main()