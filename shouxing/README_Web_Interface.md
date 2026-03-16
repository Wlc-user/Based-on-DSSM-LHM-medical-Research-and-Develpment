# 分子相似性检索系统 - 第二阶段: Web界面系统

## 📋 项目概述

基于Flask框架构建的分子相似性检索Web应用，提供直观的分子搜索和虚拟筛选界面。

## 🚀 核心功能

### 🔬 分子相似性搜索
- **实时SMILES验证**: 输入分子SMILES字符串，实时验证有效性并显示分子性质
- **多尺度指纹融合**: 结合ECFP2、ECFP4、ECFP6指纹进行相似性计算
- **LSH加速搜索**: 使用局部敏感哈希算法大幅提升搜索速度
- **智能结果排序**: 基于加权相似度进行结果排序和展示

### 🌐 Web界面特性
- **响应式设计**: 适配桌面和移动设备
- **交互式搜索**: 实时参数调整和结果预览
- **分子结构可视化**: 自动生成分子结构图像
- **结果导出**: 支持CSV格式导出搜索结果

### 📊 数据统计
- **分子库概览**: 显示当前分子库的统计信息
- **性能监控**: 实时显示搜索性能指标
- **数据源分析**: 展示不同数据源的分布情况

## 🛠️ 技术架构

```
Web界面系统
├── Flask应用 (web_interface.py)
├── HTML模板 (templates/)
│   ├── index.html      # 主页
│   └── search.html     # 搜索页面
├── RESTful API
│   ├── /api/search          # 分子搜索
│   ├── /api/validate_smiles # SMILES验证
│   ├── /api/library_stats   # 库统计
│   └── /api/export_results  # 结果导出
└── 集成系统
    ├── enhanced_virtual_screening.py  # 核心搜索引擎
    └── integrated_molecule_library.csv # 分子数据库
```

## 📦 安装依赖

```bash
# 安装Web框架依赖
pip install flask flask-cors

# 验证安装
python -c "import flask, flask_cors; print('✅ 依赖安装成功')"
```

## 🚀 启动服务

### 方法1: 直接运行
```bash
cd e:\pyspace\shouxing
python web_interface.py
```

### 方法2: 后台运行
```bash
# Windows PowerShell
Start-Process python -ArgumentList "web_interface.py" -WorkingDirectory "e:\pyspace\shouxing"

# 或使用 nohup (如果安装了)
python web_interface.py &
```

## 🌐 访问界面

启动服务后，在浏览器中访问：

- **主页**: http://localhost:5000/
- **搜索页面**: http://localhost:5000/search

## 📖 使用指南

### 1. 分子搜索
1. 在搜索页面输入分子SMILES字符串
2. 调整相似度阈值 (0.1-1.0)
3. 设置返回结果数量 (5-20)
4. 选择是否使用LSH加速
5. 点击"开始搜索"

### 2. SMILES验证
- 输入SMILES时会实时验证有效性
- 显示分子量、化学式、原子数等性质
- 无效SMILES会显示错误提示

### 3. 结果查看
- 搜索结果按相似度降序排列
- 显示分子结构图像、相似度分数、类别等信息
- 支持结果导出为CSV文件

## 🔧 API接口

### POST /api/search
执行分子相似性搜索

**请求参数:**
```json
{
  "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
  "threshold": 0.3,
  "top_k": 10,
  "use_lsh": true,
  "multiscale_weighting": [0.2, 0.6, 0.2]
}
```

**响应格式:**
```json
{
  "query": {
    "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "mw": 180.04,
    "formula": "C9H8O4"
  },
  "parameters": {...},
  "performance": {
    "elapsed_time": 0.015,
    "throughput": 66.7
  },
  "results": [...]
}
```

### POST /api/validate_smiles
验证SMILES字符串有效性

### GET /api/library_stats
获取分子库统计信息

### POST /api/export_results
导出搜索结果

## 📈 性能指标

基于当前测试结果：

- **搜索速度**: ~100个/秒 (使用LSH加速)
- **响应时间**: <50ms (典型查询)
- **分子库规模**: 60个化合物
- **支持格式**: SMILES字符串
- **相似度算法**: 多尺度Morgan指纹 + Tanimoto系数

## 🧪 测试验证

运行功能演示脚本：

```bash
python web_demo.py
```

该脚本会自动测试所有API接口并验证功能完整性。

## 🔄 集成说明

Web界面系统完全集成第一阶段的核心功能：

- ✅ 多尺度指纹融合 (ECFP2/4/6)
- ✅ LSH近似最近邻搜索
- ✅ 增强数据集 (PubChem + ZINC)
- ✅ 性能基准测试
- ✅ 结果可视化

## 🎯 下一步计划

第三阶段开发重点：

1. **高级Web功能**
   - 用户会话管理
   - 批量上传处理
   - 高级筛选条件

2. **API扩展**
   - RESTful API完善
   - GraphQL接口
   - 第三方集成

3. **生产部署**
   - Docker容器化
   - 云服务部署
   - 负载均衡

## 📝 注意事项

1. **端口占用**: 默认使用5000端口，确保该端口未被占用
2. **内存使用**: 大型分子库可能需要较多内存
3. **网络安全**: 生产环境建议配置HTTPS和认证
4. **数据更新**: 分子库可通过数据集成脚本定期更新

## 🆘 故障排除

### 常见问题

1. **端口被占用**
   ```bash
   # 检查端口占用
   netstat -ano | findstr :5000
   # 修改端口
   python web_interface.py  # 在代码中修改port参数
   ```

2. **依赖缺失**
   ```bash
   pip install -r requirements.txt
   ```

3. **浏览器缓存**
   - 硬刷新页面 (Ctrl+F5)
   - 清除浏览器缓存

## 📄 许可证

本项目采用MIT许可证。详见LICENSE文件。

## 🤝 贡献

欢迎提交Issue和Pull Request来改进系统功能。

---

**开发时间**: 2026年3月
**技术栈**: Python + Flask + RDKit + Bootstrap
**系统状态**: ✅ 第二阶段完成，Web界面运行正常