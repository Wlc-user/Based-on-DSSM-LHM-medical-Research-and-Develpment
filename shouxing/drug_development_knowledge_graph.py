#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
医药研发知识图谱系统
构建药物研发的核心关联网络，包含毒理学分析和前后置关系推理

知识图谱包含:
- 药物研发流程图谱
- 毒理学分析网络
- 临床试验关联
- 药物-靶点-疾病关联
- ADMET性质预测
- 药物相互作用网络
"""

import json
import networkx as nx
from collections import defaultdict
from chemistry_knowledge_graph import ChemistryKnowledgeGraph
from chemistry_interpreter import ChemicalFormulaInterpreter

class DrugDevelopmentKnowledgeGraph:
    """医药研发知识图谱"""

    def __init__(self):
        self.graph = nx.DiGraph()
        self.chemistry_kg = ChemistryKnowledgeGraph()
        self.interpreter = ChemicalFormulaInterpreter()
        self._build_drug_development_graph()

    def _build_drug_development_graph(self):
        """构建医药研发知识图谱"""
        print("🏥 构建医药研发知识图谱...")

        # 1. 药物研发流程
        self._build_drug_discovery_pipeline()

        # 2. 毒理学分析网络
        self._build_toxicology_network()

        # 3. 临床试验关联
        self._build_clinical_trial_network()

        # 4. 药物-靶点-疾病关联
        self._build_drug_target_disease_network()

        # 5. ADMET性质预测
        self._build_admet_network()

        # 6. 药物相互作用
        self._build_drug_interaction_network()

        print(f"✅ 医药研发知识图谱构建完成: {len(self.graph.nodes)} 个节点, {len(self.graph.edges)} 条边")

    def _build_drug_discovery_pipeline(self):
        """构建药物研发流程图谱"""
        print("   🔬 构建药物研发流程图谱...")

        # 研发阶段
        stages = {
            'hit_identification': {
                'name': '先导化合物识别', 'phase': 1,
                'description': '通过虚拟筛选、高通量筛选等方法发现活性化合物',
                'success_rate': 0.001, 'time': '3-6个月', 'cost': '50-200万美元'
            },
            'hit_to_lead': {
                'name': '先导化合物优化', 'phase': 2,
                'description': '优化先导化合物的药效、选择性和ADMET性质',
                'success_rate': 0.1, 'time': '6-12个月', 'cost': '200-500万美元'
            },
            'lead_optimization': {
                'name': '候选药物优化', 'phase': 3,
                'description': '进一步优化药代动力学和安全性',
                'success_rate': 0.3, 'time': '12-24个月', 'cost': '500-2000万美元'
            },
            'preclinical': {
                'name': '临床前研究', 'phase': 4,
                'description': '体内外药效学、毒理学和药代动力学研究',
                'success_rate': 0.7, 'time': '12-24个月', 'cost': '1000-5000万美元'
            },
            'phase1': {
                'name': 'I期临床试验', 'phase': 5,
                'description': '健康志愿者安全性、耐受性和药代动力学研究',
                'success_rate': 0.6, 'time': '6-12个月', 'cost': '1000-2000万美元'
            },
            'phase2': {
                'name': 'II期临床试验', 'phase': 6,
                'description': '患者疗效和安全性评估',
                'success_rate': 0.3, 'time': '12-24个月', 'cost': '2000-5000万美元'
            },
            'phase3': {
                'name': 'III期临床试验', 'phase': 7,
                'description': '大规模疗效和安全性验证',
                'success_rate': 0.5, 'time': '24-36个月', 'cost': '5000-2亿美元'
            },
            'nda_review': {
                'name': '新药申请审评', 'phase': 8,
                'description': '监管机构审评和批准',
                'success_rate': 0.8, 'time': '6-12个月', 'cost': '500-1000万美元'
            },
            'post_marketing': {
                'name': '上市后监测', 'phase': 9,
                'description': '上市后安全性监测和真实世界证据收集',
                'success_rate': 1.0, 'time': '持续', 'cost': '每年500-2000万美元'
            }
        }

        # 添加阶段节点
        for stage_id, data in stages.items():
            self.graph.add_node(f"stage_{stage_id}",
                              type="development_stage",
                              stage_id=stage_id,
                              **data)

        # 添加前后置关系
        stage_sequence = ['hit_identification', 'hit_to_lead', 'lead_optimization',
                         'preclinical', 'phase1', 'phase2', 'phase3',
                         'nda_review', 'post_marketing']

        for i in range(len(stage_sequence) - 1):
            current = f"stage_{stage_sequence[i]}"
            next_stage = f"stage_{stage_sequence[i+1]}"
            self.graph.add_edge(current, next_stage,
                              relation="leads_to",
                              description="研发流程顺序")

        # 添加关键决策点
        decision_points = {
            'go_no_go_1': {'name': 'GO/NO-GO决策点1', 'stage': 'hit_to_lead', 'criteria': ['活性', '选择性', '专利性']},
            'go_no_go_2': {'name': 'GO/NO-GO决策点2', 'stage': 'lead_optimization', 'criteria': ['ADMET', '安全性', '药效']},
            'go_no_go_3': {'name': 'GO/NO-GO决策点3', 'stage': 'preclinical', 'criteria': ['毒理学', '药代动力学', '临床潜力']},
            'go_no_go_4': {'name': 'GO/NO-GO决策点4', 'stage': 'phase2', 'criteria': ['疗效', '安全性', '商业潜力']}
        }

        for dp_id, data in decision_points.items():
            self.graph.add_node(f"decision_{dp_id}",
                              type="decision_point",
                              dp_id=dp_id,
                              **data)
            # 连接到对应阶段
            self.graph.add_edge(f"stage_{data['stage']}", f"decision_{dp_id}",
                              relation="has_decision")

    def _build_toxicology_network(self):
        """构建毒理学分析网络"""
        print("   ☠️ 构建毒理学分析网络...")

        # 毒理学评估类型
        toxicology_types = {
            'genotoxicity': {
                'name': '遗传毒性', 'description': 'DNA损伤、染色体畸变、基因突变',
                'assays': ['Ames试验', '染色体畸变试验', '小鼠淋巴瘤试验'],
                'regulatory': 'ICH S2(R1)', 'critical': True
            },
            'carcinogenicity': {
                'name': '致癌性', 'description': '长期暴露下的肿瘤形成风险',
                'assays': ['转基因小鼠模型', '大鼠长期毒性试验'],
                'regulatory': 'ICH S1', 'critical': True
            },
            'reproductive_toxicity': {
                'name': '生殖毒性', 'description': '对生育能力和后代的影响',
                'assays': ['生殖毒性试验', '胚胎-胎仔毒性试验'],
                'regulatory': 'ICH S5(R3)', 'critical': True
            },
            'developmental_toxicity': {
                'name': '发育毒性', 'description': '对胚胎和胎仔发育的影响',
                'assays': ['胚胎发育毒性试验'],
                'regulatory': 'ICH S5(R3)', 'critical': True
            },
            'cardiovascular_toxicity': {
                'name': '心血管毒性', 'description': '对心脏和血管系统的损害',
                'assays': ['心电图监测', '血压测量', '心脏组织病理'],
                'regulatory': 'ICH S7A', 'critical': True
            },
            'hepatotoxicity': {
                'name': '肝毒性', 'description': '对肝脏的损害',
                'assays': ['肝酶检测', '肝组织病理', '胆汁淤积评估'],
                'regulatory': 'ICH S3A', 'critical': True
            },
            'nephrotoxicity': {
                'name': '肾毒性', 'description': '对肾脏的损害',
                'assays': ['肾功能检测', '肾组织病理', '尿液分析'],
                'regulatory': 'ICH S3A', 'critical': True
            },
            'neurotoxicity': {
                'name': '神经毒性', 'description': '对神经系统的损害',
                'assays': ['行为测试', '神经病理检查', '神经递质分析'],
                'regulatory': 'ICH S7', 'critical': True
            },
            'immunotoxicity': {
                'name': '免疫毒性', 'description': '对免疫系统的损害',
                'assays': ['免疫器官病理', '免疫细胞计数', '抗体产生评估'],
                'regulatory': 'ICH S8', 'critical': True
            },
            'phototoxicity': {
                'name': '光毒性', 'description': '光照诱导的毒性反应',
                'assays': ['3T3中性红摄取试验', '体内光毒性试验'],
                'regulatory': 'ICH S10', 'critical': False
            }
        }

        # 添加毒理学节点
        for tox_id, data in toxicology_types.items():
            self.graph.add_node(f"tox_{tox_id}",
                              type="toxicology",
                              tox_id=tox_id,
                              **data)

        # 毒理学评估策略
        strategies = {
            'in_silico': {'name': '计算机预测', 'methods': ['QSAR', '分子对接', '机器学习'], 'cost': '低', 'time': '快'},
            'in_vitro': {'name': '体外试验', 'methods': ['细胞毒性', '酶抑制', '转基因细胞'], 'cost': '中', 'time': '中'},
            'in_vivo': {'name': '体内试验', 'methods': ['急性毒性', '亚慢性毒性', '慢性毒性'], 'cost': '高', 'time': '慢'},
            'clinical_monitoring': {'name': '临床监测', 'methods': ['不良反应报告', '生物标志物'], 'cost': '中', 'time': '持续'}
        }

        for strat_id, data in strategies.items():
            self.graph.add_node(f"strategy_{strat_id}",
                              type="toxicology_strategy",
                              strat_id=strat_id,
                              **data)

        # 连接毒理学类型到评估策略
        tox_strategy_relations = [
            ('genotoxicity', 'in_silico', 'primary'),
            ('genotoxicity', 'in_vitro', 'primary'),
            ('carcinogenicity', 'in_vivo', 'primary'),
            ('reproductive_toxicity', 'in_vivo', 'primary'),
            ('developmental_toxicity', 'in_vivo', 'primary'),
            ('cardiovascular_toxicity', 'in_vitro', 'primary'),
            ('cardiovascular_toxicity', 'in_vivo', 'secondary'),
            ('hepatotoxicity', 'in_vitro', 'primary'),
            ('hepatotoxicity', 'in_vivo', 'secondary'),
            ('nephrotoxicity', 'in_vitro', 'primary'),
            ('nephrotoxicity', 'in_vivo', 'secondary'),
            ('neurotoxicity', 'in_vivo', 'primary'),
            ('immunotoxicity', 'in_vivo', 'primary'),
            ('phototoxicity', 'in_vitro', 'primary')
        ]

        for tox, strat, priority in tox_strategy_relations:
            self.graph.add_edge(f"tox_{tox}", f"strategy_{strat}",
                              relation="assessed_by",
                              priority=priority)

    def _build_clinical_trial_network(self):
        """构建临床试验关联网络"""
        print("   🏥 构建临床试验关联网络...")

        # 临床试验设计要素
        trial_elements = {
            'study_design': {
                'name': '试验设计', 'category': 'design',
                'elements': ['随机化', '盲法', '对照组', '样本量计算']
            },
            'endpoints': {
                'name': '终点指标', 'category': 'outcome',
                'elements': ['主要终点', '次要终点', '安全性终点', '探索性终点']
            },
            'patient_population': {
                'name': '受试者人群', 'category': 'population',
                'elements': ['纳入标准', '排除标准', '分层因素', '代表性']
            },
            'statistical_analysis': {
                'name': '统计分析', 'category': 'analysis',
                'elements': ['假设检验', '功效分析', '多重性校正', '亚组分析']
            },
            'regulatory_compliance': {
                'name': '法规合规', 'category': 'regulatory',
                'elements': ['伦理审查', '知情同意', 'GCP合规', '数据管理']
            }
        }

        # 添加临床试验节点
        for elem_id, data in trial_elements.items():
            self.graph.add_node(f"trial_{elem_id}",
                              type="clinical_trial_element",
                              elem_id=elem_id,
                              **data)

        # 试验失败原因
        failure_reasons = {
            'efficacy_failure': {'name': '疗效不足', 'frequency': 0.4, 'phase': 'II/III'},
            'safety_failure': {'name': '安全性问题', 'frequency': 0.3, 'phase': 'I-III'},
            'commercial_failure': {'name': '商业化失败', 'frequency': 0.2, 'phase': 'III'},
            'regulatory_failure': {'name': '监管失败', 'frequency': 0.1, 'phase': 'III'}
        }

        for reason_id, data in failure_reasons.items():
            self.graph.add_node(f"failure_{reason_id}",
                              type="trial_failure_reason",
                              reason_id=reason_id,
                              **data)

    def _build_drug_target_disease_network(self):
        """构建药物-靶点-疾病关联网络"""
        print("   🎯 构建药物-靶点-疾病关联网络...")

        # 重要靶点
        targets = {
            'COX1': {'name': '环氧合酶1', 'type': '酶', 'family': '环氧合酶', 'diseases': ['炎症', '疼痛']},
            'COX2': {'name': '环氧合酶2', 'type': '酶', 'family': '环氧合酶', 'diseases': ['炎症', '疼痛', '癌症']},
            'A2A': {'name': '腺苷A2A受体', 'type': 'GPCR', 'family': '嘌呤受体', 'diseases': ['帕金森病', '睡眠障碍']},
            'EGFR': {'name': '表皮生长因子受体', 'type': 'RTK', 'family': '受体酪氨酸激酶', 'diseases': ['肺癌', '乳腺癌']},
            'VEGFR': {'name': '血管内皮生长因子受体', 'type': 'RTK', 'family': '受体酪氨酸激酶', 'diseases': ['癌症', 'AMD']},
            'PPARγ': {'name': '过氧化物酶体增殖物激活受体γ', 'type': '核受体', 'family': '核受体', 'diseases': ['糖尿病', '代谢综合征']},
            'HMGCS': {'name': 'HMG-CoA合成酶', 'type': '酶', 'family': '胆固醇合成酶', 'diseases': ['高胆固醇血症']},
            'ACE': {'name': '血管紧张素转换酶', 'type': '酶', 'family': '金属蛋白酶', 'diseases': ['高血压', '心力衰竭']}
        }

        # 添加靶点节点
        for target_id, data in targets.items():
            self.graph.add_node(f"target_{target_id}",
                              target_id=target_id,
                              **data)

        # 疾病节点
        diseases = {
            'inflammation': {'name': '炎症', 'category': '免疫系统', 'prevalence': '高'},
            'pain': {'name': '疼痛', 'category': '神经系统', 'prevalence': '高'},
            'cancer': {'name': '癌症', 'category': '肿瘤', 'prevalence': '中'},
            'diabetes': {'name': '糖尿病', 'category': '代谢疾病', 'prevalence': '高'},
            'hypertension': {'name': '高血压', 'category': '心血管疾病', 'prevalence': '高'},
            'parkinsons': {'name': '帕金森病', 'category': '神经退行性疾病', 'prevalence': '低'},
            'amd': {'name': '年龄相关性黄斑变性', 'category': '眼科疾病', 'prevalence': '中'}
        }

        for disease_id, data in diseases.items():
            self.graph.add_node(f"disease_{disease_id}",
                              disease_id=disease_id,
                              **data)

        # 连接靶点到疾病
        target_disease_relations = [
            ('COX1', 'inflammation', '治疗'),
            ('COX1', 'pain', '治疗'),
            ('COX2', 'inflammation', '治疗'),
            ('COX2', 'pain', '治疗'),
            ('COX2', 'cancer', '治疗'),
            ('A2A', 'parkinsons', '治疗'),
            ('EGFR', 'cancer', '治疗'),
            ('VEGFR', 'cancer', '治疗'),
            ('VEGFR', 'amd', '治疗'),
            ('PPARγ', 'diabetes', '治疗'),
            ('HMGCS', 'hypertension', '预防'),
            ('ACE', 'hypertension', '治疗')
        ]

        for target, disease, relation_type in target_disease_relations:
            self.graph.add_edge(f"target_{target}", f"disease_{disease}",
                              relation=relation_type,
                              type="target_disease")

    def _build_admet_network(self):
        """构建ADMET性质预测网络"""
        print("   💊 构建ADMET性质预测网络...")

        # ADMET性质
        admet_properties = {
            'absorption': {
                'name': '吸收', 'abbrev': 'A',
                'description': '药物从给药部位进入体循环的能力',
                'key_factors': ['溶解度', '通透性', '稳定性'],
                'prediction_methods': ['Caco-2', 'PAMPA', '计算模型']
            },
            'distribution': {
                'name': '分布', 'abbrev': 'D',
                'description': '药物在体内各组织和器官的分布',
                'key_factors': ['血浆蛋白结合', '组织亲和性', '血脑屏障通透'],
                'prediction_methods': ['血浆蛋白结合率', '组织分布模型']
            },
            'metabolism': {
                'name': '代谢', 'abbrev': 'M',
                'description': '药物在体内的生物转化过程',
                'key_factors': ['CYP酶', 'UGT酶', '其他代谢酶'],
                'prediction_methods': ['体外酶抑制', '微粒体代谢', '计算模型']
            },
            'excretion': {
                'name': '排泄', 'abbrev': 'E',
                'description': '药物及其代谢物从体内的清除',
                'key_factors': ['肾清除', '胆汁排泄', '半衰期'],
                'prediction_methods': ['肾清除模型', '胆汁排泄研究']
            },
            'toxicity': {
                'name': '毒性', 'abbrev': 'T',
                'description': '药物潜在的毒性反应',
                'key_factors': ['剂量依赖性', '靶器官毒性', '免疫毒性'],
                'prediction_methods': ['体外毒性试验', '体内毒性试验', '计算毒性预测']
            }
        }

        # 添加ADMET节点
        for prop_id, data in admet_properties.items():
            self.graph.add_node(f"admet_{prop_id}",
                              prop_id=prop_id,
                              **data)

        # ADMET前后置关系
        admet_sequence = ['absorption', 'distribution', 'metabolism', 'excretion', 'toxicity']
        for i in range(len(admet_sequence) - 1):
            current = f"admet_{admet_sequence[i]}"
            next_prop = f"admet_{admet_sequence[i+1]}"
            self.graph.add_edge(current, next_prop,
                              relation="influences",
                              description="ADMET性质间的相互影响")

    def _build_drug_interaction_network(self):
        """构建药物相互作用网络"""
        print("   ⚠️ 构建药物相互作用网络...")

        # 药物相互作用类型
        interaction_types = {
            'pharmacokinetic': {
                'name': '药代动力学相互作用',
                'subtypes': ['CYP抑制', 'CYP诱导', '转运体抑制', '蛋白结合位移'],
                'severity': '中-高'
            },
            'pharmacodynamic': {
                'name': '药效学相互作用',
                'subtypes': ['协同作用', '拮抗作用', '相加作用'],
                'severity': '中'
            },
            'toxicity_interaction': {
                'name': '毒性相互作用',
                'subtypes': ['协同毒性', '肾毒性增强', '肝毒性增强'],
                'severity': '高'
            }
        }

        # 添加相互作用节点
        for int_id, data in interaction_types.items():
            self.graph.add_node(f"interaction_{int_id}",
                              int_id=int_id,
                              **data)

    def analyze_drug_development_path(self, molecule_smiles, target_disease=None):
        """分析药物研发路径"""
        try:
            # 基本分子分析
            mol_analysis = self.chemistry_kg.explain_molecule_context(molecule_smiles)

            # 毒理学风险评估
            tox_risks = self._assess_toxicology_risks(mol_analysis)

            # ADMET预测
            admet_profile = self._predict_admet_properties(mol_analysis)

            # 研发路径推荐
            development_path = self._recommend_development_path(mol_analysis, tox_risks, admet_profile)

            # 如果指定疾病，添加靶点分析
            target_analysis = None
            if target_disease:
                target_analysis = self._analyze_target_suitability(mol_analysis, target_disease)

            return {
                'molecular_analysis': mol_analysis,
                'toxicology_risks': tox_risks,
                'admet_profile': admet_profile,
                'development_path': development_path,
                'target_analysis': target_analysis
            }

        except Exception as e:
            return {'error': str(e)}

    def _assess_toxicology_risks(self, mol_analysis):
        """评估毒理学风险"""
        risks = {}

        # 基于结构特征的毒性预测
        functional_groups = [fg['type'] for fg in mol_analysis.get('functional_groups', [])]

        # 芳香胺类化合物 - 遗传毒性风险
        if 'amine' in functional_groups and 'aromatic' in functional_groups:
            risks['genotoxicity'] = {
                'level': '高',
                'reason': '芳香胺类化合物可能具有遗传毒性',
                'recommendations': ['进行Ames试验', '评估DNA反应性']
            }

        # 含有多个卤素 - 潜在毒性
        halogen_count = sum(1 for fg in mol_analysis.get('functional_groups', [])
                           if fg['type'] == 'halogen')
        if halogen_count > 2:
            risks['organotoxicity'] = {
                'level': '中',
                'reason': '多卤代化合物可能具有器官毒性',
                'recommendations': ['肝肾功能监测', '长期毒性试验']
            }

        # 分子量过大 - 可能影响代谢和排泄
        mw = mol_analysis.get('basic_info', {}).get('molecular_weight', 0)
        if mw > 500:
            risks['metabolism'] = {
                'level': '中',
                'reason': '分子量较大，可能影响代谢和排泄',
                'recommendations': ['代谢稳定性评估', '肾清除率监测']
            }

        return risks

    def _predict_admet_properties(self, mol_analysis):
        """预测ADMET性质"""
        profile = {}

        # 基于分子性质的简单预测
        mw = mol_analysis.get('basic_info', {}).get('molecular_weight', 0)
        functional_groups = mol_analysis.get('functional_groups', [])

        # 吸收预测
        profile['absorption'] = {
            'prediction': '中' if 200 < mw < 500 else '低' if mw >= 500 else '高',
            'factors': [f'分子量: {mw} Da']
        }

        # 分布预测
        profile['distribution'] = {
            'prediction': '中',
            'factors': ['血浆蛋白结合预测需要实验数据']
        }

        # 代谢预测
        profile['metabolism'] = {
            'prediction': '中',
            'factors': ['CYP酶相互作用预测需要实验数据']
        }

        # 排泄预测
        profile['excretion'] = {
            'prediction': '肾排泄为主' if mw < 300 else '肝胆排泄为主',
            'factors': [f'分子量驱动的排泄途径: {mw} Da']
        }

        return profile

    def _recommend_development_path(self, mol_analysis, tox_risks, admet_profile):
        """推荐研发路径"""
        path = {
            'recommended_starting_point': 'hit_to_lead',
            'critical_checkpoints': [],
            'timeline_estimate': '24-36个月',
            'estimated_cost': '2000-5000万美元'
        }

        # 基于风险调整路径
        if tox_risks:
            path['critical_checkpoints'].append({
                'stage': 'preclinical',
                'focus': '毒理学评估',
                'reason': f'发现{len(tox_risks)}个毒性风险因素'
            })

        # 基于ADMET调整
        absorption_pred = admet_profile.get('absorption', {}).get('prediction', '中')
        if absorption_pred == '低':
            path['critical_checkpoints'].append({
                'stage': 'lead_optimization',
                'focus': '吸收改善',
                'reason': '预测吸收性较差'
            })

        return path

    def _analyze_target_suitability(self, mol_analysis, disease):
        """分析靶点适用性"""
        # 查找相关靶点
        related_targets = []
        for node, data in self.graph.nodes(data=True):
            if data.get('type') == 'drug_target':
                if disease in data.get('diseases', []):
                    related_targets.append({
                        'target': data['target_id'],
                        'name': data['name'],
                        'type': data['type']
                    })

        return {
            'disease': disease,
            'potential_targets': related_targets,
            'suitability_score': len(related_targets) / 10  # 简单评分
        }

    def export_drug_development_graph(self, filename='drug_development_knowledge_graph.json'):
        """导出医药研发知识图谱"""
        graph_data = {
            'nodes': [],
            'edges': [],
            'metadata': {
                'total_nodes': len(self.graph.nodes),
                'total_edges': len(self.graph.edges),
                'node_types': {},
                'created_date': '2026-03-16'
            }
        }

        # 统计节点类型
        for node, data in self.graph.nodes(data=True):
            node_type = data.get('type', 'unknown')
            graph_data['metadata']['node_types'][node_type] = \
                graph_data['metadata']['node_types'].get(node_type, 0) + 1

        # 导出节点
        for node, data in self.graph.nodes(data=True):
            graph_data['nodes'].append({
                'id': node,
                'type': data.get('type'),
                'data': data
            })

        # 导出边
        for source, target, data in self.graph.edges(data=True):
            graph_data['edges'].append({
                'source': source,
                'target': target,
                'relation': data.get('relation'),
                'data': data
            })

        with open(filename, 'w', encoding='utf-8') as f:
            json.dump(graph_data, f, ensure_ascii=False, indent=2)

        print(f"✅ 医药研发知识图谱已导出到: {filename}")

def create_drug_development_analyzer():
    """创建医药研发分析器"""
    analyzer = DrugDevelopmentKnowledgeGraph()

    # 导出知识图谱
    analyzer.export_drug_development_graph()

    return analyzer

if __name__ == "__main__":
    # 创建医药研发知识图谱
    print("🏥 构建医药研发知识图谱系统")
    print("=" * 50)

    analyzer = create_drug_development_analyzer()

    # 测试药物研发路径分析
    print("\n🔬 测试药物研发路径分析")
    print("-" * 40)

    test_drugs = [
        ("阿司匹林", "CC(=O)OC1=CC=CC=C1C(=O)O", "inflammation"),
        ("咖啡因", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", None),
        ("布洛芬", "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", "pain")
    ]

    for name, smiles, disease in test_drugs:
        print(f"\n分析药物: {name}")
        if disease:
            print(f"目标疾病: {disease}")

        analysis = analyzer.analyze_drug_development_path(smiles, disease)

        if 'error' not in analysis:
            print(f"  毒理学风险: {len(analysis['toxicology_risks'])} 项")
            print(f"  ADMET预测: {len(analysis['admet_profile'])} 项")
            print(f"  研发路径: {analysis['development_path']['timeline_estimate']}")

            if analysis['target_analysis']:
                targets = analysis['target_analysis']['potential_targets']
                print(f"  潜在靶点: {len(targets)} 个")
                if targets:
                    print(f"    示例: {targets[0]['name']}")
        else:
            print(f"  分析失败: {analysis['error']}")

    print("\n✅ 医药研发知识图谱系统构建完成!")
    print("📁 输出文件:")
    print("   • drug_development_knowledge_graph.json - 医药研发知识图谱数据")