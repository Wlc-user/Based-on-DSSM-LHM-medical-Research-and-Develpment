import requests
import pubchempy as pcp
from chembl_webresource_client.new_client import new_client
import pandas as pd
import os
from rdkit import Chem
from rdkit.Chem import AllChem

class DataCrawler:
    def __init__(self, data_dir='data'):
        self.data_dir = data_dir
        os.makedirs(data_dir, exist_ok=True)

import requests
import pubchempy as pcp
from chembl_webresource_client.new_client import new_client
import pandas as pd
import os
from rdkit import Chem
from rdkit.Chem import AllChem

class DataCrawler:
    def __init__(self, data_dir='data'):
        self.data_dir = data_dir
        os.makedirs(data_dir, exist_ok=True)

    def crawl_pubchem_bulk(self, limit=500):
        """批量爬取PubChem数据"""
        print(f"批量爬取PubChem数据: {limit} 个分子")
        # 获取一些常见的化合物类别
        queries = ['organic compounds', 'drug', 'natural product', 'small molecule']
        all_compounds = []
        
        for query in queries:
            try:
                compounds = pcp.get_compounds(query, 'name', list_return='flat', limit=limit//len(queries))
                all_compounds.extend(compounds)
                print(f"查询 '{query}' 获得 {len(compounds)} 个化合物")
            except Exception as e:
                print(f"查询 '{query}' 失败: {e}")
                continue
        
        # 去重
        seen_cids = set()
        unique_compounds = []
        for comp in all_compounds:
            if comp.cid and comp.cid not in seen_cids:
                seen_cids.add(comp.cid)
                unique_compounds.append(comp)
        
        data = []
        for comp in unique_compounds[:limit]:
            if comp.cid:
                # 使用新的API属性
                smiles = getattr(comp, 'isomeric_smiles', None) or getattr(comp, 'canonical_smiles', None)
                if not smiles:
                    continue
                    
                mol = Chem.MolFromSmiles(smiles)
                if not mol:
                    continue
                    
                chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
                data.append({
                    'cid': comp.cid,
                    'name': comp.iupac_name or comp.synonyms[0] if comp.synonyms else '',
                    'smiles': smiles,
                    'formula': comp.molecular_formula,
                    'mw': comp.molecular_weight,
                    'chiral_centers': len(chiral_centers),
                    'source': 'PubChem'
                })
        
        df = pd.DataFrame(data)
        df.to_csv(os.path.join(self.data_dir, 'pubchem_bulk_data.csv'), index=False)
        print(f"保存了 {len(df)} 个分子数据")
        return df

    def crawl_chembl(self, query, limit=100):
        """爬取ChEMBL数据"""
        print(f"爬取ChEMBL数据: {query}")
        molecule = new_client.molecule
        results = molecule.filter(pref_name__icontains=query)[:limit]
        data = []
        for mol in results:
            smiles = mol['molecule_structures']['canonical_smiles'] if mol.get('molecule_structures') else None
            rdkit_mol = Chem.MolFromSmiles(smiles) if smiles else None
            chiral_centers = Chem.FindMolChiralCenters(rdkit_mol, includeUnassigned=True) if rdkit_mol else []
            data.append({
                'chembl_id': mol['molecule_chembl_id'],
                'name': mol['pref_name'],
                'smiles': smiles,
                'formula': mol.get('full_molformula'),
                'mw': mol.get('full_mwt'),
                'chiral_centers': len(chiral_centers),
                'source': 'ChEMBL'
            })
        df = pd.DataFrame(data)
        df.to_csv(os.path.join(self.data_dir, 'chembl_data.csv'), index=False)
        return df

    def combine_data(self):
        """合并PubChem和ChEMBL数据"""
        pubchem_file = os.path.join(self.data_dir, 'pubchem_data.csv')
        chembl_file = os.path.join(self.data_dir, 'chembl_data.csv')
        
        combined_data = []
        if os.path.exists(pubchem_file):
            combined_data.append(pd.read_csv(pubchem_file))
        if os.path.exists(chembl_file):
            combined_data.append(pd.read_csv(chembl_file))
        
        if combined_data:
            df = pd.concat(combined_data, ignore_index=True)
            df.to_csv(os.path.join(self.data_dir, 'combined_data.csv'), index=False)
            return df
        return pd.DataFrame()