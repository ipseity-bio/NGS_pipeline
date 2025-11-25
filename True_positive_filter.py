import os
import pandas as pd

input_dir = "./" 
output_dir = "./"


os.makedirs(output_dir, exist_ok=True)

for file_name in os.listdir(input_dir):
    if file_name.endswith("coverage_summary.csv"):
        coverage_file = os.path.join(input_dir, file_name)

data_1 = pd.read_csv(coverage_file, sep=",")

def gt_to_zyg(gt):
    if pd.isna(gt):
        return "-"
    gt = str(gt).strip()

    if gt in ["./.", ".|."]:
        return "-"
    gt = gt.replace("|", "/")

    alleles = gt.split("/")
    if len(alleles) != 2 or "." in alleles:
        return "-"

    a, b = alleles

    # homozygous reference
    if a == "0" and b == "0":
        return "HOM_REF"
    # homozygous alternate 
    if a == b and a != "0":
        return "HOM_ALT"
    # heterozygous
    if a != b:
        return "HET"

    return "-"

# Process each CSV file in the directory
for file_name in os.listdir(input_dir):
    if file_name.endswith("raw_output.csv"):
        input_file = os.path.join(input_dir, file_name)
        base_name, _ = os.path.splitext(file_name) 
        output_file = os.path.join(output_dir, f"{base_name}_filtered_tp.csv") 
        report_output_file = os.path.join(output_dir, f"{base_name}_filtered_tp_report.csv")
        
        data = pd.read_csv(input_file, sep=",")
        
        zyg_missing = data["ZYG"].astype(str).str.strip().isin(["-", ".", "nan", "NaN", ""])
        data.loc[zyg_missing, "ZYG"] = data.loc[zyg_missing, "Genotype (GT)"].apply(gt_to_zyg)
        
        data['MAX_AF (Population)'] = pd.to_numeric(data['MAX_AF (Population)'].replace('-', None), errors='coerce')
        
        data['Avg Read Depth'] = data['Depth of Coverage (DP)']
        data['Gene'] = data['GeneSymbol']
        data['Transcript'] = data['HGVSc'].str.extract(r'^(NM_\d+\.\d+)')
        data['Variant'] = data['HGVSc'].str.extract(r'(c\.\S+)')[0] + data['HGVSp'].str.extract(r'(p\.\S+)')[0].radd(' ; ').fillna('')
        data['Inheritance'] = data['ZYG']
        data['Phenotype'] = data['PhenotypeList'].str.split('|').apply(lambda x: next((p for p in x if p not in ['not provided', 'not specified']), None))
        data['Classification'] = '(' + data['CLIN_SIG'] + ' ; ' + data['ClinicalSignificance (ClinVar)'] + ')'
        data['Allele State'] = data['ZYG']
        data['Allelic Read Depths'] = 'Ref(' + data['REF_ALLELE'] +'), Alt('+ data['Allele'] + ')' + ' VAF:' + (data['VAF (Sample)']*100).astype(str) + '%'
        data['Genomic Position'] = data['HGVSg']     
        data['Variant Frequency'] = data['MAX_AF (Population)'].apply(lambda x: 'Not identified in a large population studies' if pd.isna(x) or x== '' else (str(x*100) + '% Max frequency observed in Annotated 1000 genomes/ESP/gnomAD.'))
 
        filtered_data = data[( ( ( (data['CLIN_SIG'] == 'pathogenic') | (data['CLIN_SIG'] == 'pathogenic/likely_pathogenic') ) | ( (data['ClinicalSignificance (ClinVar)'] == 'Pathogenic') | (data['ClinicalSignificance (ClinVar)'] == 'Pathogenic/Likely pathogenic') ) ) & 
              ((data['MAX_AF (Population)'].isna()) | (data['MAX_AF (Population)'] == 0.001) | 
               (data['MAX_AF (Population)'] < 0.015)))]

        report_out = filtered_data[[
            'Avg Read Depth', 
            'Gene', 
            'Transcript',                
            'Variant',                   
            'Inheritance',               
            'Phenotype',                               
            'Classification',            
            'Location', 
            'Allele State',              
            'Allelic Read Depths',       
            'Genomic Position',
            'Variant Frequency'
        ]].copy()
        
        accession=base_name.split("_raw")[0]

        value_MD = data_1.loc[data_1['Sample'].str.contains(accession, na=False), 'Mean depth of coverage'].values[0]
        report_out['Avg Read Depth'] = value_MD       
        report_out['Avg Read Depth']=report_out['Avg Read Depth'].astype(str)+'X'

        value_P30 = data_1.loc[data_1['Sample'].str.contains(accession,na=False),'Percentage of bases covered at 30X'].values[0]
        report_out['Panel Coverage']=value_P30
        report_out['Panel Coverage']=report_out['Panel Coverage'].astype(str)+'%'
        report_out=report_out[[report_out.columns[-1]]+list(report_out.columns[:-1])]
        report_out = report_out.drop_duplicates(subset=['Gene', 'Variant', 'Variant Frequency','Phenotype'])

        filtered_data.to_csv(output_file, sep=",", index=False)
        report_out = report_out.drop_duplicates(subset=['Gene','Transcript','Variant'])
        report_out.to_csv(report_output_file, sep=",", index=False)
        print(f"Filtered file saved: {output_file}")