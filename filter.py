import io
import os
import pandas as pd
import glob
import numpy as np

import re
import pandas as pd
import numpy as np

def read_vcf_no_meta(path: str) -> pd.DataFrame:
    header = None
    rows = []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            if raw.startswith("##"):
                continue
            if raw.startswith("#CHROM"):
                header = raw.lstrip("#").rstrip("\n").split("\t")
                continue
            if raw.strip() == "":
                continue
            rows.append(raw.rstrip("\n").split("\t"))
    if header is None:
        raise ValueError(f"No #CHROM header found in {path}")
    df = pd.DataFrame(rows, columns=header, dtype=str)
    df = df.rename(columns={"CHROM": "CHROM", "#CHROM": "CHROM"})
    return df

def parse_sample_by_format(df: pd.DataFrame) -> pd.DataFrame:
    # take the last column as the sample column
    sample_col = df.columns[-1]

    fmt_keys = df["FORMAT"].fillna("").str.split(":")
    sample_vals = df[sample_col].fillna("").str.split(":")

    parsed = []
    for keys, vals in zip(fmt_keys, sample_vals):
        d = {}
        for k, v in zip(keys, vals):
            d[k] = v
        parsed.append(d)

    parsed_df = pd.DataFrame(parsed)

    for col in ["GT", "DP", "GQ", "AD", "AO", "RO"]:
        if col not in parsed_df.columns:
            parsed_df[col] = np.nan

    out = pd.concat([df.reset_index(drop=True), parsed_df], axis=1)

    ad_alt = out["AD"].astype(str).str.split(",").str[1]
    out["AD_ALT"] = pd.to_numeric(ad_alt, errors="coerce")
    out["DP_num"] = pd.to_numeric(out["DP"], errors="coerce")
    out["VAF_fb"] = out["AD_ALT"] / out["DP_num"]

    return out


def load_freebayes_metrics(freebayes_path: str) -> pd.DataFrame:
    fb = read_vcf_no_meta(freebayes_path)
    fb = parse_sample_by_format(fb)  # no sample_col arg now

    # normalize types for merge
    fb["POS"] = pd.to_numeric(fb["POS"], errors="coerce")

    metrics = fb[["CHROM", "POS", "REF", "ALT", "GT", "DP", "GQ", "VAF_fb"]].copy()
    metrics = metrics.rename(
        columns={
            "GT": "GT_fb",
            "DP": "DP_fb",
            "GQ": "GQ_fb"
        }
    )
    return metrics


###
import numpy as np
import pandas as pd

def add_lookup_pos(df: pd.DataFrame,
                   location_col="Location",
                   allele_col="Allele",
                   out_col="POS_LOOKUP") -> pd.DataFrame:
    start_pos = (
        df[location_col]
        .astype(str)
        .str.split(r"[:|-]", regex=True)
        .str[1]
    )
    start_pos = pd.to_numeric(start_pos, errors="coerce")

    is_del = df[allele_col].astype(str).str.strip().eq("-")

    df[out_col] = np.where(is_del, start_pos - 1, start_pos)
    return df

  
def add_pos_candidates(df, location_col="Location"):
    start_pos = (
        df[location_col]
        .astype(str)
        .str.split(r"[:|-]", regex=True)
        .str[1]
    )
    start_pos = pd.to_numeric(start_pos, errors="coerce")

    df["POS_START"] = start_pos
    df["POS_ANCHOR"] = start_pos - 1
    return df

  
def load_freebayes_full(freebayes_path: str) -> pd.DataFrame:
    fb = read_vcf_no_meta(freebayes_path)
    #fb = parse_sample_by_format(fb, sample_col="sample")
    fb = parse_sample_by_format(fb)

    fb["POS"] = pd.to_numeric(fb["POS"], errors="coerce")
    fb["ALT_LIST"] = fb["ALT"].astype(str).str.split(",")
    fb["AD_LIST"] = fb["AD"].astype(str).str.split(",")  

    fb["AD_LIST"] = fb["AD_LIST"].apply(
        lambda xs: [pd.to_numeric(x, errors="coerce") for x in xs] if isinstance(xs, list) else []
    )

    fb = fb.rename(columns={"GT": "GT_fb", "DP": "DP_fb", "GQ": "GQ_fb"})
    return fb[["CHROM", "POS", "REF", "ALT", "ALT_LIST", "AD_LIST", "GT_fb", "DP_fb", "GQ_fb"]]

import io, re
import pandas as pd

def read_vep_tab(path: str) -> pd.DataFrame:
    header_cols = None
    data_rows = []

    with open(path, 'r', encoding='utf-8', errors='replace') as fh:
        for raw in fh:

            if raw.lstrip().startswith('##'):
                continue

            line = raw.lstrip().rstrip('\n')

            if header_cols is None and line.startswith('#'):

                header_cols = re.split(r'\t+', line[1:].strip())
                continue

            if line == '':
                continue  


            fields = line.split('\t')
            data_rows.append(fields)

    if header_cols is None:
        raise ValueError(f"No VEP header ('#Uploaded_variation ...') found in {path}")


    n = len(header_cols)
    fixed_rows = []
    for r in data_rows:
        if len(r) < n:
            r = r + [''] * (n - len(r))    
        elif len(r) > n:
            r = r[:n]               
        fixed_rows.append(r)

    df = pd.DataFrame(fixed_rows, columns=header_cols, dtype=str)
    return df
  
def process_vcf(
        file_pattern: str, 
        variant_summary_path: str
    ):
    files = glob.glob(file_pattern)
    y=pd.read_csv(variant_summary_path, sep='\t')

    for file in files:
        import numpy as np
      
        vcf = read_vep_tab(file)
        if '#CHROM' in vcf.columns:
            vcf = vcf.rename(columns={'#CHROM': 'CHROM'})
        print(file)

        vcf['Existing_variation'] = vcf['Existing_variation'].str.split(',')
        vcf = vcf.explode('Existing_variation')
        
     
        vcf = vcf.dropna(subset=['Existing_variation'])
        vcf['Existing_variation'] = vcf['Existing_variation'].astype(str)
        


        x = vcf[vcf['Existing_variation'].str.startswith('rs')]
        x["rs_id"]=x['Existing_variation'].str.replace('rs', '')

        x['rs_id'] = pd.to_numeric(x['rs_id'], errors='coerce', downcast='integer')


        result_df = pd.merge(x, y, left_on='rs_id', right_on='RS# (dbSNP)', how='inner')
        

        a = result_df[result_df['Assembly'].str.startswith('GRCh37')]
        b = a.drop_duplicates().reset_index(drop=True)

        b['CLIN_SIG'] = b['CLIN_SIG'].str.split(',')
        c = b.explode('CLIN_SIG')

        e = c.drop_duplicates().reset_index(drop=True)

        final = e.drop_duplicates(subset=['SPDI', 'CLIN_SIG', 'REF_ALLELE', 'Allele', 
                                          'Existing_variation', 'CLIN_SIG', 'PhenotypeList'], 
                                          keep='first').reset_index(drop=True)
                                          

        import numpy as np       
        import os

        fname = os.path.basename(file)

        if fname.endswith('_merged_vep.vcf'):
            fil = fname[:-len('_merged_vep.vcf')]       
        elif fname.endswith('_vep.vcf'):
            fil = fname[:-len('_vep.vcf')]             
        else:
            fil = fname.split('_vep')[0]              

        partner_name = fil + '_variants.vcf'
        partner_path = os.path.join(os.path.dirname(file), partner_name)

        with open(partner_path, 'r') as f1:
            lines_1 = [l for l in f1 if not l.startswith('##')]
        


        haplocall=pd.read_csv(io.StringIO(''.join(lines_1)),dtype=str,sep='\t').rename(columns={'#CHROM': 'CHROM'})
        
        haplocall = pd.read_csv(
            io.StringIO(''.join(lines_1)),
            dtype=str,
            sep='\t'
        ).rename(columns={'#CHROM': 'CHROM'})

        # take the last column as the sample genotype column
        sample_col = haplocall.columns[-1]

        
        final = add_pos_candidates(final, location_col="Location")
        haplocall["POS"] = pd.to_numeric(haplocall["POS"], errors="coerce")

        m1 = final.merge(haplocall, left_on="POS_START", right_on="POS", how="left")
        need2 = m1[sample_col].isna()

        m2 = final.loc[need2].merge(haplocall, left_on="POS_ANCHOR", right_on="POS", how="left")
        m1.loc[need2, haplocall.columns] = m2[haplocall.columns].values

        DP_GQ_Zygo = m1

        DP_GQ_Zygo['Genotype (GT)'] = DP_GQ_Zygo[sample_col].str.split(':').str[0]
        DP_GQ_Zygo['Depth of Coverage (DP)'] = DP_GQ_Zygo[sample_col].str.split(':').str[2]
        DP_GQ_Zygo['Genotype Quality (GQ)'] = DP_GQ_Zygo[sample_col].str.split(':').str[3]
        DP_GQ_Zygo['AD_ALT'] = DP_GQ_Zygo[sample_col].str.split(':').str[1].str.split(',').str[1]



        DP_GQ_Zygo['AD_ALT'] = pd.to_numeric(DP_GQ_Zygo['AD_ALT'], errors='coerce')
        DP_GQ_Zygo['Depth of Coverage (DP)'] = pd.to_numeric(DP_GQ_Zygo['Depth of Coverage (DP)'], errors='coerce')

        DP_GQ_Zygo['VAF (Sample)'] = (DP_GQ_Zygo['AD_ALT'] / DP_GQ_Zygo['Depth of Coverage (DP)'])#*100

        freebayes_name = fil + "_freebayes.vcf"
        freebayes_path = os.path.join(os.path.dirname(file), freebayes_name)

        if os.path.exists(freebayes_path):
            fb_full = load_freebayes_full(freebayes_path)

            DP_GQ_Zygo["CHROM"] = DP_GQ_Zygo["Location"].astype(str).str.split(r"[:|-]").str[0]
            DP_GQ_Zygo = add_pos_candidates(DP_GQ_Zygo, location_col="Location")


            m1 = DP_GQ_Zygo.merge(
                fb_full,
                left_on=["CHROM", "POS_START"],
                right_on=["CHROM", "POS"],
                how="left",
                suffixes=("", "_fb")
            )

            need2 = m1["GT_fb"].isna()
            m2 = m1.loc[need2].drop(columns=fb_full.columns.difference(["CHROM"]))  
            m2 = m2.merge(
                fb_full,
                left_on=["CHROM", "POS_ANCHOR"],
                right_on=["CHROM", "POS"],
                how="left",
                suffixes=("", "_fb")
            )

            m1.loc[need2, fb_full.columns] = m2[fb_full.columns].values
            DP_GQ_Zygo = m1

            def pick_alt_idx(allele, alt_list):
                if not isinstance(alt_list, list) or len(alt_list) == 0:
                    return np.nan
                allele = str(allele)


                for i, a in enumerate(alt_list):
                    if allele == a:
                        return i + 1


                for i, a in enumerate(alt_list):
                    if allele in a or a in allele:
                        return i + 1


                lens = [abs(len(a) - len(allele)) for a in alt_list]
                return int(np.argmin(lens)) + 1

            DP_GQ_Zygo["ALT_IDX"] = DP_GQ_Zygo.apply(
                lambda r: pick_alt_idx(r["Allele"], r["ALT_LIST"]),
                axis=1
            )

            def get_ad_alt(ad_list, idx):
                if not isinstance(ad_list, list) or pd.isna(idx):
                    return np.nan
                idx = int(idx)
                return ad_list[idx] if idx < len(ad_list) else np.nan

            DP_GQ_Zygo["AD_ALT_fb"] = DP_GQ_Zygo.apply(
                lambda r: get_ad_alt(r["AD_LIST"], r["ALT_IDX"]),
                axis=1
            )

            DP_GQ_Zygo["VAF_fb"] = (
                pd.to_numeric(DP_GQ_Zygo["AD_ALT_fb"], errors="coerce") /
                pd.to_numeric(DP_GQ_Zygo["DP_fb"], errors="coerce")
            )


            DP_GQ_Zygo["Genotype (GT)"] = DP_GQ_Zygo["Genotype (GT)"].where(
                DP_GQ_Zygo["Genotype (GT)"].notna() & (DP_GQ_Zygo["Genotype (GT)"] != ""),
                DP_GQ_Zygo["GT_fb"]
            )

            DP_GQ_Zygo["Depth of Coverage (DP)"] = DP_GQ_Zygo["Depth of Coverage (DP)"].where(
                DP_GQ_Zygo["Depth of Coverage (DP)"].notna(),
                pd.to_numeric(DP_GQ_Zygo["DP_fb"], errors="coerce")
            )

            DP_GQ_Zygo["Genotype Quality (GQ)"] = DP_GQ_Zygo["Genotype Quality (GQ)"].where(
                DP_GQ_Zygo["Genotype Quality (GQ)"].notna(),
                pd.to_numeric(DP_GQ_Zygo["GQ_fb"], errors="coerce")
            )

            DP_GQ_Zygo["VAF (Sample)"] = DP_GQ_Zygo["VAF (Sample)"].where(
                DP_GQ_Zygo["VAF (Sample)"].notna(),
                DP_GQ_Zygo["VAF_fb"]
            )

        else:
            print(f"[INFO] FreeBayes file not found: {freebayes_path}")

        out = DP_GQ_Zygo[['Location', 'REF_ALLELE', 'Allele', 'BIOTYPE', 'Existing_variation','SYMBOL', 'GeneSymbol','CLIN_SIG','ClinicalSignificance','SPDI','HGVSg','HGVSc','HGVSp','HGVS_OFFSET','SIFT', 'PolyPhen', 'VAF (Sample)','MAX_AF', 'Genotype (GT)','ZYG','Depth of Coverage (DP)','Genotype Quality (GQ)','PhenotypeIDS','PhenotypeList']]

        out = out.drop_duplicates().reset_index(drop=True)
        out.rename(columns={'ClinicalSignificance': 'ClinicalSignificance (ClinVar)', 'MAX_AF' : 'MAX_AF (Population)'}, inplace=True)


        out_1 = out.drop_duplicates(subset=out.columns.difference(['PhenotypeIDS','PhenotypeList','ClinicalSignificance (ClinVar)']))
        out_1 = out
        out_1.to_csv(fil+'_raw_output.csv')



        out_f1=out[out["CLIN_SIG"].isin(["pathogenic", "pathogenic/established_risk_allele","likely_pathogenic","pathogenic/likely_pathogenic"])]

        out_f1 = out_f1.drop_duplicates(subset=out_f1.columns.difference(['CLIN_SIG','ClinicalSignificance (ClinVar)','PhenotypeIDS','PhenotypeList']))

        out_f1.to_csv(fil+'_output_p.csv', index=False)
