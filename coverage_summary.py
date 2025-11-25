import os
import pandas as pd
import glob

def process_coverage_summary(absolute_path, report_folder):

    folder_name = os.path.basename(os.path.dirname(absolute_path))

    output_file_name = f'{folder_name}_coverage_summary.csv'

    coverage_summary_path = os.path.join(report_folder, output_file_name)
    coverage_summary = pd.read_csv(coverage_summary_path, delimiter=',') 


    for index, row in coverage_summary.iterrows():
        sample_name = row['Sample']
        bam_name = str(row['Sample']).strip()
        
        if bam_name.endswith("_recal.bam"):
            sample_id = bam_name[:-len("_recal.bam")]
        else:
          	sample_id = bam_name
        print (sample_id)
        
        raw_report_path = os.path.join(report_folder, f'{sample_id}_raw_output_filtered_tp_report.csv')
        
        if os.path.exists(raw_report_path):
            raw_report = pd.read_csv(raw_report_path)
            
            if raw_report['Gene'].isna().all() and raw_report['Transcript'].isna().all():
                coverage_summary.at[index, 'Status'] = 'Negative'
            else:
                coverage_summary.at[index, 'Status'] = 'Positive'
                

            if row['Percentage of bases covered at 30X'] < 30:
                coverage_summary.at[index, 'Status'] = 'Qc Fail'

        else:
            print(f"Warning: {raw_report_path} not found.")
            coverage_summary.at[index, 'Status'] = 'Qc Fail'
            
            if sample_name.startswith('NTC'):
                if row['Percentage of bases covered at 30X'] < 1:
                    coverage_summary.at[index, 'Status'] = 'Pass'
                else:
                    coverage_summary.at[index, 'Status'] = 'Fail'
            else:
                coverage_summary.at[index, 'Status'] = 'Qc Fail'
                                  
            
        if sample_name.startswith('PC'):
            if row['Percentage of bases covered at 30X'] > 30:
                coverage_summary.at[index, 'Status'] = 'Pass'
            else:
                coverage_summary.at[index, 'Status'] = 'Fail'
    
    coverage_summary.to_csv(coverage_summary_path, index=False)
    print(f"Updated coverage summary saved to: {coverage_summary_path}")

file_paths = glob.glob('./*_coverage.txt')
absolute_path = os.path.abspath(file_paths[0])
report_folder = './'

process_coverage_summary(absolute_path, report_folder)
