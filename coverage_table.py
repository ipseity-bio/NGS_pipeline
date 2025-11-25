import pandas as pd
import glob
import os

data = []

file_paths = glob.glob('./*_coverage.txt')

for file_path in file_paths:
    with open(file_path, 'r') as file:
        content = file.readlines()

        file_name_1 = file_path.split('/')[-1] 
        file_name=str(file_name_1.split("_all_coverage")[0])
        coverage_1X = float(content[0].split(': ')[1].strip().replace('%', ''))
        coverage_10X = float(content[1].split(': ')[1].strip().replace('%', ''))
        coverage_20X = float(content[2].split(': ')[1].strip().replace('%', ''))
        coverage_30X = float(content[3].split(': ')[1].strip().replace('%', ''))
        coverage_40X = float(content[4].split(': ')[1].strip().replace('%', ''))
        coverage_50X = float(content[5].split(': ')[1].strip().replace('%', ''))
        coverage_100X = float(content[6].split(': ')[1].strip().replace('%', ''))
        coverage_500X = float(content[7].split(': ')[1].strip().replace('%', ''))
        mean_depth = float(content[8].split(': ')[1].strip())


        data.append([file_name, coverage_1X, coverage_10X, coverage_20X, coverage_30X, coverage_40X, coverage_50X, coverage_100X, coverage_500X, mean_depth])


columns = ['Sample', 'Percentage of bases covered at 1X', 'Percentage of bases covered at 10X', 'Percentage of bases covered at 20X',
           'Percentage of bases covered at 30X', 'Percentage of bases covered at 40X', 'Percentage of bases covered at 50X',
           'Percentage of bases covered at 100X', 'Percentage of bases covered at 500X', 'Mean depth of coverage']
df = pd.DataFrame(data, columns=columns)

df['Sort Key'] = df['Sample'].str.extract('S(\d+)', expand=False).astype(int)
df = df.sort_values(by='Sort Key')
df = df.drop('Sort Key', axis=1)


absolute_path = os.path.abspath(file_paths[0])


folder_name = os.path.basename(os.path.dirname(absolute_path))

output_file_name = f'{folder_name}_coverage_summary.csv'
df = df.fillna(0)
df.to_csv(output_file_name, index=False)

