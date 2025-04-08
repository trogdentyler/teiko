import pandas as pd
import matplotlib.pyplot as plt
import os

from scipy import stats


cell_count = pd.read_csv('cell-count.csv')

# P.1)
# Please return an output file in CSV format with cell count and
# relative frequency of each population of each sample per line.

cell_names = ['b_cell', 'cd8_t_cell', 'cd4_t_cell', 'nk_cell', 'monocyte']
csv_files = [f'{name}.csv' for name in cell_names]
output_columns = ['sample', 'total_count', 'population', 'count', 'percentage']

for file, name in zip(csv_files, cell_names):
    
    file_path = os.path.join('./csvs/', file)
    if not os.path.exists(file_path):
        os.makedirs(os.path.dirname(file_path), exist_ok=True)

    df = pd.DataFrame(columns=output_columns)

    df['sample'] = cell_count['sample']
    df['total_count'] = cell_count[cell_names].sum(axis=1)
    df['population'] = name
    df['count'] = cell_count[name]
    df['percentage'] = cell_count[name].div(df['total_count'], axis=0)

    df.to_csv(file_path, index=False)

# P.2) (a)
# For each immune cell population, please generate a boxplot of the
# population relative frequencies comparing responders versus non-responders.

treated_mel = cell_count[(cell_count['treatment'] == 'tr1') & \
                         (cell_count['condition'] == 'melanoma') & \
                         (cell_count['sample_type'] == 'PBMC')]

resp_y = treated_mel[treated_mel['response'] == 'y']
resp_n = treated_mel[treated_mel['response'] == 'n']

test_dict = {}
os.makedirs('./plots', exist_ok=True)
for name in cell_names:
    rel_f_y = resp_y[name].div(resp_y[cell_names].sum(axis=1), axis=0)
    rel_f_n = resp_n[name].div(resp_n[cell_names].sum(axis=1), axis=0)

    test_dict[name] = (rel_f_y, rel_f_n)

    plt.boxplot([rel_f_y, rel_f_n], tick_labels=['Responders', 'Non-responders'])
    plt.title(f'Cell Type: {name}')
    plt.ylabel('Relative Frequency')

    plt.savefig(f'./plots/{name}_boxplot.png')
    plt.close()

# P.2) (b)
# Which cell populations are significantly different in relative frequencies
# between responders and non-responders? Include statistics to support your conclusion.
alpha = 0.05
for name in cell_names:
    rel_f_y, rel_f_n = test_dict[name]
    t, p = stats.ttest_ind(rel_f_y, rel_f_n)
    
    message = f'Significant difference found for {name} between responders and non-responders.\n' if p < alpha else f'No significant difference for {name} between responders and non-responders.\n'
    
    print(f'Cell Population: {name}')
    print(f'T-test statistic: {t:.4f}, p-value: {p:.4f}')
    print(message, '\n')
