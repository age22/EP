import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import logging

df = pd.read_csv('data/Output/plink/2-POP/PCA_results.eigenvec', delim_whitespace=True, header=None)

cols = ['FID', 'IID']
for i in range(1, 21):
    cols.append(f'PC{i}')

df.columns = cols

colorcode = pd.read_csv('data/Output/plink/2-POP/popfile.txt', delim_whitespace=True)


df = pd.merge(df, colorcode, on=['FID', 'IID'])


ax = sns.lmplot('PC1', # Horizontal axis
           'PC2', # Vertical axis
           hue = 'SUPERPOP',  # color variable
           data=df, # Data source
           fit_reg=False, # Don't fix a regression line
           height = 10,
           aspect =2 ) # height and dimension

plt.title('PCA: Projection on 1000Genomes')
# Set x-axis label
plt.xlabel(f'PC1')
# Set y-axis label
plt.ylabel('PC2')

plt.savefig('data/PCA.png', dpi=240, bbox_inches='tight')


def label_point(x, y, val, ax):
    a = pd.concat({'x': x, 'y': y, 'val': val}, axis=1)
    for i, point in a.iterrows():
        ax.text(point['x']+.0001, point['y'], str(point['val']))

label_point(df[df['SUPERPOP'] == 'EP']['PC1'], df[df['SUPERPOP'] == 'EP']['PC2'], df[df['SUPERPOP'] == 'EP']['IID'], plt.gca())

plt.savefig('data/PCA_labels.png', dpi=240, bbox_inches='tight')
