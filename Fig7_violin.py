import pandas as pd
import numpy as np
import seaborn as sns
import statsmodels.api as sm
import matplotlib.pyplot as plt
lowess = sm.nonparametric.lowess

from plot_violin import plot_violin_and_other,plot_se_nstar


plt.rcParams.update({'font.size':13})
plt.rcParams.update({'axes.labelpad':9})
plt.rc('font', family = 'Times New Roman')
plt.rc('legend', fontsize=11)
N=100000



days = [50,100,150,200] # days must be evenly distributed
vl_data = pd.read_csv('.//p_start0.005/vl_nstar7_round_2delay.csv', usecols=[str(day) for day in days])
se_data = pd.read_csv('.//p_start0.005/se_nstar7_round_2delay.csv')
nstar_data = pd.read_csv('.//p_start0.005/exp_nstar7_round_2delay.csv')
# data preprocess: calculate nstar
se_data['n_star'] = nstar_data['n_star']
se_data['n_star'].fillna(1, inplace=True)
se_data['sed_nstar'] = 0




for i in se_data.index:
    se_data.loc[i, 'sed_nstar'] = se_data.loc[i, 'sed' + str(int(se_data.loc[i, 'n_star']))]

fig = plt.figure(1, figsize=(15,6))
ax1 = fig.add_subplot(2,2,1)
epi_data = pd.read_csv('.//p_start0.005/exp_nstar7_round_2delay.csv')
epi_data['TQR'] = epi_data.TP.cumsum()/(N-epi_data.S)

# epi_data = epi_data.loc(epi_data.TP>0)
ax1,_ = plot_violin_and_other(ax1, days, vl_data, epi_data, nstar_data, y2_col = 'TQR', y2_name='Testing–isolation ratio', LOD=3, title='(C) Pooled PCR screening ')



# se_data['day'] = se_data.index

# ax1 = plot_violin(ax1, days, vl_data, se_data, nstar_data, y2_col = 'sed_nstar', y2_name='Testing–isolation ratio')
# import pdb; pdb.set_trace()



ax2 = fig.add_subplot(2,2,2)
ax2,_ = plot_se_nstar(ax2, se_data, nstar_data ,cols=['sed_nstar'],max_day=max(days)+50, min_day=28, with_label=False, title='(D) Optimal pool sizes and sensitivity of pooled PCR screening')


vl_data = pd.read_csv('.//p_start0.005/vl_antigen_individual_round_14day.csv', usecols=[str(day) for day in days])
ax3 = fig.add_subplot(2,2,3)

epi_data = pd.read_csv('.//p_start0.005/exp_antigen_individual_round_14day.csv')
epi_data['TQR'] = epi_data.TP.cumsum()/(N-epi_data.S)
ax3,_ = plot_violin_and_other(ax3, days, vl_data, epi_data, nstar_data, y2_col = 'TQR', y2_name='Testing–isolation ratio', LOD=5, title='(E) Antigen screening (every 14 days) ')



days = [50,70, 90, 110,] # days must be evenly distributed


vl_data = pd.read_csv('.//p_start0.005/vl_antigen_individual_round_3day.csv', usecols=[str(day) for day in days])
epi_data = pd.read_csv('.//p_start0.005/exp_antigen_individual_round_3day.csv')
# data preprocess: calculate nstar
epi_data['TQR'] = epi_data.TP.cumsum()/(N-epi_data.S)
ax4 = fig.add_subplot(2,2,4)
# days = [50,100,150,200] # days must be evenly distributed
ax4,_ = plot_violin_and_other(ax4, days, vl_data, epi_data, nstar_data, y2_col = 'TQR', y2_name='Testing–isolation ratio', LOD=5, title='(F) Antigen screening (every 3 days) ')




fig.text(0.5, 0.05,'Time (days)', ha='center',fontsize=15)
plt.subplots_adjust(wspace=0.27, hspace=0.4)
# plt.show()
plt.savefig('.//p_start0.005/Fig7_violin.pdf', bbox_inches='tight')


