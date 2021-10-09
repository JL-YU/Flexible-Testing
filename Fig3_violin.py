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


p_start = 0.001
fig = plt.figure(1, figsize=(15,6))
ax1 = fig.add_subplot(2,2,1)
days = [20,50,80,110]
vl_data = pd.read_csv('.//p_start'+str(p_start)+'/'+'vl_nstar7_round_2delay.csv', usecols=[str(day) for day in days])
se_data = pd.read_csv('.//p_start'+str(p_start)+'/'+'se_nstar7_round_2delay.csv')
nstar_data = pd.read_csv('.//p_start'+str(p_start)+'/'+'exp_nstar7_round_2delay.csv')
# data preprocess: calculate nstar
se_data['n_star'] = nstar_data['n_star']
se_data['n_star'].fillna(1, inplace=True)
se_data['sed_nstar'] = 0

for i in se_data.index:
    se_data.loc[i, 'sed_nstar'] = se_data.loc[i, 'sed' + str(int(se_data.loc[i, 'n_star']))]
ax1,_ = plot_violin_and_other(ax1, days, vl_data, se_data, nstar_data, y2_col = 'sed_nstar', y2_name='Test sensitivity')


ax3 = fig.add_subplot(2,2,2)
days = [50,100,150,200] # days must be evenly distributed
vl_data = pd.read_csv('.//p_start'+str(p_start)+'/'+'vl_individual_round_1delay.csv', usecols=[str(day) for day in days])
se_data = pd.read_csv('.//p_start'+str(p_start)+'/'+'se_individual_round_1delay.csv')
# data preprocess: calculate nstar
se_data['n_star'] = 1
se_data['sed_nstar'] = se_data['sed1']
ax3,_ = plot_violin_and_other(ax3, days, vl_data, se_data, nstar_data, y2_col = 'sed_nstar', y2_name='Test sensitivity',LOD=3, title='(B) Individual PCR screening')


ax2 = fig.add_subplot(2,2,3)
days = [20,30,40,50]
vl_data = pd.read_csv('.//p_start'+str(p_start)+'/'+'vl_antigen_individual_round_3day.csv', usecols=[str(day) for day in days])
se_data = pd.read_csv('.//p_start'+str(p_start)+'/'+'se_antigen_individual_round_3day.csv')
# data preprocess: calculate nstar
se_data['n_star'] = 1
se_data['sed_nstar'] = se_data['sed1']
ax2,_ = plot_violin_and_other(ax2, days, vl_data, se_data, nstar_data, y2_col = 'sed_nstar', y2_name='Test sensitivity',LOD=5, title='(C) Antigen screeening (every 3 days)')


ax4 = fig.add_subplot(2,2,4)
days = [75,150,225,300] # days must be evenly distributed
vl_data = pd.read_csv('.//p_start'+str(p_start)+'/'+'vl_antigen_individual_round_14day.csv', usecols=[str(day) for day in days])
se_data = pd.read_csv('.//p_start'+str(p_start)+'/'+'se_antigen_individual_round_14day.csv')
# data preprocess: calculate nstar
se_data['n_star'] = 1
se_data['sed_nstar'] = se_data['sed1']
ax4,_ = plot_violin_and_other(ax4, days, vl_data, se_data, nstar_data, y2_col = 'sed_nstar', y2_name='Test sensitivity', LOD=5, title='(D) Antigen screeening (every 14 days)')



fig.text(0.5, 0.05,'Time (days)', ha='center',fontsize=15)
plt.subplots_adjust(wspace=0.27, hspace=0.4)
plt.savefig('.//p_start'+str(p_start)+'/'+'violin_'+str(p_start)+'.pdf', bbox_inches='tight')
