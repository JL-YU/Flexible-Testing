"""
do lowess and draw se/cpr and n_star graphs
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
lowess = sm.nonparametric.lowess
import seaborn as sns

from sir_simulation_for_se import get_n_star,lowess_data,lowess_data_time
# from plot_violin import plot_violin, plot_violin_and_other

plt.rcParams.update({'font.size':13})
plt.rcParams.update({'axes.labelpad':9})
plt.rc('font', family = 'Times New Roman')
plt.rc('legend', fontsize=11)

def plot_violin(ax, days, vl_data, se_data, nstar_data, y2_col = 'sed_nstar', y2_name='Test sensitivity', title='(A) Population-level viral load distribution of pooled PCR screening',LOD=3):
    # y2_col can be {'sed_nstar', 'sed1','sed10',....,'sep1','sep10',...}
    vl_data.rename({str(day): day for day in days}, axis='columns', inplace=True)
    vl_data.replace(0, np.nan, inplace=True)
    ax = sns.violinplot(data=vl_data,  color="skyblue", inner='quartile',ax=ax, cut=0,)
    ax = sns.swarmplot(data=vl_data.sample(n=3000), color='black', size= 3, edgecolor='gray', ax=ax)
    ax.set_ylabel('Log10 viral load')
    ax.set_ylim(0,15)
    ax.tick_params(axis='y')
    ax.set_title(title)
    ax.axhline( 3, color = 'purple',ls='--')
    ax.axhline( 5, color = 'magenta',ls= '-.')


def lowess_data_time(n_list, df_se, suffix='sep'):
    # fit curve se
    for n in n_list:
        x = df_se['time'].values
        y = df_se[suffix + str(n)].values
        z = lowess(y, x, return_sorted=False, frac=1. / 4)
        z[z > 1.] = 1.
        df_se[suffix + str(n) + '_lws'] = z
    return df_se


n_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30]
cmap= {1:'purple', 2:'deeppink',3:'crimson', 4:'orange', 5:'gold',
            6:'yellow',7:'greenyellow',8:'lawngreen', 9:'green', 10:'lightseagreen', 15:'cyan',20: 'teal',
            25:'skyblue',30: 'blue'}



days = [25,50,75,100,125,150,175,200] # days must be evenly distributed
vl_data = pd.read_csv('.//p_start0.001/vl_no_testing.csv', usecols=[str(day) for day in days])
se_data = pd.read_csv('.//p_start0.001/se_no_testing.csv')
# data preprocess: calculate nstar
nstar_data = pd.read_csv('.//p_start0.001/exp_nstar7_round_2delay.csv')
se_data['n_star'] = 1
se_data['sed_nstar'] = se_data['sed1']


# lowess regression of time & sensitivity
pcr_all = pd.read_csv('.//no_test/100_pcr_all.csv')
pcr_all = pcr_all.loc[pcr_all.p>=0.0001,]
pcr_all.sort_values('time', inplace=True)
se_all_time_sep = lowess_data_time(n_list, pcr_all, 'sep')
se_all_time_sep.to_csv('.//no_test/100_pcr_sep_time.csv')

se_all_time_sed = lowess_data_time(n_list, pcr_all, 'sed')
se_all_time_sed.to_csv('.//no_test/100_pcr_sed_time.csv')


# lowess regression of time & sensitivity
antigen_all = pd.read_csv('.//no_test/100_antigen_all.csv')
n_list_antigen= [1]
antigen_all.sort_values('time', inplace=True)
antigen_all_time_sep = lowess_data_time(n_list_antigen, antigen_all, 'sep')
antigen_all_time_sep.to_csv('.//no_test/100_pcr_sep_time.csv')



fig = plt.figure(1, figsize=(15, 6))
axvl = fig.add_subplot(2,1,1)
axvl.broken_barh([(25, 100), (100, 200), ], (0, 15),  facecolors=('yellow', 'green',),alpha = 0.05)
axvl = plot_violin(axvl, days, vl_data, se_data, nstar_data, y2_col = 'sed_nstar', y2_name='Test sensitivity', title= '(A) Viral load distribution without testing-isolation policy ')

ax1 = fig.add_subplot(2,2,3)
ax1.broken_barh([(0, 100), (100, 200), ], (0, 1),  facecolors=('yellow', 'green',),alpha =0.05)
for n in n_list:  ax1.plot(se_all_time_sep.time, se_all_time_sep['sep' + str(n) + '_lws'], label='PCR n=' + str(n), color=cmap[n])
ax1.set_ylim([0,1])
ax1.set_xlim([0,200])
ax1.set_title('(B) Sensitivity for pooled samples (first stage)')
ax1.set_ylabel('Sensitivity')

ax2 = fig.add_subplot(2,2,4)
ax2.broken_barh([(0, 100), (100, 200), ], (0, 1), facecolors=('yellow', 'green',),alpha = 0.05)
for n in n_list:  ax2.plot(se_all_time_sed.time, se_all_time_sed['sed' + str(n) + '_lws'], label='PCR n=' + str(n), color=cmap[n])
plt.plot(antigen_all_time_sep.time, antigen_all_time_sep.sep1_lws,label = 'Antigen',color = 'magenta')
ax2.legend()
ax2.set_ylim([0,1])
ax2.set_xlim([0,200])
ax2.set_ylabel('Sensitivity')
ax2.set_title('(C) Overall sensitivity for two-stage Dorfman pooling')

fig.text(0.5, 0.02,'Time (days)', ha='center',fontsize=15)
plt.subplots_adjust(wspace=0.27, hspace=0.4)
ax2.legend(ncol=1,bbox_to_anchor=(1,2),loc='upper left')
plt.savefig('.//no_test/Fig2.pdf', bbox_inches='tight')


