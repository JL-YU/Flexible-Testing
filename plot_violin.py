import pandas as pd
import numpy as np
import seaborn as sns
import statsmodels.api as sm
import matplotlib.pyplot as plt
lowess = sm.nonparametric.lowess

plt.rcParams.update({'font.size':13})
plt.rcParams.update({'axes.labelpad':9})
plt.rc('font', family = 'Times New Roman')
plt.rc('legend', fontsize=11)


def plot_violin_and_other(ax, days, vl_data, se_data, nstar_data, y2_col = 'sed_nstar', y2_name='Testing–isolation ratio', title='(A) Population-level viral load distribution',LOD=3):
    # y2_col can be {'sed_nstar', 'sed1','sed10',....,'sep1','sep10',...}
    vl_data.rename({str(day): day for day in days}, axis='columns', inplace=True)
    vl_data.replace(0, np.nan, inplace=True)

    ax = sns.violinplot(data=vl_data, color="skyblue", inner='quartile',ax=ax, cut=0,)
    ax = sns.swarmplot(data=vl_data.sample(n=3000), color='black', size= 3,edgecolor='gray', ax=ax)

    ax.set_ylabel('Log10 viral load')
    ax.tick_params(axis='y')
    ax.set_title(title)
    ax.set_ylim(0,15)

    if LOD == 3:
        color = 'purple'
        ls = 'dotted'
    else:
        color = 'magenta'
        ls = 'dashdot'
    ax.axhline( LOD, color = color, ls = ls)

    ax2 = ax.twinx()
    data_y = se_data.loc[0:days[-1], y2_col]
    data_x = 1 / days[-1] * (len(days) - 1) * np.arange(days[-1] + 1)

    if y2_name == 'Testing–isolation ratio':
        cl = 'orange'
    else:
        cl = 'red'

    ax2.scatter(data_x, data_y, s = 3, color = cl)
    ax2.set_ylabel(y2_name)
    ax2.tick_params(axis='y', labelcolor = cl)
    ax2.set_ylim(0,1)
    return ax, ax2





def plot_violin(ax, days, vl_data, se_data, nstar_data, y2_col = 'sed_nstar', y2_name='Testing–isolation ratio', title='(A) Pooled PCR screening',LOD=3):
    # y2_col can be {'sed_nstar', 'sed1','sed10',....,'sep1','sep10',...}
    vl_data.rename({str(day): day for day in days}, axis='columns', inplace=True)
    vl_data.replace(0, np.nan, inplace=True)

    ax = sns.violinplot(data=vl_data,  color="skyblue", inner='quartile',ax=ax, cut=0,)
    ax = sns.swarmplot(data=vl_data.sample(n=3000), color='black', size= 3, edgecolor='gray', ax=ax)

    ax.set_ylabel('Log10 viral load')
    ax.set_ylim(0,15)
    ax.tick_params(axis='y')
    ax.set_title(title)
    ax.axhline( 3, color = 'green')

    return ax



def plot_se_nstar(ax, se_data, nstar_data ,cols=['sed_nstar'],max_day=250, min_day=30, with_label=False, title='(B) Optimal pool sizes and sensitivities of pooled PCR screening'):
    # cols will be data used for se plot, should be a list like ['sed1','sed2',...]
    # max_day and min_day gives the range
    days = np.arange(min_day,max_day+1)

    ax.plot(days, nstar_data.loc[min_day:max_day,'n_star'],color='blue')
    ax.set_ylim(1,32)
    ax.set_xlim(0,max_day)
    ax.set_ylabel('Optimal pool size')
    ax.tick_params(axis='y',labelcolor='blue')

    ax2 = ax.twinx()

    for col_name in cols:
        if with_label:
            if col_name=='sed_nstar':
                ax2.scatter(days,se_data.loc[min_day:max_day,col_name],label='n=$n^*$', s= 3, color = 'red')
            else:
                ax2.scatter(days,se_data.loc[min_day:max_day, col_name], label='n='+col_name[3:], s= 3, color = 'red')
        else:
            ax2.scatter(days, se_data.loc[min_day:max_day, col_name], s= 3, color = 'red')
    ax2.set_title(title)
    ax2.set_ylim(0,1)

    ax2.set_ylabel('Sensitivity of Dorfman pooling')
    ax2.tick_params(axis='y',labelcolor='red')

    return ax,ax2
