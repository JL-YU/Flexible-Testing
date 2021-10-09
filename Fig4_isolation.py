import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 13})
plt.rcParams.update({'axes.labelpad': 11})
plt.rc('font', family='Times New Roman')
plt.rc('legend', fontsize=11)

def plot_for_each_var(exp_names, data_list, variable_interested, x_var_interested, exp_interested,
                      subplot_titles,
                      subplot_x_names,
                      subplot_y_names,
                      markers,
                      xlims,
                      ylims):
    # each subplot is a var
    fig = plt.figure(figsize=(15, 9))
    subplot_number = len(variable_interested)
    width = 2
    height = int(np.ceil(subplot_number / width))
    sub_list = []
    for var_i, var in enumerate(variable_interested):
        sub_list.append(fig.add_subplot(height, width, var_i + 1))
        for exp_index, exp_i in enumerate(exp_interested):
            exp_name = exp_names[exp_i]
            if x_var_interested[exp_index] is None:
                sub_list[var_i].plot(data_list[exp_i][var],
                                     markers[exp_index], label=exp_name)
            else:
                sub_list[var_i].plot(data_list[exp_i][x_var_interested[exp_index]], data_list[exp_i][var], markers[exp_index],
                                     label=exp_name)
            sub_list[var_i].set_xlabel(subplot_x_names[var_i])
            sub_list[var_i].set_ylabel(subplot_y_names[var_i])
            if xlims[var_i] is None:
                sub_list[var_i].autoscale(axis='x')
            else:
                sub_list[var_i].set_xlim(xlims[var_i])
            if ylims[var_i] is None:
                sub_list[var_i].autoscale(axis='y')
            else:
                sub_list[var_i].set_ylim(ylims[var_i])
            sub_list[var_i].set_title(subplot_titles[var_i])
        if var_i == 0:
            sub_list[var_i].legend()
    return fig, sub_list


def plot_for_each_exp(exp_names, data_list, variable_interested, x_var_interested, exp_interested,
                      subplot_titles,
                      subplot_x_names,
                      subplot_y_names,
                      markers,
                      xlims,
                      ylims,label_list):
    # each subplot is a var
    fig = plt.figure(figsize=(15, 6))
    subplot_number = len(exp_interested)
    width = 2
    height = int(np.ceil(subplot_number / width))
    sub_list = []
    for exp_index, exp_i in enumerate(exp_interested):
        sub_list.append(fig.add_subplot(height, width, exp_index + 1))
        for var_i, var in enumerate(variable_interested):
            var_name = var
            if x_var_interested[var_i] is None:
                sub_list[exp_index].plot(data_list[exp_i][var],
                                         markers[exp_i], label=label_list[var_i])
            else:
                sub_list[exp_index].plot(data_list[exp_i][x_var_interested[var_i]], data_list[exp_i][var],
                                         markers[var_i], label=label_list[var_i])
            sub_list[exp_index].set_xlabel(subplot_x_names[exp_index])
            sub_list[exp_index].set_ylabel(subplot_y_names[exp_index])
            if xlims[exp_index] is None:
                sub_list[exp_index].autoscale(axis='x')
            else:
                sub_list[exp_index].set_xlim(xlims[exp_index])
            if ylims[exp_index] is None:
                sub_list[exp_index].autoscale(axis='y')
            else:
                sub_list[exp_index].set_ylim(ylims[exp_index])
            sub_list[exp_index].set_title(subplot_titles[exp_index])
        if exp_index == 1:
            sub_list[exp_index].legend()
        plt.subplots_adjust(wspace=0.27, hspace=0.4)
    return fig, sub_list




################## change file name here #################################################################
p_start = 0.001
exp_names = ['PCR_pooled test', 'PCR Individual test', 'Antigen test(every 3 day)','Antigen test(every 14 day)', 'No test']
file_names = ['exp_nstar7_round_2delay.csv',
              'exp_individual_round_1delay.csv',
              'exp_antigen_individual_round_3day.csv',
              'exp_antigen_individual_round_14day.csv',
              'exp_no_testing.csv',
              ]

data_list = [pd.read_csv('.//p_start'+str(p_start)+'/'+ name) for name in file_names]


for df in data_list:
    df['total_I'] = df['I'] + df['Q'] + df['SQ']
    df['cum_Infetion_frac'] = (df['I'] + df['Q'] + df['SQ'] + df['R']) / (df['I'] + df['Q'] + df['SQ'] + df['R'] + df['S'])
    df['total_I_frac'] = (df['I'] + df['Q'] + df['SQ'] + df['R']) / (df['I'] + df['Q'] + df['SQ'] + df['R'] + df['S'])
    df['I_frac'] = df['I']
    df['Q_frac'] = df['Q']
    df['SQ_frac'] = df['SQ']
    df['R_frac'] = df['R'] / (df['I'] + df['Q'] + df['SQ'] + df['R'] + df['S'])
    df['total_test_cumsum'] = df['number_of_total_tests'].cumsum()
    df['FN_CS'] = df['FN'].cumsum()
    df['FP_CS'] = df['FP'].cumsum()


fig_exp, sub_list_exp = plot_for_each_exp(exp_names='Infection',
                                          data_list=data_list,
                                          variable_interested=['I', 'Q', 'SQ'],
                                          # len(variable_interesed) =number of lines on each fig
                                          markers=['r:', 'm-.', 'k--'],  # len(markers) =number of lines on each fig,
                                          label_list=[ 'Infected indivs. (mixing in population)','Infected indivs. (isolated due to test)','Infected indivs. (isolated due to symptoms)'],
                                          x_var_interested=[ 'day', 'day','day','day'],
                                          # len(exp_interested) = number of figs
                                          exp_interested=[0, 1, 2, 3],
                                          subplot_titles=['(C) PCR pooled screening', '(D) PCR individual screening', '(E) Rapid antigen screening (every 3 days)', '(F) Rapid antigen screening (every 14 days)'],
                                          subplot_x_names=[None, None, None, None],
                                          subplot_y_names=['Count', 'Count', 'Count', 'Count'],
                                          xlims=[[-5,350], [-5,350], [-5,350], [-5,350]],
                                          # ylims=[[0.5,17000], [0.5,17000], [0.5,17000], [0.5,17000]]，
                                          ylims=[None, None, None, None])
fig_exp.text(0.5, 0.05,'Time (days)', ha='center',fontsize=15)
plt.savefig('.//p_start'+str(p_start)+'/'+'isolation_'+str(p_start)+'.pdf', bbox_inches='tight')

