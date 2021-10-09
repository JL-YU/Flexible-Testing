import pandas as pd
import matplotlib.pyplot as plt
import multiprocessing
cpu_count = multiprocessing.cpu_count()



plt.rcParams.update({'font.size':13})
plt.rcParams.update({'axes.labelpad':9})
plt.rc('font', family = 'Times New Roman')
plt.rc('legend', fontsize=11)

p_start = 0.001
N=100000


SIM_i_alpha_1 = pd.read_csv('.//p_start'+str(p_start)+'/'+'exp_individual_round_1delay.csv')
SIM_g_alpha_2 = pd.read_csv('.//p_start'+str(p_start)+'/'+'exp_nstar7_round_2delay.csv')
SIM_0_alpha = pd.read_csv('.//p_start'+str(p_start)+'/'+'exp_no_testing.csv')
SIM_antigen_alpha_3 = pd.read_csv('.//p_start'+str(p_start)+'/'+'exp_antigen_individual_round_3day.csv')
SIM_antigen_alpha_14 = pd.read_csv('.//p_start'+str(p_start)+'/'+'exp_antigen_individual_round_14day.csv')


fig = plt.figure(figsize=(15, 3))

ax_1 = fig.add_subplot(1, 2, 1)
ax_2 = fig.add_subplot(1, 2, 2)

ax_1.plot(SIM_g_alpha_2.day, 1 - SIM_g_alpha_2.S/N, label='PCR pooled screening', linestyle='dashed')
ax_1.plot(SIM_i_alpha_1.day, 1 - SIM_i_alpha_1.S/N, label='PCR pooled screening', linestyle='dotted')
ax_1.plot(SIM_0_alpha.day, 1 - SIM_0_alpha.S/N, label='No screening', linestyle='solid')
ax_1.plot(SIM_antigen_alpha_3.day, 1 - SIM_antigen_alpha_3.S/N, label='Antigen Individual screening', linestyle='dashdot',color ='magenta')
ax_1.plot(SIM_antigen_alpha_14.day, 1 - SIM_antigen_alpha_14.S/N, label='Antigen Individual screening', linestyle='dashdot',color ='red')
ax_1.set_title('(A) Cumulative infection rates')
ax_1.set_ylabel('Fraction')

ax_2.plot(SIM_g_alpha_2.day,(SIM_g_alpha_2.I+SIM_g_alpha_2.Q+SIM_g_alpha_2.SQ)/N, label='PCR pooled screening', linestyle='dashed')
ax_2.plot(SIM_i_alpha_1.day, (SIM_i_alpha_1.I+SIM_i_alpha_1.Q+SIM_i_alpha_1.SQ)/N, label='PCR individual screening', linestyle='dotted')
ax_2.plot(SIM_0_alpha.day, (SIM_0_alpha.I+SIM_0_alpha.Q+SIM_0_alpha.SQ)/N, label='No screening', linestyle='solid')
ax_2.plot(SIM_antigen_alpha_3.day, (SIM_antigen_alpha_3.I+SIM_antigen_alpha_3.Q+SIM_antigen_alpha_3.SQ)/N, label='Antigen screening (every 3 days)', linestyle='dashdot',color ='magenta')
ax_2.plot(SIM_antigen_alpha_14.day, (SIM_antigen_alpha_14.I+SIM_antigen_alpha_14.Q+SIM_antigen_alpha_14.SQ)/N, label='Antigen screening (every 14 days)', linestyle='dashdot',color ='red')

ax_2.set_ylabel('Fraction')
ax_2.set_title('(B) Infection rates')

ax_1.set_xlim(-1, 390)
ax_2.set_xlim(-1, 390)

ax_2.legend(ncol=1)
plt.subplots_adjust(wspace=0.27, hspace=0.4)

plt.savefig('.//p_start'+str(p_start)+'/'+'epi_'+str(p_start)+'.pdf', bbox_inches='tight')
