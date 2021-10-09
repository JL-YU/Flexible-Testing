import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':13})
plt.rcParams.update({'axes.labelpad':9})
plt.rc('font', family = 'Times New Roman')
plt.rc('legend', fontsize=11)

N=100000
p_start = 0.001


SIM_0_alpha = pd.read_csv('.//p_start'+str(p_start)+'/'+'exp_no_testing.csv')
SIM_g_alpha_2 = pd.read_csv('.//p_start'+str(p_start)+'/'+'exp_nstar7_round_2delay.csv')
SIM_i_alpha_1 = pd.read_csv('.//p_start'+str(p_start)+'/'+'exp_individual_round_1delay.csv')
SIM_antigen_alpha_3 = pd.read_csv('.//p_start'+str(p_start)+'/'+'exp_antigen_individual_round_3day.csv')
SIM_antigen_alpha_14 = pd.read_csv('.//p_start'+str(p_start)+'/'+'exp_antigen_individual_round_14day.csv')

SIM_0_alpha = SIM_0_alpha.loc[SIM_0_alpha.total_tested_individual>0]
SIM_g_alpha_2 = SIM_g_alpha_2.loc[SIM_g_alpha_2.total_tested_individual>0]
SIM_i_alpha_1 = SIM_i_alpha_1.loc[SIM_i_alpha_1.total_tested_individual>0]
SIM_antigen_alpha_3 = SIM_antigen_alpha_3.loc[SIM_antigen_alpha_3.total_tested_individual>0]
SIM_antigen_alpha_14 = SIM_antigen_alpha_14.loc[SIM_antigen_alpha_14.total_tested_individual>0]

fig = plt.figure(figsize=(15, 9))
ax_1 = fig.add_subplot(3, 2, 1)
ax_2 = fig.add_subplot(3, 2, 2)
ax_3 = fig.add_subplot(3, 2, 3)
ax_4 = fig.add_subplot(3, 2, 4)
ax_5 = fig.add_subplot(3, 2, 5)
ax_6 = fig.add_subplot(3, 2, 6)


ax_1.plot(SIM_g_alpha_2.loc[SIM_g_alpha_2.n_star >1,'day'], SIM_g_alpha_2.loc[SIM_g_alpha_2.n_star>1,'n_star'], label='PCR pooled screening', linestyle='dashed')
ax_1.plot(SIM_i_alpha_1.day, SIM_i_alpha_1.n_star, label='PCR individual screening', linestyle='dotted')
ax_1.plot(SIM_antigen_alpha_3.day, SIM_antigen_alpha_3.n_star, label='Antigen screening (every 3 days)', linestyle='dashdot',color ='magenta')
ax_1.plot(SIM_antigen_alpha_14.day, SIM_antigen_alpha_14.n_star, label='Antigen screening(every 14 days)', linestyle='dashdot',color ='red')
ax_1.set_ylabel('Count')
ax_1.set_ylim(-1, 33)
ax_1.set_title('(A) Optimal pool sizes')

ax_2.plot(SIM_g_alpha_2.day, SIM_g_alpha_2.total_tested_individual.cumsum()/N, label='PCR pooled screening', linestyle='dashed')
ax_2.plot(SIM_i_alpha_1.day, SIM_i_alpha_1.total_tested_individual.cumsum()/N, label='PCR individual screening', linestyle='dotted')
ax_2.plot(SIM_antigen_alpha_3.day, SIM_antigen_alpha_3.total_tested_individual.cumsum()/N, label='Antigen screening (every 3 days)', linestyle='dashdot',color ='magenta')
ax_2.plot(SIM_antigen_alpha_14.day, SIM_antigen_alpha_14.total_tested_individual.cumsum()/N, label='Antigen screening(every 14 days)', linestyle='dashdot',color ='red')
ax_2.set_ylabel('Count')
ax_2.set_title('(B) Total number of screening rounds')


ax_3.plot(SIM_g_alpha_2.day, SIM_g_alpha_2.FN.cumsum(), label='PCR pooled screening', linestyle='dashed')
ax_3.plot(SIM_i_alpha_1.day, SIM_i_alpha_1.FN.cumsum(), label='PCR individual screening', linestyle='dotted')
ax_3.plot(SIM_antigen_alpha_3.day, SIM_antigen_alpha_3.FN.cumsum(), label='Antigen screening (every 3 days)', linestyle='dashdot',color ='magenta')
ax_3.plot(SIM_antigen_alpha_14.day, SIM_antigen_alpha_14.FN.cumsum(), label='Antigen screening (every 14 days)', linestyle='dashdot',color ='red')
ax_3.set_ylabel('Count')
ax_3.set_title('(C) Total number of false negative results')

ax_4.plot(SIM_g_alpha_2.day, SIM_g_alpha_2.FP.cumsum(), label='FT_TC', linestyle='dashed')
ax_4.plot(SIM_i_alpha_1.day, SIM_i_alpha_1.FP.cumsum(), label='IT_TC', linestyle='dotted')
ax_4.plot(SIM_antigen_alpha_3.day, SIM_antigen_alpha_3.FP.cumsum(), label='IT_TC', linestyle='dashdot',color ='magenta')
ax_4.plot(SIM_antigen_alpha_14.day, SIM_antigen_alpha_14.FP.cumsum(), label='IT_TC', linestyle='dashdot',color ='red')
ax_4.set_ylabel('Count')
ax_4.set_title('(D) Total number of false positive results')

ax_5.scatter(SIM_g_alpha_2.day, SIM_g_alpha_2.NPV, label='PCR pooled screening',s=3, marker = "s" )
ax_5.scatter(SIM_i_alpha_1.day, SIM_i_alpha_1.NPV, label='PCR individual screening',s=3,)
ax_5.scatter(SIM_antigen_alpha_3.day, SIM_antigen_alpha_3.NPV, label='Antigen screening (every 3 days)',s=3,color ='magenta',marker = "<")
ax_5.scatter(SIM_antigen_alpha_14.day, SIM_antigen_alpha_14.NPV, label='Antigen screening(every 14 days)',s=3,color ='red',marker = "d")
ax_5.set_ylabel('Probability')
ax_5.set_title('(E) Negative predictive values')

# ax_6.scatter(SIM_g_alpha_1.day, SIM_g_alpha_1.PPV, label='FT',s=3)
ax_6.scatter(SIM_g_alpha_2.day, SIM_g_alpha_2.PPV, label='FT',s=3, marker = "s")
ax_6.scatter(SIM_i_alpha_1.day, SIM_i_alpha_1.PPV, label='IT',s=3)
ax_6.scatter(SIM_antigen_alpha_3.day, SIM_antigen_alpha_3.PPV, label='IT',s=3,color ='magenta',marker = "<")
ax_6.scatter(SIM_antigen_alpha_14.day, SIM_antigen_alpha_14.PPV, label='IT',s=3,color ='red',marker = "d")
ax_6.set_ylabel('Probability')
ax_6.set_title('(F) Positive predictive values')


ax_1.set_xlim(-1, 390)
ax_2.set_xlim(-1, 390)
ax_3.set_xlim(-1, 390)
ax_4.set_xlim(-1, 390)
ax_5.set_xlim(-1, 390)
ax_6.set_xlim(-1, 390)


ax_5.legend(ncol=1)
ax_3.legend(ncol=1)
fig.text(0.5, 0.05,'Time (days)', ha='center',fontsize=15 )
plt.subplots_adjust(wspace=0.27, hspace=0.4)
plt.savefig('.//p_start'+str(p_start)+'/'+'test_'+str(p_start)+'.pdf', bbox_inches='tight')
