# File: customized_oncoplot_3.py

import matplotlib.pyplot as plt
from fuc import common, pymaf
common.load_dataset('tcga-laml')
mf = pymaf.MafFrame.from_file('~/fuc-data/tcga-laml/tcga_laml.maf.gz')
af = common.AnnFrame.from_file('~/fuc-data/tcga-laml/tcga_laml_annot.tsv', sample_col=0)
af.df['days_to_last_followup'] = common.convert_num2cat(af.df['days_to_last_followup'])

# Filter the MafFrame.
mf = mf.filter_annot(af, "Overall_Survival_Status == 1")

# Define the shared variables.
count=10
figsize=(15, 10)
label_fontsize=13
ticklabels_fontsize=12
legend_fontsize=12

# Create the figure. Getting the right height ratios can be tricky and often requires a trial-and-error process.
fig, axes = plt.subplots(6, 2, figsize=figsize, gridspec_kw={'height_ratios': [1, 10, 1, 1, 1, 3.5], 'width_ratios': [10, 1]})
[[ax1, ax2], [ax3, ax4], [ax5, ax6], [ax7, ax8], [ax9, ax10], [ax11, ax12]] = axes
fig.suptitle('customized_oncoplot_3.py', fontsize=20)

# Create the TMB plot.
samples = list(mf.matrix_waterfall(count).columns)
mf.plot_tmb(ax=ax1, samples=samples)
ax1.set_xlabel('')
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax1.spines['bottom'].set_visible(False)
ax1.set_ylabel('TMB', fontsize=label_fontsize)
ax1.set_yticks([0, mf.matrix_tmb().sum(axis=1).max()])
ax1.tick_params(axis='y', which='major', labelsize=ticklabels_fontsize)

ax2.remove()

# Create the waterfall plot.
mf.plot_waterfall(count=count, ax=ax3, linewidths=1, samples=samples)
ax3.set_xlabel('')
ax3.tick_params(axis='y', which='major', labelrotation=0, labelsize=ticklabels_fontsize)

# Create the genes plot.
mf.plot_genes(count=count, ax=ax4, mode='samples', width=0.95)
ax4.spines['right'].set_visible(False)
ax4.spines['left'].set_visible(False)
ax4.spines['top'].set_visible(False)
ax4.set_yticks([])
ax4.set_xlabel('Samples', fontsize=label_fontsize)
ax4.set_xticks([0, mf.matrix_genes(count=10, mode='samples').sum(axis=1).max()])
ax4.set_ylim(-0.5, count-0.5)
ax4.tick_params(axis='x', which='major', labelsize=ticklabels_fontsize)

# Create the annotation plot for 'FAB_classification'.
_, handles1 = af.plot_annot('FAB_classification', samples=samples, ax=ax5, colors='Dark2', xticklabels=False)
ax5.set_ylabel('')
ax6.remove()

# Create the annotation plot for 'days_to_last_followup'.
_, handles2 = af.plot_annot('days_to_last_followup', samples=samples, ax=ax7, colors='viridis', sequential=True, xticklabels=False)
ax7.set_ylabel('')

ax8.remove()

# Create the annotation plot for 'Overall_Survival_Status'.
_, handles3 = af.plot_annot('Overall_Survival_Status', samples=samples, ax=ax9, colors='Pastel1', xticklabels=False)
ax9.set_xlabel('Samples', fontsize=label_fontsize)
ax9.set_ylabel('')

ax10.remove()

# Create the legends. Getting the right legend locations can be tricky and often requires a trial-and-error process.
handles4 = common.legend_handles(pymaf.NONSYN_NAMES + ['Multi_Hit'], pymaf.NONSYN_COLORS + ['k'])
leg4 = ax11.legend(handles=handles4, loc=(0, 0), title='Variant_Classification', ncol=2, fontsize=legend_fontsize, title_fontsize=legend_fontsize)
leg1 = ax11.legend(handles=handles1, loc=(0.43, 0), title='FAB_classification', ncol=2, fontsize=legend_fontsize, title_fontsize=legend_fontsize)
leg2 = ax11.legend(handles=handles2, loc=(0.62, 0), title='days_to_last_followup', fontsize=legend_fontsize, title_fontsize=legend_fontsize)
leg3 = ax11.legend(handles=handles3, loc=(0.82, 0), title='Overall_Survival_Status', fontsize=legend_fontsize, title_fontsize=legend_fontsize)
ax11.add_artist(leg1)
ax11.add_artist(leg2)
ax11.add_artist(leg3)
ax11.add_artist(leg4)
ax11.axis('off')

ax12.remove()

plt.tight_layout()
plt.subplots_adjust(wspace=0.01, hspace=0.01)
plt.savefig('customized_oncoplot_3.png')
