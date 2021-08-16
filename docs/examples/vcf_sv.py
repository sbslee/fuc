from fuc import pyvcf, common
import matplotlib.pyplot as plt
import seaborn as sns

sns.set()

cyp2d6_starts = [42522500,42522852,42523448,42523843,42524175,42524785,42525034,42525739,42526613]
cyp2d6_ends = [42522754,42522994,42523636,42523985,42524352,42524946,42525187,42525911,42526883]
cyp2d7_starts = [42536213,42536565,42537161,42537543,42537877,42538479,42538728,42539410,42540284]
cyp2d7_ends = [42536467,42536707,42537349,42537685,42538054,42538640,42538881,42539582,42540576]

common.load_dataset('pyvcf')
vcf_file = '~/fuc-data/pyvcf/getrm-cyp2d6-vdr.vcf'
vf = pyvcf.VcfFrame.from_file(vcf_file)

fig, axes = plt.subplots(2, 3, figsize=(18, 12))

[[ax1, ax2, ax3], [ax4, ax5, ax6]] = axes

vf.plot_region('NA18973', ax=ax1, color='tab:green')
vf.plot_region('HG00276', ax=ax2, color='tab:green')
vf.plot_region('NA19109', ax=ax3, color='tab:green')

vf.plot_region('NA18973', ax=ax4, k='#AD_FRAC_REF', label='REF')
vf.plot_region('NA18973', ax=ax4, k='#AD_FRAC_ALT', label='ALT')
vf.plot_region('HG00276', ax=ax5, k='#AD_FRAC_REF')
vf.plot_region('HG00276', ax=ax5, k='#AD_FRAC_ALT')
vf.plot_region('NA19109', ax=ax6, k='#AD_FRAC_REF')
vf.plot_region('NA19109', ax=ax6, k='#AD_FRAC_ALT')

ax1.set_title('NA18973 (no structural variation)', fontsize=25)
ax2.set_title('NA10831 (CYP2D6 deletion)', fontsize=25)
ax3.set_title('NA19109 (CYP2D6 duplication)', fontsize=25)

ax4.set_ylabel('Allele fraction')
ax4.legend(loc='upper left', fontsize=20, markerscale=2)

for ax in [ax1, ax2, ax3]:
    ax.set_xlabel('')
    ax.set_xticklabels([])
    ax.set_ylim([-15, 120])
    common.plot_exons(cyp2d6_starts, cyp2d6_ends, name='CYP2D6', offset=8, y=-5, height=5, ax=ax)
    common.plot_exons(cyp2d7_starts, cyp2d7_ends, name='CYP2D7', offset=8, y=-5, height=5, ax=ax)

for ax in [ax2, ax3, ax5, ax6]:
    ax.set_ylabel('')
    ax.set_yticklabels([])

for ax in [ax1, ax4, ax5, ax6]:
    ax.xaxis.label.set_size(20)
    ax.yaxis.label.set_size(20)
    ax.tick_params(axis='both', which='major', labelsize=15)

for ax in [ax4, ax5, ax6]:
    ax.set_ylim([-0.15, 1.1])
    common.plot_exons(cyp2d6_starts, cyp2d6_ends, name='CYP2D6', offset=0.075, y=-0.05, height=0.05, ax=ax)
    common.plot_exons(cyp2d7_starts, cyp2d7_ends, name='CYP2D7', offset=0.075, y=-0.05, height=0.05, ax=ax)
    
plt.tight_layout()
plt.savefig('vcf_sv.png')
