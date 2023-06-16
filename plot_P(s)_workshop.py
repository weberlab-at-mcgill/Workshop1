"""
@author: lucasphilipp
a modification of:
https://cooltools.readthedocs.io/en/latest/notebooks/contacts_vs_distance.html
"""
import matplotlib.pyplot as plt
import pandas as pd
import cooler
import cooltools
import numpy as np

##############################################################################

#input the resolution of the HiC matrix in bp:
resolution  = 5000
#input the name of your organism here
organism = 'organism'
#choose a file path to save the data, end with a /
filepath = '/Users/lucasphilipp/Downloads/'

#copy path name of the .mcool file for your organism, at the end of this string add "::resolutions/5000"
clr = cooler.Cooler('/Users/lucasphilipp/Downloads/HiC Across The Tree of Life/Acropora millepora stony coral branch fragments/GSM5182734_amil_sf_1.1_HiC.mcool::resolutions/5000')

df = pd.DataFrame({'chrom': clr.chromnames,
'start': 0,
'end': clr.chromsizes.values,
'name': clr.chromnames}
)

cvd_smooth = cooltools.expected_cis(
clr=clr,
view_df=df,
smooth=True,
aggregate_smoothed=True,
nproc=8 #if you do not have multiple cores available, set to 1
)

#columns are chr, chr, s, Ps_smoothed
#print(cvd_smooth_minutum.values[:, [0,1,2,8]])
#P(s) curve for raw counts: count.avg, P(s) curve for normalized counts balanced.avg

#choose a file path to save the P(s) curve as a .csv file
pd.DataFrame(cvd_smooth.values[:, [0,1,2,8]]).to_csv(filepath + '_' + organism + '_P(s).csv', sep='\t')

cvd_smooth['s_bp'] = cvd_smooth['dist'] * resolution
cvd_smooth['balanced.avg.smoothed'].loc[cvd_smooth['dist'] < 2] = np.nan

#this loop displays the raw P(s) curves and smoothed P(s) curves for each chromosome
for region in df['name']:
    f, ax = plt.subplots(1,1)
    ax.loglog(
    cvd_smooth['s_bp'].loc[cvd_smooth['region1']==region],
    cvd_smooth['balanced.avg'].loc[cvd_smooth['region1']==region],
    cvd_smooth['s_bp'].loc[cvd_smooth['region1']==region],
    cvd_smooth['balanced.avg.smoothed'].loc[cvd_smooth['region1']==region]
    )
    ax.set(
    xlabel='s, bp',
    ylabel='P(s)') 
    ax.set_aspect(1.0)
    ax.grid(lw=0.5)
    plt.savefig(filepath + "HIC_SCAFFOLD_{0}".format(region))
    plt.show()
    
#this loop displays the smoothed P(s) curve and its log derivative: d(log(P(s)))/d(log(s)), both P(s) and s are smoothed with averaging done over logarithmically spaced bins
f, (ax, ax2) = plt.subplots(1,2, figsize=(12, 12))
for region in df['name']:
    ax.loglog(
    cvd_smooth['s_bp'].loc[cvd_smooth['region1']==region],
    cvd_smooth['balanced.avg.smoothed'].loc[cvd_smooth['region1']==region],
    alpha=0.3, color='orange', label = organism
    )
handles, labels = ax.get_legend_handles_labels()
display = (0,len(clr.chromnames))
ax.legend([handle for i,handle in enumerate(handles) if i in display],
      [label for i,label in enumerate(labels) if i in display], loc = 'best')
ax.set(
xlabel='s, bp',
ylabel='P(s)')
ax.grid(lw=0.5)
ax.set_aspect(1.465)

plt.xlim( (10**3,10**8) )
plt.ylim( (10**-6,10**-1) )

for region in df['name']:
    y=cvd_smooth['balanced.avg.smoothed'].loc[cvd_smooth['region1']==region]
    x=cvd_smooth['s_bp'].loc[cvd_smooth['region1']==region]
    logy=np.log(y)
    logx=np.log(x)
    grad=np.gradient(logy,logx)
    # Calculate derivative in log-log space
    ax2.semilogx(
    cvd_smooth['s_bp'].loc[cvd_smooth['region1']==region],
    grad, alpha=0.3, color='orange', label = organism
    )
handles, labels = ax.get_legend_handles_labels()
display = (0,len(clr.chromnames))   
ax2.legend([handle for i,handle in enumerate(handles) if i in display],
      [label for i,label in enumerate(labels) if i in display], loc = 'best')
ax2.set(
xlabel='s, bp',
ylabel='slope')
ax2.grid(lw=0.5)
ax2.set_aspect(1.0)

plt.xlim( (10**3,10**8) )
plt.ylim( (-3,2) )
plt.show()


    
    