#!/usr/bin/env python3
"""Create Figure 5 from the reproducible BFAST C3D1 landmark analysis.

All data are read from the original local CSVs and used only in memory. The
output figure contains aggregate curves and risk counts, not subject IDs.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from statsmodels.duration.survfunc import SurvfuncRight
from bfast_c3d1_landmark_analysis import prepare_data, estimate

C={'ALK cleared':'#24679A', 'Detectable ALK':'#B84C54'}
LABELS={'ALK cleared':'ALK fusion cleared at C3D1 (n=66)',
        'Detectable ALK':'ALK fusion detectable at C3D1 (n=5)'}

def make_plot(data,stats,out):
    plt.rcParams.update({'font.size':10,'font.family':'DejaVu Sans','axes.linewidth':0.8,
                         'savefig.facecolor':'white'})
    fig=plt.figure(figsize=(10,7.5),facecolor='white')
    ax=fig.add_axes([0.20,0.41,0.76,0.50])
    risk=fig.add_axes([0.20,0.105,0.76,0.225],sharex=ax)
    max_t=max(float(data['post_landmark_months'].max()),54.)
    xcap=float(np.ceil(max_t/6.)*6.)
    ticks=np.arange(0,xcap+0.01,12.)
    for group in ('ALK cleared','Detectable ALK'):
        cur=data.loc[data['clearance_group']==group]
        x=cur['post_landmark_months'].to_numpy()
        e=cur['pfs_event'].to_numpy()
        sf=SurvfuncRight(x,e)
        # Extend KM survival step through final censoring time. Previous plot
        # stopped at the final event even when patients remained under observation.
        times=np.r_[0,sf.surv_times,float(np.max(x))]
        probs=np.r_[1.,sf.surv_prob,sf.surv_prob[-1]]
        ax.step(times,probs,where='post',lw=2.4,color=C[group],label=LABELS[group])
        for t in x[e==0]:
            idx=np.searchsorted(sf.surv_times,t,side='right')-1
            survival=float(sf.surv_prob[idx]) if idx>=0 else 1.0
            ax.plot([t],[survival],marker='|',ms=7,mew=1.1,color=C[group],linestyle='None')
    ax.set_xlim(0,xcap)
    ax.set_ylim(0,1.02)
    ax.set_yticks(np.arange(0,1.01,0.2))
    ax.tick_params(axis='x',which='both',labelbottom=False,length=3)
    ax.set_ylabel('Progression-free survival probability',labelpad=11,fontsize=11)
    ax.grid(axis='y',color='#D8DDE3',alpha=.75,lw=.65)
    ax.set_axisbelow(True)
    ax.legend(loc='upper right',fontsize=10,frameon=False,handlelength=2.6)
    ax.set_title('C3D1 landmark progression-free survival',fontsize=13,fontweight='bold',loc='left',pad=16)
    ax.text(.02,.04, 'Unadjusted HR 0.17 (95% CI 0.06–0.45); log-rank P<0.001',
            transform=ax.transAxes,color='#30343A',fontsize=9.5,ha='left',va='bottom')

    risk.axis('off')
    risk.set_xlim(0,xcap)
    risk.text(-.18,.96,'Number at risk',fontweight='bold',fontsize=10.3,
              transform=risk.transAxes,ha='left',va='center')
    for i,group in enumerate(('ALK cleared','Detectable ALK')):
        cur=data.loc[data['clearance_group']==group]
        y=.66-i*.32
        risk.text(-.018,y,'Cleared' if group=='ALK cleared' else 'Detectable',
                  color=C[group],ha='right',va='center',transform=risk.transAxes,fontsize=10)
        for t in ticks:
            n=int((cur['post_landmark_months'] >= t).sum())
            risk.text(t,y,str(n),ha='center',va='center',fontsize=10)
    for t in ticks:
        risk.text(t,.065,str(int(t)),ha='center',va='top',fontsize=10)
    risk.text(.5,-.19,'Months after C3D1 landmark',ha='center',va='top',fontsize=11,transform=risk.transAxes)
    fig.savefig(out,dpi=300,bbox_inches='tight',pad_inches=.12)
    plt.close(fig)

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--data-dir',type=Path,default=Path('/mnt/data'))
    parser.add_argument('--output',type=Path,default=Path('/mnt/data/Figure_5_BFAST_C3D1_landmark_PFS.png'))
    args=parser.parse_args()
    data,audit=prepare_data(args.data_dir/'adtte.csv',
                           args.data_dir/'bfast_oSpl_cohortA.withDBIannotations.BlC3Td.csv',
                           args.data_dir/'bfast_oapl_cohortA.withDBIannotations.clean.BlC3Td.csv',56)
    stats=estimate(data)
    assert audit['eligible_landmark_n']==71 and audit['eligible_PFS_events']==47
    assert abs(stats['hazard_ratio']-.1676)<.001
    args.output.parent.mkdir(parents=True,exist_ok=True)
    make_plot(data,stats,args.output)
    print(args.output, args.output.stat().st_size)
