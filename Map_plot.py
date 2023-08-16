import numpy as np
import os,sys,argparse
import pandas as pd
import time
import glob
import matplotlib as plt

def map_plot(RPKM_file, output):
    df_RPM_filtered_length_filtered_CDS = pd.read_csv(RPKM_file, delim_whitespace=True)
    fig,ax = plt.subplots()
    test = df_RPM_filtered_length_filtered_CDS[(df_RPM_filtered_length_filtered_CDS['Rest_clip_RPM_average'] > 1) & (df_RPM_filtered_length_filtered_CDS['Tg1_clip_RPM_average'] > 1)]
    test['clip_FC_percent'] = (test['clipseq_FC'] - test['clipseq_FC'].min()) / (test['clipseq_FC'].max() - test['clipseq_FC'].min())
    test.plot(kind='scatter', x='TotalRNA_FC', y='polysome_FC', figsize=(15,15), c='clip_FC_percent',cmap='coolwarm', alpha=0.9, ax=ax)
    cbar = ax.collections[0].colorbar
    cbar.set_ticks([0, 0.5, 1])
    cbar.set_ticklabels(['Low\nBinding', 'Average\nBinding', 'High\nBinding'])
    cbar.set_label('')
    plt.axhline(y=0, color='grey', linestyle='-')
    plt.axhline(y=0.58, color='grey', linestyle='-')
    plt.axhline(y=-0.58, color='grey', linestyle='-')
    plt.axvline(x=0, color='grey', linestyle='-')
    plt.axvline(x=0.58, color='grey', linestyle='-')
    plt.axvline(x=-0.58, color='grey', linestyle='-')
    plt.plot([-4,4],[-4,4], 'grey', linewidth=1)
    plt.plot([-3.42,4],[-4,3.42], 'grey', linewidth=1)
    plt.plot([-4,3.42],[-3.42,4], 'grey', linewidth=1)
    plt.ylim(-4, 4)
    plt.xlim(-4, 4)
    ax.set_xticks([-4,-3,-2, -1, 0, 1, 2, 3, 4])
    plt.rc('xtick',labelsize=16)
    plt.rc('ytick',labelsize=16)
    ax.set_ylabel('Log2 fold change of polysome associated (RPM)', fontsize=22)
    ax.set_xlabel('Log2 fold change of RNAseq (RPM)', fontsize=22)
    plt.title('polysome vs RNAseq RPM\n', fontsize=28)
    fig = ax.get_figure()
    fig.savefig(output, dpi=100)def getArgs():
    parser = argparse.ArgumentParser('python')
    parser.add_argument('-RPKM_file', required=True)
    parser.add_argument('-output', required=True)
    return parser.parse_args()

if __name__ == "__main__":
    args = getArgs()
    start = time.time()
    filename = map_plot(args.RPKM_file, args.output)
    end = time.time()
    print ('time elapsed:' + str(end - start))
