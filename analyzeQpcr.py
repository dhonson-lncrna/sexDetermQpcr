import polars as pl
import numpy as np
import argparse

from datetime import datetime
from collections import defaultdict
from itertools import chain, combinations

import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 8
import seaborn as sns

def collect_cts(layout,
                results):
    '''
    Extracts Ct values based on a sample/primer layout file and a qPCR 
    results output.

    Parameters
    ----------
    layout : str
        Path to layout excel file containing Sample and Primer tabs.
    results: str
        Path to results excel file output from $$$ machine.

    Output
    ------
    A dictionary with the structure ct_dict['sample']['primer'] = [ct1, ct2, ...]

    '''
    df = pl.read_excel(layout,
                   sheet_id=0, infer_schema_length=0,
                   columns=[str(i) for i in range(1,25)])
    for k,v in df.items():
        df[k] = v.fill_null('empty')
    
    ct_vals = pl.read_excel(results,
                            sheet_name='Results',
                            columns=['Well Position','CT'],
                            read_options={'header_row':21}
                           )['CT'].to_numpy()
    
    ct_dict = defaultdict(lambda: defaultdict(list))
    
    s_flat = list(chain.from_iterable([i for i in df['Samples'].iter_rows()]))
    p_flat = list(chain.from_iterable([i for i in df['Primers'].iter_rows()]))
    
    s_vals = []
    p_vals = []
    
    for s,p,c in zip(s_flat,p_flat,ct_vals):
        s_vals.append(s)
        p_vals.append(p)
        if s == 'empty':
            continue
        elif (s == 'empty') != (p == 'empty'):
            raise ValueError(f'''Mismatch in sample and primer configuration. 
                Sample was {s} and primer was {p}.''')
        else:
            try:
                ct_dict[s][p].append(float(c))
            except:
                ct_dict[s][p].append(c)
    
    s_vals = set(s_vals)
    p_vals = set(p_vals)
    
    return ct_dict

def analyzer(ct_dict,
             thresh = 0.1,
             prefix = None):
    '''
    Averages Ct values, filters out poor results (high deviation, Ct > 30), 
    and performs ∆∆Ct to determine sex.

    Parameters
    ----------
    ct_dict : dict
        Output of collect_cts().
    thresh : float
        Maximum acceptable standard deviation of Ct replicates. Default: 0.1.
    prefix : str
        Prefix for _failedWells.log and _results.csv output files. If None, 
        uses current date. Default: None.

    Output
    ------
    A polars dataframe containing results.
    '''
    if prefix:
        pass
    else:
        prefix = datetime.now().strftime("%Y%m%d")
    
    logname = f'{prefix}_failedWells.log'
    outfile = f'{prefix}_results.csv'

    log = open(logname,'w')
    mean_dict = defaultdict(defaultdict)
    
    for squid in ct_dict.keys():
        mean_squid = defaultdict(lambda: defaultdict(list))
        
        for prim, ct in ct_dict[squid].items():
            if prim == 'empty':
                continue
            std = np.std(ct)
            if std > thresh:
                diffs = [abs(np.mean(ct) - i) for i in ct]
                ct_new = [v for i,v in enumerate(ct) if i != np.argmax(diffs)]
                std = np.std(ct_new)
                if std > thresh:
                    log.write(f'Failed stdev for {squid} {prim}\n')
                    continue
                else:
                    mean = np.mean(ct_new)
            else:
                mean = np.mean(ct)
            if mean > 30:
                log.write(f'Mean Ct > 30 for {squid} {prim}\n')
                continue
            if 'Z' in prim:
                mean_squid['chrZ']['mean'].append(mean)
                mean_squid['chrZ']['std'].append(std)
            else:
                mean_squid['chrA']['mean'].append(mean)
                mean_squid['chrA']['std'].append(std)
    
        # Check for failures
        if len(mean_squid.keys()) == 1:
            log.write(f'{squid} had only values for {mean_squid.keys()}\n')
            continue
        if not all([len(i['mean']) >= 2 for i in mean_squid.values()]):
            log.write(f'Too many Ct values failed for {squid}\n')
            continue
        
        for prim, val_dict in mean_squid.items():
            mean_all = np.mean(val_dict['mean'])
            std_all = np.sqrt(np.sum([i**2 for i in val_dict['std']])) / len(val_dict['std'])
            mean_dict[squid][prim] = (mean_all,std_all)
    log.close()
    
    delta_ls = []
    for squid, pdict in mean_dict.items():
        delta = pdict['chrZ'][0] - pdict['chrA'][0]
        if delta >= 1:
            sex = 'f'
        elif delta >= 0.5:
            sex = 'ambig'
        else:
            sex = 'm'
        std = 0.5 * np.sqrt(pdict['chrZ'][1]**2 + pdict['chrA'][1]**2)
        delta_ls.append((squid,delta,std,sex))
    
    df = pl.DataFrame(delta_ls, 
                      schema=['squid','mean','stdev','sex'],
                      orient='row').sort('squid')
    df = df.sort('squid')
    df.write_csv(outfile)

    return df

def plot_results(df,
                 prefix = None):
    '''
    Plots results of analysis.

    Parameters
    ----------
    df : a polars dataframe
        Output of analyzer().
    prefix : str
        Prefix for _qpcrResults.png output file. If None, uses current date.
        Default: None.
    '''
    if prefix:
        figname = f'{prefix}_qpcrResults.png'
    else:
        figname = f'{datetime.now().strftime("%Y%m%d")}_qpcrResults.png'
    fig, ax = plt.subplots(1,1,
                           figsize=(2,0.25*len(df)),
                           layout='constrained')
    
    cmap={'m':'#1B76B7',
          'f':'#F58817'}
    for row in df.iter_rows(named=True):
        ax.errorbar(row['mean'], row['squid'],
                    xerr=row['stdev'],
                    fmt='o',
                    color=cmap[row['sex']],
                    ecolor='gray',
                    elinewidth=1,
                    capsize=1.75,
                    capthick=0.75,
                    markersize=2)
    
    ax.margins(y=0.1)
    ax.axvline(1,color='gray',linewidth=1,zorder=0,ls=':')
    ax.set_xticks(np.arange(0.5,2.5,0.5))
    ax.set_xlabel('∆Ct\n(chrZ - chrA)')
    
    plt.savefig(figname,format='png',dpi=600)
    plt.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Analyze qPCR results to determine sex via ∆Ct (chrZ - chrA).'
    )
    parser.add_argument(
        'layout',
        help='Path to the layout Excel file containing Sample and Primer tabs.'
    )
    parser.add_argument(
        'results',
        help='Path to the results Excel file output from the qPCR machine.'
    )
    parser.add_argument(
        '-t', '--thresh',
        type=float,
        default=0.1,
        help='Maximum acceptable standard deviation of Ct replicates. '
             'Wells exceeding this threshold will have their worst outlier removed '
             'before re-evaluation. (default: 0.1)'
    )
    parser.add_argument(
        '-p', '--prefix',
        default=None,
        help='Prefix for output files (<prefix>_failedWells.log, <prefix>_results.csv, '
             '<prefix>_qpcrResults.png). Defaults to the current date (YYYYMMDD).'
    )
    args = parser.parse_args()
 
    ct_dict = collect_cts(args.layout, args.results)
    df = analyzer(ct_dict, thresh=args.thresh, prefix=args.prefix)
    plot_results(df, prefix=args.prefix)