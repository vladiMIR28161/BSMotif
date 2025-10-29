#!/usr/bin/env python3

import pandas as pd
import numpy as np
import re
import warnings
from copy import copy

warnings.filterwarnings("ignore", category=RuntimeWarning)

results_df = pd.DataFrame(columns=['Branch', 'List', 'Similarity'], dtype=str)

def hierarchical_classification_tf (df, query_col, target_col, next_query_col, next_target_col):
    
    df['Branch'] = 0
    
    # combine both superclass columns and get unique values
    unique_levels = pd.concat([df[query_col], df[target_col]]).dropna().unique()
    
    # sort by numeric value inside {}
    unique_levels = sorted(
        unique_levels,
        key=lambda x: tuple(map(int, x.split('{')[1].strip('}').split('.'))))
    
    for lvl in unique_levels:
        
        code_lvl = re.findall(r"\{.*?\}", lvl)[0]
        
        sub_df = df.loc[(df[query_col] == lvl) & (df[target_col] == lvl)].reset_index(drop=True)
        if sub_df.empty:
            continue
            
        next_ids = pd.concat([sub_df[next_query_col], sub_df[next_target_col]]).dropna().unique()
        next_ids = sorted(
            next_ids,
            key=lambda x: tuple(map(int, x.split('{')[1].strip('}').split('.'))))
        
        # calculate pairwise medians
        pairs = []
        for i, id_1 in enumerate(next_ids):
            for id_2 in next_ids[i:]:
                mask = sub_df.loc[(((sub_df[next_query_col] == id_1) & (sub_df[next_target_col] == id_2)) |
                                   ((sub_df[next_query_col] == id_2) & (sub_df[next_target_col] == id_1)))]
                med = np.median(mask.Score_TF.unique()).round(3)
                pairs.append([id_1, id_2, round(med, 3)])
                
        # filter low-median self-similarites
        if next_query_col != 'Query_gene':
            black_list = [x[0] for x in pairs if x[0] == x[1] and x[2] < 3]
            pairs = [x for x in pairs if x[0] not in black_list and x[1] not in black_list]

        list_level = list({item for row in pairs for item in row if not isinstance(item, (int, float))})
        list_level = sorted(list_level, key=lambda x: tuple(map(int, x.split('{')[1].strip('}').split('.'))))
        
        branch = {}
        N = 0
        for x in list_level:
            N += 1
            branch[N] = []
            branch[N].append(x)
        
        # definition of branches
        flag = 0
        while flag == 0:
            max_score = 3
            N = 0
            total_x = '0'
            for x in branch.values():
                N += 1
                N_2 = N
                for x_2 in list(branch.values())[N_2:]:
                    N_2 += 1
                    mask = sub_df[(sub_df[next_query_col].isin(x) | sub_df[next_target_col].isin(x)) & 
                                  (sub_df[next_query_col].isin(x_2) | sub_df[next_target_col].isin(x_2))].reset_index(drop=True)
                    med = np.median(mask.Score_TF.unique()).round(3)
                    if med >= max_score:
                        max_score = med
                        total_x = x
                        total_x_2 = x_2
            if total_x == '0':
                flag = 1
            else:
                for key, value in branch.items():
                    if value == total_x:
                        total_key_1 = key
                    elif value == total_x_2:
                        total_key_2 = key
                branch[total_key_1] = branch[total_key_1] + branch[total_key_2]
                del branch[total_key_2]
                branch = {i+1: v for i, v in enumerate(branch.values())}
        
        print (branch)
        global results_df
        
        if next_query_col != 'Query_gene':        
            if len(branch) == 0:
                print ('Ветвей нет!')
            else:
                for br in branch.items():
                    df_new = df[(df[next_query_col].isin(br[1]) & df[next_target_col].isin(br[1]))].reset_index(drop=True)
                    med = np.median(df_new.Score_TF.unique()).round(3)
                    if np.isnan(med):
                        med = '#'
                    print(br, med)
                    results_df.loc[len(results_df)] = [code_lvl + ' ' + str(br[0]), ', '.join(br[1]), str(med)]
                    df = pd.concat([df, df_new, df_new]).drop_duplicates(keep=False)
                print ('\n')
        
        else:    
            if len(branch) == 0:
                print('Ветвей нет')
            else:
                count = 0
                for br in branch.items():
                    if len (br[1]) == 1:
                        for br2 in br[1]:
                            df_new = df[(df[next_query_col].isin([br2]) & df[next_target_col].isin([br2]))].reset_index(drop=True)
                            med = np.median(df_new.Log10pvalue.unique()).round(3)
                            if np.isnan(med):
                                med = '#'
                            count += 1
                            print(br2, med)
                            df = pd.concat([df, df_new, df_new]).drop_duplicates(keep=False)
                            results_df.loc[len(results_df)] = [code_lvl + ' ' + str(count), br2, str(med)]
                    else:
                        br3 = copy(br[1])
                        for br2 in br3:
                            df_new = df[(df[next_query_col].isin([br2]) & df[next_target_col].isin([br2]))].reset_index(drop=True)
                            med = np.median(df_new.Log10pvalue.unique()).round(3)
                            if med < 3:
                                if np.isnan(med):
                                    med = '#'
                                count += 1
                                print(br2, med)
                                results_df.loc[len(results_df)] = [code_lvl + ' ' + str(count), br2, str(med)]
                                br[1].remove(br2)
                        df_new = df[(df[next_query_col].isin(br[1]) & df[next_target_col].isin(br[1]))].reset_index(drop=True)
                        med = np.median(df_new.Log10pvalue.unique()).round(3)
                        if med >= 3:
                            if np.isnan(med):
                                med = '#'
                            count += 1
                            print(*br[1], med)
                            df = pd.concat([df, df_new, df_new]).drop_duplicates(keep=False)
                            results_df.loc[len(results_df)] = [code_lvl + ' ' + str(count), ', '.join(br[1]), str(med)]
                        else:
                            for br2 in br[1]:
                                df_new = df[(df[next_query_col].isin([br2]) & df[next_target_col].isin([br2]))].reset_index(drop=True)
                                med = np.median(df_new.Log10pvalue.unique()).round(3)
                                if np.isnan(med):
                                    med = '#'
                                count += 1
                                print(br2, med)
                                df = pd.concat([df, df_new, df_new]).drop_duplicates(keep=False)
                                results_df.loc[len(results_df)] = [code_lvl + ' ' + str(count), br2, str(med)]
                print ('\n')
        
    return df