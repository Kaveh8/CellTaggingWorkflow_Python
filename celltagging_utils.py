import numpy as np
import pandas as pd
from scipy.spatial.distance import jaccard
from tqdm.notebook import tqdm
import networkx as nx
from itertools import product, chain
import matplotlib.pyplot as plt
import seaborn as sns
from random import sample
import re
from pyvis.network import Network
import datetime
#********************************************************************************
def single_cell_data_binarization(celltag_dat, tag_cutoff):
    
    celltags = celltag_dat.copy()
    celltags[celltags < tag_cutoff] = 0
    celltags[celltags > 0] = 1

    return celltags
#********************************************************************************
def get_stats(df):
    df1 = pd.DataFrame({'Cell_UMI_Counts':[],'CellTags_per_Cell':[]})
    df2 = pd.DataFrame({'CellTag_UMI_Counts':[],'Cells_per_CellTag':[]})
    df1['Cell_UMI_Counts'] = df.sum(axis=1).astype(int)
    df1['CellTags_per_Cell'] = df.gt(0).sum(axis=1).astype(int)
    df2['CellTag_UMI_Counts'] = df.T.sum(axis=1).astype(int)
    df2['Cells_per_CellTag'] = df.T.gt(0).sum(axis=1).astype(int)
    return pd.concat([df1['Cell_UMI_Counts'].describe(), 
                      df1['CellTags_per_Cell'].describe(), 
                      df2['CellTag_UMI_Counts'].describe(), 
                      df2['Cells_per_CellTag'].describe()], axis=1)
#********************************************************************************                     
def single_cell_data_whitelist(celltag_dat, whitels_cell_tag_file):
    
    celltags = celltag_dat.copy().T
    cell_names = celltag_dat.index
    tag_names = celltag_dat.columns
    whitelist_names = pd.read_csv(whitels_cell_tag_file, header=0, dtype=str).iloc[:, 0]
    
    whitelist_intersection = whitelist_names[whitelist_names.isin(celltags.index)]
    celltags_whitelisted = celltags.loc[whitelist_intersection, :].T
    
    return celltags_whitelisted
#********************************************************************************
def metric_based_filtering(whitelisted_celltag_data, cutoff, comparison="less"): 
    celltags_per_cell_whitelisted_pf = pd.DataFrame(whitelisted_celltag_data.sum(axis=1), columns=["Sum"])
    if comparison == "less":
        cell_filter = celltags_per_cell_whitelisted_pf[celltags_per_cell_whitelisted_pf["Sum"] < (cutoff + 1)]
    else:
        cell_filter = celltags_per_cell_whitelisted_pf[celltags_per_cell_whitelisted_pf["Sum"] > (cutoff - 1)]

    cell_bc_filter = cell_filter.index
    celltags_whitelisted = whitelisted_celltag_data.T
    celltags_whitelisted_new = celltags_whitelisted[cell_bc_filter].T

    return celltags_whitelisted_new
    
    
#********************************************************************************    
def jaccard_analysis(whitelisted_celltag_data, save_mtx=True, id=""):

    num_rows = whitelisted_celltag_data.shape[0]
    jaccard_similarities = np.zeros((num_rows, num_rows))
    
    for i in tqdm(range(num_rows), desc="Calculating Jaccard Similarities"):
        for j in range(i + 1, num_rows):
            similarity = 1 - jaccard(whitelisted_celltag_data.iloc[i], whitelisted_celltag_data.iloc[j])
            jaccard_similarities[i, j] = similarity
            jaccard_similarities[j, i] = similarity
            
    np.fill_diagonal(jaccard_similarities, 1)  
    Jac = pd.DataFrame(jaccard_similarities, index=whitelisted_celltag_data.index, columns=whitelisted_celltag_data.index)
    
    if save_mtx:
        Jac.to_csv(f"{id}_Jaccard_mtx.csv", index_label="Cell_BC")

    return Jac
#********************************************************************************    
def find_values_above_cutoff(dataframe, cutoff):
    result = []
    for col in dataframe.columns:
        rows = dataframe.index[dataframe[col] > cutoff].tolist()
        for row in rows:
            result.append((row, col))
    return result
#********************************************************************************    
def clone_calling(jaccard_df, output_file, correlation_cutoff):

    edges = find_values_above_cutoff(jaccard_df, correlation_cutoff)
    
    G = nx.Graph()
    G.add_edges_from(edges)

    clones = [clone for clone in nx.connected_components(G) if len(clone) > 1]
    
    # Create clone membership table
    clone_membership = []
    for clone_id, clone in enumerate(clones, start=1):
        for cell in clone:
            clone_membership.append({'clone_id': clone_id, 'cell_barcode': cell})
    clone_df = pd.DataFrame(clone_membership)
    
    clone_df.to_csv(output_file, index=False, quoting=0)
    
    clone_sizes = clone_df['clone_id'].value_counts().reset_index()
    clone_sizes.columns = ['Clone_ID', 'Frequency']
    
    return clone_df, clone_sizes
    
#********************************************************************************   
def convert_cell_tag_matrix_to_link_list(celltag_data):
    
    print("Preprocessing data..")
    # Remove cells that do not have any celltag
    celltag_data.dropna(how='all', inplace=True)
    print("Cells that have CellTagV1:", celltag_data['CellTagV1'].notna().sum())
    print("Cells that have CellTagV2:", celltag_data['CellTagV2'].notna().sum())
    print("Cells that have CellTagV3:", celltag_data['CellTagV3'].notna().sum())
    celltag_data.fillna("e", inplace=True)
    
    def find_root(cell_id, tag):
        tagid = celltag_data.at[cell_id, tag]
        tmp = pd.DataFrame([f"{tag}_{tagid}", cell_id, tag]).T
        tmp.columns = ["source", "target", "tag"]
        return tmp
#------------------------------- 
    print("find connection between [celltag -> cells]...") 
    # 2.1 find connection between "celltag" -> "cells"
    tags = ["CellTagV3", "CellTagV2", "CellTagV1"]
    linkList = pd.DataFrame(columns=["source", "target", "tag"])
    cell_ids = celltag_data.index
    for tag in tags:
        subcells = celltag_data[celltag_data[tag] != "e"]    
        tmp = pd.concat([find_root(i, tag) for i in subcells.index], ignore_index=True)
        linkList = pd.concat([linkList, tmp], ignore_index=True)
#------------------------------- 
    print("find hidden links [CellTagV2 -> CellTagV3], or [CellTagV1 -> CellTagV3]...") 
    # 2.2 hidden link ["CellTagV2" -> "CellTagV3"], or ["CellTagV1" -> "CellTagV3"]
    hiddenlink_D13 = pd.DataFrame(columns=["source", "target", "tag"])
    for i in celltag_data['CellTagV3'].unique():
        if pd.isna(i) or i == 'e':
            continue
        
        sub_cells = celltag_data[celltag_data['CellTagV3'] == i]
        
        for prev_tag_name in ['CellTagV2', 'CellTagV1']:
            prev_tag = sub_cells[prev_tag_name]
            prev_tag = prev_tag[prev_tag != 'e']
            
            if not prev_tag.empty:
                prev_tag = prev_tag.mode()[0]
                tmp = pd.DataFrame([[f"{prev_tag_name}_{prev_tag}", f"CellTagV3_{i}", prev_tag_name]], columns=["source", "target", "tag"])
                hiddenlink_D13 = pd.concat([hiddenlink_D13, tmp], ignore_index=True)
                break
#------------------------------- 
    print("find hidden links [CellTagV1 -> CellTagV2]...") 
    # 2.3 hidden link ["CellTagV1" -> "CellTagV2"]
    hiddenlink_D3 = pd.DataFrame(columns=["source", "target", "tag"])
    for i in celltag_data['CellTagV2'].unique():
        if pd.isna(i) or i == 'e':
            continue
        
        sub_cells = celltag_data[celltag_data['CellTagV2'] == i]
        
        prev_tag_name = 'CellTagV1'
        prev_tag = sub_cells[prev_tag_name]
        prev_tag = prev_tag[prev_tag != 'e']
        
        if not prev_tag.empty:
            prev_tag = prev_tag.mode()[0]
            tmp = pd.DataFrame([[f"{prev_tag_name}_{prev_tag}", f"CellTagV2_{i}", prev_tag_name]], columns=["source", "target", "tag"])
            hiddenlink_D3 = pd.concat([hiddenlink_D3, tmp], ignore_index=True)
#------------------------------- 
    def modify_cell_name(linkList):
        linkList['target_unmodified'] = linkList['target']
        node_cell = linkList['target'].str.contains('-')
        linkList.loc[node_cell, 'target'] += '_' + linkList.loc[node_cell, 'tag'].str.split('g').str[1]
        return linkList
    
    # 2.4 integrating all links
    linkList = pd.concat([linkList, hiddenlink_D3, hiddenlink_D13], ignore_index=True)
    linkList = modify_cell_name(linkList)
    
    print("finished")
    return linkList
#********************************************************************************
def get_nodes_from_link_list(link_list):  
    nodes = pd.unique(link_list[['target', 'source']].values.ravel('K'))
    nodes = pd.DataFrame(nodes, columns=['nodes'])
    
    def get_type(each_node):
        is_cell = sum(tag in each_node['nodes'] for tag in ["CellTagV1", "CellTagV2", "CellTagV3"]) == 0
        if is_cell:
            return "cell"
        else:
            return "clone"
    
    def get_tag_id(each_node):
        if each_node['type']=="cell":
            return link_list.loc[link_list['target'] == each_node['nodes'], 'tag'].values[0]
        else:
            return each_node['nodes'].split("_")[0]
        
    def get_unmodified_name(each_node):
        if each_node['type']=="cell":
            return link_list.loc[link_list['target'] == each_node['nodes'], 'target_unmodified'].values[0]
        else:
            return each_node['nodes']
        
    def get_clone_no(each_node):
        if each_node['type']=="cell":
            return link_list.loc[link_list['target'] == each_node['nodes'], 'source'].values[0]
        else:
            return each_node['nodes']
    
    nodes['type'] = nodes.apply(get_type, axis=1)
    nodes['tag'] = nodes.apply(get_tag_id, axis=1)
    nodes['node_name_unmodified'] = nodes.apply(get_unmodified_name, axis=1)
    nodes['clone'] = nodes.apply(get_clone_no, axis=1)
    
    return nodes
#********************************************************************************
def add_data_to_nodes(Nodes, additional_data):
    return pd.merge(Nodes, additional_data, left_on='node_name_unmodified', right_index=True, how='left')

#********************************************************************************
def return_directly_connected_nodes(node, link_list):
    tmp_link = link_list[link_list['source'].isin(node)]
    tmp_link2 = link_list[link_list['target'].isin(node)]
    
    tmp_nodes = set(tmp_link['target']).union(set(tmp_link2['source'])).union(set(node))
    return list(tmp_nodes)

#********************************************************************************

#********************************************************************************
def get_subnet(ref_nodes, link_list, nodes):
    connected_nods = return_directly_connected_nodes(ref_nodes, link_list)
    sub_link = link_list[link_list['source'].isin(connected_nods) | link_list['target'].isin(connected_nods)]
    sub_nodes = nodes[nodes['nodes'].isin(connected_nods)]
    
    return sub_link, sub_nodes
    