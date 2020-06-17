#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.lines import Line2D

embed_file = sys.argv[1]
embedd = pd.read_csv(embed_file, sep='\t', comment=None)
legend_annot = embedd[['Color', 'Label']].drop_duplicates().dropna()


lines = [Line2D([0], [0], color=c, linewidth=0, marker='o') for c in legend_annot['Color'].values]
labels = legend_annot['Label'].values

if sys.argv[3] == 'slingshot':
    embedd_orig = embedd[embedd['Type'].values == 'CB']
    embedd_tree = embedd[embedd['Type'].values != 'CB']
    
    fig = plt.figure(figsize=[10, 10])
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(embedd_orig['Dim.1'].values, embedd_orig['Dim.2'].values, embedd_orig['Dim.3'].values, c=embedd_orig['Color'].values, s=8, marker='.')

    temp_vals = embedd_tree.Type.unique()
    
    for val in temp_vals:
        local_embedd_tree = embedd_tree[embedd_tree['Type'].values == val]
        print(local_embedd_tree.head())
        ax.plot(local_embedd_tree['Dim.1'].values, local_embedd_tree['Dim.2'].values, local_embedd_tree['Dim.3'].values, c='black')
    ax.set_xlabel('Latent Dim 1')
    ax.set_ylabel('Latent Dim 2')
    ax.set_zlabel('Latent Dim 3')
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.legend(lines, labels, loc='upper right')
    plt.savefig(sys.argv[2], dpi=600)

else:
    fig = plt.figure(figsize=[10, 10])
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(embedd['Dim.1'].values, embedd['Dim.2'].values, embedd['Dim.3'].values, c=embedd['Color'].values, s=8, marker='.')
    ax.set_xlabel('Latent Dim 1')
    ax.set_ylabel('Latent Dim 2')
    ax.set_zlabel('Latent Dim 3')
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.legend(lines, labels, loc='upper right')
    plt.savefig(sys.argv[2], dpi=600)
