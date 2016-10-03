#!/usr/bin/env python

from __future__ import division
import pandas
import numpy
import os, sys
import trajutils
import matplotlib.pyplot as plt
#Now add the following code in case plotting is done on remote servers to avoid backend x-server related issue
plt.switch_backend('agg')

def plotSIFT(avg_sift=None, residue_labels=None,
             sift_types = ['Apolar','Aro_F2F','Aro_E2F','Hbond_ProD','Hbond_ProA','Elec_ProP','Elec_ProN'],
             color_list = ['r','b','g','y','c','m','k'], filename='AVG_SIFT_BAR.png'):
    '''Bar plot of Avg. contacts coming from SIFT calculation'''
    #Prepare the dataset 'sift_data' for pandas dataframe
    #dataset: is a dict (for residues with non-zero contacts)
    #dict keys 'Residue_Name', followed by entries in sift_types in order
    #NOTE: ENTRIES in sift_types MUST BE IN SAME ORDER AS ORIGINAL SIFT BITS
    #Values are list of residue names (Ex: Met124), and avg. bits for each sift type for corresponding residue
    sift_data = {}
    bit_length = len(sift_types)
    for index in xrange(0, len(avg_sift),bit_length):
        res_name = residue_labels[index].split('_')[0]
        #If residue has non-zero contacts; make a dataset entry
        if numpy.sum(avg_sift[index:index+bit_length]):# and res_name[0:3].upper() != 'GLY':
            if not sift_data.has_key('Residue_Name'):
                sift_data['Residue_Name'] = [res_name]
            else:
                sift_data['Residue_Name'].append(res_name)
            for count, elem in enumerate(sift_types):
                if not sift_data.has_key(elem):
                    sift_data[elem] = [avg_sift[index+count]]
                else:
                    sift_data[elem].append(avg_sift[index+count])
    df = pandas.DataFrame(sift_data, columns = ['Residue_Name']+sift_types)
        
    # Create the general plot and the "subplots" i.e. the bars
    fontS=18
    f, ax1 = plt.subplots(1, figsize=(18,8))
    f.subplots_adjust(bottom=0.2)
    # Set the bar width
    bar_width = 0.75
    # positions of the left bar-boundaries
    bar_l = [i+1 for i in range(len(df['Residue_Name']))]
    # positions of the x-axis ticks (center of the bars as bar labels)
    tick_pos = [i+(bar_width/2) for i in bar_l]
    bottom = []
    for sift_cat, use_col in zip(sift_types, color_list):
        # Create a bar plot, in position bar_1
        if len(bottom) == 0:
            bottom = df[sift_cat]
            ax1.bar(bar_l, df[sift_cat], width=bar_width, label=sift_cat, alpha=0.5, color=use_col, edgecolor=use_col)
        else:
            ax1.bar(bar_l, df[sift_cat], width=bar_width, label=sift_cat, alpha=0.5, color=use_col, bottom = bottom, edgecolor=use_col)
            bottom = bottom + df[sift_cat]
    # set the x ticks with names
    plt.xticks(tick_pos, df['Residue_Name'], rotation="vertical", fontsize=fontS)
    # Set the label and legends
    ax1.set_ylabel("Contact Fraction", fontsize=fontS)
    plt.yticks([0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0], fontsize=fontS)
    ax1.set_xlabel("Residues", fontsize=fontS)
    plt.legend(loc='upper left')
    # Set a buffer around the edge
    plt.xlim([min(tick_pos)-(bar_width/2), max(tick_pos)+(bar_width/2)])
    # Save fig
    plt.savefig(filename, dpi=300)
    plt.close('all')

def plot_CMDS(distance_matrix,labels, marker_size_clust=5, marker_size_noise=2, marker_clust='o', marker_noise='.',save_path=None,filename='CMDS_Embed.png'):
    '''
    Calculates and plots a classical mds embedding of the clusters
    based on the distance matrix (can be either similarity/dissimilarity matrix) of the trajectory frames
    '''
    if save_path == None:
        save_path = os.path.abspath(os.path.dirname(sys.argv[0]))
    #Do Classical MDS
    pos, eigen_val=trajutils.cmdscale(distance_matrix)
    fig = plt.figure()
    #Plot noise first
    unique_labels=sorted(set(labels))
    colors = plt.cm.Spectral(numpy.linspace(0, 1, len(unique_labels)))
    for k, col in zip(unique_labels, colors):
        msize=marker_size_clust
        marker=marker_clust
        if k == -1:
            col = 'lightgray'
            msize=marker_size_noise
            marker=marker_noise
        class_member_mask = (labels == k)
        xy=pos[class_member_mask,:]
        plt.plot(xy[:,0],xy[:,1],marker,markersize=msize,markerfacecolor=col, label='Cluster_%d'%k)
    plt.legend(loc='best', fontsize=10)
    plt.savefig(os.path.join(save_path,filename))
    plt.close('all')
    
def plot_2D(mat_2D=None, outfile='Plot_2D.png', f_size=(4.0,4.0), vmin=None, vmax=None, xlab='', ylab='', cbar_lab='',
            font_lab= 3.0, font_tick=1.5, xtick_lab=[], ytick_lab=[], cmap='jet', interpolation='nearest'):
    if vmin == None:
        vmin = numpy.amin(mat_2D)
    if vmax == None:
        vmax = numpy.amax(mat_2D)
        
    f1=plt.figure(figsize=f_size)
    #f1.subplots_adjust(left=0.18,bottom=0.15)
    ax = plt.gca()
    #cmap= colors.ListedColormap(['#FFFFFF', '#FF99FF','#FF00FF','#00FFFF','#FFD700','#FF0000','#0000CD','#006400'])
    #bounds=[0.0,0.2,0.4,0.6,0.8,1.0]
    #bounds=[0.0,0.2,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
    #norm = colors.BoundaryNorm(bounds, cmap.N)

    img= plt.imshow(mat_2D,interpolation=interpolation,origin='lower', vmin=vmin,vmax=vmax, cmap=cmap)
    plt.xlabel(xlab, fontsize = font_lab)
    plt.ylabel(ylab, fontsize = font_lab)
    #ax.set_xticks(range(len(mat_2D)))
    #ax.set_yticks(range(len(mat_2D)))
    if xtick_lab:
        ax.set_xticklabels(xtick_lab, fontsize=font_tick)
    if ytick_lab:
        ax.set_yticklabels(ytick_lab, fontsize=font_tick)

    cbar=plt.colorbar(img)
    cbar.set_label(cbar_lab)
    #for t in cbar.ax.get_yticklabels():
        #t.set_fontsize(5)
    #plot lines defining different target boundaries 
    
    plt.savefig(outfile, dpi=300)
    plt.close('all')
