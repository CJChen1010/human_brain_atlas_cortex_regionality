import os
from os import path
from code_lib.Segmentation_220209 import Segmentor2D
from code_lib.Assignment_201020 import *
from skimage.io import imread
import numpy as np, pandas as pd
from skimage.segmentation import expand_labels

def mask2centroid(maskImg):
    centroids = []
    for i in range(1, maskImg.max() + 1):
        xs, ys = np.where(maskImg == i)
        xc, yc = xs.mean().astype(int), ys.mean().astype(int)
        centroids.append((xc, yc))
    return np.array(centroids)

def mask2centroid_parallel(rng, mimg):
    cent = []
    for i in rng:
        xs, ys = np.where(mimg == i)
        xc, yc = xs.mean().astype(int), ys.mean().astype(int)
        cent.append((xc, yc))    
    return np.array(cent)

def plotRolonies2d(rolonyDf, nucLabels, coords = ['y', 'x'], label_name = 'nucleus_label', backgroudImg = None, ax = None, backgroundAlpha = 0.5):
    # myCmap = np.random.rand(np.max(nucLabels) + 1, 4)
    # myCmap[:, -1] = 1
    # myCmap[0] = (0, 0, 0, 1)
    # myCmap = ListedColormap(myCmap)
    # plt.style.use('dark_background')
    if ax is None:
        fig, ax = plt.subplots(nrows = 1, figsize = (18, 11))
    
    boundaries = nucLabels.copy()
    print("finding boundaries")
    boundaries[~find_boundaries(nucLabels)] = 0


    if not backgroudImg is None:
        ax.imshow(backgroudImg, alpha = backgroundAlpha, cmap = 'gray')
    
    # ax.imshow(mask2rgb(boundaries, myCmap), alpha = 0.7)#, vmin = 0, vmax = myCmap.N)
    # print(boundaries)
    ax.imshow(boundaries, cmap = myCmap, alpha=0.7)

    # ax.imshow(nucLabels, alpha = 0.5, cmap = myCmap, vmin = 0, vmax = myCmap.N)

    for i, rol in rolonyDf.iterrows():
        circ = plt.Circle((rol[coords[0]], rol[coords[1]]), 2 * rol['radius'], 
                          linewidth = 0.7, fill = False, alpha = 0.8, 
                          color = myCmap(rol[label_name]))
        ax.add_patch(circ)
        # if (i %100) == 0:
        #     print(i)
    plt.tight_layout()
    plt.show() 


nuc_path = '../2_Registered_NOR6/stitched/MIP_draq5_ch00.tif'

saving_path = '../4_CellAssignment'

bcmag = 'bcmag0.23'

if not path.exists(saving_path):
    os.makedirs(saving_path)
    
spot_file = '../3_Decoded/output_Starfish_normLGs2_1100_dist1/{}/all_spots_empty6_bc0.7_filtered3.0.tsv'.format(bcmag)

nuc_img = imread(nuc_path)

# segmenting the nuclear image
# segmentor = Segmentor2D()
# mask = segmentor.segment([nuc_img], diameters = 30, 
#                          out_files = [path.join(saving_path, 'segmentation_mask.npy')])[0]
mask = np.load(path.join(saving_path, 'segmentation_mask.npy'))

### Plotting the segmentation mask.
## mask=mask[0]
myCmap = np.random.rand(np.max(mask) + 1, 4)
myCmap[:, -1] = 1
myCmap[0] = (0, 0, 0, 1)
myCmap = ListedColormap(myCmap)

plt.figure(figsize = (int(mask.shape[0]/600), int(mask.shape[1]/600)))
plt.imshow(mask, cmap = myCmap)
plt.savefig(os.path.join(saving_path, 'mask_nucleus.png'), dpi = 500, bbox_inches='tight')

# # Rolony assignment
# spot_df = pd.read_csv(spot_file, index_col=0, sep = '\t')
# assigner = RolonyAssigner(nucleiImg=mask, rolonyDf=spot_df, axes = ['yg', 'xg'])
# labels, dists = assigner.getResults()

# spot_df['nucleus_label'] = labels
# spot_df['dist2nucleus'] = np.round(dists, 2)
# spot_df = spot_df.sort_values('nucleus_label', ignore_index = True)
# spot_df.to_csv(path.join(saving_path, 'spots_assigned.tsv'), sep = '\t', index = False, float_format='%.3f')

# # Making the cell by gene matrix
# # spot_df = pd.read_csv('../4_CellAssignment_1100_bcmag2.0_filtered0.5_ef0.175/spots_assigned.tsv', index_col=0, sep = '\t')
# nuc_gene_df = spot_df[['nucleus_label', 'gene']].groupby(by = ['nucleus_label', 'gene'], as_index = False).size()
# nuc_gene_df = nuc_gene_df.reset_index().pivot(index = 'nucleus_label', columns = 'gene', values = 'size').fillna(0).astype(int)
# nuc_gene_df.to_csv(path.join(saving_path, 'nucleus-by-gene.tsv'), sep = '\t')

# # # finding the nuclei centroids
# centroids = mask2centroid(mask)
# centroid_df = pd.DataFrame({'nucleus_label' : np.arange(1, mask.max() + 1), 
#                             'centroid_x' : centroids[:, 0], 'centroid_y' : centroids[:, 1]})
# centroid_df.to_csv(path.join(saving_path, 'nucleus_locations.tsv'), sep = '\t', index = False)

# # plotting assigned rolonies
# spot_df=pd.read_csv('../4_CellAssignment_troubleshoot/spots_assigned.tsv')
# # plotting assigned rolonies
# print("plotting assigned rolonies")
# fig = plt.figure(figsize = (int(mask.shape[0]/600), int(mask.shape[1]/600)))
# ax = fig.gca()
# plotRolonies2d(spot_df, mask, coords = ['xg', 'yg'], label_name='nucleus_label', ax = ax, backgroudImg=nuc_img, backgroundAlpha=0.6)
# fig.savefig(path.join(saving_path, 'assigned_rolonies.png'),
#             transparent = False, dpi = 500, bbox_inches='tight')
# print("plotting assigned rolonies done")

# # # # plotting the nuclei with their label
# fig = plt.figure(figsize = (int(mask.shape[0]/600), int(mask.shape[1]/600)))
# ax = fig.gca()
# ax.imshow(nuc_img, cmap='gray')
# ax.scatter(centroids[:, 1], centroids[:, 0], s = 1, c='red')
# for i in range(centroids.shape[0]):
#     ax.text(centroids[i, 1], centroids[i, 0], str(i), fontsize = 5, c = 'orange')
# fig.savefig(path.join(saving_path, 'nuclei_map.png'),
#             transparent = True, dpi = 500, bbox_inches='tight')

# expand nuclei mask by 15um------------------------------------------------------------
# nuclei_mask = np.load(path.join(saving_path, 'segmentation_mask.npy'))
expanded = expand_labels(mask, distance=53)

### Plotting the expanded segmentation mask.
# myCmap = np.random.rand(np.max(expanded) + 1, 4)
# myCmap[:, -1] = 1
# myCmap[0] = (0, 0, 0, 1)
# myCmap = ListedColormap(myCmap)

plt.figure(figsize = (int(expanded.shape[0]/600), int(expanded.shape[1]/600)))
plt.imshow(expanded, cmap = myCmap)
plt.savefig(os.path.join(saving_path, 'expanded_nucleus_mask.png'), dpi = 500, bbox_inches='tight')


# Assign rolonies and filter out rolonies with dists>53 pixels
spot_df = pd.read_csv(spot_file, index_col=0, sep = '\t')
assigner = RolonyAssigner(nucleiImg=mask, rolonyDf=spot_df, axes = ['yg', 'xg'])
labels, dists = assigner.getResults()

spot_df=spot_df.loc[dists<=53]
spot_df['nucleus_label'] = labels[dists<=53]
spot_df['dist2nucleus'] = np.round(dists[dists<=53], 2)
spot_df = spot_df.sort_values('nucleus_label', ignore_index = True)
spot_df.to_csv(path.join(saving_path, 'spots_assigned_expanded_filtered.tsv'), sep = '\t', index = False, float_format='%.3f')

# Making the nucleus by gene matrix
# spot_df = pd.read_csv('../4_CellAssignment_1100_bcmag2.0_filtered0.5_ef0.175/spots_assigned.tsv', index_col=0, sep = '\t')
nuc_gene_df = spot_df[['nucleus_label', 'gene']].groupby(by = ['nucleus_label', 'gene'], as_index = False).size()
nuc_gene_df = nuc_gene_df.reset_index().pivot(index = 'nucleus_label', columns = 'gene', values = 'size').fillna(0).astype(int)
nuc_gene_df.to_csv(path.join(saving_path, 'expanded_nucleus-by-gene.tsv'), sep = '\t')

# plotting assigned rolonies
print("plotting assigned rolonies")
fig = plt.figure(figsize = (int(expanded.shape[0]/600), int(expanded.shape[1]/600)))
ax = fig.gca()
plotRolonies2d(spot_df, expanded, coords = ['xg', 'yg'], label_name='nucleus_label', ax = ax, backgroudImg=nuc_img, backgroundAlpha=0.6)
fig.savefig(path.join(saving_path, 'assigned_rolonies_expanded_filtered.png'),
            transparent = False, dpi = 500, bbox_inches='tight')
print("plotting assigned rolonies done")
