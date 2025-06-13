import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os, sys

import starfish
from starfish import Experiment
from starfish import display
from starfish import data, FieldOfView
from starfish.types import Features, Axes

from starfish import IntensityTable

from starfish.image import Filter
from starfish.spots import DetectPixels

from datetime import datetime
from multiprocessing import Pool 	# Kian: added 210602
import functools	# Kian: added 210602


def DARTFISH_pipeline(fov, codebook, magnitude_threshold):
	imgs = fov.get_image(starfish.FieldOfView.PRIMARY_IMAGES)

	gauss_filt = Filter.Laplace(2, True)
	gauss_imgs = gauss_filt.run(imgs)
	
	sc_filt = Filter.Clip(p_max=100, expand_dynamic_range=True)
	norm_imgs = sc_filt.run(gauss_imgs)

	# z_filt = Filter.ZeroByChannelMagnitude(thresh=.25, normalize=True)
	# filtered_imgs = z_filt.run(norm_imgs)
	
#	filtered_imgs = filtered_imgs.apply(lambda x: 255 / 50 * np.clip(x, 0, 50/255)) # 211019


	def compute_magnitudes(stack, norm_order=2):
		pixel_intensities = IntensityTable.from_image_stack(stack)
		feature_traces = pixel_intensities.stack(traces=(Axes.CH.value, Axes.ROUND.value))
		norm = np.linalg.norm(feature_traces.values, ord=norm_order, axis=1)
		return norm

	mags = compute_magnitudes(norm_imgs)
	
	# how much magnitude should a barcode have for it to be considered by decoding? this was set by looking at
	# the plot above
#	magnitude_threshold = 0.2
	# how big do we expect our spots to me, min/max size. this was set to be equivalent to the parameters
	# determined by the Zhang lab.
	area_threshold = (1,100)
	# how close, in euclidean space, should the pixel barcode be to the nearest barcode it was called to?
	# here, I set this to be a large number, so I can inspect the distribution of decoded distances below
	# distance_threshold = 3

	# psd = DetectPixels.PixelSpotDecoder(
		# codebook=codebook,
		# metric='euclidean',
		# distance_threshold=distance_threshold,
		# magnitude_threshold=magnitude_threshold,
		# min_area=area_threshold[0],
		# max_area=area_threshold[1]
	# )

	# initial_spot_intensities, results = psd.run(filtered_imgs)
	
	# spots_df = initial_spot_intensities.to_features_dataframe()
	# spots_df['area'] = np.pi*spots_df['radius']**2
	# spots_df = spots_df.loc[spots_df[Features.PASSES_THRESHOLDS]]
	
	distance_threshold = 1
	
	psd = DetectPixels.PixelSpotDecoder(
		codebook=codebook,
		metric='euclidean',
		distance_threshold=distance_threshold,
		magnitude_threshold=magnitude_threshold,
		min_area=area_threshold[0],
		max_area=area_threshold[1]
	)

	spot_intensities, results = psd.run(norm_imgs)
	spot_intensities = IntensityTable(spot_intensities.where(spot_intensities[Features.PASSES_THRESHOLDS], drop=True))
	# reshape the spot intensity table into a RxC barcode vector
	pixel_traces = spot_intensities.stack(traces=(Axes.ROUND.value, Axes.CH.value))

	# extract dataframe from spot intensity table for indexing purposes
	pixel_traces_df = pixel_traces.to_features_dataframe()
	pixel_traces_df['area'] = np.pi*pixel_traces_df.radius**2
	return pixel_traces_df, mags


def process_experiment(experiment: starfish.Experiment, output_dir, magnitude_threshold):
	decoded_intensities = {}
	regions = {}
	count = 0
	df_pipe_partial = functools.partial(DARTFISH_pipeline, codebook=experiment.codebook, magnitude_threshold = magnitude_threshold)
	# fovs = [fov for (name_, fov) in experiment.items()]
	# names = [name_ for (name_, fov) in experiment.items()]
	names, fovs = zip(*experiment.items())
	with Pool(8) as p:
		traces_dfs, mags = zip(*p.map(df_pipe_partial, list(fovs)))
	print(traces_dfs)
	# for i, (name_, fov) in enumerate(experiment.items()):
	# 	print(datetime.now().strftime('%Y-%d-%m_%H:%M:%S: Started Processing FOV {:02d} with Barcode Magnitude threshold {}'.format(count,magnitude_threshold)))
	# 	df_pipe_partial = functools.partial(DARTFISH_pipeline, codebook=expriment.codebook, magnitude_threshold = magnitude_threshold, 
	# 		normalize=normalize)
	# 	pixel_traces_df, mags = DARTFISH_pipeline(fov, experiment.codebook, magnitude_threshold, normalize)
	for name, trace in zip(names, traces_dfs):
		newName = "FOV" + name[4:]
		trace.to_csv(os.path.join(output_dir,'starfish_table_bcmag_{0}_{1}'.format(magnitude_threshold,newName) + '.csv'))
	# 	print(datetime.now().strftime('%Y-%d-%m_%H:%M:%S: Finished Processing FOV {:02d} with Barcode Magnitude threshold {}'.format(count,magnitude_threshold)))
	# 	count += 1
		#decoded_intensities[name_] = decoded
		#regions[name_] = segmentation_results
#	return decoded_intensities, regions


data_dir = "../3_Decoded/data_Starfish"
output_dir = "../3_Decoded/output_Starfish_normLGs2_1100_dist1"
if not os.path.exists(output_dir):
	os.makedirs(output_dir)

currentTime = datetime.now() 
reportFile = os.path.join(output_dir, currentTime.strftime("%Y-%d-%m_%H:%M_starfish.log"))
sys.stdout = open(reportFile, 'x') # redirecting the stdout to the log file


exp = Experiment.from_json(os.path.join(data_dir,"experiment.json"))

magnitude_thresholds = [0.23] #,1.7
for magnitude_threshold in magnitude_thresholds:
	output_path = os.path.join(output_dir,"bcmag{}".format(magnitude_threshold))
	if not os.path.exists(output_path):
		os.makedirs(output_path)
	print(datetime.now().strftime('%Y-%d-%m_%H:%M:%S: Started Processing Experiment with Barcode Magnitude threshold ' + str(magnitude_threshold)))
	sys.stdout.flush()
	process_experiment(exp, output_path, magnitude_threshold)
	print(datetime.now().strftime('%Y-%d-%m_%H:%M:%S: Finished Processing Experiment with Barcode Magnitude threshold ' + str(magnitude_threshold)))
	sys.stdout.flush()
'''
magnitude_thresholds = [1.0]
for magnitude_threshold in magnitude_thresholds:
	output_path = os.path.join(output_dir,"bcmag{}".format(magnitude_threshold))
	if not os.path.exists(output_path):
		os.makedirs(output_path)
	print(datetime.now().strftime('%Y-%d-%m_%H:%M:%S: Started Processing Experiment with Barcode Magnitude threshold ' + str(magnitude_threshold)))
	sys.stdout.flush()
	process_experiment(exp, output_path, magnitude_threshold,normalize=True)
	print(datetime.now().strftime('%Y-%d-%m_%H:%M:%S: Finished Processing Experiment with Barcode Magnitude threshold ' + str(magnitude_threshold)))
	sys.stdout.flush()
'''
# magnitude_thresholds = [0.3]
# for magnitude_threshold in magnitude_thresholds:
# 	output_path = os.path.join(output_dir,"bcmag{}".format(magnitude_threshold))
# 	if not os.path.exists(output_path):
# 		os.makedirs(output_path)
# 	print(datetime.now().strftime('%Y-%d-%m_%H:%M:%S: Started Processing Experiment with Barcode Magnitude threshold ' + str(magnitude_threshold)))
# 	sys.stdout.flush()
# 	process_experiment(exp, output_path, magnitude_threshold,normalize=False)
# 	print(datetime.now().strftime('%Y-%d-%m_%H:%M:%S: Finished Processing Experiment with Barcode Magnitude threshold ' + str(magnitude_threshold)))
# 	sys.stdout.flush()

# magnitude_thresholds = [0.9]#[0.1,0.2,0.4,0.5,0.6,0.7,0.8,0.9]
# for magnitude_threshold in magnitude_thresholds:
# 	output_path = os.path.join(output_dir,"bcmag{}".format(magnitude_threshold))
# 	if not os.path.exists(output_path):
# 		os.makedirs(output_path)
# 	print(datetime.now().strftime('%Y-%d-%m_%H:%M:%S: Started Processing Experiment with Barcode Magnitude threshold ' + str(magnitude_threshold)))
# 	sys.stdout.flush()
# 	process_experiment(exp, output_path, magnitude_threshold,normalize=False)
# 	print(datetime.now().strftime('%Y-%d-%m_%H:%M:%S: Finished Processing Experiment with Barcode Magnitude threshold ' + str(magnitude_threshold)))
# 	sys.stdout.flush()
'''
magnitude_thresholds = [1.4, 1.7, 2.0]
for magnitude_threshold in magnitude_thresholds:
	output_path = os.path.join(output_dir,"bcmag{}".format(magnitude_threshold))
	if not os.path.exists(output_path):
		os.makedirs(output_path)
	print(datetime.now().strftime('%Y-%d-%m_%H:%M:%S: Started Processing Experiment with Barcode Magnitude threshold ' + str(magnitude_threshold)))
	sys.stdout.flush()
	process_experiment(exp, output_path, magnitude_threshold,normalize=True)
	print(datetime.now().strftime('%Y-%d-%m_%H:%M:%S: Finished Processing Experiment with Barcode Magnitude threshold ' + str(magnitude_threshold)))
	sys.stdout.flush()	
'''
sys.stdout = sys.__stdout__ # restoring the stdout pipe to normal	
