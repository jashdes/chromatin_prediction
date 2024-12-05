from __future__ import print_function, unicode_literals, absolute_import, division
import numpy as np
from csbdeep.utils import download_and_extract_zip_file, axes_dict, plot_some, plot_history
from csbdeep.utils.tf import limit_gpu_memory
from csbdeep.io import load_training_data
from csbdeep.models import Config, CARE
from csbdeep.data import RawData, create_patches
from csbdeep.utils import Path, download_and_extract_zip_file, plot_some
from csbdeep.io import load_training_data, save_tiff_imagej_compatible
from csbdeep.models import CARE
from tifffile import imread, imsave
from os import listdir, makedirs, path, remove
from os.path import isfile, join, dirname, isdir, basename, abspath
from glob import glob
import re
import shutil
import matlab.engine
import json
import random
from n2v.models import N2VConfig, N2V
from n2v.utils.n2v_utils import manipulate_val_data
from n2v.internals.N2V_DataGenerator import N2V_DataGenerator
import tensorflow as tf
from bm3d import bm3d
import subprocess

from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense

# Warmup; to ensure keras is functional
x_input = np.array([[1,2,3,4,5]])
y_input = np.array([[10]])
model = Sequential()
model.add(Dense(units=32, activation="tanh", input_dim=x_input.shape[1], kernel_initializer='random_normal'))
model.add(Dense(units=1, kernel_initializer='random_normal'))
model.compile(loss='mse', optimizer='sgd', metrics=['accuracy'])
history = model.fit(x_input, y_input, epochs=10, batch_size=32)


eng = matlab.engine.start_matlab()

trackdir = abspath("MATLAB/tracker")
helpers = abspath("MATLAB/helpers")

imageJ = None 
ndsafir = None
with open('env.json') as json_file: 
	environment = json.load(json_file)
	imageJ = environment["imageJ"]
	ndsafir = environment["ndsafir"]


def predictBM3D(images, outdir):
	for image in images:
		x = imread(image)
		xhat = bm3d(x, np.max(x) * 0.02)
		save_tiff_imagej_compatible(join(outdir, basename(image)), xhat[...,None], axes='YXC')

def predictNDSAFIR(images, outdir, time=False):
	if time:
		image = images[0]
		args = [ndsafir, "-i", image,"-2dt", "true", "-o", join(outdir,basename(image) + ".ndsafir.tif")]
		subprocess.run(args, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
	else:
		for image in images:
			args = [ndsafir, "-i", image, "-o", join(outdir,basename(image))]
			subprocess.run(args, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

def trainN2V(trainimgs, evalstring="n2v",psize=64, prefix = None, structmask = None):
	if prefix is not None:
		basedir = join(prefix,"models")
	else:
		basedir = "models"
	if not isfile(join(basedir,evalstring,"weights_last.h5")):
		datagen = N2V_DataGenerator()
		X = None
		X_val = None
		for trainimg in trainimgs:
			print("Processing " + trainimg)
			img = datagen.load_imgs([trainimg])[0]
			imgs_train = img #img[:,:int(max(psize,img.shape[1]*.9))]
			#imgs_vali = img[:,int(max(psize,img.shape[1]*.9)):]
			patch_shape = (psize,psize)
			cX = datagen.generate_patches(imgs_train, shape=patch_shape)
			#cX_val = datagen.generate_patches(imgs_vali, shape=patch_shape)
			
			if X is None:
				X = cX
			#	X_val = cX_val	
			else:
				X = np.append(X,cX, axis=0)
			#	X_val = np.append(X_val,cX_val, axis=0)
		if structmask is not None:
			perc_pix = 100/(psize*psize) 
		else:
			perc_pix = 0.198

		config = N2VConfig(X, unet_kern_size=3, 
                   train_steps_per_epoch=int(int(X.shape[0]*0.9)/128), train_epochs=100, train_loss='mse', batch_norm=True, 
                   train_batch_size=128, n2v_perc_pix=perc_pix, n2v_patch_shape=(psize, psize), 
                   n2v_manipulator='uniform_withCP', n2v_neighborhood_radius=5, structN2Vmask = structmask)
		
		vars(config)
		model=N2V(config,evalstring,basedir=basedir)
		np.random.shuffle(X)
		history = model.train(X[:int(X.shape[0]*.9)],X[int(X.shape[0]*.9):])

def predictN2V(images, outdir, evalstring="n2v", prefix=None):
	if prefix is not None:
		basedir = join(prefix,"models")
	else:
		basedir = "models"
	print("Loading model from " + join(basedir, evalstring))
	model = N2V(config=None, name=evalstring, basedir=basedir)
	model.config.axes='YXC'		
	model.load_weights('weights_last.h5')

	for image in images:
		x = imread(image)[...,None]
		xhat = model.predict(x.astype("float32"), axes='YXC')
		save_tiff_imagej_compatible(join(outdir, basename(image)), xhat, axes='YXC')

def trainCARE(trainsource, evalstring="care",psize=64, npatches=200, trainN = None, prefix = None):
	if prefix is not None:
		basedir = join(prefix,"models")
	else:
		basedir = "models"
	if not isfile(join(basedir,evalstring,"weights_best.h5")):
		if not isdir(join(trainsource,'s')):
			makedirs(join(trainsource,'s'), exist_ok=True)
			makedirs(join(trainsource,'t'), exist_ok=True)
			cells = [d for d in listdir(trainsource) if not isfile(join(trainsource, d))]
			for cell in cells:
				noisy = glob(join(trainsource,cell,"ns","*.tif"))
				for f in noisy:
					shutil.copy(f, join(trainsource,"s",basename(f)) )
				gt = glob(join(trainsource,cell,"gt","*.tif"))
				for f in gt:
					if not basename(f) == "cropframe.tif":
						shutil.copy(f, join(trainsource,"t",basename(f)) )
		images = [f for f in listdir(join(trainsource,'s')) if isfile(join(trainsource,'s', f))]
		random.shuffle(images)

		for d in ["s","t"]:
			shutil.rmtree(join("data",d))	
			makedirs(join("data",d), exist_ok=True)

		for i,image in enumerate(images):
			if trainN is not None and (i == trainN):
				break 
			m = re.search('(.*)_ns-(\d*).tif',image)
			bname = m.group(1) 
			frame = m.group(2)
			destname = bname + frame + ".tif"
			filename = bname + "_ns-" + frame + ".tif"
			shutil.copy(join(trainsource,'s', filename), join("data","s",destname))
			filename = bname + "_gt-" + frame + ".tif"
			shutil.copy(join(trainsource,'t', filename), join("data","t",destname))

		raw_data = RawData.from_folder(basepath='data', source_dirs=['s'], target_dir='t', axes='YX')
		X, Y, axes = create_patches(raw_data, patch_size=(psize,psize), n_patches_per_image=npatches, save_file="data/train.npz")
		(X,Y), (X_val,Y_val), axes = load_training_data('data/train.npz', validation_split=0.1, verbose=True)
		c = axes_dict(axes)['C']
		n_channel_in, n_channel_out = X.shape[c], Y.shape[c]
		config = Config(axes, n_channel_in, n_channel_out, probabilistic=False, train_steps_per_epoch=max(1,int(X.shape[0]/npatches)))
		vars(config)
		model = CARE(config, evalstring, basedir=basedir)
		history = model.train(X,Y, validation_data=(X_val,Y_val))
		#model.export_TF()

def predictCARE(images, outdir, evalstring="n2v", prefix=None):
	if prefix is not None:
		basedir = join(prefix,"models")
	else:
		basedir = "models"
	print("Loading model from " + join(basedir, evalstring))
	model = CARE(None,  evalstring, basedir=basedir)
	axes = 'YX'
	for image in images:
		x = imread(image)
		xhat = model.predict(x, axes).astype(x.dtype)
		save_tiff_imagej_compatible(join(outdir, basename(image)), xhat, axes)

def track(datadir, outputdir):
	# Convert denoised images to trackable stacks
	eng.cd(helpers)
	eng.addpath("functions",nargout=0)
	eng.addpath("workflows",nargout=0)
	eng.postdenoise(datadir,"templates/7x7template.tif", nargout=0)

	# Move trackable stacks to tracker and apply StackReg
	# Note that this requires ij.jar and mij.jar in the MATLAB java path
	stacks = [f for f in listdir(datadir) if isfile(join(datadir, f))]

	eng.addpath(imageJ ,nargout=0)

	
	# clean up previous data
	mout = join(trackdir, "mout", "BATCH_outputSS")
	files = [f for f in listdir(mout) if isfile(join(mout, f))]
	for f in files:
		remove(join(mout, f))

	binput = join(trackdir, "BATCH_inputData")
	files = [f for f in listdir(binput) if isfile(join(binput, f))]
	for f in files:
		remove(join(binput, f))
	
	# register new images
	for stack in stacks:
		eng.stackRegPath(join(datadir,stack),join(trackdir, "BATCH_inputData", stack), nargout=0)

	eng.addpath(trackdir,nargout=0)
	# Track stacks
	eng.Main_Batch_7x7f(binput,join(trackdir, "mout"), "templates/7x7template.tif", 9.26,0.31, 0, nargout=0)
	
	# Convert tracks to .txt format
	track = join(outputdir,"txt")
	makedirs(track, exist_ok=True)
	files = [f for f in listdir(track) if isfile(join(track, f))]
	for f in files:
		remove(join(track, f))
	mtrack = join(outputdir,"mat")
	makedirs(mtrack, exist_ok=True)
	files = [f for f in listdir(mtrack) if isfile(join(mtrack, f))]
	for f in files:
		remove(join(mtrack, f))
	files = [f for f in listdir(mout) if isfile(join(mout, f))]
	for f in files:
		shutil.copy(join(mout,f), join(mtrack,f))
		
		
	
	eng.mat2text(mout,track, nargout=0)
	
	# Calculate metrics
	total = eng.stackRestorationMetrics(track, nargout=1)
	return total

def compareMethods(experimentdir, methods):
	experimentdir = abspath(experimentdir)
	# create output directories
	for d in ["data", "models", "out", "track"]:
		makedirs(join(experimentdir,d), exist_ok=True)
	resfile = join(experimentdir,"results.json")
	if isfile(resfile):
		with open(resfile) as json_file: 
		    results = json.load(json_file)
	else:
		results = {}
	validationdir = join(experimentdir,"data")
	cells = [d for d in listdir(validationdir) if not isfile(join(validationdir, d))]

	for traintype in methods:
		# Return previous result if it exists
		if traintype in results:
			continue

		# Remove previous output
		oldfiles = [d for d in listdir(validationdir) if isfile(join(validationdir, d))]
		for f in oldfiles:
			remove(join(validationdir,f))
		for cell in cells:
			olddenoised = glob(join(validationdir,cell,"dn","*.tif"))
			for f in olddenoised:
				remove(f)
		
		for cell in cells:
			noisy = glob(join(validationdir,cell,"ns","*.tif"))
			if traintype == "gt":
				modelstring = "gt"
				gt = glob(join(validationdir,cell,"gt","*.tif"))
				for img in gt:
					shutil.copy(img, join(validationdir,cell,"dn",basename(img)) )
				continue

			if traintype == "bm3d":
				modelstring = "bm3d"
				ns = glob(join(validationdir,cell,"ns","*.tif"))
				predictBM3D(ns,join(validationdir,cell,"dn"))
				continue
						
			if traintype == "ndsafir_t":
				modelstring = "ndsafir_t"
				ns = glob(join(validationdir,cell,"*_ns.tif"))
				predictNDSAFIR(ns,join(validationdir,cell), True)
				# Convert denoised images to trackable stacks
				eng.cd(helpers)
				eng.addpath("functions",nargout=0)
				eng.addpath("workflows",nargout=0)
				dnspath = join(validationdir,cell,cell+"_ns.tif.ndsafir.tif")
				eng.explodeNDSAFIR(dnspath, join(validationdir,cell,"dn"), cell+"_ns",nargout=0)

				continue

			if traintype == "ndsafir":
				modelstring = "ndsafir"
				ns = glob(join(validationdir,cell,"ns","*.tif"))
				predictNDSAFIR(ns,join(validationdir,cell,"dn"))
				continue

			if traintype == "ns":
				modelstring = "ns"
				ns = glob(join(validationdir,cell,"ns","*.tif"))
				for img in ns:
					shutil.copy(img, join(validationdir,cell,"dn",basename(img)) )
				continue

			tsplit = traintype.split('-')
			if tsplit[0] == "care":
				modelstring = traintype
				if len(tsplit) > 1:
					ntrain = int(tsplit[1])
				else:
					ntrain = None 
				trainCARE(join(experimentdir,"train"), modelstring, 28, 64, ntrain, experimentdir)
				predictCARE(noisy, join(validationdir,cell,"dn") , modelstring, experimentdir)
				continue


			if traintype == "n2v_each":
				noisy_train = glob(join(validationdir,cell,"ns","*.tif"))
				for trainimg in noisy_train:
					targetname = basename(trainimg)
					modelstring = cell + "_" + traintype + "_" + targetname 
						
					# Train network with specified parameter
					trainN2V([trainimg],modelstring,64, experimentdir)

					# Use network to denoise images
					noisy = glob(join(validationdir,cell,"ns",targetname))
					predictN2V(noisy, join(validationdir,cell,"dn"), modelstring, experimentdir)
			else:
				modelstring = cell + "_" + traintype
				if traintype == "n2v_first":
					trainimg = glob(join(validationdir,cell,"ns","*-1.tif"))
				if traintype == "n2v_full":
					trainimg = glob(join(validationdir,cell,"ns_uncut","*.tif"))
					trainimg = trainimg[:40]
				elif traintype == "n2v" or "structn2v" in traintype:
					trainimg = glob(join(validationdir,cell,"ns","*.tif"))
					
				# Train network with specified parameter
				if "structn2v" in traintype:
					msplit = traintype.split('-')
					masky = int(msplit[1])
					maskx = int(msplit[2])
					structmask = np.zeros((maskx*2+1,masky*2+1), dtype=int)
					structmask[maskx,:] = 1
					structmask[:,masky] = 1
					trainN2V(trainimg,modelstring, 64, experimentdir, structmask.tolist())
				else:
					trainN2V(trainimg,modelstring, 64, experimentdir)

				# Use network to denoise images
				noisy = glob(join(validationdir,cell,"ns","*.tif"))
				predictN2V(noisy, join(validationdir,cell,"dn"), modelstring, experimentdir)

		# Track denoised images
		results[traintype] = track(join(validationdir),join(experimentdir,"track", traintype))

		# Save output 
		dst = join(experimentdir,"out",traintype)
		makedirs(dst, exist_ok=True)
		out = glob(join(trackdir,"BATCH_inputData","*_dn.tif"))
		for img in out:
			shutil.copy(img, dst)

		with open(resfile, 'w') as fp:
			json.dump(results, fp)
	return results

exps = [join("live","3ms"),join("live","10ms"),join("fixed","3ms"),join("fixed","10ms")]

# Denoise and track

for ex in exps:
	masks = []
	for xmask in [0,3,8,16]:
		for ymask in [0,3,8,16]:
			if xmask == 0 and ymask == 0:
				continue
			masks.append("structn2v-" + str(xmask) + "-" + str(ymask))
	compareMethods(join("experiments",ex),["gt","ndsafir","bm3d","ns", "n2v", *masks,"care"])

# Calculate statistics

eng.cd(helpers)
eng.addpath("functions", nargout=0)
eng.addpath("workflows", nargout=0)

for ex in exps:
	if "fixed" in ex:
		eng.exportDenoise(abspath(join("experiments",ex)),abspath(join("experiments",ex)), True, nargout=0)
	else:
		eng.exportDenoise(abspath(join("experiments",ex)), nargout=0)

	
eng.quit()

