This script was used to evaluate the degree to which select image restoration methods improve particle tracking outcomes. The associated manuscript can be found [here](https://doi.org/10.1091/mbc.E20-11-0689). The datasets used can be accessed [here](https://osf.io/59s4k/?view_only=bc7991629cc24ccc93856a689dd7e3c4)

## Installation

### Python 3.6 environment

The following Python packages are required to run this script:
    
    numpy==1.17.2
    tensorflow-gpu==1.12.0
    Keras==2.2.4
    csbdeep==0.4.0
    n2v==0.2.1
    bm3d==3.0.7
    tifffile==2020.9.3
    
### MATLAB

While the denoising methods we have used are available for Python, the tracking methods are written in MATLAB. To run them automatically, we require [MATLAB R2019b](https://www.mathworks.com/products/new_products/release2019b.html), the [MATLAB Image Processing Toolbox](https://www.mathworks.com/products/image.html) and the [MATLAB Engine API for Python](https://www.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html).
    
### ImageJ

We use an ImageJ plugin called StackReg for image stack registration. Please install [ImageJ](https://imagej.nih.gov/ij/), the [StackReg](http://bigwww.epfl.ch/thevenaz/stackreg/) plugin, and [MIJ](http://bigwww.epfl.ch/sage/soft/mij/) to run ImageJ within MATLAB.

### ND-SAFIR

[ND-SAFIR](https://team.inria.fr/serpico/software/nd-safir/) is available as a binary file. Please download it and indicate the binary's location in `env.json` as outlined below.

## Configuration

Please modify `env.json` so that:

 * `imageJ` contains the fully specified path to the `scripts` folder of your ImageJ installation, which should countain `Miji.m`
 * `ndsafir` contains the fully specified path of the ND-SAFIR binary you have downloaded above.
 
 If set up correctly, your `env.json` may look like this:
 
```
{
  "imageJ" : "/home/user/Documents/Fiji.app/scripts",
  "ndsafir" : "/home/user/Documents/ndsafir/bin/ndsafir"
}
```

## Running

When run, `denoise.py` will apply the denoising methods and calculate statistics to determine performance. 

## Dataset organization

Different replicates and imaging conditions are organized into `experiments`, which is reflected in the dataset directory structure:

```
experiments
├── fixed
│   ├── 10ms
│   └── 3ms
└── live
    ├── 10ms
    └── 3ms
```

Each `experiment` follows the below structure:

```
.
├── data
│   ├── cell_0
│   ├── ...
│   └── cell_n
├── models
├── out
├── track
└── train
    ├── s
    └── t
```

* `data` contains the test set images for each cell imaged.
* trained model weights are stored in `models` to be re-used
* images resulting from denoising are written to `out` for visual inspection
* raw tracking output is written to `track`
* `train` contains `s`ource and `t`arget image files for training supervised methods such as `CARE`
* Final performance statistics are written as `.csv` files at the root of each `experiment`

Finally, experiment data are organized as follows:

```
.
├── dn
├── gt
├── ns
└── ns_uncut
```

* `dn` is used to temporarily store denoised frames
* `gt` contains the ground-truth images
* `ns` contains the corresponding noisy input images
* `ns_uncut` contains the same noisy input images without cropping
