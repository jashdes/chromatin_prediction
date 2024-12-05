# depth-regression

## Environment

```
tensorflow-gpu==1.12.0
Keras==2.2.4
tifffile==2020.9.3
matplotlib==3.1.3
numpy==1.17.2
scikit-image==0.16.2
```

Also ensure your environment has a `Tk` with version at least 8.6 to use the matplotlib interactive UI.

## Data annotation

Input data should be formatted as a `.tiff` image stack, with one depth per image in the stack. To make an image usable in training, we need to specify where in 3-space the beads are within the images. To do so, we create an annotation file with the same name as the to-be-annotated image stack and `.json` appended using the following format:

```
{
    "XY": [
        [
            783.6214387055974,
            155.49871425990221
        ],
        [
            550.5631391921149,
            265.49904260195973
        ],
        [
            365.3718290527076,
            410.73668623882196
        ]
    ],
    "Z": {
        "7": "-1.5",
        "8": "-1.25",
        "9": "-1",
        "10": "-0.75",
        "11": "-0.5",
        "12": "-0.25",
        "13": "0",
        "14": "0.25",
        "15": "0.5",
        "16": "0.75",
        "17": "1",
        "18": "1.25",
        "19": "1.5",
    }
}
```

`XY` is an array of X and Y coordinate tuples; since these two axes don't change as we vary the depth, they are the same for all images.

`Z` is a dictionary that maps the image indices to z-depths: here, images 8 through 20 in the stack are used and go from a depth of -1.5 to 1.5; any other images are ignored.

The `data_labeling` Jupyter notebook shows how such an annotation file is created interactively: first, the source image stack is selected; then the user is prompted for the z-depth of each image. To use an image without the DOE to annotate X and Y positions of spots interactively, a depth of `p` is assigned. To ignore an image, `s` is assigned. Once depths and X and Y locations have been fully specified, the annotation file is written.

## Training

Example data can be downloaded [here](https://drive.google.com/file/d/1QeIpFKwxHj0NYl0gO4r8QlQStq2TjMQ4/view?usp=sharing).
The `train+test` Jupyter notebook shows how the network is both trained and validated on these labeled datasets. The training image generator will use any `.tiff` stack that has an annotation file (`.json` appeneded to its name) in the same directory. To augment the data, we introduce random X and Y shifts, but no rotation or scaling as this would impact the DOE pattern itself. Performance is validated via a separate dataset and average error is calculated.
