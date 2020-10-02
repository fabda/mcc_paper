# Code Repository

## Paper with no name yet

### Author list

Preprint : no preprint yet

doi: no doi yet

## 1 - Setup

First, you need to setup a Python Anaconda environment (https://www.anaconda.com/products/individual) or miniconda (Python 3.5 or higher), and clone this repository on your computer.

Then create a virtual environment and process to the python requirements install by typing into the cloned repository on your computer

```shell
pip install -r requirements.txt
```

## 2 - Video Refocusing and Stabilization

In order to overcome the loss of focus events, artefactual movements and tissue deformation which made the detailed analysis of cell migration and their path tracking impossible, a semi automatic video refocus and stabilization procedure was used.

```shell
$ python refocus.py -h
usage: refocus.py [-h] [--stack_name stack_name]
                  [--tracking_file tracking_file]

Refocus an ImageJ TIFF hyperstack

optional arguments:
  -h, --help            show this help message and exit
  --stack_name stack_name
                        path to the TIFF hyperstack
  --tracking_file tracking_file
                        path to the CSV file containing the tracking
```

shell command to start refocus.py

```shell
$ python refocus.py --stack_name=path_to_the_stack.tif --tracking_file=path_to_the_tracking_file.csv
```

Plugin Fiji "Image Stabilizer": https://imagej.net/Image_Stabilizer


Supplemental Video: https://amubox.univ-amu.fr/s/WxWYmd3f62xw9Bd



## 3 - Tracking and Analysis

For each refocused and stabilized video, MCCs centroids coordinates have been manually tracked using the "Manual Tracking" Fiji plugin and their position corrected using an one-dimensional  Gaussian filter (sigma=10). A set of "behavioral descriptors" are then calculated to help quantify individual cell movements and interactions throughout a whole movie using an homemade Python script

```shell
$ python analysis.py -h

usage: analysis.py [-h] [--tracking_file tracking_file] [--sigma sigma]
                   [--window window] [--pixmicron pixmicron]
                   [--framesec framesec] [--result_file result_file]

Generate a set of Behavioral Descriptors from a Fiji Tracking File

optional arguments:
  -h, --help            show this help message and exit
  --tracking_file tracking_file
                        path to the tracking CSV file
  --sigma sigma         sigma value for 1-d gaussian filter
  --window window       window size for direction calculation
  --pixmicron pixmicron
                        1 pixel equals ? µm
  --framesec framesec   1 frame equals ? s
  --result_file result_file
                        CSV result filename
```
An example of shell command to start analysis.py using "trackingfile.csv" as input

```shell
$ python analysis.py --sigma=10 --window=5 --pixmicron=0.4 --framesec=60 --tracking_file="../trackfile.csv" --result_file="result.xls"

```


## 4 - License

This code repository is under [GPL v3.0 License](https://github.com/fabda/mcc_paper/blob/master/LICENSE)
