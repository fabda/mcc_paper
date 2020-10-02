# Code Repository

## Paper with no name yet

### Author list

Preprint : no preprint yet

doi: no doi yet

## 1 - Setup

TODO

Prepare requirements.txt

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

TODO

## 4 - License

This code repository is under [GPL v3.0 License](https://github.com/fabda/mcc_paper/blob/master/LICENSE)
