# Script          : refocus.py
# Code from Paper : ?
# Author          : Fabrice Daian - 09/2020
# License         : GPL v3.0
# Github repo     : https://github.com/fabda/mcc_paper

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from skimage import io
from tifffile import imsave
import argparse
import sys
import os


def read_imagej_manualtracking(file,delimiter=","):
    """
        Read an imageJ Manual Tracking CSV file and create a Panda DataFrame
        containing two columns X and Y
        Arguments:
            - file : CSV filename
            - delimiter: delimiter used to read the CSV file (default is ",")
        Return:
            - X coordinates
            - Y coordinates
            - Pandas dataframe containing the two columns X and Y.
    """
    df = pd.read_csv(file,delimiter=delimiter)
    X = np.array(df["X"])
    Y = np.array(df["Y"])
    return X,Y,df


def read_imagej_hyperstack(file):
    """
        Read an ImageJ Hyperstack and return the red and green channel as
        separated image stack
        Argument:
            - file: TIFF Hyper stack filename
        Return:
            Red and Green channel stacks
    """
    im0      = io.imread(file)
    im_red   = im0[:,0,:,:]
    im_green = im0[:,1,:,:]
    return im_red,im_green


def get_drift(X,Y):
    """
        Quantify the horizontal and vertical drift using a set of X,Y points
        Arguments:
            - X,Y : points
        Return:
            - the drift calculated for each frame
    """

    drift = np.transpose(np.array([np.diff(X),np.diff(Y)]))
    drift[:,0]=np.cumsum(drift[:,0])
    drift[:,1]=np.cumsum(drift[:,1])

    return drift


def refocus_one_channel(stack,drift,filename_save=""):
    """
        Refocus an one channel from ImageJ Hyperstackby by using a drift value (returned from get_drift() function)
        Arguments:
            - stack         : One color channel stack (returned by read_imagej_hyperstack)
            - X,y           : set of manually tracked point coordinates (returned by read_imagej_manualtracking)
            - filename_save : name of the refocused stack filename (TIFF file). Empty by default
                              meaning no file is saved. Otherwise the file is saved using this filename_save parameter
        Return;
            - the refocused stack

    """
    A = []
    A.append(stack[0])
    for indix in range(1,stack.shape[0]):
        if indix>drift.shape[0]:
            break
        tmp = np.copy(stack[indix])
        tmp = np.roll(stack[indix],-drift[indix-1][0],axis=1)
        tmp = np.roll(tmp,-drift[indix-1][1],axis=0)
        A.append(tmp)
    A = np.array(A)
    if len(filename_save)>0:
        imsave(filename_save, A)
    return A


def define_arguments():
    """
        ArgumentParser to retrieve and check the validity of the script arguments
        Returns:
            -  stack_name   : path to the TIFF stack
            - tracking_file : path to the CSV file containing the tracking points
    """

    my_parser = argparse.ArgumentParser(description='Refocus an ImageJ TIFF hyperstack')

    my_parser.add_argument('--stack_name',
                           metavar='stack_name',
                           type=str,
                           help='path to the TIFF hyperstack')

    my_parser.add_argument('--tracking_file',
                           metavar='tracking_file',
                           type=str,
                           help='path to the CSV file containing the tracking')

    args = my_parser.parse_args()


    if not os.path.isdir(args.stack_name):
        print('The path specified for the stack_name does not exist')
        sys.exit()

    if not os.path.isdir(args.tracking_file):
        print('The path specified for the tracking_file does not exist')
        sys.exit()

    return args.stack_name,args.tracking_file


#  -- Main --
if __name__ == "__main__":

    # Get file paths as script as script arguments
    stack_name,tracking_file = define_arguments()

    # Read TIFF Stack and return separate Green and Red Channel
    R, G    = read_imagej_hyperstack(stack_name)

    # Read CSV Tracking gile and return the set of points
    X,Y,df  = read_imagej_manualtracking(tracking_file)

    # Calculate the drift
    drift   = get_drift(X,Y)

    # We use the same drift to refocus both Red and Green channel
    refocus_one_channel(R, drift, "00A.tif")
    refocus_one_channel(G, drift, "00B.tif")
