# Script          : analysis.py
# Code from Paper : ?
# Author          : Fabrice Daian - 09/2020
# License         : GPL v3.0
# Github repo     : https://github.com/fabda/mcc_paper

import numpy as np
import pandas as pd
from scipy import interpolate
from scipy.signal import savgol_filter
from scipy.ndimage import gaussian_filter1d
import math
import argparse
import os
import sys

def read_imagej_object_tracking(file,delimiter="\t"):
    """
        Read an imageJ Manual Object Tracking CSV file and create a Panda DataFrame
        Arguments:
            - file : CSV filename
            - delimiter: delimiter used to read the CSV file (default is "\t")
        Return:
            - Pandas dataframe containing the CSV file.
    """
    df = pd.read_csv(file,delimiter=delimiter)
    return df

def get_object_X_Y(object_id,df):
    """
        Extract X and Y coordinates from a given object from a given CSV file
        Arguments:
            - object_id : the unique tracked object into the Pandas dataframe
            - df        : the pandas dataframe containing the CSV file
        Return:
            - slices    : Every slice number where the object is tracked
            - X         : X coordinates of the object
            - Y         : Y coordinates of the object
    """
    # Read
    tmp = df[df["Track"]==object_id]
    slices = np.array(tmp["Slice"])
    X = np.array(tmp["X"])
    Y = np.array(tmp["Y"])
    return slices,X,Y

def get_distances(X,Y):
    """
        Calculate the distance between two consecutive positions given by X and Y
        Arguments:
            - X  : X coordinates of the object
            - Y  : Y coordinates of the object
        Return:
            - D  : euclidean distances between two consecutive positions
    """
    # Distance
    D=[0]
    for i in range(1,X.shape[0]):
        D.append(np.sqrt((X[i]-X[i-1])**2 + (Y[i]-Y[i-1])**2))
    D=np.array(D)
    return D

def get_cumulative_distances(X,Y):
    """
        Calculate the cumulative distance between two consecutive positions given by X and Y
        Arguments:
            - X  : X coordinates of the object
            - Y  : Y coordinates of the object
        Return:
            - csD  : Cumulative euclidean distances between two consecutive positions

    """
    D = get_distances(X,Y)
    return np.cumsum(D)

def get_orientation_toward_origin(X,Y):
    """
        Calculate the cumulative distance between two consecutive positions given by X and Y
        Arguments:
            - X       : X coordinates of the object
            - Y       : Y coordinates of the object
        Return:
            - angles  : angles between the coordinates and the origin

    """
    angles=[]
    for i in range(X.shape[0]):
        angle = (180*math.atan2(Y[i],X[i]))/np.pi
        angles.append(angle)
    angles = np.array(angles)
    return angles

def get_direction(X,Y,window=1):
    """
        Calculate the cumulative distance between two consecutive positions given by X and Y
        Arguments:
            - X           : X coordinates of the object
            - Y           : Y coordinates of the object
            - window      : window size of the direction calculation (default window=1 means two consecutive position)
        Return:
            - directions  : directions of the vector calculated from two X and Y coordinate in the window

    """
    directions=[]
    for i in range(window,X.shape[0]):
        direction = (180*math.atan2(Y[i]-Y[i-window],X[i]-X[i-window]))/np.pi
        directions.append(direction)
    directions = np.array(directions)
    return directions


def calculate_object_instant_velocity(object_id,df,sigmaG=10):

    """
        Calculate the instant veolocity of an object at every frame
        Arguments:
            - object_id : unique identifier of a tracked object into the CSV file
            - df        : Pandas dataframe containing the CSV file
            - sigmaG    : sigma  of the 1-d gaussian filter apply on both X and Y positions (default sigmaG = 10)
        Return:
            - slices    : Every slice number where the object is tracked
            - Vins      : Instant velocity calculated at every frame for the object
    """

    slices,X,Y = get_object_X_Y(object_id,df)

    # 1-d gaussian filter
    X = gaussian_filter1d(X,sigmaG)
    Y = gaussian_filter1d(Y,sigmaG)

    csD = get_cumulative_distances(X,Y)

    # Instant velocity
    Vins=[0]
    for i in range(1,csD.shape[0]-1):
        #Vins.append(np.sum(D[i:i+3])/2)
        Vins.append((csD[i+1]-csD[i-1])/2)
    Vins.append(0)
    Vins=np.array(Vins)

    # Return data
    return slices,Vins

def calculate_object_travelled_distance(object_id,df,sigmaG=10):
    """
        Calculate the distance travelled by an object at every frame
        Arguments:
            - object_id : unique identifier of a tracked object into the CSV file
            - df        : Pandas dataframe containing the CSV file
            - sigmaG    : sigma  of the 1-d gaussian filter apply on both X and Y positions (default sigmaG = 10)
        Return:
            - slices    : Every slice number where the object is tracked
            - D         : Distance travelled calculated at every frame for the object
    """
    slices,X,Y = get_object_X_Y(object_id,df)

    # 1-d gaussian filter
    X = gaussian_filter1d(X,sigmaG)
    Y = gaussian_filter1d(Y,sigmaG)

    D = get_distances(X,Y)

    return slices,D

def calculate_object_cumulative_travelled_distance(object_id,df,sigmaG=10):
    """
        Calculate the distance travelled by an object at every frame
        Arguments:
            - object_id : unique identifier of a tracked object into the CSV file
            - df        : Pandas dataframe containing the CSV file
            - sigmaG    : sigma  of the 1-d gaussian filter apply on both X and Y positions (default sigmaG = 10)
        Return:
            - slices    : Every slice number where the object is tracked
            - csD       : Cumulative distance travelled calculated at every frame for the object
    """
    slices,X,Y = get_object_X_Y(object_id,df)

    # 1-d gaussian filter
    X = gaussian_filter1d(X,sigmaG)
    Y = gaussian_filter1d(Y,sigmaG)

    csD = get_cumulative_distances(X,Y)

    return slices,csD


def calculate_object_orientation_toward_origin(object_id,df,sigmaG=10):
    """
        Calculate the angle between the object and the origin
        Arguments:
            - object_id : unique identifier of a tracked object into the CSV file
            - df        : Pandas dataframe containing the CSV file
            - sigmaG    : sigma  of the 1-d gaussian filter apply on both X and Y positions (default sigmaG = 10)
        Return:
            - slices    : Every slice number where the object is tracked
            - angles    : angles between the object and the origin at every frame
    """
    slices,X,Y = get_object_X_Y(object_id,df)

    # 1-d gaussian filter
    X = gaussian_filter1d(X,sigmaG)
    Y = gaussian_filter1d(Y,sigmaG)

    angles = get_orientation_toward_origin(X,Y)

    return slices,angles

def calculate_object_direction(object_id,df,sigmaG=10,window=1):
    """
        Calculate the angle between the object and the origin
        Arguments:
            - object_id : unique identifier of a tracked object into the CSV file
            - df        : Pandas dataframe containing the CSV file
            - sigmaG    : sigma  of the 1-d gaussian filter apply on both X and Y positions (default sigmaG = 10)
            - window      : window size of the direction calculation (default window=1 means two consecutive position)
        Return:
            - slices    : Every slice number where the object is tracked
            - directions  : directions of the vector calculated from two positions coordinate in the window
    """
    slices,X,Y = get_object_X_Y(object_id,df)

    # 1-d gaussian filter
    X = gaussian_filter1d(X,sigmaG)
    Y = gaussian_filter1d(Y,sigmaG)

    angles = get_direction(X,Y,window=window)

    # Window correction
    angles_corrected=np.zeros_like(slices)
    angles_corrected[0:angles.shape[0]]=angles

    return slices,angles_corrected




def generate_instant_velocity_result(df,sigmaG=10):
    """
        Concatenate every individual object instant velocity into ONE Pandas Dataframe
        Arguments:
            - df        : Pandas dataframe containing the CSV file
            - sigmaG    : sigma  of the 1-d gaussian filter apply on both objects X and Y positions (default sigmaG = 10)
        Return:
            - df_result : Pandas dataframe containing the concatenation of every individual object instant velocity
    """
    df_result = pd.DataFrame()

    for object_id in np.unique(np.array(df["Track"])):
        slices,speed = calculate_object_instant_velocity(object_id,df,sigmaG=sigmaG)

        obj = (np.zeros_like(slices)+object_id).astype('int')
        df0 = pd.DataFrame()
        df0["Cell"]=obj
        df0["Frames"]=slices
        df0["Instant_speed"]=speed

        df_result = pd.concat([df_result,df0])

    return df_result


def generate_distance_result(df,sigmaG=10):
    """
        Concatenate every individual object distance into ONE Pandas Dataframe
        Arguments:
            - df        : Pandas dataframe containing the CSV file
            - sigmaG    : sigma  of the 1-d gaussian filter apply on both objects X and Y positions (default sigmaG = 10)
        Return:
            - df_result : Pandas dataframe containing the concatenation of every object consecutive distance
    """
    df_result = pd.DataFrame()

    for object_id in np.unique(np.array(df["Track"])):
        slices,D = calculate_object_travelled_distance(object_id,df,sigmaG=sigmaG)

        obj = (np.zeros_like(slices)+object_id).astype('int')
        df0 = pd.DataFrame()
        df0["Cell"]=obj
        df0["Frames"]=slices
        df0["Distance"]=D

        df_result = pd.concat([df_result,df0])

    return df_result

def generate_cumulative_distance_result(df,sigmaG=10):
    """
        Concatenate every individual object cumulative distance into ONE Pandas Dataframe
        Arguments:
            - df        : Pandas dataframe containing the CSV file
            - sigmaG    : sigma  of the 1-d gaussian filter apply on both objects X and Y positions (default sigmaG = 10)
        Return:
            - df_result : Pandas dataframe containing the concatenation of every object cumulative distance
    """
    df_result = pd.DataFrame()

    for object_id in np.unique(np.array(df["Track"])):
        slices,csD = calculate_object_cumulative_travelled_distance(object_id,df,sigmaG=sigmaG)

        obj = (np.zeros_like(slices)+object_id).astype('int')
        df0 = pd.DataFrame()
        df0["Cell"]=obj
        df0["Frames"]=slices
        df0["Cumulative_distance"]=csD

        df_result = pd.concat([df_result,df0])

    return df_result

def generate_orientation_toward_origin_result(df,sigmaG=10):
    """
        Concatenate every individual object orientation toward origin into ONE Pandas Dataframe
        Arguments:
            - df        : Pandas dataframe containing the CSV file
            - sigmaG    : sigma  of the 1-d gaussian filter apply on both objects X and Y positions (default sigmaG = 10)
        Return:
            - df_result : Pandas dataframe containing the concatenation of every object orientation toward origin
    """
    df_result = pd.DataFrame()

    for object_id in np.unique(np.array(df["Track"])):
        slices,orientations = calculate_object_orientation_toward_origin(object_id,df,sigmaG=sigmaG)

        obj = (np.zeros_like(slices)+object_id).astype('int')
        df0 = pd.DataFrame()
        df0["Cell"]=obj
        df0["Frames"]=slices
        df0["Orientation_toward_origin"]=orientations

        df_result = pd.concat([df_result,df0])

    return df_result


def generate_direction_result(df,sigmaG=10,window=1):
    """
        Concatenate every individual object directions into ONE Pandas Dataframe
        Arguments:
            - df        : Pandas dataframe containing the CSV file
            - sigmaG    : sigma  of the 1-d gaussian filter apply on both objects X and Y positions (default sigmaG = 10)
            - window    : window size of the direction calculation (default window=1 means two consecutive position)
        Return:
            - df_result : Pandas dataframe containing the concatenation of every object directions
    """
    df_result = pd.DataFrame()

    for object_id in np.unique(np.array(df["Track"])):
        slices,directions = calculate_object_direction(object_id,df,sigmaG=sigmaG,window=window)

        obj = (np.zeros_like(slices)+object_id).astype('int')
        df0 = pd.DataFrame()
        df0["Cell"]=obj
        df0["Frames"]=slices
        df0["Direction"]=directions

        df_result = pd.concat([df_result,df0])

    return df_result



def generate_result(features,filename = None):
    """
        Generate the complete Pandas DataFrame file containing every individual analysis dataframe
        provided into the features list of pandas dataframe
        Arguments:
            - features : a list containing individual pandas dataframe result (ex: instant velocity, distance ...)
            - filename : a XLS filename to write the results (optional) (default is None if the result is not exported on an XLS file)
        Return :
            - DF       : Complete Pandas Dataframe gathering all the analysis
    """
    DF=pd.DataFrame()
    for df0 in features:
        if DF.shape[0]==0:
            DF = df0
        else:
            DF = pd.merge(DF, df0, on=["Cell","Frames"])

    if filename is not None:
        pd.to_excel(filename)
        print(filename , " has been saved successfully !")

    return DF


def define_arguments():
    """
        ArgumentParser to retrieve and check the validity of the script arguments
        Returns:
            - stack_name    : path to the TIFF stack
            - tracking_file : path to the CSV file containing the tracking points
    """

    my_parser = argparse.ArgumentParser(description='Refocus an ImageJ TIFF hyperstack')

    my_parser.add_argument('--tracking_file',
                           metavar='tracking_file',
                           type=str,
                           help='path to the tracking CSV file')

    my_parser.add_argument('--sigma',
                           metavar='sigma',
                           type=int,
                           help='sigma value for 1-d gaussian filter')

    my_parser.add_argument('--window',
                           metavar='window',
                           type=int,
                           help='window size for direction calculation')

    my_parser.add_argument('--pixmicron',
                           metavar='pixmicron',
                           type=float,
                           help='1 pixel equals ? µm')

    my_parser.add_argument('--framesec',
                           metavar='framesec',
                           type=int,
                           help='1 frame equals ? s')


    my_parser.add_argument('--result_file',
                           metavar='result_file',
                           type=str,
                           help='CSV result filename')

    args = my_parser.parse_args()


    if not os.path.isfile(args.tracking_file):
        print('The path specified for the tracking_file does not exist')
        sys.exit()

    return args.tracking_file, args.sigma,args.window, args.pixmicron, args.framesec, args.result_file



# -- main --
if __name__ == "__main__":

        """
            Examples of use:

            $> python analysis.py
                        --sigma=10
                        --window=5
                        --pixmicron=0.4
                        --framesec=60
                        --tracking_file="../Results from STABILIZED_MAX_20180517_40X_pa-Tub-LA-GFP__mRFP_st14-2 in µm per sec.csv"
                        --result_file=result.csv

        """

        # Get arguments
        tracking_file, sigmaG, window, pixmicron, framesec, result_file = define_arguments()

        df           = read_imagej_object_tracking(tracking_file,delimiter=",")

        df_speed     = generate_instant_velocity_result(df,sigmaG=sigmaG)
        df_dist      = generate_distance_result(df,sigmaG=sigmaG)
        df_cumdist   = generate_cumulative_distance_result(df,sigmaG=sigmaG)
        df_orient    = generate_orientation_toward_origin_result(df,sigmaG=sigmaG)
        df_direction = generate_direction_result(df,sigmaG=sigmaG,window=window)

        res          = generate_result([df_speed, df_dist, df_cumdist,df_orient,df_direction])

        print(res)






#
# X,Y,df = read_imagej_manualtracking(filename,delimiter=",")
#
# print(X.shape)
# print(Y.shape)
# print(df.shape)
#
# # --- CALCUL DE L'ORIENTATION
# # track  = 1
# # sigmaG = 10
# # tmp = df[df["Track"]==track]
# # X = np.array(tmp["X"])
# # Y = np.array(tmp["Y"])
# # X = gaussian_filter1d(X,sigmaG)
# # Y = gaussian_filter1d(Y,sigmaG)
# #
# # ANGLES=[]
# # for i in range(10,X.shape[0]):
# #
# #     x1 = X[i-10]
# #     x2 = X[i]
# #     y1 = Y[i-10]
# #     y2 = Y[i]
# #     # x1 = 0
# #     # x2 = X[i]
# #     # y1 = 0
# #     # y2 = Y[i]
# #
# #     angle = (180*math.atan2(y2-y1,x2-x1))/np.pi
# #     ANGLES.append(angle)
# #
# # ANGLES = np.array(ANGLES)
# #
# # print(ANGLES)
# #
# #
# # import matplotlib.pyplot as plt
# #
# # #plt.plot(np.abs(ANGLES))
# # # plt.plot(ANGLES)
# # # plt.show()
# #
# # D,Vins = instspeed(track,df,sigmaG)
# #
# # #plt.plot(D)
# # #plt.show()
# #
# # # plt.plot(Vins)
# # # plt.show()
# #
# #
# # fig,ax = plt.subplots()
# # ax.plot(ANGLES,c="r")
# # ax.set_xlabel("Frames",fontsize=14)
# # ax.set_ylabel("Angles in degree",color="red",fontsize=14)
# #
# # ax2=ax.twinx()
# # ax2.plot(Vins,c="b")
# # ax2.set_ylabel("Instant Speed",color="blue",fontsize=14)
# # plt.show()
#
#
# # ANGLE ENTRE POINT ET DISTANCE
# # sigmaG = 10
# #
# # track  = 1
# # tmp = df[df["Track"]==track]
# # X1 = gaussian_filter1d(np.array(tmp["X"]),sigmaG)
# # Y1 = gaussian_filter1d(np.array(tmp["Y"]),sigmaG)
# # D1,Vins1 = instspeed(track,df,sigmaG)
# #
# # track  = 2
# # tmp = df[df["Track"]==track]
# # X2 = gaussian_filter1d(np.array(tmp["X"]),sigmaG)
# # Y2 = gaussian_filter1d(np.array(tmp["Y"]),sigmaG)
# # D2,Vins2 = instspeed(track,df,sigmaG)
# #
# # print(X1.shape)
# # print(X2.shape)
# #
# # ANGLES=[]
# # DIST = []
# # FRAMES=[]
# # for i in range(X1.shape[0]):
# #     x1 = X1[i]
# #     x2 = X2[i]
# #     y1 = Y1[i]
# #     y2 = Y2[i]
# #     angle = (180*math.atan2(y2-y1,x2-x1))/np.pi
# #     # ANGLES.append(np.abs(angle))
# #     ANGLES.append(angle)
# #     FRAMES.append(i)
# #
# # ANGLES=np.array(ANGLES)
# #
# # for i in range(X1.shape[0]):
# #     DIST.append(np.sqrt((X1[i]-X2[i])**2 + (Y1[i]-Y2[i])**2))
# # DIST=np.array(DIST)
# #
# # import matplotlib.pyplot as plt
# #
# # plt.plot(Vins1,c="r")
# # plt.plot(Vins2,c="b")
# # plt.show()
# #
# #
# # plt.plot(ANGLES)
# # plt.show()
# #
# # plt.plot(DIST)
# # plt.show()
# #
# # plt.scatter(ANGLES,DIST,s=5,c=FRAMES)
# # plt.xlabel("angles")
# # plt.ylabel("distance")
# # plt.colorbar()
# # plt.show()
#
#
# # ANGLE ENTRE POINT ET ORIGINE ET CORRELATION AVEC LA DISTANCE ENTRE LES POINTS
#
#
# sigmaG = 10
#
# #Recap all position
# # import matplotlib.pyplot as plt
# # for i in range(1,3):
# #     track  = i
# #     tmp = df[df["Track"]==track]
# #     X1 = gaussian_filter1d(np.array(tmp["X"]),sigmaG)
# #     Y1 = gaussian_filter1d(np.array(tmp["Y"]),sigmaG)
# #     D1,Vins1 = instspeed(track,df,sigmaG)
# #     plt.plot(X1,Y1)
# #     plt.text(X1[0],Y1[0],str(i))
# # plt.show()
#
#
# def get_data(track,df,sigmaG):
#     tmp = df[df["Track"]==track]
#     X = gaussian_filter1d(np.array(tmp["X"]),sigmaG)
#     Y = gaussian_filter1d(np.array(tmp["Y"]),sigmaG)
#     return X,Y
#
#
# def orientation_origin(track,df,sigmaG):
#     X,Y = get_data(track,df,sigmaG)
#     ANGLES=[]
#     for i in range(X.shape[0]):
#         x1 = 0
#         x2 = X[i]
#         y1 = 0
#         y2 = Y[i]
#         angle = (180*math.atan2(y2-y1,x2-x1))/np.pi
#         ANGLES.append(angle)
#     ANGLES = np.array(ANGLES)
#     return ANGLES
#
# def direction(track,df,sigmaG,step=5):
#     X,Y = get_data(track,df,sigmaG)
#     ANGLES=[]
#     for i in range(step,X.shape[0]):
#         x1 = X[i-step]
#         x2 = X[i]
#         y1 = Y[i-step]
#         y2 = Y[i]
#         angle = (180*math.atan2(y2-y1,x2-x1))/np.pi
#         ANGLES.append(angle)
#     ANGLES = np.array(ANGLES)
#     return ANGLES
#
#
# def euclidean_distance(track1,track2,df,sigmaG):
#     X1,Y1 = get_data(track1,df,sigmaG)
#     X2,Y2 = get_data(track2,df,sigmaG)
#     DIST=[]
#     for i in range(X1.shape[0]):
#         DIST.append(np.sqrt((X1[i]-X2[i])**2 + (Y1[i]-Y2[i])**2))
#     DIST=np.array(DIST)
#     return DIST
#
# ANGLES1 = orientation_origin(2,df,sigmaG)
# ANGLES2 = orientation_origin(6,df,sigmaG)
#
# D1,Vins1 = instspeed(2,df,sigmaG)
# D2,Vins2 = instspeed(6,df,sigmaG)
#
# DIST     = euclidean_distance(2,6,df,sigmaG)
# print(DIST)
#
#
# import matplotlib.pyplot as plt
#
# fig,ax = plt.subplots()
# #ax.plot(ANGLES1,"r")
# #ax.plot(ANGLES2,"b")
# ax.plot(np.abs(ANGLES1-ANGLES2),c="k")
# ax.set_ylabel("Difference between cell Angles (from origin)",color="k")
# ax2=ax.twinx()
# ax2.plot(DIST,c="g")
# ax2.set_ylabel("Distance between cells",color="g")
# ax2.set_ylim([0,np.max(DIST)+10])
# plt.xlabel("Frames")
# plt.show()
#
#
# plt.plot(Vins1,c="r")
# plt.plot(Vins2,c='b')
# plt.show()
#
#
# DIRECTION1 = direction(2,df,sigmaG,step=10)
# DIRECTION2 = direction(6,df,sigmaG,step=10)
# plt.plot(DIRECTION1,c="r")
# plt.plot(DIRECTION2,c="b")
# plt.xlabel("Frames")
# plt.ylabel("Direction in degree (10 frames)")
# plt.show()
#
#
#
# X1,Y1 = get_data(2,df,sigmaG)
# X2,Y2 = get_data(6,df,sigmaG)
#
#
# for i in np.arange(0,Vins1.shape[0],10):
#     plt.plot(X1,Y1,c="r")
#     plt.plot(X2,Y2,c="b")
#     plt.scatter(X1[i],Y1[i],c="red",s=5)
#     plt.scatter(X2[i],Y2[i],c="blue",s=5)
#     plt.plot([X1[i],X2[i]],[Y1[i],Y2[i]],c="gray")
# plt.show()
#
#
# plt.plot(Vins1[10:Vins1.shape[0]],DIRECTION1)
# plt.show()
#
# plt.plot(Vins2[10:Vins2.shape[0]],DIRECTION2)
# plt.show()
#
#
#
# import sys
# sys.exit()
#
#
#
# #Paramètrage en fonction du film
#
# pixmicron = 0.4   # 1 pixel = 0.4µm
# sigmaG = 10 #Gaussian smoothing to correct manual tracking effect
#
# # if file0 == "Results from STABILIZED_SUM_20190117_40X_paTubLA-GFP_mRFP_MO-SCF-blue_st14-1_w1561_t1-1 in µm per sec":
# #     frame_per_second = 1/30  #meaning 1 frame = 120s (bizarre mais c'est comme ça)
# # elif file0 == "Results from STABILIZED_MAX_20180214_60X_moKIT_s+l_10ng_pa-tub-LA-GFP_st14-2 in µm per sec":
# #     frame_per_second = 1/120 #meand 1 frame = 30s
# #     pixmicron = 0.266667
# # elif file0 == "Results from STABILIZED_MAX_20180214_60X_moKIT_s+l_10ng_pa-tub-LA-GFP_st14-1 in µm per sec":
# #     frame_per_second = 1/120 #meand 1 frame = 30s
# #     pixmicron = 0.266667
# # elif file0 == "Results from STABILIZED_MAX_20180214_60X_moKIT_s+l_10ng_pa-tub-LA-GFP_st14-2-rfp-1 in µm per sec":
# #     frame_per_second = 1/120 #meand 1 frame = 30s
# #     pixmicron = 0.266667
# # else:
# #     frame_per_second = 1/60  #meaning 1 frame = 60s
#
# frame_per_second = 1/60  #meaning 1 frame = 60s
#
# print("FPS utilisé:",frame_per_second)
# print("Pixmicron utilisé:",pixmicron)
# print("Sigma Gaussian utilisé:", sigmaG)
#
# df_result = pd.DataFrame()
# for i in list(np.unique(df["Track"])):
#     D,Vins = instspeed(i,df,sigmaG)
#     df0=pd.DataFrame()
#     df0 = generate_result_csv(i,D,Vins,pixmicron,frame_per_second)
#     #print(df0.head())
#     if df_result.shape[0]==0:
#         df_result = df0.copy()
#     else:
#         df_result = pd.concat([df_result, df0], axis=0)
# #     gengraph_distance(D,i,pixmicron,frame_per_second)
# #     gengraph_instspeed(Vins,i,pixmicron,frame_per_second)
#
#
# # df_result.to_excel("ANALYZE_FROM" + file0 + ".xls")
