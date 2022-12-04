import os
import shutil
import pandas as pd


"""
@ Date: 15 November 2022
@ DateMod: 02 December 2022
@ author: AMAL JOSEPH VARGHESE
@ email: amaljova@gmail.com
@ github: https://github.com/amaljova


This script copies DICOM data to a newly created, properly arrange folder.
Run CollectMetaData.py and generate the required metadata csv file before running it.

"""


def makeFolders(dir_name):
    '''To create a directory/directories if it/they are basent'''
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

def getpath(x):
    return x.split("\\")[-1]

def copyFiles(metadata_file,dest_dir):
    datasheet = pd.read_csv(metadata_file)
    for i, row in datasheet.iterrows():
                try:
                    # -------------------data Copying---------------------------------
                    dest_dir = os.path.join(dest_dir, row['Patient_ID'])
                    dest_dir = os.path.join(dest_dir, row['Study_Instance_UID'])
                    dest_dir = os.path.join(dest_dir, row['Series_Instance_UID'])
                    makeFolders(dest_dir)
                    dest_file = os.path.join(dest_dir, getpath(row['File_Path']))
                    shutil.copy(row['File_Path'], dest_file)
                except Exception as e:
                    print(e)

if __name__ == "__main__":

    # Edit here

    metadata_file = "0_outfile.csv"
    dest_dir = "dest_folder"

    print("Started.")
    copyFiles(metadata_file,dest_dir)
    print("Done.")