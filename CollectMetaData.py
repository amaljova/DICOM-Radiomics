'''
* @ Date: 06 December 2021
* @ DateMod: 02 December 2022
* @ author: Amal Joseph Varghese
* @ email: amaljova@gmail.com
* @ github: https://github.com/amaljova
* @ last modified: 21 July 2022

This script reads all DICOM files in a directory/ directories and collects their metadata into .csv and  .pkl files.
Edit the final part before you run.


This scripts deals with RTSTRUCT and CT modalities.
Edit to include MR or PET or even to find patients of interest.
'''
# =========================================Need Not Modify Block=======================================
import multiprocessing
import pandas as pd
from pydicom import dcmread
import os
import pickle as pkl
# import threading
import multiprocessing

# From StackOverfloe: To measure time
from datetime import datetime
import time
_start_time = time.time()
now = datetime.now()
print("Current Time =", now.strftime("%H:%M:%S"))
def tic():
    global _start_time
    _start_time = time.time()
def tac():
    t_sec = round(time.time() - _start_time)
    (t_min, t_sec) = divmod(t_sec, 60)
    (t_hour, t_min) = divmod(t_min, 60)
    print('Time passed: {}hour:{}min:{}sec'.format(t_hour, t_min, t_sec))

# ========================================= CODE ======================================================


def getReferencedCTUID(header):
    if len(list(header[0x3006, 0x10])) > 0:
        refFrameOfRef = (header[0x3006, 0x10])[0]
        if len(list(refFrameOfRef[0x3006, 0x0012])) > 0:
            rtRefStudy = (refFrameOfRef[0x3006, 0x0012])[0]
            if len(list(rtRefStudy[0x3006, 0x14])) > 0:
                rtRefSerie = (rtRefStudy[0x3006, 0x14])[0]
                return rtRefSerie[0x20, 0xe].value
    return None


def getData(base=''):
    data_list=[]
    for root, dirs, files in os.walk(base):
        for file in files:
            try:
                f_path = os.path.join(root, file)
                header = dcmread(f_path)
                pat_id = None  # Patient ID
                pat_id = header[(0x10, 0x20)].value  # Patient ID
                pat_Name = None  # Patient Name
                pat_Name = header[(0x10, 0x10)].value  # Patient Name
                modality = None  # Modality
                modality = header[(0x08, 0x60)].value  # Modality
                # Default Variable values
                SliceThickness = None  # Slice Thickness
                stu_inst_UID = None  # Study Instance UID
                ser_inst_UID = None  # Series Instance UID
                StudyDescription = None  # StudyDescription
                SeriesDescription = None  # SeriesDescription
                SOPClassUID = None  # SOP Class UID
                SOPInstanceUID = None  # SOP Instance UID
                StructureSetLabel = None  # Structure Set Label
                ReferencedCTUID = None  # ReferencedCTUID
                ROIs = None  # ROIs

                # _______Operations on file_____________________
                try:
                    # Study Instance UID
                    try:
                        stu_inst_UID = header[(0x20, 0x0d)].value
                    except Exception as e:
                        pass
                        # print(e, "stu_inst_UID", "Patient: ", pat_id,
                        #       "Modality: ", modality)
                    # Series Instance UID
                    try:
                        ser_inst_UID = header[(0x20, 0x0e)].value
                    except Exception as e:
                        pass
                        # print(e, "ser_inst_UID", "Patient: ", pat_id,
                        #       "Modality: ", modality)
                    # StudyDescription
                    try:
                        StudyDescription = header[(0x08, 0x1030)].value
                    except Exception as e:
                        pass
                        # print(e, "StudyDescription", "Patient: ", pat_id,
                        #       "Modality: ", modality)
                    # SeriesDescription
                    try:
                        SeriesDescription = header[(0x08, 0x103e)].value
                    except Exception as e:
                        pass
                        # print(e, "SeriesDescription", "Patient: ", pat_id,
                        #       "Modality: ", modality)
                    # SOP Class UID
                    try:
                        SOPClassUID = header[(0x08, 0x16)].value
                    except Exception as e:
                        pass
                        # print(e, "SOPClassUID", "Patient: ", pat_id,
                        #       "Modality: ", modality)
                    # SOP Instance UID
                    try:
                        SOPInstanceUID = header[(0x08, 0x18)].value
                    except Exception as e:
                        pass
                        # print(e, "SOPInstanceUID", "Patient: ", pat_id,
                        #       "Modality: ", modality)

                    # Structure Set Label
                    try:
                        StructureSetLabel = header[(
                            0x3006, 0x02)].value
                    except Exception as e:
                        pass
                        # print(e, 'StructureSetLabel', 'Patient: ',
                        #       pat_id, ' Modality: ', modality)

                    # Slice Thickness
                    try:
                        SliceThickness = header[(0x18, 0x50)].value
                    except Exception as e:
                        pass
                        # print(e, "Slice Thickness", "Patient: ",
                        #       pat_id, "Modality: ", modality)
                    
                    #ReferencedCTUID
                    try:
                        ReferencedCTUID = getReferencedCTUID(header)
                    except Exception as e:
                        pass
                        # print(e, ' Pixel Spacing', 'Patient: ',
                        #       pat_id, ' Modality: ', modality)
                    #ROIs
                    try:
                        ROIs = [
                            i[(0x3006, 0x26)].value for i in header[(0x3006, 0x20)].value]
                    except Exception as e:
                        pass
                        # print(e, ' Pixel Spacing', 'Patient: ',
                        #       pat_id, ' Modality: ', modality)

                    data_list.append((
                        pat_id,
                        pat_Name,
                        stu_inst_UID,
                        StudyDescription,
                        ser_inst_UID,
                        SeriesDescription,
                        ReferencedCTUID,
                        StructureSetLabel,
                        modality,
                        SOPClassUID,
                        SOPInstanceUID,
                        SliceThickness,
                        f_path,
                        ROIs
                    ))
                # print('Data Entered.: ', pat_id)
                except Exception as e:
                    print("File level Error")
                    print(e)
                    print('Patient: ', pat_id, ' Modality: ', modality)
                    continue
            except Exception as e:
                print(e, "\nfile is not Readable")
                continue
    return data_list




def makeDataBse(path,index):
    columns = [
        'Patient_ID',
        'pat_Name',
        'Study_Instance_UID',
        'Study_Description',
        'Series_Instance_UID',
        'Series_Description',
        'ReferencedSeriesUID',
        'Structure_Set_Label',
        'Modality',
        'SOP_Class_UID',
        'SOP_Instance_UID',
        'Slice_Thickness',
        'File_Path',
        'ROIs'
    ]

    print(f"Started: {path}")
    data_list = getData(path)
    print(f"Finished: {path}")

    data = pd.DataFrame(data_list, index=None, columns=columns)
    data.to_csv(f"{index}_outfile.csv", index=None)
    print('CSV_file --done!')
    data = {'columns':columns,"data":data_list}
    with open(f"{index}_output.pkl",'wb') as f:
        pkl.dump(data,f)

# =========================================FIXME Block=======================================


if __name__ == '__main__':
    tic()

    # modify this list to include all your source directories.
    source = [
        "absolute_path_to_the_target_folder"
    ]

    thrs = dict()
    for i_num,sr in enumerate(source):
        thrs[i_num] = multiprocessing.Process(target=makeDataBse, args=(sr,i_num))

    
    for k in thrs:
        thrs[k].start()
        print("started thread: ",k)
    for k in thrs:
        thrs[k].join()
        print("join thread: ",k)

    print('Done!')
    tac()
