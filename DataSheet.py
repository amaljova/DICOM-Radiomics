'''
* @ Date: 06 December 2021
* @ author: Amal Joseph Varghese
* @ email: amaljova@gmail.com
* @ github: https://github.com/amaljova
* @ last modified: 06 December 2021


Make a datasheet of all DICOM files in a folder or in many folders.
This scripts deals with RTSTRUCT and CT modalities.
Edit this to include MR or PET or even to find patients of interest.
'''
# =========================================Need Not Modify Block=======================================

import pandas as pd
from pydicom import dcmread
import os

# -----------------------------------------------------------------------------------------------------

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
    if len(list(header[0x3006,0x10])) > 0:
        refFrameOfRef = (header[0x3006,0x10])[0]
        if len(list(refFrameOfRef[0x3006, 0x0012])) > 0:
            rtRefStudy = (refFrameOfRef[0x3006, 0x0012])[0]
            if len(list(rtRefStudy[0x3006,0x14])) > 0:
                rtRefSerie = (rtRefStudy[0x3006,0x14])[0]
                return rtRefSerie[0x20,0xe].value
    return None

def theFunction(base='', modalities=[], data_list=''):
    for root, dirs, files in os.walk(base):
        for file in files:
            try:
                f_path = os.path.join(root, file)
                header = dcmread(f_path)
                pat_id = None  # Patient ID
                pat_id = header[(0x10, 0x20)].value  # Patient ID
                modality = None  # Modality
                modality = header[(0x08, 0x60)].value  # Modality
                if modality in modalities:
                    # Initialize Some Variables
                    Manufacturer = None  # Manufacturer
                    ManufacturersModelName = None  # Manufacturer's Model Name
                    SliceThickness = None  # Slice Thickness
                    Kvp = None  # KVP
                    XRayTubeCurrent = None  # X-Ray Tube Current
                    Exposure = None  # Exposure
                    ConvolutionKernel = None  # Convolution Kernel
                    PixelSpacingx = None  # Pixel Spacingx
                    PixelSpacingy = None  # Pixel Spacingy
                    stu_inst_UID = None  # Study Instance UID
                    ser_inst_UID = None  # Series Instance UID
                    StudyDescription = None  # StudyDescription
                    SeriesDescription = None  # SeriesDescription
                    SOPClassUID = None  # SOP Class UID
                    SOPInstanceUID = None  # SOP Instance UID
                    StructureSetLabel = None  # Structure Set Label
                    ReferencedCTUID = None # ReferencedCTUID

                    # _______Operations on file_____________________
                    try:

                        # Study Instance UID
                        try:
                            stu_inst_UID = header[(0x20, 0x0d)].value
                        except Exception as e:
                            print(e, "stu_inst_UID", "Patient: ", pat_id,
                                  "Modality: ", modality)
                        # Series Instance UID
                        try:
                            ser_inst_UID = header[(0x20, 0x0e)].value
                        except Exception as e:
                            print(e, "ser_inst_UID", "Patient: ", pat_id,
                                  "Modality: ", modality)
                        # StudyDescription
                        try:
                            StudyDescription = header[(0x08, 0x1030)].value
                        except Exception as e:
                            print(e, "StudyDescription", "Patient: ", pat_id,
                                  "Modality: ", modality)
                        # SeriesDescription
                        try:
                            SeriesDescription = header[(0x08, 0x103e)].value
                        except Exception as e:
                            print(e, "SeriesDescription", "Patient: ", pat_id,
                                  "Modality: ", modality)
                        # SOP Class UID
                        try:
                            SOPClassUID = header[(0x08, 0x16)].value
                        except Exception as e:
                            print(e, "SOPClassUID", "Patient: ", pat_id,
                                  "Modality: ", modality)
                        # SOP Instance UID
                        try:
                            SOPInstanceUID = header[(0x08, 0x18)].value
                        except Exception as e:
                            print(e, "SOPInstanceUID", "Patient: ", pat_id,
                                  "Modality: ", modality)

                        if modality == "RTSTRUCT":
                            # Structure Set Label
                            try:
                                StructureSetLabel = header[(
                                    0x3006, 0x02)].value
                            except Exception as e:
                                print(e, 'StructureSetLabel', 'Patient: ',
                                      pat_id, ' Modality: ', modality)

                        if modality == "CT":
                            # Manufacturer
                            try:
                                Manufacturer = header[(0x08, 0x70)].value
                            except Exception as e:
                                print(e, "Manufacturer", "Patient: ",
                                      pat_id, "Modality: ", modality)
                            # Manufacturer's Model Name
                            try:
                                ManufacturersModelName = header[(
                                    0x08, 0x1090)].value
                            except Exception as e:
                                print(e, "Manufacturer's Model Name",
                                      "Patient: ", pat_id, "Modality: ", modality)
                            # Slice Thickness
                            try:
                                SliceThickness = header[(0x18, 0x50)].value
                            except Exception as e:
                                print(e, "Slice Thickness", "Patient: ",
                                      pat_id, "Modality: ", modality)
                            # KVP
                            try:
                                Kvp = header[(0x18, 0x60)].value
                            except Exception as e:
                                print(e, "KVP", "Patient: ", pat_id,
                                      "Modality: ", modality)
                            # X-Ray Tube Current
                            try:
                                XRayTubeCurrent = header[(0x18, 0x1151)].value
                            except Exception as e:
                                print(e, ' X-Ray Tube Current', 'Patient: ',
                                      pat_id, ' Modality: ', modality)
                            # Exposure
                            try:
                                Exposure = header[(0x18, 0x1152)].value
                            except Exception as e:
                                print(e, ' Exposure', 'Patient: ',
                                      pat_id, ' Modality: ', modality)
                            # Convolution Kernel
                            try:
                                ConvolutionKernel = header[(
                                    0x18, 0x1210)].value
                            except Exception as e:
                                print(e, ' Convolution Kernel', 'Patient: ',
                                      pat_id, ' Modality: ', modality)
                            # Pixel Spacing
                            try:
                                PixelSpacingx = header[(0x28, 0x30)].value[0]
                                PixelSpacingy = header[(0x28, 0x30)].value[1]
                            except Exception as e:
                                print(e, ' Pixel Spacing', 'Patient: ',
                                      pat_id, ' Modality: ', modality)

                        if modality == "RTSTRUCT":
                            ReferencedCTUID = getReferencedCTUID(header) 

                        data_list.append([pat_id,
                                        stu_inst_UID,
                                        StudyDescription,
                                        ser_inst_UID,
                                        SeriesDescription,
                                        ReferencedCTUID,
                                        StructureSetLabel,
                                        modality,
                                        SOPClassUID,
                                        SOPInstanceUID,
                                        Manufacturer,
                                        ManufacturersModelName,
                                        SliceThickness,
                                        Kvp,
                                        XRayTubeCurrent,
                                        Exposure,
                                        ConvolutionKernel,
                                        PixelSpacingx,
                                        PixelSpacingy,
                                        f_path])
                        print('Data Entered.: ', pat_id)
                    except Exception as e:
                        print("File level Error")
                        print(e)
                        print('Patient: ', pat_id, ' Modality: ', modality)
                        continue
            except Exception as e:
                print(e, "\nfile is not Readable")
                continue

    return data_list

def makeDataBse(source='', modalities=[]):
    data_list = [[
        'Patient_ID',
        'Study_Instance_UID',
        'Study_Description',
        'Series_Instance_UID',
        'Series_Description',
        'ReferencedCTUID'
        'Structure_Set_Label',
        'Modality',
        'SOP_Class_UID',
        'SOP_Instance_UID',
        'Manufacturer',
        'Manufacturers_Model_Name',
        'Slice_Thickness',
        'KVP',
        'X-Ray_Tube_Current',
        'Exposure',
        'Convolution_Kernel',
        'Pixel_Spacing_x',
        'Pixel_Spacing_y',
        'File_Path']]

    # makeFolders(destination)
    if type(source) == str:
        data_list = theFunction(source, modalities, data_list)
    elif type(source) == list:
        for path in source:
            print(f"Started: {path}")
            data_list = theFunction(path, modalities, data_list)
            print(f"Finished: {path}")
    data  = pd.DataFrame(data_list,index=None, columns=data_list[0])
    data.to_csv("outfile.csv",index=None)
    print('CSV_file --done!')
    data.to_pickle("outfile.pkl")
    print('pkl_file --done!')


# =========================================FIXME Block=======================================

# source = []

source = '$path'
modalities = ["RTSTRUCT", "CT"]

# ====================================Need Not Modify Block==================================
if __name__ == '__main__':
    tic()
    makeDataBse(
        source=source,
        modalities=modalities
    )
    print('Done!')
    tac()
