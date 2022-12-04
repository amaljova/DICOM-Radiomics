"""
@ Date: 15 November 2022
@ DateMod: 02 December 2022
@ author: AMAL JOSEPH VARGHESE
@ email: amaljova@gmail.com
@ github: https://github.com/amaljova



This Script contains snippets from https://github.com/zhenweishi/O-RAW (The Class DicomDatabase) and https://github.com/zhenweishi/Py-rex
"""

import radiomics
from radiomics import featureextractor
import pydicom
import os
import numpy as np
from skimage import draw
import SimpleITK as sitk
import re
import glob
import pandas as pd


# -------------------------------------------------------------------------------------------
# DICOM DATABASE CLASS
# from https://github.com/zhenweishi/O-RAW
class DicomDatabase:
    def __init__(self):
        self.patient = dict()

#============================    MODIFIED !    ===================================   
# Modified by Amal (https://github.com/amaljova)
    """
    Use this when using a csv file with all needed paths
    """

    def parseFiles(self, files):
        for file_path in files:
            if(file_path.endswith(".dcm") or file_path.endswith(".DCM")):
                dcmHeader = pydicom.dcmread(file_path)
                patientId = dcmHeader[0x10,0x20].value
                patient = self.getOrCreatePatient(patientId)
                patient.addFile(file_path, dcmHeader)


#==================================================================================

    def parseFolder(self, folderPath):
        for root, subdirs, files in os.walk(folderPath):
            for filename in files:
                file_path = os.path.join(root, filename)
                if(file_path.endswith(".dcm") or file_path.endswith(".DCM")):
                    dcmHeader = pydicom.dcmread(file_path)
                    patientId = dcmHeader[0x10,0x20].value
                    patient = self.getOrCreatePatient(patientId)
                    patient.addFile(file_path, dcmHeader)
        print("done")

#==================================================================================

    def getOrCreatePatient(self, patientId):
        if not (patientId in self.patient):
            self.patient[patientId] = Patient()
        return self.patient[patientId]
    
    def countPatients(self):
        return len(self.patient)
    
    def getPatient(self, patientId):
        return self.patient[patientId]
    
    def getPatientIds(self):
        return self.patient.keys()
    
    def doesPatientExist(self, patientId):
        # return self.patient.has_key(patientId)
        return patientId in self.patient
class Patient:
    def __init__(self):
        self.ct = dict()
        self.rtstruct = dict()

    def addFile(self, filePath, dcmHeader):
        modality = dcmHeader[0x8,0x60].value
        sopInstanceUid = dcmHeader[0x8,0x18].value
        seriesInstanceUid = dcmHeader[0x20,0xe].value
        if(modality == "CT") or (modality == "PT") or (modality == "MR"):
            if not (seriesInstanceUid in self.ct):
                self.ct[seriesInstanceUid] = CT()
            myCT = self.ct[seriesInstanceUid]
            myCT.addCtSlice(filePath)
        if(modality == "RTSTRUCT"):
            struct = RTStruct(filePath)
            self.rtstruct[sopInstanceUid] = struct
    
    def countCTScans(self):
        return len(self.ct)
        
    def countRTStructs(self):
        return len(self.rtstruct)
    
    def getCTScan(self, seriesInstanceUid):
        if seriesInstanceUid is not None:
            if self.doesCTExist(seriesInstanceUid):
                return self.ct[seriesInstanceUid]
        return None
    def getRTStruct(self, sopInstanceUid):
        return self.rtstruct[sopInstanceUid]
    
    def getCTScans(self):
        return self.ct.keys()
    def getRTStructs(self):
        return self.rtstruct.keys()
    
    def doesCTExist(self, seriesInstanceUid):
        # return self.ct.has_key(seriesInstanceUid)
        return seriesInstanceUid in self.ct
    def doesRTStructExist(self, sopInstanceUid):
        # return self.rtstruct.has_key(sopInstanceUid)
        return sopInstanceUid in self.rtstruct
    
    def getCTForRTStruct(self, rtStruct):
        if rtStruct.getReferencedCTUID() is not None:
            return self.getCTScan(rtStruct.getReferencedCTUID())
        else:
            return None

class CT:
    def __init__(self):
        self.filePath = list()
    def addCtSlice(self, filePath):
        self.filePath.append(filePath)
    def getSlices(self):
        return self.filePath
    def getSliceCount(self):
        return len(self.filePath)
    def getSliceHeader(self, index):
        return pydicom.dcmread(self.filePath[index])

class RTStruct:
    def __init__(self, filePath):
        self.filePath = filePath
    def getHeader(self):
        return pydicom.dcmread(self.filePath)
    def getReferencedCTUID(self):
        dcmHeader = self.getHeader()
        if len(list(dcmHeader[0x3006,0x10])) > 0:
            refFrameOfRef = (dcmHeader[0x3006,0x10])[0]
            if len(list(refFrameOfRef[0x3006, 0x0012])) > 0:
                rtRefStudy = (refFrameOfRef[0x3006, 0x0012])[0]
                if len(list(rtRefStudy[0x3006,0x14])) > 0:
                    rtRefSerie = (rtRefStudy[0x3006,0x14])[0]
                    return rtRefSerie[0x20,0xe].value
        return None
    def getFileLocation(self):
        return self.filePath


# ----------------------------------------------------------------------------------------------------------------
# To get the image and mask array of ROI of interest, from dicom files.
# MODIFIED PYREX from GITHUB
# from https://github.com/zhenweishi/Py-rex


# module PyrexReader:
def match_ROIid(rtstruct_path,ROI_name): # Match ROI id in RTSTURCT to a given ROI name in the parameter file
    mask_vol = Read_RTSTRUCT(rtstruct_path)
    M= mask_vol[0]
    for i in range(len(M.StructureSetROISequence)):
        if str(ROI_name)==M.StructureSetROISequence[i].ROIName:
            ROI_number = M.StructureSetROISequence[i].ROINumber
            break
    for ROI_id in range(len(M.StructureSetROISequence)):
        if ROI_number == M.ROIContourSequence[ROI_id].ReferencedROINumber:
            break
    return ROI_id

def Read_scan(myCT): # Read scans under the specified path 
    scan = [myCT.getSliceHeader(i) for i in range(myCT.getSliceCount())]
    try:
        scan.sort(key = lambda x: int(x.ImagePositionPatient[2])) # sort slices based on Z coordinate
    except:
        print('AttributeError: Cannot read scans')
    return scan

def Read_RTSTRUCT(myStruct): # Read RTSTRUCT under the specified path
    try:
        rt = [myStruct.getHeader()]
    except:
        print('AttributeError: Cannot read RTSTRUCT')
    return rt

def poly2mask(vertex_row_coords, vertex_col_coords, shape): # Mask interpolation
    fill_row_coords, fill_col_coords = draw.polygon(vertex_row_coords, vertex_col_coords, shape)
    mask = np.zeros(shape, dtype=bool)
    mask[fill_row_coords, fill_col_coords] = True
    return mask

def get_pixels_hu(scans): # convert to Hounsfield Unit (HU) by multiplying rescale slope and adding intercept
    image = np.stack([s.pixel_array for s in scans])
    image = image.astype(np.int16) #convert to int16
    # the code below checks if the image has slope and intercept
    # since MRI images often do not provide these
    try:
        intercept = scans[0].RescaleIntercept
        slope = scans[0].RescaleSlope
    except AttributeError:
        pass
    else:
        if slope != 1:
            image = slope * image.astype(np.float64)
            image = image.astype(np.int16)        
        image += np.int16(intercept)    
    return np.array(image, dtype=np.int16)

def Img_Bimask(img_path,rtstruct_path,ROI_name): # generating image array and binary mask
    print('Generating binary mask based on ROI: %s ......' % ROI_name)
    img_vol = Read_scan(img_path)
    mask_vol=Read_RTSTRUCT(rtstruct_path)
    IM=img_vol[0] # Slices usually have the same basic information including slice size, patient position, etc.
    IM_P=get_pixels_hu(img_vol)
    M=mask_vol[0]
    num_slice=len(img_vol)
    mask=np.zeros([num_slice, IM.Rows, IM.Columns],dtype=np.uint8)
    xres=np.array(IM.PixelSpacing[0])
    yres=np.array(IM.PixelSpacing[1])
    slice_thickness=np.abs(img_vol[1].ImagePositionPatient[2]-img_vol[0].ImagePositionPatient[2])

    ROI_id = match_ROIid(rtstruct_path,ROI_name)
    #Check DICOM file Modality
    # Modified by Amal (https://github.com/amaljova)
    if IM.Modality == 'CT' or IM.Modality == 'PT':
        for k in range(len(M.ROIContourSequence[ROI_id].ContourSequence)):
            Cpostion_rt = M.ROIContourSequence[ROI_id].ContourSequence[k].ContourData[2]
            for i in range(num_slice):
                if np.int64(Cpostion_rt) == np.int64(img_vol[i].ImagePositionPatient[2]): # match the binary mask and the corresponding slice
                    sliceOK = i
                    break
            x=[]
            y=[]
            z=[]
            m=M.ROIContourSequence[ROI_id].ContourSequence[k].ContourData
            for i in range(0,len(m),3):
                x.append(m[i+1])
                y.append(m[i+0])
                z.append(m[i+2])
            x=np.array(x)
            y=np.array(y)
            z=np.array(z)
            x-= IM.ImagePositionPatient[1]
            y-= IM.ImagePositionPatient[0]
            z-= IM.ImagePositionPatient[2]
            pts = np.zeros([len(x),3])
            pts[:,0] = x
            pts[:,1] = y
            pts[:,2] = z
            a=0
            b=1
            p1 = xres
            p2 = yres
            m=np.zeros([2,2])
            m[0,0]=img_vol[sliceOK].ImageOrientationPatient[a]*p1
            m[0,1]=img_vol[sliceOK].ImageOrientationPatient[a+3]*p2
            m[1,0]=img_vol[sliceOK].ImageOrientationPatient[b]*p1
            m[1,1]=img_vol[sliceOK].ImageOrientationPatient[b+3]*p2
          # Transform points from reference frame to image coordinates           
            m_inv=np.linalg.inv(m)
            pts = (np.matmul((m_inv),(pts[:,[a,b]]).T)).T
            mask[sliceOK,:,:] = np.logical_or(mask[sliceOK,:,:],poly2mask(pts[:,0],pts[:,1],[IM_P.shape[1],IM_P.shape[2]]))
    elif IM.Modality == 'MR':
        slice_0 = img_vol[0]
        slice_n = img_vol[-1]

		 # the screen coordinates, including the slice number can then be computed 
		  # using the inverse of this matrix
        transform_matrix = np.r_[slice_0.ImageOrientationPatient[3:], 0, slice_0.ImageOrientationPatient[:3], 0, 0, 0, 0, 0, 1, 1, 1, 1].reshape(4, 4).T # yeah that's ugly but I didn't have enough time to make anything nicer
        T_0 = np.array(slice_0.ImagePositionPatient)
        T_n = np.array(slice_n.ImagePositionPatient)
        col_2 = (T_0 - T_n) / (1 - len(img_vol))
        pix_s = slice_0.PixelSpacing
        transform_matrix[:, -1] = np.r_[T_0, 1] 
        transform_matrix[:, 2] = np.r_[col_2, 0] 
        transform_matrix[:, 0] *= pix_s[1]
        transform_matrix[:, 1] *= pix_s[0]
        
        transform_matrix = np.linalg.inv(transform_matrix)
        for s in M.ROIContourSequence[ROI_id].ContourSequence:
            Cpostion_rt = np.r_[s.ContourData[:3], 1]
            roi_slice_nb = int(transform_matrix.dot(Cpostion_rt)[2])
            for i in range(num_slice):
                print(roi_slice_nb, i)
                if roi_slice_nb == i:
                    sliceOK = i
                    break
            x=[]
            y=[]
            z=[]
            m=s.ContourData
            for i in range(0,len(m),3):
                x.append(m[i+1])
                y.append(m[i+0])
                z.append(m[i+2])
            x=np.array(x)
            y=np.array(y)
            z=np.array(z)
            x-= IM.ImagePositionPatient[1]
            y-= IM.ImagePositionPatient[0]
            z-= IM.ImagePositionPatient[2]
            pts = np.zeros([len(x),3])
            pts[:,0] = x
            pts[:,1] = y
            pts[:,2] = z
            a=0
            b=1
            p1 = xres
            p2 = yres
            m=np.zeros([2,2])
            m[0,0]=img_vol[sliceOK].ImageOrientationPatient[a]*p1
            m[0,1]=img_vol[sliceOK].ImageOrientationPatient[a+3]*p2
            m[1,0]=img_vol[sliceOK].ImageOrientationPatient[b]*p1
            m[1,1]=img_vol[sliceOK].ImageOrientationPatient[b+3]*p2
            # Transform points from reference frame to image coordinates
            m_inv=np.linalg.inv(m)
            pts = (np.matmul((m_inv),(pts[:,[a,b]]).T)).T
            mask[sliceOK,:,:] = np.logical_or(mask[sliceOK,:,:],poly2mask(pts[:,0],pts[:,1],[IM_P.shape[1],IM_P.shape[2]]))
    
    # The pixel intensity values are normalized to range [0 255] using linear translation   
    IM_P=IM_P.astype(np.float32)
    # IM_P = (IM_P-np.min(IM_P))*255/(np.max(IM_P)-np.min(IM_P))  

    Img=sitk.GetImageFromArray(IM_P) # convert image_array to image
    Mask=sitk.GetImageFromArray(mask)
    # try:
    #     origin = IM.GetOrigin()
    # except:
    #     origin = (0.0, 0.0, 0.0)

    # Set voxel spacing [[pixel spacing_x, pixel spacing_y, slice thickness]
    #slice_thickness = IM.SliceThickness
    Img.SetSpacing([np.float64(xres),np.float64(yres),np.float64(slice_thickness)])
    Mask.SetSpacing([np.float64(xres),np.float64(yres),np.float64(slice_thickness)])
    # sitk.WriteImage(Img,'./data/test_image.nrrd',True)
    # sitk.WriteImage(Mask,'./data/test_label.nrrd',True)
    return Img, Mask




# ---------------------------------------------------------------------------------------

# Extracting Radiomic Features
def extractFeatures(dicomDb, roi, parameter, out_name):
    """
    This Fuction extracts radiomic features automatically.
    """
    # create an empty dataframe to collect results
    result = pd.DataFrame()
    # Get Patients
    for ptid in dicomDb.getPatientIds():
        print("staring with Patient %s" % (ptid))
        myPatient = dicomDb.getPatient(ptid)
        # get RTSTRUCTs
        for myStructUID in myPatient.getRTStructs():
            print("Starting with RTStruct %s" % myStructUID)
            myStruct = myPatient.getRTStruct(myStructUID)
            # Get the related CT
            myCT = myPatient.getCTForRTStruct(myStruct)
            if myCT == None:
                continue
            # Collect all ROI names
            roi_list= [i[(0x3006, 0x26)].value for i in myStruct.getHeader()[(0x3006, 0x20)].value]
            # Take each ROI name
            for roi_name in roi_list:
                # Check if it is GTV or not.
                if re.search(roi,roi_name):
                    # Convert DICOM images into Raw data and Mask
                    Img, Mask = Img_Bimask(myCT,myStruct,roi_name)
                    # Initialize feature extractor object using parameter file
                    extractor = featureextractor.RadiomicsFeatureExtractor(parameter)
                    # Extract and collect features, add patient ID, structure set UID and Contour name to it.
                    featureVector = pd.Series({"patient": ptid, "structUID": myStructUID, "contour": roi_name})
                    featureVector = featureVector.append(pd.Series(extractor.execute(Img, Mask)))
                    # Collect the results in a dataframe
                    result = result.append(featureVector.to_frame().T)
                    # simply set these variables to empty list 
                    Img, Mask = [],[]
    # Save Results in CSV file.
    result.to_csv(out_name,index=False)
    print(f"Created {out_name}")

def makeFolders(dir_name):
    '''To create a directory/directories if it/they are basent'''
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)


def execute(data_dir, roi, parameter_list,outDir,data_name):
    # create dicom databese from the given folder
    dicomDb = DicomDatabase()
    print("Instantiated Dicomdatabase & Started parsing folders...")
    # walk over all files in folder, and index in the database
    dicomDb.parseFolder(data_dir)
    print("Parsing is done!")
    for param in parameter_list:
        print(f"Using {param}")
        out_name = os.path.join(outDir, f'{data_name}_{((param.split("/")[-1]).split(".")[0])}_Out.csv')
        extractFeatures(dicomDb, roi, param, out_name)


if __name__ == "__main__":

    # EDIT HERE

    data_dir = "path_to_dataset"

    # ROI of interest (regular expresion)
    roi = '[Gg][Tt][Vv]'

    # Edit list to include more parameter files (if interested.)
    parameter_list = [
        "Params/Pyradiomics_Params.yaml"
    ]

    outDir = "results_folder" # where you want to save csv files with extracted features
    data_name = "dataset_name" # name of the dataset, just to name output file

    makeFolders(outDir)
    execute(data_dir, roi, parameter_list, outDir, data_name)


