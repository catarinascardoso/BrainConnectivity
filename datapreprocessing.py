import nibabel as nib
import numpy as np
import os
from glob import glob
import fnmatch
import matplotlib.pyplot as plt

DATAPATH = '../data/hcp/'
MASKPATH = '../data/mask/'

def searchData(dataPath):
    scansPath = glob(DATAPATH+'*.gz', recursive=True)
    return scansPath

def readScans(scansPath):
    print("Number of scans read:", len(scansPath))
    print()
    for file in scansPath:
        path, filename = os.path.split(file)
        print(filename)

def loadMask(maskPath):
    maskData = nib.load(MASKPATH+'brainmask_fs.2.nii.gz').get_data()
    mask = np.asarray(maskData)
    return mask

def checkMaskShape(mask):
    print("Mask's shape: ", mask.shape)

def visualizeMask(mask):
    mask=np.rot90(mask)
    f, ax = plt.subplots(5,9,figsize=(8,7))
    for i in range(0,45):
        ax[i//9,i%9].imshow(mask[:,:,i*2], cmap='gray')
        ax[i//9,i%9].axis('off')
        ax[i//9,i%9].set_title('Slice {}'.format(i*2), fontsize=10)
    f.suptitle('Mascara')
    plt.show()
    
def coordsArray(mask):
    coords=[] 
    for x in range(mask.shape[0]):
        for y in range(mask.shape[1]):
            for z in range(mask.shape[2]):
                if mask[x][y][z]==1:
                    coords.append((x,y,z))
    return coords

def checkCoords(coords, c):
    print('Number of coordinates within the brain: %d' % (len(coords)))
    print('Coordinate # %d within the brain mask: [%d, %d, %d] '% (c,coords[c][0],coords[c][1],coords[c][2]))
    
def loadScan(scansPath, scanID):
    scanData = nib.load(scansPath[scanID]).get_data()
    scan = np.asarray(scanData)
    return scan
    
def checkScanShape(scan):
    scanShape = scan.shape
    print("Scan's shape: ", scan.shape)
    return scanShape

def visualizeBrain(scan, t):
    scan = np.rot90(scan)
    f,ax = plt.subplots(5,9,figsize=(8,7))
    for i in range(0,45):
        ax[i//9,i%9].imshow(scan[:,:,i*2,t], cmap='gray')
        ax[i//9,i%9].axis('off')
        ax[i//9,i%9].set_title('Slice {}'.format(i*2), fontsize=10)
    f.suptitle('Connectome instant = %d' % (t))
    plt.show()
    
def createNumpyArray(coords, scanShape):
    timeseriesNo = scanShape[3] 
    scanData = np.zeros((len(coords),timeseriesNo),dtype=np.float32)
    print("scanData com shape: ", scanData.shape)
    return scanData, timeseriesNo

def saveNumpyDataPerSubject(scansPath, coords, scanData):
    scansNo = len(scansPath)
    scanID = 0
    for file in scansPath:
        data = nib.load(file).get_data()
        for c in range(len(coords)): 
            x, y, z = coords[c] 
            scanData[c] = data[x][y][z]

        fileID = scansPath[scanID]
        ID = fileID[17:30]
        print("Case " +str(scanID+1)+" built out of "+str(scansNo)+". Shape: " + str(scanData.shape))
        print()

        np.save('../npdata/'+str(ID)+'.npy', scanData)
        scanID += 1
    return scanData
    
def createNumpyArrayForLineGraphs(scansPath, timeseriesNo, coords):
    twoScans = scansPath[181:183]
    graphData = np.zeros((timeseriesNo,len(coords)),dtype=np.float32)
    print("graphData com shape: ", graphData.shape)
    return twoScans, graphData
    
def plotLineGraphs(twoScans, graphData, timeseriesNo, coords, lt, ht, lc, hc):
    scanID = 0
    for file in twoScans:
        data = nib.load(file).get_data() 
        for t in range(timeseriesNo): 
            for c in range(len(coords)):
                x, y, z = coords[c]
                graphData[t][c] = data[x][y][z][t]
        
        scanID += 1
    return graphData
        
def visualizeTimeSeries(graphData, lt, ht, lc, hc, i):
    plt.title("Subject's Time Series", fontsize = 12)
    plt.plot(graphData[lt:ht,lc:hc:i], linewidth = 1)
    plt.show()
    print()
    print("Shape: ", graphData[lt:ht,lc:hc:i].shape)
    
def checkMinMax(graphData):
    m = np.amin(graphData)
    M = np.amax(graphData)
    print("Minimo: ", m)
    print("Maximo: ", M)
    
def loadNumPyArrays(normalizedDataPath):
    allfiles = os.listdir(path)
    return allfiles

def MinMaxNormalization(allfiles):
    for file in allfiles:
        data = np.load(dataPath+"/"+file)
        dataMin = data.min()
        dataMax = data.max()
        print("Max: %f  Min: %f "%(dataMax,dataMin)
        normdata = (data - dataMin)/(dataMax-dataMin)
        print("NormMax: %f  NormMin: %f "%(normdata.max(),normdata.min()))
        np.save("../npdatanorm/"+file, data)