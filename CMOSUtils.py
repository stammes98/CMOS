# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 09:40:44 2019

@author: stamm
"""

from pyueye import ueye
import sys
import numpy as np
import datetime
import astropy.io.fits as pyfits
import time
import math
import glob
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import curve_fit
import gc
import os
import psutil
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)

badShapes = [4112, 2048, 528384, 1052672, 8392704, 4097, 4098, 4100, 4104, 4120, 4128, 5120, 6148, 4256, 4608, 4736, 5152, 13312, 20480, 20992, 36864, 135172, 
             135174, 25169920, 16781312, 8394752, 135176, 405505, 530432, 544768, 659456, 1054912, 1060864, 2101248, 2103488, 2109440, 2109824, 4198528, 4200448, 
             8394944, 8400896, 4102, 4108, 4225, 5248, 6145, 6156, 6160, 6656, 6672, 12289, 12291, 12292, 12294, 12296, 12321, 12448, 12673, 12704, 13568, 20608, 
             22656, 22720, 36992, 37888, 45056, 45440, 135169, 135680, 136192, 136224, 202756, 202760, 203264, 398368, 405508, 405512, 405536, 406528, 661504, 
             1054720, 1085440, 2103296, 2109568, 3152000, 4198400, 4200640, 4206592, 4206720, 6295552, 8392832, 8401280, 12587008, 16781440, 16783424, 16789888,
             16912384, 16914432, 17305600, 176128, 4099, 4240, 4624, 6152, 6288, 46080, 135704, 137217, 143392, 202758, 202768, 405506, 405516, 405520, 1052800, 
             3158400, 12587136, 12595200, 16789504, 135180, 398336, 2103424, 1183744, 1191936, 1454080, 4129, 4481, 6800, 46208, 135696, 137728, 143366, 202753,
             3149824, 16920576, 25170048, 12290, 4131, 4632, 4800, 5280, 6146, 6147, 6150, 6168, 6273, 6352, 6664, 6848, 12304, 12320, 12417, 12545, 12705, 
             13344, 13472, 13696, 21120, 22528, 38016, 46464, 72192, 135170, 135171, 135192, 135200, 137220, 137224, 137240, 143362, 143376, 143384, 144384, 144416,
             151552, 152064, 153600, 154112, 167936, 168960, 177152, 202754, 202776, 219136, 219648, 282624, 397315, 397316, 397320, 405507, 439296, 528512, 530624, 
             544896, 546816, 675840, 727040, 743424, 1060992, 1061248, 1185792, 2101376, 2101440, 3182592, 3152064, 4200512, 6297600, 6297792, 6303744, 6304128,
             8401024, 8460288, 12589056, 12589184, 16783360, 16783552, 16789632, 17307648, 17321984, 25178112, 45184]

forceHot = [[261, 824], [249, 295], [1080, 1713], [1049, 605], [341, 290], [396, 830], [601, 820], [1051, 1852], [584, 909],
            [243, 688], [441, 1323], [731, 1371], [513, 1309], [9, 1132], [961, 984], [651, 310], [6, 850]]

def aduToeV(adu, m, b):
    return ((m * adu ) + b) * 1000

def smoothing(data, pts):
    box = np.ones(pts)/pts
    smooth = np.convolve(data, box, mode='same')
    return smooth

def mean_2d(img) :
  """
  Find the mean of a 2d image
  Paramaters:
      img - the image
      
  """
  #img is a 2-d array, need to change to 1-d to do calculations
  
  nx, ny = img.shape # find the size of the array
  img1 = np.reshape(img, nx*ny) # change the shape to be 1d
  return np.mean(img1)

# find the median of a 2d image
def median_2d(img) :
  """
  Find the median of a 2d image
  Paramaters:
    img - the image
      
  """
  # img is a 2-d array, need to change to 1-d to do calculations
  nx, ny = img.shape # find the size of the array
  img1 = np.reshape(img, nx*ny) # change the shape to be 1d
  return np.median(img1)


# find the standard deviation of a 2d image
def std_2d(img) :
  """
  Find the std of a 2d image
  Paramaters:
      img - the image
      
  """
  # img is a 2-d array, need to change to 1-d to do calculations
  nx, ny = img.shape # find the size of the array
  img1 = np.reshape(img, nx*ny) # change the shape to be 1d
  return np.std(img1)


def imStats(img):
    """
    Prints the mean, median and std of a 2d image
    Paramaters:
        img - the image
    """
    print("Image statistics:")
    print("Mean: " + str(mean_2d(img)))
    print("Median: " + str(median_2d(img)))
    print("Std: " + str(std_2d(img)))
    print()

def getImStats(img):
    """
    Returns the mean, median and std of a 2d image as a tuple
    Paramaters:
        img - the image
    """
    mean = mean_2d(img)
    med = median_2d(img)
    std = std_2d(img)
    return [mean, med, std]

def avgBias(biasPath):
    """
    Averages out a set of bias frames and returns it
    Paramaters:
        biasPath - path to the bias images, include file name but not ?.FTS
    """
    if (os.path.exists(biasPath + "avgBias.FTS")):
        print("Average bias already compiled")
        f = pyfits.open(biasPath + "avgBias.FTS")
        dat = f[0].data
        f.close()
        return dat
    
    files = glob.glob(biasPath + "?.FTS") + glob.glob(biasPath + "??.FTS") + glob.glob(biasPath + "???.FTS")
    shaping = pyfits.open(files[0])
    index = 0
    if(len(shaping) == 1):
        index = 0
    else:
        index = 1
    shape = shaping[index].data
    shaping.close()
    medBias = np.empty_like(shape)
    totBias = np.empty_like(medBias)
    for i in range(0, len(files)):
        data = pyfits.open(files[i])
        imgData = data[index].data
        totBias += imgData
        data.close()
    medBias = totBias / len(files)
    medBias = medBias.astype(int)
    print("Bias mean: " + str(mean_2d(medBias)))
    print("Bias frame compiled from " + str(len(files)) + " frames")
    pyfits.writeto(biasPath + 'avgBias.FTS', medBias, None, overwrite=True)
    return medBias
    
def avgDark(darkPath, verbose=True):
    if (os.path.exists(darkPath + "avgDark.FTS")):
        print("Average dark already compiled.")
        f = pyfits.open(darkPath + "avgDark.FTS")
        dat = f[0].data
        f.close()
        return dat
        
    """
    Averages out a set of dark frames
    Paramaters:
        darkPath - path to the bias images, include file name but not ?.FTS
    """
    files = glob.glob(darkPath + "?.FTS") + glob.glob(darkPath + "??.FTS") + glob.glob(darkPath + "???.FTS")
    shaping = pyfits.open(files[0])
    index = 0
    if(len(shaping) == 1):
        index = 0
    else:
        index = 1
    shape = shaping[index].data
    shaping.close()
    medDark = np.empty_like(shape)
    totDark = np.empty_like(medDark)
    for i in range(0, len(files)):
        data = pyfits.open(files[i])
        imgData = data[index].data
        totDark += imgData
        data.close()
    medDark = totDark / len(files)
    medDark = medDark.astype(int)
    if(verbose):
        print("Dark mean: " + str(mean_2d(medDark)))
        print("Dark frame compiled from " + str(len(files)) + " frames")
    pyfits.writeto(darkPath + 'avgDark.FTS', medDark, None, overwrite=True)

    return medDark
    
def getFilesForBrute(filePath, bar=True):
    fs = glob.glob(filePath + "?.FTS") + glob.glob(filePath + "??.FTS") + glob.glob(filePath + "???.FTS") + glob.glob(filePath + "????.FTS") +  glob.glob(filePath + "?????.FTS")
    ps = []
    for i in range(0, len(fs)):
        prog = (i/len(fs)) * 100
        if(bar):
            printProgressBar(prog, 100, prefix = 'Loading images:', suffix='Complete', length=50)
        dat = pyfits.open(fs[i])
        if(len(dat) == 1):
            ps.append(dat[0].data)
        else:        
            ps.append(dat[1].data)
        dat.close()
    print()
    return ps
    
def fileAnalysisBrute(phots, biasPath, darkPath, thres, thresp, savePath, x=2, y=2, autoThres=False):
    t = thres
    tp = thresp
    """
    Performs the analysis of the photos at filePath
    Returns a tuple of event intensity, x location and y location lists
    Paramaters:
        filePath - list of the photos (must be photos taken with C++ code, not Python)
        biasPath - path to the bias frames (don't include .FTS in path)
        darkPath - path to the dark frames (don't include .FTS in path)
        thres - threshold intensity for an x-ray to be detected
        thresp - threshold intensity for neighboring pixel to x-ray to be counted
        savePath - where to save the .CSV with the lists (don't include .CSV in path)
    """
    ep, ex, ey = [], [], []

    # set image size
    nx, ny = 1936, 1096
    # set integration region
    dx, dy = x, y
    biasFrame = avgBias(biasPath)
    darkFrame = avgDark(darkPath)
    print("Performing analysis of " + str(len(phots)) + " photos")
    
    if autoThres:
        tp = 2.5 * std_2d(biasFrame + darkFrame)
        t = 5 * std_2d(biasFrame + darkFrame)
    print("thres set to " + str(t))
    print("thresp set to " + str(tp))
    print("dx set to " + str(dx))
    print("dy set to " + str(dy))
    hotPixels = getHotPix(darkPath)
    #end len(files)
    for c in range(0, len(phots)):
        phot = phots[c] - (biasFrame + darkFrame)
        phot[phot < 0] = 0
        for ii in range(0, len(hotPixels)):
            phot[hotPixels[ii][0]][hotPixels[ii][1]] = 0
        #print(phot)
        nx, ny = phot.shape # find size of image
        f1 = np.reshape(phot, nx*ny) # change the shape to be 1d
        q = np.argsort(f1)[::-1]
        #print(f1[q[0]])
        #print(q[0])
        # go through pixels with large signal
        j = 0
        above_thres = True
        while above_thres :
            i = q[j]
            #print(j, q[j], f1[i])
            j += 1
            if f1[i] >= t : # pixel is above threshold
              #print("i: " + str(i) + " j: " + str(j))
              x = (i % 1936) + 1
              y = math.floor((i / 1936) + 1)
              yR = (i % 1936)
              xR = math.floor((i / 1936))
              #print("Checking pixel " + str(x) + ", " + str(y) + " with value " + str(phot[xR][yR]))
              # i = ny*x + y # index into 1d array
              # check if pixel is far enough from edge of chip
              if (xR > dx) and (xR < nx-dx-1) and (yR > dy) and (yR < ny-dy-1) :
                #p = sum(f1[i+di])
                p = 0
                for xi in range(xR-dx, xR+dx+1) :
                  for yi in range(yR-dy, yR+dy+1):
                      if phot[xi, yi] > tp : p += phot[xi,yi]
                      phot[xi, yi] = 0 # zero out pixels in this event
                #if (p > 106):
                #print(str(p) + " pulse at " + str(x) + ", " + str(y) + " on image " + str(files[c]) + " - c: " + str(c))                    
                #print("pulse at " + str(y) + ", " + str(x) + " - " + str(p))
                
                if (p > 0):
                    ep.append(p)
                    ex.append(x)
                    ey.append(y)

                #print (x, y, p)
                #f1[i+di] = 0 # zero out pixels in this event
            # if we reach a pixel below threshold than has not been zeroed, then we are done
            elif f1[i] > 0 :
              above_thres = False # stop looping
        prog = (c/len(phots)) * 100
        printProgressBar(prog, 100, prefix = 'Progress:', suffix='Complete', length=50)
          
    print()
    ep = np.array(ep)
    ex = np.array(ex)
    ey = np.array(ey)
    fi = open(savePath + '.csv', 'w')
    fi.write('ep, ex, ey\n')
    for cntr in range(0, len(ep)):
        fi.write(str(ep[cntr]) + ', ' + str(ex[cntr]) + ', ' + str(ey[cntr]) + '\n')
    fi.close()
    print(str(len(ep)) + " events written")
    return [ep, ex, ey]


def fileAnalysis(filePath, biasPath, darkPath, thres, thresp, savePath, x=2, y=2, autoThres=False):
    t = thres
    tp = thresp
    """
    Performs the analysis of the photos at filePath
    Returns a tuple of event intensity, x location and y location lists
    Paramaters:
        filePath - path to the photos (must be photos taken with C/C++ code, not Python)
        biasPath - path to the bias frames (don't include .FTS in path)
        darkPath - path to the dark frames (don't include .FTS in path)
        thres - threshold intensity for an x-ray to be detected
        thresp - threshold intensity for neighboring pixel to x-ray to be counted
        savePath - where to save the .CSV with the lists (don't include .CSV in path)
    """
    ep, ex, ey = [], [], []
    shp = []
    # set image size
    nx, ny = 1936, 1096
    # set integration region
    dx, dy = x, y
    biasFrame = avgBias(biasPath)
    darkFrame = avgDark(darkPath)
    darkFrame[darkFrame > t] = 0
    files = glob.glob(filePath + "?.FTS") + glob.glob(filePath + "??.FTS") + glob.glob(filePath + "???.FTS") + glob.glob(filePath + "????.FTS") + glob.glob(filePath + "?????.FTS")
    print("Performing analysis of " + str(len(files)) + " photos")
    if autoThres:
        tp = 2.5 * std_2d(biasFrame + darkFrame)
        t = 5 * std_2d(biasFrame + darkFrame)
    print("thres set to " + str(t))
    print("thresp set to " + str(tp))
    print("dx set to " + str(dx))
    print("dy set to " + str(dy))
    hotPixels = getHotPix(darkPath)
    #for ii in range(0, len(hotPixels)):
        #print("x: " + str(hotPixels[ii][1]) + " y: " + str(hotPixels[ii][0]))
    #end len(files)
    for c in range(0, len(files)):
        data = pyfits.open(files[c])
        if(len(data) == 1):
            imData = data[0].data
        else:
            imData = data[1].data
        data.close()
        phot = imData - (biasFrame + darkFrame)
        phot[phot < 0] = 0
        for ii in range(0, len(hotPixels)):
            phot[hotPixels[ii][0]][hotPixels[ii][1]] = 0
            
        #print(phot)
        nx, ny = phot.shape # find size of image
        f1 = np.reshape(phot, nx*ny) # change the shape to be 1d
        q = np.argsort(f1)[::-1]
        #print(f1[q[0]])
        #print(q[0])
        # go through pixels with large signal
        j = 0
        above_thres = True
        while above_thres :
            i = q[j]
            #print(j, q[j], f1[i])
            j += 1
            if f1[i] >= t : # pixel is above threshold
              #print("i: " + str(i) + " j: " + str(j))
              x = (i % 1936) + 1
              y = math.floor((i / 1936) + 1)
              yR = (i % 1936)
              xR = math.floor((i / 1936))
              #print("Checking pixel " + str(x) + ", " + str(y) + " with value " + str(phot[xR][yR]))
              # i = ny*x + y # index into 1d array
              # check if pixel is far enough from edge of chip
              if (xR > dx) and (xR < nx-dx-1) and (yR > dy) and (yR < ny-dy-1) :
                #p = sum(f1[i+di])
                p = 0
                shape = ""
                for xi in range(xR-dx, xR+dx+1) :
                  for yi in range(yR-dy, yR+dy+1):
                      if (phot[xi, yi] > tp):
                          shape = shape + "1"
                          p += phot[xi,yi]
                      else:
                          shape = shape + "0"
                      phot[xi, yi] = 0 # zero out pixels in this event
                #if (p > 106):
                    #print(str(p) + " pulse at " + str(x) + ", " + str(y) + " on image " + str(files[c]) + " - c: " + str(c))                    
                #print("pulse at " + str(y) + ", " + str(x) + " - " + str(p))
                
                if (p > 0):
                    ep.append(p)
                    ex.append(x)
                    ey.append(y)
                    shp.append(int(shape, 2))
                    #if (p > 3000):
                        #print("Potential cosmic ray of ADU " + str(p) + " on image " + str(files[c]))
                #print (x, y, p)
                #f1[i+di] = 0 # zero out pixels in this event
            # if we reach a pixel below threshold than has not been zeroed, then we are done
            elif f1[i] > 0 :
              above_thres = False # stop looping
        prog = (c/len(files)) * 100
        printProgressBar(prog, 100, prefix = 'Progress:', suffix='Complete', length=50)
          
    ep = np.array(ep)
    ex = np.array(ex)
    ey = np.array(ey)
    fi = open(savePath + '.csv', 'w')
    fi.write('ep, ex, ey\n')
    for cntr in range(0, len(ep)):
        fi.write(str(ep[cntr]) + ', ' + str(ex[cntr]) + ', ' + str(ey[cntr]) + ':' + str(shp[cntr]) + "\n")
    fi.close()
    print(str(len(ep)) + " events written")
    return [ep, ex, ey]

def getHotPix(darkPath, thr=9, rep=5):
    lpath = darkPath + "hotT" + str(thr) + "R" + str(rep) + ".csv"
    if(os.path.exists(lpath)):
        print("Hot pixel list already compiled")
        (x, y) = np.loadtxt(lpath, delimiter=',', unpack=True)
        pos = []
        for i in range(0, len(x)):
            pos.append(x[i], y[i])
        return pos
    
    ps = getFilesForBrute(darkPath, bar=False)
    print("Parsing " + str(len(ps)) + " files for hot pixels bounded by " + str(thr) + " and multiplicity " + str(rep))
    hots = []
    pos = []
    for phot in ps:
        nx, ny = phot.shape # find size of image
        f1 = np.reshape(phot, nx*ny) # change the shape to be 1d
        q = np.argsort(f1)[::-1]
        hot = True
        i = 0
        while (hot):
            if(f1[q[i]] > thr):
                i += 1
                hots.append(q[i])
            else:
                hot = False

    h, c = np.unique(hots, return_counts=True)
    for i in range(0, len(h)):
        y = (h[i] % 1936)
        x = math.floor((h[i] / 1936))
        if(c[i] > rep):
            pos.append([x,y])
            #print("(" + str(x) + ", " + str(y) + "): " + str(c[i]))
            
    for i in range(0, len(forceHot)):
        pos.append(forceHot[i])
    print(str(len(pos)) + " pixels found")
    print("Loss: %.2f" % ((len(pos)/(1096 * 1936))*100))
    f = open(darkPath + "hots.csv", 'w')
    for i in range(0, len(pos)):
        f.write(str(pos[i][0]) + ", " + str(pos[i][1]) + "\n")
    f.close()
    return pos



def spectraString(thres, thresp, dx, dy, ep):
    '''
    Converts a spectrum and its paramaters into a string
    See readSpectraString to undo this process 
    Paramaters:
        thres - the thres of the spectrum
        thresp - the thresp of the spectrum
        dx - x intergration size
        dy - y intergraiton size
        ep - list of the spectrum events
    '''
    lngStr = str(thres) + ", " + str(thresp) + ", " + str(dx) + ", " + str(dy) + ", "
    for i in range(0, len(ep)):
        if(i == 0):
            lngStr = lngStr + str(int(ep[i]))
        else:
            lngStr = lngStr + "-" + (str(int(ep[i])))

    return lngStr + "\n"

def readSpectraString(specStr):
    '''
    Given a string representing a spectrum, this will return a tuple of all the important values
    Returns [thres, thresp, dx, dy, ep]
    See spectraString to get aforementioned string
    Paramaters:
        specStr - the string representing the spectrum
    '''
    splitBits = specStr.split(',')
    thres = int(splitBits[0])
    thresp = int(splitBits[1])
    dx = int(splitBits[2])
    dy = int(splitBits[3])
    bigList = splitBits[4].split('-')
    ep = []
    for i in range(0, len(bigList)):
        vals = bigList[i].split(':')
        ep.append(int(vals[0]))
    return [thres, thresp, dx, dy, ep]
    
    
def addToMultiFile(ep, thres, thresp, dx, dy, fileName, checkForDuplicate=False):
    '''
    Finds the multispectrum CSV file located at fileName (include the .CSV in the path) and adds the specified spectrum to the file
    Paramaters:
        ep - the events of the spectrum
        thres - what the thres was set to
        thresp - what the thresp was set to
        dx - x intergration size
        dy - y intergration size
        fileName - name of the multispectrum CSV file
        checkForDuplicate - default false - compares the given spectrum with spectra within the file, if they're identical the spectrum won't be written
        
    '''
    if (os.path.exists(fileName)):
        f = open(fileName, 'r')
        oldFile = []
        for line in f:
            oldFile.append(line)
        #print(oldFile)
        f.close()
        if(checkForDuplicate):
            for i in range(0, len(oldFile)):
                if (spectraString(thres, thresp, dx, dy, ep) == oldFile[i]):
                    return
        dat = oldFile
        dat.append(spectraString(thres, thresp, dx, dy, ep))
    else:
        dat = spectraString(thres, thresp, dx, dy, ep)
    f = open(fileName, 'w')
    for i in range(0, len(dat)):
        f.write(dat[i])
    f.close()


def histogram(ep, saveImg=False, saveName="default.png"):
    """
    Display a histogram of the events
    Paramaters:
        ep - event intensity list
    """
    plt.ion()
    plt.figure('Event histogram')
    plt.clf()
    unis = np.unique(ep)
    b = int(max(unis) - min(unis))
    plt.hist(ep, bins=b, range=[min(unis),max(unis)])
    if(saveImg):
        plt.savefig(saveName)
    plt.ioff()
    plt.show()

def spectrum(ep, m, b, saveImg=False, saveName="default.png", tickSpace=400):
    """
    Produces a spectrum intensity graph
    Paramaters:
        ep - event intensity list
        m - slope of ADU to eV fit
        b - intercept of ADU to eV fit
    """
    unis, cnts = np.unique(ep, return_counts=True)
    unisList = dict(zip(unis, cnts))
    plt.ion()
    plt.figure('Spectrum', figsize=(18,10))
    plt.clf()
    lists = sorted(unisList.items()) # sorted by key, return a list of tuples
    x, y = zip(*lists) # unpack a list of pairs into two tuples
    x = np.asarray(x)
    xdat = aduToeV(x, m, b)
    y = np.asarray(y)
    ydat = smoothing(y, 1)
    #plt.rc('xtick', labelsize=8) 
    plt.xlabel("Energy (eV)")
    plt.ylabel("Intensity")
    plt.plot(xdat, ydat)
    #plt.axvline(537)
    plt.gca().xaxis.set_major_locator(MultipleLocator(200))
    plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%d'))
    plt.gca().xaxis.set_minor_locator(MultipleLocator(50))
    plt.xticks(np.arange(100, max(xdat), tickSpace))
    if(saveImg):
        plt.savefig(saveName)
    plt.ioff()
    plt.show()

def compSpectrum(ep1, ep2, m, b):
    """
    Produces two spectrum intensity graphs
    Paramaters:
        ep1 - event intensity list 1
        ep2 - event intensity list 2
        m - slope of ADU to eV fit
        b - intercept of ADU to eV fit
    """
    unis1, cnts1 = np.unique(ep1, return_counts=True)
    unis2, cnts2 = np.unique(ep2, return_counts=True)
    unisList1 = dict(zip(unis1, cnts1))
    unisList2 = dict(zip(unis2, cnts2))
    plt.ion()
    plt.figure('Spectrum')
    plt.clf()
    lists1 = sorted(unisList1.items()) # sorted by key, return a list of tuples
    lists2 = sorted(unisList2.items())
    x1, y1 = zip(*lists1) # unpack a list of pairs into two tuples
    x2, y2 = zip(*lists2)
    x1 = np.asarray(x1)
    x2 = np.asarray(x2)
    xdat1 = aduToeV(x1, m, b)
    xdat2 = aduToeV(x2, m, b)
    y1 = np.asarray(y1)
    y2 = np.asarray(y2)
    ydat1 = smoothing(y1, 2)
    ydat2 = smoothing(y2, 2)
    plt.xlabel("Energy (ev)")
    plt.ylabel("Intensity")
    plt.plot(xdat1, ydat1, label='Spectrum 1')
    plt.plot(xdat2, ydat2, label='Spectrum 2')
    plt.legend()
    plt.xticks(np.arange(min(xdat1), max(xdat1)+1, 100.0))
    plt.ioff()
    plt.show()


def spectrumFromFile(fileName, m, b):
    """
    Produces a spectrum intensity graph
    Paramaters:
        fileName - name of the file with events (don't include .CSV)
    """
    (ep, ex, ey) = np.loadtxt(fileName + '.csv', delimiter=',', unpack=True, skiprows=1)    
    spectrum(ep, m, b)


def findADUtoeV(ep, ev, plot):
    """
    Finds the relationship between pulse height and associated energy
    Returns a tuple of m, b, r value, p value, and std of the fit
    Paramaters:
        ep - event peaks intensity list
        ev - corresponding energy to each event (eV)
        plot - boolean: do you want a plot of the relationship or not
    """
    slope, intercept, r_value, p_value, std_err = stats.linregress(ep,ev)
    if (plot):
        plt.ion()
        plt.figure("Energy versus ADU")
        graphSpace = np.linspace(min(ep), max(ep), 500)
        plt.plot(ev, ep, 'o', label='Data points')
        plt.plot(((graphSpace * slope) + intercept), graphSpace, label='Best fit')
        plt.xlabel("Energy (eV)")
        plt.ylabel("Pulse height (ADU)")
        plt.legend()
        plt.ioff()
        plt.show()
    return [slope, intercept, r_value, p_value, std_err]

def basicADUtoeV(ep, ev, plot):
    """
    Finds the relationship between pulse height and associated energy
    Returns just the slope and intercept, for more on the fit use findADUtoeV()
    Paramaters:
        ep - event peaks intensity list
        ev - corresponding energy to each event (keV)
        plot - boolean: do you want a plot of the relationship or not
    """
    (slope, intercept, r_value, p_value, std_err) = findADUtoeV(ep, ev, plot)
    return [slope, intercept]

def eventsFromFile(fileName, forceNoShape=False):
    """
    Returns the tuple of ep, ex and ey from a file name
    Paramaters:
        fileName - name of file with events (don't include .CSV)
    """
    fP = fileName + ".csv"
    ep, ex, ey, sh = [], [], [], []
    hasShape = False
    with open(fP) as file:
        line = file.readline()
        while line:
            l1 = line
            bits = l1.split(":")
            vals = bits[0].split(",")
            if (str(line) != "ep, ex, ey\n"):
                if (len(bits) == 1):
                    ep.append(int(vals[0]))
                    ex.append(int(vals[1]))
                    ey.append(int(vals[2]))
                else:
                    hasShape = True
                    ep.append(int(vals[0]))
                    ex.append(int(vals[1]))
                    ey.append(int(vals[2]))
                    sh.append(int(bits[1]))
            line = file.readline()
    ep = np.array(ep)
    ex = np.array(ex)
    ey = np.array(ey)
    if (hasShape and (forceNoShape == False)):
        sh = np.array(sh)
        return [ep, ex, ey, sh]
    else:
        return [ep, ex, ey]

def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '='):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar =  fill * filledLength +'>' + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
    # Print New Line on Complete
    if iteration == total: 
        print()
        print()
        
def func_gaussc(x, norm, cent, width, cons):
    return norm*np.exp(-(x-cent)**2/(2*width**2)) + cons
gaussc = lambda p, x: p[0]*np.exp(-(x-p[1])**2/(2*p[2]**2)) + p[3]
        
def fitgaus(ph, cent, fwhm, plot=0) :
    # ph is event pulse height data
    # cent, fwhm are guesses at Gaussian centroid, fwhm
    # histogram the data
    bins = int((max(ph)-min(ph))/1.0)
    (counts, edges) = np.histogram(ph, bins=bins)
    n = len(counts)
    chan = 0.5*(edges[0:n]+edges[1:n+1])
    # guess at peak parameters
    tot = sum(abs(ph-cent) < fwhm) 
    sig = fwhm/2.35
    # fit a gaussian+constant to the data, use Gehrels variance function for error
    pinit = [tot/(sig*np.sqrt(2*np.pi)), cent, sig, 0.0] # use calculated statistics as first guess for fit  
    cerr = 1+np.sqrt(counts+0.75)
    p,t = curve_fit(func_gaussc, chan, counts, p0=pinit, sigma=cerr)
    #print ('Gaussian fit = ', p)
    if (plot > 0) :
        plt.figure(plot)
        plt.clf()
        plt.errorbar(chan, counts, cerr, fmt='.', capsize=0)
        plt.xlabel('Pulse height (ADU)')
        plt.ylabel('Number events/bin')
        #plt.title(filename)
        #s1 = 'Mean = %.1f' % average(v)
        #s2 = 'FWHM = %.1f' % 2.35(std(v))
        #s4 = 'Counts = %i' % len(v)
        #plt.text((max(chan)-min(chan))*0.75+min(chan), max(counts)*0.8, s1+'\n'+s2+'\n'+s3+'\n'+s4)
        plt.plot(chan, gaussc(p, chan))
        plt.show()
    return p, t

def energyRes(E, N):
    """
    Returns energy resolution as a function of energy E and noise N
    """
    return (2.355 * 3.66) * np.sqrt(N**2 + ((3.66 * E)/3.66))

def printEventShape(shape):
    b = np.binary_repr(shape)
    padding = 25 - len(b)
    b = ("0" * padding) + b
    for i in range(0, 5):
        for j in range(0, 5):
            px = 5 * i + j
            if(px == 12):
                sys.stdout.write("X")
            else:
                sys.stdout.write(b[px])
        sys.stdout.write("\n")
    sys.stdout.flush()

def fileAnalysisDebug(filePath, biasPath, darkPath, thres, thresp, savePath, x=2, y=2, autoThres=False, returnShape=False, searchIndex=1, hThr=9, hRep=4):
    t = thres
    tp = thresp
    """
    Performs the analysis of the photos at filePath
    Returns a tuple of event intensity, x location and y location lists
    Paramaters:
        filePath - path to the photos (must be photos taken with C/C++ code, not Python)
        biasPath - path to the bias frames (don't include .FTS in path)
        darkPath - path to the dark frames (don't include .FTS in path)
        thres - threshold intensity for an x-ray to be detected
        thresp - threshold intensity for neighboring pixel to x-ray to be counted
        savePath - where to save the .CSV with the lists (don't include .CSV in path)
    """
    ep, ex, ey = [], [], []
    shp = []
    # set image size
    nx, ny = 1936, 1096
    # set integration region
    dx, dy = x, y
    biasFrame = avgBias(biasPath)
    darkFrame = avgDark(darkPath)
    darkFrame[darkFrame > t] = 0
    files = glob.glob(filePath + "?.FTS") + glob.glob(filePath + "??.FTS") + glob.glob(filePath + "???.FTS")+ glob.glob(filePath + "????.FTS") + glob.glob(filePath + "?????.FTS")
    print("Performing analysis of " + str(len(files)) + " photos")
    if autoThres:
        tp = 2.5 * std_2d(biasFrame + darkFrame)
        t = 5 * std_2d(biasFrame + darkFrame)
    print("thres set to " + str(t))
    print("thresp set to " + str(tp))
    print("dx set to " + str(dx))
    print("dy set to " + str(dy))
    hotPixels = getHotPix(darkPath, thr=hThr, rep=hRep)
    #for ii in range(0, len(hotPixels)):
        #print("x: " + str(hotPixels[ii][1]) + " y: " + str(hotPixels[ii][0]))
    #end len(files)
    problemCount = 0
    prx, pry = [], []
    problemSets = []
    for c in range(0, len(files)):
        data = pyfits.open(files[c])
        if(len(data) == 1):
            imData = data[0].data
        else:
            imData = data[1].data
        data.close()
        phot = imData - (biasFrame + darkFrame)
        phot[phot < 0] = 0
        for ii in range(0, len(hotPixels)):
            phot[hotPixels[ii][0]][hotPixels[ii][1]] = 0
            
        #print(phot)
        nx, ny = phot.shape # find size of image
        f1 = np.reshape(phot, nx*ny) # change the shape to be 1d
        q = np.argsort(f1)[::-1]
        #print(f1[q[0]])
        #print(q[0])
        # go through pixels with large signal
        j = 0
        above_thres = True
        while above_thres :
            i = q[j]
            #print(j, q[j], f1[i])
            j += 1
            if f1[i] >= t : # pixel is above threshold
              #print("i: " + str(i) + " j: " + str(j))
              x = (i % 1936) + 1
              y = math.floor((i / 1936) + 1)
              yR = (i % 1936)
              xR = math.floor((i / 1936))
              #print("Checking pixel " + str(x) + ", " + str(y) + " with value " + str(phot[xR][yR]))
              # i = ny*x + y # index into 1d array
              # check if pixel is far enough from edge of chip
              if (xR > dx) and (xR < nx-dx-1) and (yR > dy) and (yR < ny-dy-1) :
                #p = sum(f1[i+di])
                p = 0
                area = phot[(xR-dx):(xR+dx+1), (yR-dy):(yR+dy+1)]
                pAdjust, s = pixelAnalysis(area, tp)
                for xi in range(xR-dx, xR+dx+1) :
                  for yi in range(yR-dy, yR+dy+1):
                      if (phot[xi, yi] > tp):
                          p += phot[xi,yi]
                          if(pAdjust == 0):
                              phot[xi, yi] = 0 # zero out pixels in this event
                          elif(xi-xR+dx != 0 and xi-xR+dx != 4 and yi-yR+dy != 0 and yi-yR+dy != 4):
                              phot[xi, yi] = 0
                if(p < 90):
                    #print("Problematic event of %d at %d, %d on file %s. Adjust: %d" % (p, xR, yR, files[c][45:], pAdjust))
                    problemCount += 1
                    if([xR, yR] in problemSets):
                        print("Repeated problem pixel at %d, %d" % (xR, yR))
                    prx.append(xR)
                    pry.append(yR)
                    problemSets.append([xR, yR])
                if(p == 3821):
                    print("That weird pulse at %d %d on image %s" % (xR, yR, files[c]))
                if (p > 0):
                    ep.append(p - pAdjust)
                    ex.append(x)
                    ey.append(y)
                    shp.append(s)
                    #if (p > 3000):
                        #print("Potential cosmic ray of ADU " + str(p) + " on image " + str(files[c]))
                #print (x, y, p)
                #f1[i+di] = 0 # zero out pixels in this event
            # if we reach a pixel below threshold than has not been zeroed, then we are done
            elif f1[i] > 0 :
              above_thres = False # stop looping
        prog = (c/len(files)) * 100
        printProgressBar(prog, 100, prefix = 'Progress:', suffix='Complete', length=50)
          
    ep = np.array(ep)
    ex = np.array(ex)
    ey = np.array(ey)
    shp = np.array(shp)
    fi = open(savePath + '.csv', 'w')
    fi.write('ep, ex, ey\n')
    for cntr in range(0, len(ep)):
        fi.write(str(ep[cntr]) + ', ' + str(ex[cntr]) + ', ' + str(ey[cntr]) + ':' + str(shp[cntr]) + "\n")
    fi.close()
    '''
    print(str(len(ep)) + " events written")
    print(str(problemCount) + " total problematic events")
    badPointZip = dict(zip(prx, pry))
    badSort = sorted(badPointZip.items(), key=lambda kv: kv[0])[::-1]
    print(len(badSort))
    for a in range(0, len(badSort)):
        print(badSort[a])
    '''
    if(returnShape):
        return [ep, ex, ey, shp]
    else:
        return [ep, ex, ey]
    
def pixelAnalysis(area, tp):
    shp = ""
    truePulse = 0
    smallPulse = 0
    for i in range(4, -1, -1):
        for j in range(0, 5):
            if(area[i][j] > tp):
                shp += "1"
                truePulse += area[i][j]
            else:
                shp += "0"
    s = int(shp, 2)
    newShp = ""
    for t in range(0, len(badShapes)):
        if (str(s) == str(badShapes[t])):
            for i in range(4, -1, -1):
                for j in range(0, 5):
                    if(i == 0 or i == 4 or j == 0 or j == 4):
                        newShp += "0"
                    elif(area[i][j] > tp):
                        newShp += "1"
                        smallPulse += area[i][j]
                    else:
                        newShp += "0"
    if(newShp == ''):
        return [0, s]
    else:
        nS = int(newShp, 2)
        return [(truePulse - smallPulse), nS]
      
      
      
def lorentzfit(data,amp,cent,fwhm):
    m = 0.0037168339915133394
    b = 0.002116795017471418
    bins = int((max(data) - min(data))/2)
    (counts, edges) = np.histogram(data, bins=bins)
    n = len(counts)
    chan = 0.5 * (edges[0:n] + edges[1:n + 1])
    chan = aduToeV(chan, m, b)
    sig = fwhm / 2
    cerr = 1 + np.sqrt(counts + 0.75)
    cerr = cerr[:-1]
    chan = chan[:-1]
    counts = counts[:-1]
    lorentzian = lmfit.models.LorentzianModel(prefix = 'l1_')
    pars = lorentzian.make_params()
    pars['l1_amplitude'].set(value = amp,min = 0)
    pars['l1_center'].set(value=cent, min=min(chan),max = max(chan))
    pars['l1_sigma'].set(value=sig, min=0)
    out = lorentzian.fit(counts, pars, x=chan,weights = 1/cerr)
    r_chisq = out.redchi
    l1_norm = out.best_values['l1_amplitude']
    l1_cent = out.best_values['l1_center']
    l1_sigma = out.best_values['l1_sigma']
    l1_fwhm = 2*l1_sigma
    covar = out.covar
    #err_norm = np.sqrt(abs(covar[0][0]))
    err_cent = np.sqrt(abs(covar[1][1]))
    err_sigma = np.sqrt(abs(covar[2][2]))
    err_fwhm = np.sqrt((err_sigma*2)**2)
    height = 0.3183099 * l1_norm / max(2.220446049250313e-16, l1_sigma)
    #err_h = height*0.3183099* np.sqrt(((err_norm)/(l1_norm))**2 + (
            #((err_sigma) / (l1_sigma))**2))
    print(out.fit_report(min_correl=0.5))
    fig, axes = plt.subplots(1, figsize=(10, 8))
    textstr = '\n'.join((
        r'Counts= %.2f' % (height),
        r'Centroid = %.2f +- %.2f eV' % (l1_cent,err_cent),
        r'FWHM = %.2f +- %.2f eV' % (l1_fwhm,err_fwhm),
        r'$\chi^2/dov $= %.2f' % (r_chisq)))
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    axes.text(.22, .89, textstr, fontsize=8, ha='center', va='center', transform=axes.transAxes, bbox=props)
    axes.errorbar(chan, counts, cerr, fmt='.')
    axes.plot(chan, out.best_fit, 'r-', label='best fit')
    axes.legend(loc='best')
    return (height,l1_cent,l1_fwhm,r_chisq)


def voigtfit(data,amplitude,cent,fwhm):
    m = 0.0037168339915133394
    b = 0.002116795017471418
    bins = int((max(data) - min(data))/2)
    (counts, edges) = np.histogram(data, bins=bins)
    n = len(counts)
    chan = 0.5 * (edges[0:n] + edges[1:n + 1])
    chan = aduToeV(chan, m, b)
    sig = fwhm / 2
    cerr = 1 + np.sqrt(counts + 0.75)
    cerr = cerr[:-1]
    chan = chan[:-1]
    counts = counts[:-1]
    lorentzian = lmfit.models.VoigtModel(prefix = 'v1_')
    pars = lorentzian.make_params()
    pars['v1_amplitude'].set(value= amplitude,min = 0)
    pars['v1_center'].set(value=cent, min=min(chan),max = max(chan))
    pars['v1_gamma'].set(value = sig)
    pars['v1_sigma'].set(value=sig, min=0)
    out = lorentzian.fit(counts, pars, x=chan,weights = 1/cerr)
    r_chisq = out.redchi
    v1_norm = out.best_values['v1_amplitude']
    v1_cent = out.best_values['v1_center']
    v1_sigma = out.best_values['v1_sigma']
    v1_gamma = out.best_values['v1_gamma']
    v1_fwhm = 1.0692*v1_gamma+np.sqrt(0.8664*v1_gamma**2+5.545083*v1_sigma**2)
    covar = out.covar
    #err_norm = np.sqrt(abs(covar[0][0]))
    err_cent = np.sqrt(abs(covar[1][1]))
    err_sigma = np.sqrt(abs(covar[2][2]))
    #err_fwhm = np.sqrt((err_sigma*2)**2)
    height = (v1_norm/(max(2.220446049250313e-16, v1_sigma*np.sqrt(2*np.pi))))*special.wofz((1j*v1_gamma)/(max(2.220446049250313e-16, v1_sigma*np.sqrt(2)))).real
    print(out.fit_report(min_correl=0.5))
    fig, axes = plt.subplots(1, figsize=(10, 8))
    textstr = '\n'.join((
        r'Height= %.2f' % (height),
        r'Centroid = %.2f +- %.2f eV' % (v1_cent,err_cent),
        r'FWHM = %.2f eV' % (v1_fwhm),
        r'$\chi^2/dov $= %.2f' % (r_chisq)))
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    axes.text(.22, .89, textstr, fontsize=8, ha='center', va='center', transform=axes.transAxes, bbox=props)
    axes.errorbar(chan, counts, cerr, fmt='.')
    axes.plot(chan, out.best_fit, 'r-', label='best fit')
    axes.legend(loc='best')
    return (height,v1_cent,v1_fwhm,r_chisq)
