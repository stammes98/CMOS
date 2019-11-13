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


def takePhotos(nPhotos, sleep, exposure, fileName, FPS=1, pixelClock=20, gain=30):
    gc.set_threshold(200,5,5)
    #process = psutil.Process(os.getpid())
    """
    Takes photos and saves them as FTS files
    Paramaters:
        nPhotos - number of photos (int)
        sleep - how long to sleep between taking photos
        exposure - exposure time of photo (ms)
        FPS - frames per second
        pixelClock - pixel clock setting (int)
        gain - hardware gain (0-100)
        fileName - path to file. do not include the .FTS ending
    """
    hCam = ueye.HIDS(0)             #0: first available camera;  1-254: The camera with the specified camera ID
    sInfo = ueye.SENSORINFO()
    cInfo = ueye.CAMINFO()
    pcImageMemory = ueye.c_mem_p()
    MemID = ueye.int()
    rectAOI = ueye.IS_RECT()
    pitch = ueye.INT()
    nBitsPerPixel = ueye.INT(8)    #24: bits per pixel for color mode; take 8 bits per pixel for monochrome
    m_nColorMode = ueye.INT()		# Y8/RGB16/RGB24/REG32
    bytes_per_pixel = int(nBitsPerPixel / 8)
    print("START")
    print()
    
    # Starts the driver and establishes the connection to the camera
    nRet = ueye.is_InitCamera(hCam, None)
    if nRet != ueye.IS_SUCCESS:
        print("is_InitCamera ERROR")
        sys.exit("Init Error. Is the camera plugged in?")
    
    # Reads out the data hard-coded in the non-volatile camera memory and writes it to the data structure that cInfo points to
    nRet = ueye.is_GetCameraInfo(hCam, cInfo)
    if nRet != ueye.IS_SUCCESS:
        print("is_GetCameraInfo ERROR")
        #exit(0)
    
    
    # You can query additional information about the sensor type used in the camera
    nRet = ueye.is_GetSensorInfo(hCam, sInfo)
    if nRet != ueye.IS_SUCCESS:
        print("is_GetSensorInfo ERROR")
    
    nRet = ueye.is_ResetToDefault(hCam)
    if nRet != ueye.IS_SUCCESS:
        print("is_ResetToDefault ERROR")
    
    # Set display mode to DIB
    nRet = ueye.is_SetDisplayMode(hCam, ueye.IS_SET_DM_DIB)
    
    # Set the right color mode
    if int.from_bytes(sInfo.nColorMode.value, byteorder='big') == ueye.IS_COLORMODE_BAYER:
        # setup the color depth to the current windows setting
        ueye.is_GetColorDepth(hCam, nBitsPerPixel, m_nColorMode)
        bytes_per_pixel = int(nBitsPerPixel / 8)
        print("IS_COLORMODE_BAYER: ", )
        print("\tm_nColorMode: \t\t", m_nColorMode)
        print("\tnBitsPerPixel: \t\t", nBitsPerPixel)
        print("\tbytes_per_pixel: \t\t", bytes_per_pixel)
        print()
    
    elif int.from_bytes(sInfo.nColorMode.value, byteorder='big') == ueye.IS_COLORMODE_CBYCRY:
        # for color camera models use RGB32 mode
        m_nColorMode = ueye.IS_CM_BGRA8_PACKED
        nBitsPerPixel = ueye.INT(32)
        bytes_per_pixel = int(nBitsPerPixel / 8)
        print("IS_COLORMODE_CBYCRY: ", )
        print("\tm_nColorMode: \t\t", m_nColorMode)
        print("\tnBitsPerPixel: \t\t", nBitsPerPixel)
        print("\tbytes_per_pixel: \t\t", bytes_per_pixel)
        print()
    
    elif int.from_bytes(sInfo.nColorMode.value, byteorder='big') == ueye.IS_COLORMODE_MONOCHROME:
        # for color camera models use RGB32 mode
        m_nColorMode = ueye.IS_CM_MONO12
        nBitsPerPixel = ueye.INT(12)
        bytes_per_pixel = int(nBitsPerPixel / 8)
        print("IS_COLORMODE_MONOCHROME: ", )
        print("\tm_nColorMode: \t\t", m_nColorMode)
        print("\tnBitsPerPixel: \t\t", nBitsPerPixel)
        print("\tbytes_per_pixel: \t\t", bytes_per_pixel)
        print()
    
    else:
        # for monochrome camera models use Y8 mode
        m_nColorMode = ueye.IS_CM_MONO8
        nBitsPerPixel = ueye.INT(8)
        bytes_per_pixel = int(nBitsPerPixel / 8)
        print("else")
    
    # Can be used to set the size and position of an "area of interest"(AOI) within an image
    nRet = ueye.is_AOI(hCam, ueye.IS_AOI_IMAGE_GET_AOI, rectAOI, ueye.sizeof(rectAOI))
    if nRet != ueye.IS_SUCCESS:
        print("is_AOI ERROR")
    
    width = rectAOI.s32Width
    height = rectAOI.s32Height
    
    
    # Prints out some information about the camera and the sensor
    print("Camera model:\t\t", sInfo.strSensorName.decode('utf-8'))
    print("Camera serial no.:\t", cInfo.SerNo.decode('utf-8'))
    print("Maximum image width:\t", width)
    print("Maximum image height:\t", height)
    print()
    
    
    nRet = ueye.is_AllocImageMem(hCam, width, height, nBitsPerPixel, pcImageMemory, MemID)
    if nRet != ueye.IS_SUCCESS:
        print("is_AllocImageMem ERROR")
    else:
        # Makes the specified image memory the active memory
        nRet = ueye.is_SetImageMem(hCam, pcImageMemory, MemID)
        if nRet != ueye.IS_SUCCESS:
            print("is_SetImageMem ERROR")
        else:
            # Set the desired color mode
            nRet = ueye.is_SetColorMode(hCam, m_nColorMode)
    
    # Enables the queue mode for existing image memory sequences
    nRet = ueye.is_InquireImageMem(hCam, pcImageMemory, MemID, width, height, nBitsPerPixel, pitch)
    if nRet != ueye.IS_SUCCESS:
        print("is_InquireImageMem ERROR")
    else:
        print("Memory queued")
        
    
    nRet = 0
    nCurrentState = ueye.UINT(0)
    nRet = ueye.is_IO(hCam, ueye.IS_IO_CMD_LED_GET_STATE, nCurrentState, ueye.sizeof(nCurrentState))
    nCurrentState = ueye.UINT(3)
    nRet = ueye.is_IO(hCam, ueye.IS_IO_CMD_LED_SET_STATE, nCurrentState, ueye.sizeof(nCurrentState))
    if (nRet == ueye.IS_SUCCESS):
        print("LED disabled")
    
    fps = ueye.DOUBLE(FPS)
    newFps = ueye.double(0.0)
    nRet = ueye.is_SetFrameRate(hCam, fps, newFps)
    print("Current FPS: " + str(newFps))
    
    
    pixClock = ueye.INT(pixelClock)
    nRet = ueye.is_PixelClock(hCam, ueye.IS_PIXELCLOCK_CMD_SET, pixClock, ueye.sizeof(pixClock))
    print("Pixel clock (MHz): " + str(pixClock))
    
    
    if (exposure > 999):
        long = ueye.UINT(1)
        nRet = ueye.is_Exposure(hCam, ueye.IS_EXPOSURE_CMD_SET_LONG_EXPOSURE_ENABLE, long, ueye.sizeof(long))
    
    hardGain = ueye.INT(gain)
    nRet = ueye.is_SetHardwareGain(hCam, hardGain, hardGain, hardGain, hardGain)
    print("Hardware gain setting: " + str(hardGain))
    
    expo = ueye.DOUBLE(exposure)
    nRet = ueye.is_Exposure(hCam, ueye.IS_EXPOSURE_CMD_SET_EXPOSURE, expo, ueye.sizeof(expo))
    print("Exposure time (ms): " + str(expo))
    
    time.sleep((int((expo/1000))))
    for i in range(0,nPhotos + 1):
        try:
            time.sleep(sleep / 2)
            nRet = ueye.is_FreezeVideo(hCam, ueye.IS_DONT_WAIT)
            time.sleep((int((expo/1000))))
            if (i == 0):
                continue
            print("Taking image #" + str(i))
            if(nRet == ueye.IS_SUCCESS):
                # In order to display the image in an OpenCV window we need to...
                # ...extract the data of our image memory
                nRet = 0
                dInfo = ueye.IS_DEVICE_INFO()
                nRet = ueye.is_DeviceInfo(hCam, ueye.IS_DEVICE_INFO_CMD_GET_DEVICE_INFO, dInfo, ueye.sizeof(dInfo))
                temp = int(dInfo.infoDevHeartbeat.wTemperature)
                tempBin = np.binary_repr(temp, width=15)
                tempStr = str(tempBin)[::-1]
                tempDec = tempStr[0:3][::-1]
                tempPreDec = tempStr[4:10][::-1]
                tempSign = tempStr[14]
                tempPreDecReal = int(tempPreDec, 2)
                tempDecReal = int(tempDec, 2)
                if (int(tempSign) == 0):
                    tempReal = str(tempPreDecReal) + '.' + str(tempDecReal)
                else:
                    tempReal = '-' + str(tempPreDecReal) + '.' + str(tempDecReal)
                print("Camera temp is " + tempReal + "*C")
                temp = None
                tempBin = None
                tempStr = None
                tempDec = None
                tempPreDec = None
                tempSign = None
                tempPreDecReal = None
                tempDecReal = None
                dInfo = None
    
                #if (float(tempReal) >= 35.1):
                #    sys.exit("Camera is overheating")
                array = ueye.get_data(pcImageMemory, width, height, nBitsPerPixel, pitch, copy=False)
        
                bytes_per_pixel = 2
                # ...reshape it in an numpy array...
                print (array.shape)
                frame = np.reshape(array,(height.value, width.value, bytes_per_pixel))
                print(frame.shape)
                frame1 = frame[:,:,0] + (256 * frame[:,:,1])
                print(frame1.shape)
                imStats(frame1)
        
                
                head = pyfits.Header()
                head.append(('Time', str(datetime.datetime.now()),'YEAR-MONTH-DAY HOUR:MINUTE:SECOND'))
                head.append(('Exposure', round(float(expo), 4), 'Exposure time in ms'))
                head.append(('FPS', round(float(newFps),4), 'Frames / second'))
                head.append(('Clock', pixelClock, 'Pixel clock in MHz'))
                head.append(('Sleep', sleep, 'Seconds between each photo taken'))
                head.append(('Number', (i), 'Which photo in the batch is this one'))
                if int.from_bytes(sInfo.nColorMode.value, byteorder='big') == ueye.IS_COLORMODE_BAYER:
                    head.append(('Mode', 'IS_COLORMODE_BAYER', 'Camera Mode'))
                elif int.from_bytes(sInfo.nColorMode.value, byteorder='big') == ueye.IS_COLORMODE_CBYCRY:
                    # for color camera models use RGB32 mode
                    head.append(('Mode', 'IS_COLORMODE_CBYCRY', 'Camera Mode'))
                elif int.from_bytes(sInfo.nColorMode.value, byteorder='big') == ueye.IS_COLORMODE_MONOCHROME:
                    # for color camera models use RGB32 mode
                    head.append(('Mode', 'IS_COLORMODE_MONOCHROME', 'Camera Mode'))
                else:
                    head.append(('Mode', 'OTHER', 'Camera Mode'))
                head.append(('Model', sInfo.strSensorName.decode('utf-8'), 'Camera Model'))
                head.append(('Serial', cInfo.SerNo.decode('utf-8'), 'Camera serial number'))
                head.append(('Width', int(width), 'Maximum width of image'))
                head.append(('Height', int(height), 'Maximum height of image'))
                head.append(('BPS', int(nBitsPerPixel), 'Bits per pixel'))
                head.append(('Temp', float(tempReal), 'Temperature of camera when image was taken (C)'))
                head.append(('Author', 'Steve Tammes', 'Author of the code to take these images'))
                head.append(('Gain', int(gain), 'Gain'))
                date = str(datetime.datetime.now().year) + '-' + str(datetime.datetime.now().month) + '-' + str(datetime.datetime.now().day) + ' ' + str(datetime.datetime.now().hour) + '-' + str(datetime.datetime.now().minute) + '-' + str(datetime.datetime.now().second)
                head.append(('Time', date, 'When the image was taken'))
                pyfits.writeto(fileName + str(i) + '.FTS', frame1, head, overwrite=True)
                frame = None
                frame1 = None
                array = None
                gc.collect()
                """
                print (str(delItems) + " items cleared from memory")
                print("Garbage counts: " + str(gc.get_count()))
                print("img refrence count: " + str(sys.getrefcount(frame1)))
                print("Memory usage: " + str(process.memory_info().rss))
                print()
                """
                time.sleep(sleep/2)
        except (KeyboardInterrupt):
            break
            print("Cleaning up")
            return
            
    
    ueye.is_FreeImageMem(hCam, pcImageMemory, MemID)
    
    # Disables the hCam camera handle and releases the data structures and memory areas taken up by the uEye camera
    ueye.is_ExitCamera(hCam)
    
    print()
    print("END")


def avgBias(biasPath):
    """
    Averages out a set of bias frames and returns it
    Paramaters:
        biasPath - path to the bias images, include file name but not ?.FTS
    """
    files = glob.glob(biasPath + "?.FTS") + glob.glob(biasPath + "??.FTS") + glob.glob(biasPath + "???.FTS")
    shaping = pyfits.open(files[0])
    shape = shaping[0].data
    shaping.close()
    medBias = np.empty_like(shape)
    totBias = np.empty_like(medBias)
    for i in range(0, len(files)):
        data = pyfits.open(files[i])
        imgData = data[0].data
        totBias += imgData
        data.close()
    medBias = totBias / len(files)
    medBias = medBias.astype(int)
    print("Bias frame compiled")
    return medBias
    
def avgDark(darkPath):
    """
    Averages out a set of dark frames
    Paramaters:
        darkPath - path to the bias images, include file name but not ?.FTS
    """
    files = glob.glob(darkPath + "?.FTS") + glob.glob(darkPath + "??.FTS") + glob.glob(darkPath + "???.FTS")
    shaping = pyfits.open(files[0])
    shape = shaping[0].data
    shaping.close()
    medDark = np.empty_like(shape)
    totDark = np.empty_like(medDark)
    for i in range(0, len(files)):
        data = pyfits.open(files[i])
        imgData = data[0].data
        totDark += imgData
        data.close()
    medDark = totDark / len(files)
    medDark = medDark.astype(int)
    print("Dark frame compiled")
    return medDark


def activeAnalysis(filePath, biasPath, darkPath, thres, thresp, savePath, x=2, y=2, autoThres=False):
    t = thres
    tp = thresp
    """
    Performs the analysis of the photos at filePath
    Returns a tuple of event intensity, x location and y location lists
    Paramaters:
        filePath - path to the photos (must be photos taken with C++ code, not Python)
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
    darkFrame[darkFrame > t] = 0
    files = glob.glob(filePath + "*.FTS")
    print("Performing analysis of " + str(len(files)) + " photos")
    if autoThres:
        tp = 2.5 * std_2d(biasFrame + darkFrame)
        t = 5 * std_2d(biasFrame + darkFrame)
    print("thres set to " + str(t))
    print("thresp set to " + str(tp))
    print("dx set to " + str(dx))
    print("dy set to " + str(dy))
    #end len(files)
    for c in range(0, len(files)):
        data = pyfits.open(files[c])
        imData = data[0].data
        data.close()
        phot = imData - (biasFrame + darkFrame)
        phot[phot < 0] = 0
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
              x = int(math.floor(i/1936) +1)
              y = (i % 1936) + 1
              # i = ny*x + y # index into 1d array
              # check if pixel is far enough from edge of chip
              if (x > dx) and (x < nx-dx-1) and (y > dy) and (y < ny-dy-1) :
                #p = sum(f1[i+di])
                p = 0
                for xi in range(x-dx, x+dx+1) :
                  for yi in range(y-dy, y+dy+1):
                  	if phot[xi, yi] > tp : p += phot[xi,yi]
                  	phot[xi, yi] = 0 # zero out pixels in this event
                #if (p > 490 and p < 515):
                    #print(str(p) + " pulse at " + str(y) + ", " + str(x) + " on image " + str(files[c]) + " - c: " + str(c))
                    
                #print("pulse at " + str(y) + ", " + str(x) + " - " + str(p))
                ep.append(p)
                ex.append(y)
                ey.append(x)
                #print (x, y, p)
                #f1[i+di] = 0 # zero out pixels in this event
            # if we reach a pixel below threshold than has not been zeroed, then we are done
            elif f1[i] > 0 :
              above_thres = False # stop looping
        prog = (c/len(files)) * 100
        printProgressBar(prog, 100, prefix = 'Progress:', suffix='Complete', length=50)
          
    print("100% complete")      
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


#DON'T USE I'M NOT DONE WITH THIS ONE YET
def pythonAnalysis(filePath, biasPath, darkPath, thres, thresp, savePath, x=1, y=1, autoThres=False):
    t = thres
    tp = thresp
    """
    Performs the analysis of the photos listed at the filePath
    Returns a tuple of event intensity, x location and y location lists
    Paramaters:
        filePath - path to the image files (must be images taken with Python code, not C++)
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
    darkFrame[darkFrame > t] = 0
    files = glob.glob(filePath + "*.FTS")
    print("Performing analysis of " + str(len(files)) + " photos")
    if autoThres:
        tp = 2.5 * std_2d(biasFrame + darkFrame)
        t = 5 * std_2d(biasFrame + darkFrame)
    print("thresp set to " + str(thresp))
    print("thres set to " + str(thres))
    #end len(files)
    for c in range(0, len(files)):
        data = pyfits.open(files[c])
        imData = data[0].data
        data.close()
        phot = imData - (biasFrame + darkFrame)
        phot[phot < 0] = 0
        #print(phot)
        nx, ny = phot.shape # find size of image
        f1 = np.reshape(phot, nx*ny) # change the shape to be 1d
        q = np.argsort(f1)[::-1]
        print(f1[q[0]])
        print(q[0])
        # go through pixels with large signal
        j = 0
        above_thres = True
        while above_thres :
            i = q[j]
            #print(j, q[j], f1[i])
            j += 1
            if f1[i] >= t : # pixel is above threshold
              print("i: " + str(i) + " j: " + str(j))
              x = int(math.floor(i/1096) + 1)
              y = (i % 1096)
              # i = ny*x + y # index into 1d array
              # check if pixel is far enough from edge of chip
              if (x > dx) and (x < nx-dx-1) and (y > dy) and (y < ny-dy-1) :
                #p = sum(f1[i+di])
                p = 0
                for xi in range(x-dx, x+dx+1) :
                  for yi in range(y-dy, y+dy+1):
                  	if phot[yi-1, xi-1] > tp : p += phot[yi-1,xi-1]
                  	phot[yi-1, xi-1] = 0 # zero out pixels in this event
                if (p > 490 and p < 515):
                    print(str(p) + " pulse at " + str(x) + ", " + str(y) + " on image " + str(files[c]) + " - C: " + str(c))
                    
                print("pulse at " + str(x) + ", " + str(y) + " - " + str(p))
                ep.append(p)
                ex.append(x)
                ey.append(y)
                #print (x, y, p)
                #f1[i+di] = 0 # zero out pixels in this event
            # if we reach a pixel below threshold than has not been zeroed, then we are done
            elif f1[i] > 0 :
              above_thres = False # stop looping
          
    print("100% complete")      
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



def fileAnalysis(photoPath, biasPath, darkPath, thres, thresp, savePath, dx=1, dy=1, autoThres=False):
    """
    Performs the analysis of a list of photos saved as FTS files
    Paramaters:
        photoPath - path to the photos (don't include .FTS in path)
        biasPath - path to the bias frames (don't include .FTS in path)
        darkPath - path to the dark frames (don't include .FTS in path)
        thres - threshold intensity for an x-ray to be detected
        thresp - threshold intensity for neighboring pixel to x-ray to be counted
    """
        
    return activeAnalysis(photoPath, biasPath, darkPath, thres, thresp, savePath, dx, dy, autoThres)


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
    plt.hist(ep, bins=b, range=[0,max(unis)])
    if(saveImg):
        plt.savefig(saveName)
    plt.ioff()
    plt.show()

def spectrum(ep, m, b, saveImg=False, saveName="default.png"):
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
    ydat = smoothing(y, 2)
    #plt.rc('xtick', labelsize=8) 
    plt.xlabel("Energy (eV)")
    plt.ylabel("Intensity")
    plt.plot(xdat, ydat)
    plt.gca().xaxis.set_major_locator(MultipleLocator(200))
    plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%d'))
    plt.gca().xaxis.set_minor_locator(MultipleLocator(50))
    plt.xticks(np.arange(100, max(xdat), 200.0))
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
    plt.xticks(np.arange(min(xdat1), max(xdat1)+1, 80.0))
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

def eventsFromFile(fileName):
    """
    Returns the tuple of ep, ex and ey from a file name
    Paramaters:
        fileName - name of file with events (don't include .CSV)
    """
    (ep, ex, ey) = np.loadtxt(fileName + '.csv', delimiter=',', unpack=True, skiprows=1)    
    return [ep, ex, ey]

def readMultimeter(fileName):
    """
    Reads the .csv file output by the multimeter program, returns a tuple of the datetime of the reading, the value of the reading and a string of the unit (ex: mA dc)
    Paramaters:
        fileName - path to the file (don't include .CSV)
    """
    date, reading, unit = [], [], []
    f = open(fileName + '.csv')
    for line in f:
        if (line.startswith("Time")):
            continue
        elif (line.startswith("Date:")):
            continue
        else:
            s1 = line.split(',')
            for i in range(0, len(s1)):
                hour, minute, sec = s1[0].split(':')
                year, month, day = s1[1].split('/')
                year = "20" + year
                time = datetime.datetime(int(year), int(month), int(day), int(hour), int(minute), int(sec))
                reading.append(float(s1[2].replace(" ", "")))
                unit.append(s1[3])
                date.append(time)
    f.close()
    return [date, reading, unit]

def plotMultimeter(date, reading, unit):
    """
    Plots the output of a multimeter's readings
    Paramaters:
        date - list of datetime elements when each reading was taken
        reading - list of the readings from the multimeter
        unit - list of the unit of each reading from the multimeter
    """
    graphTime = []
    graphTime.append(0)
    for i in range(1, len(date)):
        diff = date[i] - date[0]
        diff = diff.total_seconds()
        graphTime.append(diff)
    
    plt.ion()
    plt.figure("Current vs time")
    plt.xlabel("Seconds")
    plt.ylabel(unit[0])
    plt.clf()
    plt.plot(graphTime, reading)
    plt.ioff()    
    plt.show()
    
def plotMultimeterVsPhotos(date, reading, unit, fileName):
    """
    Plots the output of a multimeter's readings with vertical lines to indicate when the photos in the batch were taken
    Paramaters:
        date - list of datetime elements when each reading was taken
        reading - list of the readings from the multimeter
        unit - list of the unit of each reading from the multimeter
        fileName - path to the images taken during the multimeter readings (don't include .FTS in path)
    """
    graphTime = []
    graphTime.append(0)
    for i in range(1, len(date)):
        diff = date[i] - date[0]
        diff = diff.total_seconds()
        graphTime.append(diff)
    plt.ion()
    plt.figure("Current vs time")
    plt.clf()
    plt.xlabel("Seconds")
    plt.ylabel(unit[0])
    plt.plot(graphTime, reading)
    pics = glob.glob(fileName + "?.FTS") + glob.glob(fileName + "??.FTS") + glob.glob(fileName + "???.FTS")  + glob.glob(fileName + "????.FTS")
    print("Comparing multimeter to " + str(len(pics)) + " photos")
    picTimes = []
    for j in range(0, len(pics)):
        f = pyfits.open(pics[j])
        head = f[0].header
        fitsDate = head['Time']
        calender, time = fitsDate.split(" ")
        hour, minute, second = time.split(":")
        second = second[0:2]
        year, month, day = calender.split("-")
        photoTime = datetime.datetime(int(year), int(month), int(day), int(hour), int(minute), int(second))
        deltaTime = photoTime - date[0]
        deltaTime = deltaTime.total_seconds()
        picTimes.append(deltaTime)
        plt.axvline(deltaTime, color='k', linewidth='0.2')
        f.close()
    plt.ioff()
    plt.show()

def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█'):
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
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
    # Print New Line on Complete
    if iteration == total: 
        print()    
        
def func_gaussc(x, norm, cent, width, cons):
    return norm*np.exp(-(x-cent)**2/(2*width**2)) + cons
gaussc = lambda p, x: p[0]*np.exp(-(x-p[1])**2/(2*p[2]**2)) + p[3]
        
def fitgaus(ph, cent, fwhm, plot=0) :
    # ph is event pulse height data
    # cent, fwhm are guesses at Gaussian centroid, fwhm
    # histogram the data
    bins = int((max(ph)-min(ph))/4.0)
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
    print ('Gaussian fit = ', p)
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
