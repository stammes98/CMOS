# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 09:40:44 2019

@author: stamm
"""

import sys
import lmfit
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
import time
import gc
from scipy.special import erf
from scipy import special


def aduToeV(adu, m, b):
    return ((m * adu) + b) * 1000


def smoothing(data, pts):
    box = np.ones(pts) / pts
    smooth = np.convolve(data, box, mode='same')
    return smooth


def mean_2d(img):
    """
    Find the mean of a 2d image
    Paramaters:
        img - the image

    """
    # img is a 2-d array, need to change to 1-d to do calculations

    nx, ny = img.shape  # find the size of the array
    img1 = np.reshape(img, nx * ny)  # change the shape to be 1d
    return np.mean(img1)


# find the median of a 2d image
def median_2d(img):
    """
    Find the median of a 2d image
    Paramaters:
      img - the image

    """
    # img is a 2-d array, need to change to 1-d to do calculations
    nx, ny = img.shape  # find the size of the array
    img1 = np.reshape(img, nx * ny)  # change the shape to be 1d
    return np.median(img1)


# find the standard deviation of a 2d image
def std_2d(img):
    """
    Find the std of a 2d image
    Paramaters:
        img - the image

    """
    # img is a 2-d array, need to change to 1-d to do calculations
    nx, ny = img.shape  # find the size of the array
    img1 = np.reshape(img, nx * ny)  # change the shape to be 1d
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


"""
def takePhotos(nPhotos, sleep, exposure, fileName, FPS=1, pixelClock=20, gain=30):
    gc.set_threshold(200,5,5)
    #process = psutil.Process(os.getpid())
    Takes photos and saves them as FTS files
    Paramaters:
        nPhotos - number of photos (int)
        sleep - how long to sleep between taking photos
        exposure - exposure time of photo (ms)
        FPS - frames per second
        pixelClock - pixel clock setting (int)
        gain - hardware gain (0-100)
        fileName - path to file. do not include the .FTS ending
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
                print (str(delItems) + " items cleared from memory")
                print("Garbage counts: " + str(gc.get_count()))
                print("img refrence count: " + str(sys.getrefcount(frame1)))
                print("Memory usage: " + str(process.memory_info().rss))
                print()
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
    """





def avgBias(biasPath):
    """
    Averages out a set of bias frames and returns it
    Paramaters:
        biasPath - path to the bias images, include file name but not ?.FTS
    """
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
    return medBias


def avgDark(darkPath):
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
    print("Dark mean: " + str(mean_2d(medDark)))
    print("Dark frame compiled from " + str(len(files)) + " frames")
    return medDark


def getFilesForBrute(filePath, bar=True):
	fs = glob.glob(filePath + "*")
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



def fileAnalysisBrute(phots, biasPath, darkPath, thres, thresp, x=2, y=2, autoThres=False):
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
    # end len(files)
    for c in range(0, len(phots)):
        phot = phots[c] - (biasFrame + darkFrame)
        phot[phot < 0] = 0
        # print(phot)
        nx, ny = phot.shape  # find size of image
        f1 = np.reshape(phot, nx * ny)  # change the shape to be 1d
        q = np.argsort(f1)[::-1]
        # print(f1[q[0]])
        # print(q[0])
        # go through pixels with large signal
        j = 0
        above_thres = True
        while above_thres:
            i = q[j]
            # print(j, q[j], f1[i])
            j += 1
            if f1[i] >= t:  # pixel is above threshold
                # print("i: " + str(i) + " j: " + str(j))
                x = int(math.floor(i / 1936) + 1)
                y = (i % 1936) + 1
                # i = ny*x + y # index into 1d array
                # check if pixel is far enough from edge of chip
                if (x > dx) and (x < nx - dx - 1) and (y > dy) and (y < ny - dy - 1):
                    # p = sum(f1[i+di])
                    p = 0
                    for xi in range(x - dx, x + dx + 1):
                        for yi in range(y - dy, y + dy + 1):
                            if phot[xi, yi] > tp: p += phot[xi, yi]
                            phot[xi, yi] = 0  # zero out pixels in this event
                    # if (p > 490 and p < 515):
                    # print(str(p) + " pulse at " + str(y) + ", " + str(x) + " on image " + str(files[c]) + " - c: " + str(c))

                    # print("pulse at " + str(y) + ", " + str(x) + " - " + str(p))
                    ep.append(p)
                    ex.append(y)
                    ey.append(x)
                    # print (x, y, p)
                    # f1[i+di] = 0 # zero out pixels in this event
            # if we reach a pixel below threshold than has not been zeroed, then we are done
            elif f1[i] > 0:
                above_thres = False  # stop looping
        prog = (c / len(phots)) * 100
        printProgressBar(prog, 100, prefix='Progress:', suffix='Complete', length=50)

    print("100% complete")
    ep = np.array(ep)
    ex = np.array(ex)
    ey = np.array(ey)
    # fi = open(savePath + '.csv', 'w')
    # fi.write('ep, ex, ey\n')
    # for cntr in range(0, len(ep)):
    #   fi.write(str(ep[cntr]) + ', ' + str(ex[cntr]) + ', ' + str(ey[cntr]) + '\n')
    # fi.close()
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

    # set image size
    nx, ny = 1936, 1096
    # set integration region
    dx, dy = x, y
    biasFrame = avgBias(biasPath)
    darkFrame = avgDark(darkPath)
    darkFrame[darkFrame > t] = 0
    files = glob.glob(filePath + "?.FTS") + glob.glob(filePath + "??.FTS") + glob.glob(
        filePath + "???.FTS") + glob.glob(filePath + "????.FTS") + glob.glob(filePath + "?????.FTS")
    print("Performing analysis of " + str(len(files)) + " photos")
    if autoThres:
        tp = 2.5 * std_2d(biasFrame + darkFrame)
        t = 5 * std_2d(biasFrame + darkFrame)
    print("thres set to " + str(t))
    print("thresp set to " + str(tp))
    print("dx set to " + str(dx))
    print("dy set to " + str(dy))
    hotPixels = getHotPix(darkPath)
    # for ii in range(0, len(hotPixels)):
    # print("x: " + str(hotPixels[ii][1]) + " y: " + str(hotPixels[ii][0]))
    # end len(files)
    for c in range(0, len(files)):
        data = pyfits.open(files[c])
        if (len(data) == 1):
            imData = data[0].data
        else:
            imData = data[1].data
        data.close()
        phot = imData - (biasFrame + darkFrame)
        phot[phot < 0] = 0
        for ii in range(0, len(hotPixels)):
            phot[hotPixels[ii][0]][hotPixels[ii][1]] = 0

        # print(phot)
        nx, ny = phot.shape  # find size of image
        f1 = np.reshape(phot, nx * ny)  # change the shape to be 1d
        q = np.argsort(f1)[::-1]
        # print(f1[q[0]])
        # print(q[0])
        # go through pixels with large signal
        j = 0
        above_thres = True
        while above_thres:
            i = q[j]
            # print(j, q[j], f1[i])
            j += 1
            if f1[i] >= t:  # pixel is above threshold
                # print("i: " + str(i) + " j: " + str(j))
                x = (i % 1936) + 1
                y = math.floor((i / 1936) + 1)
                yR = (i % 1936)
                xR = math.floor((i / 1936))
                # print("Checking pixel " + str(x) + ", " + str(y) + " with value " + str(phot[xR][yR]))
                # i = ny*x + y # index into 1d array
                # check if pixel is far enough from edge of chip
                if (xR > dx) and (xR < nx - dx - 1) and (yR > dy) and (yR < ny - dy - 1):
                    # p = sum(f1[i+di])
                    p = 0
                    for xi in range(xR - dx, xR + dx + 1):
                        for yi in range(yR - dy, yR + dy + 1):
                            if phot[xi, yi] > tp: p += phot[xi, yi]
                            phot[xi, yi] = 0  # zero out pixels in this event
                    # if (p > 106):
                    # print(str(p) + " pulse at " + str(x) + ", " + str(y) + " on image " + str(files[c]) + " - c: " + str(c))
                    # print("pulse at " + str(y) + ", " + str(x) + " - " + str(p))

                    if (p > 0):
                        ep.append(p)
                        ex.append(x)
                        ey.append(y)
                    # print (x, y, p)
                    # f1[i+di] = 0 # zero out pixels in this event
            # if we reach a pixel below threshold than has not been zeroed, then we are done
            elif f1[i] > 0:
                above_thres = False  # stop looping
        prog = (c / len(files)) * 100
        printProgressBar(prog, 100, prefix='Progress:', suffix='Complete', length=50)

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


def fileAnalysisDebug(filePath, biasPath, darkPath, thres, thresp, savePath, x=2, y=2, autoThres=False):
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

    # set image size
    nx, ny = 1936, 1096
    # set integration region
    dx, dy = x, y
    biasFrame = avgBias(biasPath)
    darkFrame = avgDark(darkPath)
    darkFrame[darkFrame > t] = 0
    files = glob.glob(filePath + "?.FTS") + glob.glob(filePath + "??.FTS") + glob.glob(
        filePath + "???.FTS") + glob.glob(filePath + "????.FTS") + glob.glob(filePath + "?????.FTS")
    print("Performing analysis of " + str(len(files)) + " photos")
    if autoThres:
        tp = 2.5 * std_2d(biasFrame + darkFrame)
        t = 5 * std_2d(biasFrame + darkFrame)
    print("thres set to " + str(t))
    print("thresp set to " + str(tp))
    print("dx set to " + str(dx))
    print("dy set to " + str(dy))
    #hotPixels = getHotPix(darkPath)
    # for ii in range(0, len(hotPixels)):
    # print("x: " + str(hotPixels[ii][1]) + " y: " + str(hotPixels[ii][0]))
    # end len(files)
    for c in range(0, len(files)):
        data = pyfits.open(files[c])
        if (len(data) == 1):
            imData = data[0].data
        else:
            imData = data[1].data
        data.close()
        phot = imData - (biasFrame + darkFrame)
        phot[phot < 0] = 0
        ##for ii in range(0, len(hotPixels)):
          ##  phot[hotPixels[ii][0]][hotPixels[ii][1]] = 0

        # print(phot)
        nx, ny = phot.shape  # find size of image
        f1 = np.reshape(phot, nx * ny)  # change the shape to be 1d
        q = np.argsort(f1)[::-1]
        # print(f1[q[0]])
        # print(q[0])
        # go through pixels with large signal
        j = 0
        above_thres = True
        while above_thres:
            i = q[j]
            # print(j, q[j], f1[i])
            j += 1
            if f1[i] >= t:  # pixel is above threshold
                # print("i: " + str(i) + " j: " + str(j))
                x = (i % 1936) + 1
                y = math.floor((i / 1936) + 1)
                # i = ny*x + y # index into 1d array
                # check if pixel is far enough from edge of chip
                if (x > dx) and (x < nx - dx - 1) and (y > dy) and (y < ny - dy - 1):
                    # p = sum(f1[i+di])
                    p = 0
                    for xi in range(x - dx, x + dx + 1):
                        for yi in range(y - dy, y + dy + 1):
                            if phot[xi, yi] > tp: p += phot[xi, yi]
                            phot[xi, yi] = 0  # zero out pixels in this event
                    # if (p > 106):
                    # print(str(p) + " pulse at " + str(x) + ", " + str(y) + " on image " + str(files[c]) + " - c: " + str(c))
                    # print("pulse at " + str(y) + ", " + str(x) + " - " + str(p))

                    ep.append(p)
                    ex.append(y)
                    ey.append(x)
                    # print (x, y, p)
                    # f1[i+di] = 0 # zero out pixels in this event
            # if we reach a pixel below threshold than has not been zeroed, then we are done
            elif f1[i] > 0:
                above_thres = False  # stop looping
        prog = (c / len(files)) * 100
        printProgressBar(prog, 100, prefix='Progress:', suffix='Complete', length=50)

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


def getHotPix(darkPath):
    ps = getFilesForBrute(darkPath, bar=False)
    print("Parsing " + str(len(ps)) + " files for hot pixels")
    hots = []
    pos = []
    for phot in ps:
        nx, ny = phot.shape # find size of image
        f1 = np.reshape(phot, nx*ny) # change the shape to be 1d
        q = np.argsort(f1)[::-1]
        for i in range(0, 50):
            if(f1[q[i]] > 25):
                hots.append(q[i])

    h = np.unique(hots)
    for i in range(0, len(h)):
        y = (h[i] % 1936)
        x = math.floor((h[i] / 1936))
        pos.append([x,y])
    print(str(len(h)) + " pixels found")
    return pos


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
    nx = 1936
    ny = 1096
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
    # end len(files)
    for c in range(0, len(files)):
        data = pyfits.open(files[c])
        imData = data[0].data
        data.close()
        phot = imData - (biasFrame + darkFrame)
        phot[phot < 0] = 0
        # print(phot)
        nx, ny = phot.shape  # find size of image
        f1 = np.reshape(phot, nx * ny)  # change the shape to be 1d
        q = np.argsort(f1)[::-1]
        print(f1[q[0]])
        print(q[0])
        # go through pixels with large signal
        j = 0
        above_thres = True
        while above_thres:
            i = q[j]
            # print(j, q[j], f1[i])
            j += 1
            if f1[i] >= t:  # pixel is above threshold
                print("i: " + str(i) + " j: " + str(j))
                x = int(i / ny)
                y = i % ny
                # i = ny*x + y # index into 1d array
                # check if pixel is far enough from edge of chip
                if (x > dx) and (x < nx - dx - 1) and (y > dy) and (y < ny - dy - 1):
                    # p = sum(f1[i+di])
                    p = 0
                    for xi in range(x - dx, x + dx + 1):
                        for yi in range(y - dy, y + dy + 1):
                            if phot[yi - 1, xi - 1] > tp: p += phot[yi - 1, xi - 1]
                            phot[yi - 1, xi - 1] = 0  # zero out pixels in this event
                    if (p > 490 and p < 515):
                        print(str(p) + " pulse at " + str(x) + ", " + str(y) + " on image " + str(
                            files[c]) + " - C: " + str(c))

                    print("pulse at " + str(x) + ", " + str(y) + " - " + str(p))
                    ep.append(p)
                    ex.append(x)
                    ey.append(y)
                    # print (x, y, p)
                    # f1[i+di] = 0 # zero out pixels in this event
            # if we reach a pixel below threshold than has not been zeroed, then we are done
            elif f1[i] > 0:
                above_thres = False  # stop looping

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
        if (i == 0):
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
        # print(oldFile)
        f.close()
        if (checkForDuplicate):
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
    plt.xlabel('Pulse Height (ADU)')
    plt.ylabel('Intensity')
    plt.hist(ep, bins=b, range=[0, max(unis)])
    if (saveImg):
        plt.savefig(saveName)
    plt.ioff()
    plt.show()


def spectrum(ep, m, b, saveImg=False, saveName="default.png", title=False):
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
    plt.figure('Spectrum', figsize=(18, 10))
    plt.clf()
    lists = sorted(unisList.items())  # sorted by key, return a list of tuples
    x, y = zip(*lists)  # unpack a list of pairs into two tuples
    x = np.asarray(x)
    xdat = aduToeV(x, m, b)
    y = np.asarray(y)
    ydat = smoothing(y, 2)
    # plt.rc('xtick', labelsize=8)
    plt.xlabel("Energy (eV)")
    plt.ylabel("Intensity")
    plt.plot(xdat, ydat)
    #plt.gca().xaxis.set_major_locator(MultipleLocator(200))
    #plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%d'))
    #plt.gca().xaxis.set_minor_locator(MultipleLocator(50))
    #plt.xticks(np.arange(100, max(xdat), 200.0))
    if (saveImg):
        plt.savefig(saveName)
    plt.ioff()
    plt.show()
    if (title):
        label = input('Title: ')
        plt.title(label)
    return (xdat, y)


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
    lists1 = sorted(unisList1.items())  # sorted by key, return a list of tuples
    lists2 = sorted(unisList2.items())
    x1, y1 = zip(*lists1)  # unpack a list of pairs into two tuples
    x2, y2 = zip(*lists2)
    x1 = np.asarray(x1)
    x2 = np.asarray(x2)
    xdat1 = aduToeV(x1, m, b)
    xdat2 = aduToeV(x2, m, b)
    y1 = np.asarray(y1)
    y2 = np.asarray(y2)
    ydat1 = smoothing(y1, 2)
    ydat2 = smoothing(y2, 2)
    plt.xlabel("Energy (eV)")
    plt.ylabel("Intensity")
    plt.plot(xdat1, ydat1, label='Spectrum 1')
    plt.plot(xdat2, ydat2, label='Spectrum 2')
    plt.legend()
    plt.xticks(np.arange(min(xdat1), max(xdat1) + 1, 200.0))
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

def linear_fit(x,y,yerr):
    # Sums we will use to find A,B that minimize chi-squared
    sumy = np.sum(y/(yerr**2))
    sum1 = np.sum(1/(yerr**2))
    sumx = np.sum(x/(yerr**2))
    sumxy = np.sum((x*y)/(yerr**2))
    sumx2 = np.sum((x**2)/(yerr**2))
    delta = ((sum1*sumx2) - (sumx**2))
    # Fit value for A
    A = ((sumx2*sumy) - (sumx*sumxy))*(1/delta)
    # Fit value for B
    B = ((sum1*sumxy) - (sumx*sumy))*(1/delta)
    # Error on A and B
    error_A = np.sqrt((sumx)*(1/delta))
    error_B = np.sqrt((sum1)*(1/delta))
    # Chi_Squared value
    chi_squared = np.sum(((y-A-(B*x))/yerr)**2)
    # In the linear fit we have two parameters
    num_parameters = 2
    # Calculate the num of data points
    num_points = np.size(x)
    # Degrees of freedom
    deg_freedom = num_points - num_parameters
    # Make an array for [A,B]
    ab_array = np.array([A,B])
    # Make an array for error on [A,B]
    error_array = np.array([error_A,error_B])
    return(ab_array,error_array,chi_squared,deg_freedom)

# x_name and y_name will be the name of the axes of the plot


def findadutoev2(ev, evp, err):
    fit = linear_fit(ev, evp, err)
    print(fit)
    intercept = fit[0][0]  # The A value that minimizes Chi-Squared
    slope = fit[0][1]  # The B value that minimizes Chi-Squared
    intercept_Error = fit[1][0]  # Error on A
    slope_Error = fit[1][1]  # Error on B
    Chi_Squared = fit[2]  # Chi-Squared val
    print(slope_Error)  # ue
    DOF = fit[3]  # DOF on fit
    p = 1 - stats.chi2.cdf(Chi_Squared, DOF)  # Probability of Chi-Squared
    print("The best fit parameters are:\nA = {} +/- {}\nB = {} +/- {}".format(intercept, intercept_Error, slope,
                                                                              slope_Error))
    print("The chi-squared value is {}".format(np.round(Chi_Squared, 3)))
    print("There are {} degrees of freedom".format(DOF))
    print("The fit probability is {}".format(np.round(p, 2)))
    fig = plt.figure("Fig1")  # Create a window
    ax = fig.add_subplot(111)
    # Title the plot
    ax.set_title("Energy Calibration")
    # Plot the data points w/ error
    ax.errorbar(ev, evp, err, fmt='o', markersize=3, color='k')
    # Plot both axes with the input names
    ax.set_xlabel('Pulse height (ADU)')
    ax.set_ylabel('Energy (keV')
    textstr = '\n'.join((
        r'Slope = %.2f +- %.2f' % (slope, slope_Error), r'Intercept = %.2f +- %.2f' % (intercept, intercept_Error)))
    # these are matplotlib.patch.Patch properties
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

    # place a text box in upper left in axes coords
    ax.text(.19, .87, textstr, fontsize=10, ha='center', va='center', transform=ax.transAxes, bbox=props)
    # Call linear fit on input data
    x_axis = np.linspace(0, max(ev), np.size(ev))  # X axis of plot
    # Plot the linear model with the best fit parameters
    ax.plot(x_axis, slope * x_axis + intercept, label='Y = {}x + {}'.format(np.round(slope, 5), np.round(intercept, 5)),
            color='r')
    # plt.legend()    # Show legend
    # plt.grid(True)    # Show grid lines
    # plt.figure("Residuals Plot")   # Make new window for residual plot
    # plt.clf()              # Clear window
    # Residuals = (y-(B*x+A))/err      # Calculate the residual values
    # plt.scatter(x_axis,Residuals, color = 'k')  # Scatter plot residuals vs x
    # plt.axhline(0, color='r',linestyle = 'dashed')  # Plot a line at y = 0 for reference
    # plt.xlabel(str(x_name))   # Label x axes
    # plt.ylabel("Residuals")  # Label y axes
    # Title the plot
    # plt.title("Residuals vs {}".format(str(x_name)))
    return(slope,slope_Error,intercept,intercept_Error)

def findADUtoeV(ep, ev, plot):
    """
    Finds the relationship between pulse height and associated energy
    Returns a tuple of m, b, r value, p value, and std of the fit
    Paramaters:
        ep - event peaks intensity list
        ev - corresponding energy to each event (eV)
        plot - boolean: do you want a plot of the relationship or not
    """
    slope, intercept, r_value, p_value, std_err = stats.linregress(ep, ev)
    if (plot):
        plt.ion()
        fig = plt.figure("Energy versus ADU")
        ax = fig.add_subplot(1,1,1)

        graphSpace = np.linspace(0, max(ep), 500)
        ax.plot(ev, ep,marker =  'o',label='Data points')
        ax.plot(((graphSpace * slope) + intercept), graphSpace, label='Best fit')
        ax.set_xlabel("Energy (keV)")
        ax.set_ylabel("Pulse height (ADU)")
        textstr = '\n'.join((
            r'Slope = %.4f (ADU/keV)' % (slope),
            r'Intercept = %.4f ADU' % (intercept),
            r'r value  = %.5f' % (r_value), r'std error = %.4f' % (std_err)))
        # these are matplotlib.patch.Patch properties
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

        # place a text box in upper left in axes coords
        ax.text(.23, .87, textstr, fontsize=10, ha='center', va='center', transform=ax.transAxes, bbox=props)
        ax.plot()
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
    pics = glob.glob(fileName + "?.FTS") + glob.glob(fileName + "??.FTS") + glob.glob(fileName + "???.FTS") + glob.glob(
        fileName + "????.FTS")
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


def printProgressBar(iteration, total, prefix='', suffix='', decimals=1, length=100, fill='█'):
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
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end='\r')
    # Print New Line on Complete
    if iteration == total:
        print()


def func_gaussc(x, norm, cent, width, cons):
    return norm * np.exp(-(x - cent) ** 2 / (2 * width ** 2)) + cons



def index_of(arrval, value):
    """Return index of array *at or below* value."""
    if value < min(arrval):
        return 0
    return max(np.where(arrval <= value)[0])




def gausfitc(data,amp,cent,fwhm):
    m = 0.0037168339915133394
    b = 0.002116795017471418
    bins = int((max(data) - min(data))/2)
    (counts, edges) = np.histogram(data, bins=bins)
    n = len(counts)
    chan = 0.5 * (edges[0:n] + edges[1:n + 1])
    chan = aduToeV(chan, m, b)
    sig = fwhm / 2.35
    cerr = 1 + np.sqrt(counts + 0.75)
    cerr = cerr[:-1]
    chan = chan[:-1]
    counts = counts[:-1]
    gauss = lmfit.Model(func_gaussc,prefix= 'g1_')
    pars = gauss.make_params()
    pars['g1_norm'].set(value=amp,min = 0)
    pars['g1_cent'].set(value=cent, min=min(chan),max = max(chan))
    pars['g1_width'].set(value=sig, min=0)
    pars['g1_cons'].set(value = 0)
    init = gauss.eval(pars, x=chan)
    out = gauss.fit(counts, pars, x=chan,weights = 1/cerr)
    r_chisq = out.redchi
    g1_norm = out.best_values['g1_norm']
    g1_cent = out.best_values['g1_cent']
    g1_sigma = out.best_values['g1_width']
    g1_cons = out.best_values['g1_cons']
    g1_fwhm = 2.35 * g1_sigma
    covar = out.covar
    err_norm = np.sqrt(abs(covar[0][0]))
    err_cent = np.sqrt(abs(covar[1][1]))
    err_sigma = np.sqrt(abs(covar[2][2]))
    err_cons = np.sqrt(abs(covar[3][3]))
    err_fwhm = np.sqrt((err_sigma*2)**2.35)
    height = g1_norm + g1_cons
    err_h = np.sqrt(err_norm**2 + err_cons**2)
    print(out.fit_report(min_correl=0.5))
    fig, axes = plt.subplots(1, figsize=(10, 8))

    textstr = '\n'.join((
        r'Height = %.2f +- %.2f' % (height,err_h),
        r'Centroid = %.2f +- %.2f eV' % (g1_cent,err_cent),
        r'FWHM = %.2f +- %.2f eV' % (g1_fwhm,err_fwhm),
        r'Constant = %.2f +- %.2f eV' % (g1_cons, err_cons),
        r'$\chi^2/dov $= %.2f' % (r_chisq)))
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    axes.text(.22, .89, textstr, fontsize=8, ha='center', va='center', transform=axes.transAxes, bbox=props)
    axes.errorbar(chan, counts, cerr, fmt='.')
    axes.plot(chan, out.best_fit, 'r-', label='best fit')
    axes.legend(loc='best')
    return (height,g1_cent,g1_fwhm,g1_cons,r_chisq)


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

