import glob
import astropy.io.fits as pyfits
import numpy as np
import os
import math

badShapes = [4112, 2048, 528384, 1052672, 8392704, 4097, 4098, 4100, 4104, 4120, 4128, 5120, 6148, 4256, 4608, 4736, 5152, 13312, 20480, 20992, 36864, 135172, 135174, 25169920, 16781312, 8394752, 135176, 405505, 530432, 544768, 659456, 1054912, 1060864, 2101248, 2103488, 2109440, 2109824, 4198528, 4200448, 8394944, 8400896, 4102, 4108, 4225, 5248, 6145, 6156, 6160, 6656, 6672, 12289, 12291, 12292, 12294, 12296, 12321, 12448, 12673, 12704, 13568, 20608, 22656, 22720, 36992, 37888, 45056, 45440, 135169, 135680, 136192, 136224, 202756, 202760, 203264, 398368, 405508, 405512, 405536, 406528, 661504, 1054720, 1085440, 2103296, 2109568, 3152000, 4198400, 4200640, 4206592, 4206720, 6295552, 8392832, 8401280, 12587008, 16781440, 16783424, 16789888, 16912384, 16914432, 17305600, 176128, 4099, 4240, 4624, 6152, 6288, 46080, 135704, 137217, 143392, 202758, 202768, 405506, 405516, 405520, 1052800, 3158400, 12587136, 12595200, 16789504, 135180, 398336, 2103424, 1183744, 1191936, 1454080, 4129, 4481, 6800, 46208, 135696, 137728, 143366, 202753, 3149824, 16920576, 25170048, 12290, 4131, 4632, 4800, 5280, 6146, 6147, 6150, 6168, 6273, 6352, 6664, 6848, 12304, 12320, 12417, 12545, 12705, 13344, 13472, 13696, 21120, 22528, 38016, 46464, 72192, 135170, 135171, 135192, 135200, 137220, 137224, 137240, 143362, 143376, 143384, 144384, 144416,151552, 152064, 153600, 154112, 167936, 168960, 177152, 202754, 202776, 219136, 219648, 282624, 397315, 397316, 397320, 405507, 439296, 528512, 530624, 544896, 546816, 675840, 727040, 743424, 1060992, 1061248, 1185792, 2101376, 2101440, 3182592, 3152064, 4200512, 6297600, 6297792, 6303744, 6304128,8401024, 8460288, 12589056, 12589184, 16783360, 16783552, 16789632, 17307648, 17321984, 25178112, 45184]


def makeDarkFrame(darkpath):
	files = glob.glob(darkpath + "?.FTS") + glob.glob(darkpath + "??.FTS") + glob.glob(darkpath + "???.FTS")
	if (len(files) == 0):
		print("No files found at path!")
		return None

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
	print("Dark frame compiled from " + str(len(files)) + " frames")
	pyfits.writeto(darkpath + 'avgDark.FTS', medBias, None, overwrite=True)
	return medBias

def makeBiasFrame(biaspath):
	files = glob.glob(biaspath + "?.FTS") + glob.glob(biaspath + "??.FTS") + glob.glob(biaspath + "???.FTS")
	if (len(files) == 0):
		print("No files found at path!")
		return None

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
	print("Bias frame compiled from " + str(len(files)) + " frames")
	pyfits.writeto(biaspath + 'avgBias.FTS', medBias, None, overwrite=True)
	return medBias

def getImg(path):
	dat = pyfits.open(path)
	if (len(dat) == 1):
		index = 0
	else:
		index = 1
	img = dat[index].data
	dat.close()
	return img

def getHotPix(darkPath, thr=9, rep=5, argonne=False):
	lpath = darkPath + "hotT" + str(thr) + "R" + str(rep) + ".csv"
	if(os.path.exists(lpath)):
		print("Hot pixel list already compiled")
		(x, y) = np.loadtxt(lpath, delimiter=',', unpack=True)
		pos = []
		for i in range(0, len(x)):
			pos.append([x[i], y[i]])
		return pos

	files = glob.glob(darkPath + "?.FTS") + glob.glob(darkPath + "??.FTS") + glob.glob(darkPath + "???.FTS")
	print("Parsing " + str(len(files)) + " files for hot pixels bounded by " + str(thr) + " and multiplicity " + str(rep))
	hots = []
	pos = []
	for f in files:
		phot = getImg(f)
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
	print(str(len(pos)) + " pixels found")
	print("Loss: %.2f" % ((len(pos)/(1096 * 1936))*100))

	f = open(lpath, 'w')
	for i in range(0, len(pos)):
		f.write(str(pos[i][0]) + ", " + str(pos[i][1]) + "\n")
	f.close()
	return pos

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
