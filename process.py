import argparse
import astropy.io.fits as pyfits
import os
import CMOSUtilsLive as cu
import numpy as np
import math

forceHot = [[261, 824], [249, 295], [1080, 1713], [1049, 605], [341, 290], [396, 830], [601, 820], [1051, 1852], [584, 909], [243, 688], [441, 1323], [731, 1371], [513, 1309], [9, 1132], [961, 984], [651, 310], [6, 850]]

parser = argparse.ArgumentParser()
parser.add_argument('fileName', nargs=1)
parser.add_argument('overWrite', nargs=1)
parser.add_argument('outPath', nargs=1)
parser.add_argument('dName', nargs=1)
parser.add_argument('bName', nargs=1)
args = parser.parse_args()

#print(str(args.fileName[0]))
#print(str(args.overWrite[0]))
#print(str(args.dName[0]))
#print(str(args.bName[0]))
overWrite = False
if(str(args.overWrite[0]) == "True"):
	overWrite = True


darkPath = str(args.dName[0])
biasPath = str(args.bName[0])

if(overWrite):
	if (os.path.exists(darkPath + "avgDark.FTS")):
		print("Dark frame found")
		ds = cu.getImg(darkPath + "avgDark.FTS")
	else:
		print("Dark frame not found")
		ds = cu.makeDarkFrame(darkPath)
	if (os.path.exists(biasPath + "avgBias.FTS")):
		print("Bias frame found")
		bs = cu.getImg(biasPath + "avgBias.FTS")
	else:
		print("Bias frame not found")
		bs = cu.makeBiasFrame(biasPath)
	if (type(ds) == type(None) or type(bs) == type(None)):
		print("You need to take the proper bias and dark frames before using the live camera function")
else:
	ds = cu.getImg(darkPath + "avgDark.FTS")
	bs = cu.getImg(biasPath + "avgBias.FTS")

img = cu.getImg(str(args.fileName[0]))
t = 20
tp = 10
dx, dy = 2, 2
ep, ex, ey = [], [], []
shp = []
nx, ny = 1936, 1096
hotPixels = cu.getHotPix(darkPath, thr=8, rep=10, argonne=False)

img = img - ds - bs
img[img < 0] = 0

for ii in range(0, len(hotPixels)):
	img[int(hotPixels[ii][0])][int(hotPixels[ii][1])] = 0
for ii in range(0, len(forceHot)):
	img[forceHot[ii][0]][forceHot[ii][1]] = 0

f1 = np.reshape(img, nx*ny)
q = np.argsort(f1)[::-1]
j = 0
above_thres = True
while above_thres:
	i = q[j]
	j += 1
	if f1[i] >= t:
		x = (i % 1936) + 1
		y = math.floor((i / 1936) + 1)
		xR = math.floor(i / 1936)
		yR = (i % 1936)
		if (xR > dx) and (xR < nx-dx-1) and (yR > dy) and (yR < ny-dy-1):
			p = 0
			area = img[(xR-dx):(xR+dx+1), (yR-dy):(yR+dy+1)]
			pAdjust, s = cu.pixelAnalysis(area, tp)
			for xi in range(xR-dx, xR+dx+1):
				for yi in range(yR-dy, yR+dy+1):
					if(img[xi, yi] > tp):
						p += img[xi, yi]
						if (pAdjust == 0):
							img[xi, yi] = 0
						elif(xi-xR+dx != 0 and xi-xR+dx != 4 and yi-yR+dy != 0 and yi-yR+dy != 4):
							img[xi, yi] = 0
			if (p > 0):
				ep.append(p - pAdjust)
				ex.append(x)
				ey.append(y)
				shp.append(s)
		elif f1[i] > 0:
			above_thres = False
print(ep)
if(overWrite):
	f = open(str(args.outPath[0]), 'w')
	f.write('ep, ex, ey\n')
	ep = np.array(ep)
	ex = np.array(ex)
	ey = np.array(ey)
	shp = np.array(shp)
	for i in range(0, len(ep)):
		f.write(str(ep[i]) + ", " + str(ex[i]) + ", " + str(ey[i]) + ":" + str(shp[i]) + "\n")
	f.close()
else:
	f = open(str(args.outPath[0]), 'a')
	ep = np.array(ep)
	ex = np.array(ex)
	ey = np.array(ey)
	shp = np.array(shp)
	for i in range(0, len(ep)):
		f.write(str(ep[i]) + ", " + str(ex[i]) + ", " + str(ey[i]) + ":" + str(shp[i]) + "\n")
	f.close()
