import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import style
import time
import numpy as np
import astropy.io.fits as pyfits
from scipy.stats import gaussian_kde
import sys
from matplotlib import cm


def aduToeV(adu, m, b):
    return ((m * adu) + b) * 1000


def getImg(path):
    dat = pyfits.open(path)
    if (len(dat) == 1):
        index = 0
    else:
        index = 1
    img = dat[index].data
    dat.close()
    return img


# No spaces, no string
filepath = input('CSV Filepath: ')
show_images = input('Show Images? ')

try:
    if int(show_images) == 1:
        image_path = input('Image Filepath: ')
except ValueError:
    print('Enter 0 for no images or 1 for images')
    print('Exiting Script')
    sys.exit()

fig = plt.figure('Live Plotting', figsize=(16, 8))
ax1 = fig.add_subplot(2, 2, 1)
ax2 = fig.add_subplot(2, 2, 2)
ax3 = fig.add_subplot(2, 2, 3)
ax4 = fig.add_subplot(2, 2, 4)


def animate_intensity(i):
    """Live Plot of Energy Histogram
    """
    # Energy calibration
    m = 0.0037168339915133394
    b = 0.002116795017471418
    # Open the CSV file
    graph_data = open(filepath, 'r').read()
    lines = graph_data.split('\n')
    ep_list = []
    ex_list = []
    ey_list = []
    # Loop through each line of CSV
    for line in lines:
        # Blocks bug when last line has 0 length
        if len(line) > 1:
            image, ep, ep_x, ep_y = line.split(',')
            # Ignore the first line that says ep,ex,ey....
            try:
                ep_x = float(ep_x)
                ep_y = float(ep_y)
                ep = float(ep)
                ep_list.append(ep)
                ex_list.append(ep_x)
                ey_list.append(ep_y)
            except ValueError:
                print('skipping iteration')
                pass
    # Create the full path to the image
    ax1.clear();
    ax2.clear();
    ax3.clear();
    ax4.clear()
    # Create spectrum for subplot 1, but only works after first iteration
    if len(ep_list) > 1:
        bins = int((max(ep_list) - min(ep_list)) / 2)
        (counts, edges) = np.histogram(ep_list, bins=bins)
        n = len(counts)
        chan = 0.5 * (edges[0:n] + edges[1:n + 1])
        chan = aduToeV(chan, m, b)
        cerr = 1 + np.sqrt(counts + 0.75)
        ax1.errorbar(chan, counts, cerr, fmt='.')
    ax1.set_title('Energy vs Counts')
    ax1.set_xlabel('Energy eV')
    ax1.set_ylabel('Counts/Bin')

    ep_arr = np.array(ep_list)
    ex_arr = np.array(ex_list)
    ey_arr = np.array(ey_list)
    idx = ep_arr.argsort()[::-1]
    sort_ep = ep_arr[idx]
    sort_ex = ex_arr[idx]
    sort_ey = ey_arr[idx]
    ax2.set_xlim([0, 1936])
    ax2.set_ylim([0, 1096])
    d = ax2.scatter(sort_ex, sort_ey, c=aduToeV(sort_ep,m,b), s=25, cmap=cm.jet)
    fig.colorbar(d,ax=ax2)
    if int(show_images) == 1:
        img = getImg(image_path + '/' + image)
        ax3.imshow(img, cmap='gray')
        ax3.set_xlabel('{}'.format(image))
    if int(show_images) != 1:
        ax3.hist(ep_list, bins='auto')

    ax4.set_xlabel('X pixels')
    ax4.set_ylabel('Y pixels')
    ax4.set_title('Event Location Density')
    ax4.set_xlim([0, 1936])
    ax4.set_ylim([0, 1096])
    # Create density plot of events
    xy = np.vstack([ex_list, ey_list])
    z = gaussian_kde(xy)(xy)
    ax4.scatter(ex_list, ey_list, c=z, s=50, edgecolor='')
    ax4.grid(True)
    return ()


# This function does all the live plotting by repeatedly calling animate_intensity based on an interval (ms) you choose

ani = animation.FuncAnimation(fig, animate_intensity, interval=500)
plt.show()

# Testing workflows


# for i in range(len(ep)):
#     f = open('test_file.csv', 'a+')
#     f.write('{}, {}, {}\n'.format(1,2,ep[i]))
#     print('{}, {}, {}\n'.format(1,2,ep[i]))
#     f.close()
#     time.sleep(.7)


# for i in range(len(epc)):
#     f = open('test_file.csv', 'a+')
#     f.write('Neg20DegTest{}.FTS,{}, {}, {}\n'.format(i,epc[i],np.random.randint(0,1936),np.random.randint(0,1096)))
#     print('Neg20DegTest{}.FTS,{}, {}, {}\n'.format(i, epc[i], np.random.randint(0, 1936), np.random.randint(0, 1096)))
#     f.close()
#     time.sleep(1.5)
