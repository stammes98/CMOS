# CMOS
To use CMOSUtils, include CMOSUtils.h into your main.cpp file. Make sure it's the last import
To use the camera function, call the function doLiveCamera();
It takes the following paramaters
int nPhotos - how many photos you want to take
double exposure - exposure of the photo in ms
string imgPath - path to save the image, if saveImages is true.
bool saveImages - if true, the image will be saved to imgPath along with temp.FTS (which will be overwritten by the next photo)
string darkPath - if you're doing live processing, this will be the path to your dark frames (ex /cImages/Dark/DarkFrame)
string biasPath - same as above but for the bias frames
string outPath - path to the .CSV where each processed frame's ep, ex and ey will be written to
bool liveProccess - determines if you will process each image after taking it.

Software flow:
Program init -> doLiveCamera init -> take photo -> save photo to temp.FTS -> save photo to imgPath if saveImages is true -> call process.py if liveProccess is true -> update event list at outPath -> take photo -> repeat

To run the program:
simply type "make" and then "make run"
If you encounter any errors in the compiling of the program, type "make clean" and then try again
If you made any updates to the .cpp or .h files, type "make" again for the program to reflect those changes. This is not needed for making edits to the python files

All the constants (t, tp, etc) are defined in process.py rather than CMOSUtilsLive.py
