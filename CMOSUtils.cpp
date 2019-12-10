#include <stdio.h>
#include <stdlib.h>
#include <ueye.h>
#include <fitsio.h>
#include <string>
#include "CMOSUtils.h"

void camStop(HIDS hCam, char* message);
void doSleep(int ms);

using namespace std;

int doLiveCamera(int nPhotos, double exposure, string imgPath, bool saveImages, string darkPath, string biasPath, string outPath, bool liveProccess) {
	/*
		for(int i = 0; i < nPhotos; i++) {
			//printf("%d\n", i);
			string args = darkPath;
			args = args + string(" ") + biasPath;
			string command = "python process.py ";
			if(i == 0) {
				command = command  + imgPath + string(" True ") + outPath + " " + args;
			} else {
				command = command + imgPath +  string(" False ") + outPath + " " +  args;
			}
			//printf("%s\n", command.c_str());
			system(command.c_str());
		}
		*/

		HIDS hCam = 0;
		printf("BEGIN CAMERA ROUTINE\n");
		char* pcImgMem;
		int memID;
		SENSORINFO sInfo;
		CAMINFO cInfo;
		IS_RECT rectAOI;
		INT nBitsPerPixel = 8;
		INT m_nColorMode;
		int bytesPerPixel = (int)nBitsPerPixel / 8;
		int gain = 30;
		int nRet = is_InitCamera(&hCam, NULL);
		if (nRet == IS_SUCCESS) {
			printf("Camera init.\n");
		} else {
			camStop(hCam, "General camera init error");
			return nRet;
		}
		nRet = is_GetCameraInfo(hCam, &cInfo);
		if (nRet == IS_SUCCESS) {
			printf("Camera info acquired.\n");
		} else {
			camStop(hCam, "Error acquiring camera info.");
			return nRet;
		}
		nRet = is_GetSensorInfo(hCam, &sInfo);
		if (nRet == IS_SUCCESS) {
			printf("Sensor info acquired.\n");
		} else {
			camStop(hCam, "Error acquiring sensor info.");
			return nRet;
		}

		nRet = is_ResetToDefault(hCam);
		if (nRet == IS_SUCCESS) {
			printf("Settings cleared and ready to be assigned.\n");
		} else {
			camStop(hCam, "Error cleaning settings");
			return nRet;
		}
		nRet = is_SetDisplayMode(hCam, IS_SET_DM_DIB);
		if (nRet == IS_SUCCESS) {
			printf("Display mode set.\n");
		} else {
			camStop(hCam, "Error setting display mode.");
			return nRet;
		}
		if ((int)sInfo.nColorMode != IS_COLORMODE_MONOCHROME) {
			camStop(hCam, "Invalid sensor colormode detected.");
			return 1;
		} else {
			m_nColorMode = IS_CM_MONO12;
			bytesPerPixel = 2;
			nBitsPerPixel = 12;
			printf("IS_COLORMODE_MONOCHROME:\n");
			printf("\tm_nColorMode: \t\t%i\n", m_nColorMode);
			printf("\tnBitsPerPixel: \t\t%i\n", nBitsPerPixel);
			printf("\tBytesPerPixel: \t\t%i\n", bytesPerPixel);
		}
		nRet = is_AOI(hCam, IS_AOI_IMAGE_GET_AOI, &rectAOI, sizeof(rectAOI));
		if (nRet == IS_SUCCESS) {
			printf("AOI settings retrieved.\n");
		} else {
			camStop(hCam, "AOI settings error.");
			return nRet;
		}
		int width = rectAOI.s32Width;
		int height = rectAOI.s32Height;

		printf("Camera model: \t\t%s\n", sInfo.strSensorName);
		printf("Camera serial no: \t%s\n", cInfo.SerNo);
		printf("Max image width: \t%i\n", width);
		printf("Max image height: \t%i\n", height);

		nRet = is_AllocImageMem(hCam, width, height, nBitsPerPixel, &pcImgMem, &memID);
		if (nRet == IS_SUCCESS) {
			printf("Memory successfully alloc'd.\n");
			nRet = is_SetImageMem(hCam, pcImgMem, memID);
			if (nRet == IS_SUCCESS) {
				printf("Memory successfully set.\n");
				nRet = is_SetColorMode(hCam, m_nColorMode);
				if (nRet == IS_SUCCESS) {
					printf("Colormode applied.\n");
				} else {
					camStop(hCam, "Error setting colormode.");
					return nRet;
				}
			} else {
				camStop(hCam, "is_SetImageMem error");
				return nRet;
			}
		} else {
			camStop(hCam, "is_AllocImageMem error");
			return nRet;
		}

		UINT ledOff = 3;
		nRet = is_IO(hCam, IS_IO_CMD_LED_SET_STATE, &ledOff, sizeof(ledOff));
		if (nRet == IS_SUCCESS) {
			printf("LED disabled.\n");
		} else {
			camStop(hCam, "Error disabling LED.");
			return nRet;
		}
		printf("\nCamera prep finished. Applying settings.\n");
		int pixelClock = 20;
		nRet = is_PixelClock(hCam, IS_PIXELCLOCK_CMD_SET, &pixelClock, sizeof(pixelClock));
		if (nRet == IS_SUCCESS) {
			printf("Pixel clock set to %i\n", pixelClock);
		} else {
			camStop(hCam, "Error setting pixel clock.");
			return nRet;
		}
		UINT longExpo = 1;
		nRet = is_Exposure(hCam, IS_EXPOSURE_CMD_SET_LONG_EXPOSURE_ENABLE, &longExpo, sizeof(longExpo));
		if (nRet == IS_SUCCESS) {
			printf("Long exposures enabled\n");
		} else {
			camStop(hCam, "Error setting long exposure.");
			return nRet;
		}

		double FPS = 1;
		nRet = is_SetFrameRate(hCam, FPS, &FPS);
		if (nRet == IS_SUCCESS) {
			printf("FPS set to %f\n", FPS);
		} else {
			camStop(hCam, "Error setting FPS.");
			return nRet;
		}

		nRet = is_SetHardwareGain(hCam, gain, gain, gain, gain);
		if (nRet == IS_SUCCESS) {
			printf("Gain set to %i\n", gain);
		} else {
			camStop(hCam, "Error setting gain.");
			return nRet;
		}
		nRet = is_Exposure(hCam, IS_EXPOSURE_CMD_SET_EXPOSURE, &exposure, sizeof(exposure));
		if (nRet == IS_SUCCESS) {
			printf("Exposure set to %f\n", exposure);
		} else {
			camStop(hCam, "Error setting exposure");
			return nRet;
		}
		printf("Settings applied.\n\n");
		nRet = 0;
		for (int i = 0; i < nPhotos + 1; ++i) {
			nRet = is_FreezeVideo(hCam, IS_DONT_WAIT);
			if (i == 0) {
				doSleep(exposure);
				continue;
			} else if (nRet == 140) {
				doSleep(exposure);
				i--;
				continue;
			} else if (nRet == IS_SUCCESS) {
				printf("Image %i taken\n", i);
				VOID *pMem;
				is_GetImageMem(hCam, &pMem);
				if(saveImages) {
					string fName = outPath + to_string(i) + string(".FTS");
					fitsfile *fptr;
					int status = 0;
					long fpixel = 1, naxis = 2, nelements;
					long naxes[2] = { 1936, 1096 };
					fits_create_file(&fptr, fName.c_str(), &status);
					fits_set_compression_type(fptr, RICE_1, &status);
					fits_create_img(fptr, SHORT_IMG, naxis, naxes, &status);
					nelements = naxes[0] * naxes[1];
					fits_write_date(fptr, &status);
					fits_write_img(fptr, TSHORT, fpixel, nelements, pMem, &status);
					fits_close_file(fptr, &status);
					fits_report_error(stderr, status);
				}
				fitsfile *fptr;
				int status = 0;
				long fpixel = 1, naxis = 2, nelements;
				long naxes[2] = { 1936, 1096 };
				fits_create_file(&fptr, "!temp.FTS", &status);
				fits_set_compression_type(fptr, RICE_1, &status);
				fits_create_img(fptr, SHORT_IMG, naxis, naxes, &status);
				nelements = naxes[0] * naxes[1];
				fits_write_date(fptr, &status);
				fits_write_img(fptr, TSHORT, fpixel, nelements, pMem, &status);
				fits_close_file(fptr, &status);
				fits_report_error(stderr, status);
				//doProcess(tp, t);
				if (i == 1) {
					system("python process.py temp.FTS True");
				} else {
					system("python process.py temp.FTS False");
				}
				doSleep(exposure);
			} else {
				printf("Error %d\n", nRet);
				is_FreeImageMem(hCam, pcImgMem, memID);
				camStop(hCam, "Error taking image.");
				return nRet;
			}
		}
		camStop(hCam, "Camera finished.");

		return 0;
}


void camStop(HIDS hCam, char* message) {
	is_ExitCamera(hCam);
	printf(message);
	printf("\n");
	printf("Camera stopped.\n");
}

void doSleep(int ms) {

}
