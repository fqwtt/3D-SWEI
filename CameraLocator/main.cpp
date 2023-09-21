#include "CamControl.h"
#include "Locator.h"
#include <thread>
#include <shared_mutex>

using std::shared_mutex;
using std::mutex;
using std::thread;

shared_mutex sm_pData;
queue<unsigned char*> pDatas;
MV_FRAME_OUT_INFO_EX g_pFrameInfo;
bool ispFrameInfoWritten = false;

shared_mutex sm_USImages;
vector<cv::Mat> USImages; //saved in memory until whole program ends

shared_mutex sm_Rt;
queue<cv::Mat> R_w2cs, t_w2cs;

bool ExitFlag = false;
mutex mut_Exit;

bool LocateExitFlag = false; //After global exit flag is set true, the locating thread needs some time to exit, thus
	//the locating thread needs an independent exit flag (I don't know why, but it works)
mutex mut_LocateExit;

/*****ScreenCapturer Variables*******/
WindowInfo window;
int tl_x, tl_y;
int cap_width, cap_height;
string windowName = "Flash_64";

/******Camera Control Paramaters******/
const float exposureTime = 50000.000; //us
const int gain = 3000;
float fps = 5;

/******Locator Paramaters******/
#define USE_ARUCO 1
#define SET_MODE 1 //0: Set board pattern mannually, 1: reading the board pattern from intrinsFile
string dataPath = "F:\\Users\\jxy\\recon0907\\"; //Path for storing images and output texts
string intrinsFile = dataPath + "CamParams_0804.xml"; //xml or yaml file generated from opencv calibration sample

bool InitScreenCapturer()
{
	/****** GDI+ Init ******/
	Gdiplus::GdiplusStartupInput gdiplusStartupInput;
	ULONG_PTR token;
	Gdiplus::GdiplusStartup(&token, &gdiplusStartupInput, NULL);
	//if (ScreenCapturer::Init(windowName, window, tl_x, tl_y, cap_width, cap_height) == false)
	if (ScreenCapturer::Init(windowName, window, cap_width, cap_height) == false)
	{
		cout << "InitScreenCapturer() failed! Please open the " + windowName + "window first!" << endl;
		return false;
	}
	return true;
}

void Locate(CamController* camcontroller)
{
	string saveRtPath = dataPath + "fileRt.txt";
	if (remove(saveRtPath.c_str()) == 0)
		cout << saveRtPath << "already exists. Previous file has been deleted." << endl;

	/****** Locator Init ******/
#if USE_ARUCO
	//Locator locator = Locator(dataPath, intrinsFile, camcontroller);
	Locator locator = Locator(dataPath, intrinsFile, camcontroller, 0.3, 1); //Not resizing ArUco
#else
	int segWidth = 250;
	int segHeight = 250;
	float minArea = 3000;
	float maxArea = 30000;
	Locator locator = Locator(dataPath, intrinsFile, camcontroller, segWidth, segHeight, minArea, maxArea);
#endif // USE_ARUCO

#if SET_MODE
	if (!locator.init())
	{
		camcontroller->CloseCam();
		throw(-1);
	}
#else //Set board pattern mannually, instead of reading from intrinsFile. 
	int patternHeight = 11;
	int patternWidth = 4;
	float intervalSize_width = 5;
	float intervalSize_height = 5;
	if (!locator.init(patternHeight, patternWidth, intervalSize_width, intervalSize_height))
	{
		camcontroller->CloseCam();
		throw(-1);
	}
#endif

	locator.run();
}

int main()
{
	if (!InitScreenCapturer())
		return -1;

	//Let the camera run. It uses callback mechanism, thus won't block the following codes.
	CamController camcontroller(exposureTime, gain, fps);
	if (!camcontroller.StartCam())
		return -1;
#if 1
	thread t_locate(Locate, &camcontroller);
	t_locate.join();

#endif

	saveImages(dataPath + "USIMG\\", "USImage");

	return 0;
}