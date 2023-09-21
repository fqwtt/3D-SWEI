#pragma once
#include <stdio.h>
#include <Windows.h>
#include <conio.h>
#include "MvCameraControl.h"
#include <opencv2/opencv.hpp>
#include <queue>
#include <thread>
#include <shared_mutex>
#include <vector>
#include "ScreenCapturer.h"

using std::vector;
using std::queue;
using std::shared_mutex;

void __stdcall ImageCallBackEx(unsigned char* pData, MV_FRAME_OUT_INFO_EX* pFrameInfo, void* pUser);
bool Convert2Mat(MV_FRAME_OUT_INFO_EX* pstImageInfo, unsigned char* pData, cv::Mat& mat);

class CamController {
public:
    CamController(const float& exposureTime, const int& gain, const float& fps);

    bool StartCam();

    void CloseCam();

private:
	bool PrintInterfaceInfo(MV_GENTL_IF_INFO* pInterfaceInfo);

	bool PrintDeviceInfo(MV_GENTL_DEV_INFO* pDeviceInfo);

private:
    int nRet = MV_OK;
    void* handle = NULL;
    MV_GENTL_IF_INFO_LIST myInterfaceList;
    MV_GENTL_IF_INFO* pInterfaceInfo = NULL;
    MV_GENTL_DEV_INFO_LIST myDeviceList;
    MV_GENTL_DEV_INFO* pDeviceInfo = NULL;

    const float _exposureTime; //us
    const int _gain;
    const float _fps;
};







