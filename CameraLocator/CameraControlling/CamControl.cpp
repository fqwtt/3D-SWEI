#include "CamControl.h"

extern shared_mutex sm_pData;
extern queue<unsigned char*> pDatas;
extern MV_FRAME_OUT_INFO_EX g_pFrameInfo;
extern bool ispFrameInfoWritten;

extern WindowInfo window;
extern int tl_x, tl_y, cap_width, cap_height;
extern shared_mutex sm_USImages;
extern vector<cv::Mat> USImages;

int frameCount = 1;

CamController::CamController(const float& exposureTime, const int& gain, const float& fps) :
	_exposureTime(exposureTime), _gain(gain), _fps(fps) {}

bool CamController::StartCam()
{
	bool isStart = false;
	do
	{
		// Enum interfaces
		memset(&myInterfaceList, 0, sizeof(MV_GENTL_IF_INFO_LIST));
		const char* GenTLPath = "C:/Program Files (x86)/Common Files/MVS/Runtime/Win64_x64/MvFGProducerCXP.cti";
		nRet = MV_CC_EnumInterfacesByGenTL(&myInterfaceList, GenTLPath);
		if (MV_OK != nRet)
		{
			printf("Enum Interface fail! nRet [0x%x]\n", nRet);
			break;
		}

		if (myInterfaceList.nInterfaceNum > 0)
		{
			for (unsigned int i = 0; i < myInterfaceList.nInterfaceNum; i++)
			{
				printf("[interface %d]:\n", i);
				pInterfaceInfo = myInterfaceList.pIFInfo[i];
				if (NULL == pInterfaceInfo)
				{
					break;
				}
				PrintInterfaceInfo(pInterfaceInfo);
			}
		}
		else
		{
			printf("\nFind No Interfaces!\n");
			break;
		}

		// Enum devices
		memset(&myDeviceList, 0, sizeof(MV_GENTL_DEV_INFO_LIST));
		nRet = MV_CC_EnumDevicesByGenTL(pInterfaceInfo, &myDeviceList);
		if (MV_OK != nRet)
		{
			printf("\nEnum device fail! nRet [0x%x]\n", nRet);
			break;
		}

		if (myDeviceList.nDeviceNum > 0)
		{
			for (unsigned int i = 0; i < myDeviceList.nDeviceNum; i++)
			{
				printf("\n[device %d]:\n", i);
				pDeviceInfo = myDeviceList.pDeviceInfo[i];
				if (NULL == pDeviceInfo)
				{
					break;
				}
				PrintDeviceInfo(pDeviceInfo);
			}
		}
		else
		{
			printf("\nFind No Devices!\n");
			break;
		}

		printf("\nPlease Input camera index:");
		unsigned int nIndex = 0;
		scanf_s("%d", &nIndex);

		if (nIndex >= myDeviceList.nDeviceNum)
		{
			printf("Input error!\n");
			break;
		}

		// Select device and create handle
		nRet = MV_CC_CreateHandleByGenTL(&handle, myDeviceList.pDeviceInfo[nIndex]);
		if (MV_OK != nRet)
		{
			printf("Create Handle fail! nRet [0x%x]\n", nRet);
			break;
		}

		// Open the specified device
		nRet = MV_CC_OpenDevice(handle);
		if (MV_OK != nRet)
		{
			printf("Open Device fail! nRet [0x%x]\n", nRet);
			break;
		}

		// Set trigger mode as off
		nRet = MV_CC_SetEnumValue(handle, "TriggerMode", MV_TRIGGER_MODE_OFF);
		if (MV_OK != nRet)
		{
			printf("Set Trigger Mode fail! nRet [0x%x]\n", nRet);
			break;
		}

		//Set the gain value
		nRet = MV_CC_SetEnumValue(handle, "PreampGain", _gain);
		if (MV_OK != nRet)
		{
			printf("Set Gain fail! nRet [0x%x]\n", nRet);
			break;
		}

		//Set exposure time
		nRet = MV_CC_SetFloatValue(handle, "ExposureTime", _exposureTime);
		if (MV_OK != nRet)
		{
			printf("Set Exposure Time fail! nRet [0x%x]\n", nRet);
			break;
		}

		//Set frame rate
		nRet = MV_CC_SetFloatValue(handle, "AcquisitionFrameRate", _fps);
		if (MV_OK != nRet)
		{
			printf("Set Acquisition Frame Rate fail! nRet [0x%x]\n", nRet);
			break;
		}

		// Register image callback
		nRet = MV_CC_RegisterImageCallBackEx(handle, ImageCallBackEx, handle);
		if (MV_OK != nRet)
		{
			printf("Register Image CallBack fail! nRet [0x%x]\n", nRet);
			break;
		}

		// Start grab image
		nRet = MV_CC_StartGrabbing(handle);
		if (MV_OK != nRet)
		{
			printf("Start Grabbing fail! nRet [0x%x]\n", nRet);
			break;
		}
		printf("\nStartCam succeeds!\n");
		isStart = true;
	} while (0);

	if(!isStart)
		printf("\nStartCam fails!\n");

	return isStart;
}

void CamController::CloseCam()
{
	do {
		// Stop grab image
		nRet = MV_CC_StopGrabbing(handle);
		if (MV_OK != nRet)
		{
			printf("Stop Grabbing fail! nRet [0x%x]\n", nRet);
			break;
		}
		else
			printf("Stop Grabbing succeeds!\n");

		//取消注册后，回调函数的内存空间会被清空，如果最后一次回调函数还没有执行完，就会产生异常，这里等待一会儿
		Sleep(500);
		// Unregister image callback
		nRet = MV_CC_RegisterImageCallBackEx(handle, NULL, NULL);
		if (MV_OK != nRet)
		{
			printf("Unregister Image CallBack fail! nRet [0x%x]\n", nRet);
			break;
		}

		// Close device
		nRet = MV_CC_CloseDevice(handle);
		if (MV_OK != nRet)
		{
			printf("Close Device fail! nRet [0x%x]\n", nRet);
			break;
		}

		// Destroy handle
		nRet = MV_CC_DestroyHandle(handle);
		if (MV_OK != nRet)
		{
			printf("Destroy Handle fail! nRet [0x%x]\n", nRet);
			break;
		}
	} while (0);

	if (nRet != MV_OK)
	{
		if (handle != NULL)
		{
			MV_CC_DestroyHandle(handle);
			handle = NULL;
		}
	}
	printf("\nCloseCam succeeds!\n");
}

void __stdcall ImageCallBackEx(unsigned char* pData, MV_FRAME_OUT_INFO_EX* pFrameInfo, void* pUser)
{
	//Aquire correspond US image
	//cv::Mat USImage = ScreenCapturer::Capture(window, tl_x, tl_y, cap_width, cap_height);
	cv::Mat USImage = ScreenCapturer::Capture(window, cap_width, cap_height);

	//Push to USImages queue
	sm_USImages.lock();
	USImages.push_back(USImage);
	sm_USImages.unlock();

	//必须先push超声图，再push相机图，不然可能定位结束还没有push进去超声图，导致重建失败
	
	//Acquired one frame
	unsigned char* new_pData = (unsigned char*)calloc(pFrameInfo->nHeight * pFrameInfo->nWidth, sizeof(unsigned char));
	memcpy_s(new_pData, pFrameInfo->nHeight * pFrameInfo->nWidth * sizeof(unsigned char), 
		pData, pFrameInfo->nHeight * pFrameInfo->nWidth * sizeof(unsigned char));

	sm_pData.lock();
	pDatas.push(new_pData);
	if (!ispFrameInfoWritten)
	{
		g_pFrameInfo = *pFrameInfo;
		ispFrameInfoWritten = true;
	}
	sm_pData.unlock();
}

bool CamController::PrintInterfaceInfo(MV_GENTL_IF_INFO* pInterfaceInfo)
{
	if (NULL == pInterfaceInfo)
	{
		printf("The Pointer of pInterfaceInfo is NULL!\n");
		return false;
	}
	else
	{
		printf("chInterfaceID: [%s]\n", pInterfaceInfo->chInterfaceID);
		printf("chTLType: [%s]\n", pInterfaceInfo->chTLType);
		printf("chDisplayName: [%s]\n", pInterfaceInfo->chDisplayName);
	}

	return true;
}

bool CamController::PrintDeviceInfo(MV_GENTL_DEV_INFO* pDeviceInfo)
{
	if (NULL == pDeviceInfo)
	{
		printf("The Pointer of pDeviceInfo is NULL!\n");
		return false;
	}
	else
	{
		printf("chInterfaceID: [%s]\n", pDeviceInfo->chInterfaceID);
		printf("chDeviceID: [%s]\n", pDeviceInfo->chDeviceID);
		printf("chVendorName: [%s]\n", pDeviceInfo->chVendorName);
		printf("chModelName: [%s]\n", pDeviceInfo->chModelName);
		printf("chTLType: [%s]\n", pDeviceInfo->chTLType);
		printf("chDisplayName: [%s]\n", pDeviceInfo->chDisplayName);
		printf("chSerialNumber: [%s]\n", pDeviceInfo->chSerialNumber);
		printf("chDeviceVersion: [%s]\n", pDeviceInfo->chDeviceVersion);
	}

	return true;
}

// convert data stream in Mat format
bool Convert2Mat(MV_FRAME_OUT_INFO_EX* pstImageInfo, unsigned char* pData, cv::Mat& mat)
{
	if (pstImageInfo->enPixelType == PixelType_Gvsp_Mono8)
	{
		mat = cv::Mat(pstImageInfo->nHeight, pstImageInfo->nWidth, CV_8UC1, pData).clone();
		return true;
	}
	else
	{
		printf("unsupported pixel format\n");
		return false;
	}
}
