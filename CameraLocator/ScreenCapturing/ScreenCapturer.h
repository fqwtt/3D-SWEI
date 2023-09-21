#pragma once
#include <Windows.h>
#include <windowsx.h>

#include <iostream>
#include <algorithm>
#include <string>
#include <vector>
#include <conio.h>
#include <atlbase.h> // 字符串转换用到
#include <functional>
#include <memory>
#include <numeric>

using std::vector;

namespace Gdiplus
{
	using std::min;
	using std::max;
};

#include <Gdiplus.h> // 保存图片用到了GDI+
#include <dwmapi.h>
#include <opencv2/opencv.hpp>

#pragma comment(lib, "gdiplus.lib") // 保存图片需要
#pragma comment(lib, "Dwmapi.lib")  // 判断是否是隐形窗口以及获取窗口大小会用到

// 为了将屏幕和窗口进行统一,因此使用了结构体
struct WindowInfo
{
	HWND hwnd; /* 为空表示屏幕截图 */
	std::string desc; // 窗口标题
	RECT rect{ 0,0,0,0 }; /* hwnd不为空时,此参数无效 */
};

static int GetEncoderClsid(const WCHAR* format, CLSID* pClsid)
{
	UINT  num = 0;          // number of image encoders
	UINT  size = 0;         // size of the image encoder array in bytes

	Gdiplus::ImageCodecInfo* pImageCodecInfo = NULL;

	Gdiplus::GetImageEncodersSize(&num, &size);
	if (size == 0)
		return -1;  // Failure

	pImageCodecInfo = (Gdiplus::ImageCodecInfo*)(malloc(size));
	if (pImageCodecInfo == NULL)
		return -1;  // Failure

	Gdiplus::GetImageEncoders(num, size, pImageCodecInfo);

	for (UINT j = 0; j < num; ++j) {
		if (wcscmp(pImageCodecInfo[j].MimeType, format) == 0) {
			*pClsid = pImageCodecInfo[j].Clsid;
			free(pImageCodecInfo);
			return j;  // Success
		}
	}

	free(pImageCodecInfo);
	return -1;  // Failure
}

// 将bitmap对象保存为bmp图片
static bool SaveBitmapAsBmp(const std::shared_ptr<Gdiplus::Bitmap>& bitmap, const std::string& filename)
{
	if (bitmap == nullptr) return false;
	CLSID bmp_clsid;
	GetEncoderClsid(L"image/bmp", &bmp_clsid);
	Gdiplus::Status ok = bitmap->Save(CA2T(filename.c_str(), CP_ACP), &bmp_clsid, nullptr);
	return ok == Gdiplus::Status::Ok;
}

class Enumerator
{
public:
	using EnumCallback = std::function<void(const WindowInfo&)>;

	static bool EnumMonitor(EnumCallback callback)
	{
		// 调用Win32Api进行显示器遍历
		return ::EnumDisplayMonitors(NULL, NULL, MonitorEnumProc, (LPARAM)&callback);
	}

	static bool EnumWindow(EnumCallback callback)
	{
		// 调用Win32Api进行窗口遍历
		return ::EnumWindows(EnumWindowsProc, (LPARAM)&callback);
	}

private:
	static BOOL CALLBACK EnumWindowsProc(HWND hwnd, LPARAM lParam)
	{
		//::GetParent获取的有可能是所有者窗口,因此使用GetAncestor获取父窗口句柄
		HWND parent = ::GetAncestor(hwnd, GA_PARENT);
		HWND desktop = ::GetDesktopWindow(); // 获取桌面的句柄
		TCHAR szTitle[MAX_PATH] = { 0 };
		::GetWindowText(hwnd, szTitle, MAX_PATH); // 获取标题

		// 排除image_show窗口
		if (_tcscmp(szTitle, _T("image_show")) == 0) return TRUE;

		// 排除父窗口不是桌面的
		if (parent != nullptr && parent != desktop) return TRUE;

		// 排除标题为空的
		if (wcscmp(szTitle, L"") == 0) return TRUE;

		// 排除最小化窗口(因为获取最小化窗口的区域数据是不对的,因此也没办法进行截图等操作)
		if (::IsIconic(hwnd)) return TRUE;

		// 排除不可见窗口,被其他窗口遮挡的情况是可见的
		if (!::IsWindowVisible(hwnd)) return TRUE;

		// 排除对用户隐形的窗口,参考[https://docs.microsoft.com/en-us/windows/win32/api/dwmapi/ne-dwmapi-dwmwindowattribute]
		DWORD flag = 0;
		DwmGetWindowAttribute(hwnd, DWMWA_CLOAKED, &flag, sizeof(flag));
		if (flag) return TRUE;

		if (lParam) {
			WindowInfo wnd_info{ hwnd,(LPCSTR)CT2A(szTitle, CP_ACP) };
			EnumCallback* callback_ptr = reinterpret_cast<EnumCallback*>(lParam);
			callback_ptr->operator()(wnd_info);
		}
		return TRUE;
	}

	static BOOL CALLBACK MonitorEnumProc(HMONITOR hMonitor, HDC hdcMonitor, LPRECT lprcMonitor, LPARAM dwData)
	{
		MONITORINFOEX mi;
		mi.cbSize = sizeof(MONITORINFOEX);
		GetMonitorInfo(hMonitor, &mi);
		if (dwData) {
			std::string device_name = (LPCSTR)CT2A(mi.szDevice, CP_ACP);
			if (mi.dwFlags == MONITORINFOF_PRIMARY) device_name += "(Primary)"; // 主显示器,可根据需要进行操作
			WindowInfo wnd_info{ nullptr, device_name, mi.rcMonitor };

			EnumCallback* callback = reinterpret_cast<EnumCallback*>(dwData);
			(*callback)(wnd_info);
		}
		return TRUE;
	}
};

class WindowCapture
{
public:
	using BitmapPtr = std::shared_ptr<Gdiplus::Bitmap>;

	static BitmapPtr Capture(const WindowInfo& wnd_info)
	{
		HDC hWndDC = GetWindowDC(wnd_info.hwnd);
		RECT capture_rect{ 0,0,0,0 }; // 最终要截取的区域
		RECT wnd_rect; // 窗口区域
		RECT real_rect; // 真实的窗口区域,实际上也不是百分百准确

		if (wnd_info.hwnd) {
			::GetWindowRect(wnd_info.hwnd, &wnd_rect);
			DwmGetWindowAttribute(wnd_info.hwnd, DWMWINDOWATTRIBUTE::DWMWA_EXTENDED_FRAME_BOUNDS, &real_rect, sizeof(RECT));
			int offset_left = real_rect.left - wnd_rect.left;
			int offset_top = real_rect.top - wnd_rect.top;
			capture_rect = RECT{ offset_left,offset_top,real_rect.right - real_rect.left + offset_left,real_rect.bottom - real_rect.top + offset_top };
		}
		else {
			capture_rect = wnd_info.rect;
		}

		int width = capture_rect.right - capture_rect.left;
		int height = capture_rect.bottom - capture_rect.top;

		HDC hMemDC = CreateCompatibleDC(hWndDC);
		HBITMAP hBitmap = CreateCompatibleBitmap(hWndDC, width, height);
		SelectObject(hMemDC, hBitmap);

		BitmapPtr bitmap;
		// 获取指定区域的rgb数据
		bool ok = BitBlt(hMemDC, 0, 0, width, height, hWndDC, capture_rect.left, capture_rect.top, SRCCOPY);
		// hBitmap就是得到的图片对象,转GDI的Bitmap进行保存
		if (ok) bitmap = std::make_shared<Gdiplus::Bitmap>(hBitmap, nullptr);

		DeleteDC(hWndDC);
		DeleteDC(hMemDC);
		DeleteObject(hBitmap);

		return bitmap;
	}
};
#if 1
namespace ScreenCapturer
{
	static cv::Mat cropUSImage(const cv::Mat& src, const int& width, const int& height)
	{
		//L7-3
		const int tl_x = 43;
		const int tl_y = 92;
		cv::Rect rect = cv::Rect(tl_x, tl_y, width, height);
		return src(rect).clone();
	}

	static bool Init(const std::string& WindowName, WindowInfo& windowInfo, int& width, int& height)
	{
		std::vector<WindowInfo> window_vec; // 用来保存窗口信息
		// 枚举显示器
		Enumerator::EnumMonitor([&window_vec](const WindowInfo& wnd_info)
			{
				window_vec.push_back(wnd_info);
			});
		// 计算生成所有屏幕加在一起的区域大小
		if (window_vec.size() > 1) { // 也可大于1,这样只有一个显示器时不会显示全屏选项
			int width = 0, height = 0;
			for (const auto& wnd_info : window_vec) {
				width += wnd_info.rect.right - wnd_info.rect.left;
				int h = wnd_info.rect.bottom - wnd_info.rect.top;
				if (h > height) height = h; // 高度可能不一样,需要以最高的为准
			}
			WindowInfo wnd_info{ nullptr, "FullScreen", { 0, 0, width, height} };
			window_vec.push_back(wnd_info);
		}
		// 枚举窗口
		Enumerator::EnumWindow([&window_vec](const WindowInfo& wnd_info)
			{
				window_vec.push_back(wnd_info);
			});

		bool isWindowExist = false;
		for (const auto& window : window_vec) {
			//cout << window.desc << endl;
			if (window.desc == WindowName)
			{
				windowInfo = window;
				isWindowExist = true;
			}
		}
		if (!isWindowExist)
			return false;

		auto bitmap = WindowCapture::Capture(windowInfo);
		if (!bitmap)
		{
			std::cout << "Failed to get bitmap!" << std::endl;
			return false;
		}

		SaveBitmapAsBmp(bitmap, "C:\\Users\\USImage_temp.bmp");
		cv::Mat temp = cv::imread("C:\\Users\\USImage_temp.bmp", cv::IMREAD_GRAYSCALE);
		//img = convertBitmap2Mat(bitmap);
		if (temp.empty())
			return false;
#if 0
		//need to be changed
		const int tl_x = 43;
		const int tl_y = 92;

		cv::Rect rect1 = cv::Rect(tl_x, tl_y, temp.cols - tl_x, temp.rows - tl_y);
		cv::Mat mid = temp(rect1).clone();

		int r = mid.rows - 1; // the last row
		// Scan each row in reverse order
		// Find the 1st row of all white pixels
		for (; r > 0; --r)
		{
			bool isAllWhite = true;
			for (int c = 0; c < mid.cols; ++c)
			{
				int v = *(mid.ptr(r, 1));
				if (v < 190)
				{
					isAllWhite = false;
					break;
				}
			}

			if (isAllWhite)
				break;
		}
		if (r < 1)
		{
			std::cout << "Failed to find all white row in image!" << std::endl;
			return false;
		}

		// Find the next row of >=100 continuous black pixels
		for (--r; r > 0; --r)
		{
			int row_total = 0;
			for (int c = 0; c < 100; ++c)
			{
				int v = *(mid.ptr(r, c));
				row_total += v;
			}

			if (row_total < 2000)
				break;
		}
		if (r < 1)
		{
			std::cout << "Failed to find all black row in image!" << std::endl;
			return false;
		}

		// Find the corner in this row
		int i = 0;
		for (; i < mid.cols; ++i)
		{
			int v = *(mid.ptr(r, i));
			if (v > 190)
				break;
		}
		if (i == mid.cols)
		{
			std::cout << "Failed to find the corner in image!" << std::endl;
			return false;
		}

		int br_x = i - 1;
		int br_y = r;

		width = br_x;
		height = br_y;
		cv::Mat out = temp(cv::Rect(tl_x, tl_y, width, height));
#endif
		return true;
	}

	static cv::Mat Capture(const WindowInfo& window, const int& cap_width, const int& cap_height)
	{
		auto bitmap = WindowCapture::Capture(window);
		if (bitmap)
		{
			SaveBitmapAsBmp(bitmap, "C:\\Users\\USImage_temp.bmp");
			cv::Mat temp = cv::imread("C:\\Users\\USImage_temp.bmp", cv::IMREAD_COLOR);
			//img = convertBitmap2Mat(bitmap);
			if (temp.empty())
				return cv::Mat();
			else
				//image = temp;
				//return cropUSImage(temp, cap_width, cap_height);
				return temp;
		}
	}
}
#else
namespace ScreenCapturer
{
	static cv::Mat cropUSImage(const cv::Mat& src, const int& tl_x, const int& tl_y, 
		const int& width, const int& height)
	{
		cv::Rect rect = cv::Rect(tl_x, tl_y, width, height);
		return src(rect).clone();
	}

	static void extractRect(cv::Mat& imageSource, int& tl_x, int& tl_y, int& width, int& height)
	{
		//cv::imwrite("window.bmp", source);
		using namespace cv;
		int min_x = 0, max_x = 0, min_y = 0, max_y = 0;
		int a = 0;

		Mat image, edge;
		GaussianBlur(imageSource, image, Size(3, 3), 0);
		Canny(image, edge, 100, 250, 3);
		vector<vector<Point>> contours;
		vector<Vec4i> hierarchy;
		findContours(edge, contours, hierarchy, RETR_TREE, CHAIN_APPROX_NONE, Point());
		for (int i = 0; i < contours.size(); i++)
		{
			//contours[i]代表的是第i个轮廓，contours[i].size()代表的是第i个轮廓上所有的像素点数
			if ((1000 < contours[i].size()) && (hierarchy[i][2] != -1))
			{
				a = a + 1;
				vector<int> x(contours[i].size()), y(contours[i].size());
				for (int j = 0; j < contours[i].size(); j++)
				{
					x[j] = contours[i][j].x;
					y[j] = contours[i][j].y;
				}
				max_x = *max_element(x.begin(), x.end());
				min_x = *min_element(x.begin(), x.end());
				max_y = *max_element(y.begin(), y.end());
				min_y = *min_element(y.begin(), y.end());
			}
			if (max_x != 0 && a == 3)break;
		}

		tl_x = min_x + 1;
		tl_y = min_y + 1;
		width = max_x - min_x - 1;
		height = max_y - min_y;

		cv::Mat out = imageSource(Range(min_y, max_y), Range(min_x, max_x));
		return;
	}

	static bool Init(const std::string& WindowName, WindowInfo& windowInfo, int& tl_x, int& tl_y, 
		int& width, int& height)
	{
		std::vector<WindowInfo> window_vec; // 用来保存窗口信息
		// 枚举显示器
		Enumerator::EnumMonitor([&window_vec](const WindowInfo& wnd_info)
			{
				window_vec.push_back(wnd_info);
			});
		// 计算生成所有屏幕加在一起的区域大小
		if (window_vec.size() > 1) { // 也可大于1,这样只有一个显示器时不会显示全屏选项
			int width = 0, height = 0;
			for (const auto& wnd_info : window_vec) {
				width += wnd_info.rect.right - wnd_info.rect.left;
				int h = wnd_info.rect.bottom - wnd_info.rect.top;
				if (h > height) height = h; // 高度可能不一样,需要以最高的为准
			}
			WindowInfo wnd_info{ nullptr, "FullScreen", { 0, 0, width, height} };
			window_vec.push_back(wnd_info);
		}
		// 枚举窗口
		Enumerator::EnumWindow([&window_vec](const WindowInfo& wnd_info)
			{
				window_vec.push_back(wnd_info);
			});

		bool isWindowExist = false;
		for (const auto& window : window_vec) {
			//cout << window.desc << endl;
			if (window.desc == WindowName)
			{
				windowInfo = window;
				isWindowExist = true;
			}
		}
		if (!isWindowExist)
			return false;

		auto bitmap = WindowCapture::Capture(windowInfo);
		if (!bitmap)
		{
			std::cout << "Failed to get bitmap!" << std::endl;
			return false;
		}

		SaveBitmapAsBmp(bitmap, "C:\\Users\\USImage_temp.bmp");
		cv::Mat source = cv::imread("C:\\Users\\USImage_temp.bmp", cv::IMREAD_GRAYSCALE);
		//img = convertBitmap2Mat(bitmap);
		if (source.empty())
			return false;

		extractRect(source, tl_x, tl_y, width, height);

		return true;
	}

	static cv::Mat Capture(const WindowInfo& window, const int& tl_x, const int& tl_y, 
		const int& cap_width, const int& cap_height)
	{
		auto bitmap = WindowCapture::Capture(window);
		if (bitmap)
		{
			SaveBitmapAsBmp(bitmap, "C:\\Users\\USImage_temp.bmp");
			cv::Mat temp = cv::imread("C:\\Users\\USImage_temp.bmp", cv::IMREAD_GRAYSCALE);
			//img = convertBitmap2Mat(bitmap);
			if (temp.empty())
				return cv::Mat();
			else
				//image = temp;
				return cropUSImage(temp, tl_x, tl_y, cap_width, cap_height);
		}
	}
}
#endif