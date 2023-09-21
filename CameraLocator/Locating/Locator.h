#pragma once
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <vector>
#include <queue>
#include <opencv2/opencv.hpp>
#include <thread>
#include <shared_mutex>
#include <chrono>
#include <opencv2/aruco.hpp>
#include <omp.h>
#include "CamControl.h"

using std::queue;
using std::shared_mutex;
using std::thread;
using std::vector;
using std::string;
using std::queue;
using std::cout;
using std::endl;

constexpr auto PI = 3.14159265358979323846264338327950288;

class Locator {

public:
	//Constructor using aruco
	//Parameters:
	//		dataPath -- Path for storing images and output texts
	//		intrinsFile -- xml or yaml file generated from opencv calibration sample
	//		resizeFactor -- how much the camera frame will be resized for initial frame detection
	//		ArUcoFactor -- how much the camera frame will be resized for ArUco detection
	Locator(const string& dataPath, const string& intrinsFile, CamController* camcontroller, 
		double globalResizeFactor = 0.3, double ArUcoResizeFactor = 0.3);

	//Constructor not using aruco. Set segment size and SimpleBlobDetector params explicitly
	Locator(const string& dataPath, const string& intrinsFile, CamController* camcontroller, const int& segWidth,
		const int& segHeight, const float& minArea, const float& maxArea, double resizeFactor = 0.3, 
		double ArUcoResizeFactor = 0.5);

	//Initiation Function
	bool init();
	//Reload of initiation method.  Set board parameters mannually instead of reading from intrinsFile
	bool init(int patternHeight, int patternWidth, float intervalSize_width, float intervalSize_height);

	void run();

private:
	bool detectArUco(const cv::Mat& inputImage);

	//Detect the first frame to provide location of circles
	bool detectFirstFrame(cv::Mat& inputImage);

	//Locate circles and write their coordinates into file
	void locate(cv::Mat& inputImage);

	//Return the four outer-corners' coordinates.  Four ArUco markers are defined in the following order
	//  M3														M1
	//      *		*		*		*		*		*		*
	//			*		*		*		*		*		*
	//      *		*		*		*		*		*		*
	//			*		*		*		*		*		*
	//  M4														M2
	vector<cv::Point2f> detectArucoCorners(const cv::Mat& src);

	//Set parameters of SimpleblobDetector and segment size, based on location of aruco markers.
	void setParams(cv::SimpleBlobDetector::Params& params, const vector<cv::Point2f>& corners);

	//Inter-frame detection of circle centers
	void interFrameDetect(vector<cv::Point2f>& centers, bool& isFound, const cv::Mat& image_orig);

	//Circle center detection based on previous frame of image
	bool detectCircleInSeg(const int& seg_x, const int& seg_y, const cv::Mat& image_seg, cv::Mat& image_border, cv::Point2f& center);

	//Calculate world coordinates( board plane is z = 0 plane )
	static void calcWorldCords(cv::Size boardSize, float squareSize_width, float squareSize_height, 
		vector<cv::Point3f>& centers);

	//Calculate circle centers in camera coordinates
	void calcCamPoints(vector<cv::Point2f>& centers, const int& count, vector<cv::Point3f>& camPoints);

	//Draw circle centers.
	void drawCenters(cv::Mat image, const vector<cv::Point2f>& centers);

private:
	
	int m_count = 0; //frame counter

	bool m_detectArUco;

	string m_dataPath;
	string m_intrinsFile;

	double m_globalResizeFactor; //how much the camera frame will be resized for initial frame detection
	double m_ArUcoResizeFactor;

	int m_imageHeight;
	int m_imageWidth;

	int m_patternHeight;
	int m_patternWidth;
	float m_intervalSize_width; //half of the distance between two adjacent circles on the same row
	float m_intervalSize_height;

	int m_segWidth, m_segHeight, m_segSize;
	cv::SimpleBlobDetector::Params m_SBD_Params;

	cv::Mat m_cameraMatrix, m_distCoeffs;
	cv::Mat m_rvec, m_tvec;
	vector<cv::Mat> m_rvecs, m_tvecs;
	vector<double> m_reprojectionError;
	cv::Mat m_RotationM;
	vector<cv::Point2f> centers; //image coordinates of circle centers
public:
	CamController* m_camcontroller; //when exiting, shut off the camera
};
/****** End of Locator ******/

//Used in an independent thread. For rendering the camera frame.
void renderCameraFrame(CamController* camcontroller);

void writeRt(vector<double> R, vector<double> t, const string& path);

void writeRt(cv::Mat R, cv::Mat t, const string& path);

vector<double> Mat2vector(const cv::Mat& mat);

//Save screen captured images
void saveImages(const string& path, const string& name);

cv::Mat vector2Mat(const vector<double>& vec);

