#include "Locator.h"

extern shared_mutex sm_pData;
extern queue<unsigned char*> pDatas;
extern MV_FRAME_OUT_INFO_EX g_pFrameInfo;
extern bool ispFrameInfoWritten;

extern vector<cv::Mat> USImages;

extern shared_mutex sm_Rt;
extern queue<cv::Mat> R_w2cs, t_w2cs;

shared_mutex sm_camRender;

extern bool ExitFlag;
extern std::mutex mut_Exit;

extern bool LocateExitFlag;
extern std::mutex mut_LocateExit;

bool CamRenderedFlag = false;
cv::Mat camRenderImage; //storing the image rendered in the camera frame window

extern float fps;

Locator::Locator(const string& dataPath, const string& intrinsFile, CamController* camcontroller, 
	double globalResizeFactor, double ArUcoResizeFactor) :
	m_dataPath(dataPath), m_intrinsFile(intrinsFile), m_camcontroller(camcontroller)
	, m_detectArUco(true), m_segHeight(0), m_segWidth(0),
	m_segSize(0), m_patternHeight(0), m_patternWidth(0), m_intervalSize_width(0), m_intervalSize_height(0),
	m_globalResizeFactor(globalResizeFactor), m_ArUcoResizeFactor(ArUcoResizeFactor) {}

Locator::Locator(const string& dataPath, const string& intrinsFile, CamController* camcontroller,
	const int& segWidth, const int& segHeight, const float& minArea, const float& maxArea, double globalResizeFactor, 
	double ArUcoResizeFactor) : m_dataPath(dataPath), m_intrinsFile(intrinsFile), m_camcontroller(camcontroller),
	m_detectArUco(false), m_segHeight(segHeight), m_segWidth(segWidth), m_segSize(0),
	m_patternHeight(0), m_patternWidth(0), m_intervalSize_width(0), m_intervalSize_height(0),
	m_globalResizeFactor(globalResizeFactor), m_ArUcoResizeFactor(ArUcoResizeFactor)
{
	m_SBD_Params.minArea = minArea;
	m_SBD_Params.maxArea = maxArea;
}

bool Locator::init() {
	//Read board pattern and camera intrinsics from intrinsFile
	cv::FileStorage fs(m_intrinsFile, cv::FileStorage::READ);
	fs["image_width"] >> m_imageWidth;
	fs["image_height"] >> m_imageHeight;
	fs["board_width"] >> m_patternWidth;
	fs["board_height"] >> m_patternHeight; 
	fs["square_size_width"] >> m_intervalSize_width;
	fs["square_size_height"] >> m_intervalSize_height;
	fs["camera_matrix"] >> m_cameraMatrix;
	fs["distortion_coefficients"] >> m_distCoeffs;
	cout << "Camera intrinsics£º" << endl << m_cameraMatrix << endl;
	cout << "Deformation coefficient: " << endl << m_distCoeffs << endl;

	if (m_cameraMatrix.empty() || m_distCoeffs.empty())
	{
		std::cerr << "Intrinsics empty! Please check intrinsics data file!" << endl;
		return false;
	}

	//Delete out-dated output file
	char outfile[100];
	strcpy_s(outfile, strlen((m_dataPath + "fileRt.txt").c_str()) + 1, (m_dataPath + "fileRt.txt").c_str());
	if (remove(outfile) == 0)
		cout << "'fileRt.txt' already exists.  Previous file has been deleted." << endl;

	cout << "Locator initiation succeeds.\n" << endl;

	return true;
}

bool Locator::init(int patternHeight, int patternWidth, float intervalSize_width, float intervalSize_height)
{
	//Read camera intrinsics from intrinsFile.  Set board parameters mannually.
	m_patternHeight = patternHeight;
	m_patternWidth = patternWidth;
	m_intervalSize_width = intervalSize_width;
	m_intervalSize_height = intervalSize_height;

	cv::FileStorage fs(m_intrinsFile, cv::FileStorage::READ);
	fs["image_width"] >> m_imageWidth;
	fs["image_height"] >> m_imageHeight;
	fs["camera_matrix"] >> m_cameraMatrix;
	fs["distortion_coefficients"] >> m_distCoeffs;
	cout << "Camera intrinsics£º" << endl << m_cameraMatrix << endl;
	cout << "Deformation coefficient: " << endl << m_distCoeffs << endl;

	if (m_cameraMatrix.empty() || m_distCoeffs.empty())
	{
		std::cerr << "Intrinsics empty! Please check intrinsics data file!" << endl;
		return false;
	}

	//Delete out-dated output file
	char outfile[100];
	strcpy_s(outfile, strlen((m_dataPath + "fileRt.txt").c_str()) + 1, (m_dataPath + "fileRt.txt").c_str());
	if (remove(outfile) == 0)
		cout << "'fileRt.txt' already exists.  Previous file has been deleted." << endl;

	cout << "Locator initiation succeeds.\n" << endl;

	return true;
}

void Locator::run()
{
	int sleep_time = 10;

	while (1)
	{
		/*ÅÐ¶ÏÍË³ö*/
		{
			std::lock_guard<std::mutex> lk_Exit(mut_Exit);
			if (ExitFlag)
			{
				sm_pData.lock_shared();
				bool ispDataEmpty = pDatas.empty();
				sm_pData.unlock_shared();
				if (ispDataEmpty)
				{
					std::lock_guard<std::mutex> lk_LocateExit(mut_LocateExit);
					LocateExitFlag = true;
					return;
				}
			}
		}

		cv::Mat frame0; //Camera original output
		sm_pData.lock();
		bool ispDataEmpty = pDatas.empty();
		if (!ispDataEmpty)
		{
			Convert2Mat(&g_pFrameInfo, pDatas.front(), frame0);
			free(pDatas.front());
			pDatas.pop();
			//Every time the camera gets a frame, a pair of R t will be written (if detection fails, write -9999)
			locate(frame0);
			sm_pData.unlock();
		}
		else
		{
			sm_pData.unlock();
			std::this_thread::sleep_for(std::chrono::milliseconds(sleep_time));
			continue;
		}
	}
}

void Locator::locate(cv::Mat& inputImage)
{
	time_t t_start, t_end;

	if (!CamRenderedFlag)
	{
		inputImage.copyTo(camRenderImage);
		CamRenderedFlag = true;
		thread t_showImage(renderCameraFrame, m_camcontroller);
		t_showImage.detach();
	}

	if (m_detectArUco) //if using ArUco, we use the first frame to set parameters
	{
		if(detectArUco(inputImage))
			m_detectArUco = false;
		vector<double> R(9, -9999);
		vector<double>	t(3, -9999);
		writeRt(R, t, m_dataPath);
		return;
	}

	m_count++;
	bool isFound = true;

	//For frame number 2 or larger
	if (m_count > 1) {
		t_start = clock();

		interFrameDetect(centers, isFound, inputImage);

		if (isFound) {
			drawCenters(inputImage, centers);
		}
		//Inter-frame detection fails.  Take this frame as initial frame.
		else {
			m_count = 1;
			cout << "Detection using the segmented image failed.  Using whole image instead." << endl;
		}

		sm_camRender.lock();
		inputImage.copyTo(camRenderImage);
		sm_camRender.unlock();

		t_end = clock();
		if(isFound)
			cout << "Time Consumption (inter-frame): " << t_end - t_start << endl << endl;
	}//if (count > 1)

	//For frame number 1
	if (m_count == 1) {
		t_start = clock();
		if (!detectFirstFrame(inputImage))
		{
			t_end = clock();
			cout << "Initial frame detection failed." << endl << endl;
		}
		else
		{
			t_end = clock();
			cout << "Time Consumption (initial-frame): " << t_end - t_start << endl << endl;
		}
		sm_camRender.lock();
		inputImage.copyTo(camRenderImage);
		sm_camRender.unlock();
		vector<double> R(9, -9999);
		vector<double>	t(3, -9999);
		writeRt(R, t, m_dataPath);
		return;
	}//if (count == 1)

	vector<cv::Point3f> camPoints; //circle centers in camera coordinates
	calcCamPoints(centers, m_count, camPoints);

	cv::Mat R_w2c;
	m_RotationM.copyTo(R_w2c);
	cv::Mat t_w2c;
	m_tvec.copyTo(t_w2c);

	//Write R and t to txt
	writeRt(R_w2c, t_w2c, m_dataPath);
}

bool Locator::detectArUco(const cv::Mat& inputImage)
{
	static vector<cv::Point2f> corners; //the four corners of the board marked by four aruco markers

	cv::Mat inputImage_resized;
	cv::resize(inputImage, inputImage_resized, cv::Size(), m_ArUcoResizeFactor, m_ArUcoResizeFactor, cv::INTER_AREA);

	//Detect the board corners. Set blobDetector m_SBD_Params accordingly.
	corners.clear();
	corners = detectArucoCorners(inputImage_resized);
	if (corners.size() != 4) //board out of FOV
	{
		cout << "Detect ArUco failed." << endl << endl;
		return false; //load next image
	}
	//Cut original image according to board corners to accelerate
	else // (corners.size() == 4)
	{
		for (auto& c : corners)
		{
			c.x /= m_ArUcoResizeFactor;
			c.y /= m_ArUcoResizeFactor;
		}
		setParams(m_SBD_Params, corners);

		cout << "Detect ArUco succeeded. Set SBD params complete." << endl;
		return true;
	}
}

bool Locator::detectFirstFrame(cv::Mat& inputImage)
{
	bool isFound = true;

	cv::Mat inputImage_resized;
	cv::resize(inputImage, inputImage_resized, cv::Size(), m_globalResizeFactor, m_globalResizeFactor, cv::INTER_AREA);

	m_count = 0;
	cv::Ptr<cv::FeatureDetector> ptrBlobDetector = cv::SimpleBlobDetector::create(m_SBD_Params);

	isFound = cv::findCirclesGrid(inputImage_resized, cv::Size(m_patternWidth, m_patternHeight), centers,
		cv::CALIB_CB_ASYMMETRIC_GRID, ptrBlobDetector);

	if (isFound && centers.size() == m_patternHeight * m_patternWidth) {
		for (size_t i = 0; i < centers.size(); i++) {
			centers[i].x = centers[i].x / m_globalResizeFactor;
			centers[i].y = centers[i].y / m_globalResizeFactor;
		}

		cv::drawChessboardCorners(inputImage, cv::Size(m_patternWidth, m_patternHeight), centers, 1);

		//Only when global first frame detection succeeds, will m_count be increased to 1, otherwise m_count is still 0
		m_count = 1;

		cout << "\nGlobal board detection succeeded." << endl;
		return true;
	}
	//!isFound
	else {
		cout << "\nGlobal board detection failed." << endl;
		return false; //load next image
	}
}

vector<cv::Point2f> Locator::detectArucoCorners(const cv::Mat& src)
{
	std::vector<int> markerIds;
	std::vector<std::vector<cv::Point2f>> markerCorners, rejectedCandidates;
	cv::Ptr<cv::aruco::DetectorParameters> parameters = cv::aruco::DetectorParameters::create();
	//parameters->adaptiveThreshWinSizeMax = 100;
	//parameters->adaptiveThreshWinSizeStep = 5;
	//parameters->maxMarkerPerimeterRate = 0.0;
	parameters->minMarkerPerimeterRate = 0.003;
	//parameters->cornerRefinementMethod = aruco::CORNER_REFINE_SUBPIX;

	//pixels changed from 6 to 4
	cv::Ptr<cv::aruco::Dictionary> dictionary = cv::aruco::getPredefinedDictionary(cv::aruco::DICT_4X4_50);
	//cv::Ptr<cv::aruco::Dictionary> dictionary = cv::aruco::generateCustomDictionary(4, 6);
	cv::aruco::detectMarkers(src, dictionary, markerCorners, markerIds, parameters, rejectedCandidates);
	if (markerCorners.size() != 4) {
		if (markerCorners.size() == 0)
			cout << "Detected NO ArUco corner!" << endl;
		else
			cout << "Detected " << markerCorners.size() << " ArUco corners!. Board is out of vision!" << endl;
		return vector<cv::Point2f>();
	}
	std::vector<cv::Point2f> corners;
	corners.resize(4);
	int idx;
	for (idx = 0; idx < 4; ++idx) {
		if (markerIds[idx] == 0)
			corners[3] = markerCorners[idx][2];
		else if (markerIds[idx] == 1)
			corners[1] = markerCorners[idx][1];
		else if (markerIds[idx] == 2)
			corners[2] = markerCorners[idx][3];
		else if (markerIds[idx] == 3)
			corners[0] = markerCorners[idx][0];
	}
	return corners;
}

void Locator::setParams(cv::SimpleBlobDetector::Params& params, const vector<cv::Point2f>& corners) {
	//Unit: mm
	int boardWidth = 60;
	int boardHeight = 85;
	int whiteOffset = 6.5;
	int circleSize = 4;

	float len1 = sqrt((corners[0].x - corners[1].x) * (corners[0].x - corners[1].x)
		+ (corners[0].y - corners[1].y) * (corners[0].y - corners[1].y)) / (boardWidth - 2 * whiteOffset) * circleSize;
	float len2 = sqrt((corners[0].x - corners[2].x) * (corners[0].x - corners[2].x) 
		+ (corners[0].y - corners[2].y) * (corners[0].y - corners[2].y)) / (boardHeight - 2 * whiteOffset) * circleSize;

	params.maxArea = 1.5 * len1 * len2 * m_globalResizeFactor * m_globalResizeFactor;
	params.minArea = 0.5 * len1 * len2 * m_globalResizeFactor * m_globalResizeFactor;

	m_segSize = 2 * (len1 > len2 ? len1 : len2);
	m_segHeight = m_segSize;
	m_segWidth = m_segSize;
}

void Locator::interFrameDetect(vector<cv::Point2f>& centers, bool& isFound, const cv::Mat& image_orig)
{
	//Split the image into segments as many as the size of circles.  Detect circle in every single segment.
	vector<cv::Mat> image_seg(centers.size());
	vector<cv::Mat> image_border(centers.size());
	vector<bool> seg_found(centers.size());
	vector<int> seg_x(centers.size()); //top-left corner of segment
	vector<int> seg_y(centers.size()); //top-left corner of segment
	vector<int> seg_size(centers.size(), m_segSize); //top-left corner of segment

//Traverse the position of each center of the previous frame of image
	int threadNum = 3;
#pragma omp parallel num_threads(threadNum)
	{
		int threadID = omp_get_thread_num();
		int sectionLength = centers.size() / threadNum + 1;
		int sectionStart = sectionLength * threadID;
		int sectionlEnd = sectionStart + sectionLength;
		if (threadID == threadNum - 1)
			sectionlEnd = centers.size();
		for (int i = sectionStart; i < sectionlEnd; i++)
		{
			seg_x[i] = centers[i].x - seg_size[i] / 2;
			seg_y[i] = centers[i].y - seg_size[i] / 2;

			//Check if the board meets the edge of FOV
			if ((seg_x[i] - seg_size[i] < 0) || (seg_y[i] - seg_size[i] < 0) ||
				(seg_x[i] + seg_size[i] > image_orig.cols) || (seg_y[i] + seg_size[i] > image_orig.rows)) {
				cout << "The board is near the image border, please move it to the center!" << endl;
				isFound = false;
			}

			if (isFound)
			{
				cv::Rect seg(seg_x[i], seg_y[i], m_segSize, m_segSize);
				image_orig(seg).copyTo(image_seg[i]);

				seg_found[i] = detectCircleInSeg(seg_x[i], seg_y[i], image_seg[i], image_border[i], centers[i]);
			}
		}
	}
#pragma omp barrier

	if (isFound == false)
		return;

	for (int j = 0; j < centers.size(); j++)
	{
		if (seg_found[j] == false)
		{
			isFound = false;
			return;
		}
	}
	isFound = true;
}

bool Locator::detectCircleInSeg(const int& seg_x, const int& seg_y, const cv::Mat& image_seg, cv::Mat& image_border, cv::Point2f& center)
{
	cv::cuda::GpuMat image_gpu(image_seg);

	//Threshold
	cv::Mat image_binary;
	int thresh = cv::threshold(image_seg, image_binary, 0, 255, cv::THRESH_BINARY | cv::THRESH_OTSU);
	cv::cuda::threshold(image_gpu, image_gpu, thresh, 255, cv::THRESH_BINARY);

	//Morphology erosion
	cv::Mat mKernel = cv::getStructuringElement(cv::MORPH_RECT, cv::Size(3, 3));
	auto morph_filter = cv::cuda::createMorphologyFilter(cv::MORPH_ERODE, image_seg.type(), mKernel);
	morph_filter->apply(image_gpu, image_gpu);
	cv::Mat image_eroded(image_gpu);

	//Extract circle border and computed center
	image_border = image_binary - image_eroded;
	cv::Mat m_labels, m_stats, m_centroids;
	cv::connectedComponentsWithStats(image_border, m_labels, m_stats, m_centroids);

	bool isFound = true;

	__int32 maxArea = 0;
	int maxAreaIdx = 1;

	for (size_t j = 1; j < m_stats.rows; j++) {
		if (maxArea < m_stats.at<__int32>(j, 4)) {
			maxArea = m_stats.at<__int32>(j, 4);
			maxAreaIdx = j;
		}
	}

	for (size_t m = 0; m < m_segWidth; m++) {
		for (size_t n = 0; n < m_segHeight; n++) {
			if (m_labels.at<__int32>(m, n) == maxAreaIdx) {
				//The border out of FOV
				if ((m < 3) || (m > m_segWidth - 3) || (n < 3) || (n > m_segHeight - 3)) {
					isFound = false;
					break;
				}
			}
			else
				image_border.at<uchar>(m, n) = 0;
		}
	}
	center.x = m_centroids.at<double>(maxAreaIdx, 0) + seg_x;
	center.y = m_centroids.at<double>(maxAreaIdx, 1) + seg_y;

	return isFound;
}

void Locator::calcWorldCords(cv::Size boardSize, float squareSize_width, float squareSize_height, 
	vector<cv::Point3f>& centers)
{
	centers.clear();
	for (int i = 0; i < boardSize.height; i++)
		for (int j = 0; j < boardSize.width; j++)
			centers.push_back(cv::Point3f(float((2 * j + i % 2) * squareSize_width), 
				float(i * squareSize_height), 0));
}

void Locator::calcCamPoints(vector<cv::Point2f>& centers, const int& count, vector<cv::Point3f>& camPoints)
{
	//Calculate world coordinates( Take the plane where the calibration board is located as the z = 0 plane )
	vector<cv::Point3f> WorldPoints;
	calcWorldCords(cv::Size(m_patternWidth, m_patternHeight), 
		m_intervalSize_width, m_intervalSize_height, WorldPoints);

	//Solve Pnp equation to get R, t and reprojection error
	if (count == 2) {
		solvePnPGeneric(WorldPoints, centers, m_cameraMatrix, m_distCoeffs, m_rvecs, m_tvecs,
			false, cv::SOLVEPNP_ITERATIVE, cv::noArray(), cv::noArray(), m_reprojectionError);
	}
	//Using previous frame to accelerate
	else {
		solvePnPGeneric(WorldPoints, centers, m_cameraMatrix, m_distCoeffs, m_rvecs, m_tvecs,
			true, cv::SOLVEPNP_ITERATIVE, m_rvec, m_tvec, m_reprojectionError);
	}
	//solvePnpGeneric returns only one set of solution in SOLVEPNP_ITERATIVE method.  So rvecs, tvecs, reprojectionError are all of size one.
	m_rvec = m_rvecs[0];
	m_tvec = m_tvecs[0];
	cv::Rodrigues(m_rvec, m_RotationM); //rotation vector transforms to rotation matrix

	//Calculate camera coordinates
	camPoints.resize(WorldPoints.size());
	for (int i = 0; i < WorldPoints.size(); i++) {
		camPoints[i].x = m_RotationM.at<double>(0, 0) * WorldPoints[i].x
			+ m_RotationM.at<double>(0, 1) * WorldPoints[i].y
			+ m_RotationM.at<double>(0, 2) * WorldPoints[i].z
			+ m_tvec.at<double>(0, 0);
		camPoints[i].y = m_RotationM.at<double>(1, 0) * WorldPoints[i].x
			+ m_RotationM.at<double>(1, 1) * WorldPoints[i].y
			+ m_RotationM.at<double>(1, 2) * WorldPoints[i].z
			+ m_tvec.at<double>(1, 0);
		camPoints[i].z = m_RotationM.at<double>(2, 0) * WorldPoints[i].x
			+ m_RotationM.at<double>(2, 1) * WorldPoints[i].y
			+ m_RotationM.at<double>(2, 2) * WorldPoints[i].z
			+ m_tvec.at<double>(2, 0);
	}
}

void Locator::drawCenters(cv::Mat image, const vector<cv::Point2f>& centers) {
	for (size_t i = 0; i < centers.size(); i++)
		cv::circle(image, centers[i], 8, cv::Scalar(255, 255, 255), -1);
}

void renderCameraFrame(CamController* camcontroller)
{
	cv::namedWindow("image_show", cv::WINDOW_NORMAL);
	cv::resizeWindow("image_show", camRenderImage.cols / 8, camRenderImage.rows / 8);

	do {
		{
			std::lock_guard<std::mutex> lk_Exit(mut_Exit);
			if (cv::waitKey(10) == 27 || ExitFlag == true)
			{
				camcontroller->CloseCam();
				ExitFlag = true;
				return;
			}
		}
		
		sm_camRender.lock_shared();
		cv::Mat temp = camRenderImage.clone(); //Need to be cloned because camRenderImage may be changed when rendering
		sm_camRender.unlock_shared();

		cv::imshow("image_show", temp);
	} while (cv::waitKey(1));
}

vector<double> Mat2vector(const cv::Mat& mat)
{
	vector<double> vec;
	for (int i = 0; i < mat.rows; ++i) {
		for (int j = 0; j < mat.cols; ++j) {
			vec.push_back(mat.at<double>(i, j));
		}
	}
	return vec;
}

void saveImages(const string& path, const string& name)
{
	cout << "Saving captured images......" << endl;
	int count = 1;
	for (auto& image : USImages)
		cv::imwrite(path + name + std::to_string(count++) + ".bmp", image);
	cout << "Save captured images done." << endl;
}

cv::Mat vector2Mat(const vector<double>& vec)
{
	cv::Mat mat;
	if (vec.size() == 3)
	{
		mat.create(3, 1, CV_64FC1);
	}
	else if (vec.size() == 9)
		mat.create(3, 3, CV_64FC1);
	else
	{
		cout << "In vector2Mat: wrong input size!" << endl;
		return cv::Mat();
	}
	for (int i = 0; i < vec.size(); ++i)
	{
		mat.at<double>(i) = vec[i];
	}
	return mat;
}

void writeRt(vector<double> R, vector<double> t, const string& path)
{
	cv::Mat Rmat = vector2Mat(R);
	cv::Mat tmat = vector2Mat(t);
	sm_Rt.lock();
	R_w2cs.push(Rmat);
	t_w2cs.push(tmat);
	sm_Rt.unlock();

	std::ofstream fileRt;
	fileRt.open(path + "fileRt.txt", std::ios::out | std::ios::app);

	if (fileRt.is_open())
	{
		for (int i = 0; i < R.size(); i++)
			fileRt << R[i] << " ";
		fileRt << "\n";
		for (int i = 0; i < t.size(); i++)
			fileRt << t[i] << " ";
		fileRt << "\n";
		fileRt.close();
	}
	else
		cout << "Can't write fileRt.txt" << endl;

	return;
}
void writeRt(cv::Mat R, cv::Mat t, const string& path)
{
	sm_Rt.lock();
	R_w2cs.push(R);
	t_w2cs.push(t);
	sm_Rt.unlock();

	std::ofstream fileRt;
	fileRt.open(path + "fileRt.txt", std::ios::out | std::ios::app);

	if (fileRt.is_open())
	{
		for (int i = 0; i < R.rows; i++) {
			for (int j = 0; j < R.cols; ++j) {
				fileRt << R.at<double>(i, j) << " ";
			}
		}
		fileRt << "\n";
		for (int i = 0; i < t.rows; i++)
			fileRt << t.at<double>(i) << " ";
		fileRt << "\n";
		fileRt.close();
	}
	else
		cout << "Can't write fileRt.txt" << endl;

	return;
}

