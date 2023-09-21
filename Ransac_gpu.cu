#include "mex.h"
#include "gpu/mxGPUArray.h"

using namespace std;
void __global__ calculateImRe(double* idata,double* qdata,double* Im1,double* Re1,double* Im2,double* Re2,
        int depth,int width,int na)
{
    int id = (blockIdx.x * blockDim.x + threadIdx.x);
	if (id < depth * width * (na - 1))
	{   
		Im1[id] = qdata[id] * idata[id + depth * width] - \
			idata[id] * qdata[id + depth * width];
		Re1[id] = idata[id] * idata[id + depth * width] + \
			qdata[id] * qdata[id + depth * width];
	}
	if (id < (depth - 1) * width * na)
	{
        int i = id % (depth - 1);
		int j = id / (depth - 1);
        id= j * depth + i;
		Im2[id] = qdata[id] * idata[id + 1] - \
			idata[id] * qdata[id + 1];
		Re2[id] = idata[id] * idata[id + 1] + \
			qdata[id] * qdata[id + 1];
	}
}

__global__ void caculate_dis_withM(double* Im1, double* Re1, double* Im2, double* Re2,
	double* Im1_M, double* Re1_M, double* Im2_M, double* Re2_M, int depth, int width, int na, int M,  int N)
{
	int id = (blockIdx.x * blockDim.x + threadIdx.x);
	if (0 <= id && id < (depth - M + 1) * width * (na - 1))
	{
		int i = id % (depth - M + 1);
		int j = id / (depth - M + 1);
        double Im1_sum = 0.0;
		double Re1_sum = 0.0;
		for (int m = 0; m < M; ++m)
		{
            Im1_sum += Im1[j * depth + i + m];
            Re1_sum += Re1[j * depth + i + m];
		}
        Im1_M[j * depth + i] = Im1_sum;
        Re1_M[j * depth + i] = Re1_sum;
    }
    
    id = id - (depth - M + 1) * width * (na - 1);
    if (0 <= id && id < (depth - M + 1) * width * na)
	{
        int i = id % (depth - M + 1);
		int j = id / (depth - M + 1);
        double Im2_sum = 0.0;
		double Re2_sum = 0.0;
		for (int m = 0; m < M - 1; ++m)
		{
            Im2_sum += Im2[j * depth + i + m];
            Re2_sum += Re2[j * depth + i + m];
		}
        Im2_M[j * depth + i] = Im2_sum;
        Re2_M[j * depth + i] = Re2_sum;
	}
}
__global__ void caculate_dis_withN(double* Im1_N, double* Re1_N, double* Im2_N, double* Re2_N,
	double* Im1_M, double* Re1_M, double* Im2_M, double* Re2_M, int depth, int width, int na, int M,  int N)
{
	int id = (blockIdx.x * blockDim.x + threadIdx.x);
    
	if (0 <= id && id < (depth - M + 1) * width * (na - N + 1 ))
	{
		int i = id % (depth - M + 1);
		int j = id / (depth - M + 1) % width;
        int k = id / (depth - M + 1) / width;
		id = k * depth * width + j * depth + i;
        double Im1_sum = 0.0;
		double Re1_sum = 0.0;
		for (int n = 0; n < N - 1; ++n)
		{
            Im1_sum += Im1_M[id + n * depth * width];
            Re1_sum += Re1_M[id + n * depth * width];
		}   
        Im1_N[id] = Im1_sum;
        Re1_N[id] = Re1_sum;
    }
    id = id - (depth - M + 1) * width * (na - N +1);
    if (0 <= id && id < (depth - M + 1) * width * (na - N + 1))
	{
        int i = id % (depth - M + 1);
		int j = id / (depth - M + 1) % width;
        int k = id / (depth - M + 1) / width;
        id = k * depth * width + j * depth + i;
        double Im2_sum = 0.0;
		double Re2_sum = 0.0;
		for (int n = 0; n < N; ++n)
		{
            Im2_sum += Im2_M[id + n * depth * width];
            Re2_sum += Re2_M[id + n * depth * width];
		}
        Im2_N[id] = Im2_sum;
        Re2_N[id] = Re2_sum;
    }
}
__global__ void caculate_dis_without_conv(double* Dis, double* Im1, double* Re1, double* Im2, double* Re2,
	int depth, int width, int na, int M,  int N, double c, double fc)
{
	int id = (blockIdx.x * blockDim.x + threadIdx.x);
	if (id < (depth - M + 1) * width * (na - N + 1 ))
	{
		int i = id % (depth - M + 1);
		int j = id / (depth - M + 1) % width;
        int k = id / (depth - M + 1) / width;
        int id1 = k * depth * width + j * depth + i;
		Dis[id] = atan2(Im1[id1], Re1[id1]) / (1 + atan2(Im2[id1], Re2[id1]) / (2.0 * 3.1416)) * c / (4.0 * 3.1416 * fc);
    }
}
__global__ void caculate_dis_with_conv(double* Dis_k2D, double* Dis,
	int depth, int width, int na)
{
	double kernel2D[5][5] = {
		{ 0.0073, 0.0208, 0.0294, 0.0208, 0.0073 },
		{ 0.0208, 0.0589, 0.0833, 0.0589, 0.0208 },
		{ 0.0294, 0.0833, 0.1179, 0.0833, 0.0294 },
		{ 0.0208, 0.0589, 0.0833, 0.0589, 0.0208 },
		{ 0.0073, 0.0208, 0.0294, 0.0208, 0.0073 }
	};
	int id = (blockIdx.x * blockDim.x + threadIdx.x);
	if (id < depth * width * na)
	{
		int i = id % depth;
		int j = id / depth % width;
		double sum = 0.0;
		for (int m = -2; m < 3; ++m)
		{
			for (int n = -2; n < 3; ++n)
			{
				if (i + m >= 0 && i + m < depth && j + n >= 0 && j + n < width)
				{
					sum += kernel2D[m + 2][n + 2] * Dis[id + n * depth + m];
				}
			}
		}
		Dis_k2D[id] = sum;
	}
}
__global__ void caculate_ave_dis_in_axial(double* Dis_aver, double* Dis_k2D,
    int averN, int depth, int width, int na)
{
	int id = (blockIdx.x * blockDim.x + threadIdx.x);
	if (id < depth * width * na)
	{
		int i = id % depth;
		double sum = 0.0;
		for (int m = (-averN + 1) / 2; m <= averN / 2; ++m)
		{
			if (i + m >= 0 && i + m < depth)
			{
				sum += Dis_k2D[id + m];
			}
		}
		Dis_aver[id] = sum / averN;
	}
}

__global__ void caculate_min_index(double* minIndex, double* Dis,
	int depth, int width, int na)
{
	int id = (blockIdx.x * blockDim.x + threadIdx.x);
	if (id < depth * width)
	{
		double minValue = 0;
		double minInd = 0;
		for (int k = 0; k < na; ++k)
		{
			if (Dis[id + k * depth * width] < minValue)
			{
				minValue = Dis[id + k * depth * width];
				minInd = (double)k;
			}
		}
		minIndex[id] = minInd + 1;
	}
}

// __device__ double polyfit_wt(double* points, int pointNum)
// {   
//     double t1=0.0;
//     double t2=0.0;
//     double t3=0.0;
//     double t4=0.0;
//     for (int i = 0; i < pointNum; ++i)
//     {
//         t1 += points[2 * i] * points[2 * i];
//         t2 += points[2 * i];
//         t3 += points[2 * i + 1] * points[2 * i];
//         t4 += points[2 * i + 1];
//     }
//     return (t3 * pointNum - t2 * t4) / (pointNum * t1 - t2 * t2);
// }

// __device__ double fitFunck(double x1, double y1, double x2, double y2)
// {
//     return (y2 - y1) / (x2 - x1);
// }
// 
// __device__ double fitFuncb(double x1, double y1, double x2, double y2)
// {
//     return (x2 * y1 - x1 * y2) / (x2 - x1);
// }

// __device__ int computeLoopNumber(double confidence, int pointNum, int inlierNum)
// {
//     double inlierProbability = (inlierNum / pointNum) * (inlierNum / pointNum);
//     double conf = 0.01 * confidence;
//     double num = log10(1 - conf);
//     double den = log10(1 - inlierProbability);
//     return (int)ceil(num / den);
// }

__device__ double msac_wt(double* points, int pointNum, double* randpoint, double maxDistance)
{   
    int numTrials  = 50;
    double bestDis = maxDistance * pointNum;
    double bestk = 0.0;
    double bestb = 0.0;
    double accDis = 0.0;
    double x1 = 0.0;
    double x2 = 0.0;
    double y1 = 0.0;
    double y2 = 0.0;
    double k  = 0.0;
    double b  = 0.0;
    double loss = 0.0;
    for (int j = 0; j < numTrials; ++j)
    {
        x1 = randpoint[j * 2];
        x2 = randpoint[j * 2 + 1];
        y1 = points[(int)x1];
        y2 = points[(int)x2];
        k  = (y2 - y1) / (x2 - x1);
        b  = (x2 * y1 - x1 * y2) / (x2 - x1);
        accDis = 0.0;
        for(int i = 0; i < pointNum; ++i)
        {
            loss=fabs((points[i]-k*i-b)/sqrt(k*k+b*b));
            accDis+=min(loss,maxDistance);
        }
        if (accDis < bestDis)
        {
            bestDis = accDis;
            bestk = k;
            bestb = b;
        }
    }
    double t1=0.0;
    double t2=0.0;
    double t3=0.0;
    double t4=0.0;
    int inlierNum = 0;
    for (int i = 0; i < pointNum; ++i)
    {
        loss=fabs((points[i]-bestk*i-bestb)/sqrt(bestk*bestk+bestb*bestb));
        if (loss < maxDistance)
        {
            inlierNum += 1;
            t1 += i*i;
            t2 += i;
            t3 += points[i] * i;
            t4 += points[i];
        }
    }
    return (t3 * inlierNum - t2 * t4) / (inlierNum * t1 - t2 * t2);
}


__global__ void caculate_Yang(double* Yang, double* minInd, double* randpoint, int pointNum,
	int Ledge, int Redge, int depth, double unitDistance, double unitTime)
{
    int id = (blockIdx.x * blockDim.x + threadIdx.x);
    if (id >= depth * Ledge && id < depth * Redge)
    {
        double* points = new double[pointNum];
        for (int i = 0; i < pointNum; ++i)
        {
            points[i]=minInd[id+(i-(pointNum-1)/2)*depth];
        }
        Yang[id]=pow(unitDistance/(unitTime * msac_wt(points,pointNum,randpoint,0.5)),2) * 3000;
        delete []points;
    }
}

void mexFunction(int nlhs, mxArray *plhs[],int nrhs,const mxArray *prhs[])
{   
    
    mxInitGPU();
    int const threadsPerBlock = 1024;
    mxGPUArray const *idata      = mxGPUCreateFromMxArray(prhs[0]);
    mxGPUArray const *qdata      = mxGPUCreateFromMxArray(prhs[1]);
    mxGPUArray const *randpoint  = mxGPUCreateFromMxArray(prhs[2]);
    double *parameter            = mxGetPr(prhs[3]);
    
    double* d_idata    = (double*)(mxGPUGetDataReadOnly(idata));
    double* d_qdata    = (double*)(mxGPUGetDataReadOnly(qdata));
    double* d_randpoint= (double*)(mxGPUGetDataReadOnly(randpoint));
    mxClassID      cid = mxGPUGetClassID(idata);
    mxComplexity   ccx = mxGPUGetComplexity(idata);
    const mwSize* dims = mxGPUGetDimensions(idata);
    int          depth = (int)dims[0];
    int          width = (int)dims[1];
    int             na = (int)dims[2];
    int              M = (int)parameter[0];
    int              N = (int)parameter[1];
    double           c = (double)parameter[2];
    double          fc = (double)parameter[3];
    int          averN = (int)parameter[4];
    int       pointNum = (int)parameter[5];
    int          Ledge = (int)parameter[6];
    int          Redge = (int)parameter[7];
    double unitDistance = (double)parameter[8];
    double     unitTime = (double)parameter[9];
    mwSize* dim1s=new mwSize[3];
    dim1s[0]=dims[0]-M+1;dim1s[1]=dims[1];dim1s[2]=dims[2]-N+1;
    mwSize* dim2s=new mwSize[2];
    dim2s[0]=dims[0]-M+1;dim2s[1]=dims[1];
    mxGPUArray* Im1 = mxGPUCreateGPUArray(3,dims,cid,ccx,MX_GPU_DO_NOT_INITIALIZE);
    mxGPUArray* Re1 = mxGPUCreateGPUArray(3,dims,cid,ccx,MX_GPU_DO_NOT_INITIALIZE);
    mxGPUArray* Im2 = mxGPUCreateGPUArray(3,dims,cid,ccx,MX_GPU_DO_NOT_INITIALIZE);
    mxGPUArray* Re2 = mxGPUCreateGPUArray(3,dims,cid,ccx,MX_GPU_DO_NOT_INITIALIZE);
    mxGPUArray* Im1_M = mxGPUCreateGPUArray(3,dims,cid,ccx,MX_GPU_DO_NOT_INITIALIZE);
    mxGPUArray* Re1_M = mxGPUCreateGPUArray(3,dims,cid,ccx,MX_GPU_DO_NOT_INITIALIZE);
    mxGPUArray* Im2_M = mxGPUCreateGPUArray(3,dims,cid,ccx,MX_GPU_DO_NOT_INITIALIZE);
    mxGPUArray* Re2_M = mxGPUCreateGPUArray(3,dims,cid,ccx,MX_GPU_DO_NOT_INITIALIZE);
    mxGPUArray* Im1_N = mxGPUCreateGPUArray(3,dims,cid,ccx,MX_GPU_DO_NOT_INITIALIZE);
    mxGPUArray* Re1_N = mxGPUCreateGPUArray(3,dims,cid,ccx,MX_GPU_DO_NOT_INITIALIZE);
    mxGPUArray* Im2_N = mxGPUCreateGPUArray(3,dims,cid,ccx,MX_GPU_DO_NOT_INITIALIZE);
    mxGPUArray* Re2_N = mxGPUCreateGPUArray(3,dims,cid,ccx,MX_GPU_DO_NOT_INITIALIZE);
    mxGPUArray* Dis = mxGPUCreateGPUArray(3,dim1s,cid,ccx,MX_GPU_DO_NOT_INITIALIZE);
    mxGPUArray* Dis_k2D = mxGPUCreateGPUArray(3,dim1s,cid,ccx,MX_GPU_DO_NOT_INITIALIZE);
    mxGPUArray* Dis_aver = mxGPUCreateGPUArray(3,dim1s,cid,ccx,MX_GPU_DO_NOT_INITIALIZE);
    mxGPUArray* minInd = mxGPUCreateGPUArray(2,dim2s,cid,ccx,MX_GPU_DO_NOT_INITIALIZE);
    mxGPUArray* Yang = mxGPUCreateGPUArray(2,dim2s,cid,ccx,MX_GPU_DO_NOT_INITIALIZE);
    double* d_Im1 = (double *)(mxGPUGetData(Im1));
    double* d_Re1 = (double *)(mxGPUGetData(Re1));
    double* d_Im2 = (double *)(mxGPUGetData(Im2));
    double* d_Re2 = (double *)(mxGPUGetData(Re2));
    double* d_Im1_M = (double *)(mxGPUGetData(Im1_M));
    double* d_Re1_M = (double *)(mxGPUGetData(Re1_M));
    double* d_Im2_M = (double *)(mxGPUGetData(Im2_M));
    double* d_Re2_M = (double *)(mxGPUGetData(Re2_M));
    double* d_Im1_N = (double *)(mxGPUGetData(Im1_N));
    double* d_Re1_N = (double *)(mxGPUGetData(Re1_N));
    double* d_Im2_N = (double *)(mxGPUGetData(Im2_N));
    double* d_Re2_N = (double *)(mxGPUGetData(Re2_N));
    double* d_Dis = (double *)(mxGPUGetData(Dis));
    double* d_Dis_k2D = (double *)(mxGPUGetData(Dis_k2D));
    double* d_Dis_aver = (double *)(mxGPUGetData(Dis_aver));
    double* d_minInd = (double *)(mxGPUGetData(minInd));
    double* d_Yang   = (double *)(mxGPUGetData(Yang));
    int blocksPerGrid = (mxGPUGetNumberOfElements (idata)-1) / threadsPerBlock + 1;
    calculateImRe<<<blocksPerGrid, threadsPerBlock>>>(d_idata,d_qdata,d_Im1,d_Re1,d_Im2,d_Re2,depth,width,na);
    caculate_dis_withM <<<2*blocksPerGrid, threadsPerBlock >>>(d_Im1, d_Re1, d_Im2, d_Re2, \
            d_Im1_M, d_Re1_M, d_Im2_M, d_Re2_M, depth, width, na, M, N);
    caculate_dis_withN <<<2*blocksPerGrid, threadsPerBlock >>>(d_Im1_N, d_Re1_N, d_Im2_N, d_Re2_N, \
            d_Im1_M, d_Re1_M, d_Im2_M, d_Re2_M, depth, width, na, M, N);
    caculate_dis_without_conv <<<blocksPerGrid, threadsPerBlock >>> (d_Dis, d_Im1_N, d_Re1_N, d_Im2_N, d_Re2_N, \
            depth, width, na, M, N, c, fc);
	caculate_dis_with_conv <<<blocksPerGrid, threadsPerBlock >>> (d_Dis_k2D, d_Dis, depth - M + 1, width, na - N + 1);
    caculate_ave_dis_in_axial <<<blocksPerGrid, threadsPerBlock >>> (d_Dis_aver, d_Dis_k2D, averN, depth - M + 1, width, na - N + 1);
    caculate_min_index <<<blocksPerGrid, threadsPerBlock >>> (d_minInd, d_Dis_aver, depth - M + 1, width, na - N + 1);
    caculate_Yang <<<blocksPerGrid, threadsPerBlock >>> (d_Yang, d_minInd, d_randpoint, pointNum, Ledge, Redge, depth - M + 1, unitDistance, unitTime);
    plhs[0] = mxGPUCreateMxArrayOnGPU(Yang);
    mxGPUDestroyGPUArray(idata);
    mxGPUDestroyGPUArray(qdata);
    mxGPUDestroyGPUArray(randpoint);
    mxGPUDestroyGPUArray(Im1);
    mxGPUDestroyGPUArray(Re1);
    mxGPUDestroyGPUArray(Im2);
    mxGPUDestroyGPUArray(Re2);
    mxGPUDestroyGPUArray(Im1_M);
    mxGPUDestroyGPUArray(Re1_M);
    mxGPUDestroyGPUArray(Im2_M);
    mxGPUDestroyGPUArray(Re2_M);
    mxGPUDestroyGPUArray(Im1_N);
    mxGPUDestroyGPUArray(Re1_N);
    mxGPUDestroyGPUArray(Im2_N);
    mxGPUDestroyGPUArray(Re2_N);
    mxGPUDestroyGPUArray(Dis);
    mxGPUDestroyGPUArray(Dis_k2D);
    mxGPUDestroyGPUArray(Dis_aver);
    delete []dim1s;
    delete []dim2s;
}