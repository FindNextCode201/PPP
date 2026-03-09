#include "SWAS.h"
#include <Eigen/Dense>

// 曲线拟合函数,输出结果是从低阶开始：y=a0+a1x+a2x^2
Eigen::VectorXd curveFitting(const Eigen::VectorXd& x, const Eigen::VectorXd& y, int degree)
{
	int n = x.size();
	Eigen::MatrixXd A(n, degree + 1);

	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j <= degree; ++j)
		{
			A(i, j) = pow(x(i), j);
		}
	}

	Eigen::VectorXd result = A.householderQr().solve(y);
	return result;
}

//存最近的B2b信息
extern void RecordRAC(double *dR, double *dA, double *dC, double *dT,ssr_t b2bs, int ne, gtime_t ts)
{
	int i, j;
	if (ne < MAXPRENUM) {
		dR[ne] = b2bs.deph[0];
		dA[ne] = b2bs.deph[1];
		dC[ne] = b2bs.deph[2];
		dT[ne] = timediff(b2bs.t0[0], ts);		
	}
	else {
		for (i = 0; i < MAXPRENUM - 1; i++) {
			dR[i] = dR[i + 1];
			dA[i] = dA[i + 1];
			dC[i] = dC[i + 1];
			dT[i] = dT[i + 1];
		}
		dR[i] = b2bs.deph[0];
		dA[i] = b2bs.deph[1];
		dC[i] = b2bs.deph[2];
		dT[i] = timediff(b2bs.t0[0], ts);
	}
}

//存最近的B2b信息
extern void RecordCLK(double* dT, double* dClk, ssr_t b2bs, int ne, gtime_t ts)
{
	int i, j;
	if (ne < MAXCLKNUM) {
		dClk[ne] = b2bs.dclk[0];
		dT[ne] = timediff(b2bs.t0[1], ts);
	}
	else {
		for (i = 0; i < MAXCLKNUM - 1; i++) {
			dT[i] = dT[i + 1];
			dClk[i] = dClk[i + 1];
		}
		dT[i] = timediff(b2bs.t0[1], ts);
		dClk[i] = b2bs.dclk[0];
	}

}


//
extern double CalFit(const double *coff, double dxt) {
	double result = 0;
	for (int i = 0; i <= DEGREE; i++)
	{
		result += coff[i] * pow(dxt, i);
	}
	return result;
}
extern double CalFitClk(const double* coff, double dxt) {
	double result = 0;
	for (int i = 0; i <= DEGREE_CLK; i++)
	{
		result += coff[i] * pow(dxt, i);
	}
	return result;
}

extern void GetCoff(const double* arr_x, const double* arr_y, double *coff,int degree)
{
	//数组转矩阵
	Eigen::VectorXd x(MAXPRENUM), y(MAXPRENUM);
	for (int i = 0; i < MAXPRENUM; i++)
	{
		x[i] = arr_x[i];
		y[i] = arr_y[i];
	}
	Eigen::VectorXd result = curveFitting(x, y, degree);
	for (int i = 0; i < degree+1; i++) coff[i] = result[i];
	
}

extern void GetCoffClk(const double* arr_x, const double* arr_y, double* coff, int degree)
{
	//数组转矩阵
	Eigen::VectorXd x(MAXCLKNUM), y(MAXCLKNUM);
	for (int i = 0; i < MAXCLKNUM; i++)
	{
		x[i] = arr_x[i];
		y[i] = arr_y[i];
	}
	Eigen::VectorXd result = curveFitting(x, y, degree);
	for (int i = 0; i < degree + 1; i++) coff[i] = result[i];

}