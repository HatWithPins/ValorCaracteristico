#include "MetodoQR.h"
#include <fstream>
#include <string>

using namespace std;

MetodoQR::MetodoQR(int N, double lambda)
{
	_N = N;
	_lambda = lambda;
	autovalores = new double[N];
}

MetodoQR::~MetodoQR()
{
	delete[] autovalores;
}

void MetodoQR::Solve()
{
	int M = 400;
	int n = _N - 1;
	double TOL = 0.0000001;
	double SHIFT = 0.0;
	double b;
	double c;
	double d;
	double mu_1;
	double mu_2;
	double sigma;

	double* x = new double[_N];
	double* y = new double[_N];
	double* z = new double[_N];
	double* d_n = new double[_N];
	double* sigma_n = new double[_N];
	double* c_n = new double[_N];
	double* q = new double[_N];

	double* b_n = new double[_N];
	double* a = new double[_N];

	for (int i = 0; i < _N; i++)
	{
		b_n[i] = _lambda;
		a[i] = 1 - 2 * _lambda;
	}

	for (int k = 0; k < M; k++)
	{
		if (b_n[n] <= TOL)
		{
			autovalores[n] = a[n] + SHIFT;
			n--;
		}
		if (b_n[1] <= TOL)
		{
			autovalores[0] = a[0] + SHIFT;
			n--;
			for (int j = 0; j < n; j++)
			{
				a[j] = a[j + 1];
				b_n[j+1] = b_n[j + 2];
			}
		}

		if (n == 0)
		{
			autovalores[0] = a[0] + SHIFT;
			break;
		}

		b = -(a[n - 1] + a[n]);
		c = a[n] * a[n - 1] - b_n[n] * b_n[n];
		d = sqrt(b * b - 4 * c);

		if (b > 0)
		{
			mu_1 = -2 * c / (b + d);
			mu_2 = -(b + d) / 2;
		}
		else
		{
			mu_1 = (d - b) / 2;
			mu_2 = 2 * c / (d - b);
		}

		if (n == 1)
		{
			autovalores[0] = mu_1 + SHIFT;
			autovalores[1] = mu_2 + SHIFT;
			break;
		}

		sigma = mu_1 * (abs(mu_1 - a[n]) < abs(mu_2 - a[n])) + mu_2 * (abs(mu_1 - a[n]) > abs(mu_2 - a[n]));
		SHIFT += sigma;

		for (int j = 0; j <= n; j++)
		{
			d_n[j] = a[j] - sigma;
		}

		x[0] = d_n[0];
		y[0] = b_n[1];

		for (int j = 1; j <= n; j++)
		{
			z[j-1] = sqrt(x[j-1] * x[j-1] + b_n[j] * b_n[j]);
			c_n[j] = x[j-1] / z[j-1];
			sigma_n[j] = b_n[j] / z[j-1];
			q[j-1] = c_n[j] * y[j-1] + sigma_n[j] * d_n[j];
			x[j] = -sigma_n[j] * y[j - 1] + c_n[j] * d_n[j];

			if (j != n)
			{
				y[j] = c_n[j] * b_n[j + 1];
			}
		}

		z[n] = x[n];
		a[0] = sigma_n[1] * q[0] + c_n[1] * z[0];
		b_n[1] = sigma_n[1] * z[1];

		for (int j = 1; j <= n ; j++)
		{
			a[j] = sigma_n[j + 1] * q[j] + c_n[j]*c_n[j + 1] * z[j];
			b_n[j + 1] = sigma_n[j + 1] * z[j + 1];
		}

		a[n] = c_n[n] * z[n];
	}
}

double* MetodoQR::GetSolution()
{
	return autovalores;
}

void MetodoQR::WriteSolution()
{

}