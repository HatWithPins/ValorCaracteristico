#include "MetodoPotencias.h"
#include <fstream>
#include <string>

using namespace std;

MetodoPotencias::MetodoPotencias(int N, double lambda)
{
	_N = N;
	_lambda = lambda;
	x = new double[_N];

	for (int i = 0; i < _N; i++)
	{
		x[i] = i;
	}

	mu = 0.0;
	err = 0.0;
}

MetodoPotencias::~MetodoPotencias()
{
	delete[] x;
}


void MetodoPotencias::Solve()
{
	int M = 400;
	double TOL = 0.0000001;
	double* y = new double[_N];

	int p = 0;
	double x_p = x[p];
	double y_p;

	for (int i = 1; i < _N; i++)
	{
		if (abs(x_p) < abs(x[i]))
		{
			p = i;
			x_p = x[p];
		}
	}


	for (int i = 0; i < _N; i++)
	{
		x[i] = x[i]/x_p;
	}

	for (int k = 0; k < M; k++)
	{
		y[0] = (1 - 2 * _lambda) * x[0] + _lambda * x[1];
		y[_N - 1] = (1 - 2 * _lambda) * x[_N - 1] + _lambda * x[_N - 2];
		for (int j = 1; j <= _N - 2; j++)
		{
			y[j] = _lambda*x[j - 1] + (1 - 2 * _lambda) * x[j] + _lambda * x[j + 1];
		}

		p = 0;
		y_p = y[p];
		for (int j = 1; j < _N; j++)
		{
			if (abs(y_p) < abs(y[j]))
			{
				p = j;
				y_p = y[p];
			}
		}

		mu = y_p;
		if (abs(mu) < TOL) break;
		
		err = abs(x[0] - y[0] / y_p);
		for (int j = 1; j < _N; j++)
		{
			if (err < abs(x[j] - y[j] / y_p))
			{
				err = abs(x[j] - y[j] / y_p);
			}
		}

		for (int j = 0; j < _N; j++)
		{
			x[j] = y[j] / y_p;
		}

		if (err < TOL) break;
	}
}

double MetodoPotencias::GetSolution()
{
	return mu;
}

double MetodoPotencias::GetError()
{
	return err;
}

void MetodoPotencias::WriteSolution()
{
	ofstream file{ "results/metodo_Potencias-lambda-" + to_string(_lambda) + "-N-" + to_string(_N) + ".csv" };

	file << "x_0=(";

	for (int i = 0; i < _N; i++)
	{
		if (i < _N - 1)
		{
			file << to_string(x[i]) + ",";
		}
		else
		{
			file << to_string(x[i]) + ")\n";
		}
	}


	file << "lambda,err\n";

	for (int i = 0; i < _N; i++)
	{
		file << to_string(mu) + "," + to_string(err) + "\n";
	}

	file.close();
}
