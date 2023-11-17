#include "MetodoQR.h"
#include <fstream>
#include <string>
#include <cmath>

using namespace std;

//Constructor
MetodoQR::MetodoQR(int N, double lambda)
{
	_N = N;
	_lambda = lambda;
	autovalores = new double[N];
	//Aunque b_n en realidad tendría 1 elemento menos, para simplificar el trabajo con los índices, se crea de tamaño
	//N y se ignora el índice 0.
	//Although b_n actually has one element less, in order to ease indexing, b_n has size N but we ignore the first element.
	b_n = new double[_N];
}

//Destructor
MetodoQR::~MetodoQR()
{
	delete[] autovalores;
}

//Implementación del método QR.
//QR method implementation.
void MetodoQR::Solve()
{
	int M = 1000;
	int n = _N - 1;
	double TOL = 0.00000001;
	double SHIFT = 0.0;
	double b;
	double c;
	double d;
	double mu_1;
	double mu_2;
	double sigma;

	//Variables auxiliares para calcular las matrices Q y R.
	//Helping variables to calcultate Q and R matrices.
	double* x = new double[_N];
	double* y = new double[_N];
	double* z = new double[_N];
	double* d_n = new double[_N];
	double* sigma_n = new double[_N];
	double* c_n = new double[_N];
	double* q = new double[_N];

	double* a = new double[_N];

	//Primera asignación.
	//Initial assigment.
	for (int i = 0; i < _N; i++)
	{
		b_n[i] = _lambda;
		a[i] = 1 - 2 * _lambda;
	}

	for (int k = 0; k < M; k++)
	{
		//Comprobando si se ha resuelto el problema para el último b_n.
		//Checking if problem is solved for latest b_n.
		if (std::abs(b_n[n]) <= TOL)
		{
			autovalores[n] = a[n] + SHIFT;
			n--;
		}
		if (std::abs(b_n[1]) <= TOL)
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

		//Seguimos resolviendo para b_n.
		// Solving for b_n. 
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

		//Comprobamos si el problema está resuelto por completo.
		//Checking if problem is completely solved.
		if (n == 1)
		{
			autovalores[0] = mu_1 + SHIFT;
			autovalores[1] = mu_2 + SHIFT;
			break;
		}

		//Los siguientes pasos son el cálculo de las matrices Q y R y su aplicación.
		//Next steps are for calculating Q and R matrices and their application.

		//Cálculo del desplazamiento.
		//Fancy way to calculate SHIFT.
		sigma = mu_1 * (std::abs(mu_1 - a[n]) < std::abs(mu_2 - a[n])) + mu_2 * (std::abs(mu_1 - a[n]) > std::abs(mu_2 - a[n]));
		SHIFT += sigma;

		for (int j = 0; j <= n; j++)
		{
			d_n[j] = a[j] - sigma;
		}

		x[0] = d_n[0];
		y[0] = b_n[1];

		//Cálculo de PA y R.
		//PA and R calculation.
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


		//Cálculo de la siguiente iteración de A.
		//Next iteration of A.
		for (int j = 1; j <= n ; j++)
		{
			a[j] = sigma_n[j + 1] * q[j] + c_n[j]*c_n[j + 1] * z[j];
			b_n[j + 1] = sigma_n[j + 1] * z[j + 1];
		}

		a[n] = c_n[n] * z[n];
	}
}

//Para devolver los autovalores.
//To return eigenvalues.
double* MetodoQR::GetSolution()
{
	return autovalores;
}

//Para escribir la solución en un CSV.
//To write solution to CSV.
void MetodoQR::WriteSolution()
{
	ofstream latex_file{ "results/metodo_QR-lambda-" + to_string(_lambda) + "-N-" + to_string(_N) + ".txt" };
	ofstream csv_file{ "results/metodo_QR-lambda-" + to_string(_lambda) + "-N-" + to_string(_N) + ".csv" };

	latex_file << "$\\lambda=$";
	csv_file << "lambda,b_n\n";

	for (int i = 0; i < _N; i++)
	{
		csv_file << to_string(autovalores[i]) + "," + to_string(b_n[i]) +"\n";
		if (i != _N - 1)
		{
			latex_file << to_string(autovalores[i]) + ", ";
		}
		else
		{
			latex_file << to_string(autovalores[i]) + "\n";
		}
		
	}

	csv_file.close();
	latex_file << "$b_{n}=$";
	for (int i = 1; i < _N; i++)
	{
		if (i != _N - 1)
		{
			latex_file << to_string(b_n[i]) + ", ";
		}
		else
		{
			latex_file << to_string(b_n[i]) + "\n";
		}

	}

	latex_file.close();
}