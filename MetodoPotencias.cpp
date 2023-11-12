#include "MetodoPotencias.h"
#include <fstream>
#include <string>

using namespace std;

//Constructor
MetodoPotencias::MetodoPotencias(int N, double lambda)
{
	_N = N;
	_lambda = lambda;
	x = new double[_N];

	//Vector x de prueba. No hay un método para elegir uno concreto. Lo suyo sería poder especificar esto sin tener que recompilar.
	//Test x vector. There is no criteria to pick a particular one. Ideally, this should be specified without recompiling.
	for (int i = 0; i < _N; i++)
	{
		x[i] = i + 1;
	}

	mu = 0.0;
	err = 0.0;
}

//Destructor
MetodoPotencias::~MetodoPotencias()
{
	delete[] x;
}


//Implementación del método.
//Method implementation.
void MetodoPotencias::Solve()
{
	int M = 1000;
	double TOL = 0.00000001;
	double* y = new double[_N];

	int p = 0;
	double x_p = x[p];
	double y_p;

	//Cálculo de la norma infinita de x.
	//Computation of infinity norm of x.
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

	//Aplicación del método.
	//Method application.
	for (int k = 0; k < M; k++)
	{
		//y = Ax
		y[0] = (1 - 2 * _lambda) * x[0] + _lambda * x[1];
		y[_N - 1] = (1 - 2 * _lambda) * x[_N - 1] + _lambda * x[_N - 2];
		for (int j = 1; j <= _N - 2; j++)
		{
			y[j] = _lambda*x[j - 1] + (1 - 2 * _lambda) * x[j] + _lambda * x[j + 1];
		}

		//Cálculo de la norma infinita de y.
		//Computation of infinity norm of y.
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
		//Si el autovalor está por debajo de la tolerancia, asumimos que es 0.0 y acabamos.
		//If eigenvalue is below tolerance, we assume it's 0.0 and end.
		if (abs(mu) < TOL) {
			mu = 0.0;
			break;
		}
		
		//Cálculo del error estimado.
		//Estimating error.
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

//Para devolver la solución.
//To return solution.
double MetodoPotencias::GetSolution()
{
	return mu;
}

//Para devolver el error estimado.
//To return estimated error.
double MetodoPotencias::GetError()
{
	return err;
}

//Para escribir el resultado y el vector utilizado a un fichero.
//To write results and used vector to a file.
void MetodoPotencias::WriteSolution()
{
	ofstream file{ "results/metodo_Potencias-lambda-" + to_string(_lambda) + "-N-" + to_string(_N) + ".txt" };
	
	file << "$\textbf{x}_0=($";

	for (int i = 0; i < _N; i++)
	{
		if (i < _N - 1)
		{
			file << to_string(x[i]) + ", ";
		}
		else
		{
			file << to_string(x[i]) + "$)$\n";
		}
	}


	file << "lambda,err\n";
	file << to_string(mu) + "," + to_string(err) + "\n";

	file.close();
}
