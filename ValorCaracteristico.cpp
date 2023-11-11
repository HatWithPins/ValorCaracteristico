#include "ValorCaracteristico.h"

using namespace std;

int main(int argc, char** argv)
{
	string argument;
	size_t pos;
	size_t check;
	double lambda;
	int N;
	int expectedArguments = 3;
	vector<string> expectedArgumentsList = {"N=", "lambda="};
	vector<string> receivedArguments;

	//Se comprueba el número de parámetros recibidos. 
	//Checking the number of received parameters.
	if (argc != expectedArguments)
	{
		cout << "Error, expected " << expectedArguments - 1 << " as much, but received " << argc - 1 << endl;
		return -1;
	}

	//Se comprueban los nombres de los parámetros.
	//Checking if parameter names are correct.
	try
	{
		for (int i = 1; i < argc; i++)
		{
			argument = argv[i];
			check = argument.find(expectedArgumentsList[i - 1]);
			if (check < 0 || check > argument.size())
			{
				vector<string> exceptionVector = { argument, expectedArgumentsList[i - 1] };
				throw(exceptionVector);
			}
			pos = argument.find("=");
			receivedArguments.push_back(argument.substr(pos + 1));
		}
	}
	catch (vector<string> errorVector)
	{
		cout << "Error, expected " << errorVector[1] << "something, but received " << errorVector[0] << endl;
		return -1;
	}

	//Se castean los parámetros a sus respectivas variables.
	//Casting parameters to their variables.
	try
	{
		lambda = stod(receivedArguments[1]);
		N = stoi(receivedArguments[0]);
	}
	catch (const std::exception& e)
	{
		std::cerr << e.what() << '\n';
		return -1;
	}

	MetodoQR metodoQR(N, lambda);
	metodoQR.Solve();
	double* lambdas = metodoQR.GetSolution();
	for (int i = 0; i < N; i++)
	{
		cout << "Lambda_" << i << " = " << lambdas[i] << ", ";
	}

	MetodoPotencias metodoPotencias(N, lambda);
	metodoPotencias.Solve();
	cout << "\n" << "mu = " << metodoPotencias.GetSolution() << ", err = " << metodoPotencias.GetError() << endl;

	return 0;
}
