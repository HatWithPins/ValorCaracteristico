class MetodoQR
{
    private:
        //Las variables del problema.
        //Variables of the problem.
        int _N;
        double _lambda;
        double* autovalores;
        double* b_n;

    public:
        MetodoQR(int N, double lambda);
        ~MetodoQR();
        void Solve();
        double* GetSolution();
        void WriteSolution();
};