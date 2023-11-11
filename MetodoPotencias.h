class MetodoPotencias
{
    private:
        //Las variables del problema.
        //Variables of the problem.
        int _N;
        double _lambda;
        double* x;
        double mu;
        double err;

    public:
        MetodoPotencias(int N, double lambda);
        ~MetodoPotencias();
        void Solve();
        double GetSolution();
        double GetError();
        void WriteSolution();
};