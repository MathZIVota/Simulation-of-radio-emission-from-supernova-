#include"Tube.hpp"
#include"AnalyticSolver.hpp"

int main(){
    cout << "#0\n";
    Tube tub(501,400);
    // p u rho
    tub.Initial_data({1, 0, 1},{0.1, 0, 0.125}, 1.4);
    cout << "#1\n";
    Solver solver;
    cout << "#2\n";
    solver.Solve(tub);
    cout << "#3\n";
    tub.UploadToFile();
    cout << "#4\n";

    cout << " -----------Finish_Calculation-------------- " << endl;
    return 0;
}

