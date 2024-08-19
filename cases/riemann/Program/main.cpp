#include"Tube.h"
#include"AnalyticSolver.h"

int main(){
    cout << "#0\n";
    Tube tub(501,400);
    
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

