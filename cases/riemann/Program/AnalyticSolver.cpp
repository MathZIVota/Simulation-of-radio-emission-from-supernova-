#include"AnalyticSolver.hpp"
#include"Nonlinear_sol_Newton.hpp"

void Solver::Calc_1(size_t k, size_t i, Tube &tube){
    tube.u[k][i] = 2/(tube.gamma+1)*(tube.c_left+tube.x[i]/tube.t[k]); 
    tube.rho[k][i] = rho_left*pow(1-(tube.gamma-1)*tube.u[k][i]/(2*tube.c_left), 2/(tube.gamma-1)); 
    tube.p[k][i] = p_left*pow(1-(tube.gamma-1)*tube.u[k][i]/(2*tube.c_left), 2*tube.gamma/(tube.gamma-1));
}

void Solver::Calc_2(size_t k, size_t i, Tube &tube){
    tube.rho[k][i] = rho_left*pow(p_3/p_left, 1/tube.gamma);
    tube.p[k][i] = p_3;
    tube.u[k][i] = u_3;
}

void Solver::Calc_u_p(){
    tuple<double, double> result = solve_nonlinear_shock();
    this->u_3 = get<0>(result); // Получаем u из кортежа
    this->p_3 = get<1>(result); // Получаем p из кортежа
    cout << u_3 << " " << p_3 << endl;
    u_3 = 0.9275;
    p_3 = 0.3031;
}

void Solver::Calc_3(size_t k, size_t i, Tube &tube){
    u_3 = 0.9275;
    p_3 = 0.3031;
    
    tube.rho[k][i] =  rho_right*((tube.gamma-1)*p_right+(tube.gamma+1)*p_3)/((tube.gamma+1)*p_right+(tube.gamma-1)*p_3);
    tube.p[k][i] = p_3;
    tube.u[k][i] = u_3;
}

double Solver::Calc_Temperature(size_t k, size_t i, Tube &tube){
    return tube.p[k][i]/tube.rho[k][i];    
}

double Solver::Calc_Entropy(size_t k, size_t i, Tube &tube){
    return log(tube.p[k][i]/pow(tube.rho[k][i],tube.gamma));
}
/*
double Solver::Calc_Energy(size_t k, size_t i, Tube &tube){
    return 0;
}
*/
void Solver::Solve_inTimeMoment(size_t k, Tube &tube){ // k = time iteration 
    //double c3 = sqrt(tube.gamma*p_3/rho_3);
    double a3 = tube.c_left - (tube.gamma+1)*u_3/2; 

    int j = 0;
    while(tube.x[j] < -tube.c_left*tube.t[k]) j++;
    tube.index[0] = j;

    while(tube.x[j] < -a3*tube.t[k]) j++;
    tube.index[1] = j;


    while(tube.x[j] < u_3*tube.t[k]) j++;
    tube.index[2] = j;

    while(tube.x[j] < tube.D*tube.t[k]) j++;            
    tube.index[3] = j;

    if(tube.index[3] >= Nx) tube.index[3] = Nx-1;
    if(tube.index[2] >= Nx) tube.index[2] = Nx-1;
    if(tube.index[1] >= Nx) tube.index[1] = Nx-1;
    if(tube.index[0] < 0) tube.index[0] = 0;
    
    
    for(int i=0; i<Nx; i++){
        if(i < tube.index[0]){
            tube.rho[k][i] = rho_left;
            tube.p[k][i] = p_left;
            tube.u[k][i] = u_left;
        }
        else if(i < tube.index[1]){
            Calc_1(k, i, tube);
        }else if(i < tube.index[2]){
            Calc_2(k, i, tube);
        }else if(i < tube.index[3]){
            Calc_3(k, i, tube);
            tube.D = (p_3-p_right)/(rho_right*u_3);
        }else{
            tube.rho[k][i] = rho_right;
            tube.p[k][i] = p_right;
            tube.u[k][i] = u_right;
        }
        tube.T[k][i] = Calc_Temperature(k,i,tube);
        tube.S[k][i] = Calc_Entropy(k,i,tube);
    }
}

void Solver::Solve(Tube &tube){
    Nx = tube.Nx;
    Nt = tube.Nt;

    tube.getInitial_P(p_left, p_right);
    tube.getInitial_Rho(rho_left, rho_right);
    tube.getInitial_U(u_left, u_right);
    
    tube.M = tube.D/tube.c_right;
    
    for(size_t k = 1; k < (size_t)Nt; k++){
        cout << k;
        tube.addNullVector();
        Solve_inTimeMoment(k, tube);
        cout << "$" << endl;
    }
}
    