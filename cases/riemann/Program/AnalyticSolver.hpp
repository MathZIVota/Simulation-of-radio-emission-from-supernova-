#pragma once
#include<iostream>
#include<cmath>
#include<vector>
#include <fstream>
using namespace std;
#include "Tube.hpp"

class Solver{
public:
    int Nx, Nt;
    double p_right=0, rho_right=0, u_right=0;
    double p_3=0, rho_3=0, u_3=0; 
    double p_2=0, rho_2=0, u_2=0;
    double p_1=0, rho_1=0, u_1=0;
    double p_left=0, rho_left=0, u_left=0;
    

    void Solve(Tube &);
    
    void Calc_1(size_t, size_t, Tube &);
    
    void Calc_2(size_t, size_t, Tube &);

    void Calc_3(size_t, size_t, Tube &);

    void Calc_u_p();

    double Calc_Temperature(size_t, size_t, Tube &);

    double Calc_Entropy(size_t, size_t, Tube &);

    double Calc_Energy(size_t, size_t, Tube &);

    void Solve_inTimeMoment(size_t, Tube &);
};