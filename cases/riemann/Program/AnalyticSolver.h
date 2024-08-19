#pragma once
#include<iostream>
#include<cmath>
#include<vector>
#include <fstream>
using namespace std;
#include "Tube.h"

class Solver{
public:
    int Nx, Nt;
    double p_right, rho_right, u_right;
    double p_3, rho_3, u_3=0; 
    double p_2, rho_2, u_2=0;
    double p_1, rho_1, u_1=0;
    double p_left, rho_left, u_left;
    

    void Solve(Tube &);
    
    void Calc_1(size_t, size_t, Tube &);
    
    void Calc_2(size_t, size_t, Tube &);

    void Calc_3(size_t, size_t, Tube &);

    double Calc_Temperature(size_t, size_t, Tube &);

    double Calc_Entropy(size_t, size_t, Tube &);

    double Calc_Energy(size_t, size_t, Tube &);

    void Solve_inTimeMoment(size_t, Tube &);
};