#pragma once
#include<iostream>
#include<cmath>
#include<vector>
#include <fstream>
using namespace std;

class Tube{
public:
    int Nx;
    double x_l, x_m, x_r;
    vector<int> index = vector<int>(4,0);
    double dx;
    
    int Nt;
    double t_max;
    double dt;

    vector<double> t, x;
    vector<vector<double>> p,u,rho,T,E,h,S; //параметры газа

    double D, M; 
    double c_left,c_right;
    double c3, Uc; //скорости звука слева и справа
    double u0; //3006317067;
    double gamma;//cкорость ударной волны
    
    void Initial_data(vector<double>, vector<double>, double);

    void addNullVector();

    Tube(int, int);

    void getInitial_P(double &, double &);

    void getInitial_Rho(double &, double &);

    void getInitial_U(double &, double &);

    void UploadToFile();
};
