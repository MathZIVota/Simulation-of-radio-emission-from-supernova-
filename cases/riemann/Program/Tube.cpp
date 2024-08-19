#include"Tube.h"

void Tube::Initial_data(vector<double> data_left = {1, 0, 1}, vector<double> data_right = {0.125, 0, 0.1}, double Gamma = 1.4){
    //data = {p, u, rho}
    gamma = Gamma;
    for(int i=0; i<Nx; i++){
        if(x[i]<=0){
            p[0][i] = data_left[0];
            u[0][i] = data_left[1];
            rho[0][i] = data_left[2];
        }else{
            p[0][i] = data_right[0];
            u[0][i] = data_right[1];
            rho[0][i] = data_right[2];
        }
        double R = 1; //8.31;
        T[0][i] = p[0][i]/rho[0][i]/R;
        S[0][i] =log(p[0][i]/pow(rho[0][i],gamma));
    }
    c_right = sqrt(gamma*data_right[0]/data_right[2]);
    c_left = sqrt(gamma*data_left[0]/data_left[2]);
}

void Tube::addNullVector(){
    vector<double> tmp(Nx,{0});
    p.push_back(tmp);
    u.push_back(tmp);
    rho.push_back(tmp);
    T.push_back(tmp);
    E.push_back(tmp);
    h.push_back(tmp);
    S.push_back(tmp);
}


Tube::Tube(int num_x, int num_t):
        x_l(-1), x_m(0), x_r(1),
        t_max(1)
{
    u0 = 1.254032;
    Nx = num_x;
    Nt = num_t;
    dt = t_max/Nt;
    dx = (x_r-x_l)/(Nx-1);

    for (int i=0; i<Nt; i++){
        t.push_back((double)i*dt);
    }
    for (double i=0; i<Nx; i++) {
        double xi = x_l+i*dx;
        x.push_back(xi);
    }
    index = vector<int>(4, Nx/2);
    addNullVector();
    Initial_data();
}

void Tube::getInitial_P(double &p4, double &p1){
    p4 = p[0][Nx/2];
    p1 = p[0][Nx/2+1];
}

void Tube::getInitial_Rho(double &rho4, double &rho1){
    rho4 = rho[0][Nx/2];
    rho1 = rho[0][Nx/2+1];
}

void Tube::getInitial_U(double &u4, double &u1){
    u4 = u[0][Nx/2];
    u1 = u[0][Nx/2+1];
}