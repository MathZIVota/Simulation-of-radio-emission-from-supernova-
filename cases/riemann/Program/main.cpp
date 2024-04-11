#include<iostream>
#include<cmath>
#include<vector>
#include <fstream>
using namespace std;

class Tube{
public:
    class Geometry
    {
    public:
        double xl, xr, xi; 
        Geometry(double in_xl, double in_xr, double in_xi){
            xl = in_xl;
            xr = in_xr;
            xi = in_xi;
        }
    };

    class State
    {
    public:
        double pressure, density, velocity;
        State(double p, double rho, double v){
            velocity = v;
            density = rho;
            pressure = p;
        }
    };
    
    Geometry geometry{-1,1,0};
    State left_state{0,0,0}, right_state{0,0,0};
    double gamma = 1.4;
    Tube(){};

    Tube(Geometry geom, State l_state, State r_state){
        geometry.xl = geom.xl;
        geometry.xr = geom.xr;
        geometry.xi = geom.xi;
        left_state.density = l_state.density;
        left_state.velocity = l_state.velocity;
        left_state.pressure = l_state.pressure;
        right_state.density = r_state.density;
        right_state.velocity = r_state.velocity;
        right_state.pressure = r_state.pressure;
    };

    Tube& operator=(const Tube& t){
        if(&t != this)
        {
            geometry.xi = t.geometry.xi;
            geometry.xl = t.geometry.xl;
            geometry.xr = t.geometry.xr;

            left_state.density = t.left_state.density;
            left_state.pressure =  t.left_state.pressure;
            left_state.velocity =  t.left_state.velocity;

            right_state.density =  t.right_state.density;
            right_state.pressure =  t.right_state.pressure;
            right_state.velocity =  t.right_state.velocity;
        }
        return *this;
    }

    void print() {
        std::cout << "Geometry: " << geometry.xl << " " << geometry.xi << " " << geometry.xr << std::endl;
        std::cout << "Left: " << left_state.pressure << " " << left_state.density  << " " << left_state.velocity << "\n";
        std::cout << "Right: " << right_state.pressure << " " << right_state.density << " " << right_state.velocity << "\n";
    }

};

class Solver{
public:
    Tube tube;

    Solver(Tube t){
        tube = t;
    }
    


    double sound_speed(){
        return sqrt(tube.gamma * tube.left_state.pressure/tube.left_state.density);
    }
};
// sound speed: 1.1832159566199232
// double a0 = 1.1832159566199232;

class Сond
{
public:
    double b[3];
    void Set(double rho, double rho_u, double rho_E){
        b[0] = rho;
        b[1] = rho_u;
        b[2] = rho_E;
    }
    Сond& operator=(const Сond& data) {
        if(&data != this){
            for(int i=0; i<3; i++) b[i] = data.b[i];
        }
        return *this;
    }
    Сond operator+(const Сond& data) const{
        //Сond tmp;
        //tmp.Set();
        //for(size_t i=0; i<3; i++) tmp.b[i] = b[i]+data[i];
        return Сond{data.b[0]+b[0], data.b[1]+b[1],data.b[2]+b[2]};
    }  

    Сond operator*(const double a){
        return Сond{a*b[0], a*b[1], a*b[2]};
    }

    void print(){
        cout << b[0] << " " << b[1] << " " << b[2] << endl;
    }
};

void Func(vector<Сond> A){
    cout << "Func: " << A.size() << endl;
    for(int i=0; i<A.size(); i++){
        A[i].print();
        cout << endl;
    }   
}

void Func1(vector<Сond> &A){
    cout << "Func1: " << A.size() << endl;
    for(int i=0; i<A.size(); i++){
        A[i].Set(i/10,i/2,i/5);
    }   
}

int main(){
    Tube tube({-1,1,0},{1,1,0},{0.125, 0.1, 0});
    tube.print();
    Solver solver(tube);
    int Nx = 2000;
    double L = 1, t_max = 1, c = 1;
    double x_l = -1, x_m = 0, x_r = 1;
    double dx = (x_r-x_l)/(Nx-1);
    int Nt = 350;
    double dt = t_max/Nt;
    
    //double dt = T/Nt;
    //int Nt = round(T/dt) + 1;
    
    vector<double> t,x;
    cout << "t("<< t_max << " "<< Nt << " " << dt <<"):";

    // Разбиение времени
    for (int i=0; i<Nt; i++){
        t.push_back((double)i*dt);
        cout << i*dt << " ";
    }
    cout << endl;

    //Разбиение координаты
    cout << "x("<< L << " "<< Nx << " " << dx <<"):";
    for (double i=0; i<Nx; i++) {
        double xi = x_l+i*dx;
        x.push_back(xi);
        cout << xi << " ";
    }
    cout << endl << endl;
    

    vector<vector<double>> p,u,rho,T,e,E,h,S;
    double p4 = 1, p1 = 0.1, u4 = 0, u1 = 0, rho4 = 1, rho1 = 0.1, T1 = 1.4, T4 = 1;
    double gamma = 1.4;
    
    //cout << "sound_speed:" <<  a0 << endl;
    
    cout << "---" << endl;
    double D = 0, M = 0, u_2 = 0, p_2 = 0, rho_2 = 0, u_3=0, p_3=0, rho_3=0;
    double c1 = sqrt(gamma*p1/rho1);
    cout << "x size:"<< Nx << " "<< endl;
    double eps = 1e-15;
    D = rho4/rho1*(p1-p4)/(rho1-rho4);
    // cout << D << " - " << c1  << " /" << x[Nx/2] << endl;
    
    int i0 = 0, i1=Nx/2, i2=Nx/2, i3=Nx/2, i4=Nx/2, i5 = Nx-1, i_initial = 200;
    p_3 = p4;
    rho_3 = rho4;
    //cout << t[200] << endl;
    cout << " -----------Begin_Calculation-------------- " << endl;
    for (size_t k = 0; k < Nt; k++)
    {
        p.push_back(x);
        u.push_back(x);
        rho.push_back(x);
        //e.push_back(x);
        //E.push_back(x);
        //h.push_back(x);
        T.push_back(x);
        S.push_back(x);
        //M.push_back(x);
        for(size_t i = 0; i < Nx; i++){
            p[k][i] = p4;
            u[k][i] = u4;
            rho[k][i] = rho4;
            T[k][i] = T4;
            if(x[i] > 0) {
                p[k][i] = p1;
                u[k][i] = u1;
                rho[k][i] = rho1;
                T[k][i] = T1;
            } 
            M = 1.6;

        }

        if(k > 0){
            //cout << "! " << k << "/ " << p_3 << ", " << p1 << "/ " << rho_3 << ", " << rho1 << "/" << endl;//<< p[0][199] << "," << p[0][200] << "," << p[0][201]  << endl;
            if(k > 1) D = rho_3*u_3/(rho_3-rho1);//rho_3/rho1*(p_3-p1)/(rho_3-rho1); ////////
            else D = rho4/rho1*(p1-p4)/(rho1-rho4);//M*c1;
            
            cout << "! " << D << " " << D*dt << endl;
            cout << "! " << c1 << " " << c1*dt << endl;

            int j = i3;
            while(x[j] < x[i3]+c1*dt) j++;
            i3 = j;
            j = i4;
            while(x[j] < x[i4]+D*dt) j++;            
            i4 = j;
            

            if(i4 > Nx-1) i4 = Nx-1;
            if(i3 > Nx-1) i3 = Nx-1;

            //p_3 = p1*pow(1-(gamma-1)*u_3/(2*c1),2*gamma/(gamma-1));
            p_3 = p1*(2*gamma*M*M-(gamma-1))/(gamma+1);
            rho_3 = rho1*(gamma+1)*M*M/((gamma-1)*M*M+2);
            //rho_3 = rho1*((gamma-1)*p4+(gamma+1)*p_3)/((gamma-1)*p_3+(gamma+1)*p4);
           // u_3 = c1*2*(M-1/M)/(gamma+1);
             u_3 = sqrt(2*pow((p_3-p1),2)/(rho1*((gamma-1)*p1+(gamma+1)*p_3)));
            //u_3 = 2*c1*(M-1/M)/(gamma+1);
            //u_3 = 1 - 2*c1/(gamma-1)*pow(p_3/p4, (gamma-1)/2/gamma);
            //u_3 = sqrt(2/rho1)*(p_3-p1)/sqrt((gamma-1)*p1+(gamma+1)*p_3);

            p_2 = p_3;
            u_2 = u_3;
            rho_2 = rho4*pow(p_3/p4,1/gamma);

            cout << "|" << i2 << " " << i3 << " " << i4  << " " << u_3<< endl;
            cout << x[i2] << " " << x[i3] << endl;
            for(int i=i2; i<i3;i++){
                rho[k][i] = rho_2;
                p[k][i] = p_2;
                u[k][i] = u_2;
            }
            cout << x[i3] << " " << x[i4] << endl;
            for(int i=i3; i<i4 ; i++){
                rho[k][i] = rho_3;
                p[k][i] = p_3;
                u[k][i] = u_3;
            }
            
            M = D/c1;
        }
        for(int i=0; i<Nx; i++){
            S[k][i] = log(p[k][i]/pow(rho[k][i], gamma));//+S[0][i];
        }
        cout << "M = " << M  << " " << p_3/p4 << endl;
    }
    //cout << dx << " " << Nx << " | " << dt << " " << Nt << endl;
    cout << " -----------End_Calculation-------------- " << endl;

    //**********************************************************//
     //*********************************************************//
    ofstream out1,out2;
    out2.open("x_t_diagram.txt");
    if(out2.is_open()){
        out2 << c1 << " " << D << endl;
        out2 <<  0 << " " << x[Nx-1] << endl;
    }
    out2.close();

    out1.open("result_rho.txt");
    if(out1.is_open()){
        for (int i=0; i<Nx; i++){
            out1 << x[i] << " ";
        }
        out1 << endl;
        for (int i=0; i<Nt; i++){
            out1 << t[i] << " ";
            for(int j=0; j<Nx; j++) out1  << rho[i][j] << " ";
            out1 << endl;
        }
    }
    out1.close();
    cout << solver.sound_speed() << endl;

    
    out1.open("result_p.txt");
    if(out1.is_open()){
        for (int i=0; i<Nx; i++){
            out1 << x[i] << " ";
        }
        out1 << endl;
        for (int i=0; i<Nt; i++){
            out1 << t[i] << " ";
            for(int j=0; j<Nx; j++) out1  << p[i][j] << " ";
            out1 << endl;
        }
    }
    out1.close();

    out1.open("result_u.txt");
    if(out1.is_open()){
        for (int i=0; i<Nx; i++){
            out1 << x[i] << " ";
        }
        out1 << endl;
        for (int i=0; i<Nt; i++){
            out1 << t[i] << " ";
            for(int j=0; j<Nx; j++) out1  << u[i][j] << " ";
            out1 << endl;
        }
    }
    out1.close();
     out1.open("result_S.txt");
    if(out1.is_open()){
        for (int i=0; i<Nx; i++){
            out1 << x[i] << " ";
        }
        out1 << endl;
        for (int i=0; i<Nt; i++){
            out1 << t[i] << " ";
            for(int j=0; j<Nx; j++) out1  << S[i][j] << " ";
            out1 << endl;
        }
    }
    out1.close();


   // cout << " -----------Begin_Plotting------------- " << endl;
   // system("plotting.exe");
   // cout << " -----------End_Plotting------------- " << endl;
    return 0;
}

