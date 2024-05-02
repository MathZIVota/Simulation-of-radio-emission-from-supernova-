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

class TUBE{
    private:
    void Set_vector_Value(vector<vector<double>> &a,double a4, double a1, int j){
        vector<double> tmp(Nx);
        a.push_back(tmp);
        for(size_t i=0; i<Nx; i++){
            a[j][i] = a4;
            if(x[i] > 0) a[j][i] = a1;
        }
    }
    
    class Geometry
    {
    public:
        double xl, xr, xi;
        double x1, x2, x3, x4; 
        Geometry(double in_xl=-1, double in_xr=1, double in_xi=0){
            xl = in_xl;
            xr = in_xr;
            xi = in_xi;
            x1 = xi;
            x2 = xi;
            x3 = xi;
            x4 = xi;
        }
    };

    public:
    int Nx = 0;
    int Nt = 0;
    vector<double> x, t;
    vector<vector<double>> p,u,rho,T,e,E,h,S;
    Geometry geom;
    void Set_Geometry(int N1 = 501, int N2 = 300, double xl = -1, double xr = 1){
        Nx = N1;
        Nt = N2;
        geom.xi = 0;
        geom.xl = xl;
        geom.xr = xr;

        double dx = (xr-xl)/(Nx-1);
        for (double i=0; i<Nx; i++) {
            double xi = xl+i*dx;
            x.push_back(xi);
        }

        double dt = 1./(Nt);
        for (int i=0; i<Nt; i++){
            t.push_back((double)i*dt);
        }
    };
    
    void Set_Initial_Values(vector<double> initial_value = {1,1,0,0.1,0.125,0}){
        Set_vector_Value(p, initial_value[0], initial_value[3], 0);
        Set_vector_Value(rho, initial_value[1], initial_value[4], 0);
        Set_vector_Value(u, initial_value[2], initial_value[5], 0);
    };

    void Print(){
        for(int i=0; i<Nx; i++){
            //printf("", x[i]);
            if(i % 50 == 0) printf("%lf:%lf %lf %lf", p[0][i], rho[0][i], u[0][i]);
            cout << endl;
        }
    };

};


class SOLVER{
public:
    int Nt;
    vector<double> t;
    TUBE tube;
    // i1 - передний фронт волны разрежения, i2 - задний, i3 - контактный разрыв, i4 - ударная волна
    int i1,i2,i3,i4; 

    SOLVER(){
        Nt = 200;
        double dt = 1./Nt;
        for (int i=0; i<Nt; i++){
            t.push_back((double)i*dt);
        }
        tube.Set_Geometry();
        tube.Set_Initial_Values();
        i1 = tube.geom.xi;
        i2 = tube.geom.xi;
        i3 = tube.geom.xi;
        i4 = tube.geom.xi;
    }
};

int main(){
    int Nx = 301;
    int Nt = 350;

    double t_max = 1;
    double x_l = -1, x_m = 0, x_r = 1;
    double dx = (x_r-x_l)/(Nx-1);
    
    double dt = t_max/Nt;

    vector<double> t,x;
    cout << "t("<< t_max << " "<< Nt << " " << dt <<"):";

    // Разбиение времени
    for (int i=0; i<Nt; i++){
        t.push_back((double)i*dt);
    //    cout << i*dt << " ";
    }
    cout << endl;

    // Разбиение координаты
    cout << "x("<< 2 << " "<< Nx << " " << dx <<"):";
    for (double i=0; i<Nx; i++) {
        double xi = x_l+i*dx;
        x.push_back(xi);
    //    cout << xi << " ";
    }
    cout << endl << endl;
    
    // Основные величины
    vector<vector<double>> p,u,rho,T,e,E,h,S;
    double p4 = 1, p1 = 0.125, u4 = 0, u1 = 0, rho4 = 1, rho1 = 0.1;
    double R = 8.31, T1 = p1/rho1, T4 = p4/rho4;
    double gamma = 1.4;    
    double eps = 1e-15;
    
    // Расчеты    
    cout << " -----------Start_Calculation-------------- " << endl;
    int i0 = 0, i1=Nx/2, i2=Nx/2, i3=Nx/2, i4=Nx/2, i5 = Nx-1, i_initial = 200;
    double M = 0, u_2 = 0, p_2 = 0, rho_2 = 0, u_3=0, p_3=0, rho_3=0;
    double T_2 = 0, T_3 = 0;
    double D = rho4/rho1*(p1-p4)/(rho1-rho4); //скорость ударной волны
    double c1 = sqrt(gamma*p1/rho1), Uc = 0, c3 = 0; // скорость звука справа
    double c4 = sqrt(gamma*p4/rho4), a = 0; // скорость звука слева
    p_3 = p4;
    rho_3 = rho4;
    double A4 = 1;
    
    vector<double> emp(Nx,{0});
    for (size_t k = 0; k < Nt; k++)
    {

        p.push_back(emp);
        u.push_back(emp);
        rho.push_back(emp);
        //e.push_back(x);
        //E.push_back(x);
        //h.push_back(x);
        T.push_back(emp);
        S.push_back(emp);
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
        }
        if(k==0){    
            M = 1.6;
            D = M*c1;
            u_3 = u1;
            p_3 = p1;
            rho_3 = rho1;
            //D = rho4/rho1*(p1-p4)/(rho1-rho4);
            Uc = 2*c1*(M-1/M)/(gamma+1);
        }
        else if(k > 0){
            c3 = (gamma+1)*Uc/2-c4;
            cout << "! " << k  << ":" << Uc << " " << D << endl; 
            //<< "/ " << p_3 << ", " << p1 << "/ " << rho_3 << ", " << rho1 << "/" << endl;//<< p[0][199] << "," << p[0][200] << "," << p[0][201]  << endl;
            //if(k < 1) D = (p_3-p1)/(rho1*u_3);//rho_3/rho1*(p_3-p1)/(rho_3-rho1); ////////
            //else D = rho4/rho1*(p1-p4)/(rho1-rho4);//M*c1;
            
            //D = rho_3/rho1*(p1-p_3)/(rho1-rho_3);
            //cout << "! " << D << " " << D*dt << endl;
            //cout << "! " << c1 << " " << c1*dt << endl;
            int j = i0;
            while(x[j] < -c4*t[k]) j++;
            i1 = j;
            j = i2;
            while(x[j] < c3*t[k]) j++;
            i2 = j;
            cout << "i:" << i3;
            j = i3;
            while(x[j] < Uc*t[k]) j++;
            i3 = j;
            cout << " " << i3 << "|" << i4 << " ";
            j = i4;
            while(x[j] < D*t[k]) j++;            
            i4 = j;
            cout << i4 << endl;

            if(i4 >= Nx) i4 = Nx-1;
            if(i3 >= Nx) i3 = Nx-1;
            if(i2 >= Nx) i2 = Nx-1;
            if(i1 < 0) i1 = 0;

            p_3 = p1*(2*gamma*M*M-(gamma-1))/(gamma+1);
            rho_3 = rho1*(gamma+1)*M*M/((gamma-1)*M*M+2);
            u_3 = c1*2*(M-1/M)/(gamma+1);
            T_3 = p_3/rho_3;

            p_2 = p_3;
            u_2 = u_3;
            rho_2 = rho4*pow(p_3/p4,1/gamma);
            T_2 = p_2/rho_2;

            
            //rho_3 = rho1*((gamma-1)*p4+(gamma+1)*p_3)/((gamma-1)*p_3+(gamma+1)*p4);
            //p_3 = p4*pow(1-(gamma-1)*u_3/(2*c1), 2*gamma/(gamma-1));
            //u_3 = sqrt(2*pow((p_3-p1),2)/(rho1*((gamma-1)*p1+(gamma+1)*p_3)));
            //p_3 = p1*pow(1-(gamma-1)*u_3/(2*c1),2*gamma/(gamma-1));
            //T_3 = T1*(2*gamma*M*M-(gamma-1))/(gamma+1)*(2/((gamma+1)*M*M)+(gamma-1)/(gamma+1));
            //u_3 = 2*c1*(M-1/M)/(gamma+1);
            //u_3 = 1 - 2*c1/(gamma-1)*pow(p_3/p4, (gamma-1)/2/gamma);
            //u_3 = sqrt(2/rho1)*(p_3-p1)/sqrt((gamma-1)*p1+(gamma+1)*p_3);

            
            //cout << "index:" << i2 << " " << i3 << " " << i4 << endl;
            //cout << "rho:" << rho_2 << " " << rho_3 << " " << rho4 << endl;
            

            for(int i=i1; i<Nx;i++){
                a = (-(gamma-1)*x[i]/t[k]+2*c4)/(gamma+1);

                if(i<i2){
                    rho[k][i] = 0;
                    p[k][i] = pow(pow(a*a/gamma,gamma)/A4,1/(gamma-1));
                    u[k][i] = (2*x[i]/t[k]+2*c4)/(gamma+1);
                    T[k][i] = 0;
                }
                if(i>=i2 && i < i3){
                    //if(k%50==0) cout << "i2-i3\n";
                    rho[k][i] = rho_2;
                    p[k][i] = p_2;
                    u[k][i] = u_2;
                    T[k][i] = T_2;
                }
                if(i >= i3 && i < i4){
                    //if(k%50==0) cout << "i3-i4\n";
                    rho[k][i] = rho_3;
                    p[k][i] = p_3;
                    u[k][i] = u_3;
                    T[k][i] = T_3;
                }
            }
        }
        for(int i=0; i<Nx; i++){
            S[k][i] = log(p[k][i]/pow(rho[k][i], gamma));//+S[0][i];
        }
        M = D/c1;
        
        //D = rho_3/rho1*(p_3-p1)/(rho_3-rho1);
        //cout << "M = " << M  << " " << p_3/p4 << endl;
    }
    //cout << dx << " " << Nx << " | " << dt << " " << Nt << endl;
    cout << " -----------Finish_Calculation-------------- " << endl;

    //**********************************************************//
    //**********Вывод в файл**********************************//

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
    //cout << solver.sound_speed() << endl;

    
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

    out1.close();
    out1.open("result_T.txt");
    if(out1.is_open()){
        for (int i=0; i<Nx; i++){
            out1 << x[i] << " ";
        }
        out1 << endl;
        for (int i=0; i<Nt; i++){
            out1 << t[i] << " ";
            for(int j=0; j<Nx; j++) out1  << T[i][j] << " ";
            out1 << endl;
        }
    }
    out1.close();

   // cout << " -----------Begin_Plotting------------- " << endl;
   // system("plotting.exe");
   // cout << " -----------End_Plotting------------- " << endl;
    return 0;
}

