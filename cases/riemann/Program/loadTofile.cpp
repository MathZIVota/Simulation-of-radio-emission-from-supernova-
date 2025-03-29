#include "Tube.hpp"
void Tube::UploadToFile(){
        ofstream out;
        out.open("result_rho.txt");
        if(out.is_open()){
            for (int i=0; i<Nx; i++){
                out <<  x[i] << " ";
            }
            out << endl;
            for (int i=0; i<Nt; i++){
                out <<  t[i] << " ";
                for(int j=0; j<Nx; j++) out  <<  rho[i][j] << " ";
                out << endl;
            }
        }
        out.close();
        
        out.open("result_p.txt");
        if(out.is_open()){
            for (int i=0; i<Nx; i++){
                out <<  x[i] << " ";
            }
            out << endl;
            for (int i=0; i<Nt; i++){
                out <<  t[i] << " ";
                for(int j=0; j<Nx; j++) out  <<  p[i][j] << " ";
                out << endl;
            }
        }
        out.close();

        out.open("result_u.txt");
        if(out.is_open()){
            for (int i=0; i<Nx; i++){
                out <<  x[i] << " ";
            }
            out << endl;
            for (int i=0; i<Nt; i++){
                out <<  t[i] << " ";
                for(int j=0; j<Nx; j++) out  <<  u[i][j] << " ";
                out << endl;
            }
        }
        out.close();

        out.open("result_S.txt");
        if(out.is_open()){
            for (int i=0; i<Nx; i++){
                out <<  x[i] << " ";
            }
            out << endl;
            for (int i=0; i<Nt; i++){
                out <<  t[i] << " ";
                for(int j=0; j<Nx; j++) out  <<  S[i][j] << " ";
                out << endl;
            }
        }
        out.close();

        out.open("result_T.txt");
        if(out.is_open()){
            for (int i=0; i<Nx; i++){
                out <<  x[i] << " ";
            }
            out << endl;
            for (int i=0; i<Nt; i++){
                out <<  t[i] << " ";
                for(int j=0; j<Nx; j++) out  <<  T[i][j] << " ";
                out << endl;
            }
        }
        out.close();
    }