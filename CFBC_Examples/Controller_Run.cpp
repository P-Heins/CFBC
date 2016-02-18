//
//  Controller_Run.cpp
//
//
//  Created by Peter Heins on 27/05/2015.
//  Copyright 2015 University of Sheffield. All rights reserved.


#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "channelflow/flowfield.h"
#include "channelflow/periodicfunc.h"
#include "channelflow/dns.h"
#include "channelflow/utilfuncs.h"
#include "channelflow/controller.h"
#include "channelflow/turbstats.h"

using namespace std;
using namespace channelflow;

int main(){
    
    
    
    
    bool ReStart    = false;
    bool Controlled = true;
    
    
    
    // Num of data points in x, y and z directions
    const int Nx=80;
    const int Ny=151;
    const int Nz=80;
    
    // Size of domain
    const Real Lx= 2*pi;
    const Real a=-1.0;
    const Real b= 1.0;
    const Real Lz= 2*pi;
    
    // Define flow parameters
    const Real Reynolds = 2230.0;
    const Real nu = 1/Reynolds; // Kinematic viscosity
    
    // Variable time-stepping parameters
    const Real CFLmin= 0.30;
    const Real CFLmax= 0.70;
    const Real dtmax=  0.010; // <<<< CONTROLLER time interval
    const Real dtmin=  0.00010;
    const Real dt0=    dtmax;
    const bool variable_dt = true;
    const Real T0 = 0.0; // Start time
    const Real T1 = 10.0; // End time
    const Real dT = 0.050; // dtmax; //
    
    // Save parameters
    Real SaveInt = T1-T0;
    int nSave = iround(SaveInt/dT);
    
    DNSFlags flags; // DNSFlags = class, flags = object
    flags.baseflow      = Parabolic;
    flags.nonlinearity  = SkewSymmetric;
    flags.timestepping  = SBDF3;
    flags.initstepping  = CNRK2;
    flags.dealiasing    = NoDealiasing;
    flags.constraint    = BulkVelocity;
    flags.Ubulk         = 2.0/3.0;
    
    
    cout << setprecision(8); // Sets number of output decimal places
    
    Vector x = periodicpoints(Nx, Lx); // x = object, periodicpoints = fn
    Vector y = chebypoints(Ny,a,b); // y = object, chebypoints = fn
    Vector z = periodicpoints(Nz, Lz); // z = object, pp = fn
    
    
    // FlowField = class
    const int kxmax = 4;
    const int kzmax = 4;
    const Real decay = 0.7;
    const Real perturbMag = 0.50;
    
    FlowField u(Nx,Ny,Nz,3,Lx,Lz,a,b); // Velocity flowfield (u=object)
    FlowField q(Nx,Ny,Nz,1,Lx,Lz,a,b); // Pressure flowfield (q=object)
    u.addPerturbations(kxmax,kzmax,1,decay); // Add pertubations to base flow
    u *= perturbMag/L2Norm(u);
    
    
    
    FlowField F(u);
    
    
    
    cout << "Optimising FFTW..." << flush;
    fftw_loadwisdom();
    u.optimizeFFTW();
    fftw_savewisdom();
    cout << "Done" << endl;
    
    // Construct DNS
    cout << "Constructing DNS..." << flush;
    // TimeStep = class, dt = object
    TimeStep dt(dt0, dtmin, dtmax, dT, CFLmin, CFLmax, variable_dt);
    // DNS = class, dns = object
    DNS dns(u, nu, dt, flags, T0, Controlled);
    cout << "Done" << endl;
    
    
    FlowField tmp(Nx,Ny,Nz,6,Lx,Lz,a,b);
    ChebyTransform trans(Ny);
    
    ChebyCoeff Ubase_Stat(Ny,a,b,Physical);
    for (int i=0; i<Ny; ++i)
        Ubase_Stat[i] = 1 - square(y[i]);
    //Ubase.save("Ubase");
    Ubase_Stat.makeSpectral(trans);
    
    TurbStats stats(Ubase_Stat, nu);
    
    
    
    // Make data directory
    mkdir("data-walllaw");
    
    
    
    //*******---CONTROLLER---*******************************
    // ******************************************************
    FlowField BC(Nx,2,Nz,3,Lx,Lz,a,b,Physical,Physical);
    
    
    for (int ny=0;ny<2;++ny){
        for (int nx=0;nx<Nx;++nx){
            for (int nz=0;nz<Nz;++nz){
                for (int i=0;i<3;++i){
                    BC(nx,ny,nz,i) = 0.0;
                }
            }
        }
    }
    BC.makeState(Spectral,Physical);
    
    
    // Area of wavenumber space controlled
    int minKx = 0;
    int maxKx = 0;
    int minKz = 0;
    int maxKz = 10;
    
    int maxMx = (maxKx-minKx)+1;
    int maxMz = (maxKz-minKz)+1;
    
    int uvw = 1; // u=0, v=1, w=2
    int NumInputs = 8; // inputs to controller i.e. no. of flow measurements - 2*(lower/upper)*(real/imag)
    const char* SIFile = "Mult_Control_Mat/StateInfo.bin"; //state info file
    
    double tau = 0.010; // actuator time constant
    
    bool Spectral_states = false;
    
    Controller Hinf(u.Ny(),minKx,maxKx,minKz,maxKz,uvw,SIFile,NumInputs,tau,Spectral_states);
    
    double*** CStateMat = Hinf.ConStates(); // controller state array
    
    
    // 3D input/output array
    int NIO = maxMx*maxMz; // no. of kx,kz pairs
    double *** IO; // input/output data array
    IO = new double**[2];
    
    for (int i=0;i<2;++i){
        IO[i] = new double*[3+NumInputs];
        for (int j=0;j<3+NumInputs;++j){
            IO[i][j] = new double[NIO];
        }
    }
    
    for (int i=0;i<2;++i){
        for (int j=0;j<3+NumInputs;++j){
            for (int k=0;k<NIO;++k){
                IO[i][j][k] = 0.0;
            }
        }
    }
    
    if (ReStart){
        // Load controller states
        cout << "Loading Controller States: ";
        double** CStateInf = Hinf.CStateInfo();
        fstream ConStatesfile("Controller_States.asc", ios_base::in|ios_base::out);
        ConStatesfile.seekg(0);
        for (int i=0;i<(maxKx-minKx+1);++i){
            for (int j=0;j<(maxKz-minKz+1);++j){
                if (i==0 && j==0){
                    continue;
                }
                else {
                    //cout << CStateInf[i][j] << endl;
                    for (int k=0;k<CStateInf[i][j];++k){
                        ConStatesfile >> CStateMat[i][j][k];
                        
                    }
                }
                
            }
        }
        ConStatesfile.close();
        cout << "Done" << endl;
        
        // Load BCs
        cout << "loading BCs:";
        fstream BCfile("BCs.asc", ios_base::in|ios_base::out);
        ConStatesfile.seekg(0);
        for (int UL=0;UL<2;++UL){
            for (int imxc=0;imxc<maxMx;++imxc){
                int ikx = Hinf.kx_c(imxc);
                int imx = u.mx(ikx);
                for (int imz=0;imz<maxMz;++imz){
                    BCfile >> BC.cmplx(imx,UL,imz,uvw);
                }
            }
            
        }
        BCfile.close();
        cout << "Done" << endl;
        
    }
    
    
    //*******************************************************
    
    
    
    // Data output to files
    mkdir("Input_Output");
    
    
    ofstream Outfile;
    Outfile.open("Input_Output/Shear_Stresses.asc");
    ofstream Infile;
    Infile.open("Input_Output/Control_Signals.asc");
    
    ofstream Dragfile;
    Dragfile.open("Input_Output/Drag.asc");
    ofstream Ustarfile;
    Ustarfile.open("Input_Output/Ustar.asc");
    
    ofstream Energyfile;
    Energyfile.open("Input_Output/Energy.asc");
    
    ofstream EnergyMxMz_File;
    EnergyMxMz_File.open("Input_Output/EnergyMxMz.asc");
    
    
    
    Outfile << "Sim Output / Controller Input \nTime  kx   kz  Str_R_Up  Str_I_Up   Str_R_Low   Str_I_Low Span_R_Up  Span_I_Up   Span_R_Low   Span_I_Low  \n";
    Infile << "Sim Input / Controller Output  \nTime  kx   kz  R_Up  I_Up   R_Low   I_Low\n";
    
    
    Dragfile << "Time  Drag_Lower  Drag_Upper\n";
    Ustarfile << "Time   Ustar\n";
    
    Energyfile << "Time  Perturbation Energy\n";
    
    EnergyMxMz_File << "Time   kx   kz   Energy\n";
    
    
    int KxC_1 = 0; // For printing BC at this pair to screen
    int KzC_1 = 1;
    int mxC_1 = BC.mx(KxC_1);
    int mzC_1 = BC.mz(KzC_1);
    
    int KxC_2 = 0; // For printing BC at this pair to screen
    int KzC_2 = 2;
    int mxC_2 = BC.mx(KxC_2);
    int mzC_2 = BC.mz(KzC_2);
    
    
    int it=0; // time-step counter
    
    
    
    if (Controlled){
        dns.reset_dtIH(dt,BC);
    }
    
    
    
    
    
    
    // Time-stepping loop - Controller and sim
    for (double t=T0; t<T1+(dT/2); t += dT) {
        
        Real RUpBC_1 = real(BC.cmplx(mxC_1,1,mzC_1,1)); // to print to screen
        Real IUpBC_1 = imag(BC.cmplx(mxC_1,1,mzC_1,1));
        Real RLowBC_1 = real(BC.cmplx(mxC_1,0,mzC_1,1));
        Real ILowBC_1 = imag(BC.cmplx(mxC_1,0,mzC_1,1));
        
        Real RUpBC_2 = real(BC.cmplx(mxC_2,1,mzC_2,1)); // to print to screen
        Real IUpBC_2 = imag(BC.cmplx(mxC_2,1,mzC_2,1));
        Real RLowBC_2 = real(BC.cmplx(mxC_2,0,mzC_2,1));
        Real ILowBC_2 = imag(BC.cmplx(mxC_2,0,mzC_2,1));
        
        ChebyCoeff u00 = Re(u.profile(0,0,0));
        ChebyCoeff du00dy = diff(u00);
        Real drag_L = nu*du00dy.eval_a();
        Real drag_U = nu*du00dy.eval_b();
        
        Real Energy = L2Norm2(u);
        
        cout << "          t == " << t << endl;
        cout << "         dt == " << dt << endl;
        cout << "        CFL == " << dns.CFL() << endl;
        cout << " L2Norm2(u) == " << Energy << endl;
        cout << "divNorm2(u) == " << divNorm(u)/L2Norm(u) << endl;
        cout << "      Ubulk == " << dns.Ubulk() << endl;
        cout << "      ubulk == " << Re(u.profile(0,0,0)).mean()/2 << endl;
        cout << "       dPdx == " << dns.dPdx() << endl;
        cout << "       drag == " << 0.5*(abs(drag_L)+abs(drag_U)) << endl;
        cout << "Kx="+i2s(KxC_1)+", Kz="+i2s(KzC_1)+" Upper BC ==" << RUpBC_1 << "+(" << IUpBC_1 << ")i" << endl;
        cout << "Kx="+i2s(KxC_1)+", Kz="+i2s(KzC_1)+" Lower BC ==" << RLowBC_1 << "+(" << ILowBC_1 << ")i" << endl;
        cout << "Kx="+i2s(KxC_2)+", Kz="+i2s(KzC_2)+" Upper BC ==" << RUpBC_2 << "+(" << IUpBC_2 << ")i" << endl;
        cout << "Kx="+i2s(KxC_2)+", Kz="+i2s(KzC_2)+" Lower BC ==" << RLowBC_2 << "+(" << ILowBC_2 << ")i" << endl;
        
        //--------------------------------------------
        //--------------------------------------------
        // Wall-law and turb stats
        u00.makePhysical();
        u00.save("data-walllaw/u00");
        stats.addData(u,tmp);
        cout << "centerline Re = " << stats.centerlineReynolds() << endl;
        cout << " parabolic Re = " << stats.parabolicReynolds() << endl;
        cout << "      bulk Re = " << stats.bulkReynolds() << endl;
        cout << "        ustar = " << stats.ustar() << endl;
        
        stats.msave("data-walllaw/uu");
        stats.msave("data-walllaw/uustar", true);
        
        Real ustar = stats.ustar();
        Vector yp = stats.yplus();
        yp.save("data-walllaw/yp");
        
        ChebyCoeff Umean = stats.U();
        Umean.makeSpectral(trans);
        ChebyCoeff Umeany = diff(Umean);
        Umean.makePhysical(trans);
        Umeany.makePhysical(trans);
        Umean.save("data-walllaw/Umean");
        Umeany.save("data-walllaw/Umeany");
        
        Umean /= ustar;
        Umean.save("data-walllaw/Uplus");
        
        ChebyCoeff ubase = stats.ubase();
        ubase.save("data-walllaw/ubase");
        
        ChebyCoeff uv = stats.uv();
        uv.save("data-walllaw/uv");
        save(ustar, "data-walllaw/ustar");
        save(nu, "data-walllaw/nu");
        
        Ustarfile << t << "  " << ustar << endl;
        //-------------------------------------
        //-------------------------------------
        
        
        
        
        
        if (it == nSave) {
            cout << "Saving flowfields..." << endl;
            u.save("u"+i2s(iround(t)));
            q.save("q"+i2s(iround(t)));
            cout << "done" << endl;
            
            cout << "Saving Controller States and BCs:";
            
            // Save Controller States to file
            double** CStateInf = Hinf.CStateInfo();
            ofstream ConStatesfile;
            ConStatesfile.open("Controller_States.asc");
            for (int i=0;i<(maxKx-minKx+1);++i){
                for (int j=0;j<(maxKz-minKz+1);++j){
                    for (int k=0;k<CStateInf[i][j];++k){
                        if (i==0 && j==0){
                            continue;
                        }
                        else {
                            ConStatesfile << CStateMat[i][j][k] << endl;
                        }
                    }
                }
            }
            ConStatesfile.close();
            
            // Save BCs
            ofstream BCfile;
            BCfile.open("BCs.asc");
            for (int UL=0;UL<2;++UL){
                for (int imxc=0;imxc<maxMx;++imxc){
                    int ikx = Hinf.kx_c(imxc);
                    int imx = u.mx(ikx);
                    for (int imz=0;imz<maxMz;++imz){
                        BCfile << BC.cmplx(imx,UL,imz,uvw) << endl;
                    }
                }
                
            }
            BCfile.close();
            cout << " Done" << endl;
            
        }
        
        
               
        
        Dragfile << t << " " << drag_L << " " << drag_U << endl;
        Energyfile << t << " " << Energy << endl;
        //***************************************************
             
              
        
        
        // Record energy at each mx,mz
        for (int imx=0;imx<1;++imx){
            for (int imz=0;imz<u.Mz();++imz){
                EnergyMxMz_File << t << " " << u.kx(imx) << " " << u.kz(imz) << " " << u.energy(imx,imz,true) << endl;
            }
        }
        
        
        
        if (dt.adjust(dns.CFL())) {
            cerr << "resetting dt to " << dt << ", new CFL ==" << dt.CFL() << endl;
            if (Controlled){
                dns.reset_dtIH(dt,BC); // Controlled
            }
            else {
                dns.reset_dt(dt); // Not controlled
            }
            
        }
        
        if (dns.CFL()>1.0){
            cout << "Error: CFL>1" << endl;
            return(0);
        }
        
        
        if (Controlled){
            dns.advance_inhom_CON(Hinf,u,q,BC,F,dt.n(),IO,CStateMat); // Controlled
        }
        else {
            dns.advance(u,q,dt.n()); // Not controlled
        }
        
        
        // Print IO data to file
        for (int k=0;k<NIO;++k){
            for (int j=0;j<3+NumInputs;++j){
                if (j<3+NumInputs-1){
                    Outfile << IO[0][j][k] << " ";
                }
                else {
                    Outfile << IO[0][j][k] << endl;
                }
            }                  
        }
        
        for (int k=0;k<NIO;++k){
            for (int j=0;j<3+4;++j){
                if (j<7-1){
                    Infile << IO[1][j][k] << " ";
                }
                else {
                    Infile << IO[1][j][k] << endl;
                }               
            }                  
        }
        
        
        
        
        it +=1; // timestep index
        cout << endl;
        
        
        
    }
    cout << "Done!" << endl;
    
    
    Outfile.close();
    Infile.close();
    Dragfile.close();
    Ustarfile.close();
    Energyfile.close();
    EnergyMxMz_File.close();
    
    
}      




