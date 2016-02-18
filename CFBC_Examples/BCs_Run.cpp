//
//  BC_Run.cpp
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
    const int Nx=100;
    const int Ny=81;
    const int Nz=100;
    
    // Size of domain
    const Real Lx= 2*pi;
    const Real a=-1.0;
    const Real b= 1.0;
    const Real Lz= 2*pi;
    
    // Define flow parameters
    const Real Reynolds = 2000.0;
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
    const Real decay = 0.5;
    const Real perturbMag = 0.10;
    
    FlowField u(Nx,Ny,Nz,3,Lx,Lz,a,b); // Velocity flowfield (u=object)
    FlowField q(Nx,Ny,Nz,1,Lx,Lz,a,b); // Pressure flowfield (q=object)
    u.addPerturbations(kxmax,kzmax,1,decay); // Add pertubations to base flow
    u *= perturbMag/L2Norm(u);
    
    
    
    
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
    
    
    
    //*******---BCs ---*******************************
    // ***********************************************
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
    
    
    
    if (ReStart){
        
        // Load BCs
        cout << "loading BCs:";
        fstream BCfile("BCs.asc", ios_base::in|ios_base::out);
        
        for (int UL=0;UL<2;++UL){
            for (int imx=0;imx<u.Mx();++imx){
                for (int imz=0;imz<u.Mz();++imz){
                    for (int uvw=0;uvw<3;++uvw){
                        BCfile >> BC.cmplx(imx,UL,imz,uvw);
                    }
                }
            }
            
        }
        BCfile.close();
        cout << "Done" << endl;
        
    }
    
    
    //*******************************************************
    
    
    
    
    
    int it=0; // time-step counter
    
    
    
    if (Controlled){
        dns.reset_dtIH(dt,BC);
    }
    
    
    
    
    
    
    // Time-stepping loop - Controller and sim
    for (double t=T0; t<T1+(dT/2); t += dT) {
        
        
        //*****************************************************
        // Assign BCs to v-component on kx=0,kz=1 at both walls
        //*****************************************************
        BC.cmplx(0,0,1,1) = 0.01*sin(t);
        BC.cmplx(0,1,1,1) = 0.01*sin(t);
        
        
        
        
        
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
        
        cout << "   BC upper == " << BC.cmplx(0,1,1,1) << endl;
        cout << "   BC lower == " << BC.cmplx(0,0,1,1) << endl;
        
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
        
        
        //-------------------------------------
        //-------------------------------------
        
        
        
        
        
        if (it == nSave) {
            cout << "Saving flowfields..." << endl;
            u.save("u"+i2s(iround(t)));
            q.save("q"+i2s(iround(t)));
            cout << "done" << endl;
            
            cout << "Saving BCs:";
            
            // Save BCs
            ofstream BCfile;
            BCfile.open("BCs.asc");
            for (int UL=0;UL<2;++UL){
                for (int imx=0;imx<u.Mx();++imx){
                    
                    for (int imz=0;imz<u.Mz();++imz){
                        for (int uvw=0;uvw<3;++uvw){
                            BCfile << BC.cmplx(imx,UL,imz,uvw) << endl;
                        }
                    }
                }
                
            }
            BCfile.close();
            cout << " Done" << endl;
            
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
            dns.advance_inhom(u,q,BC,dt.n()); // Controlled
        }
        else {
            dns.advance(u,q,dt.n()); // Not controlled
        }
        
        
        
        
        
        
        it +=1; // timestep index
        cout << endl;
        
        
        
    }
    cout << "Done!" << endl;
    
    
    
    
}      




