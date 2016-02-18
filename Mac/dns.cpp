/* dns.cpp: time-integration classes for spectral Navier-Stokes simulation
 * Channelflow-1.0
 *
 * Copyright (C) 2001-2007  John F. Gibson
 *
 * Center for Nonlinear Science
 * School of Physics
 * Georgia Institute of Technology
 * Atlanta, GA 30332-0430
 * 404 385 2509
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, U
 */
//#include <Accelerate/Accelerate.h>
#include "channelflow/dns.h"
#include "channelflow/utilfuncs.h"
//#include "channelflow/controller.h"
#include <limits.h>
#include <complex.h>

//#include "orrsommfunc.h"
//#include <fstream> // tmp debugging need

using namespace std;
using namespace channelflow;

namespace channelflow {
    
    const Real EPSILON = 1e-11;
    
    //NonlinearMethod convection = Convection;
    
    // This used to be in diffops.cpp, but its logic is closely tied to the
    // configuration of the DNS algorithms that use it, so I moved it here.
    
    void navierstokesNL(const FlowField& u, const FlowField& ubase,
                        const ChebyCoeff& Ubase, FlowField& f, FlowField& tmp,
                        FlowField& tmp2, NonlinearMethod& method) {
        
        if (method == LinearAboutField)
            linearAboutFieldNL(u, ubase, Ubase, f, tmp, tmp2);
        else
            navierstokesNL(u, Ubase, f, tmp, method);
    }
    
    void navierstokesNL(const FlowField& u_, const ChebyCoeff& Ubase,
                        FlowField& f, FlowField& tmp, NonlinearMethod& method) {
        
        assert(u_.xzstate() == Spectral && u_.ystate() == Spectral);
        assert(Ubase.state() == Spectral);
        
        if (method == LinearAboutProfile)
            linearizedNL(u_, Ubase,f);
        else {
            FlowField& u = (FlowField&) u_;
            u += Ubase;
            switch (method) {
                case Rotational:
                    rotationalNL(u,f,tmp);
                    break;
                case Convection:
                    convectionNL(u,f,tmp);
                    break;
                case LinearAboutField:
                    cerr << "error in navierstokesNL : LinearAboutField requires a different\n"
                    << "function signature.\nPlease submit a bug report.";
                    exit(1);
                    break;
                case SkewSymmetric:
                    skewsymmetricNL(u,f,tmp);
                    break;
                case Divergence:
                    divergenceNL(u,f,tmp);
                    break;
                case Alternating:
                    divergenceNL(u,f,tmp);
                    method = Alternating_;
                    break;
                case Alternating_:
                    convectionNL(u,f,tmp);
                    method = Alternating;
                    break;
                default:
                    cferror("navierstokesNL(method, u,U,f,tmp) : unknown method");
            }
            u -= Ubase;
        }
    }
    
    ostream& operator<<(ostream& os, const TimeStep& dt) {
        os << "{dt=" << dt.dt()
        << ", n=" << dt.n()
        << ", dT=" << dt.dT()
        << ", N=" << dt.N()
        << ", dtmin=" << dt.dtmin()
        << ", dtmax=" << dt.dtmax()
        << ", CFLmin=" << dt.CFLmin()
        << ", CFL=" << dt.CFL()
        << ", CFLmax=" << dt.CFLmax()
        << ", variable=" << dt.variable()
        << "}";
        return os;
    }
    
    TimeStep::TimeStep()
    :
    n_(0),
    N_(0),
    dt_(0),
    dtmin_(0),
    dtmax_(0),
    dT_(0.0),
    T_(0.0),
    CFLmin_(0.0),
    CFL_(0.0), // will take on meaninful value after first adjust
    CFLmax_(0.0),
    variable_(false)
    {}
    
    
    TimeStep::TimeStep(Real dt, Real dtmin, Real dtmax, Real dT,
                       Real CFLmin, Real CFLmax, bool variable)
    :
    n_(0),
    N_(0),
    dt_(dt),
    dtmin_(dtmin),
    dtmax_(dtmax),
    dT_(dT),
    T_(0.0),
    CFLmin_(CFLmin),
    CFL_((CFLmax+CFLmin)/2), // will take on meaningful value after first adjust
    CFLmax_(CFLmax),
    variable_(variable)
    {
        
        if (dtmin < 0 || dt < dtmin || dtmax < dt) {
            cerr << "error in TimeStep::TimeStep(dt, dtmin, dtmax, dT, CFLmin, CFLmax, variable) :\n"
            << "condition 0 <= dtmin <= dt <= dtmax does not hold" << endl;
            exit(1);
        }
        if (CFLmin < 0 || CFLmax < CFLmin) {
            cerr << "error in TimeStep::TimeStep(dt, dtmin, dtmax, dT, CFLmin, CFLmax, variable) :\n"
            << "condition 0 <= CFLmin <= CFLmax does not hold" << endl;
            exit(1);
        }
        if (dT < dtmin) {
            cerr << "error in TimeStep::TimeStep(dt, dtmin, dtmax, dT, CFLmin, CFLmax, variable) :\n"
            << "dT < dtmin" << endl;
            exit(1);
        }
        
        // Adjust dt to be integer divisor of dT. At this point we have 0 <= dtmin <= dt and dtmin <= dT
        n_ = Greater(iround(dT/dt), 1);
        dt_ = dT_/n_;                    // 0 <= dt <= dT and dtmin <= dT
        
        // Bump up or down to get within dtmin, dtmax
        while (dt_ < dtmin_ && n_ >= 2 && dT_ != 0) {
            dt_ = dT_/--n_;                 // guaranteed to terminate at  dtmin <= dt == dT
        }
        
        
        while (dt_ > dtmax_ && n_ <= INT_MAX && dT_ != 0) {
            dt_ = dT_/++n_;                 // guaranteed to terminate at  dt == dT/INT_MAX
        }
        assert(dt_>0 && dt_<=dT);
        assert(dt_>=dtmin && dt<=dtmax);
    }
    
    
    
    
    
    
    // relations
    // n*dt = dT
    bool TimeStep::adjust(Real CFL, bool verbose, ostream& os) {
        CFL_ = CFL;
        if (!variable_ && CFLmin_ <= CFL && CFL <= CFLmax_)
            return false;
        else
            return adjustToMiddle(CFL, verbose, os);
    }
    
    
    bool TimeStep::adjustToMiddle(Real CFL, bool verbose, ostream& os) {
        
        if (dtmin_ == dtmax_ || dT_ == 0.0)
            return false;
        
        // New update algorithm puts CFL at midpoint btwn bounds
        // Aim for      CFL' == (CFLmax+CFLmin)/2
        // Change is    CFL' == CFL * dt'/dt
        // (CFLmax+CFLmin)/2 == CFL * n/n'      since dt=dT/n
        // So             n' == 2 n CFL/(CFLmax+CFLmin)
        //
        int n  = Greater(iround(2*n_*CFL/(CFLmax_ + CFLmin_)), 1);
        
        /*
        //if (n % 2 != 0) {
        //    ++n;            
        //}
        
        if (n==3){
            n = 4;
        }
        if (n>5 && n<8){
            n = 8;
        }
        if (n==9){
            n = 10;
        }
        if (n>10 && n<16){
            n = 16;
        }
        if (n>16){
            n = 20;
        }
        */
        
        Real dt = dT_/n;
               
        // Bump dt up or down to get within [dtmin, dtmax]
        while (dt < dtmin_ && dt < dT_)
            dt = dT_/--n;                 // guaranteed to terminate at  dtmin <= dt == dT
        while (dt > dtmax_ && n <= INT_MAX)
            dt = dT_/++n;                 // guaranteed to terminate at  dtmin <= dt == dT
        
        CFL *= dt/dt_;
        
        // Check to see if adjustment took dt out of range
        if (verbose && (CFL > CFLmax_ || CFL < CFLmin_)) {
            os << "TimeStep::adjust(CFL) : dt " << (CFL > CFLmax_ ? "bottomed" : "topped") << " out at\n"
            << " dt  == " << dt << endl
            << " CFL == " << CFL  << endl
            << " n   == " << n <<endl;
        }
        
        // If final choice for n differs from original n_, reset internal values
        bool adjustment = (n == n_) ? false : true;
        if (adjustment && verbose) {
            os << "TimeStep::adjust(CFL) { " << endl;
            os << "   n : " << n_ << " -> " << n << endl;
            os << "  dt : " << dt_ << " -> " << dt << endl;
            os << " CFL : " << CFL_ << " -> " << CFL << endl;
            os << "}" << endl;
            n_ = n;
            dt_ = dt;
            CFL_ = CFL;
        }
        return adjustment;
    }
    
    
    bool TimeStep::adjust_for_T(Real T, bool verbose, ostream& os) {
        T_ = T;
        if (T < 0) {
            cerr << "TimeStep::adjust_for_T : can't integrate backwards in time.\n"
            << "Exiting." << endl;
            exit(1);
        }
        if (T == 0) {
            bool adjustment = (dt_ == 0) ? false : true;
            dt_ = 0;
            n_  = 0;
            dT_ = 0;
            T_ = 0;
            return adjustment;
        }
        int N = Greater(iround(T/dT_), 1);
        Real dT = T/N;
        int n = Greater(iround(dT/dt_), 1);
        Real dt = dT/n;
        
        while (dt < dtmin_ && n>2 && dT != 0)
            dt = dT/--n;
        while (dt > dtmax_ &&  n <= INT_MAX && dT != 0)
            dt = dT/++n;
        
        Real CFL = dt*CFL_/dt_;
        
        bool adjustment = (dt == dt_) ? false : true;
        if (adjustment  && verbose) {
            os << "TimeStep::adjust_for_T(Real T) { " << endl;
            os << "   T : " << T << endl;
            os << "  dT : " << dT_ << " -> " << dT << endl;
            os << "  dt : " << dt_ << " -> " << dt << endl;
            os << "  n  : " << n_  << " -> " << n << endl;
            os << "  N  : " << N_  << " -> " << N << endl;
            os << " CFL : " << CFL_ << " -> " << (dt*CFL_)/dt_ << endl;
            os << "}" << endl;
        }
        n_ = n;
        N_ = N;
        dt_ = dt;
        dT_ = dT;
        CFL_ = CFL;
        
        return adjustment;
    }
    
    int  TimeStep::n() const {return n_;}
    int  TimeStep::N() const {return N_;}
    Real TimeStep::dt() const {return dt_;}
    Real TimeStep::dT() const {return dT_;}
    Real TimeStep::T() const {return T_;}
    Real TimeStep::dtmin() const {return dtmin_;}
    Real TimeStep::dtmax() const {return dtmax_;}
    Real TimeStep::CFL() const {return CFL_;}
    Real TimeStep::CFLmin() const {return CFLmin_;}
    Real TimeStep::CFLmax() const {return CFLmax_;}
    bool TimeStep::variable() const {return variable_;}
    TimeStep::operator Real() const {return dT_/n_;}
    
    //====================================================================
    DNS::DNS()
    :
    main_algorithm_(0),
    init_algorithm_(0)
    {}
    
    DNS::DNS(const FlowField& u, const ChebyCoeff& Ubase, Real nu, Real dt,
             const DNSFlags& flags, Real t, bool controlled)
    :
    main_algorithm_(0),
    init_algorithm_(0)
    {
        main_algorithm_ = newAlgorithm(u, Ubase, nu, dt, flags, t, controlled);
        
        if (!main_algorithm_->full() &&
            flags.initstepping != flags.timestepping) {
            DNSFlags initflags = flags;
            initflags.timestepping = flags.initstepping;
            init_algorithm_ = newAlgorithm(u,Ubase,nu,dt,initflags,t, controlled);
            
            // Safety check
            if (init_algorithm_->Ninitsteps() != 0)
                cerr << "DNS::DNS(u, Ubase, nu, dt, flags, t) :\n" << flags.initstepping
                << " can't initialize " << flags.timestepping
                << " since it needs initialization itself.\n";
        }
    }
    /*
    DNS::DNS(const FlowField& u, Real nu, Real dt, const DNSFlags& flags, Real t)
    :
    main_algorithm_(0),
    init_algorithm_(0)
    {
        ChebyCoeff Ubase(u.Ny(), u.a(), u.b(), Spectral);
        switch(flags.baseflow) {
            case PlaneCouette :
                Ubase[1] = 1;     // Ubase = y;
                break;
            case Parabolic :
                Ubase[0] =  0.5;  // Ubase = 1 - y^2;
                Ubase[2] = -0.5;
                break;
            case Zero:
            default:
                break;
        }
        
        main_algorithm_ = newAlgorithm(u, Ubase, nu, dt, flags, t);
        
        if (!main_algorithm_->full() &&
            flags.initstepping != flags.timestepping) {
            DNSFlags initflags = flags;
            initflags.timestepping = flags.initstepping;
            init_algorithm_ = newAlgorithm(u,Ubase,nu,dt,initflags,t);
            
            // Safety check
            if (init_algorithm_->Ninitsteps() != 0)
                cerr << "DNS::DNS(u, nu, dt, flags, t) :\n" << flags.initstepping
                << " can't initialize " << flags.timestepping
                << " since it needs initialization itself.\n";
        }
    }
    */
    
    DNS::DNS(const DNS& dns)
    :
    main_algorithm_(dns.main_algorithm_ ? dns.main_algorithm_->clone() : 0),
    init_algorithm_(dns.init_algorithm_ ? dns.init_algorithm_->clone() : 0)
    {}
    
    
    
    ////////////////////////////////////////////////////
    // Added by Peter H 10/13
    DNS::DNS(const FlowField& u, Real nu, Real dt, const DNSFlags& flags, Real t, bool controlled)
    :
    main_algorithm_(0),
    init_algorithm_(0)
    {
        
        
        ChebyCoeff Ubase(u.Ny(), u.a(), u.b(), Spectral);
        switch(flags.baseflow) {
            case PlaneCouette :
                Ubase[1] = 1;     // Ubase = y;
                break;
            case Parabolic :
                Ubase[0] =  0.5;  // Ubase = 1 - y^2;
                Ubase[2] = -0.5;
                break;
            case Zero:
            default:
                break;
        }
        
        main_algorithm_ = newAlgorithm(u, Ubase, nu, dt, flags, t, controlled);
        
        if (!main_algorithm_->full() &&
            flags.initstepping != flags.timestepping) {
            DNSFlags initflags = flags;
            initflags.timestepping = flags.initstepping;
            init_algorithm_ = newAlgorithm(u, Ubase, nu, dt, initflags, t, controlled);
            
            // Safety check
            if (init_algorithm_->Ninitsteps() != 0)
                cerr << "DNS::DNS(u, nu, dt, flags, t) :\n" << flags.initstepping
                << " can't initialize " << flags.timestepping
                << " since it needs initialization itself.\n";
        }
    }

    
    ////////////////////////////////////////////////////
    
    
    
    /*
    DNSAlgorithm* DNS::newAlgorithm(const FlowField& u, const ChebyCoeff& Ubase,
                                    Real nu,Real dt, const DNSFlags& flags,Real t){
        DNSAlgorithm* alg = 0;
        switch (flags.timestepping) {
            case CNFE1:
            case SBDF1:
            case SBDF2:
            case SBDF3:
            case SBDF4:
                alg = new MultistepDNS(u, Ubase, nu, dt, flags, t);
                break;
            case CNRK2:
                alg = new RungeKuttaDNS(u, Ubase, nu, dt, flags, t);
                break;
            case SMRK2:
            case CNAB2:
                alg = new CNABstyleDNS(u, Ubase, nu, dt, flags, t);
                break;
            default:
                cerr << "DNS::newAlgorithm : algorithm " << flags.timestepping
                << " is unimplemented" << endl;
        }
        return alg;
    }
    */
    
    
    ////////////////////////////////////////////////////
    // Added by Peter H 10/13
    DNSAlgorithm* DNS::newAlgorithm(const FlowField& u, const ChebyCoeff& Ubase,
                                    Real nu,Real dt, const DNSFlags& flags,Real t,bool controlled){
        DNSAlgorithm* alg = 0;
        switch (flags.timestepping) {
            case CNFE1:
            case SBDF1:
            case SBDF2:
            case SBDF3:
            case SBDF4:
                alg = new MultistepDNS(u, Ubase, nu, dt, flags, t, controlled);
                break;
            case CNRK2:
                alg = new RungeKuttaDNS(u, Ubase, nu, dt, flags, t, controlled);
                break;
            case SMRK2:
            case CNAB2:
                alg = new CNABstyleDNS(u, Ubase, nu, dt, flags, t, controlled);
                break;
            default:
                cerr << "DNS::newAlgorithm : algorithm " << flags.timestepping
                << " is unimplemented" << endl;
        }
        return alg;
    }

    
    ////////////////////////////////////////////////////
    
    
    
    
    
    DNS::~DNS() {
        delete main_algorithm_;
        delete init_algorithm_;
    }
    
    DNS& DNS::operator=(const DNS& dns) {
        delete main_algorithm_;
        delete init_algorithm_;
        main_algorithm_ = dns.main_algorithm_ ? dns.main_algorithm_->clone() : 0;
        init_algorithm_ = dns.init_algorithm_ ? dns.init_algorithm_->clone() : 0;
        return *this;
    }
    
    void DNS::advance(FlowField& u, FlowField& q, int Nsteps) {
        assert(main_algorithm_);
        
        // Error check
        if (!main_algorithm_->full() && !init_algorithm_) {
            cerr << "DNS::advance(u,q,Nsteps) : the main algorithm is uninitialized,\n"
            << "and the initialization algorithm is not set. This should not be\n"
            << "possible. Please submit a bug report (see documentation)."
            << endl;
            exit(1);
        }
        if (!q.geomCongruent(u))
            q.resize(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b());
        
        int n=0;
        while (!main_algorithm_->full() && n<Nsteps) {
	  init_algorithm_->advance(u,q,1);
            main_algorithm_->push(u);
            if (main_algorithm_->full()) {
                delete init_algorithm_;
                init_algorithm_ = 0;
            }
            ++n;
        }
        main_algorithm_->advance(u,q,Nsteps-n);
        main_algorithm_->project();
        u.project(flags().symmetries);
        q.project(flags().symmetries);
    }
    
    //void DNS::reset() {
    //}
    
    void DNS::advance_NL(FlowField& u, FlowField& q, FlowField& F, int Nsteps) {
        assert(main_algorithm_);
        
        // Error check
        if (!main_algorithm_->full() && !init_algorithm_) {
            cerr << "DNS::advance(u,q,Nsteps) : the main algorithm is uninitialized,\n"
            << "and the initialization algorithm is not set. This should not be\n"
            << "possible. Please submit a bug report (see documentation)."
            << endl;
            exit(1);
        }
        if (!q.geomCongruent(u))
            q.resize(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b());
        
        int n=0;
        while (!main_algorithm_->full() && n<Nsteps) {
	  init_algorithm_->advance_NL(u,q,F,1);
            main_algorithm_->push(u);
            if (main_algorithm_->full()) {
                delete init_algorithm_;
                init_algorithm_ = 0;
            }
            ++n;
        }
        main_algorithm_->advance_NL(u,q,F,Nsteps-n);
        main_algorithm_->project();
        u.project(flags().symmetries);
        q.project(flags().symmetries);
    }
    /************************************************************************    
     ***********************************************************************
     
     Written by: Peter H Heins (Postgraduate Research Student) 
     
     Date: January 2012
     
     Institution: University of Sheffield (ACSE Department)
     
     Purpose: This function advances the u and q fields but incorporates inhomogeneous BCs
     
     Based on code written by Binh Lieu, University of Minnesota March 2009
     
     *********************************************************************** 
     *************************************************************************/
    
  void DNS::advance_inhom(FlowField& u, FlowField& q, FlowField& BCs, int Nsteps) {
        
        assert(main_algorithm_);
        
        // Error check
        if (!main_algorithm_->full() && !init_algorithm_) {
            cerr << "DNS::advance(u,q,Nsteps) : the main algorithm is uninitialized,\n"
            << "and the initialization algorithm is not set. This should not be\n"
            << "possible. Please submit a bug report (see documentation)."
            << endl;
            exit(1);
        }
        if (!q.geomCongruent(u))
            q.resize(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b());
        
        int n=0;
        while (!main_algorithm_->full() && n<Nsteps) {
            init_algorithm_->advance_inhom(u,q,BCs,1);
            main_algorithm_->push(u);
            if (main_algorithm_->full()) {
                delete init_algorithm_;
                init_algorithm_ = 0;
            }
            ++n;
        }
        main_algorithm_->advance_inhom(u,q,BCs,Nsteps-n);
        main_algorithm_->project();
        u.project(flags().symmetries);
        q.project(flags().symmetries);
    }
    
    
    
    
    
    
  void DNS::advance_inhom_CON(Controller& controller, FlowField& u, FlowField& q, FlowField& BCs, FlowField& F, int Nsteps, double*** IO,  double*** CStateMat) {
        
        assert(main_algorithm_);
        
        // Error check
        if (!main_algorithm_->full() && !init_algorithm_) {
            cerr << "DNS::advance(u,q,Nsteps) : the main algorithm is uninitialized,\n"
            << "and the initialization algorithm is not set. This should not be\n"
            << "possible. Please submit a bug report (see documentation)."
            << endl;
            exit(1);
        }
        if (!q.geomCongruent(u))
            q.resize(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b());
        
        int n=0;
        while (!main_algorithm_->full() && n<Nsteps) {
	  init_algorithm_->advance_inhom_CON(controller,u,q,BCs,F,1,IO,CStateMat);
            main_algorithm_->push(u);
            if (main_algorithm_->full()) {
                delete init_algorithm_;
                init_algorithm_ = 0;
            }
            ++n;
        }
        main_algorithm_->advance_inhom_CON(controller,u,q,BCs,F,Nsteps-n,IO,CStateMat);
        main_algorithm_->project();
        u.project(flags().symmetries);
        q.project(flags().symmetries);
    }

    /*
  void DNS::advance_inhom_CON_SF(Controller& controller, FlowField& u, FlowField& q, FlowField& BCs, FlowField& F, int Nsteps, double*** IO) {
        
        assert(main_algorithm_);
        
        // Error check
        if (!main_algorithm_->full() && !init_algorithm_) {
            cerr << "DNS::advance(u,q,Nsteps) : the main algorithm is uninitialized,\n"
            << "and the initialization algorithm is not set. This should not be\n"
            << "possible. Please submit a bug report (see documentation)."
            << endl;
            exit(1);
        }
        if (!q.geomCongruent(u))
            q.resize(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b());
        
        int n=0;
        while (!main_algorithm_->full() && n<Nsteps) {
	  init_algorithm_->advance_inhom_CON_SF(controller,u,q,BCs,F,1,IO);
            main_algorithm_->push(u);
            if (main_algorithm_->full()) {
                delete init_algorithm_;
                init_algorithm_ = 0;
            }
            ++n;
        }
        main_algorithm_->advance_inhom_CON_SF(controller,u,q,BCs,F,Nsteps-n,IO);
        main_algorithm_->project();
        u.project(flags().symmetries);
        q.project(flags().symmetries);
    }
*/




    /**************************************************************************
     
     ---End: Peter H Heins---
     
     **************************************************************************/  
    
    
    
    
    void DNS::project() {
        if (init_algorithm_)
            init_algorithm_->project();
        if (main_algorithm_)
            main_algorithm_->project();
    }
    
    void DNS::operator*=(const FieldSymmetry& sigma) {
        if (init_algorithm_)
            *init_algorithm_ *= sigma;
        if (main_algorithm_)
            *main_algorithm_ *= sigma;
    }
    
    void DNS::reset_dt(Real dt) {
        assert(main_algorithm_);
        main_algorithm_->reset_dt(dt);
        
        DNSFlags mainflags = main_algorithm_->flags();
        if (!main_algorithm_->full() &&
            mainflags.initstepping != mainflags.timestepping) {
            
            DNSFlags initflags = mainflags;
            initflags.timestepping = mainflags.initstepping;
            
            // An initialization algorithm needs only u's parameters at construction,
            // not its data, so we can construct it from a zero-valued u of right size
            FlowField u(main_algorithm_->Nx(), main_algorithm_->Ny(),
                        main_algorithm_->Nz(), 3,
                        main_algorithm_->Lx(), main_algorithm_->Lz(),
                        main_algorithm_->a(),  main_algorithm_->b());
            
            
            if (init_algorithm_)
                delete init_algorithm_;
            
            init_algorithm_ = newAlgorithm(u,
                                           main_algorithm_->Ubase(),
                                           main_algorithm_->nu(), dt, initflags,
                                           main_algorithm_->time(), false);
            
            // Safety check
            if (init_algorithm_->Ninitsteps() != 0)
                cerr << "DNS::DNS(u, Ubase, nu, dt, flags, t) :\n"
                << mainflags.initstepping  << " can't initialize "
                << mainflags.timestepping
                << " since it needs initialization itself.\n";
        }
        // No need for safety check on init_algorithm, has already been ok'd in ctor
    }
    
    //*******************************************************************************
    //**************** Added by P.Heins *********************************************
    //*******************************************************************************
    
    void DNS::reset_dtIH(Real dt, FlowField& BCs) {
        assert(main_algorithm_);
        main_algorithm_->reset_dtIH(dt,BCs);
        
        DNSFlags mainflags = main_algorithm_->flags();
        if (!main_algorithm_->full() &&
            mainflags.initstepping != mainflags.timestepping) {
            
            DNSFlags initflags = mainflags;
            initflags.timestepping = mainflags.initstepping;
            
            // An initialization algorithm needs only u's parameters at construction,
            // not its data, so we can construct it from a zero-valued u of right size
            FlowField u(main_algorithm_->Nx(), main_algorithm_->Ny(),
                        main_algorithm_->Nz(), 3,
                        main_algorithm_->Lx(), main_algorithm_->Lz(),
                        main_algorithm_->a(),  main_algorithm_->b());
            
            
            if (init_algorithm_)
                delete init_algorithm_;
            
            init_algorithm_ = newAlgorithm(u,
                                           main_algorithm_->Ubase(),
                                           main_algorithm_->nu(), dt, initflags,
                                           main_algorithm_->time(), true);
            
            // Safety check
            if (init_algorithm_->Ninitsteps() != 0)
                cerr << "DNS::DNS(u, Ubase, nu, dt, flags, t) :\n"
                << mainflags.initstepping  << " can't initialize "
                << mainflags.timestepping
                << " since it needs initialization itself.\n";
        }
        // No need for safety check on init_algorithm, has already been ok'd in ctor
    }

    //*********************************************************************************
    //*********************************************************************************
    
    // The mindless hassle of wrapper classes in C++ follows
    void DNS::reset_time(Real t) {
        if (init_algorithm_)
            init_algorithm_->reset_time(t);
        if (main_algorithm_)
            main_algorithm_->reset_time(t);
    }
    void DNS::reset_dPdx(Real dPdx) {
        if (init_algorithm_)
            init_algorithm_->reset_dPdx(dPdx);
        if (main_algorithm_)
            main_algorithm_->reset_dPdx(dPdx);
    }
    void DNS::reset_Ubulk(Real Ubulk) {
        if (init_algorithm_)
            init_algorithm_->reset_Ubulk(Ubulk);
        if (main_algorithm_)
            main_algorithm_->reset_Ubulk(Ubulk);
    }
    //void DNS::reset_uj(const FlowField& uj, int j) {
    //assert(main_algorithm_);
    //main_algorithm_->reset_uj(uj, j);
    //}
    bool DNS::push(const FlowField& u) {
        if (main_algorithm_)
            return main_algorithm_->push(u);
        else
            return false;
    }
    bool DNS::full() const {
        if (main_algorithm_)
            return main_algorithm_->full();
        else
            return false;
    }
    int DNS::order() const {
        if (main_algorithm_)
            return main_algorithm_->order();
        else if (init_algorithm_)
            return init_algorithm_->order();
        else
            return 0;
    }
    
    int DNS::Ninitsteps() const {
        if (main_algorithm_)
            return main_algorithm_->Ninitsteps();
        else
            return 0;
    }
    
    Real DNS::nu() const {
        if (main_algorithm_)
            return main_algorithm_->nu();
        else if (init_algorithm_)
            return init_algorithm_->nu();
        else
            return 0.0;
    }
    Real DNS::dt() const {
        if (main_algorithm_)
            return main_algorithm_->dt();
        else if (init_algorithm_)
            return init_algorithm_->dt();
        else
            return 0.0;
    }
    Real DNS::CFL() const {
        if (main_algorithm_)
            return main_algorithm_->CFL();
        else if (init_algorithm_)
            return init_algorithm_->CFL();
        else
            return 0.0;
    }
    Real DNS::time() const {
        if (main_algorithm_)
            return main_algorithm_->time();
        else if (init_algorithm_)
            return init_algorithm_->time();
        else
            return 0.0;
    }
    Real DNS::dPdx() const {
        if (main_algorithm_)
            return main_algorithm_->dPdx();
        else if (init_algorithm_)
            return init_algorithm_->dPdx();
        else
            return 0.0;
    }
    Real DNS::Ubulk() const {
        if (main_algorithm_)
            return main_algorithm_->Ubulk();
        else if (init_algorithm_)
            return init_algorithm_->Ubulk();
        else
            return 0.0;
    }
    Real DNS::dPdxRef() const {
        if (main_algorithm_)
            return main_algorithm_->dPdxRef();
        else if (init_algorithm_)
            return init_algorithm_->dPdxRef();
        else
            return 0.0;
    }
    Real DNS::UbulkRef() const {  // the bulk velocity enforced during integ.
        if (main_algorithm_)
            return main_algorithm_->UbulkRef();
        else if (init_algorithm_)
            return init_algorithm_->UbulkRef();
        else
            return 0.0;
    }
    const DNSFlags& DNS::flags() const {
        if (main_algorithm_)
            return main_algorithm_->flags();
        else if (init_algorithm_)
            return init_algorithm_->flags();
        else {
            cerr << "Error in DNS::flags(): flags are currently undefined" << endl;
            exit(1);
            return init_algorithm_->flags(); // to make compiler happy
        }
    }
    TimeStepMethod DNS::timestepping() const {
        if (main_algorithm_)
            return main_algorithm_->timestepping();
        else if (init_algorithm_)
            return init_algorithm_->timestepping();
        else
            return CNFE1;
    }
    
    void DNS::uq2p(const FlowField& u_, const FlowField& q_, FlowField& p) const {
        if (flags().nonlinearity != Rotational) {
            p = q_;
            return;
        }
        
        assert(u_.Nd() == 3);
        assert(q_.Nd() == 1);
        assert(main_algorithm_);
        FlowField& u = const_cast<FlowField&>(u_);
        FlowField& q = const_cast<FlowField&>(q_);
        ChebyCoeff& U = const_cast<ChebyCoeff&>(main_algorithm_->Ubase());
        
        fieldstate uxzstate = u.xzstate();
        fieldstate uystate = u.ystate();
        fieldstate qxzstate = q.xzstate();
        fieldstate qystate = q.ystate();
        fieldstate Ustate = U.state();
        
        u.makePhysical();
        q.makePhysical();
        U.makePhysical();
        
        int Nx=u.Nx();
        int Ny=u.Ny();
        int Nz=u.Nz();
        
        // Set p = q - 1/2 (u+U) dot (u+U)
        p = q;
        for (int ny=0; ny<Ny; ++ny) {
            Real Uny = U(ny);
            for (int nx=0; nx<Nx; ++nx)
                for (int nz=0; nz<Nz; ++nz)
                    p(nx,ny,nz,0) -= 0.5*(square(u(nx,ny,nz,0) + Uny) +
                                          square(u(nx,ny,nz,1)) +
                                          square(u(nx,ny,nz,2)));
        }
        u.makeState(uxzstate, uystate);
        U.makeState(Ustate);
        q.makeState(qxzstate, qystate);
        p.makeState(qxzstate, qystate);
    }
    
    void DNS::up2q(const FlowField& u_, const FlowField& p_, FlowField& q) const{
        if (flags().nonlinearity != Rotational) {
            q = p_;
            return;
        }
        
        assert(main_algorithm_);
        assert(u_.Nd() == 3);
        assert(p_.Nd() == 1);
        FlowField& u = const_cast<FlowField&>(u_);
        FlowField& p = const_cast<FlowField&>(p_);
        ChebyCoeff& U = const_cast<ChebyCoeff&>(main_algorithm_->Ubase());
        
        fieldstate uxzstate = u.xzstate();
        fieldstate uystate = u.ystate();
        fieldstate pxzstate = p.xzstate();
        fieldstate pystate = p.ystate();
        fieldstate Ustate = U.state();
        
        u.makePhysical();
        p.makePhysical();
        U.makePhysical();
        
        int Nx=u.Nx();
        int Ny=u.Ny();
        int Nz=u.Nz();
        
        // Set q = p + 1/2 (u+U) dot (u+U) to q
        q.makePhysical();
        for (int ny=0; ny<Ny; ++ny) {
            Real Uny = U(ny);
            for (int nx=0; nx<Nx; ++nx)
                for (int nz=0; nz<Nz; ++nz)
                    q(nx,ny,nz,0) += 0.5*(square(u(nx,ny,nz,0) + Uny) +
                                          square(u(nx,ny,nz,1)) +
                                          square(u(nx,ny,nz,2)));
        }
        u.makeState(uxzstate, uystate);
        U.makeState(Ustate);
        p.makeState(pxzstate, pystate);
        q.makeState(pxzstate, pystate);
    }
    
    void DNS::printStack() const {
        assert(main_algorithm_);
        main_algorithm_->printStack();
    }
    
    
    //====================================================================
    
    DNSAlgorithm::~DNSAlgorithm() {}
    
    DNSAlgorithm::DNSAlgorithm()
    :
    Nx_(0),
    Ny_(0),
    Nz_(0),
    Mx_(0),
    Mz_(0),
    Nyd_(0),
    kxd_max_(0),
    kzd_max_(0),
    Lx_(0),
    Lz_(0),
    a_(0),
    b_(0),
    flags_(),
    order_(0),
    Ninitsteps_(0),
    nu_(0),
    dt_(0),
    t_(0),
    cfl_(0),
    dPdxRef_(0),
    dPdxAct_(0),
    UbulkRef_(0),
    UbulkAct_(0),
    UbulkBase_(0),
    ubulkBase_(0),
    Ubase_(),
    Ubaseyy_(),
    ubase_(),
    ubtot_(),
    tmp2_(),
    //lapl_ubase_(),
    //nonl_ubase_(),
    tmp_(),
    uk_(),
    vk_(),
    wk_(),
    Pk_(),
    Pyk_(),
    Rxk_(),
    Ryk_(),
    Rzk_()
    
    {}
    
    DNSAlgorithm::DNSAlgorithm(const DNSAlgorithm& d)
    :
    Nx_(d.Nx_),
    Ny_(d.Ny_),
    Nz_(d.Nz_),
    Mx_(d.Mx_),
    Mz_(d.Mz_),
    Nyd_(d.Nyd_),
    kxd_max_(d.kxd_max_),
    kzd_max_(d.kzd_max_),
    Lx_(d.Lx_),
    Lz_(d.Lz_),
    a_(d.a_),
    b_(d.b_),
    flags_(d.flags_),
    order_(d.order_),
    Ninitsteps_(d.Ninitsteps_),
    nu_(d.nu_),
    dt_(d.dt_),
    t_(d.t_),
    cfl_(d.cfl_),
    dPdxRef_(d.dPdxRef_),
    dPdxAct_(d.dPdxAct_),
    UbulkRef_(d.UbulkRef_),
    UbulkAct_(d.UbulkAct_),
    UbulkBase_(d.UbulkBase_),
    ubulkBase_(d.ubulkBase_),
    Ubase_(d.Ubase_),
    Ubaseyy_(d.Ubaseyy_),
    ubase_(d.ubase_),
    ubtot_(d.ubtot_),
    tmp2_(d.tmp2_),
    //lapl_ubase_(d.lapl_ubase_),
    //nonl_ubase_(d.nonl_ubase_),
    tmp_(d.tmp_),
    uk_(d.uk_),
    vk_(d.vk_),
    wk_(d.wk_),
    Pk_(d.Pk_),
    Pyk_(d.Pyk_),
    Rxk_(d.Rxk_),
    Ryk_(d.Ryk_),
    Rzk_(d.Rzk_)
    {}
    
    DNSAlgorithm::DNSAlgorithm(const FlowField& u, const ChebyCoeff& Ubase,
                               Real nu, Real dt, const DNSFlags& flags, Real t0, bool controlled)
    :
    Nx_(u.Nx()),
    Ny_(u.numYmodes()),
    Nz_(u.Nz()),
    Mx_(u.numXmodes()),
    Mz_(u.numZmodes()),
    Nyd_(flags.dealias_y() ? 2*(u.numYmodes()-1)/3 + 1 : u.numYmodes()),
    kxd_max_(flags.dealias_xz() ? u.Nx()/3-1 : u.kxmax()),
    kzd_max_(flags.dealias_xz() ? u.Nz()/3-1 : u.kzmax()),
    Lx_(u.Lx()),
    Lz_(u.Lz()),
    a_(u.a()),
    b_(u.b()),
    flags_(flags),
    order_(0),
    Ninitsteps_(0),
    nu_(nu),
    dt_(dt),
    t_(t0),
    cfl_(0),
    dPdxRef_(0),
    dPdxAct_(0),
    UbulkRef_(0),
    UbulkAct_(0),
    UbulkBase_(0),
    ubulkBase_(0),
    Ubase_(Ubase),
    Ubaseyy_(),
    ubase_(),
    ubtot_(),
    tmp2_(),
    //lapl_ubase_(),
    //nonl_ubase_(),
    tmp_(),
    uk_(Ny_,a_,b_,Spectral),
    vk_(Ny_,a_,b_,Spectral),
    wk_(Ny_,a_,b_,Spectral),
    Pk_(Ny_,a_,b_,Spectral),
    Pyk_(Ny_,a_,b_,Spectral),
    Rxk_(Ny_,a_,b_,Spectral),
    Ryk_(Ny_,a_,b_,Spectral),
    Rzk_(Ny_,a_,b_,Spectral)
    {
        
       
        
        //u.makeSpectral();
        assert(u.vectorDim() == 3);
        
        // These methods require a 9d (3x3) tmp flowfield
        if (flags_.nonlinearity ==  Alternating ||
            flags_.nonlinearity ==  Alternating_ ||
            flags_.nonlinearity ==  Convection ||
            flags_.nonlinearity ==  LinearAboutField ||
            flags_.nonlinearity ==  LinearAboutProfile ||
            flags_.nonlinearity ==  Divergence ||
            flags_.nonlinearity ==  SkewSymmetric) {
            tmp_.resize(u.Nx(), u.Ny(), u.Nz(), 9, u.Lx(), u.Lz(), u.a(), u.b());
        }
        else
            tmp_.resize(u.Nx(), u.Ny(), u.Nz(), 3, u.Lx(), u.Lz(), u.a(), u.b());
        
        // Set the aliased modes to zero.
        // (Obsolete, I think, Try removing 2007-03-14).
        
        //if (flags_.dealias_xz())
        //u.zeroAliasedModes();
        
        // =====================================================================
        // Resolution of ubase, Ubase, and flags.baseflow
        // Construct Ubase specified by flags.baseflow.
        
        ChebyCoeff UbaseFlag(Ny_, a_, b_, Spectral);
        switch(flags_.baseflow) {
            case PlaneCouette :
                UbaseFlag[1] = 1;     // Ubase = y;
                break;
            case Parabolic :
                UbaseFlag[0] =  0.5;  // Ubase = 1 - y^2;
                UbaseFlag[2] = -0.5;
                break;
            case Zero:
            default:
                break;
        }
        
        if (Ubase_.N() != 0 && flags_.baseflow != Zero && L2Dist(UbaseFlag, Ubase_) > 1e-14) {
            cerr << "warning in DNSAlgorithm::DNSAlgorithm : \n"
            << "discrepancy between explicitly-set Ubase and flags.baseflow\n";
            cerr << "L2Dist(UbaseFlag, Ubase_) == " << L2Dist(UbaseFlag, Ubase_) << endl;
            cerr << "Setting Ubase to value specified by flags.baseflow." << endl;
            //Ubase_.save("Ubase");
            //UbaseFlag.save("Uflag");
            Ubase_ = UbaseFlag;
        }
        
        if (flags_.nonlinearity == LinearAboutField) {
            ubase_  = u;
            ubase_.makeSpectral();
            
            ubtot_  = ubase_;
            ubtot_ += Ubase_;
            grad(ubtot_, tmp2_);
            
            //ubtot_.makePhysical();
            //tmp2_.makePhysical();
            
            ubulkBase_ = Re(ubase_.profile(0,0,0)).mean();
            
            //lapl(ubase_, lapl_ubase_);
            //convectionNL(ubase_, nonl_ubase_, tmp_);
            //nonl_ubase_.makeSpectral();
            //lapl_ubase_.makeSpectral();
        }
        
        // Calculate Ubaseyy_ and related quantities
        UbulkBase_ = Ubase_.mean();
        ChebyCoeff Ubasey = diff(Ubase_);
        Ubaseyy_ = diff(Ubasey);
        
        cfl_ = u.CFLfactor(Ubase_);
        cfl_ *= flags_.dealias_xz() ? 2.0*pi/3.0*dt_ : pi*dt_;
        
        // Determine actual Ubulk and dPdx from initial data Ubase + u.
        ChebyCoeff u00(Ny_,a_,b_,Spectral);
        for (int ny=0; ny<Ny_; ++ny)
            u00[ny] = Re(u.cmplx(0,ny,0,0));
        ChebyCoeff du00dy = diff(u00);
        
        UbulkAct_ = UbulkBase_ + ubulkBase_ + u00.mean();
        dPdxAct_  = nu*(du00dy.eval_b() - du00dy.eval_a())/(b_-a_);
        if (Ubase_.length() != 0)
            dPdxAct_  += nu*(Ubasey.eval_b() - Ubasey.eval_a())/(b_-a_);
        
        if (flags_.constraint == BulkVelocity)
            UbulkRef_ = flags_.Ubulk;
        else {
            dPdxAct_ = flags_.dPdx;
            dPdxRef_ = flags_.dPdx;
        }
    }
    
    
    //DNSAlgorithm& DNSAlgorithm::operator=(const DNSAlgorithm& dns) {
    //  os << "No No!" << endl;
    //  exit(1);
    //}
    void DNSAlgorithm::project() {}
    void DNSAlgorithm::operator*=(const FieldSymmetry& symm) {}
    
    bool DNSAlgorithm::push(const FlowField& u) {return true;}
    bool DNSAlgorithm::full() const {return true;}
    
    //void DNSAlgorithm::reset_uj(const FlowField& uj, int j) {;}
    
    void DNSAlgorithm::reset_time(Real t) {t_=t;}
    
    void DNSAlgorithm::reset_dPdx(Real dPdx) {
        flags_.constraint = PressureGradient;
        flags_.dPdx = dPdx;
        flags_.Ubulk = 0.0;
        dPdxRef_ = dPdx;
        UbulkRef_ = 0.0;
    }
    void DNSAlgorithm::reset_Ubulk(Real Ubulk) {
        flags_.constraint = BulkVelocity;
        flags_.Ubulk = Ubulk;
        flags_.dPdx = 0.0;
        UbulkRef_ = Ubulk;
        dPdxRef_ = 0.0;
    }
    
    int DNSAlgorithm::Nx() const {return Nx_;}
    int DNSAlgorithm::Ny() const {return Ny_;}
    int DNSAlgorithm::Nz() const {return Nz_;}
    Real DNSAlgorithm::Lx() const {return Lx_;}
    Real DNSAlgorithm::Lz() const {return Lz_;}
    Real DNSAlgorithm::a() const {return a_;}
    Real DNSAlgorithm::b() const {return b_;}
    Real DNSAlgorithm::dt() const {return dt_;}
    Real DNSAlgorithm::nu() const {return nu_;}
    Real DNSAlgorithm::CFL() const {return cfl_;}
    Real DNSAlgorithm::time() const {return t_;}
    Real DNSAlgorithm::dPdx() const {return dPdxAct_;}
    Real DNSAlgorithm::dPdxRef() const {return dPdxRef_;}
    Real DNSAlgorithm::Ubulk() const {return UbulkAct_;}
    Real DNSAlgorithm::UbulkRef() const {return UbulkRef_;}
    int DNSAlgorithm::order() const {return order_;}
    int DNSAlgorithm::Ninitsteps() const {return Ninitsteps_;}
    const DNSFlags& DNSAlgorithm::flags() const {return flags_;}
    const ChebyCoeff&  DNSAlgorithm::Ubase() const {return Ubase_;}
    const FlowField&   DNSAlgorithm::ubase() const {return ubase_;}
    TimeStepMethod DNSAlgorithm::timestepping() const {return flags_.timestepping;}
    int DNSAlgorithm::kxmaxDealiased() const {return kxd_max_;}
    int DNSAlgorithm::kzmaxDealiased() const {return kzd_max_;}
    bool DNSAlgorithm::isAliasedMode(int kx, int kz) const {
        return (abs(kx) > kxd_max_ || (abs(kz) > kzd_max_)) ? true : false;
    }
    
    void DNSAlgorithm::printStack() const {
        //os << "DNSAlgorithm::printStack()" << endl;
    }
    
    // ====================================================================
    // Multistep algorithms
    MultistepDNS::MultistepDNS()
    :
    DNSAlgorithm()
    {}
    
    MultistepDNS::MultistepDNS(const MultistepDNS& dns)
    :
    DNSAlgorithm(dns),
    eta_(dns.eta_),
    alpha_(dns.alpha_),
    beta_(dns.beta_),
    u_(dns.u_),
    f_(dns.f_),
    countdown_(dns.countdown_)
    {
        // Copy tausolvers
        tausolver_ = new TauSolver*[Mx_];       // new #1
        for (int mx=0; mx<Mx_; ++mx) {
            tausolver_[mx] = new TauSolver[Mz_];  // new #2
            for (int mz=0; mz<Mz_; ++mz)
                tausolver_[mx][mz] = dns.tausolver_[mx][mz];
        }
    }
    
    MultistepDNS& MultistepDNS::operator=(const MultistepDNS& dns) {
        cerr << "MultistepDNS::operator=(const MultistepDNS& dns) unimplemented\n";
        exit(1);
    }
    
    
    /*
    MultistepDNS::MultistepDNS(const FlowField& u, const ChebyCoeff& Ubase,
                               Real nu, Real dt, const DNSFlags& flags, Real t)
    
    :
    DNSAlgorithm(u,Ubase,nu,dt,flags, t)
    {
        TimeStepMethod algorithm = flags.timestepping;
        switch (algorithm) {
            case CNFE1:
            case SBDF1:
                order_ = 1;
                eta_ = 1.0;
                alpha_.resize(order_);
                beta_.resize(order_);
                alpha_[0] = -1.0;
                beta_[0]  =  1.0;
                break;
            case SBDF2:
                order_ = 2;
                alpha_.resize(order_);
                beta_.resize(order_);
                eta_ = 1.5;
                alpha_[0] = -2.0; alpha_[1] =  0.5;
                beta_[0]  =  2.0;  beta_[1] = -1.0;
                break;
            case SBDF3:
                order_ = 3;
                alpha_.resize(order_);
                beta_.resize(order_);
                eta_ = 11.0/6.0;
                alpha_[0] = -3.0;  alpha_[1] = 1.5; alpha_[2] = -1.0/3.0;
                beta_[0]  =  3.0;   beta_[1] = -3.0; beta_[2] = 1.0;
                break;
            case SBDF4:
                order_ = 4;
                alpha_.resize(order_);
                beta_.resize(order_);
                eta_ = 25.0/12.0;
                alpha_[0] = -4.0; alpha_[1] =  3.0; alpha_[2] = -4.0/3.0; alpha_[3] = 0.25;
                beta_[0]  =  4.0;  beta_[1] = -6.0;  beta_[2] =  4.0;      beta_[3] = -1.0;
                break;
            default:
                cerr << "MultistepDNS::MultistepDNS(un,Ubase,nu,dt,flags,t0)\n"
                << "error: flags.timestepping == " << algorithm
                << "is a non-multistepping algorithm" << endl;
                exit(1);
        }
        
        // Configure tausolvers
        tausolver_ = new TauSolver*[Mx_];       // new #1
        for (int mx=0; mx<Mx_; ++mx)
            tausolver_[mx] = new TauSolver[Mz_];  // new #2
        
        reset_dt(dt_);
        
        // Initialize arrays of previous u's and f's
        FlowField tmp(u);
        tmp.setToZero();
        
        u_.resize(order_);
        f_.resize(order_);
        for (int j=0; j<order_; ++j) {
            u_[j] = tmp;
            f_[j] = tmp;
        }
        //if (order_ > 0)  // should always be true
        //u_[0] = u;
        cfl_ = u.CFLfactor(Ubase_);
        cfl_ *= flags_.dealias_xz() ? 2.0*pi/3.0*dt_ : pi*dt_;
        
        Ninitsteps_ = order_;
        countdown_ = Ninitsteps_;
    }
    */
    
    
    
    
    
    
    //////////////////////////////////////////////////////////////////
    // Added by P. Heins 10/13
    
    MultistepDNS::MultistepDNS(const FlowField& u, const ChebyCoeff& Ubase,
                               Real nu, Real dt, const DNSFlags& flags, Real t, bool controlled)
    
    :
    DNSAlgorithm(u,Ubase,nu,dt,flags, t, controlled)
    {
        TimeStepMethod algorithm = flags.timestepping;
        switch (algorithm) {
            case CNFE1:
            case SBDF1:
                order_ = 1;
                eta_ = 1.0;
                alpha_.resize(order_);
                beta_.resize(order_);
                alpha_[0] = -1.0;
                beta_[0]  =  1.0;
                break;
            case SBDF2:
                order_ = 2;
                alpha_.resize(order_);
                beta_.resize(order_);
                eta_ = 1.5;
                alpha_[0] = -2.0; alpha_[1] =  0.5;
                beta_[0]  =  2.0;  beta_[1] = -1.0;
                break;
            case SBDF3:
                order_ = 3;
                alpha_.resize(order_);
                beta_.resize(order_);
                eta_ = 11.0/6.0;
                alpha_[0] = -3.0;  alpha_[1] = 1.5; alpha_[2] = -1.0/3.0;
                beta_[0]  =  3.0;   beta_[1] = -3.0; beta_[2] = 1.0;
                break;
            case SBDF4:
                order_ = 4;
                alpha_.resize(order_);
                beta_.resize(order_);
                eta_ = 25.0/12.0;
                alpha_[0] = -4.0; alpha_[1] =  3.0; alpha_[2] = -4.0/3.0; alpha_[3] = 0.25;
                beta_[0]  =  4.0;  beta_[1] = -6.0;  beta_[2] =  4.0;      beta_[3] = -1.0;
                break;
            default:
                cerr << "MultistepDNS::MultistepDNS(un,Ubase,nu,dt,flags,t0)\n"
                << "error: flags.timestepping == " << algorithm
                << "is a non-multistepping algorithm" << endl;
                exit(1);
        }
        
        // Configure tausolvers
        tausolver_ = new TauSolver*[Mx_];       // new #1
        for (int mx=0; mx<Mx_; ++mx)
            tausolver_[mx] = new TauSolver[Mz_];  // new #2
        
        
        if (controlled){
            FlowField BC(Nx_,2,Nz_,3,Lx_,Lz_,a_,b_);
            BC.makeState(Spectral,Physical);
            reset_dtIH(dt_,BC);  
        }
        else {
            reset_dt(dt_);
        }
        
        
        // Initialize arrays of previous u's and f's
        FlowField tmp(u);
        tmp.setToZero();
        
        u_.resize(order_);
        f_.resize(order_);
        for (int j=0; j<order_; ++j) {
            u_[j] = tmp;
            f_[j] = tmp;
        }
        //if (order_ > 0)  // should always be true
        //u_[0] = u;
        cfl_ = u.CFLfactor(Ubase_);
        cfl_ *= flags_.dealias_xz() ? 2.0*pi/3.0*dt_ : pi*dt_;
        
        Ninitsteps_ = order_;
        countdown_ = Ninitsteps_;
    }
    //////////////////////////////////////////////////////////////////
    
    
    
    
    
    
    
    
    
    
    MultistepDNS::~MultistepDNS() {
        if (tausolver_) {
            for (int mx=0; mx<Mx_; ++mx)
                delete[] tausolver_[mx];  // undo new #2
            delete[] tausolver_;        // undo new #1
        }
        tausolver_ = 0;
    }
    
    DNSAlgorithm* MultistepDNS::clone() const {
        return new MultistepDNS(*this);
    }
    
    void MultistepDNS::reset_dt(Real dt) {
        
        
        cfl_ *= dt/dt_;
        //nu_ = nu;
        dt_ = dt;
        const Real c = 4.0*square(pi)*nu_;
        const int kxmax = tmp_.kxmax();
        const int kzmax = tmp_.kzmax();
        
        // This loop replaces the TauSolver objects at tausolver_[substep][mx][mz]
        // with new TauSolver objects, with the given parameters.
        for (int mx=0; mx<Mx_; ++mx) {
            int kx = tmp_.kx(mx);
            for (int mz=0; mz<Mz_; ++mz) {
                int kz = tmp_.kz(mz);
                Real lambda = eta_/dt_ + c*(square(kx/Lx_) + square(kz/Lz_));
                
                // When using dealiasing, some modes get set to zero, rather than
                // updated with momentum eqns. Don't initialize TauSolvers for these.
                if ((kx != kxmax && kz != kzmax) &&
                    (!flags_.dealias_xz() || !isAliasedMode(kx,kz)))
                    
                    tausolver_[mx][mz] = TauSolver(kx, kz, Lx_, Lz_, a_, b_, lambda,
                                                   nu_, Ny_, flags_.taucorrection);
            }
        }
        // Start from beginning on initialization
        countdown_ = Ninitsteps_;
    }
    
    //****************************************************************************
    //****************************************************************************
    void MultistepDNS::reset_dtIH(Real dt, FlowField& BCs) {
        cfl_ *= dt/dt_;
        
        
        //nu_ = nu;
        dt_ = dt;
        const Real c = 4.0*square(pi)*nu_;
        const int kxmax = tmp_.kxmax();
        const int kzmax = tmp_.kzmax();
        
        //cout << "Nyd_ = " << Nyd_ << endl;
        
        
        // This loop replaces the TauSolver objects at tausolver_[substep][mx][mz]
        // with new TauSolver objects, with the given parameters.
        for (int mx=0; mx<Mx_; ++mx) {
            int kx = tmp_.kx(mx);
            for (int mz=0; mz<Mz_; ++mz) {
                int kz = tmp_.kz(mz);
                
                Complex v_upper=BCs.cmplx(mx,1,mz,1);
                Complex v_lower=BCs.cmplx(mx,0,mz,1);
                                
                double vr_upper=real(v_upper);
                double vi_upper=imag(v_upper);
                double vr_lower=real(v_lower);
                double vi_lower=imag(v_lower);
                                
                 
                Real lambda = eta_/dt_ + c*(square(kx/Lx_) + square(kz/Lz_));
                
                // When using dealiasing, some modes get set to zero, rather than
                // updated with momentum eqns. Don't initialize TauSolvers for these.
                if ((kx != kxmax && kz != kzmax) &&
                    (!flags_.dealias_xz() || !isAliasedMode(kx,kz)))
                    
                                        
                    tausolver_[mx][mz] = TauSolver(kx,kz,Lx_,Lz_,a_,b_,lambda,nu_,Ny_, vr_lower, vi_lower, vr_upper, vi_upper, flags_.taucorrection);
                    
                    //tausolver_[mx][mz] = TauSolver(kx, kz, Lx_, Lz_, a_, b_, lambda,                                                   nu_, Nyd_, flags_.taucorrection);
            }
        }
        
                
        
        // Start from beginning on initialization
        countdown_ = Ninitsteps_;
    }

    //****************************************************************************
    //****************************************************************************
    
    // This calculation follows Peyret section 4.5.1(b) pg 131.
  void MultistepDNS::advance(FlowField& un, FlowField& qn, int Nsteps) {
        
        const int kxmax = un.kxmax();
        const int kzmax = un.kzmax();
        const int J = order_ -1 ;
        
        u_[0] = un;
        

        for (int step=0; step<Nsteps; ++step) {
            
            //*flags_.logstream << "Multistep::advance(...) step == " << step << " {" <<endl;
            //printStack();
            
            // Calculate nonlinearity
            if (order_ > 0)
                navierstokesNL(u_[0], ubase_, Ubase_, f_[0], tmp_, tmp2_, flags_.nonlinearity);
            
	   

            // Update each Fourier mode with time-stepping algorithm
            for (int mx=0; mx<Mx_; ++mx) {
                const int kx = un.kx(mx);
                
                for (int mz=0; mz<Mz_; ++mz) {
                    const int kz = un.kz(mz);
                    
                    // Zero out the aliased modes and break to next kx,kz
                    if ((kx == kxmax || kz == kzmax) ||
                        (flags_.dealias_xz() && isAliasedMode(kx,kz))) {
                        for (int ny=0; ny<Ny_; ++ny) {
                            u_[0].cmplx(mx,ny,mz,0) = 0.0;
                            u_[0].cmplx(mx,ny,mz,1) = 0.0;
                            u_[0].cmplx(mx,ny,mz,2) = 0.0;
                            qn.cmplx(mx,ny,mz,0) = 0.0;
                        }
                        break;
                    }
                    
                    // For nonaliased modes
                    Rxk_.setToZero();
                    Ryk_.setToZero();
                    Rzk_.setToZero();
                    
                    // Add up multistepping terms of linear and nonlinear terms
                    for (int j=0; j<order_; ++j) {
                        const Real a = -alpha_[j]/dt_;
                        const Real b = -beta_[j];
                        for (int ny=0; ny<Ny_; ++ny) {
                            Rxk_.add(ny, a*u_[j].cmplx(mx,ny,mz,0)+b*f_[j].cmplx(mx,ny,mz,0));
                            Ryk_.add(ny, a*u_[j].cmplx(mx,ny,mz,1)+b*f_[j].cmplx(mx,ny,mz,1));
                            Rzk_.add(ny, a*u_[j].cmplx(mx,ny,mz,2)+b*f_[j].cmplx(mx,ny,mz,2));
                        }
                    }
                    
                    // Solve the tau solutions
                    if (kx!=0 || kz!=0)
                        tausolver_[mx][mz].solve(uk_,vk_,wk_,Pk_, Rxk_,Ryk_,Rzk_);
                    
                    else { // kx,kz == 0,0
                        
                        if (Ubaseyy_.length() > 0)
                            for (int ny=0; ny<Ny_; ++ny)
                                Rxk_.re[ny] += nu_*Ubaseyy_[ny];   // Rx has addl'l term from Ubase
                        
                        if (flags_.constraint == PressureGradient) {
                            // pressure is supplied, put on RHS of tau eqn
                            Rxk_.re[0] -= dPdxRef_;
                            
                            // Solve the tau equations
                            tausolver_[mx][mz].solve(uk_, vk_, wk_, Pk_, Rxk_,Ryk_,Rzk_);
                            
                            // Bulk vel is free variable determined from soln of tau eqn
                            UbulkAct_ = UbulkBase_ + ubulkBase_ + uk_.re.mean();
                            dPdxAct_ = dPdxRef_;
                        }
                        else { // const bulk velocity
                            // bulk velocity is supplied, put on RHS of tau eqn
                            // dPdxAct (i.e. at next time step is solved for)
                            // constraint:    UbulkBase + ubulkBase + mean(u) = UbulkRef.
                            tausolver_[mx][mz].solve(uk_, vk_, wk_, Pk_, dPdxAct_,
                                                     Rxk_, Ryk_, Rzk_,
                                                     UbulkRef_ - UbulkBase_ - ubulkBase_);
                            UbulkAct_ = UbulkBase_ + ubulkBase_ + uk_.re.mean(); // should == UbulkRef_
                            //UbulkAct_ = UbulkRef_;
                        }
                    }
                    
                    // Load solutions into u_[J], the last element of u history, which is not needed anymore.
                    // Because of FFTW complex symmetries
                    // The 0,0 mode must be real.
                    // For Nx even, the kxmax,0 mode must be real
                    // For Nz even, the 0,kzmax mode must be real
                    // For Nx,Nz even, the kxmax,kzmax mode must be real
                    if ((kx == 0 && kz == 0) ||
                        (Nx_%2 == 0 && kx == kxmax && kz == 0) ||
                        (Nz_%2 == 0 && kz == kzmax && kx == 0) ||
                        (Nx_%2 == 0 && Nz_%2 == 0 && kx == kxmax && kz == kzmax)) {
                        
                        for (int ny=0; ny<Ny_; ++ny) {
                            u_[J].cmplx(mx,ny,mz,0) = Complex(Re(uk_[ny]), 0.0);
                            u_[J].cmplx(mx,ny,mz,1) = Complex(Re(vk_[ny]), 0.0);
                            u_[J].cmplx(mx,ny,mz,2) = Complex(Re(wk_[ny]), 0.0);
                            qn.cmplx(mx,ny,mz,0) = Complex(Re(Pk_[ny]), 0.0);
                        }
                    }
                    // The normal case, for general kx,kz
                    else
                        for (int ny=0; ny<Ny_; ++ny) {
                            u_[J].cmplx(mx,ny,mz,0) = uk_[ny];
                            u_[J].cmplx(mx,ny,mz,1) = vk_[ny];
                            u_[J].cmplx(mx,ny,mz,2) = wk_[ny];
                            qn.cmplx(mx,ny,mz,0) = Pk_[ny];
                        }
                    
                    // And now set the y-aliased modes to zero.
                    for (int ny=Nyd_; ny<Ny_; ++ny) {
                        u_[J].cmplx(mx,ny,mz,0) = 0.0;
                        u_[J].cmplx(mx,ny,mz,1) = 0.0;
                        u_[J].cmplx(mx,ny,mz,2) = 0.0;
                        qn.cmplx(mx,ny,mz,0) = 0.0;
                    }
                }
            }
            
            // The solution is stored in u_[J]. Shift entire u and f arrays in time
            // Ie shift u_[J] <- u_[J-1] <- ... <- u_[0] <- u_[J]
            for (int j=order_-1; j>0; --j) {
                swap(f_[j], f_[j-1]);
                swap(u_[j], u_[j-1]);
            }
            t_ += dt_;
            
            //printStack();
            //*flags_.logstream << "} Multistep::advance(...) step == " << step << " }" <<endl;
            
            if (flags_.verbosity == PrintTime || flags_.verbosity == PrintAll)
                *flags_.logstream << t_ << ' ' << flush;
            else if (flags_.verbosity == PrintTicks)
                *flags_.logstream << '.' << flush;
        }
        
        un = u_[0];
        qn.setPadded(flags_.dealias_xz());
        
        cfl_ = u_[0].CFLfactor(Ubase_);
        cfl_ *= flags_.dealias_xz() ? 2.0*pi/3.0*dt_ : pi*dt_;
        
        if (flags_.verbosity == PrintTime || flags_.verbosity == PrintAll ||
            flags_.verbosity == PrintTicks)
            *flags_.logstream << endl;
        
        return;
    }
    
    
    /************************************************************************
     ***********************************************************************
     
     Written by: Peter H Heins (Postgraduate Research Student)
     
     Date: January 2012
     
     Institution: University of Sheffield (ACSE Department)
     
     Purpose: This function implements imhomogeneous BCs into the Multistep algorithm
     
     Based on code written by Binh Lieu, University of Minnesota March 2009
     
     ***********************************************************************
     *************************************************************************/
    
  void MultistepDNS::advance_NL(FlowField& un, FlowField& qn, FlowField& Fn, int Nsteps) {
        
        const int kxmax = un.kxmax();
        const int kzmax = un.kzmax();
        const int J = order_ -1 ;
        
        u_[0] = un;
        

        for (int step=0; step<Nsteps; ++step) {
            
            //*flags_.logstream << "Multistep::advance(...) step == " << step << " {" <<endl;
            //printStack();
            
            // Calculate nonlinearity
            if (order_ > 0)
                navierstokesNL(u_[0], ubase_, Ubase_, f_[0], tmp_, tmp2_, flags_.nonlinearity);
            
	    //Fn = f_[0];

            // Update each Fourier mode with time-stepping algorithm
            for (int mx=0; mx<Mx_; ++mx) {
                const int kx = un.kx(mx);
                
                for (int mz=0; mz<Mz_; ++mz) {
                    const int kz = un.kz(mz);
                    
                    // Zero out the aliased modes and break to next kx,kz
                    if ((kx == kxmax || kz == kzmax) ||
                        (flags_.dealias_xz() && isAliasedMode(kx,kz))) {
                        for (int ny=0; ny<Ny_; ++ny) {
                            u_[0].cmplx(mx,ny,mz,0) = 0.0;
                            u_[0].cmplx(mx,ny,mz,1) = 0.0;
                            u_[0].cmplx(mx,ny,mz,2) = 0.0;
                            qn.cmplx(mx,ny,mz,0) = 0.0;
                        }
                        break;
                    }
                    
                    // For nonaliased modes
                    Rxk_.setToZero();
                    Ryk_.setToZero();
                    Rzk_.setToZero();
                    
                    // Add up multistepping terms of linear and nonlinear terms
                    for (int j=0; j<order_; ++j) {
                        const Real a = -alpha_[j]/dt_;
                        const Real b = -beta_[j];
                        for (int ny=0; ny<Ny_; ++ny) {
                            Rxk_.add(ny, a*u_[j].cmplx(mx,ny,mz,0)+b*f_[j].cmplx(mx,ny,mz,0));
                            Ryk_.add(ny, a*u_[j].cmplx(mx,ny,mz,1)+b*f_[j].cmplx(mx,ny,mz,1));
                            Rzk_.add(ny, a*u_[j].cmplx(mx,ny,mz,2)+b*f_[j].cmplx(mx,ny,mz,2));
                        }
                    }
                    
                    // Solve the tau solutions
                    if (kx!=0 || kz!=0)
                        tausolver_[mx][mz].solve(uk_,vk_,wk_,Pk_, Rxk_,Ryk_,Rzk_);
                    
                    else { // kx,kz == 0,0
                        
                        if (Ubaseyy_.length() > 0)
                            for (int ny=0; ny<Ny_; ++ny)
                                Rxk_.re[ny] += nu_*Ubaseyy_[ny];   // Rx has addl'l term from Ubase
                        
                        if (flags_.constraint == PressureGradient) {
                            // pressure is supplied, put on RHS of tau eqn
                            Rxk_.re[0] -= dPdxRef_;
                            
                            // Solve the tau equations
                            tausolver_[mx][mz].solve(uk_, vk_, wk_, Pk_, Rxk_,Ryk_,Rzk_);
                            
                            // Bulk vel is free variable determined from soln of tau eqn
                            UbulkAct_ = UbulkBase_ + ubulkBase_ + uk_.re.mean();
                            dPdxAct_ = dPdxRef_;
                        }
                        else { // const bulk velocity
                            // bulk velocity is supplied, put on RHS of tau eqn
                            // dPdxAct (i.e. at next time step is solved for)
                            // constraint:    UbulkBase + ubulkBase + mean(u) = UbulkRef.
                            tausolver_[mx][mz].solve(uk_, vk_, wk_, Pk_, dPdxAct_,
                                                     Rxk_, Ryk_, Rzk_,
                                                     UbulkRef_ - UbulkBase_ - ubulkBase_);
                            UbulkAct_ = UbulkBase_ + ubulkBase_ + uk_.re.mean(); // should == UbulkRef_
                            //UbulkAct_ = UbulkRef_;
                        }
                    }
                    
                    // Load solutions into u_[J], the last element of u history, which is not needed anymore.
                    // Because of FFTW complex symmetries
                    // The 0,0 mode must be real.
                    // For Nx even, the kxmax,0 mode must be real
                    // For Nz even, the 0,kzmax mode must be real
                    // For Nx,Nz even, the kxmax,kzmax mode must be real
                    if ((kx == 0 && kz == 0) ||
                        (Nx_%2 == 0 && kx == kxmax && kz == 0) ||
                        (Nz_%2 == 0 && kz == kzmax && kx == 0) ||
                        (Nx_%2 == 0 && Nz_%2 == 0 && kx == kxmax && kz == kzmax)) {
                        
                        for (int ny=0; ny<Ny_; ++ny) {
                            u_[J].cmplx(mx,ny,mz,0) = Complex(Re(uk_[ny]), 0.0);
                            u_[J].cmplx(mx,ny,mz,1) = Complex(Re(vk_[ny]), 0.0);
                            u_[J].cmplx(mx,ny,mz,2) = Complex(Re(wk_[ny]), 0.0);
                            qn.cmplx(mx,ny,mz,0) = Complex(Re(Pk_[ny]), 0.0);
                        }
                    }
                    // The normal case, for general kx,kz
                    else
                        for (int ny=0; ny<Ny_; ++ny) {
                            u_[J].cmplx(mx,ny,mz,0) = uk_[ny];
                            u_[J].cmplx(mx,ny,mz,1) = vk_[ny];
                            u_[J].cmplx(mx,ny,mz,2) = wk_[ny];
                            qn.cmplx(mx,ny,mz,0) = Pk_[ny];
                        }
                    
                    // And now set the y-aliased modes to zero.
                    for (int ny=Nyd_; ny<Ny_; ++ny) {
                        u_[J].cmplx(mx,ny,mz,0) = 0.0;
                        u_[J].cmplx(mx,ny,mz,1) = 0.0;
                        u_[J].cmplx(mx,ny,mz,2) = 0.0;
                        qn.cmplx(mx,ny,mz,0) = 0.0;
                    }
                }
            }
            
            // The solution is stored in u_[J]. Shift entire u and f arrays in time
            // Ie shift u_[J] <- u_[J-1] <- ... <- u_[0] <- u_[J]
            for (int j=order_-1; j>0; --j) {
                swap(f_[j], f_[j-1]);
                swap(u_[j], u_[j-1]);
            }
            t_ += dt_;
            
            //printStack();
            //*flags_.logstream << "} Multistep::advance(...) step == " << step << " }" <<endl;
            
            if (flags_.verbosity == PrintTime || flags_.verbosity == PrintAll)
                *flags_.logstream << t_ << ' ' << flush;
            else if (flags_.verbosity == PrintTicks)
                *flags_.logstream << '.' << flush;
        }
        
        un = u_[0];
	Fn = f_[0];  

        qn.setPadded(flags_.dealias_xz());
        
        cfl_ = u_[0].CFLfactor(Ubase_);
        cfl_ *= flags_.dealias_xz() ? 2.0*pi/3.0*dt_ : pi*dt_;
        
        if (flags_.verbosity == PrintTime || flags_.verbosity == PrintAll ||
            flags_.verbosity == PrintTicks)
            *flags_.logstream << endl;
        
        return;
    }
    


   
    
    void MultistepDNS::advance_inhom(FlowField& un, FlowField& qn, FlowField& BCs,
                                     int Nsteps) {
        
        const int kxmax = un.kxmax();
        const int kzmax = un.kzmax();
        const int J = order_ -1 ;
        
        u_[0] = un;
        
        
        for (int step=0; step<Nsteps; ++step) {
            
            //*flags_.logstream << "Multistep::advance(...) step == " << step << " {" <<endl;
            //printStack();
            
            // Calculate nonlinearity
            if (order_ > 0)
                navierstokesNL(u_[0], ubase_, Ubase_, f_[0], tmp_, tmp2_, flags_.nonlinearity);
            
            // Update each Fourier mode with time-stepping algorithm
            for (int mx=0; mx<Mx_; ++mx) {
                const int kx = un.kx(mx);
                
                for (int mz=0; mz<Mz_; ++mz) {
                    const int kz = un.kz(mz);
                    
                    //*************************
                    
                    Complex u_upper=BCs.cmplx(mx,1,mz,0);
                    Complex u_lower=BCs.cmplx(mx,0,mz,0);
                    Complex v_upper=BCs.cmplx(mx,1,mz,1);
                    Complex v_lower=BCs.cmplx(mx,0,mz,1);
                    Complex w_upper=BCs.cmplx(mx,1,mz,2);
                    Complex w_lower=BCs.cmplx(mx,0,mz,2);
                    
                    Real vr_upper=real(v_upper);
                    Real vi_upper=imag(v_upper);
                    Real vr_lower=real(v_lower);
                    Real vi_lower=imag(v_lower);
                    
                    const Real c = 4.0*square(pi)*nu_;
                    
                    Real lambda = eta_/dt_ + c*(square(kx/Lx_) + square(kz/Lz_));
                    
                    tausolver_[mx][mz] = TauSolver(kx,kz,Lx_,Lz_,a_,b_,lambda,nu_,Ny_, vr_lower, vi_lower, vr_upper, vi_upper, flags_.taucorrection);
                    
                    // ***************************
                    
                    // Zero out the aliased modes and break to next kx,kz
                    if ((kx == kxmax || kz == kzmax) ||
                        (flags_.dealias_xz() && isAliasedMode(kx,kz))) {
                        for (int ny=0; ny<Ny_; ++ny) {
                            u_[0].cmplx(mx,ny,mz,0) = 0.0;
                            u_[0].cmplx(mx,ny,mz,1) = 0.0;
                            u_[0].cmplx(mx,ny,mz,2) = 0.0;
                            qn.cmplx(mx,ny,mz,0) = 0.0;
                        }
                        break;
                    }
                    
                    /*
                    // Zero out y-aliased modes for all kx,kz
                    if (flags.dealias_y() || flags.dealias_xyz() {
                        for (int ny=iround(2*Nyd_)/3; ny<Nyd_; ++ny) {
                            u_[0].cmplx(mx,ny,mz,0) = 0.0;
                            u_[0].cmplx(mx,ny,mz,1) = 0.0;
                            u_[0].cmplx(mx,ny,mz,2) = 0.0;
                            qn.cmplx(mx,ny,mz,0) = 0.0;
                        }                    
                    }
                     */
                    
                    // For nonaliased modes
                    Rxk_.setToZero();
                    Ryk_.setToZero();
                    Rzk_.setToZero();
                    
                    // Add up multistepping terms of linear and nonlinear terms
                    for (int j=0; j<order_; ++j) {
                        const Real a = -alpha_[j]/dt_;
                        const Real b = -beta_[j];
                        for (int ny=0; ny<Ny_; ++ny) {
                            Rxk_.add(ny, a*u_[j].cmplx(mx,ny,mz,0)+b*f_[j].cmplx(mx,ny,mz,0));
                            Ryk_.add(ny, a*u_[j].cmplx(mx,ny,mz,1)+b*f_[j].cmplx(mx,ny,mz,1));
                            Rzk_.add(ny, a*u_[j].cmplx(mx,ny,mz,2)+b*f_[j].cmplx(mx,ny,mz,2));
                        }
                    }
                    
                    // Solve the tau solutions
                    if (kx!=0 || kz!=0)
                        tausolver_[mx][mz].solve_Inhom(uk_,vk_,wk_,Pk_, Rxk_,Ryk_,Rzk_, vr_lower, vr_upper, vi_lower, vi_upper, u_lower, u_upper, w_lower, w_upper);
                    
                    else { // kx,kz == 0,0
                        
                        if (Ubaseyy_.length() > 0)
                            for (int ny=0; ny<Ny_; ++ny)
                                Rxk_.re[ny] += nu_*Ubaseyy_[ny];   // Rx has addl'l term from Ubase
                        
                        if (flags_.constraint == PressureGradient) {
                            // pressure is supplied, put on RHS of tau eqn
                            Rxk_.re[0] -= dPdxRef_;
                            
                            // Solve the tau equations
                            //tausolver_[mx][mz].solve(uk_, vk_, wk_, Pk_, Rxk_,Ryk_,Rzk_);
                            tausolver_[mx][mz].solve_Inhom(uk_,vk_,wk_,Pk_, Rxk_,Ryk_,Rzk_, vr_lower, vr_upper, vi_lower, vi_upper, u_lower, u_upper, w_lower, w_upper);
                            
                            // Bulk vel is free variable determined from soln of tau eqn
                            UbulkAct_ = UbulkBase_ + ubulkBase_ + uk_.re.mean();
                            dPdxAct_ = dPdxRef_;
                        }
                        else { // const bulk velocity
                            // bulk velocity is supplied, put on RHS of tau eqn
                            // dPdxAct (i.e. at next time step is solved for)
                            // constraint:    UbulkBase + ubulkBase + mean(u) = UbulkRef.
                            tausolver_[mx][mz].solve_Inhom(uk_, vk_, wk_, Pk_, dPdxAct_,
                                                           Rxk_, Ryk_, Rzk_, vr_lower, vr_upper, vi_lower, vi_upper,                                                   UbulkRef_ - UbulkBase_ - ubulkBase_);
                            UbulkAct_ = UbulkBase_ + ubulkBase_ + uk_.re.mean(); // should == UbulkRef_
                            //UbulkAct_ = UbulkRef_;
                        }
                    }
                    
                    // Load solutions into u_[J], the last element of u history, which is not needed anymore.
                    // Because of FFTW complex symmetries
                    // The 0,0 mode must be real.
                    // For Nx even, the kxmax,0 mode must be real
                    // For Nz even, the 0,kzmax mode must be real
                    // For Nx,Nz even, the kxmax,kzmax mode must be real
                    if ((kx == 0 && kz == 0) ||
                        (Nx_%2 == 0 && kx == kxmax && kz == 0) ||
                        (Nz_%2 == 0 && kz == kzmax && kx == 0) ||
                        (Nx_%2 == 0 && Nz_%2 == 0 && kx == kxmax && kz == kzmax)) {
                        
                        for (int ny=0; ny<Ny_; ++ny) {
                            u_[J].cmplx(mx,ny,mz,0) = Complex(Re(uk_[ny]), 0.0);
                            u_[J].cmplx(mx,ny,mz,1) = Complex(Re(vk_[ny]), 0.0);
                            u_[J].cmplx(mx,ny,mz,2) = Complex(Re(wk_[ny]), 0.0);
                            qn.cmplx(mx,ny,mz,0) = Complex(Re(Pk_[ny]), 0.0);
                        }
                    }
                    // The normal case, for general kx,kz
                    else
                        for (int ny=0; ny<Ny_; ++ny) {
                            u_[J].cmplx(mx,ny,mz,0) = uk_[ny];
                            u_[J].cmplx(mx,ny,mz,1) = vk_[ny];
                            u_[J].cmplx(mx,ny,mz,2) = wk_[ny];
                            qn.cmplx(mx,ny,mz,0) = Pk_[ny];
                        }
                    
                    // And now set the y-aliased modes to zero.
                    for (int ny=Nyd_; ny<Ny_; ++ny) {
                        u_[J].cmplx(mx,ny,mz,0) = 0.0;
                        u_[J].cmplx(mx,ny,mz,1) = 0.0;
                        u_[J].cmplx(mx,ny,mz,2) = 0.0;
                        qn.cmplx(mx,ny,mz,0) = 0.0;
                    }
                }
            }
            
            
            
            // The solution is stored in u_[J]. Shift entire u and f arrays in time
            // Ie shift u_[J] <- u_[J-1] <- ... <- u_[0] <- u_[J]
            for (int j=order_-1; j>0; --j) {
                swap(f_[j], f_[j-1]);
                swap(u_[j], u_[j-1]);
            }
            t_ += dt_;
            
            //printStack();
            //*flags_.logstream << "} Multistep::advance(...) step == " << step << " }" <<endl;
            
            if (flags_.verbosity == PrintTime || flags_.verbosity == PrintAll)
                *flags_.logstream << t_ << ' ' << flush;
            else if (flags_.verbosity == PrintTicks)
                *flags_.logstream << '.' << flush;
        }
        
        un = u_[0];
        qn.setPadded(flags_.dealias_xz());
        
        cfl_ = u_[0].CFLfactor(Ubase_);
        cfl_ *= flags_.dealias_xz() ? 2.0*pi/3.0*dt_ : pi*dt_;
        
        if (flags_.verbosity == PrintTime || flags_.verbosity == PrintAll ||
            flags_.verbosity == PrintTicks)
            *flags_.logstream << endl;
        
        return;
    }
    
    //**************************************************************************
    
    
    
    
    
     void MultistepDNS::advance_inhom_CON(Controller& controller, FlowField& un, FlowField& qn, FlowField& BCs, FlowField& Fn,  int Nsteps, double*** IO, double*** CStateMat) {
        
        const int kxmax = un.kxmax();
        const int kzmax = un.kzmax();
        const int J = order_ -1 ;
        
        u_[0] = un;
       
                 
        
        for (int step=0; step<Nsteps; ++step) {
            
            

	  controller.advance_Con_CN(u_[0],qn,BCs,CStateMat,dt_,t_+dt_,IO);
            //cout << "dt=" << dt_ <<  " " << "t=" << t_ << " " << "BC" << BCs.cmplx(0,0,1,1) << endl;

                        
            
            //*flags_.logstream << "Multistep::advance(...) step == " << step << " {" <<endl;
            //printStack();
            
            // Calculate nonlinearity
            if (order_ > 0)
                navierstokesNL(u_[0], ubase_, Ubase_, f_[0], tmp_, tmp2_, flags_.nonlinearity);
            
	     Fn = f_[0];

            // Update each Fourier mode with time-stepping algorithm
            for (int mx=0; mx<Mx_; ++mx) {
                const int kx = un.kx(mx);
                
                for (int mz=0; mz<Mz_; ++mz) {
                    const int kz = un.kz(mz);
                    //*************************
                                       
                    
                    Complex u_upper=BCs.cmplx(mx,1,mz,0);
                    Complex u_lower=BCs.cmplx(mx,0,mz,0);
                    Complex v_upper=BCs.cmplx(mx,1,mz,1);
                    Complex v_lower=BCs.cmplx(mx,0,mz,1);
                    Complex w_upper=BCs.cmplx(mx,1,mz,2);
                    Complex w_lower=BCs.cmplx(mx,0,mz,2);
                    
                    Real vr_upper=real(v_upper);
                    Real vi_upper=imag(v_upper);
                    Real vr_lower=real(v_lower);
                    Real vi_lower=imag(v_lower);
                    
                    const Real c = 4.0*square(pi)*nu_;
                    
                    Real lambda = eta_/dt_ + c*(square(kx/Lx_) + square(kz/Lz_));
                    
                    tausolver_[mx][mz] = TauSolver(kx,kz,Lx_,Lz_,a_,b_,lambda,nu_,Ny_, vr_lower, vi_lower, vr_upper, vi_upper, flags_.taucorrection);
                    
                    // ***************************
                    
                    // Zero out the aliased modes and break to next kx,kz
                    if ((kx == kxmax || kz == kzmax) ||
                        (flags_.dealias_xz() && isAliasedMode(kx,kz))) {
                        for (int ny=0; ny<Ny_; ++ny) {
                            u_[0].cmplx(mx,ny,mz,0) = 0.0;
                            u_[0].cmplx(mx,ny,mz,1) = 0.0;
                            u_[0].cmplx(mx,ny,mz,2) = 0.0;
                            qn.cmplx(mx,ny,mz,0) = 0.0;
                        }
                        break;
                    }
                    
                    /*
                     // Zero out y-aliased modes for all kx,kz
                     if (flags.dealias_y() || flags.dealias_xyz() {
                     for (int ny=iround(2*Nyd_)/3; ny<Nyd_; ++ny) {
                     u_[0].cmplx(mx,ny,mz,0) = 0.0;
                     u_[0].cmplx(mx,ny,mz,1) = 0.0;
                     u_[0].cmplx(mx,ny,mz,2) = 0.0;
                     qn.cmplx(mx,ny,mz,0) = 0.0;
                     }                    
                     }
                     */
                    
                    // For nonaliased modes
                    Rxk_.setToZero();
                    Ryk_.setToZero();
                    Rzk_.setToZero();
                    
                    // Add up multistepping terms of linear and nonlinear terms
                    for (int j=0; j<order_; ++j) {
                        const Real a = -alpha_[j]/dt_;
                        const Real b = -beta_[j];
                        for (int ny=0; ny<Ny_; ++ny) {
                            Rxk_.add(ny, a*u_[j].cmplx(mx,ny,mz,0)+b*f_[j].cmplx(mx,ny,mz,0));
                            Ryk_.add(ny, a*u_[j].cmplx(mx,ny,mz,1)+b*f_[j].cmplx(mx,ny,mz,1));
                            Rzk_.add(ny, a*u_[j].cmplx(mx,ny,mz,2)+b*f_[j].cmplx(mx,ny,mz,2));
                        }
                    }
                    
                    // Solve the tau solutions
                    if (kx!=0 || kz!=0)
                        tausolver_[mx][mz].solve_Inhom(uk_,vk_,wk_,Pk_, Rxk_,Ryk_,Rzk_, vr_lower, vr_upper, vi_lower, vi_upper, u_lower, u_upper, w_lower, w_upper);
                    
                    else { // kx,kz == 0,0
                        
                        if (Ubaseyy_.length() > 0)
                            for (int ny=0; ny<Ny_; ++ny)
                                Rxk_.re[ny] += nu_*Ubaseyy_[ny];   // Rx has addl'l term from Ubase
                        
                        if (flags_.constraint == PressureGradient) {
                            // pressure is supplied, put on RHS of tau eqn
                            Rxk_.re[0] -= dPdxRef_;
                            
                            // Solve the tau equations
                            //tausolver_[mx][mz].solve(uk_, vk_, wk_, Pk_, Rxk_,Ryk_,Rzk_);
                            tausolver_[mx][mz].solve_Inhom(uk_,vk_,wk_,Pk_, Rxk_,Ryk_,Rzk_, vr_lower, vr_upper, vi_lower, vi_upper, u_lower, u_upper, w_lower, w_upper);
                            
                            // Bulk vel is free variable determined from soln of tau eqn
                            UbulkAct_ = UbulkBase_ + ubulkBase_ + uk_.re.mean();
                            dPdxAct_ = dPdxRef_;
                        }
                        else { // const bulk velocity
                            // bulk velocity is supplied, put on RHS of tau eqn
                            // dPdxAct (i.e. at next time step is solved for)
                            // constraint:    UbulkBase + ubulkBase + mean(u) = UbulkRef.
                            tausolver_[mx][mz].solve_Inhom(uk_, vk_, wk_, Pk_, dPdxAct_,
                                                           Rxk_, Ryk_, Rzk_, vr_lower, vr_upper, vi_lower, vi_upper,                                                   UbulkRef_ - UbulkBase_ - ubulkBase_);
                            UbulkAct_ = UbulkBase_ + ubulkBase_ + uk_.re.mean(); // should == UbulkRef_
                            //UbulkAct_ = UbulkRef_;
                        }
                    }
                    
                    // Load solutions into u_[J], the last element of u history, which is not needed anymore.
                    // Because of FFTW complex symmetries
                    // The 0,0 mode must be real.
                    // For Nx even, the kxmax,0 mode must be real
                    // For Nz even, the 0,kzmax mode must be real
                    // For Nx,Nz even, the kxmax,kzmax mode must be real
                    if ((kx == 0 && kz == 0) ||
                        (Nx_%2 == 0 && kx == kxmax && kz == 0) ||
                        (Nz_%2 == 0 && kz == kzmax && kx == 0) ||
                        (Nx_%2 == 0 && Nz_%2 == 0 && kx == kxmax && kz == kzmax)) {
                        
                        for (int ny=0; ny<Ny_; ++ny) {
                            u_[J].cmplx(mx,ny,mz,0) = Complex(Re(uk_[ny]), 0.0);
                            u_[J].cmplx(mx,ny,mz,1) = Complex(Re(vk_[ny]), 0.0);
                            u_[J].cmplx(mx,ny,mz,2) = Complex(Re(wk_[ny]), 0.0);
                            qn.cmplx(mx,ny,mz,0) = Complex(Re(Pk_[ny]), 0.0);
                        }
                    }
                    // The normal case, for general kx,kz
                    else
                        for (int ny=0; ny<Ny_; ++ny) {
                            u_[J].cmplx(mx,ny,mz,0) = uk_[ny];
                            u_[J].cmplx(mx,ny,mz,1) = vk_[ny];
                            u_[J].cmplx(mx,ny,mz,2) = wk_[ny];
                            qn.cmplx(mx,ny,mz,0) = Pk_[ny];
                        }
                    
                    // And now set the y-aliased modes to zero.
                    for (int ny=Nyd_; ny<Ny_; ++ny) {
                        u_[J].cmplx(mx,ny,mz,0) = 0.0;
                        u_[J].cmplx(mx,ny,mz,1) = 0.0;
                        u_[J].cmplx(mx,ny,mz,2) = 0.0;
                        qn.cmplx(mx,ny,mz,0) = 0.0;
                    }
                }
            }
            
            
            
            // The solution is stored in u_[J]. Shift entire u and f arrays in time
            // Ie shift u_[J] <- u_[J-1] <- ... <- u_[0] <- u_[J]
            for (int j=order_-1; j>0; --j) {
                swap(f_[j], f_[j-1]);
                swap(u_[j], u_[j-1]);
            }
            t_ += dt_;
            
            //printStack();
            //*flags_.logstream << "} Multistep::advance(...) step == " << step << " }" <<endl;
            
            if (flags_.verbosity == PrintTime || flags_.verbosity == PrintAll)
                *flags_.logstream << t_ << ' ' << flush;
            else if (flags_.verbosity == PrintTicks)
                *flags_.logstream << '.' << flush;
            
            
            
           

            //un = u_[0]; 
            qn.setPadded(flags_.dealias_xz());
            
            cfl_ = u_[0].CFLfactor(Ubase_);
            cfl_ *= flags_.dealias_xz() ? 2.0*pi/3.0*dt_ : pi*dt_;
            
            if (cfl_ > 1.0){
                cout << "Error: CFL>1" << endl;
                break;
            }

            
        }
        
        un = u_[0]; 
        
        qn.setPadded(flags_.dealias_xz());
        
        cfl_ = u_[0].CFLfactor(Ubase_);
        cfl_ *= flags_.dealias_xz() ? 2.0*pi/3.0*dt_ : pi*dt_;
        
        if (flags_.verbosity == PrintTime || flags_.verbosity == PrintAll ||
            flags_.verbosity == PrintTicks)
            *flags_.logstream << endl;
        
        return;
    }

/*
  void MultistepDNS::advance_inhom_CON_SF(Controller& controller, FlowField& un, FlowField& qn, FlowField& BCs, FlowField& Fn, int Nsteps, double*** IO) {
        
        const int kxmax = un.kxmax();
        const int kzmax = un.kzmax();
        const int J = order_ -1 ;
        
        u_[0] = un;
       
                 
        
        for (int step=0; step<Nsteps; ++step) {
            
            

            controller.advance_Con_SF(u_[0],qn,BCs,dt_,t_+dt_,IO);
            //cout << "dt=" << dt_ <<  " " << "t=" << t_ << " " << "BC" << BCs.cmplx(0,0,1,1) << endl;

                        
            
            //*flags_.logstream << "Multistep::advance(...) step == " << step << " {" <<endl;
            //printStack();
            
            // Calculate nonlinearity
            if (order_ > 0)
                navierstokesNL(u_[0], ubase_, Ubase_, f_[0], tmp_, tmp2_, flags_.nonlinearity);
            
	     Fn = f_[0];

            // Update each Fourier mode with time-stepping algorithm
            for (int mx=0; mx<Mx_; ++mx) {
                const int kx = un.kx(mx);
                
                for (int mz=0; mz<Mz_; ++mz) {
                    const int kz = un.kz(mz);
                    //*************************
                                       
                    
                    Complex u_upper=BCs.cmplx(mx,1,mz,0);
                    Complex u_lower=BCs.cmplx(mx,0,mz,0);
                    Complex v_upper=BCs.cmplx(mx,1,mz,1);
                    Complex v_lower=BCs.cmplx(mx,0,mz,1);
                    Complex w_upper=BCs.cmplx(mx,1,mz,2);
                    Complex w_lower=BCs.cmplx(mx,0,mz,2);
                    
                    Real vr_upper=real(v_upper);
                    Real vi_upper=imag(v_upper);
                    Real vr_lower=real(v_lower);
                    Real vi_lower=imag(v_lower);
                    
                    const Real c = 4.0*square(pi)*nu_;
                    
                    Real lambda = eta_/dt_ + c*(square(kx/Lx_) + square(kz/Lz_));
                    
                    tausolver_[mx][mz] = TauSolver(kx,kz,Lx_,Lz_,a_,b_,lambda,nu_,Ny_, vr_lower, vi_lower, vr_upper, vi_upper, flags_.taucorrection);
                    
                    // ***************************
                    
                    // Zero out the aliased modes and break to next kx,kz
                    if ((kx == kxmax || kz == kzmax) ||
                        (flags_.dealias_xz() && isAliasedMode(kx,kz))) {
                        for (int ny=0; ny<Ny_; ++ny) {
                            u_[0].cmplx(mx,ny,mz,0) = 0.0;
                            u_[0].cmplx(mx,ny,mz,1) = 0.0;
                            u_[0].cmplx(mx,ny,mz,2) = 0.0;
                            qn.cmplx(mx,ny,mz,0) = 0.0;
                        }
                        break;
                    }
                    
                    
                     // Zero out y-aliased modes for all kx,kz
                     if (flags.dealias_y() || flags.dealias_xyz() {
                     for (int ny=iround(2*Nyd_)/3; ny<Nyd_; ++ny) {
                     u_[0].cmplx(mx,ny,mz,0) = 0.0;
                     u_[0].cmplx(mx,ny,mz,1) = 0.0;
                     u_[0].cmplx(mx,ny,mz,2) = 0.0;
                     qn.cmplx(mx,ny,mz,0) = 0.0;
                     }                    
                     }
                     
                    
                    // For nonaliased modes
                    Rxk_.setToZero();
                    Ryk_.setToZero();
                    Rzk_.setToZero();
                    
                    // Add up multistepping terms of linear and nonlinear terms
                    for (int j=0; j<order_; ++j) {
                        const Real a = -alpha_[j]/dt_;
                        const Real b = -beta_[j];
                        for (int ny=0; ny<Ny_; ++ny) {
                            Rxk_.add(ny, a*u_[j].cmplx(mx,ny,mz,0)+b*f_[j].cmplx(mx,ny,mz,0));
                            Ryk_.add(ny, a*u_[j].cmplx(mx,ny,mz,1)+b*f_[j].cmplx(mx,ny,mz,1));
                            Rzk_.add(ny, a*u_[j].cmplx(mx,ny,mz,2)+b*f_[j].cmplx(mx,ny,mz,2));
                        }
                    }
                    
                    // Solve the tau solutions
                    if (kx!=0 || kz!=0)
                        tausolver_[mx][mz].solve_Inhom(uk_,vk_,wk_,Pk_, Rxk_,Ryk_,Rzk_, vr_lower, vr_upper, vi_lower, vi_upper, u_lower, u_upper, w_lower, w_upper);
                    
                    else { // kx,kz == 0,0
                        
                        if (Ubaseyy_.length() > 0)
                            for (int ny=0; ny<Ny_; ++ny)
                                Rxk_.re[ny] += nu_*Ubaseyy_[ny];   // Rx has addl'l term from Ubase
                        
                        if (flags_.constraint == PressureGradient) {
                            // pressure is supplied, put on RHS of tau eqn
                            Rxk_.re[0] -= dPdxRef_;
                            
                            // Solve the tau equations
                            //tausolver_[mx][mz].solve(uk_, vk_, wk_, Pk_, Rxk_,Ryk_,Rzk_);
                            tausolver_[mx][mz].solve_Inhom(uk_,vk_,wk_,Pk_, Rxk_,Ryk_,Rzk_, vr_lower, vr_upper, vi_lower, vi_upper, u_lower, u_upper, w_lower, w_upper);
                            
                            // Bulk vel is free variable determined from soln of tau eqn
                            UbulkAct_ = UbulkBase_ + ubulkBase_ + uk_.re.mean();
                            dPdxAct_ = dPdxRef_;
                        }
                        else { // const bulk velocity
                            // bulk velocity is supplied, put on RHS of tau eqn
                            // dPdxAct (i.e. at next time step is solved for)
                            // constraint:    UbulkBase + ubulkBase + mean(u) = UbulkRef.
                            tausolver_[mx][mz].solve_Inhom(uk_, vk_, wk_, Pk_, dPdxAct_,
                                                           Rxk_, Ryk_, Rzk_, vr_lower, vr_upper, vi_lower, vi_upper,                                                   UbulkRef_ - UbulkBase_ - ubulkBase_);
                            UbulkAct_ = UbulkBase_ + ubulkBase_ + uk_.re.mean(); // should == UbulkRef_
                            //UbulkAct_ = UbulkRef_;
                        }
                    }
                    
                    // Load solutions into u_[J], the last element of u history, which is not needed anymore.
                    // Because of FFTW complex symmetries
                    // The 0,0 mode must be real.
                    // For Nx even, the kxmax,0 mode must be real
                    // For Nz even, the 0,kzmax mode must be real
                    // For Nx,Nz even, the kxmax,kzmax mode must be real
                    if ((kx == 0 && kz == 0) ||
                        (Nx_%2 == 0 && kx == kxmax && kz == 0) ||
                        (Nz_%2 == 0 && kz == kzmax && kx == 0) ||
                        (Nx_%2 == 0 && Nz_%2 == 0 && kx == kxmax && kz == kzmax)) {
                        
                        for (int ny=0; ny<Ny_; ++ny) {
                            u_[J].cmplx(mx,ny,mz,0) = Complex(Re(uk_[ny]), 0.0);
                            u_[J].cmplx(mx,ny,mz,1) = Complex(Re(vk_[ny]), 0.0);
                            u_[J].cmplx(mx,ny,mz,2) = Complex(Re(wk_[ny]), 0.0);
                            qn.cmplx(mx,ny,mz,0) = Complex(Re(Pk_[ny]), 0.0);
                        }
                    }
                    // The normal case, for general kx,kz
                    else
                        for (int ny=0; ny<Ny_; ++ny) {
                            u_[J].cmplx(mx,ny,mz,0) = uk_[ny];
                            u_[J].cmplx(mx,ny,mz,1) = vk_[ny];
                            u_[J].cmplx(mx,ny,mz,2) = wk_[ny];
                            qn.cmplx(mx,ny,mz,0) = Pk_[ny];
                        }
                    
                    // And now set the y-aliased modes to zero.
                    for (int ny=Nyd_; ny<Ny_; ++ny) {
                        u_[J].cmplx(mx,ny,mz,0) = 0.0;
                        u_[J].cmplx(mx,ny,mz,1) = 0.0;
                        u_[J].cmplx(mx,ny,mz,2) = 0.0;
                        qn.cmplx(mx,ny,mz,0) = 0.0;
                    }
                }
            }
            
            
            
            // The solution is stored in u_[J]. Shift entire u and f arrays in time
            // Ie shift u_[J] <- u_[J-1] <- ... <- u_[0] <- u_[J]
            for (int j=order_-1; j>0; --j) {
                swap(f_[j], f_[j-1]);
                swap(u_[j], u_[j-1]);
            }
            t_ += dt_;
            
            //printStack();
            //*flags_.logstream << "} Multistep::advance(...) step == " << step << " }" <<endl;
            
            if (flags_.verbosity == PrintTime || flags_.verbosity == PrintAll)
                *flags_.logstream << t_ << ' ' << flush;
            else if (flags_.verbosity == PrintTicks)
                *flags_.logstream << '.' << flush;
            
            
            
           

            //un = u_[0]; 
            qn.setPadded(flags_.dealias_xz());
            
            cfl_ = u_[0].CFLfactor(Ubase_);
            cfl_ *= flags_.dealias_xz() ? 2.0*pi/3.0*dt_ : pi*dt_;
            
            if (cfl_ > 1.0){
                cout << "Error: CFL>1" << endl;
                break;
            }

            
        }
        
        un = u_[0]; 
        
        qn.setPadded(flags_.dealias_xz());
        
        cfl_ = u_[0].CFLfactor(Ubase_);
        cfl_ *= flags_.dealias_xz() ? 2.0*pi/3.0*dt_ : pi*dt_;
        
        if (flags_.verbosity == PrintTime || flags_.verbosity == PrintAll ||
            flags_.verbosity == PrintTicks)
            *flags_.logstream << endl;
        
        return;
    }

*/
    

//////////////////////////////////////////////////////////////////////////////////
    
    
    
    /**************************************************************************
     
     ---End: Peter H Heins---
     
     **************************************************************************/   
    
    
    
    
    
    
    // u1 = P(s,u0) :
    // tmp = u0;
    // tmp *= s
    // tmp += u0
    // tmp *= 0.5
    // u1 = tmp
    
    void MultistepDNS::project() {
        for (int n=0; n<u_.length(); ++n)
            u_[n].project(flags_.symmetries);
        for (int n=0; n<f_.length(); ++n)
            f_[n].project(flags_.symmetries);
    }
    
    void MultistepDNS::operator*=(const FieldSymmetry& sigma) {
        for (int n=0; n<u_.length(); ++n)
            u_[n] *= sigma;
        for (int n=0; n<f_.length(); ++n)
            f_[n] *= sigma;
    }
    
    void MultistepDNS::printStack() const {
        *flags_.logstream << "Multistep::printStack() {" << endl;
        *flags_.logstream << "        t == " << t_ << endl;
        *flags_.logstream << "countdown == " << countdown_ << endl;
        *flags_.logstream << "     full == " << full() << endl;
        
        for (int j=order_-1; j>=0; --j)
            printf("j=%2d t=%5.2f L2(uj)=%13.10f L2(fj)=%13.10f\n",
                   j, t_-j*dt_, L2Norm(u_[j]), L2Norm(f_[j]));
        *flags_.logstream << endl;
        *flags_.logstream << "}" << endl;
    }
    
    bool MultistepDNS::push(const FlowField& un) {
        
        //*flags_.logstream << "MultistepDNS::push(const FlowField& un) { " << endl;
        //printStack();
        // Let K = order-1. Arrays are then u_[0:K], f_[0:K]
        // Shift u_[K] <- u_[K-1] <- ... <- u_[0] <- un
        for (int j=order_-1; j>0; --j) {
            swap(u_[j], u_[j-1]);
            swap(f_[j], f_[j-1]);
        }
        
        if (order_ > 0) {
            u_[0] = un;
            navierstokesNL(u_[0], ubase_, Ubase_, f_[0], tmp_, tmp2_, flags_.nonlinearity);
            
            cfl_ = u_[0].CFLfactor(Ubase_);
            cfl_ *= flags_.dealias_xz() ? 2.0*pi/3.0*dt_ : pi*dt_;
        }
        
        t_ += dt_;
        --countdown_;
        //printStack();
        //*flags_.logstream << "}" << endl;
        return full();
    }
    
    bool MultistepDNS::full() const {
        return (countdown_ == 0) ? true : false;
    }
    
    /**************************************
     void MultistepDNS::reset_uj(const FlowField& uj, int j) {
     if (j<0 || j>=order_) {
     cerr << "error in MultistepDNS::reset_uj(uj, j) : j = " << j
	 << " is out of bounds" << endl;
     exit(1);
     }
     u_[j] = uj;
     navierstokesNL(u_[j], Ubase_, f_[j], tmp_, flags_.nonlinearity);
     }
     **********************************************/
    // ==============================================================
    // Runge-Kutta algorithms
    
    RungeKuttaDNS::RungeKuttaDNS()
    :
    DNSAlgorithm(),
    Qj1_(),
    Qj_()
    {}
    
    RungeKuttaDNS::RungeKuttaDNS(const RungeKuttaDNS& dns)
    :
    DNSAlgorithm(dns),
    Nsubsteps_(dns.Nsubsteps_),
    Qj1_(dns.Qj1_),
    Qj_(dns.Qj_),
    A_(dns.A_),
    B_(dns.B_),
    C_(dns.C_)
    {
        // Allocate memory for [Nsubsteps x Mx_ x Mz_] Tausolver arrays
        // and copy tausolvers from dns argument
        tausolver_ = new TauSolver**[Nsubsteps_];    // new #1
        for (int j=0; j<Nsubsteps_; ++j) {
            tausolver_[j] = new TauSolver*[Mx_];       // new #2
            for (int mx=0; mx<Mx_; ++mx) {
                tausolver_[j][mx] = new TauSolver[Mz_];  // new #3
                for (int mz=0; mz<Mz_; ++mz)
                    tausolver_[j][mx][mz] = dns.tausolver_[j][mx][mz];
                
            }
        }
    }
    
    // This algorithm is described in "Spectral Methods for Incompressible Viscous
    // Flow", Roger Peyret, Springer-Verlag Appl Math Sci series vol 148, 2002.
    // section 4.5.2.c.2 "Three-stage scheme (RK3/CN)". I use C for his B'
    
    /*
    RungeKuttaDNS::RungeKuttaDNS(const FlowField& u, const ChebyCoeff& Ubase,
                                 Real nu, Real dt, const DNSFlags& flags, Real t)
    :
    DNSAlgorithm(u,Ubase,nu,dt,flags,t),
    Qj1_(u),
    Qj_(u)
    {
        Qj1_.setToZero();
        Qj_.setToZero();
        
        TimeStepMethod algorithm = flags.timestepping;
        switch (algorithm) {
            case CNRK2:
                order_ = 2;
                Nsubsteps_ = 3;
                Ninitsteps_ = 0;
                A_.resize(Nsubsteps_);
                B_.resize(Nsubsteps_);
                C_.resize(Nsubsteps_);
                A_[0] = 0.0;     A_[1] = -5.0/9.0;  A_[2] = -153.0/128.0; // Peyret A
                B_[0] = 1.0/3.0; B_[1] = 15.0/16.0; B_[2] = 8.0/15.0;     // Peyret B
                C_[0] = 1.0/6.0; C_[1] = 5.0/24.0;  C_[2] = 1.0/8.0;      // Peyret B'
                break;
            default:
                cerr << "RungeKuttaDNS::RungeKuttaDNS(un,Ubase,nu,dt,flags,t0)\n"
                << "error: flags.timestepping == " << algorithm
                << " is a non-runge-kutta algorithm" << endl;
                exit(1);
        }
        
        // Allocate memory for [Nsubsteps x Mx_ x Mz_] Tausolver array
        tausolver_ = new TauSolver**[Nsubsteps_];           // new #1
        for (int j=0; j<Nsubsteps_; ++j) {
            tausolver_[j] = new TauSolver*[Mx_];       // new #2
            for (int mx=0; mx<Mx_; ++mx)
                tausolver_[j][mx] = new TauSolver[Mz_];  // new #3
        }
        
        
        // Could create BC flowfield of zeros here for arg in reset_dtIH(dt,BC)
        // and have a boolean variable for controlled/uncontrolled
        
        cout << "Reset_RK_dt:";
        reset_dt(dt_);
        cout << "Done" << endl;
    }
    */
    
    ////////////////////////////////////////////////////////////
    // Added by P. Heins 24/10
    RungeKuttaDNS::RungeKuttaDNS(const FlowField& u, const ChebyCoeff& Ubase,
                                 Real nu, Real dt, const DNSFlags& flags, Real t, bool controlled)
    :
    DNSAlgorithm(u,Ubase,nu,dt,flags,t,controlled),
    Qj1_(u),
    Qj_(u)
    {
        Qj1_.setToZero();
        Qj_.setToZero();
        
        TimeStepMethod algorithm = flags.timestepping;
        switch (algorithm) {
            case CNRK2:
                order_ = 2;
                Nsubsteps_ = 3;
                Ninitsteps_ = 0;
                A_.resize(Nsubsteps_);
                B_.resize(Nsubsteps_);
                C_.resize(Nsubsteps_);
                A_[0] = 0.0;     A_[1] = -5.0/9.0;  A_[2] = -153.0/128.0; // Peyret A
                B_[0] = 1.0/3.0; B_[1] = 15.0/16.0; B_[2] = 8.0/15.0;     // Peyret B
                C_[0] = 1.0/6.0; C_[1] = 5.0/24.0;  C_[2] = 1.0/8.0;      // Peyret B'
                break;
            default:
                cerr << "RungeKuttaDNS::RungeKuttaDNS(un,Ubase,nu,dt,flags,t0)\n"
                << "error: flags.timestepping == " << algorithm
                << " is a non-runge-kutta algorithm" << endl;
                exit(1);
        }
        
        // Allocate memory for [Nsubsteps x Mx_ x Mz_] Tausolver array
        tausolver_ = new TauSolver**[Nsubsteps_];           // new #1
        for (int j=0; j<Nsubsteps_; ++j) {
            tausolver_[j] = new TauSolver*[Mx_];       // new #2
            for (int mx=0; mx<Mx_; ++mx)
                tausolver_[j][mx] = new TauSolver[Mz_];  // new #3
        }
        
        
        // Could create BC flowfield of zeros here for arg in reset_dtIH(dt,BC)
        // and have a boolean variable for controlled/uncontrolled
        
        
        if (controlled){
            FlowField BC(Nx_,2,Nz_,3,Lx_,Lz_,a_,b_);
            BC.makeState(Spectral,Physical);
            reset_dtIH(dt_,BC);  
        }
        else {
            reset_dt(dt_);
        }
    }

    ////////////////////////////////////////////////////////////
    
    
    
    
    
    
    RungeKuttaDNS::~RungeKuttaDNS() {
        if (tausolver_) {
            for (int j=0; j<Nsubsteps_; ++j) {
                for (int mx=0; mx<Mx_; ++mx)
                    delete[] tausolver_[j][mx];  // undo new #3
                delete[] tausolver_[j];        // undo new #2
            }
            delete[] tausolver_;                   // undo new #1
        }
        tausolver_ = 0;
    }
    
    DNSAlgorithm* RungeKuttaDNS::clone() const {
        return new RungeKuttaDNS(*this);
    }
    
    void RungeKuttaDNS::reset_dt(Real dt) {
        
        cfl_ *= dt/dt_;
        //nu_ = nu;
        dt_ = dt;
        const Real c = 4.0*square(pi)*nu_;
        const int kxmax = tmp_.kxmax();
        const int kzmax = tmp_.kzmax();
        
        // This loop replaces the TauSolver objects at tausolver_[i][mx][mz]
        // with new TauSolver objects configured with appropriate parameters
        for (int j=0; j<Nsubsteps_; ++j) {
            for (int mx=0; mx<Mx_; ++mx) {
                int kx = tmp_.kx(mx);
                for (int mz=0; mz<Mz_; ++mz) {
                    int kz = tmp_.kz(mz);
                    Real lambda = 1.0/(C_[j]*dt_) + c*(square(kx/Lx_)+square(kz/Lz_));
                    
                    if ((kx != kxmax || kz != kzmax) &&
                        (!flags_.dealias_xz() || !isAliasedMode(kx,kz)))
                        
                        tausolver_[j][mx][mz] =
                        TauSolver(kx, kz, Lx_, Lz_, a_, b_, lambda, nu_, Ny_,
                                  flags_.taucorrection);
                }
            }
        }
    }
    
    //****************************************************************************
    //****************************************************************************
    void RungeKuttaDNS::reset_dtIH(Real dt, FlowField& BCs) {
        

        cfl_ *= dt/dt_;
        //nu_ = nu;
        dt_ = dt;
        const Real c = 4.0*square(pi)*nu_;
        const int kxmax = tmp_.kxmax();
        const int kzmax = tmp_.kzmax();
        
        
        // This loop replaces the TauSolver objects at tausolver_[i][mx][mz]
        // with new TauSolver objects configured with appropriate parameters
        for (int j=0; j<Nsubsteps_; ++j) {
            for (int mx=0; mx<Mx_; ++mx) {
                int kx = tmp_.kx(mx);
                for (int mz=0; mz<Mz_; ++mz) {
                    int kz = tmp_.kz(mz);
                    
                    
                    Complex v_upper=BCs.cmplx(mx,1,mz,1);
                    Complex v_lower=BCs.cmplx(mx,0,mz,1);
                    
                    double vr_upper=real(v_upper);
                    double vi_upper=imag(v_upper);
                    double vr_lower=real(v_lower);
                    double vi_lower=imag(v_lower);
                    
                    
                    Real lambda = 1.0/(C_[j]*dt_) + c*(square(kx/Lx_)+square(kz/Lz_));
                    
                    if ((kx != kxmax || kz != kzmax) &&
                        (!flags_.dealias_xz() || !isAliasedMode(kx,kz)))
                        
                        //tausolver_[j][mx][mz] = TauSolver(kx, kz, Lx_, Lz_, a_, b_, lambda,nu_, Nyd_, flags_.taucorrection);
                        
                        tausolver_[j][mx][mz] = TauSolver(kx,kz,Lx_,Lz_,a_,b_,lambda,nu_,Ny_, vr_lower, vi_lower, vr_upper, vi_upper, flags_.taucorrection);
                }
            }
        }
    }
    //****************************************************************************
    //****************************************************************************
    
  void RungeKuttaDNS::advance(FlowField& un, FlowField& qn, int Nsteps) {
        
        const int kxmax = un.kxmax();
        const int kzmax = un.kzmax();
        
        for (int n=0; n<Nsteps; ++n) {
            for (int j=0; j<Nsubsteps_; ++j) {
                
                FlowField& uj(un); // Store uj in un during substeps, reflect in notation
                
                // Efficient implementation of
                // Q_{j+1} = A_j Q_j + N(u_j)}  where N = -u grad u
                // Q_{j+1} = A_j Q_j - f(u_j)}  where f =  u grad u
                // Q_j = Q_{j+1}
                Qj_ *= A_[j];
                navierstokesNL(uj, ubase_, Ubase_, Qj1_, tmp_, tmp2_, flags_.nonlinearity);
                Qj_ -= Qj1_; // subtract because navierstokesNL(u) = u grad u = -N(u)
                
               

                // Update each Fourier mode with time-stepping algorithm
                for (int mx=0; mx<Mx_; ++mx) {
                    const int kx = un.kx(mx);
                    
                    for (int mz=0; mz<Mz_; ++mz) {
                        const int kz = un.kz(mz);
                        
                        // Zero out the aliased modes
                        if ((kx == kxmax || kz == kzmax) ||
                            (flags_.dealias_xz() && isAliasedMode(kx,kz))) {
                            for (int ny=0; ny<Ny_; ++ny) {
                                un.cmplx(mx,ny,mz,0) = 0.0;
                                un.cmplx(mx,ny,mz,1) = 0.0;
                                un.cmplx(mx,ny,mz,2) = 0.0;
                                qn.cmplx(mx,ny,mz,0) = 0.0;
                            }
                            break;
                        }
                        
                        Rxk_.setToZero();
                        Ryk_.setToZero();
                        Rzk_.setToZero();
                        
                        // Make the following assignments in prep for computation of RHS
                        
                        // nu (uk,vk,wk) = nu ujn(0,1,2)
                        //            Pk = qn
                        
                        // Goal is to compute
                        // R = nu uj" + [1/(Cj dt) - nu kappa2]    uj  - grad qj + Bj/Cj Qj
                        //   = nu uj" + [1/(nu Cj dt) - kappa2] nu uj  - grad qj + Bj/Cj Qj
                        
                        // Extract relevant Fourier modes of uj and qj for computations
                        for (int ny=0; ny<Ny_; ++ny) {
                            uk_.set(ny, nu_*uj.cmplx(mx,ny,mz,0));
                            vk_.set(ny, nu_*uj.cmplx(mx,ny,mz,1));
                            wk_.set(ny, nu_*uj.cmplx(mx,ny,mz,2));
                            Pk_.set(ny, qn.cmplx(mx,ny,mz,0));
                        }
                        
                        // (1) Put nu uj" into in R. (Pyk_ is used as tmp workspace)
                        diff2(uk_, Rxk_, Pyk_);
                        diff2(vk_, Ryk_, Pyk_);
                        diff2(wk_, Rzk_, Pyk_);
                        
                        // (2) Put qn' into Pyk (compute y-comp of pressure gradient).
                        diff(Pk_, Pyk_);
                        
                        // (3) Add [1/(nu Cj dt)- kappa2] nu uj - grad qj + Bj/Cj Qj to R.
                        const Real c =
                        1.0/(nu_*C_[j]*dt_) - 4*pi*pi*(square(kx/Lx_) + square(kz/Lz_));
                        const Real B_C = B_[j]/C_[j];
                        const Complex Dx = un.Dx(mx);
                        const Complex Dz = un.Dz(mz);
                        
                        for (int ny=0; ny<Ny_; ++ny) {
                            Rxk_.add(ny, c*uk_[ny] + B_C*Qj_.cmplx(mx,ny,mz,0) - Dx*Pk_[ny]);
                            Ryk_.add(ny, c*vk_[ny] + B_C*Qj_.cmplx(mx,ny,mz,1) - Pyk_[ny]);
                            Rzk_.add(ny, c*wk_[ny] + B_C*Qj_.cmplx(mx,ny,mz,2) - Dz*Pk_[ny]);
                        }
                        
                        // Do the tau solutions
                        if (kx!=0 || kz!=0)
                            tausolver_[j][mx][mz].solve(uk_,vk_,wk_,Pk_, Rxk_,Ryk_,Rzk_);
                        
                        else { // kx,kz == 0,0
                            // Rx has additional terms, nu Uyy at both t=j and t=j+1
                            const Real c = 2*nu_;
                            if (Ubaseyy_.length() > 0)
                                for (int ny=0; ny<Ny_; ++ny)
                                    Rxk_.re[ny] += c*Ubaseyy_[ny];
                            
                            if (flags_.constraint == PressureGradient) {
                                // dPdx is supplied, put dPdx at both t=j and t=j+1 on RHS
                                Rxk_.re[0] -= dPdxAct_ + dPdxRef_;
                                
                                // Solve the tau equations
                                tausolver_[j][mx][mz].solve(uk_, vk_, wk_, Pk_, Rxk_, Ryk_, Rzk_);
                                
                                // Bulk velocity is free variable on LHS solved by tau eqn
                                UbulkAct_ = UbulkBase_ + ubulkBase_ + uk_.re.mean();
                                dPdxAct_ = dPdxRef_;
                            }
                            else { // const bulk velocity
                                // Add the previous time-step's -dPdx to the RHS. The next
                                // timestep's dPdx term appears on LHS as unknown.
                                Rxk_.re[0] -= dPdxAct_;
                                
                                // Use tausolver with additional variable and constraint:
                                // free variable: dPdxAct at next time-step,
                                // constraint:    UbulkBase + ubulkBase + mean(u) = UbulkRef.
                                tausolver_[j][mx][mz].solve(uk_, vk_, wk_, Pk_, dPdxAct_,
                                                            Rxk_, Ryk_, Rzk_,
                                                            UbulkRef_ - UbulkBase_ - ubulkBase_);
                                UbulkAct_ = UbulkBase_ + ubulkBase_ + uk_.re.mean(); // should == UbulkRef_
                                //UbulkAct_ = UbulkRef_;
                            }
                            // for kx=kz=0, constant term of pressure is arbitrary 3/19/05
                            // Pk_.set(0, Complex(0.0, 0.0));
                        }
                        
                        // Load solutions back into the external 3d data arrays.
                        // Because of FFTW complex symmetries
                        // The 0,0 mode must be real.
                        // For Nx even, the kxmax,0 mode must be real
                        // For Nz even, the 0,kzmax mode must be real
                        // For Nx,Nz even, the kxmax,kzmax mode must be real
                        if ((kx == 0 && kz == 0) ||
                            (Nx_%2 == 0 && kx == kxmax && kz == 0) ||
                            (Nz_%2 == 0 && kz == kzmax && kx == 0) ||
                            (Nx_%2 == 0 && Nz_%2 == 0 && kx == kxmax && kz == kzmax)) {
                            
                            for (int ny=0; ny<Ny_; ++ny) {
                                un.cmplx(mx,ny,mz,0) = Complex(Re(uk_[ny]), 0.0);
                                un.cmplx(mx,ny,mz,1) = Complex(Re(vk_[ny]), 0.0);
                                un.cmplx(mx,ny,mz,2) = Complex(Re(wk_[ny]), 0.0);
                                qn.cmplx(mx,ny,mz,0) = Complex(Re(Pk_[ny]), 0.0);
                            }
                        }
                        // The normal case, for general kx,kz
                        else
                            for (int ny=0; ny<Ny_; ++ny) {
                                un.cmplx(mx,ny,mz,0) = uk_[ny];
                                un.cmplx(mx,ny,mz,1) = vk_[ny];
                                un.cmplx(mx,ny,mz,2) = wk_[ny];
                                qn.cmplx(mx,ny,mz,0) = Pk_[ny];
                            }
                        
                        // And now set the y-aliased modes to zero.
                        for (int ny=Nyd_; ny<Ny_; ++ny) {
                            un.cmplx(mx,ny,mz,0) = 0.0;
                            un.cmplx(mx,ny,mz,1) = 0.0;
                            un.cmplx(mx,ny,mz,2) = 0.0;
                            qn.cmplx(mx,ny,mz,0) = 0.0;
                        }
                    }
                }
            }
            t_ += dt_;
            
            if (flags_.verbosity == PrintTime || flags_.verbosity == PrintAll)
                *flags_.logstream << t_ << ' ' << flush;
            else if (flags_.verbosity == PrintTicks)
                *flags_.logstream << '.' << flush;
        }
        
        cfl_ = un.CFLfactor(Ubase_);
        cfl_ *= flags_.dealias_xz() ? 2.0*pi/3.0*dt_ : pi*dt_;
        
        // If using dealiasing, set flag in FlowField that compactifies binary IO
        un.setPadded(flags_.dealias_xz());
        qn.setPadded(flags_.dealias_xz());
        
        if (flags_.verbosity == PrintTime || flags_.verbosity == PrintAll)
            *flags_.logstream << endl;
        
        return;
    }
    
    //void RungeKuttaDNS::project() {}
    //void RungeKuttaDNS::operator*= {}
    
    void RungeKuttaDNS::advance_NL(FlowField& un, FlowField& qn, FlowField& Fn, int Nsteps) {
        
        const int kxmax = un.kxmax();
        const int kzmax = un.kzmax();
        
        for (int n=0; n<Nsteps; ++n) {
            for (int j=0; j<Nsubsteps_; ++j) {
                
                FlowField& uj(un); // Store uj in un during substeps, reflect in notation
                
                // Efficient implementation of
                // Q_{j+1} = A_j Q_j + N(u_j)}  where N = -u grad u
                // Q_{j+1} = A_j Q_j - f(u_j)}  where f =  u grad u
                // Q_j = Q_{j+1}
                Qj_ *= A_[j];
                navierstokesNL(uj, ubase_, Ubase_, Qj1_, tmp_, tmp2_, flags_.nonlinearity);
                Qj_ -= Qj1_; // subtract because navierstokesNL(u) = u grad u = -N(u)
                
                Fn = Qj1_;

                // Update each Fourier mode with time-stepping algorithm
                for (int mx=0; mx<Mx_; ++mx) {
                    const int kx = un.kx(mx);
                    
                    for (int mz=0; mz<Mz_; ++mz) {
                        const int kz = un.kz(mz);
                        
                        // Zero out the aliased modes
                        if ((kx == kxmax || kz == kzmax) ||
                            (flags_.dealias_xz() && isAliasedMode(kx,kz))) {
                            for (int ny=0; ny<Ny_; ++ny) {
                                un.cmplx(mx,ny,mz,0) = 0.0;
                                un.cmplx(mx,ny,mz,1) = 0.0;
                                un.cmplx(mx,ny,mz,2) = 0.0;
                                qn.cmplx(mx,ny,mz,0) = 0.0;
                            }
                            break;
                        }
                        
                        Rxk_.setToZero();
                        Ryk_.setToZero();
                        Rzk_.setToZero();
                        
                        // Make the following assignments in prep for computation of RHS
                        
                        // nu (uk,vk,wk) = nu ujn(0,1,2)
                        //            Pk = qn
                        
                        // Goal is to compute
                        // R = nu uj" + [1/(Cj dt) - nu kappa2]    uj  - grad qj + Bj/Cj Qj
                        //   = nu uj" + [1/(nu Cj dt) - kappa2] nu uj  - grad qj + Bj/Cj Qj
                        
                        // Extract relevant Fourier modes of uj and qj for computations
                        for (int ny=0; ny<Ny_; ++ny) {
                            uk_.set(ny, nu_*uj.cmplx(mx,ny,mz,0));
                            vk_.set(ny, nu_*uj.cmplx(mx,ny,mz,1));
                            wk_.set(ny, nu_*uj.cmplx(mx,ny,mz,2));
                            Pk_.set(ny, qn.cmplx(mx,ny,mz,0));
                        }
                        
                        // (1) Put nu uj" into in R. (Pyk_ is used as tmp workspace)
                        diff2(uk_, Rxk_, Pyk_);
                        diff2(vk_, Ryk_, Pyk_);
                        diff2(wk_, Rzk_, Pyk_);
                        
                        // (2) Put qn' into Pyk (compute y-comp of pressure gradient).
                        diff(Pk_, Pyk_);
                        
                        // (3) Add [1/(nu Cj dt)- kappa2] nu uj - grad qj + Bj/Cj Qj to R.
                        const Real c =
                        1.0/(nu_*C_[j]*dt_) - 4*pi*pi*(square(kx/Lx_) + square(kz/Lz_));
                        const Real B_C = B_[j]/C_[j];
                        const Complex Dx = un.Dx(mx);
                        const Complex Dz = un.Dz(mz);
                        
                        for (int ny=0; ny<Ny_; ++ny) {
                            Rxk_.add(ny, c*uk_[ny] + B_C*Qj_.cmplx(mx,ny,mz,0) - Dx*Pk_[ny]);
                            Ryk_.add(ny, c*vk_[ny] + B_C*Qj_.cmplx(mx,ny,mz,1) - Pyk_[ny]);
                            Rzk_.add(ny, c*wk_[ny] + B_C*Qj_.cmplx(mx,ny,mz,2) - Dz*Pk_[ny]);
                        }
                        
                        // Do the tau solutions
                        if (kx!=0 || kz!=0)
                            tausolver_[j][mx][mz].solve(uk_,vk_,wk_,Pk_, Rxk_,Ryk_,Rzk_);
                        
                        else { // kx,kz == 0,0
                            // Rx has additional terms, nu Uyy at both t=j and t=j+1
                            const Real c = 2*nu_;
                            if (Ubaseyy_.length() > 0)
                                for (int ny=0; ny<Ny_; ++ny)
                                    Rxk_.re[ny] += c*Ubaseyy_[ny];
                            
                            if (flags_.constraint == PressureGradient) {
                                // dPdx is supplied, put dPdx at both t=j and t=j+1 on RHS
                                Rxk_.re[0] -= dPdxAct_ + dPdxRef_;
                                
                                // Solve the tau equations
                                tausolver_[j][mx][mz].solve(uk_, vk_, wk_, Pk_, Rxk_, Ryk_, Rzk_);
                                
                                // Bulk velocity is free variable on LHS solved by tau eqn
                                UbulkAct_ = UbulkBase_ + ubulkBase_ + uk_.re.mean();
                                dPdxAct_ = dPdxRef_;
                            }
                            else { // const bulk velocity
                                // Add the previous time-step's -dPdx to the RHS. The next
                                // timestep's dPdx term appears on LHS as unknown.
                                Rxk_.re[0] -= dPdxAct_;
                                
                                // Use tausolver with additional variable and constraint:
                                // free variable: dPdxAct at next time-step,
                                // constraint:    UbulkBase + ubulkBase + mean(u) = UbulkRef.
                                tausolver_[j][mx][mz].solve(uk_, vk_, wk_, Pk_, dPdxAct_,
                                                            Rxk_, Ryk_, Rzk_,
                                                            UbulkRef_ - UbulkBase_ - ubulkBase_);
                                UbulkAct_ = UbulkBase_ + ubulkBase_ + uk_.re.mean(); // should == UbulkRef_
                                //UbulkAct_ = UbulkRef_;
                            }
                            // for kx=kz=0, constant term of pressure is arbitrary 3/19/05
                            // Pk_.set(0, Complex(0.0, 0.0));
                        }
                        
                        // Load solutions back into the external 3d data arrays.
                        // Because of FFTW complex symmetries
                        // The 0,0 mode must be real.
                        // For Nx even, the kxmax,0 mode must be real
                        // For Nz even, the 0,kzmax mode must be real
                        // For Nx,Nz even, the kxmax,kzmax mode must be real
                        if ((kx == 0 && kz == 0) ||
                            (Nx_%2 == 0 && kx == kxmax && kz == 0) ||
                            (Nz_%2 == 0 && kz == kzmax && kx == 0) ||
                            (Nx_%2 == 0 && Nz_%2 == 0 && kx == kxmax && kz == kzmax)) {
                            
                            for (int ny=0; ny<Ny_; ++ny) {
                                un.cmplx(mx,ny,mz,0) = Complex(Re(uk_[ny]), 0.0);
                                un.cmplx(mx,ny,mz,1) = Complex(Re(vk_[ny]), 0.0);
                                un.cmplx(mx,ny,mz,2) = Complex(Re(wk_[ny]), 0.0);
                                qn.cmplx(mx,ny,mz,0) = Complex(Re(Pk_[ny]), 0.0);
                            }
                        }
                        // The normal case, for general kx,kz
                        else
                            for (int ny=0; ny<Ny_; ++ny) {
                                un.cmplx(mx,ny,mz,0) = uk_[ny];
                                un.cmplx(mx,ny,mz,1) = vk_[ny];
                                un.cmplx(mx,ny,mz,2) = wk_[ny];
                                qn.cmplx(mx,ny,mz,0) = Pk_[ny];
                            }
                        
                        // And now set the y-aliased modes to zero.
                        for (int ny=Nyd_; ny<Ny_; ++ny) {
                            un.cmplx(mx,ny,mz,0) = 0.0;
                            un.cmplx(mx,ny,mz,1) = 0.0;
                            un.cmplx(mx,ny,mz,2) = 0.0;
                            qn.cmplx(mx,ny,mz,0) = 0.0;
                        }
                    }
                }
            }
            t_ += dt_;
            
            if (flags_.verbosity == PrintTime || flags_.verbosity == PrintAll)
                *flags_.logstream << t_ << ' ' << flush;
            else if (flags_.verbosity == PrintTicks)
                *flags_.logstream << '.' << flush;
        }
        
        cfl_ = un.CFLfactor(Ubase_);
        cfl_ *= flags_.dealias_xz() ? 2.0*pi/3.0*dt_ : pi*dt_;
        
        // If using dealiasing, set flag in FlowField that compactifies binary IO
        un.setPadded(flags_.dealias_xz());
        qn.setPadded(flags_.dealias_xz());
        
        if (flags_.verbosity == PrintTime || flags_.verbosity == PrintAll)
            *flags_.logstream << endl;
        
        return;
    }
    
    //void RungeKuttaDNS::project() {}
    //void RungeKuttaDNS::operator*= {}


    /************************************************************************    
     ***********************************************************************
     
     Written by: Peter H Heins (Postgraduate Research Student) 
     
     Date: January 2012
     
     Institution: University of Sheffield (ACSE Department)
     
     Purpose: This function implements imhomogeneous BCs into the RK algorithm, no changes to original, just need for _inhom name
     
     Based on code written by Binh Lieu, University of Minnesota March 2009
     
     *********************************************************************** 
     *************************************************************************/   
    
    
    void RungeKuttaDNS::advance_inhom(FlowField& un, FlowField& qn, FlowField& BCs, int Nsteps) {
        
        const int kxmax = un.kxmax();
        const int kzmax = un.kzmax();
        
        
        for (int n=0; n<Nsteps; ++n) {
            for (int j=0; j<Nsubsteps_; ++j) {
                
                FlowField& uj(un); // Store uj in un during substeps, reflect in notation
                
                // Efficient implementation of
                // Q_{j+1} = A_j Q_j + N(u_j)}  where N = -u grad u
                // Q_{j+1} = A_j Q_j - f(u_j)}  where f =  u grad u
                // Q_j = Q_{j+1}
                Qj_ *= A_[j];
                navierstokesNL(uj, ubase_, Ubase_, Qj1_, tmp_, tmp2_, flags_.nonlinearity);
                Qj_ -= Qj1_; // subtract because navierstokesNL(u) = u grad u = -N(u)
                
                // Update each Fourier mode with time-stepping algorithm
                for (int mx=0; mx<Mx_; ++mx) {
                    const int kx = un.kx(mx);
                    
                    for (int mz=0; mz<Mz_; ++mz) {
                        const int kz = un.kz(mz);
                        
                        Complex u_upper=BCs.cmplx(mx,1,mz,0);
                        Complex u_lower=BCs.cmplx(mx,0,mz,0);
                        Complex v_upper=BCs.cmplx(mx,1,mz,1);
                        Complex v_lower=BCs.cmplx(mx,0,mz,1);
                        Complex w_upper=BCs.cmplx(mx,1,mz,2);
                        Complex w_lower=BCs.cmplx(mx,0,mz,2);
                        
                        Real vr_upper=real(v_upper);
                        Real vi_upper=imag(v_upper);
                        Real vr_lower=real(v_lower);
                        Real vi_lower=imag(v_lower);                       
                        
                        
                        
                        // Zero out the aliased modes
                        if ((kx == kxmax || kz == kzmax) ||
                            (flags_.dealias_xz() && isAliasedMode(kx,kz))) {
                            for (int ny=0; ny<Ny_; ++ny) {
                                un.cmplx(mx,ny,mz,0) = 0.0;
                                un.cmplx(mx,ny,mz,1) = 0.0;
                                un.cmplx(mx,ny,mz,2) = 0.0;
                                qn.cmplx(mx,ny,mz,0) = 0.0;
                            }
                            break;
                        }
                        
                        Rxk_.setToZero();
                        Ryk_.setToZero();
                        Rzk_.setToZero();
                        
                        // Make the following assignments in prep for computation of RHS
                        
                        // nu (uk,vk,wk) = nu ujn(0,1,2)
                        //            Pk = qn
                        
                        // Goal is to compute
                        // R = nu uj" + [1/(Cj dt) - nu kappa2]    uj  - grad qj + Bj/Cj Qj
                        //   = nu uj" + [1/(nu Cj dt) - kappa2] nu uj  - grad qj + Bj/Cj Qj
                        
                        // Extract relevant Fourier modes of uj and qj for computations
                        for (int ny=0; ny<Ny_; ++ny) {
                            uk_.set(ny, nu_*uj.cmplx(mx,ny,mz,0));
                            vk_.set(ny, nu_*uj.cmplx(mx,ny,mz,1));
                            wk_.set(ny, nu_*uj.cmplx(mx,ny,mz,2));
                            Pk_.set(ny, qn.cmplx(mx,ny,mz,0));
                        }
                        
                        // (1) Put nu uj" into in R. (Pyk_ is used as tmp workspace)
                        diff2(uk_, Rxk_, Pyk_);
                        diff2(vk_, Ryk_, Pyk_);
                        diff2(wk_, Rzk_, Pyk_);
                        
                        // (2) Put qn' into Pyk (compute y-comp of pressure gradient).
                        diff(Pk_, Pyk_);
                        
                        // (3) Add [1/(nu Cj dt)- kappa2] nu uj - grad qj + Bj/Cj Qj to R.
                        const Real c =
                        1.0/(nu_*C_[j]*dt_) - 4*pi*pi*(square(kx/Lx_) + square(kz/Lz_));
                        const Real B_C = B_[j]/C_[j];
                        const Complex Dx = un.Dx(mx);
                        const Complex Dz = un.Dz(mz);
                        
                        for (int ny=0; ny<Ny_; ++ny) {
                            Rxk_.add(ny, c*uk_[ny] + B_C*Qj_.cmplx(mx,ny,mz,0) - Dx*Pk_[ny]);
                            Ryk_.add(ny, c*vk_[ny] + B_C*Qj_.cmplx(mx,ny,mz,1) - Pyk_[ny]);
                            Rzk_.add(ny, c*wk_[ny] + B_C*Qj_.cmplx(mx,ny,mz,2) - Dz*Pk_[ny]);
                        }
                        
                        // Do the tau solutions
                        if (kx!=0 || kz!=0)
                            //tausolver_[j][mx][mz].solve(uk_,vk_,wk_,Pk_, Rxk_,Ryk_,Rzk_);
                            tausolver_[j][mx][mz].solve_Inhom(uk_,vk_,wk_,Pk_, Rxk_,Ryk_,Rzk_, vr_lower, vr_upper, vi_lower, vi_upper, u_lower, u_upper, w_lower, w_upper);
                            
                            
                        else { // kx,kz == 0,0
                            // Rx has additional terms, nu Uyy at both t=j and t=j+1
                            const Real c = 2*nu_;
                            if (Ubaseyy_.length() > 0)
                                for (int ny=0; ny<Ny_; ++ny)
                                    Rxk_.re[ny] += c*Ubaseyy_[ny];
                            
                            if (flags_.constraint == PressureGradient) {
                                // dPdx is supplied, put dPdx at both t=j and t=j+1 on RHS
                                Rxk_.re[0] -= dPdxAct_ + dPdxRef_;
                                
                                // Solve the tau equations
                                //tausolver_[j][mx][mz].solve(uk_, vk_, wk_, Pk_, Rxk_, Ryk_, Rzk_);
                                tausolver_[j][mx][mz].solve_Inhom(uk_,vk_,wk_,Pk_, Rxk_,Ryk_,Rzk_, vr_lower, vr_upper, vi_lower, vi_upper, u_lower, u_upper, w_lower, w_upper);
                                
                                // Bulk velocity is free variable on LHS solved by tau eqn
                                UbulkAct_ = UbulkBase_ + ubulkBase_ + uk_.re.mean();
                                dPdxAct_ = dPdxRef_;
                            }
                            else { // const bulk velocity
                                // Add the previous time-step's -dPdx to the RHS. The next
                                // timestep's dPdx term appears on LHS as unknown.
                                Rxk_.re[0] -= dPdxAct_;
                                
                                // Use tausolver with additional variable and constraint:
                                // free variable: dPdxAct at next time-step,
                                // constraint:    UbulkBase + ubulkBase + mean(u) = UbulkRef.
                                //tausolver_[j][mx][mz].solve(uk_, vk_, wk_, Pk_, dPdxAct_, Rxk_, Ryk_, Rzk_, UbulkRef_ - UbulkBase_ - ubulkBase_);
                                tausolver_[j][mx][mz].solve_Inhom(uk_, vk_, wk_, Pk_, dPdxAct_, Rxk_, Ryk_, Rzk_, vr_lower, vr_upper, vi_lower, vi_upper, UbulkRef_ - UbulkBase_ - ubulkBase_);
                                
                                UbulkAct_ = UbulkBase_ + ubulkBase_ + uk_.re.mean(); // should == UbulkRef_
                                //UbulkAct_ = UbulkRef_;
                            }
                            // for kx=kz=0, constant term of pressure is arbitrary 3/19/05
                            // Pk_.set(0, Complex(0.0, 0.0));
                        }
                        
                        // Load solutions back into the external 3d data arrays.
                        // Because of FFTW complex symmetries
                        // The 0,0 mode must be real.
                        // For Nx even, the kxmax,0 mode must be real
                        // For Nz even, the 0,kzmax mode must be real
                        // For Nx,Nz even, the kxmax,kzmax mode must be real
                        if ((kx == 0 && kz == 0) ||
                            (Nx_%2 == 0 && kx == kxmax && kz == 0) ||
                            (Nz_%2 == 0 && kz == kzmax && kx == 0) ||
                            (Nx_%2 == 0 && Nz_%2 == 0 && kx == kxmax && kz == kzmax)) {
                            
                            for (int ny=0; ny<Ny_; ++ny) {
                                un.cmplx(mx,ny,mz,0) = Complex(Re(uk_[ny]), 0.0);
                                un.cmplx(mx,ny,mz,1) = Complex(Re(vk_[ny]), 0.0);
                                un.cmplx(mx,ny,mz,2) = Complex(Re(wk_[ny]), 0.0);
                                qn.cmplx(mx,ny,mz,0) = Complex(Re(Pk_[ny]), 0.0);
                            }
                        }
                        // The normal case, for general kx,kz
                        else
                            for (int ny=0; ny<Ny_; ++ny) {
                                un.cmplx(mx,ny,mz,0) = uk_[ny];
                                un.cmplx(mx,ny,mz,1) = vk_[ny];
                                un.cmplx(mx,ny,mz,2) = wk_[ny];
                                qn.cmplx(mx,ny,mz,0) = Pk_[ny];
                            }
                        
                        // And now set the y-aliased modes to zero.
                        for (int ny=Nyd_; ny<Ny_; ++ny) {
                            un.cmplx(mx,ny,mz,0) = 0.0;
                            un.cmplx(mx,ny,mz,1) = 0.0;
                            un.cmplx(mx,ny,mz,2) = 0.0;
                            qn.cmplx(mx,ny,mz,0) = 0.0;
                        }
                    }
                }
            }
            
            
            t_ += dt_;
            
            if (flags_.verbosity == PrintTime || flags_.verbosity == PrintAll)
                *flags_.logstream << t_ << ' ' << flush;
            else if (flags_.verbosity == PrintTicks)
                *flags_.logstream << '.' << flush;
        }
        
        cfl_ = un.CFLfactor(Ubase_);
        cfl_ *= flags_.dealias_xz() ? 2.0*pi/3.0*dt_ : pi*dt_;
        
        // If using dealiasing, set flag in FlowField that compactifies binary IO
        un.setPadded(flags_.dealias_xz());
        qn.setPadded(flags_.dealias_xz());
        
        if (flags_.verbosity == PrintTime || flags_.verbosity == PrintAll)
            *flags_.logstream << endl;
        
        return;
    }
    
    
    
        
    
    
    
    
    
  void RungeKuttaDNS::advance_inhom_CON(Controller& controller, FlowField& un, FlowField& qn, FlowField& BCs, FlowField& Fn, int Nsteps, double*** IO, double*** CStateMat) {
        
        const int kxmax = un.kxmax();
        const int kzmax = un.kzmax();
        
        

        
        for (int n=0; n<Nsteps; ++n) {
            
	  controller.advance_Con_CN(un,qn,BCs,CStateMat,dt_,t_,IO);
            
            //cout << "dt=" << dt_ <<  " " << "t=" << t_ << " " << "BC" << BCs.cmplx(0,0,1,1) << endl;
            
            
            for (int j=0; j<Nsubsteps_; ++j) {
                
                FlowField& uj(un); // Store uj in un during substeps, reflect in notation
                
                // Efficient implementation of
                // Q_{j+1} = A_j Q_j + N(u_j)}  where N = -u grad u
                // Q_{j+1} = A_j Q_j - f(u_j)}  where f =  u grad u
                // Q_j = Q_{j+1}
                Qj_ *= A_[j];
                navierstokesNL(uj, ubase_, Ubase_, Qj1_, tmp_, tmp2_, flags_.nonlinearity);
                Qj_ -= Qj1_; // subtract because navierstokesNL(u) = u grad u = -N(u)
                

		Fn = Qj1_;

                // Update each Fourier mode with time-stepping algorithm
                for (int mx=0; mx<Mx_; ++mx) {
                    const int kx = un.kx(mx);
                    
                    for (int mz=0; mz<Mz_; ++mz) {
                        const int kz = un.kz(mz);
                        
                        Complex u_upper=BCs.cmplx(mx,1,mz,0);
                        Complex u_lower=BCs.cmplx(mx,0,mz,0);
                        Complex v_upper=BCs.cmplx(mx,1,mz,1);
                        Complex v_lower=BCs.cmplx(mx,0,mz,1);
                        Complex w_upper=BCs.cmplx(mx,1,mz,2);
                        Complex w_lower=BCs.cmplx(mx,0,mz,2);
                        
                        Real vr_upper=real(v_upper);
                        Real vi_upper=imag(v_upper);
                        Real vr_lower=real(v_lower);
                        Real vi_lower=imag(v_lower);                       
                        
                        
                        
                        // Zero out the aliased modes
                        if ((kx == kxmax || kz == kzmax) ||
                            (flags_.dealias_xz() && isAliasedMode(kx,kz))) {
                            for (int ny=0; ny<Ny_; ++ny) {
                                un.cmplx(mx,ny,mz,0) = 0.0;
                                un.cmplx(mx,ny,mz,1) = 0.0;
                                un.cmplx(mx,ny,mz,2) = 0.0;
                                qn.cmplx(mx,ny,mz,0) = 0.0;
                            }
                            break;
                        }
                        
                        Rxk_.setToZero();
                        Ryk_.setToZero();
                        Rzk_.setToZero();
                        
                        // Make the following assignments in prep for computation of RHS
                        
                        // nu (uk,vk,wk) = nu ujn(0,1,2)
                        //            Pk = qn
                        
                        // Goal is to compute
                        // R = nu uj" + [1/(Cj dt) - nu kappa2]    uj  - grad qj + Bj/Cj Qj
                        //   = nu uj" + [1/(nu Cj dt) - kappa2] nu uj  - grad qj + Bj/Cj Qj
                        
                        // Extract relevant Fourier modes of uj and qj for computations
                        for (int ny=0; ny<Ny_; ++ny) {
                            uk_.set(ny, nu_*uj.cmplx(mx,ny,mz,0));
                            vk_.set(ny, nu_*uj.cmplx(mx,ny,mz,1));
                            wk_.set(ny, nu_*uj.cmplx(mx,ny,mz,2));
                            Pk_.set(ny, qn.cmplx(mx,ny,mz,0));
                        }
                        
                        // (1) Put nu uj" into in R. (Pyk_ is used as tmp workspace)
                        diff2(uk_, Rxk_, Pyk_);
                        diff2(vk_, Ryk_, Pyk_);
                        diff2(wk_, Rzk_, Pyk_);
                        
                        // (2) Put qn' into Pyk (compute y-comp of pressure gradient).
                        diff(Pk_, Pyk_);
                        
                        // (3) Add [1/(nu Cj dt)- kappa2] nu uj - grad qj + Bj/Cj Qj to R.
                        const Real c =
                        1.0/(nu_*C_[j]*dt_) - 4*pi*pi*(square(kx/Lx_) + square(kz/Lz_));
                        const Real B_C = B_[j]/C_[j];
                        const Complex Dx = un.Dx(mx);
                        const Complex Dz = un.Dz(mz);
                        
                        for (int ny=0; ny<Ny_; ++ny) {
                            Rxk_.add(ny, c*uk_[ny] + B_C*Qj_.cmplx(mx,ny,mz,0) - Dx*Pk_[ny]);
                            Ryk_.add(ny, c*vk_[ny] + B_C*Qj_.cmplx(mx,ny,mz,1) - Pyk_[ny]);
                            Rzk_.add(ny, c*wk_[ny] + B_C*Qj_.cmplx(mx,ny,mz,2) - Dz*Pk_[ny]);
                        }
                        
                        // Do the tau solutions
                        if (kx!=0 || kz!=0)
                            //tausolver_[j][mx][mz].solve(uk_,vk_,wk_,Pk_, Rxk_,Ryk_,Rzk_);
                            tausolver_[j][mx][mz].solve_Inhom(uk_,vk_,wk_,Pk_, Rxk_,Ryk_,Rzk_, vr_lower, vr_upper, vi_lower, vi_upper, u_lower, u_upper, w_lower, w_upper);
                        
                        
                        else { // kx,kz == 0,0
                            // Rx has additional terms, nu Uyy at both t=j and t=j+1
                            const Real c = 2*nu_;
                            if (Ubaseyy_.length() > 0)
                                for (int ny=0; ny<Ny_; ++ny)
                                    Rxk_.re[ny] += c*Ubaseyy_[ny];
                            
                            if (flags_.constraint == PressureGradient) {
                                // dPdx is supplied, put dPdx at both t=j and t=j+1 on RHS
                                Rxk_.re[0] -= dPdxAct_ + dPdxRef_;
                                
                                // Solve the tau equations
                                //tausolver_[j][mx][mz].solve(uk_, vk_, wk_, Pk_, Rxk_, Ryk_, Rzk_);
                                tausolver_[j][mx][mz].solve_Inhom(uk_,vk_,wk_,Pk_, Rxk_,Ryk_,Rzk_, vr_lower, vr_upper, vi_lower, vi_upper, u_lower, u_upper, w_lower, w_upper);
                                
                                // Bulk velocity is free variable on LHS solved by tau eqn
                                UbulkAct_ = UbulkBase_ + ubulkBase_ + uk_.re.mean();
                                dPdxAct_ = dPdxRef_;
                            }
                            else { // const bulk velocity
                                // Add the previous time-step's -dPdx to the RHS. The next
                                // timestep's dPdx term appears on LHS as unknown.
                                Rxk_.re[0] -= dPdxAct_;
                                
                                // Use tausolver with additional variable and constraint:
                                // free variable: dPdxAct at next time-step,
                                // constraint:    UbulkBase + ubulkBase + mean(u) = UbulkRef.
                                //tausolver_[j][mx][mz].solve(uk_, vk_, wk_, Pk_, dPdxAct_, Rxk_, Ryk_, Rzk_, UbulkRef_ - UbulkBase_ - ubulkBase_);
                                tausolver_[j][mx][mz].solve_Inhom(uk_, vk_, wk_, Pk_, dPdxAct_, Rxk_, Ryk_, Rzk_, vr_lower, vr_upper, vi_lower, vi_upper, UbulkRef_ - UbulkBase_ - ubulkBase_);
                                
                                UbulkAct_ = UbulkBase_ + ubulkBase_ + uk_.re.mean(); // should == UbulkRef_
                                //UbulkAct_ = UbulkRef_;
                            }
                            // for kx=kz=0, constant term of pressure is arbitrary 3/19/05
                            // Pk_.set(0, Complex(0.0, 0.0));
                        }
                        
                        // Load solutions back into the external 3d data arrays.
                        // Because of FFTW complex symmetries
                        // The 0,0 mode must be real.
                        // For Nx even, the kxmax,0 mode must be real
                        // For Nz even, the 0,kzmax mode must be real
                        // For Nx,Nz even, the kxmax,kzmax mode must be real
                        if ((kx == 0 && kz == 0) ||
                            (Nx_%2 == 0 && kx == kxmax && kz == 0) ||
                            (Nz_%2 == 0 && kz == kzmax && kx == 0) ||
                            (Nx_%2 == 0 && Nz_%2 == 0 && kx == kxmax && kz == kzmax)) {
                            
                            for (int ny=0; ny<Ny_; ++ny) {
                                un.cmplx(mx,ny,mz,0) = Complex(Re(uk_[ny]), 0.0);
                                un.cmplx(mx,ny,mz,1) = Complex(Re(vk_[ny]), 0.0);
                                un.cmplx(mx,ny,mz,2) = Complex(Re(wk_[ny]), 0.0);
                                qn.cmplx(mx,ny,mz,0) = Complex(Re(Pk_[ny]), 0.0);
                            }
                        }
                        // The normal case, for general kx,kz
                        else
                            for (int ny=0; ny<Ny_; ++ny) {
                                un.cmplx(mx,ny,mz,0) = uk_[ny];
                                un.cmplx(mx,ny,mz,1) = vk_[ny];
                                un.cmplx(mx,ny,mz,2) = wk_[ny];
                                qn.cmplx(mx,ny,mz,0) = Pk_[ny];
                            }
                        
                        // And now set the y-aliased modes to zero.
                        for (int ny=Nyd_; ny<Ny_; ++ny) {
                            un.cmplx(mx,ny,mz,0) = 0.0;
                            un.cmplx(mx,ny,mz,1) = 0.0;
                            un.cmplx(mx,ny,mz,2) = 0.0;
                            qn.cmplx(mx,ny,mz,0) = 0.0;
                        }
                    }
                }
            }
            
            
            t_ += dt_;
            
            if (flags_.verbosity == PrintTime || flags_.verbosity == PrintAll)
                *flags_.logstream << t_ << ' ' << flush;
            else if (flags_.verbosity == PrintTicks)
                *flags_.logstream << '.' << flush;
            
            
            
            //controller.advance_Con_CN(un,qn,BCs,CStateMat,dt_,t_,IO);
	    cfl_ = un.CFLfactor(Ubase_);
            cfl_ *= flags_.dealias_xz() ? 2.0*pi/3.0*dt_ : pi*dt_;
            
            if (cfl_ > 1.0){
                cout << "Error: CFL>1" << endl;
                break;
            }



            
        }
        
        cfl_ = un.CFLfactor(Ubase_);
        cfl_ *= flags_.dealias_xz() ? 2.0*pi/3.0*dt_ : pi*dt_;
        
        // If using dealiasing, set flag in FlowField that compactifies binary IO
        un.setPadded(flags_.dealias_xz());
        qn.setPadded(flags_.dealias_xz());
        
        if (flags_.verbosity == PrintTime || flags_.verbosity == PrintAll)
            *flags_.logstream << endl;
        
        return;
    }


/*

  void RungeKuttaDNS::advance_inhom_CON_SF(Controller& controller, FlowField& un, FlowField& qn, FlowField& BCs, FlowField& Fn, int Nsteps, double*** IO) {
        
        const int kxmax = un.kxmax();
        const int kzmax = un.kzmax();
        
        

        
        for (int n=0; n<Nsteps; ++n) {
            
            controller.advance_Con_SF(un,qn,BCs,dt_,t_,IO);
            
            //cout << "dt=" << dt_ <<  " " << "t=" << t_ << " " << "BC" << BCs.cmplx(0,0,1,1) << endl;
            
            
            for (int j=0; j<Nsubsteps_; ++j) {
                
                FlowField& uj(un); // Store uj in un during substeps, reflect in notation
                
                // Efficient implementation of
                // Q_{j+1} = A_j Q_j + N(u_j)}  where N = -u grad u
                // Q_{j+1} = A_j Q_j - f(u_j)}  where f =  u grad u
                // Q_j = Q_{j+1}
                Qj_ *= A_[j];
                navierstokesNL(uj, ubase_, Ubase_, Qj1_, tmp_, tmp2_, flags_.nonlinearity);
                Qj_ -= Qj1_; // subtract because navierstokesNL(u) = u grad u = -N(u)
                
		Fn = Qj1_;

                // Update each Fourier mode with time-stepping algorithm
                for (int mx=0; mx<Mx_; ++mx) {
                    const int kx = un.kx(mx);
                    
                    for (int mz=0; mz<Mz_; ++mz) {
                        const int kz = un.kz(mz);
                        
                        Complex u_upper=BCs.cmplx(mx,1,mz,0);
                        Complex u_lower=BCs.cmplx(mx,0,mz,0);
                        Complex v_upper=BCs.cmplx(mx,1,mz,1);
                        Complex v_lower=BCs.cmplx(mx,0,mz,1);
                        Complex w_upper=BCs.cmplx(mx,1,mz,2);
                        Complex w_lower=BCs.cmplx(mx,0,mz,2);
                        
                        Real vr_upper=real(v_upper);
                        Real vi_upper=imag(v_upper);
                        Real vr_lower=real(v_lower);
                        Real vi_lower=imag(v_lower);                       
                        
                        
                        
                        // Zero out the aliased modes
                        if ((kx == kxmax || kz == kzmax) ||
                            (flags_.dealias_xz() && isAliasedMode(kx,kz))) {
                            for (int ny=0; ny<Ny_; ++ny) {
                                un.cmplx(mx,ny,mz,0) = 0.0;
                                un.cmplx(mx,ny,mz,1) = 0.0;
                                un.cmplx(mx,ny,mz,2) = 0.0;
                                qn.cmplx(mx,ny,mz,0) = 0.0;
                            }
                            break;
                        }
                        
                        Rxk_.setToZero();
                        Ryk_.setToZero();
                        Rzk_.setToZero();
                        
                        // Make the following assignments in prep for computation of RHS
                        
                        // nu (uk,vk,wk) = nu ujn(0,1,2)
                        //            Pk = qn
                        
                        // Goal is to compute
                        // R = nu uj" + [1/(Cj dt) - nu kappa2]    uj  - grad qj + Bj/Cj Qj
                        //   = nu uj" + [1/(nu Cj dt) - kappa2] nu uj  - grad qj + Bj/Cj Qj
                        
                        // Extract relevant Fourier modes of uj and qj for computations
                        for (int ny=0; ny<Ny_; ++ny) {
                            uk_.set(ny, nu_*uj.cmplx(mx,ny,mz,0));
                            vk_.set(ny, nu_*uj.cmplx(mx,ny,mz,1));
                            wk_.set(ny, nu_*uj.cmplx(mx,ny,mz,2));
                            Pk_.set(ny, qn.cmplx(mx,ny,mz,0));
                        }
                        
                        // (1) Put nu uj" into in R. (Pyk_ is used as tmp workspace)
                        diff2(uk_, Rxk_, Pyk_);
                        diff2(vk_, Ryk_, Pyk_);
                        diff2(wk_, Rzk_, Pyk_);
                        
                        // (2) Put qn' into Pyk (compute y-comp of pressure gradient).
                        diff(Pk_, Pyk_);
                        
                        // (3) Add [1/(nu Cj dt)- kappa2] nu uj - grad qj + Bj/Cj Qj to R.
                        const Real c =
                        1.0/(nu_*C_[j]*dt_) - 4*pi*pi*(square(kx/Lx_) + square(kz/Lz_));
                        const Real B_C = B_[j]/C_[j];
                        const Complex Dx = un.Dx(mx);
                        const Complex Dz = un.Dz(mz);
                        
                        for (int ny=0; ny<Ny_; ++ny) {
                            Rxk_.add(ny, c*uk_[ny] + B_C*Qj_.cmplx(mx,ny,mz,0) - Dx*Pk_[ny]);
                            Ryk_.add(ny, c*vk_[ny] + B_C*Qj_.cmplx(mx,ny,mz,1) - Pyk_[ny]);
                            Rzk_.add(ny, c*wk_[ny] + B_C*Qj_.cmplx(mx,ny,mz,2) - Dz*Pk_[ny]);
                        }
                        
                        // Do the tau solutions
                        if (kx!=0 || kz!=0)
                            //tausolver_[j][mx][mz].solve(uk_,vk_,wk_,Pk_, Rxk_,Ryk_,Rzk_);
                            tausolver_[j][mx][mz].solve_Inhom(uk_,vk_,wk_,Pk_, Rxk_,Ryk_,Rzk_, vr_lower, vr_upper, vi_lower, vi_upper, u_lower, u_upper, w_lower, w_upper);
                        
                        
                        else { // kx,kz == 0,0
                            // Rx has additional terms, nu Uyy at both t=j and t=j+1
                            const Real c = 2*nu_;
                            if (Ubaseyy_.length() > 0)
                                for (int ny=0; ny<Ny_; ++ny)
                                    Rxk_.re[ny] += c*Ubaseyy_[ny];
                            
                            if (flags_.constraint == PressureGradient) {
                                // dPdx is supplied, put dPdx at both t=j and t=j+1 on RHS
                                Rxk_.re[0] -= dPdxAct_ + dPdxRef_;
                                
                                // Solve the tau equations
                                //tausolver_[j][mx][mz].solve(uk_, vk_, wk_, Pk_, Rxk_, Ryk_, Rzk_);
                                tausolver_[j][mx][mz].solve_Inhom(uk_,vk_,wk_,Pk_, Rxk_,Ryk_,Rzk_, vr_lower, vr_upper, vi_lower, vi_upper, u_lower, u_upper, w_lower, w_upper);
                                
                                // Bulk velocity is free variable on LHS solved by tau eqn
                                UbulkAct_ = UbulkBase_ + ubulkBase_ + uk_.re.mean();
                                dPdxAct_ = dPdxRef_;
                            }
                            else { // const bulk velocity
                                // Add the previous time-step's -dPdx to the RHS. The next
                                // timestep's dPdx term appears on LHS as unknown.
                                Rxk_.re[0] -= dPdxAct_;
                                
                                // Use tausolver with additional variable and constraint:
                                // free variable: dPdxAct at next time-step,
                                // constraint:    UbulkBase + ubulkBase + mean(u) = UbulkRef.
                                //tausolver_[j][mx][mz].solve(uk_, vk_, wk_, Pk_, dPdxAct_, Rxk_, Ryk_, Rzk_, UbulkRef_ - UbulkBase_ - ubulkBase_);
                                tausolver_[j][mx][mz].solve_Inhom(uk_, vk_, wk_, Pk_, dPdxAct_, Rxk_, Ryk_, Rzk_, vr_lower, vr_upper, vi_lower, vi_upper, UbulkRef_ - UbulkBase_ - ubulkBase_);
                                
                                UbulkAct_ = UbulkBase_ + ubulkBase_ + uk_.re.mean(); // should == UbulkRef_
                                //UbulkAct_ = UbulkRef_;
                            }
                            // for kx=kz=0, constant term of pressure is arbitrary 3/19/05
                            // Pk_.set(0, Complex(0.0, 0.0));
                        }
                        
                        // Load solutions back into the external 3d data arrays.
                        // Because of FFTW complex symmetries
                        // The 0,0 mode must be real.
                        // For Nx even, the kxmax,0 mode must be real
                        // For Nz even, the 0,kzmax mode must be real
                        // For Nx,Nz even, the kxmax,kzmax mode must be real
                        if ((kx == 0 && kz == 0) ||
                            (Nx_%2 == 0 && kx == kxmax && kz == 0) ||
                            (Nz_%2 == 0 && kz == kzmax && kx == 0) ||
                            (Nx_%2 == 0 && Nz_%2 == 0 && kx == kxmax && kz == kzmax)) {
                            
                            for (int ny=0; ny<Ny_; ++ny) {
                                un.cmplx(mx,ny,mz,0) = Complex(Re(uk_[ny]), 0.0);
                                un.cmplx(mx,ny,mz,1) = Complex(Re(vk_[ny]), 0.0);
                                un.cmplx(mx,ny,mz,2) = Complex(Re(wk_[ny]), 0.0);
                                qn.cmplx(mx,ny,mz,0) = Complex(Re(Pk_[ny]), 0.0);
                            }
                        }
                        // The normal case, for general kx,kz
                        else
                            for (int ny=0; ny<Ny_; ++ny) {
                                un.cmplx(mx,ny,mz,0) = uk_[ny];
                                un.cmplx(mx,ny,mz,1) = vk_[ny];
                                un.cmplx(mx,ny,mz,2) = wk_[ny];
                                qn.cmplx(mx,ny,mz,0) = Pk_[ny];
                            }
                        
                        // And now set the y-aliased modes to zero.
                        for (int ny=Nyd_; ny<Ny_; ++ny) {
                            un.cmplx(mx,ny,mz,0) = 0.0;
                            un.cmplx(mx,ny,mz,1) = 0.0;
                            un.cmplx(mx,ny,mz,2) = 0.0;
                            qn.cmplx(mx,ny,mz,0) = 0.0;
                        }
                    }
                }
            }
            
            
            t_ += dt_;
            
            if (flags_.verbosity == PrintTime || flags_.verbosity == PrintAll)
                *flags_.logstream << t_ << ' ' << flush;
            else if (flags_.verbosity == PrintTicks)
                *flags_.logstream << '.' << flush;
            
            
            
            //controller.advance_Con_CN(un,qn,BCs,CStateMat,dt_,t_,IO);
	    cfl_ = un.CFLfactor(Ubase_);
            cfl_ *= flags_.dealias_xz() ? 2.0*pi/3.0*dt_ : pi*dt_;
            
            if (cfl_ > 1.0){
                cout << "Error: CFL>1" << endl;
                break;
            }



            
        }
        
        cfl_ = un.CFLfactor(Ubase_);
        cfl_ *= flags_.dealias_xz() ? 2.0*pi/3.0*dt_ : pi*dt_;
        
        // If using dealiasing, set flag in FlowField that compactifies binary IO
        un.setPadded(flags_.dealias_xz());
        qn.setPadded(flags_.dealias_xz());
        
        if (flags_.verbosity == PrintTime || flags_.verbosity == PrintAll)
            *flags_.logstream << endl;
        
        return;
    }
    */
    

    /**************************************************************************
     
     ---End: Peter H Heins---
     
     **************************************************************************/  
    
    

    
    
    
    
    
    
    // ==============================================================
    // CNAB-style algorithms
    
    CNABstyleDNS::CNABstyleDNS()
    :
    DNSAlgorithm(),
    full_(false),
    fj1_(),
    fj_()
    {}
    
    CNABstyleDNS::CNABstyleDNS(const CNABstyleDNS& dns)
    :
    DNSAlgorithm(dns),
    Nsubsteps_(dns.Nsubsteps_),
    full_(dns.full_),
    fj1_(dns.fj1_),
    fj_(dns.fj_),
    alpha_(dns.alpha_),
    beta_(dns.beta_),
    gamma_(dns.gamma_),
    zeta_(dns.zeta_)
    {
        // Allocate memory for [Nsubsteps x Mx_ x Mz_] Tausolver arrays
        // and copy tausolvers from dns argument
        tausolver_ = new TauSolver**[Nsubsteps_];    // new #1
        for (int j=0; j<Nsubsteps_; ++j) {
            tausolver_[j] = new TauSolver*[Mx_];       // new #2
            for (int mx=0; mx<Mx_; ++mx) {
                tausolver_[j][mx] = new TauSolver[Mz_];  // new #3
                for (int mz=0; mz<Mz_; ++mz)
                    tausolver_[j][mx][mz] = dns.tausolver_[j][mx][mz];
                
            }
        }
    }
    /*
    CNABstyleDNS::CNABstyleDNS(const FlowField& u, const ChebyCoeff& Ubase,
                               Real nu, Real dt, const DNSFlags& flags, Real t)
    :
    DNSAlgorithm(u,Ubase,nu,dt,flags, t),
    full_(false),
    fj1_(u),
    fj_(u)
    {
        fj1_.setToZero();
        fj_.setToZero();
        
        TimeStepMethod algorithm = flags.timestepping;
        switch (algorithm) {
                // The classic Crank-Nicolson/Adams-bashforth algorithm
            case CNAB2:
                order_ = 2;
                Nsubsteps_ = 1;
                Ninitsteps_ = 1;
                full_ = false;
                alpha_.resize(Nsubsteps_);
                beta_.resize(Nsubsteps_);
                gamma_.resize(Nsubsteps_);
                zeta_.resize(Nsubsteps_);
                alpha_[0] = 0.5;
                beta_[0]  = 0.5;
                gamma_[0] = 1.5;
                zeta_[0]  = -0.5;
                break;
                
                // Constants taken from P.R. Spalart, R.D. Moser, M.M. Rogers,
                // Spectral methods for the Navier-Stokes equations with one infinite and
                // two periodic directions, J. Comp. Phys. 96, 297–324 (1990).
                //
            case SMRK2:
                order_ = 2;
                Nsubsteps_ = 3;
                Ninitsteps_ = 0;
                full_ = true;
                alpha_.resize(Nsubsteps_);
                beta_.resize(Nsubsteps_);
                gamma_.resize(Nsubsteps_);
                zeta_.resize(Nsubsteps_);
                alpha_[0] = 29.0/96.0;  alpha_[1] = -3.0/40.0;  alpha_[2] = 1.0/6.0;
                beta_[0]  = 37.0/160.0;  beta_[1] =  5.0/24.0;  beta_[2]  = 1.0/6.0;
                gamma_[0] = 8.0/15.0;   gamma_[1] =  5.0/12.0;  gamma_[2] = 3.0/4.0;
                zeta_[0]  = 0.0;         zeta_[1] = -17.0/60.0;  zeta_[2] = -5.0/12.0;
                break;
            default:
                cerr << "CNABstyleDNS::CNABstyleDNS(un,Ubase,nu,dt,flags,t0)\n"
                << "error: flags.timestepping == " << algorithm
                << " is not a CNAB-style algorithm." << endl;
                exit(1);
        }
        
        // Allocate memory for [Nsubsteps x Mx_ x Mz_] Tausolver array
        tausolver_ = new TauSolver**[Nsubsteps_];    // new #1
        for (int j=0; j<Nsubsteps_; ++j) {
            tausolver_[j] = new TauSolver*[Mx_];       // new #2
            for (int mx=0; mx<Mx_; ++mx)
                tausolver_[j][mx] = new TauSolver[Mz_];  // new #3
        }
        reset_dt(dt_);
    }
    */
    
    
    
    /////////////////////////////////////////////////////////////
    // Added by P. Heins 10/13
    CNABstyleDNS::CNABstyleDNS(const FlowField& u, const ChebyCoeff& Ubase,
                               Real nu, Real dt, const DNSFlags& flags, Real t, bool controlled)
    :
    DNSAlgorithm(u,Ubase,nu,dt,flags, t,controlled),
    full_(false),
    fj1_(u),
    fj_(u)
    {
        fj1_.setToZero();
        fj_.setToZero();
        
        TimeStepMethod algorithm = flags.timestepping;
        switch (algorithm) {
                // The classic Crank-Nicolson/Adams-bashforth algorithm
            case CNAB2:
                order_ = 2;
                Nsubsteps_ = 1;
                Ninitsteps_ = 1;
                full_ = false;
                alpha_.resize(Nsubsteps_);
                beta_.resize(Nsubsteps_);
                gamma_.resize(Nsubsteps_);
                zeta_.resize(Nsubsteps_);
                alpha_[0] = 0.5;
                beta_[0]  = 0.5;
                gamma_[0] = 1.5;
                zeta_[0]  = -0.5;
                break;
                
                // Constants taken from P.R. Spalart, R.D. Moser, M.M. Rogers,
                // Spectral methods for the Navier-Stokes equations with one infinite and
                // two periodic directions, J. Comp. Phys. 96, 297–324 (1990).
                //
            case SMRK2:
                order_ = 2;
                Nsubsteps_ = 3;
                Ninitsteps_ = 0;
                full_ = true;
                alpha_.resize(Nsubsteps_);
                beta_.resize(Nsubsteps_);
                gamma_.resize(Nsubsteps_);
                zeta_.resize(Nsubsteps_);
                alpha_[0] = 29.0/96.0;  alpha_[1] = -3.0/40.0;  alpha_[2] = 1.0/6.0;
                beta_[0]  = 37.0/160.0;  beta_[1] =  5.0/24.0;  beta_[2]  = 1.0/6.0;
                gamma_[0] = 8.0/15.0;   gamma_[1] =  5.0/12.0;  gamma_[2] = 3.0/4.0;
                zeta_[0]  = 0.0;         zeta_[1] = -17.0/60.0;  zeta_[2] = -5.0/12.0;
                break;
            default:
                cerr << "CNABstyleDNS::CNABstyleDNS(un,Ubase,nu,dt,flags,t0)\n"
                << "error: flags.timestepping == " << algorithm
                << " is not a CNAB-style algorithm." << endl;
                exit(1);
        }
        
        // Allocate memory for [Nsubsteps x Mx_ x Mz_] Tausolver array
        tausolver_ = new TauSolver**[Nsubsteps_];    // new #1
        for (int j=0; j<Nsubsteps_; ++j) {
            tausolver_[j] = new TauSolver*[Mx_];       // new #2
            for (int mx=0; mx<Mx_; ++mx)
                tausolver_[j][mx] = new TauSolver[Mz_];  // new #3
        }
        
        if (controlled){
            FlowField BC(Nx_,2,Nz_,3,Lx_,Lz_,a_,b_);
            BC.makeState(Spectral,Physical);
            reset_dtIH(dt_,BC);            
        }
        else {
            reset_dt(dt_);
        }

        
    }

    /////////////////////////////////////////////////////////////
    
    
    
    
    
    CNABstyleDNS::~CNABstyleDNS() {
        if (tausolver_) {
            for (int j=0; j<Nsubsteps_; ++j) {
                for (int mx=0; mx<Mx_; ++mx) {
                    delete[] tausolver_[j][mx];  // undo new #3
                    tausolver_[j][mx] = 0;
                }
                delete[] tausolver_[j];        // undo new #2
                tausolver_[j] = 0;
            }
            delete[] tausolver_;             // undo new #1
            tausolver_ = 0;
        }
    }
    
    DNSAlgorithm* CNABstyleDNS::clone() const {
        return new CNABstyleDNS(*this);
    }
    
    
    void CNABstyleDNS::reset_dt(Real dt) {
        cfl_ *= dt/dt_;
        //nu_ = nu;
        dt_ = dt;
        const Real c = 4.0*square(pi)*nu_;
        const int kxmax = tmp_.kxmax();
        const int kzmax = tmp_.kzmax();
        
        // This loop replaces the TauSolver objects at tausolver_[i][mx][mz]
        // with new TauSolver objects configured with appropriate parameters
        for (int j=0; j<Nsubsteps_; ++j) {
            for (int mx=0; mx<Mx_; ++mx) {
                int kx = tmp_.kx(mx);
                for (int mz=0; mz<Mz_; ++mz) {
                    int kz = tmp_.kz(mz);
                    
                    Real lambda = 1.0/(beta_[j]*dt_)+c*(square(kx/Lx_)+square(kz/Lz_));
                    
                    if ((kx != kxmax) && (kz != kzmax) &&
                        (!flags_.dealias_xz() || !isAliasedMode(kx,kz)))
                        
                        tausolver_[j][mx][mz] =
                        TauSolver(kx, kz, Lx_, Lz_, a_, b_, lambda, nu_, Ny_,
                                  flags_.taucorrection);
                }
            }
        }
        // For some forms of CNABstyle, need to reinitialize again
        switch (flags_.timestepping) {
            case CNAB2:
                full_ = false;
                break;
            default:
                ;
        }
    }
    
    void CNABstyleDNS::reset_dtIH(Real dt, FlowField& BCs) {
        cfl_ *= dt/dt_;
        //nu_ = nu;
        dt_ = dt;
        const Real c = 4.0*square(pi)*nu_;
        const int kxmax = tmp_.kxmax();
        const int kzmax = tmp_.kzmax();
        
        // This loop replaces the TauSolver objects at tausolver_[i][mx][mz]
        // with new TauSolver objects configured with appropriate parameters
        for (int j=0; j<Nsubsteps_; ++j) {
            for (int mx=0; mx<Mx_; ++mx) {
                int kx = tmp_.kx(mx);
                for (int mz=0; mz<Mz_; ++mz) {
                    int kz = tmp_.kz(mz);
                    
                    Complex v_upper=BCs.cmplx(mx,1,mz,1);
                    Complex v_lower=BCs.cmplx(mx,0,mz,1);
                    
                    double vr_upper=real(v_upper);
                    double vi_upper=imag(v_upper);
                    double vr_lower=real(v_lower);
                    double vi_lower=imag(v_lower);         
                    
                    
                    
                    Real lambda = 1.0/(beta_[j]*dt_)+c*(square(kx/Lx_)+square(kz/Lz_));
                    
                    if ((kx != kxmax) && (kz != kzmax) &&
                        (!flags_.dealias_xz() || !isAliasedMode(kx,kz)))
                        
                        //tausolver_[j][mx][mz] = TauSolver(kx, kz, Lx_, Lz_, a_, b_, lambda, nu_, Nyd_, flags_.taucorrection);
                        tausolver_[j][mx][mz] = TauSolver(kx,kz,Lx_,Lz_,a_,b_,lambda,nu_,Ny_, vr_lower, vi_lower, vr_upper, vi_upper, flags_.taucorrection);
                }
            }
        }
        // For some forms of CNABstyle, need to reinitialize again
        switch (flags_.timestepping) {
            case CNAB2:
                full_ = false;
                break;
            default:
                ;
        }
    }

    
    bool CNABstyleDNS::push(const FlowField& un) {
        swap(fj_, fj1_);
        navierstokesNL(un, ubase_, Ubase_, fj_, tmp_, tmp2_, flags_.nonlinearity);
        
        t_ += dt_;
        full_ = true;
        return full_;
    }
    
    void CNABstyleDNS::printStack() const {
        //os << "CNABstyleDNS::printStack() {" << endl;
        //printf("L2(fj) =%13.10f\n", 	 L2Norm(fj_));
        //printf("L2(fj1)=%13.10f\n", 	 L2Norm(fj1_));
        //os << "}" << endl;
        
    }
    
    bool CNABstyleDNS::full() const {
        return full_;
    }
    
    void CNABstyleDNS::project() {
        fj_.project(flags_.symmetries);
        //project(symm, fj1_); // no need; gets overwritten w fj_ first thing in advance
    }
    void CNABstyleDNS::operator*=(const FieldSymmetry& sigma) {
        fj_ *= sigma;
        //project(symm, fj1_); // no need; gets overwritten w fj_ first thing in advance
    }
    
  void CNABstyleDNS::advance(FlowField& un, FlowField& qn, int Nsteps) {
        const int kxmax = tmp_.kxmax();
        const int kzmax = tmp_.kzmax();
        
        for (int n=0; n<Nsteps; ++n) {
            for (int j=0; j<Nsubsteps_; ++j) {
                
                // Store substeps uj,qj in un,qn; reflect this in notation
                FlowField& uj(un);
                FlowField& qj(qn);
                
                swap(fj_, fj1_);
                navierstokesNL(un, ubase_, Ubase_, fj_, tmp_, tmp2_, flags_.nonlinearity);
                
	
		
                // Set convenience variables
                Real a_b = alpha_[j]/beta_[j];
                Real g_b = gamma_[j]/beta_[j];
                Real z_b = zeta_[j]/beta_[j];
                Real anu_b = a_b*nu_;
                Real anu = alpha_[j]*nu_;
                
                // Update each Fourier mode with time-stepping algorithm
                for (int mx=0; mx<Mx_; ++mx) {
                    const int kx = uj.kx(mx);
                    
                    for (int mz=0; mz<Mz_; ++mz) {
                        const int kz = uj.kz(mz);
                        
                        // Zero out the aliased modes
                        if ((kx == kxmax || kz == kzmax) ||
                            (flags_.dealias_xz() && isAliasedMode(kx,kz))) {
                            for (int ny=0; ny<Ny_; ++ny) {
                                uj.cmplx(mx,ny,mz,0) = 0.0;
                                uj.cmplx(mx,ny,mz,1) = 0.0;
                                uj.cmplx(mx,ny,mz,2) = 0.0;
                                qj.cmplx(mx,ny,mz,0) = 0.0;
                            }
                            break;
                        }
                        
                        Rxk_.setToZero();
                        Ryk_.setToZero();
                        Rzk_.setToZero();
                        
                        // Goal is to compute
                        // R = a/b nu uj" + [1/(b dt)- a/b nu kappa2] uj - a/b grad qj
                        //     - g/b fj - z/b fj1
                        //
                        //   = a/b nu uj" + [1/(nu a dt)- kappa2] a/b nu uj
                        //     - a/b grad qj - g/b fj - z/b fj1
                        
                        // Extract relevant Fourier modes of uj and qj for computations
                        // set (uk,vk,wk) to a/b nu uj, Pk to a/b qj
                        
                        for (int ny=0; ny<Ny_; ++ny) {
                            uk_.set(ny, anu_b*uj.cmplx(mx,ny,mz,0));
                            vk_.set(ny, anu_b*uj.cmplx(mx,ny,mz,1));
                            wk_.set(ny, anu_b*uj.cmplx(mx,ny,mz,2));
                            Pk_.set(ny, a_b*qj.cmplx(mx,ny,mz,0));
                        }
                        
                        // (1) Put a/b nu uj" into in R. (Pyk_ is used as tmp workspace)
                        diff2(uk_, Rxk_, Pyk_);
                        diff2(vk_, Ryk_, Pyk_);
                        diff2(wk_, Rzk_, Pyk_);
                        
                        // (2) Put a/b qj' into Pyk (compute y-comp of pressure gradient).
                        diff(Pk_, Pyk_);
                        
                        // (3) Add [1/(nu a dt)- kappa2] a/b nu uj - a/b grad qj - g/b fj - z/b fj1
                        // to R, completing calculation of R
                        const Real kappa2 = 4*pi*pi*(square(kx/Lx_) + square(kz/Lz_));
                        const Real c = 1.0/(anu*dt_) - kappa2;
                        const Complex Dx = un.Dx(mx);
                        const Complex Dz = un.Dz(mz);
                        for (int ny=0; ny<Ny_; ++ny) {
                            Rxk_.add(ny, c*uk_[ny] - Dx*Pk_[ny]
                                     - g_b*fj_.cmplx(mx,ny,mz,0) - z_b*fj1_.cmplx(mx,ny,mz,0));
                            Ryk_.add(ny, c*vk_[ny] - Pyk_[ny]
                                     - g_b*fj_.cmplx(mx,ny,mz,1) - z_b*fj1_.cmplx(mx,ny,mz,1));
                            Rzk_.add(ny, c*wk_[ny] - Dz*Pk_[ny]
                                     - g_b*fj_.cmplx(mx,ny,mz,2) - z_b*fj1_.cmplx(mx,ny,mz,2));
                        }
                        
                        // Solve the tau solutions
                        if (kx!=0 || kz!=0)
                            tausolver_[j][mx][mz].solve(uk_,vk_,wk_,Pk_, Rxk_,Ryk_,Rzk_);
                        
                        else { // kx,kz == 0,0
                            // Rx has additional terms, nu Uyy at both t=j and t=j+1
                            Real a_b = alpha_[j]/beta_[j];
                            Real c = nu_*(a_b+1.0);
                            if (Ubaseyy_.length() > 0)
                                for (int ny=0; ny<Ny_; ++ny)
                                    Rxk_.re[ny] += c*Ubaseyy_[ny];
                            
                            if (flags_.constraint == PressureGradient) {
                                // dPdx is supplied, put dPdx at both t=j and t=j+1 on RHS
                                Rxk_.re[0] -= a_b*dPdxAct_ + dPdxRef_;
                                
                                // Solve the tau equations
                                tausolver_[j][mx][mz].solve(uk_, vk_, wk_, Pk_, Rxk_, Ryk_, Rzk_);
                                // Bulk velocity is free variable on LHS solved by tau eqn
                                UbulkAct_ = UbulkBase_ + uk_.re.mean();
                                dPdxAct_ = dPdxRef_;
                            }
                            else { // const bulk velocity
                                // Add the previous time-step's -dPdx to the RHS. The next
                                // timestep's dPdx term appears on LHS as unknown.
                                Rxk_.re[0] -= a_b*dPdxAct_;
                                
                                // Use tausolver with additional variable and constraint:
                                // free variable: dPdxAct at next time-step,
                                // constraint:    UbulkBase + mean(u) = UbulkRef.
                                tausolver_[j][mx][mz].solve(uk_, vk_, wk_, Pk_, dPdxAct_,
                                                            Rxk_, Ryk_, Rzk_,
                                                            UbulkRef_ - UbulkBase_);
                                
                                UbulkAct_ = UbulkBase_ + uk_.re.mean(); // should == UbulkRef_
                                //UbulkAct_ = UbulkRef_;
                            }
                            // for kx=kz=0, constant term of pressure is arbitrary 3/19/05
                            // Pk_.set(0, Complex(0.0, 0.0));
                        }
                        
                        // Load solutions back into the external 3d data arrays.
                        // Because of FFTW complex symmetries
                        // The 0,0 mode must be real.
                        // For Nx even, the kxmax,0 mode must be real
                        // For Nz even, the 0,kzmax mode must be real
                        // For Nx,Nz even, the kxmax,kzmax mode must be real
                        if ((kx == 0 && kz == 0) ||
                            (Nx_%2 == 0 && kx == kxmax && kz == 0) ||
                            (Nz_%2 == 0 && kz == kzmax && kx == 0) ||
                            (Nx_%2 == 0 && Nz_%2 == 0 && kx == kxmax && kz == kzmax)) {
                            
                            for (int ny=0; ny<Ny_; ++ny) {
                                un.cmplx(mx,ny,mz,0) = Complex(Re(uk_[ny]), 0.0);
                                un.cmplx(mx,ny,mz,1) = Complex(Re(vk_[ny]), 0.0);
                                un.cmplx(mx,ny,mz,2) = Complex(Re(wk_[ny]), 0.0);
                                qn.cmplx(mx,ny,mz,0) = Complex(Re(Pk_[ny]), 0.0);
                            }
                        }
                        // The normal case, for general kx,kz
                        else
                            for (int ny=0; ny<Ny_; ++ny) {
                                un.cmplx(mx,ny,mz,0) = uk_[ny];
                                un.cmplx(mx,ny,mz,1) = vk_[ny];
                                un.cmplx(mx,ny,mz,2) = wk_[ny];
                                qn.cmplx(mx,ny,mz,0) = Pk_[ny];
                            }
                        
                        // And now set the y-aliased modes to zero.
                        for (int ny=Nyd_; ny<Ny_; ++ny) {
                            un.cmplx(mx,ny,mz,0) = 0.0;
                            un.cmplx(mx,ny,mz,1) = 0.0;
                            un.cmplx(mx,ny,mz,2) = 0.0;
                            qn.cmplx(mx,ny,mz,0) = 0.0;
                        }
                    }
                }
            }
            t_ += dt_;
            if (flags_.verbosity == PrintTime || flags_.verbosity == PrintAll)
                *flags_.logstream << t_ << ' ' << flush;
            else if (flags_.verbosity == PrintTicks)
                *flags_.logstream << '.' << flush;
        }
        
        cfl_ = un.CFLfactor(Ubase_);
        cfl_ *= flags_.dealias_xz() ? 2.0*pi/3.0*dt_ : pi*dt_;
        
        // If using dealiasing, set flag in FlowField that compactifies binary IO
        un.setPadded(flags_.dealias_xz());
        qn.setPadded(flags_.dealias_xz());
        
        if (flags_.verbosity == PrintTime || flags_.verbosity == PrintAll)
            *flags_.logstream << endl;
        
        return;
    }


    
    /************************************************************************
     ***********************************************************************
     
     Written by: Peter H Heins (Postgraduate Research Student)
     
     Date: January 2012
     
     Institution: University of Sheffield (ACSE Department)
     
     Purpose: This function implements imhomogeneous BCs into the CNAB algorithm,
     
     Based on code written by Binh Lieu, University of Minnesota March 2009
     
     ***********************************************************************
     *************************************************************************/

    
    void CNABstyleDNS::advance_NL(FlowField& un, FlowField& qn, FlowField& Fn, int Nsteps) {
        const int kxmax = tmp_.kxmax();
        const int kzmax = tmp_.kzmax();
        
        for (int n=0; n<Nsteps; ++n) {
            for (int j=0; j<Nsubsteps_; ++j) {
                
                // Store substeps uj,qj in un,qn; reflect this in notation
                FlowField& uj(un);
                FlowField& qj(qn);
                
                swap(fj_, fj1_);
                navierstokesNL(un, ubase_, Ubase_, fj_, tmp_, tmp2_, flags_.nonlinearity);
                
		Fn = fj_;
		
                // Set convenience variables
                Real a_b = alpha_[j]/beta_[j];
                Real g_b = gamma_[j]/beta_[j];
                Real z_b = zeta_[j]/beta_[j];
                Real anu_b = a_b*nu_;
                Real anu = alpha_[j]*nu_;
                
                // Update each Fourier mode with time-stepping algorithm
                for (int mx=0; mx<Mx_; ++mx) {
                    const int kx = uj.kx(mx);
                    
                    for (int mz=0; mz<Mz_; ++mz) {
                        const int kz = uj.kz(mz);
                        
                        // Zero out the aliased modes
                        if ((kx == kxmax || kz == kzmax) ||
                            (flags_.dealias_xz() && isAliasedMode(kx,kz))) {
                            for (int ny=0; ny<Ny_; ++ny) {
                                uj.cmplx(mx,ny,mz,0) = 0.0;
                                uj.cmplx(mx,ny,mz,1) = 0.0;
                                uj.cmplx(mx,ny,mz,2) = 0.0;
                                qj.cmplx(mx,ny,mz,0) = 0.0;
                            }
                            break;
                        }
                        
                        Rxk_.setToZero();
                        Ryk_.setToZero();
                        Rzk_.setToZero();
                        
                        // Goal is to compute
                        // R = a/b nu uj" + [1/(b dt)- a/b nu kappa2] uj - a/b grad qj
                        //     - g/b fj - z/b fj1
                        //
                        //   = a/b nu uj" + [1/(nu a dt)- kappa2] a/b nu uj
                        //     - a/b grad qj - g/b fj - z/b fj1
                        
                        // Extract relevant Fourier modes of uj and qj for computations
                        // set (uk,vk,wk) to a/b nu uj, Pk to a/b qj
                        
                        for (int ny=0; ny<Ny_; ++ny) {
                            uk_.set(ny, anu_b*uj.cmplx(mx,ny,mz,0));
                            vk_.set(ny, anu_b*uj.cmplx(mx,ny,mz,1));
                            wk_.set(ny, anu_b*uj.cmplx(mx,ny,mz,2));
                            Pk_.set(ny, a_b*qj.cmplx(mx,ny,mz,0));
                        }
                        
                        // (1) Put a/b nu uj" into in R. (Pyk_ is used as tmp workspace)
                        diff2(uk_, Rxk_, Pyk_);
                        diff2(vk_, Ryk_, Pyk_);
                        diff2(wk_, Rzk_, Pyk_);
                        
                        // (2) Put a/b qj' into Pyk (compute y-comp of pressure gradient).
                        diff(Pk_, Pyk_);
                        
                        // (3) Add [1/(nu a dt)- kappa2] a/b nu uj - a/b grad qj - g/b fj - z/b fj1
                        // to R, completing calculation of R
                        const Real kappa2 = 4*pi*pi*(square(kx/Lx_) + square(kz/Lz_));
                        const Real c = 1.0/(anu*dt_) - kappa2;
                        const Complex Dx = un.Dx(mx);
                        const Complex Dz = un.Dz(mz);
                        for (int ny=0; ny<Ny_; ++ny) {
                            Rxk_.add(ny, c*uk_[ny] - Dx*Pk_[ny]
                                     - g_b*fj_.cmplx(mx,ny,mz,0) - z_b*fj1_.cmplx(mx,ny,mz,0));
                            Ryk_.add(ny, c*vk_[ny] - Pyk_[ny]
                                     - g_b*fj_.cmplx(mx,ny,mz,1) - z_b*fj1_.cmplx(mx,ny,mz,1));
                            Rzk_.add(ny, c*wk_[ny] - Dz*Pk_[ny]
                                     - g_b*fj_.cmplx(mx,ny,mz,2) - z_b*fj1_.cmplx(mx,ny,mz,2));
                        }
                        
                        // Solve the tau solutions
                        if (kx!=0 || kz!=0)
                            tausolver_[j][mx][mz].solve(uk_,vk_,wk_,Pk_, Rxk_,Ryk_,Rzk_);
                        
                        else { // kx,kz == 0,0
                            // Rx has additional terms, nu Uyy at both t=j and t=j+1
                            Real a_b = alpha_[j]/beta_[j];
                            Real c = nu_*(a_b+1.0);
                            if (Ubaseyy_.length() > 0)
                                for (int ny=0; ny<Ny_; ++ny)
                                    Rxk_.re[ny] += c*Ubaseyy_[ny];
                            
                            if (flags_.constraint == PressureGradient) {
                                // dPdx is supplied, put dPdx at both t=j and t=j+1 on RHS
                                Rxk_.re[0] -= a_b*dPdxAct_ + dPdxRef_;
                                
                                // Solve the tau equations
                                tausolver_[j][mx][mz].solve(uk_, vk_, wk_, Pk_, Rxk_, Ryk_, Rzk_);
                                // Bulk velocity is free variable on LHS solved by tau eqn
                                UbulkAct_ = UbulkBase_ + uk_.re.mean();
                                dPdxAct_ = dPdxRef_;
                            }
                            else { // const bulk velocity
                                // Add the previous time-step's -dPdx to the RHS. The next
                                // timestep's dPdx term appears on LHS as unknown.
                                Rxk_.re[0] -= a_b*dPdxAct_;
                                
                                // Use tausolver with additional variable and constraint:
                                // free variable: dPdxAct at next time-step,
                                // constraint:    UbulkBase + mean(u) = UbulkRef.
                                tausolver_[j][mx][mz].solve(uk_, vk_, wk_, Pk_, dPdxAct_,
                                                            Rxk_, Ryk_, Rzk_,
                                                            UbulkRef_ - UbulkBase_);
                                
                                UbulkAct_ = UbulkBase_ + uk_.re.mean(); // should == UbulkRef_
                                //UbulkAct_ = UbulkRef_;
                            }
                            // for kx=kz=0, constant term of pressure is arbitrary 3/19/05
                            // Pk_.set(0, Complex(0.0, 0.0));
                        }
                        
                        // Load solutions back into the external 3d data arrays.
                        // Because of FFTW complex symmetries
                        // The 0,0 mode must be real.
                        // For Nx even, the kxmax,0 mode must be real
                        // For Nz even, the 0,kzmax mode must be real
                        // For Nx,Nz even, the kxmax,kzmax mode must be real
                        if ((kx == 0 && kz == 0) ||
                            (Nx_%2 == 0 && kx == kxmax && kz == 0) ||
                            (Nz_%2 == 0 && kz == kzmax && kx == 0) ||
                            (Nx_%2 == 0 && Nz_%2 == 0 && kx == kxmax && kz == kzmax)) {
                            
                            for (int ny=0; ny<Ny_; ++ny) {
                                un.cmplx(mx,ny,mz,0) = Complex(Re(uk_[ny]), 0.0);
                                un.cmplx(mx,ny,mz,1) = Complex(Re(vk_[ny]), 0.0);
                                un.cmplx(mx,ny,mz,2) = Complex(Re(wk_[ny]), 0.0);
                                qn.cmplx(mx,ny,mz,0) = Complex(Re(Pk_[ny]), 0.0);
                            }
                        }
                        // The normal case, for general kx,kz
                        else
                            for (int ny=0; ny<Ny_; ++ny) {
                                un.cmplx(mx,ny,mz,0) = uk_[ny];
                                un.cmplx(mx,ny,mz,1) = vk_[ny];
                                un.cmplx(mx,ny,mz,2) = wk_[ny];
                                qn.cmplx(mx,ny,mz,0) = Pk_[ny];
                            }
                        
                        // And now set the y-aliased modes to zero.
                        for (int ny=Nyd_; ny<Ny_; ++ny) {
                            un.cmplx(mx,ny,mz,0) = 0.0;
                            un.cmplx(mx,ny,mz,1) = 0.0;
                            un.cmplx(mx,ny,mz,2) = 0.0;
                            qn.cmplx(mx,ny,mz,0) = 0.0;
                        }
                    }
                }
            }
            t_ += dt_;
            if (flags_.verbosity == PrintTime || flags_.verbosity == PrintAll)
                *flags_.logstream << t_ << ' ' << flush;
            else if (flags_.verbosity == PrintTicks)
                *flags_.logstream << '.' << flush;
        }
        
        cfl_ = un.CFLfactor(Ubase_);
        cfl_ *= flags_.dealias_xz() ? 2.0*pi/3.0*dt_ : pi*dt_;
        
        // If using dealiasing, set flag in FlowField that compactifies binary IO
        un.setPadded(flags_.dealias_xz());
        qn.setPadded(flags_.dealias_xz());
        
        if (flags_.verbosity == PrintTime || flags_.verbosity == PrintAll)
            *flags_.logstream << endl;
        
        return;
    }

    
    
    void CNABstyleDNS::advance_inhom(FlowField& un, FlowField& qn, FlowField& BCs, int Nsteps) {
        const int kxmax = tmp_.kxmax();
        const int kzmax = tmp_.kzmax();
        
        for (int n=0; n<Nsteps; ++n) {
            for (int j=0; j<Nsubsteps_; ++j) {
                
                // Store substeps uj,qj in un,qn; reflect this in notation
                FlowField& uj(un);
                FlowField& qj(qn);
                
                swap(fj_, fj1_);
                navierstokesNL(un, ubase_, Ubase_, fj_, tmp_, tmp2_, flags_.nonlinearity);
                
                // Set convenience variables
                Real a_b = alpha_[j]/beta_[j];
                Real g_b = gamma_[j]/beta_[j];
                Real z_b = zeta_[j]/beta_[j];
                Real anu_b = a_b*nu_;
                Real anu = alpha_[j]*nu_;
                
                // Update each Fourier mode with time-stepping algorithm
                for (int mx=0; mx<Mx_; ++mx) {
                    const int kx = uj.kx(mx);
                    
                    for (int mz=0; mz<Mz_; ++mz) {
                        const int kz = uj.kz(mz);
                        
                        Complex u_upper=BCs.cmplx(mx,1,mz,0);
                        Complex u_lower=BCs.cmplx(mx,0,mz,0);
                        Complex v_upper=BCs.cmplx(mx,1,mz,1);
                        Complex v_lower=BCs.cmplx(mx,0,mz,1);
                        Complex w_upper=BCs.cmplx(mx,1,mz,2);
                        Complex w_lower=BCs.cmplx(mx,0,mz,2);
                        
                        Real vr_upper=real(v_upper);
                        Real vi_upper=imag(v_upper);
                        Real vr_lower=real(v_lower);
                        Real vi_lower=imag(v_lower);               
                       
                        
                        // Zero out the aliased modes
                        if ((kx == kxmax || kz == kzmax) ||
                            (flags_.dealias_xz() && isAliasedMode(kx,kz))) {
                            for (int ny=0; ny<Ny_; ++ny) {
                                uj.cmplx(mx,ny,mz,0) = 0.0;
                                uj.cmplx(mx,ny,mz,1) = 0.0;
                                uj.cmplx(mx,ny,mz,2) = 0.0;
                                qj.cmplx(mx,ny,mz,0) = 0.0;
                            }
                            break;
                        }
                        
                        Rxk_.setToZero();
                        Ryk_.setToZero();
                        Rzk_.setToZero();
                        
                        // Goal is to compute
                        // R = a/b nu uj" + [1/(b dt)- a/b nu kappa2] uj - a/b grad qj
                        //     - g/b fj - z/b fj1
                        //
                        //   = a/b nu uj" + [1/(nu a dt)- kappa2] a/b nu uj
                        //     - a/b grad qj - g/b fj - z/b fj1
                        
                        // Extract relevant Fourier modes of uj and qj for computations
                        // set (uk,vk,wk) to a/b nu uj, Pk to a/b qj
                        
                        for (int ny=0; ny<Ny_; ++ny) {
                            uk_.set(ny, anu_b*uj.cmplx(mx,ny,mz,0));
                            vk_.set(ny, anu_b*uj.cmplx(mx,ny,mz,1));
                            wk_.set(ny, anu_b*uj.cmplx(mx,ny,mz,2));
                            Pk_.set(ny, a_b*qj.cmplx(mx,ny,mz,0));
                        }
                        
                        // (1) Put a/b nu uj" into in R. (Pyk_ is used as tmp workspace)
                        diff2(uk_, Rxk_, Pyk_);
                        diff2(vk_, Ryk_, Pyk_);
                        diff2(wk_, Rzk_, Pyk_);
                        
                        // (2) Put a/b qj' into Pyk (compute y-comp of pressure gradient).
                        diff(Pk_, Pyk_);
                        
                        // (3) Add [1/(nu a dt)- kappa2] a/b nu uj - a/b grad qj - g/b fj - z/b fj1
                        // to R, completing calculation of R
                        const Real kappa2 = 4*pi*pi*(square(kx/Lx_) + square(kz/Lz_));
                        const Real c = 1.0/(anu*dt_) - kappa2;
                        const Complex Dx = un.Dx(mx);
                        const Complex Dz = un.Dz(mz);
                        for (int ny=0; ny<Ny_; ++ny) {
                            Rxk_.add(ny, c*uk_[ny] - Dx*Pk_[ny]
                                     - g_b*fj_.cmplx(mx,ny,mz,0) - z_b*fj1_.cmplx(mx,ny,mz,0));
                            Ryk_.add(ny, c*vk_[ny] - Pyk_[ny]
                                     - g_b*fj_.cmplx(mx,ny,mz,1) - z_b*fj1_.cmplx(mx,ny,mz,1));
                            Rzk_.add(ny, c*wk_[ny] - Dz*Pk_[ny]
                                     - g_b*fj_.cmplx(mx,ny,mz,2) - z_b*fj1_.cmplx(mx,ny,mz,2));
                        }
                        
                        // Solve the tau solutions
                        if (kx!=0 || kz!=0)
                            //tausolver_[j][mx][mz].solve(uk_,vk_,wk_,Pk_, Rxk_,Ryk_,Rzk_);
                            tausolver_[j][mx][mz].solve_Inhom(uk_,vk_,wk_,Pk_, Rxk_,Ryk_,Rzk_, vr_lower, vr_upper, vi_lower, vi_upper, u_lower, u_upper, w_lower, w_upper);
                            
                            
                        else { // kx,kz == 0,0
                            // Rx has additional terms, nu Uyy at both t=j and t=j+1
                            Real a_b = alpha_[j]/beta_[j];
                            Real c = nu_*(a_b+1.0);
                            if (Ubaseyy_.length() > 0)
                                for (int ny=0; ny<Ny_; ++ny)
                                    Rxk_.re[ny] += c*Ubaseyy_[ny];
                            
                            if (flags_.constraint == PressureGradient) {
                                // dPdx is supplied, put dPdx at both t=j and t=j+1 on RHS
                                Rxk_.re[0] -= a_b*dPdxAct_ + dPdxRef_;
                                
                                // Solve the tau equations
                                //tausolver_[j][mx][mz].solve(uk_, vk_, wk_, Pk_, Rxk_, Ryk_, Rzk_);
                                tausolver_[j][mx][mz].solve_Inhom(uk_,vk_,wk_,Pk_, Rxk_,Ryk_,Rzk_, vr_lower, vr_upper, vi_lower, vi_upper, u_lower, u_upper, w_lower, w_upper);
                                
                                
                                // Bulk velocity is free variable on LHS solved by tau eqn
                                UbulkAct_ = UbulkBase_ + uk_.re.mean();
                                dPdxAct_ = dPdxRef_;
                            }
                            else { // const bulk velocity
                                // Add the previous time-step's -dPdx to the RHS. The next
                                // timestep's dPdx term appears on LHS as unknown.
                                Rxk_.re[0] -= a_b*dPdxAct_;
                                
                                // Use tausolver with additional variable and constraint:
                                // free variable: dPdxAct at next time-step,
                                // constraint:    UbulkBase + mean(u) = UbulkRef.
                                //tausolver_[j][mx][mz].solve(uk_, vk_, wk_, Pk_, dPdxAct_, Rxk_, Ryk_, Rzk_, UbulkRef_ - UbulkBase_);
                                tausolver_[j][mx][mz].solve_Inhom(uk_, vk_, wk_, Pk_, dPdxAct_,
                                                               Rxk_, Ryk_, Rzk_, vr_lower, vr_upper, vi_lower, vi_upper,                                                   UbulkRef_ - UbulkBase_);
                                
                                UbulkAct_ = UbulkBase_ + uk_.re.mean(); // should == UbulkRef_
                                //UbulkAct_ = UbulkRef_;
                            }
                            // for kx=kz=0, constant term of pressure is arbitrary 3/19/05
                            // Pk_.set(0, Complex(0.0, 0.0));
                        }
                        
                        // Load solutions back into the external 3d data arrays.
                        // Because of FFTW complex symmetries
                        // The 0,0 mode must be real.
                        // For Nx even, the kxmax,0 mode must be real
                        // For Nz even, the 0,kzmax mode must be real
                        // For Nx,Nz even, the kxmax,kzmax mode must be real
                        if ((kx == 0 && kz == 0) ||
                            (Nx_%2 == 0 && kx == kxmax && kz == 0) ||
                            (Nz_%2 == 0 && kz == kzmax && kx == 0) ||
                            (Nx_%2 == 0 && Nz_%2 == 0 && kx == kxmax && kz == kzmax)) {
                            
                            for (int ny=0; ny<Ny_; ++ny) {
                                un.cmplx(mx,ny,mz,0) = Complex(Re(uk_[ny]), 0.0);
                                un.cmplx(mx,ny,mz,1) = Complex(Re(vk_[ny]), 0.0);
                                un.cmplx(mx,ny,mz,2) = Complex(Re(wk_[ny]), 0.0);
                                qn.cmplx(mx,ny,mz,0) = Complex(Re(Pk_[ny]), 0.0);
                            }
                        }
                        // The normal case, for general kx,kz
                        else
                            for (int ny=0; ny<Ny_; ++ny) {
                                un.cmplx(mx,ny,mz,0) = uk_[ny];
                                un.cmplx(mx,ny,mz,1) = vk_[ny];
                                un.cmplx(mx,ny,mz,2) = wk_[ny];
                                qn.cmplx(mx,ny,mz,0) = Pk_[ny];
                            }
                        
                        // And now set the y-aliased modes to zero.
                        for (int ny=Nyd_; ny<Ny_; ++ny) {
                            un.cmplx(mx,ny,mz,0) = 0.0;
                            un.cmplx(mx,ny,mz,1) = 0.0;
                            un.cmplx(mx,ny,mz,2) = 0.0;
                            qn.cmplx(mx,ny,mz,0) = 0.0;
                        }
                    }
                }
            }
            t_ += dt_;
            if (flags_.verbosity == PrintTime || flags_.verbosity == PrintAll)
                *flags_.logstream << t_ << ' ' << flush;
            else if (flags_.verbosity == PrintTicks)
                *flags_.logstream << '.' << flush;
        }
        
        cfl_ = un.CFLfactor(Ubase_);
        cfl_ *= flags_.dealias_xz() ? 2.0*pi/3.0*dt_ : pi*dt_;
        
        // If using dealiasing, set flag in FlowField that compactifies binary IO
        un.setPadded(flags_.dealias_xz());
        qn.setPadded(flags_.dealias_xz());
        
        if (flags_.verbosity == PrintTime || flags_.verbosity == PrintAll)
            *flags_.logstream << endl;
        
        return;
    }
    
    
    
   
    
    
    
    
    
  void CNABstyleDNS::advance_inhom_CON(Controller& controller, FlowField& un, FlowField& qn, FlowField& BCs, FlowField& Fn, int Nsteps, double*** IO, double*** CStateMat) {
        const int kxmax = tmp_.kxmax();
        const int kzmax = tmp_.kzmax();
        

        
        
        for (int n=0; n<Nsteps; ++n) {
            
            
	  controller.advance_Con_CN(un,qn,BCs,CStateMat,dt_,t_,IO);
                        
            
            for (int j=0; j<Nsubsteps_; ++j) {
                
                // Store substeps uj,qj in un,qn; reflect this in notation
                FlowField& uj(un);
                FlowField& qj(qn);
                
                swap(fj_, fj1_);
                navierstokesNL(un, ubase_, Ubase_, fj_, tmp_, tmp2_, flags_.nonlinearity);
                
		Fn = fj_;

                // Set convenience variables
                Real a_b = alpha_[j]/beta_[j];
                Real g_b = gamma_[j]/beta_[j];
                Real z_b = zeta_[j]/beta_[j];
                Real anu_b = a_b*nu_;
                Real anu = alpha_[j]*nu_;
                
                // Update each Fourier mode with time-stepping algorithm
                for (int mx=0; mx<Mx_; ++mx) {
                    const int kx = uj.kx(mx);
                    
                    for (int mz=0; mz<Mz_; ++mz) {
                        const int kz = uj.kz(mz);
                        
                        Complex u_upper=BCs.cmplx(mx,1,mz,0);
                        Complex u_lower=BCs.cmplx(mx,0,mz,0);
                        Complex v_upper=BCs.cmplx(mx,1,mz,1);
                        Complex v_lower=BCs.cmplx(mx,0,mz,1);
                        Complex w_upper=BCs.cmplx(mx,1,mz,2);
                        Complex w_lower=BCs.cmplx(mx,0,mz,2);
                        
                        Real vr_upper=real(v_upper);
                        Real vi_upper=imag(v_upper);
                        Real vr_lower=real(v_lower);
                        Real vi_lower=imag(v_lower);               
                        
                        
                        // Zero out the aliased modes
                        if ((kx == kxmax || kz == kzmax) ||
                            (flags_.dealias_xz() && isAliasedMode(kx,kz))) {
                            for (int ny=0; ny<Ny_; ++ny) {
                                uj.cmplx(mx,ny,mz,0) = 0.0;
                                uj.cmplx(mx,ny,mz,1) = 0.0;
                                uj.cmplx(mx,ny,mz,2) = 0.0;
                                qj.cmplx(mx,ny,mz,0) = 0.0;
                            }
                            break;
                        }
                        
                        Rxk_.setToZero();
                        Ryk_.setToZero();
                        Rzk_.setToZero();
                        
                        // Goal is to compute
                        // R = a/b nu uj" + [1/(b dt)- a/b nu kappa2] uj - a/b grad qj
                        //     - g/b fj - z/b fj1
                        //
                        //   = a/b nu uj" + [1/(nu a dt)- kappa2] a/b nu uj
                        //     - a/b grad qj - g/b fj - z/b fj1
                        
                        // Extract relevant Fourier modes of uj and qj for computations
                        // set (uk,vk,wk) to a/b nu uj, Pk to a/b qj
                        
                        for (int ny=0; ny<Ny_; ++ny) {
                            uk_.set(ny, anu_b*uj.cmplx(mx,ny,mz,0));
                            vk_.set(ny, anu_b*uj.cmplx(mx,ny,mz,1));
                            wk_.set(ny, anu_b*uj.cmplx(mx,ny,mz,2));
                            Pk_.set(ny, a_b*qj.cmplx(mx,ny,mz,0));
                        }
                        
                        // (1) Put a/b nu uj" into in R. (Pyk_ is used as tmp workspace)
                        diff2(uk_, Rxk_, Pyk_);
                        diff2(vk_, Ryk_, Pyk_);
                        diff2(wk_, Rzk_, Pyk_);
                        
                        // (2) Put a/b qj' into Pyk (compute y-comp of pressure gradient).
                        diff(Pk_, Pyk_);
                        
                        // (3) Add [1/(nu a dt)- kappa2] a/b nu uj - a/b grad qj - g/b fj - z/b fj1
                        // to R, completing calculation of R
                        const Real kappa2 = 4*pi*pi*(square(kx/Lx_) + square(kz/Lz_));
                        const Real c = 1.0/(anu*dt_) - kappa2;
                        const Complex Dx = un.Dx(mx);
                        const Complex Dz = un.Dz(mz);
                        for (int ny=0; ny<Ny_; ++ny) {
                            Rxk_.add(ny, c*uk_[ny] - Dx*Pk_[ny]
                                     - g_b*fj_.cmplx(mx,ny,mz,0) - z_b*fj1_.cmplx(mx,ny,mz,0));
                            Ryk_.add(ny, c*vk_[ny] - Pyk_[ny]
                                     - g_b*fj_.cmplx(mx,ny,mz,1) - z_b*fj1_.cmplx(mx,ny,mz,1));
                            Rzk_.add(ny, c*wk_[ny] - Dz*Pk_[ny]
                                     - g_b*fj_.cmplx(mx,ny,mz,2) - z_b*fj1_.cmplx(mx,ny,mz,2));
                        }
                        
                        // Solve the tau solutions
                        if (kx!=0 || kz!=0)
                            //tausolver_[j][mx][mz].solve(uk_,vk_,wk_,Pk_, Rxk_,Ryk_,Rzk_);
                            tausolver_[j][mx][mz].solve_Inhom(uk_,vk_,wk_,Pk_, Rxk_,Ryk_,Rzk_, vr_lower, vr_upper, vi_lower, vi_upper, u_lower, u_upper, w_lower, w_upper);
                        
                        
                        else { // kx,kz == 0,0
                            // Rx has additional terms, nu Uyy at both t=j and t=j+1
                            Real a_b = alpha_[j]/beta_[j];
                            Real c = nu_*(a_b+1.0);
                            if (Ubaseyy_.length() > 0)
                                for (int ny=0; ny<Ny_; ++ny)
                                    Rxk_.re[ny] += c*Ubaseyy_[ny];
                            
                            if (flags_.constraint == PressureGradient) {
                                // dPdx is supplied, put dPdx at both t=j and t=j+1 on RHS
                                Rxk_.re[0] -= a_b*dPdxAct_ + dPdxRef_;
                                
                                // Solve the tau equations
                                //tausolver_[j][mx][mz].solve(uk_, vk_, wk_, Pk_, Rxk_, Ryk_, Rzk_);
                                tausolver_[j][mx][mz].solve_Inhom(uk_,vk_,wk_,Pk_, Rxk_,Ryk_,Rzk_, vr_lower, vr_upper, vi_lower, vi_upper, u_lower, u_upper, w_lower, w_upper);
                                
                                
                                // Bulk velocity is free variable on LHS solved by tau eqn
                                UbulkAct_ = UbulkBase_ + uk_.re.mean();
                                dPdxAct_ = dPdxRef_;
                            }
                            else { // const bulk velocity
                                // Add the previous time-step's -dPdx to the RHS. The next
                                // timestep's dPdx term appears on LHS as unknown.
                                Rxk_.re[0] -= a_b*dPdxAct_;
                                
                                // Use tausolver with additional variable and constraint:
                                // free variable: dPdxAct at next time-step,
                                // constraint:    UbulkBase + mean(u) = UbulkRef.
                                //tausolver_[j][mx][mz].solve(uk_, vk_, wk_, Pk_, dPdxAct_, Rxk_, Ryk_, Rzk_, UbulkRef_ - UbulkBase_);
                                tausolver_[j][mx][mz].solve_Inhom(uk_, vk_, wk_, Pk_, dPdxAct_,
                                                                  Rxk_, Ryk_, Rzk_, vr_lower, vr_upper, vi_lower, vi_upper,                                                   UbulkRef_ - UbulkBase_);
                                
                                UbulkAct_ = UbulkBase_ + uk_.re.mean(); // should == UbulkRef_
                                //UbulkAct_ = UbulkRef_;
                            }
                            // for kx=kz=0, constant term of pressure is arbitrary 3/19/05
                            // Pk_.set(0, Complex(0.0, 0.0));
                        }
                        
                        // Load solutions back into the external 3d data arrays.
                        // Because of FFTW complex symmetries
                        // The 0,0 mode must be real.
                        // For Nx even, the kxmax,0 mode must be real
                        // For Nz even, the 0,kzmax mode must be real
                        // For Nx,Nz even, the kxmax,kzmax mode must be real
                        if ((kx == 0 && kz == 0) ||
                            (Nx_%2 == 0 && kx == kxmax && kz == 0) ||
                            (Nz_%2 == 0 && kz == kzmax && kx == 0) ||
                            (Nx_%2 == 0 && Nz_%2 == 0 && kx == kxmax && kz == kzmax)) {
                            
                            for (int ny=0; ny<Ny_; ++ny) {
                                un.cmplx(mx,ny,mz,0) = Complex(Re(uk_[ny]), 0.0);
                                un.cmplx(mx,ny,mz,1) = Complex(Re(vk_[ny]), 0.0);
                                un.cmplx(mx,ny,mz,2) = Complex(Re(wk_[ny]), 0.0);
                                qn.cmplx(mx,ny,mz,0) = Complex(Re(Pk_[ny]), 0.0);
                            }
                        }
                        // The normal case, for general kx,kz
                        else
                            for (int ny=0; ny<Ny_; ++ny) {
                                un.cmplx(mx,ny,mz,0) = uk_[ny];
                                un.cmplx(mx,ny,mz,1) = vk_[ny];
                                un.cmplx(mx,ny,mz,2) = wk_[ny];
                                qn.cmplx(mx,ny,mz,0) = Pk_[ny];
                            }
                        
                        // And now set the y-aliased modes to zero.
                        for (int ny=Nyd_; ny<Ny_; ++ny) {
                            un.cmplx(mx,ny,mz,0) = 0.0;
                            un.cmplx(mx,ny,mz,1) = 0.0;
                            un.cmplx(mx,ny,mz,2) = 0.0;
                            qn.cmplx(mx,ny,mz,0) = 0.0;
                        }
                    }
                }
            }
            t_ += dt_;
            if (flags_.verbosity == PrintTime || flags_.verbosity == PrintAll)
                *flags_.logstream << t_ << ' ' << flush;
            else if (flags_.verbosity == PrintTicks)
                *flags_.logstream << '.' << flush;
            
            
            //controller.advance_Con_CN(un,qn,BCs,CStateMat,dt_,t_,IO);
	    cfl_ = un.CFLfactor(Ubase_);
            cfl_ *= flags_.dealias_xz() ? 2.0*pi/3.0*dt_ : pi*dt_;
            
            if (cfl_ > 1.0){
                cout << "Error: CFL>1" << endl;
                break;
            }

            
        }
        
        cfl_ = un.CFLfactor(Ubase_);
        cfl_ *= flags_.dealias_xz() ? 2.0*pi/3.0*dt_ : pi*dt_;
        
        // If using dealiasing, set flag in FlowField that compactifies binary IO
        un.setPadded(flags_.dealias_xz());
        qn.setPadded(flags_.dealias_xz());
        
        if (flags_.verbosity == PrintTime || flags_.verbosity == PrintAll)
            *flags_.logstream << endl;
        
        return;
    }

/*
  void CNABstyleDNS::advance_inhom_CON_SF(Controller& controller, FlowField& un, FlowField& qn, FlowField& BCs, FlowField& Fn, int Nsteps, double*** IO) {
        const int kxmax = tmp_.kxmax();
        const int kzmax = tmp_.kzmax();
        

        
        
        for (int n=0; n<Nsteps; ++n) {
            
            
            controller.advance_Con_SF(un,qn,BCs,dt_,t_,IO);
                        
            
            for (int j=0; j<Nsubsteps_; ++j) {
                
                // Store substeps uj,qj in un,qn; reflect this in notation
                FlowField& uj(un);
                FlowField& qj(qn);
                
                swap(fj_, fj1_);
                navierstokesNL(un, ubase_, Ubase_, fj_, tmp_, tmp2_, flags_.nonlinearity);
                
		Fn = fj_;

                // Set convenience variables
                Real a_b = alpha_[j]/beta_[j];
                Real g_b = gamma_[j]/beta_[j];
                Real z_b = zeta_[j]/beta_[j];
                Real anu_b = a_b*nu_;
                Real anu = alpha_[j]*nu_;
                
                // Update each Fourier mode with time-stepping algorithm
                for (int mx=0; mx<Mx_; ++mx) {
                    const int kx = uj.kx(mx);
                    
                    for (int mz=0; mz<Mz_; ++mz) {
                        const int kz = uj.kz(mz);
                        
                        Complex u_upper=BCs.cmplx(mx,1,mz,0);
                        Complex u_lower=BCs.cmplx(mx,0,mz,0);
                        Complex v_upper=BCs.cmplx(mx,1,mz,1);
                        Complex v_lower=BCs.cmplx(mx,0,mz,1);
                        Complex w_upper=BCs.cmplx(mx,1,mz,2);
                        Complex w_lower=BCs.cmplx(mx,0,mz,2);
                        
                        Real vr_upper=real(v_upper);
                        Real vi_upper=imag(v_upper);
                        Real vr_lower=real(v_lower);
                        Real vi_lower=imag(v_lower);               
                        
                        
                        // Zero out the aliased modes
                        if ((kx == kxmax || kz == kzmax) ||
                            (flags_.dealias_xz() && isAliasedMode(kx,kz))) {
                            for (int ny=0; ny<Ny_; ++ny) {
                                uj.cmplx(mx,ny,mz,0) = 0.0;
                                uj.cmplx(mx,ny,mz,1) = 0.0;
                                uj.cmplx(mx,ny,mz,2) = 0.0;
                                qj.cmplx(mx,ny,mz,0) = 0.0;
                            }
                            break;
                        }
                        
                        Rxk_.setToZero();
                        Ryk_.setToZero();
                        Rzk_.setToZero();
                        
                        // Goal is to compute
                        // R = a/b nu uj" + [1/(b dt)- a/b nu kappa2] uj - a/b grad qj
                        //     - g/b fj - z/b fj1
                        //
                        //   = a/b nu uj" + [1/(nu a dt)- kappa2] a/b nu uj
                        //     - a/b grad qj - g/b fj - z/b fj1
                        
                        // Extract relevant Fourier modes of uj and qj for computations
                        // set (uk,vk,wk) to a/b nu uj, Pk to a/b qj
                        
                        for (int ny=0; ny<Ny_; ++ny) {
                            uk_.set(ny, anu_b*uj.cmplx(mx,ny,mz,0));
                            vk_.set(ny, anu_b*uj.cmplx(mx,ny,mz,1));
                            wk_.set(ny, anu_b*uj.cmplx(mx,ny,mz,2));
                            Pk_.set(ny, a_b*qj.cmplx(mx,ny,mz,0));
                        }
                        
                        // (1) Put a/b nu uj" into in R. (Pyk_ is used as tmp workspace)
                        diff2(uk_, Rxk_, Pyk_);
                        diff2(vk_, Ryk_, Pyk_);
                        diff2(wk_, Rzk_, Pyk_);
                        
                        // (2) Put a/b qj' into Pyk (compute y-comp of pressure gradient).
                        diff(Pk_, Pyk_);
                        
                        // (3) Add [1/(nu a dt)- kappa2] a/b nu uj - a/b grad qj - g/b fj - z/b fj1
                        // to R, completing calculation of R
                        const Real kappa2 = 4*pi*pi*(square(kx/Lx_) + square(kz/Lz_));
                        const Real c = 1.0/(anu*dt_) - kappa2;
                        const Complex Dx = un.Dx(mx);
                        const Complex Dz = un.Dz(mz);
                        for (int ny=0; ny<Ny_; ++ny) {
                            Rxk_.add(ny, c*uk_[ny] - Dx*Pk_[ny]
                                     - g_b*fj_.cmplx(mx,ny,mz,0) - z_b*fj1_.cmplx(mx,ny,mz,0));
                            Ryk_.add(ny, c*vk_[ny] - Pyk_[ny]
                                     - g_b*fj_.cmplx(mx,ny,mz,1) - z_b*fj1_.cmplx(mx,ny,mz,1));
                            Rzk_.add(ny, c*wk_[ny] - Dz*Pk_[ny]
                                     - g_b*fj_.cmplx(mx,ny,mz,2) - z_b*fj1_.cmplx(mx,ny,mz,2));
                        }
                        
                        // Solve the tau solutions
                        if (kx!=0 || kz!=0)
                            //tausolver_[j][mx][mz].solve(uk_,vk_,wk_,Pk_, Rxk_,Ryk_,Rzk_);
                            tausolver_[j][mx][mz].solve_Inhom(uk_,vk_,wk_,Pk_, Rxk_,Ryk_,Rzk_, vr_lower, vr_upper, vi_lower, vi_upper, u_lower, u_upper, w_lower, w_upper);
                        
                        
                        else { // kx,kz == 0,0
                            // Rx has additional terms, nu Uyy at both t=j and t=j+1
                            Real a_b = alpha_[j]/beta_[j];
                            Real c = nu_*(a_b+1.0);
                            if (Ubaseyy_.length() > 0)
                                for (int ny=0; ny<Ny_; ++ny)
                                    Rxk_.re[ny] += c*Ubaseyy_[ny];
                            
                            if (flags_.constraint == PressureGradient) {
                                // dPdx is supplied, put dPdx at both t=j and t=j+1 on RHS
                                Rxk_.re[0] -= a_b*dPdxAct_ + dPdxRef_;
                                
                                // Solve the tau equations
                                //tausolver_[j][mx][mz].solve(uk_, vk_, wk_, Pk_, Rxk_, Ryk_, Rzk_);
                                tausolver_[j][mx][mz].solve_Inhom(uk_,vk_,wk_,Pk_, Rxk_,Ryk_,Rzk_, vr_lower, vr_upper, vi_lower, vi_upper, u_lower, u_upper, w_lower, w_upper);
                                
                                
                                // Bulk velocity is free variable on LHS solved by tau eqn
                                UbulkAct_ = UbulkBase_ + uk_.re.mean();
                                dPdxAct_ = dPdxRef_;
                            }
                            else { // const bulk velocity
                                // Add the previous time-step's -dPdx to the RHS. The next
                                // timestep's dPdx term appears on LHS as unknown.
                                Rxk_.re[0] -= a_b*dPdxAct_;
                                
                                // Use tausolver with additional variable and constraint:
                                // free variable: dPdxAct at next time-step,
                                // constraint:    UbulkBase + mean(u) = UbulkRef.
                                //tausolver_[j][mx][mz].solve(uk_, vk_, wk_, Pk_, dPdxAct_, Rxk_, Ryk_, Rzk_, UbulkRef_ - UbulkBase_);
                                tausolver_[j][mx][mz].solve_Inhom(uk_, vk_, wk_, Pk_, dPdxAct_,
                                                                  Rxk_, Ryk_, Rzk_, vr_lower, vr_upper, vi_lower, vi_upper,                                                   UbulkRef_ - UbulkBase_);
                                
                                UbulkAct_ = UbulkBase_ + uk_.re.mean(); // should == UbulkRef_
                                //UbulkAct_ = UbulkRef_;
                            }
                            // for kx=kz=0, constant term of pressure is arbitrary 3/19/05
                            // Pk_.set(0, Complex(0.0, 0.0));
                        }
                        
                        // Load solutions back into the external 3d data arrays.
                        // Because of FFTW complex symmetries
                        // The 0,0 mode must be real.
                        // For Nx even, the kxmax,0 mode must be real
                        // For Nz even, the 0,kzmax mode must be real
                        // For Nx,Nz even, the kxmax,kzmax mode must be real
                        if ((kx == 0 && kz == 0) ||
                            (Nx_%2 == 0 && kx == kxmax && kz == 0) ||
                            (Nz_%2 == 0 && kz == kzmax && kx == 0) ||
                            (Nx_%2 == 0 && Nz_%2 == 0 && kx == kxmax && kz == kzmax)) {
                            
                            for (int ny=0; ny<Ny_; ++ny) {
                                un.cmplx(mx,ny,mz,0) = Complex(Re(uk_[ny]), 0.0);
                                un.cmplx(mx,ny,mz,1) = Complex(Re(vk_[ny]), 0.0);
                                un.cmplx(mx,ny,mz,2) = Complex(Re(wk_[ny]), 0.0);
                                qn.cmplx(mx,ny,mz,0) = Complex(Re(Pk_[ny]), 0.0);
                            }
                        }
                        // The normal case, for general kx,kz
                        else
                            for (int ny=0; ny<Ny_; ++ny) {
                                un.cmplx(mx,ny,mz,0) = uk_[ny];
                                un.cmplx(mx,ny,mz,1) = vk_[ny];
                                un.cmplx(mx,ny,mz,2) = wk_[ny];
                                qn.cmplx(mx,ny,mz,0) = Pk_[ny];
                            }
                        
                        // And now set the y-aliased modes to zero.
                        for (int ny=Nyd_; ny<Ny_; ++ny) {
                            un.cmplx(mx,ny,mz,0) = 0.0;
                            un.cmplx(mx,ny,mz,1) = 0.0;
                            un.cmplx(mx,ny,mz,2) = 0.0;
                            qn.cmplx(mx,ny,mz,0) = 0.0;
                        }
                    }
                }
            }
            t_ += dt_;
            if (flags_.verbosity == PrintTime || flags_.verbosity == PrintAll)
                *flags_.logstream << t_ << ' ' << flush;
            else if (flags_.verbosity == PrintTicks)
                *flags_.logstream << '.' << flush;
            
            
            //controller.advance_Con_CN(un,qn,BCs,CStateMat,dt_,t_,IO);
	    cfl_ = un.CFLfactor(Ubase_);
            cfl_ *= flags_.dealias_xz() ? 2.0*pi/3.0*dt_ : pi*dt_;
            
            if (cfl_ > 1.0){
                cout << "Error: CFL>1" << endl;
                break;
            }

            
        }
        
        cfl_ = un.CFLfactor(Ubase_);
        cfl_ *= flags_.dealias_xz() ? 2.0*pi/3.0*dt_ : pi*dt_;
        
        // If using dealiasing, set flag in FlowField that compactifies binary IO
        un.setPadded(flags_.dealias_xz());
        qn.setPadded(flags_.dealias_xz());
        
        if (flags_.verbosity == PrintTime || flags_.verbosity == PrintAll)
            *flags_.logstream << endl;
        
        return;
    }

    
    */
    
    
    /**************************************************************************
     
     ---End: Peter H Heins---
     
     **************************************************************************/ 
    
    
    
    
    
    
    
    void changeBaseFlow(const ChebyCoeff& ubase0, const FlowField& ufluc0,
                        const FlowField& q0arg,
                        const ChebyCoeff& ubase1, FlowField& u1, FlowField& q1){
        ChebyCoeff& U0 = (ChebyCoeff&) ubase0;
        fieldstate U0state = U0.state();
        
        ChebyCoeff& U1 = (ChebyCoeff&) ubase1;
        fieldstate U1state = U1.state();
        
        FlowField& u0 = (FlowField&) ufluc0;
        fieldstate u0xzstate = u0.xzstate();
        fieldstate u0ystate = u0.ystate();
        
        FlowField& q0 = (FlowField&) q0arg;
        fieldstate q0xzstate = q0.xzstate();
        fieldstate q0ystate = q0.ystate();
        
        int Nx=u0.numXgridpts();
        int Ny=u0.numYgridpts();
        int Nz=u0.numZgridpts();
        
        u1 = u0; // want u1 FPF
        u1.makeState(Spectral, Physical);
        u0.makePhysical();
        q0.makePhysical();
        q1 = q0; // want q1 physical
        
        // At this point
        // u1 == utot - U0
        // q1 == p + 1/2 u0 dot u0
        
        // Remove 1/2 u0 dot u0 from q1
        for (int ny=0; ny<Ny; ++ny)
            for (int nx=0; nx<Nx; ++nx)
                for (int nz=0; nz<Nz; ++nz)
                    q1(nx,ny,nz,0) -= 0.5*(square(u0(nx,ny,nz,0)) +
                                           square(u0(nx,ny,nz,1)) +
                                           square(u0(nx,ny,nz,2)));
        // At this point
        // u1 == utot - U0
        // q1 == p
        
        ChebyTransform t(U0.numModes());
        U0.makePhysical(t);
        U1.makePhysical(t);
        
        // Add U0-U1 to u1
        ChebyCoeff delta_U(U0);
        delta_U -= U1;
        u1 += delta_U;
        u1.makePhysical();
        
        // At this point
        // u1 == utot - U1
        // q1 == p
        
        // Add 1/2 u1 dot u1 to q1
        for (int ny=0; ny<Ny; ++ny)
            for (int nx=0; nx<Nx; ++nx)
                for (int nz=0; nz<Nz; ++nz)
                    q1(nx,ny,nz,0) += 0.5*(square(u1(nx,ny,nz,0)) +
                                           square(u1(nx,ny,nz,1)) +
                                           square(u1(nx,ny,nz,2)));
        // At this point
        // u1 == utot - U1
        // q1 == p + 1/2 u1 dot u1
        // et, voila
        
        U0.makeState(U0state,t);
        U1.makeState(U1state,t);
        u0.makeState(u0xzstate, u0ystate);
        q0.makeState(q0xzstate, q0ystate);
        u1.makeState(u0xzstate, u0ystate);
        q1.makeState(q0xzstate, q0ystate);
    }
    
    DNSFlags::DNSFlags(BaseFlow        baseflow_,
                       MeanConstraint  constraint_,
                       TimeStepMethod  timestepping_,
                       TimeStepMethod  initstepping_,
                       NonlinearMethod nonlinearity_,
                       Dealiasing      dealiasing_,
                       bool            taucorrection_,
                       Real            dPdx_,
                       Real            Ubulk_,
                       Verbosity       verbosity_,
                       ostream*        logstream_)
    :
    baseflow(baseflow_),
    constraint(constraint_),
    timestepping(timestepping_),
    initstepping(initstepping_),
    nonlinearity(nonlinearity_),
    dealiasing(dealiasing_),
    taucorrection(taucorrection_),
    dPdx(dPdx_),
    Ubulk(Ubulk_),
    verbosity(verbosity_),
    logstream(logstream_)
    {
        if (dealias_y() && (nonlinearity != Rotational)) {
            cerr << "DNSFlags::DNSFlags: DealiasY and DealiasXYZ work only with\n";
            cerr << "Rotational nonlinearity in the current version of channelflow.\n";
            cerr << "Setting nonlinearity to Rotational." << endl;
            nonlinearity = Rotational;
        }
        // maybe should print warnings about initstepping only mattering for SBDF and CNAB
    }
    
    bool DNSFlags::dealias_xz() const {
        return ((dealiasing == DealiasXZ || dealiasing == DealiasXYZ) ? true:false);
    }
    
    bool DNSFlags::dealias_y() const {
        return ((dealiasing == DealiasY || dealiasing == DealiasXYZ) ? true:false);
    }
    
    ostream& operator<<(ostream& os, Dealiasing d) {
        string s;
        switch(d) {
            case NoDealiasing: s="NoDealiasing"; break;
            case DealiasXZ: s="DealiasXZ"; break;
            case DealiasY: s="DealiasY"; break;
            case DealiasXYZ: s="DealiasXYZ"; break;
            default: s="Invalid Dealiasing value: please submit bug report";
        }
        os << s;
        return os;
    }
    
    ostream& operator<<(ostream& os, MeanConstraint m) {
        string s;
        switch(m) {
            case PressureGradient: s="PressureGradient"; break;
            case BulkVelocity: s="BulkVelocity"; break;
            default: s="Invalid MeanConstraint value: please submit bug report";
        }
        os << s;
        return os;
    }
    
    ostream& operator<<(ostream& os, TimeStepMethod t) {
        string s;
        switch(t) {
            case CNFE1: s="CNFE1"; break;
            case CNAB2: s="CNAB2"; break;
            case CNRK2: s="CNRK2"; break;
            case SMRK2: s="SMRK2"; break;
            case SBDF1: s="SBDF1"; break;
            case SBDF2: s="SBDF2"; break;
            case SBDF3: s="SBDF3"; break;
            case SBDF4: s="SBDF4"; break;
            default: s="Invalid TimeStepMethod value: please submit bug report";
        }
        os << s;
        return os;
    }
    ostream& operator<<(ostream& os, Verbosity v) {
        string s;
        switch(v) {
            case Silent:     s="Silent"; break;
            case PrintTicks: s="PrintTicks"; break;
            case PrintTime:  s="PrintTime"; break;
            case VerifyTauSolve: s="VerifyTauSolve"; break;
            case PrintAll:   s="PrintAll"; break;
            default:         s="Invalid Verbosity value: please submit bug report";
        }
        os << s;
        return os;
    }
    
    ostream& operator<<(ostream& os, const DNSFlags& flags) {
        string s(", ");
        string tau = (flags.taucorrection) ? "TauCorrection" : "NoTauCorrection";
        os << flags.baseflow << s
        << flags.constraint << s
        << flags.timestepping << s
        << flags.initstepping << s
        << flags.nonlinearity << s
        << flags.dealiasing << s
        << tau << s
        << flags.verbosity << s
        << "dPdx=="<<flags.dPdx << s
        << "Ubulk=="<<flags.Ubulk;
        return os;
    }
    
    ostream& operator<<(ostream& os, NonlinearMethod nonl) {
        string s;
        switch(nonl) {
            case Rotational: s="Rotational"; break;
            case Convection: s="Convection"; break;
            case Divergence: s="Divergence"; break;
            case SkewSymmetric: s="SkewSymmetric"; break;
            case Alternating: s="Alternating"; break;
            case Alternating_: s="Alternating_"; break;
            case LinearAboutProfile: s="LinearAboutProfile"; break;
            case LinearAboutField: s="LinearAboutField"; break;
            default: s = "Invalid NonlinearMethod: please submit bug report";
        }
        os << s;
        return os;
    }
    
    ostream& operator<<(ostream& os, BaseFlow base) {
        string s;
        switch(base) {
            case Zero:          s="Zero"; break;
            case PlaneCouette:  s="PlaneCouette"; break;
            case Parabolic:     s="Parabolic"; break;
            default: s = "Invalid BaseFlow: please submit bug report";
        }
        os << s;
        return os;
    }
    
    // Make a lowercase copy of s:
    string lowercase(const string& s) {
        char* buf = new char[s.length()];
        s.copy(buf, s.length());
        for(uint i = 0; i < s.length(); i++)
            buf[i] = tolower(buf[i]);
        string r(buf, s.length());
        delete[] buf;
        return r;
    }
    
    TimeStepMethod  s2stepmethod(const string& s_) {
        TimeStepMethod step = SBDF3;
        string s = lowercase(s_);
        if (s.find("cnfe1") != string::npos) step = CNFE1;
        else if (s.find("cnab2") != string::npos) step = CNAB2;
        else if (s.find("cnrk2") != string::npos) step = CNRK2;
        else if (s.find("smrk2") != string::npos) step = SMRK2;
        else if (s.find("sbdf1") != string::npos) step = SBDF1;
        else if (s.find("sbdf2") != string::npos) step = SBDF2;
        else if (s.find("sbdf3") != string::npos) step = SBDF3;
        else if (s.find("sbdf4") != string::npos) step = SBDF4;
        else {
            cerr << "warning : s2stepstepod(string) : unrecognized string " <<s<<endl;
            exit(1);
        }
        return step;
    }
    
    NonlinearMethod s2nonlmethod(const string& s_) {
        NonlinearMethod nonl = Rotational;
        string s = lowercase(s_);
        if (s.find("rot") != string::npos) nonl = Rotational;
        else if (s.find("conv") != string::npos) nonl = Convection;
        else if (s.find("skew") != string::npos) nonl = SkewSymmetric;
        else if (s.find("alt")  != string::npos) nonl = Alternating;
        else if (s.find("div")  != string::npos) nonl = Divergence;
        else if (s.find("field") != string::npos) nonl = LinearAboutField;
        else if (s.find("profile") != string::npos) nonl = LinearAboutProfile;
        else {
            cerr << "warning : s2nonlmethod(string) : unrecognized string " <<s<<endl;
            exit(1);
        }
        return nonl;
    }
    
    Dealiasing s2dealaising(const std::string& s_) {
        Dealiasing d = DealiasXZ;
        string s = lowercase(s_);
        if (s.find("nodealiasing")    != string::npos) d = NoDealiasing;
        else if (s.find("dealiasxz")  != string::npos) d = DealiasXZ;
        else if (s.find("dealiasy")   != string::npos) d = DealiasY;
        else if (s.find("dealiasxyz") != string::npos) d = DealiasXYZ;
        else {
            cerr << "warning : s2dealiasing(string) : unrecognized string " <<s<<endl;
            exit(1);
        }
        return d;
    }
    
    MeanConstraint s2constraint(const std::string& s_) {
        MeanConstraint m;
        string s = lowercase(s_);
        if (s.find("pressure")    != string::npos) m = PressureGradient;
        else if (s.find("velocity")  != string::npos) m = BulkVelocity;
        else {
            cerr << "warning : s2constraint(string) : unrecognized string " <<s<<endl;
            exit(1);
        }
        return m;
    }
    
    BaseFlow s2baseflow(const std::string& s_) {
        BaseFlow b = Zero;
        string s = lowercase(s_);
        if (s.find("none") != string::npos ||
            s.find("zero") != string::npos)
            b = Zero;
        else if (s.find("pcf")     != string::npos ||
                 s.find("couette") != string::npos ||
                 s.find("linear") != string::npos)
            b = PlaneCouette;
        else if (s.find("parab") != string::npos ||
                 s.find("poiss") != string::npos)
            b = Parabolic;
        else {
            cerr << "warning : s2baseflow(string) : unrecognized string " <<s<<endl;
            exit(1);
        }
        return b;
    }
    
    void uUq2p(const FlowField& u_, const ChebyCoeff& U_, const FlowField& q,
               FlowField& p, NonlinearMethod nonl_method) {
        
    }
    
    void uUp2q(const FlowField& u_, const ChebyCoeff& U_, const FlowField& p,
               FlowField& q, NonlinearMethod nonl_method) {
        
        assert(u_.Nd() == 3);
        assert(q.Nd() == 1);
        assert(p.Nd() == 1);
        
        q = p;
        if (nonl_method != Rotational)
            return;
        
        FlowField& u = const_cast<FlowField&>(u_);
        ChebyCoeff& U = const_cast<ChebyCoeff&>(U_);
        fieldstate uxzstate = u.xzstate();
        fieldstate uystate = u.ystate();
        fieldstate Ustate = U.state();
        
        u.makePhysical();
        U.makePhysical();
        q.makePhysical();
        
        int Nx=u_.Nx();
        int Ny=u_.Ny();
        int Nz=u_.Nz();
        
        // Add 1/2 u dot u to q
        for (int ny=0; ny<Ny; ++ny) {
            Real Uny = U[ny];
            for (int nx=0; nx<Nx; ++nx)
                for (int nz=0; nz<Nz; ++nz)
                    q(nx,ny,nz,0) += 0.5*(square(u(nx,ny,nz,0) + Uny) +
                                          square(u(nx,ny,nz,1)) +
                                          square(u(nx,ny,nz,2)));
        }
        u.makeState(uxzstate, uystate);
        U.makeState(Ustate);
        q.makeState(p.xzstate(), p.ystate());
        return;
    }
    
    
    // *******************************************************************************************
    // BEGIN EXPERIMENTAL CODE
    
    
    PoincareCondition::PoincareCondition() {}
    
    PlaneIntersection::~PlaneIntersection() {}
    PlaneIntersection::PlaneIntersection() {}
    PlaneIntersection::PlaneIntersection(const FlowField& ustar, const FlowField& estar)
    :
    estar_(estar),
    cstar_(L2IP(ustar, estar))
    {}
    
    Real PlaneIntersection::operator()(const FlowField& u) {
        return L2IP(u, estar_) - cstar_;
    }
    
    DragDissipation::DragDissipation() {}
    
    Real DragDissipation::operator()(const FlowField& u) {
        return forcing(u) - dissipation(u);
    };
    
    DNSPoincare::DNSPoincare()
    :
    DNS(),
    e_(),
    sigma_(),
    h_(0),
    ucrossing_(),
    pcrossing_(),
    tcrossing_(0),
    scrossing_(0),
    hcrossing_(0.0),
    hcurrent_(0.0),
    t0_(0)
    {}
    
    
    DNSPoincare::
    DNSPoincare(const FlowField& u,  PoincareCondition* h, Real nu, Real dt, const DNSFlags& flags, Real t0)
    :
    DNS(u, nu, dt, flags, t0),
    e_(),
    sigma_(),
    h_(h),
    ucrossing_(),
    pcrossing_(),
    tcrossing_(0),
    scrossing_(0),
    hcrossing_(0.0),
    hcurrent_((*h)(u)),
    t0_(t0)
    {
    }
    
    DNSPoincare::
    DNSPoincare(const FlowField& u,  const array<FlowField>& e, const array<FieldSymmetry>& sigma,
                PoincareCondition* h, Real nu, Real dt, const DNSFlags& flags, Real t0)
    :
    DNS(u, nu, dt, flags, t0),
    e_(e),
    sigma_(sigma),
    h_(h),
    ucrossing_(),
    pcrossing_(),
    tcrossing_(0),
    scrossing_(0),
    hcrossing_(0.0),
    hcurrent_((*h)(u)),
    t0_(t0)
    {
        // check that sigma[n] e[n] = -e[n]
        //*flags_.logstream << "Checking symmetries and basis sets for fundamental domain " << endl;
        //FlowField tmp;
        //*flags_.logstream << "n \t L2Norm(e[n] + s[n] e[n])/L2Norm(e[n])" << endl;
        //for (int n=0; n<e.length(); ++n) {
        //tmp = e_[n];
        //tmp += sigma_[n](e_[n]);
        //*flags_.logstream << n << '\t' << L2Norm(tmp)/L2Norm(e_[n]) << endl;
        //}
    }
    
    const FlowField& DNSPoincare::ucrossing() const {return ucrossing_;}
    const FlowField& DNSPoincare::pcrossing() const {return pcrossing_;}
    Real DNSPoincare::tcrossing() const {return tcrossing_;}
    int  DNSPoincare::scrossing() const {return scrossing_;}
    Real DNSPoincare::hcrossing() const {return hcrossing_;}
    Real DNSPoincare::hcurrent() const {return hcurrent_;}
    //Real DNSPoincare::f(const FlowField& u) const {
    //  return L2IP(u, estar_) - cstar_;
    //}
    
    bool DNSPoincare::advanceToSection(FlowField& u, FlowField& p, int nSteps, int crosssign, Real Tmin, Real epsilon) {
        
        ostream* os = flags().logstream;
        FlowField uprev(u);
        FlowField pprev(p);
        
        // Take nSteps of length dt, advancing t -> t + nSteps dt
        advance(u, p, nSteps);
        
        // Check for u(t) cross of fundamental domain boundary, map back in
        // Map uprev, pprev back, too, so that
        FieldSymmetry s, identity; // both identity at this point
        for (int n=0; n<e_.length(); ++n) {
            if (L2IP(u, e_[n]) <0) {
                // *os << "Crossed " << n << "th boundary of fundamental domain" << endl;
                s  *= sigma_[n];
            }
        }
        if (s != identity) {
            // *os << "Mapping fields back into fund. domain w s = " << s << endl;
            u *= s;
            p *= s;
            uprev *= s;
            pprev *= s;
            (*this) *= s; // maps all FlowField members of DNS
        }
        
        Real tcoarse = DNS::time();
        if (tcoarse-t0_ > Tmin) {
            //Real cstar = L2IP(ustar_, estar_);
            Real hprev = (*h_)(uprev);
            Real hcurr = (*h_)(u);
            //hcrossing_ = hcurr;
            hcurrent_  = hcurr;
            //*os << tcoarse << '\t' << c << endl;
            
            bool dhdt_pos_cross  =  (hprev<0 && 0<=hcurr) ? true : false;
            bool dhdt_neg_cross  =  (hprev>0 && 0>=hcurr) ? true : false;
            
            // If we cross the Poincare section in required direction, back up and
            // reintegrate, checking condition every integration time step dt.
            if ((crosssign > 0  && dhdt_pos_cross) ||
                (crosssign < 0  && dhdt_neg_cross) ||
                (crosssign == 0 && (dhdt_pos_cross || dhdt_neg_cross))) {
                
                //*os << "u(t) crossed Poincare section..." << endl;
                *os << (dhdt_pos_cross ? '+' : '-') << flush;
                scrossing_ = dhdt_pos_cross ? +1 : -1;
                
                // Allocate arrays to store time, velocity, and Poincare conditions
                // at three consecutive time steps for quadratic interpolation. E.g.
                // v[n],q[n],f[n] are u,p,f at three successive fine-scale timesteps
                array<Real> s(3);      // s[n] = tcoarse - dT + n dt    time-like variable
                array<Real> h(3);      // h[n] = h(s[n])                space-like variable
                array<FlowField> v(3); // v[n] = u(s[n])                velocity field
                array<FlowField> q(3); // q[n] = p(s[n])                pressure field
                
                Real dt = DNS::dt();
                Real dT = dt*nSteps;
                
                s[0] = tcoarse - dT;
                s[1] = 0.0;
                s[2] = 0.0;
                //s[3] = 0.0;
                h[0] = hprev;
                h[1] = 0.0;
                h[2] = 0.0;
                //h[3] = 0.0;
                v[0] = uprev;
                q[0] = pprev;
                
                //*os << "h:";
                //for (int i=0; i<3; ++i)
                //*os << h[i] << '\t';
                //*os << endl;
                //*os << "s:";
                //for (int i=0; i<3; ++i)
                //*os << s[i] << '\t';
                //*os << endl;
                
                DNSFlags fineflags = DNS::flags();
                fineflags.verbosity=Silent;
                //*os << "constructing DNS for fine-time integration..." << endl;
                DNS dns(v[0], DNS::nu(), DNS::dt(), fineflags, tcoarse-dT);
                //*os << "finished contructing DNS..." << endl;
                
                int count=1; // need four data points for cubic interpolation
                
                // Now take a number of small-scale (dt) time-steps until we hit the section
                // Hitting the section will be defined by 0.0 lying between d[0] and d[2]
                for (Real tfine=tcoarse-dT; tfine <= tcoarse+dt; tfine += dt) {
                    
                    //*os << "time  shifts..." << endl;
                    // Time-shift velocity and Poincare condition arrays
                    // v[2] <- v[1] <- v[0], same w d in prep for advancing v[0],q[0] under DNS
                    for (int n=2; n>0; --n) {
                        s[n] = s[n-1];
                        v[n] = v[n-1];
                        q[n] = q[n-1];
                        h[n] = h[n-1];
                    }
                    
                    //*os << "time step..." << endl;
                    dns.advance(v[0], q[0], 1); // take one step of length dt
                    h[0] = (*h_)(v[0]);
                    s[0] = tfine + dt;
                    *os << ':' << flush;
                    
                    //*os << "crossing check..." << endl;
                    // Check for Poincare section crossing in midpoint of h[0],h[1],h[2],h[3]
                    if (++count >= 3 && ((h[2]<0 && 0<=h[0]) || (h[2]>0 && 0>=h[0]))) {
                        
                        // Newton search for zero of h(s) == h(v(s)) == 0 where v(s) is
                        // quadratic interpolant of v. Interpolating s as a function of h
                        // at h==0 gives a good initial guess for s
                        
                        // Newton iteration variables
                        Real sN = polynomialInterpolate(s, h, 0);
                        Real eps = 1e-9; // used for approximation dh/ds = (h(s + eps s) - h(s))/(eps s)
                        Real hsN, hsN_ds;
                        FlowField vN;
                        //os << "Newtown iteration on interpolated Poincare crossing" << endl;
                        
                        int Newtsteps = 6;
                        for (int n=0; n<Newtsteps; ++n) {
                            
                            vN = polynomialInterpolate(v, s, sN);
                            vN.makeSpectral();
                            hsN = (*h_)(vN);
                            //*os << n << flush;
                            
                            if (abs(hsN) < epsilon/2 || n==Newtsteps-1) {
                                if (abs(hsN) < epsilon/2)
                                    *os << "|" << flush;  // signal an accurate computation of a Poincare crossing
                                else
                                    *os << "~|" << flush; // signal an inaccurate computation of a Poincare crossing
                                
                                tcrossing_ = sN;
                                ucrossing_ = vN;
                                pcrossing_ = polynomialInterpolate(q, s, sN);
                                pcrossing_.makeSpectral();
                                break;
                            }
                            else {
                                vN = polynomialInterpolate(v, s, sN + eps*sN);
                                vN.makeSpectral();
                                hsN_ds = (*h_)(vN);
                                Real dhds = (hsN_ds - hsN)/(eps*sN);
                                Real ds = -hsN/dhds;
                                sN += ds;
                                
                                //*os << "Not good enough. Taking Newton step. " << endl;
                                //*os << "dhds == " << dhds << endl;
                                //*os << "ds   == " << ds << endl;
                                //*os << "s+ds == " << sN << endl;
                                
                            }
                        }
                        
                        // output time of crossing
                        hcrossing_ = hsN;
                        //Real cross = (*h_)(ucrossing_);
                        //*os << "Estimated poincare crossing: " << endl;
                        //*os << "  h(u) == " << cross << endl;
                        //*os << "  time == " << tcrossing_ << endl;
                        
                        return true;
                    }
                }
                *os << "Strange.... the large-scale steps crossed the Poincare section,\n";
                *os << "but when we went back and looked with finer-scale steps, there\n";
                *os << "was no crossing. Exiting." << endl;
                exit(1);
            }
        }
        return false; // didn't cross Poincare section
    }
    // END EXPERIMENTAL CODE
    // *******************************************************************************************
    
    
} //namespace channelflow
