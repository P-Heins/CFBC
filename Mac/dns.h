// dns.h: time-integration class for spectral Navier-Stokes simulator
// channelflow-1.3 www.channelflow.org
// copyright (C) 2001-2009 John F. Gibson, license declaration at end of file

// DNS is a class for integrating Navier-Stokes equation.
// DNSFlags is used to specify the integration parameters of DNS.
// TimeStep manages variable time-stepping, adjusting dt to keep CFL in range


#ifndef CHANNELFLOW_DNS_H
#define CHANNELFLOW_DNS_H

//#include <Accelerate/Accelerate.h>
#include "channelflow/controller.h"
#include "channelflow/mathdefs.h"
#include "channelflow/vector.h"
#include "channelflow/chebyshev.h"
#include "channelflow/array.h"
#include "channelflow/flowfield.h"
#include "channelflow/diffops.h"
#include "channelflow/tausolver.h"
#include "channelflow/symmetry.h"
#include <fftw3.h>

namespace channelflow {

// Enum types for specifying the behavior of DNS, fields of DNSFlags.
enum BaseFlow        {Zero, PlaneCouette, Parabolic};
enum MeanConstraint  {PressureGradient, BulkVelocity};
enum TimeStepMethod  {CNFE1, CNAB2, CNRK2, SMRK2, SBDF1, SBDF2, SBDF3, SBDF4};
enum NonlinearMethod {Rotational, Convection, Divergence, SkewSymmetric,
		      Alternating, Alternating_, LinearAboutField, LinearAboutProfile};
enum Dealiasing      {NoDealiasing, DealiasXZ, DealiasY, DealiasXYZ};
enum Verbosity       {Silent, PrintTicks, PrintTime, VerifyTauSolve, PrintAll};

// CNFE1  == Crank-Nicolson Forward-Euler order 1   (no init steps needed)
// CNAB2  == Crank-Nicolson Adams-Bashforth order 2 (needs 1 init steps)
// SMRK2  == Spalart,Moser,R?, Runge-Kutta order 2  (no init steps needed)
// SBDFn  == Semimplicit backwards-differentiation order n (needs n-1 init steps)
// See channelflow manual for detailed description of these algorithms
// Note: CNFE1 and SBDF1 are the same algorithm.

TimeStepMethod  s2stepmethod(const std::string& s);
NonlinearMethod s2nonlmethod(const std::string& s);
Dealiasing      s2dealaising(const std::string& s);
MeanConstraint  s2constraint(const std::string& s);
BaseFlow        s2baseflow(const std::string& s);

void  navierstokesNL(const FlowField& u, const ChebyCoeff& Ubase,
		     FlowField& f, FlowField& tmp, NonlinearMethod& method);

/*****************************************************************
void  navierstokesNL(const FlowField& u, const FlowField& ubase,
		     const ChebyCoeff& Ubase, const FlowField& ubtot,
		     const FlowField& grad_ubtot, FlowField& f,
		     FlowField& tmp, NonlinearMethod& method);
********************************************/

void  navierstokesNL(const FlowField& u, const FlowField& ubase,
		     const ChebyCoeff& Ubase, FlowField& f, FlowField& tmp,
		     FlowField& tmp2, NonlinearMethod& method);

// Specify the behavior of NSIntegrators by setting fields of DNSFlags.
class DNSFlags {
public:
  //       Option type     Option name      Default value
  DNSFlags(BaseFlow        baseflow       = Zero,
	   MeanConstraint  constraint     = PressureGradient,
	   TimeStepMethod  timestepping   = SBDF3,
	   TimeStepMethod  initstepping   = SMRK2,
	   NonlinearMethod nonlinearity   = Rotational,
	   Dealiasing      dealiasing     = DealiasXZ,
	   bool            taucorrection  = true,
	   Real            dPdx           = 0.0,
	   Real            Ubulk          = 0.0,
	   Verbosity       verbosity      = PrintTicks,
	   std::ostream*   logstream      = &std::cout);

  BaseFlow        baseflow;     // utot = u + Ubase(y) ex, Ubase incorps BCs
  MeanConstraint  constraint;   // Enforce const press grad or const bulk vel
  TimeStepMethod  timestepping; // Time-stepping algorithm
  TimeStepMethod  initstepping; // Algorithm for initializing multistep methods
  NonlinearMethod nonlinearity; // Method of calculating nonlinearity of NS eqn
  Dealiasing      dealiasing;   // Use 3/2 rule to eliminate aliasing
  bool            taucorrection;// Remove divergence caused by discretization
  Real            dPdx;         // Constraint value
  Real            Ubulk;        // Constraint value
  Verbosity       verbosity;    // Print diagnostics, times, ticks, or nothing
  std::ostream*   logstream;           // stream for output

  array<FieldSymmetry> symmetries; //restrict u(t) to these symmetries

  bool dealias_xz() const;
  bool dealias_y() const;
};


// TimeStep keeps dt between dtmin and dtmax, and CFL between CFLminand CFLmax,
// in such a way that dt*n = dT for some integer n. That's useful if you
// want to plot/save data at regular dT intervals, but use a variable timestep
// dt for efficiency. You can mandate a fixed timestep by setting dtmin==dtmax.
// For example of use, see example codes.

class TimeStep {
public:
  TimeStep();
  TimeStep(Real dt, Real dtmin, Real dtmax, Real dT, Real CFLmin,
	   Real CFLmax, bool variable=true);

  // If variable, adjust dt to keep CFLmin<=CFL<=CFLmax (soft),
  // and dtmin<=dt<=dtmax (hard). Returns true if dt changes, false otherwise
  bool adjust(Real CFL, bool verbose=true, std::ostream& os=std::cout);
  bool adjustToMiddle(Real CFL, bool verbose=true, std::ostream& os=std::cout);
  bool adjust_for_T(Real T, bool verbose=true, std::ostream& os=std::cout); // tweak dT and dt to fit T exactly


  int  n() const;        // n*dt == dT
  int  N() const;        // N*dT == T
  Real dt() const;       // integration timestep
  Real dtmin() const;
  Real dtmax() const;
  Real dT() const;       // plot/CFL-check interval
  Real T() const;        // total integration time
  Real CFL() const;
  Real CFLmin() const;
  Real CFLmax() const;
  bool variable() const;
  operator Real() const; // same as dt()

private:
  int n_;
  int N_;
  Real dt_;
  Real dtmin_;  //
  Real dtmax_;  //
  Real dT_;     // dT_ == n_*dt_, plot interval
  Real T_;      // T_  == N_*dt_, total integration time
  Real CFLmin_;
  Real CFL_;
  Real CFLmax_;
  bool variable_;
};

std::ostream& operator<<(std::ostream& os, const TimeStep& ts);

class DNSAlgorithm;

// DNS is a wrapper class for DNSAlgorithms. It's the main class for
// integrating the Navier-Stokes equations in top-level programs.
// Specify the integration algorithm and other parameters in the DNSFlags.
// If you like, you can construct and use specific DNS algorithms like
// MultiStepDNS in top-level programs --any class derived from DNSAlgorithm.
// Look in example codes for examples of initialization and use.

class DNS {
public:
  DNS();
  DNS(const DNS& dns);

  //DNS(const FlowField& u, Real nu, Real dt, const DNSFlags& flags, Real t=0.0);
    
  DNS(const FlowField& u, Real nu, Real dt, const DNSFlags& flags, Real t=0.0, bool controlled=false);


  DNS(const FlowField& u, const ChebyCoeff& Ubase,
      Real nu, Real dt, const DNSFlags& flags, Real t=0.0, bool controlled=false);

  ~DNS();

  DNS& operator=(const DNS& dns);

  void advance(FlowField& u, FlowField& q, int nSteps=1);
  
     
  //*****Added by Peter H Heins******************
  void advance_NL(FlowField& u, FlowField& q, FlowField& F, int nSteps=1);

  void advance_inhom_CON(Controller& cont, FlowField& u, FlowField& q, FlowField& BCs, FlowField& F, int nSteps=1, double*** IO=0, double*** CStateMat=0);
  void advance_inhom(FlowField& u, FlowField& q, FlowField& BCs, int nSteps=1);
  //void advance_inhom_CON_SF(Controller& cont, FlowField& u, FlowField& q, FlowField& BCs, FlowField& F, int nSteps=1, double*** IO=0);
  //*********************************************
    
    
  void project();      // Project onto symmetric subspace
  void operator *= (const FieldSymmetry& symm);  // Apply symmetry to internal fields

  // Convert potentially fake pressure q and true pressure p, back and forth.
  void up2q(const FlowField& u, const FlowField& p, FlowField& q) const;
  void uq2p(const FlowField& u, const FlowField& q, FlowField& p) const;

  //void reset();                  // flush state, prepare for new integration
  void reset_dt(Real dt);
    
  void reset_dtIH(Real dt, FlowField& BCs);  // ***Added by Peter H Heins***
  
  void reset_time(Real t);
  void reset_dPdx(Real dPdx);    // change dPdx and enforce const dPdx
  void reset_Ubulk(Real Ubulk);  // change Ubulk and enforce const Ubulk


  //void reset_uj(const FlowField& uj, int j);  // set u[j]=u(t-j*dt)
  bool push(const FlowField& u); // push into u[j] stack, true when full, t+=dt
  bool full() const;             // pushed enough init data into u[j]?

  int order() const;             // err should scale as dt^order
  int Ninitsteps() const;        // number of steps needed to initialize

  Real nu() const;
  Real dt() const;
  Real CFL() const;
  Real time() const;
  Real dPdx() const;      // the mean pressure gradient at the current time
  Real Ubulk() const;     // the actual bulk velocity at the current time
  Real dPdxRef() const;   // the mean press grad enforced during integration
  Real UbulkRef() const;  // the bulk velocity enforced during integ.

  const DNSFlags& flags() const;
  TimeStepMethod timestepping() const;

  void printStack() const;

 private:
  DNSAlgorithm* main_algorithm_;        // same as
  DNSAlgorithm* init_algorithm_;

  //DNSAlgorithm* newAlgorithm(const FlowField& u, const ChebyCoeff& Ubase, Real nu,
			    // Real dt, const DNSFlags& flags, Real t);
    
  DNSAlgorithm* newAlgorithm(const FlowField& u, const ChebyCoeff& Ubase, Real nu,
                               Real dt, const DNSFlags& flags, Real t, bool controlled);  
};


// DNSAlgorithm is a base class for classes representing time-stepping
// algorithms for the Navier-Stokes equations, using a Fourier x Chebyshev
// x Fourier FlowField for spatial discretization and finite-differencing
// and tau method for temporal discretization.

class DNSAlgorithm {
public:
  DNSAlgorithm();
  DNSAlgorithm(const DNSAlgorithm& dns);
  //DNSAlgorithm(const FlowField& u, const ChebyCoeff& Ubase, Real nu, Real dt,
//	       const DNSFlags& flags, Real t=0);
   
  DNSAlgorithm(const FlowField& u, const ChebyCoeff& Ubase, Real nu, Real dt,
                 const DNSFlags& flags, Real t=0, bool controlled=false);   


  virtual ~DNSAlgorithm();
  //DNSAlgorithm& operator=(const DNSAlgorithm& dns);

  virtual void advance(FlowField& u, FlowField& q, int nSteps=1) = 0;

  //***********Added by Peter H Heins*******************
  virtual void advance_NL(FlowField& u, FlowField& q, FlowField& F, int nSteps=1) = 0;

  virtual void advance_inhom_CON(Controller& cont, FlowField& u, FlowField& q, FlowField& BCs, FlowField& F, int nSteps=1, double*** IO=0, double*** CStateMat=0) = 0;  
  virtual void advance_inhom(FlowField& u, FlowField& q, FlowField& BCs, int nSteps=1) = 0; 
  //virtual void advance_inhom_CON_SF(Controller& cont, FlowField& u, FlowField& q, FlowField& BCs, FlowField& F, int nSteps=1, double*** IO=0) = 0;
    //****************************************************  
    
  virtual void project(); // project onto symm subspace (a member of flags)
  virtual void operator *= (const FieldSymmetry& symm);   // apply symmetry operator


  //UNIMPLEMENTED
  // Convert potentially fake pressure q and true pressure p, back and forth.
  //  void up2q(const FlowField& u, const FlowField& p, FlowField& q) const;
  //  void uq2p(const FlowField& u, const FlowField& q, FlowField& p) const;


  //virtual void reset();             // flush state, prepare for new integration
  virtual void reset_dt(Real dt) = 0;     // somewhat expensive
  
  virtual void reset_dtIH(Real dt, FlowField& BCs) =0; //**Added by PHH**  
    
  virtual bool push(const FlowField& u);  // push u onto u[j] stack, t += dt
  virtual bool full() const;              // have enough init data?


  void reset_time(Real t);
  void reset_dPdx(Real dPdx);    // change dPdx and enforce const dPdx
  void reset_Ubulk(Real Ubulk);  // change Ubulk and enforce const Ubulk

  int order() const;             // err should scale as dt^order
  int Ninitsteps() const;        // number of steps needed to initialize

  int Nx() const;
  int Ny() const;
  int Nz() const;

  // These methods are not implemented.
  //  int Mx() const;
  //  int My() const;
  //  int Mz() const;

  Real Lx() const;
  Real Lz() const;
  Real a() const;
  Real b() const;
  Real nu() const;
  Real dt() const;
  Real CFL() const;
  Real time() const;
  Real dPdx() const;      // the mean pressure gradient at the current time
  Real Ubulk() const;     // the actual bulk velocity at the current time
  Real dPdxRef() const;   // the mean press grad enforced during integration
  Real UbulkRef() const;  // the bulk velocity enforced during integ.

  const DNSFlags& flags() const;
  const ChebyCoeff& Ubase() const;
  const FlowField& ubase() const;
  TimeStepMethod timestepping() const;

  virtual DNSAlgorithm* clone() const = 0;  // new copy of *this

  virtual void printStack() const;

protected:
  // Spatial parameters
  int Nx_;      // number of X gridpoints
  int Ny_;      // number of Chebyshev T(y) modes
  int Nz_;      // number of Z gridpoints
  int Mx_;      // number of X modes
  int Mz_;      // number of Z modes
  int Nyd_;     // number of dealiased Chebyshev T(y) modes
  int kxd_max_; // maximum value of kx among dealiased modes
  int kzd_max_; // maximum value of kz among dealiased modes
  Real Lx_;
  Real Lz_;
  Real a_;
  Real b_;

  // Temporal integration parameters
  DNSFlags flags_; // User-defined integration parameters
  int order_;
  int Ninitsteps_; // number of initialization steps required
  Real nu_;
  Real dt_;
  Real t_;         // time in convective units
  Real cfl_;       // CFL number
  Real dPdxRef_;   // Enforced mean pressure gradient (0.0 if unused).
  Real dPdxAct_;   // Actual mean pressure gradient at previous timestep.
  Real UbulkRef_;  // Enforced total bulk velocity (0.0 if unused).
  Real UbulkAct_;  // Actual total bulk velocity bulk obtained.
  Real UbulkBase_; // Bulk velocity of Ubase
  Real ubulkBase_; // bulk velocity of ubase

  ChebyCoeff Ubase_;   // baseflow physical
  ChebyCoeff Ubaseyy_; // baseflow'' physical

  FlowField ubase_;    // non-null only if linearizing about genl field...
  FlowField ubtot_;    // (ubase+Ubase)
  FlowField tmp2_;     // grad(ubase+Ubase)
  FlowField tmp_;      // tmp space for nonlinearity calculation

  // These variables are used as temp storage when solving indpt tau problems.
  ComplexChebyCoeff uk_;   // profile of u_{kx,kz} (y) at t = n dt
  ComplexChebyCoeff vk_;
  ComplexChebyCoeff wk_;
  ComplexChebyCoeff Pk_;   // profile of P_{kx,kz} (y)
  ComplexChebyCoeff Pyk_;  // profile of dP_{kx,kz}/dy (y)
  ComplexChebyCoeff Rxk_;
  ComplexChebyCoeff Ryk_;
  ComplexChebyCoeff Rzk_;

  int kxmaxDealiased() const;
  int kzmaxDealiased() const;
  bool isAliasedMode(int kx, int kz) const;

  void init(FlowField& u); // common constructor code
};


// Multistep algorithms, all using Backwards Differentiation: SBDFk
// Based on Peyret section. The order is set by flags.timestepping.
class MultistepDNS : public DNSAlgorithm {
public:
  MultistepDNS();
  MultistepDNS(const MultistepDNS& dns);
//MultistepDNS(const FlowField& u, const ChebyCoeff& Ubase, Real nu, Real dt,
//	       const DNSFlags& flags, Real t=0);
    
  MultistepDNS(const FlowField& u, const ChebyCoeff& Ubase, Real nu, Real dt,
                 const DNSFlags& flags, Real t=0, bool controlled=false);  
    
  ~MultistepDNS();

  MultistepDNS& operator=(const MultistepDNS& dns);

  virtual void advance(FlowField& u, FlowField& q, int nSteps=1);
  
  //***********Added by Peter H Heins************************
  virtual void advance_NL(FlowField& u, FlowField& q, FlowField& F, int nSteps=1);

  virtual void advance_inhom_CON(Controller& cont, FlowField& u, FlowField& q, FlowField& BCs, FlowField& F, int nSteps=1, double*** IO=0, double*** CStateMat=0);  
  virtual void advance_inhom(FlowField& u, FlowField& q, FlowField& BCs, int nSteps=1);
  //virtual void advance_inhom_CON_SF(Controller& cont, FlowField& u, FlowField& q, FlowField& BCs, FlowField& F, int nSteps=1, double*** IO=0);
  //*********************************************************  
    
  virtual void project();
  virtual void operator *= (const FieldSymmetry& symm);
  virtual void reset_dt(Real dt);
    
  virtual void reset_dtIH(Real dt, FlowField& BCs); //**Added by PHH** 
    
  virtual bool push(const FlowField& u); // for initialization
  virtual bool full() const;             // have enough init data?
  //virtual void reset();       // flush state, prepare for new integration


  virtual DNSAlgorithm* clone() const;  // new copy of *this

  virtual void printStack() const;

protected:
  Real eta_;
  array<Real> alpha_;
  array<Real> beta_;
  array<FlowField> u_;  // u[j] == u at t-j*dt for multistep algorithms
  array<FlowField> f_;  // f[j] == f at t-j*dt for multistep algorithms

  TauSolver** tausolver_;  // 2d array of tausolvers, indexed by [mx][mz]

  int countdown_;
};

// CNRK2 and hopefully another. Based on algorithm in Peyret pg 149
class RungeKuttaDNS : public DNSAlgorithm {
public:
  RungeKuttaDNS();
  RungeKuttaDNS(const RungeKuttaDNS& dns);
  //RungeKuttaDNS(const FlowField& u, const ChebyCoeff& Ubase, Real nu, Real dt,
	//	const DNSFlags& flags, Real t=0);
    
  RungeKuttaDNS(const FlowField& u, const ChebyCoeff& Ubase, Real nu, Real dt,
                  const DNSFlags& flags, Real t=0, bool controlled=false);  
    
  ~RungeKuttaDNS();

  RungeKuttaDNS& operator=(const RungeKuttaDNS& dns);

  virtual void advance(FlowField& u, FlowField& q, int nSteps=1);
    
  //***********Added by Peter H Heins************************
  virtual void advance_NL(FlowField& u, FlowField& q, FlowField& F, int nSteps=1);

  virtual void advance_inhom_CON(Controller& cont, FlowField& u, FlowField& q, FlowField& BCs, FlowField& F, int nSteps=1, double*** IO=0, double*** CStateMat=0);
  virtual void advance_inhom(FlowField& u, FlowField& q, FlowField& BCs, int nSteps=1);
  //virtual void advance_inhom_CON_SF(Controller& cont, FlowField& u, FlowField& q, FlowField& BCs, FlowField& F, int nSteps=1, double*** IO=0);
  
  //*********************************************************  
    
  virtual void reset_dt(Real dt);

  virtual void reset_dtIH(Real dt, FlowField& BCs); //**Added by PHH**     
    
  // next few functions are no-ops so base class defns suffice
  // virtual void project();
  // virtual void operator *= (const FieldSymmetry& symm);
  // virtual void reset();  // flush state, prepare for new integration

  virtual DNSAlgorithm* clone() const;  // new copy of *this
protected:
  int Nsubsteps_;
  FlowField Qj1_;  // Q_{j-1} (Q at previous substep)
  FlowField Qj_;   // Q_j     (Q at current  substep)
  array<Real> A_;  // Q_{j+1} = A_j Q_j + N(u_j)
  array<Real> B_;  // u_{j+1} = u_j + dt B_j Q_j + dt C_j (L u_j + L u_{j+1})
  array<Real> C_;

  TauSolver*** tausolver_; // 3d array indexed by [step][mx][mz]
};

// A generalization of CNAB2 with substeps. Implements CNAB2 and SMRK2
class CNABstyleDNS : public DNSAlgorithm {
public:
  CNABstyleDNS();
  CNABstyleDNS(const CNABstyleDNS& dns);
  //CNABstyleDNS(const FlowField& u, const ChebyCoeff& Ubase, Real nu, Real dt,
	//       const DNSFlags& flags, Real t=0);
    
    
  CNABstyleDNS(const FlowField& u, const ChebyCoeff& Ubase, Real nu, Real dt,
                 const DNSFlags& flags, Real t=0, bool controlled=false);  
    
  ~CNABstyleDNS();

  CNABstyleDNS& operator=(const CNABstyleDNS& dns);

  virtual void advance(FlowField& u, FlowField& q, int nSteps=1);
    
  //***********Added by Peter H Heins************************
  virtual void advance_NL(FlowField& u, FlowField& q, FlowField& F, int nSteps=1);

  virtual void advance_inhom_CON(Controller& cont, FlowField& u, FlowField& q, FlowField& BCs, FlowField& F, int nSteps=1, double*** IO=0, double*** CStateMat=0);
  virtual void advance_inhom(FlowField& u, FlowField& q, FlowField& BCs, int nSteps=1);
  //virtual void advance_inhom_CON_SF(Controller& cont, FlowField& u, FlowField& q, FlowField& BCs, FlowField& F, int nSteps=1, double*** IO=0);
  //*********************************************************
    
  virtual void project();
  virtual void operator *= (const FieldSymmetry& symm);   // apply symmetry operator
  virtual void reset_dt(Real dt);
  virtual void reset_dtIH(Real dt, FlowField& BCs); //**Added by PHH**    
    
  virtual bool push(const FlowField& u);  // push u onto u[j] stack, t += dt
  virtual bool full() const;              // have enough init data?
  virtual void printStack() const;


private:
  int Nsubsteps_;
  bool full_;
  FlowField fj1_;      // f_{j-1} (f at previous substep)
  FlowField fj_;       // f_j     (f at current  substep)
  array<Real> alpha_;  // u_{j+1} = u_j + dt L (alpha_j u_j + beta_j u_{j+1})
  array<Real> beta_;    //           + dt gamma_j N(u_j) + dt zeta N(u_{j-1})
  array<Real> gamma_;
  array<Real> zeta_;

  TauSolver*** tausolver_; // 3d array indexed by [i][mx][mz]

  virtual DNSAlgorithm* clone() const;    // new copy of *this
};

std::ostream& operator<<(std::ostream& os, BaseFlow b);
std::ostream& operator<<(std::ostream& os, MeanConstraint m);
std::ostream& operator<<(std::ostream& os, NonlinearMethod n);
std::ostream& operator<<(std::ostream& os, TimeStepMethod t);
std::ostream& operator<<(std::ostream& os, Verbosity v);
std::ostream& operator<<(std::ostream& os, Dealiasing d);
std::ostream& operator<<(std::ostream& os, const DNSFlags& flags);



// Given a baseflow, fluctation, modified pressure triple (U0,u0,q0) and
// a new baseflow U1, compute new fluctuation u1 and modified pressure q1.
//     utot == U0 + u0 == U1 + u1
// pressure == q0 - 1/2 || u0 ||^2 == q1 - 1/2 || u1 ||^2
void changeBaseFlow(const ChebyCoeff& U0, const FlowField& u0, const FlowField& q0,
		    const ChebyCoeff& U1, FlowField& u1, FlowField& q1);




// **********************************************************************************
// BEGIN EXPERIMENTAL CODE: DNS that integrates to a Poincare section and maps back to a
// fundamental domain via symmetries, whenever certain boundaries are crossed.


// DNSPoincare is a class for integrating u to a Poincare section and mapping u back
// into a fundamental domain of a discrete symmetry group. The Poincare intersections
// are well-tested. The fundamental domain stuff is not. For the Poincare section,
// the stopping condition is a geometric condition rather than a stopping time.
// As of now there are two forms for the Poincare condtion, I-D=0 (DragDissipation)
// or (u(t) - ustar, estar) = 0 (PlaneIntersection). The fundamental domain is defined as
// (u(t), e[n]) >= 0. e[n] is antisymmetric under symmetry sigma[n], so that when
// (u(t), e[n]) <  0, we can get back to the fundamental domain by applying sigma[n],
// since (sigma[n] u(t), e[n]) = (u(t), sigma[n] e[n]) = (u(t), -e[n]) > 0.

// advanceToSection should be used this way
// 1. Start with some initial condition (u,q,t), t arbitrary
// 2. Make repeated calls to advanceToSection(u,q,nsteps,crosssign,eps). Each call will
//    advance (u,q,t) by dT = nSteps*dt and map u(t) back into the fundamental domain
//    should it leave. The sign argument determines which kinds of crossings will
//    return: sign<0 => only dh/dt < 0,
//            sign=0 => either direction
//            sign>0 => only dh/dt > 0,
// 3. Check the return value of advanceToSection. It will be
//      FALSE if u(t) does not cross the section during the the advancement and
//      TRUE  if u(t) does cross the section
// 4. When the return value is TRUE, you can then access the values of (u,q,t)
//    at the crossing through ucrossing(), pcrossing(), and tcrossing().
// 5. The signcrossing() member function returns the sign of dh/dt at h==0
// 6. Continue on with more calls to advanceToSection to find the next intersection,
//    if you like.
// The integration to the section is done in multiple steps and multiple calls
// to advanceToSection so that the field can be saved, projected, etc.
// over the course of integration by the caller.

class PoincareCondition {
public:
  PoincareCondition();
  virtual Real operator()(const FlowField& u) = 0;
};

// Section defined by (u - ustar, estar) == (u, estar) - (ustar, estar)
class PlaneIntersection : public PoincareCondition {
public:
  ~PlaneIntersection();
  PlaneIntersection();
  PlaneIntersection(const FlowField& ustar, const FlowField& estar);
  Real operator()(const FlowField& u);
private:
  FlowField estar_; // A normal that defines orientation of section
  Real      cstar_; // L2IP(ustar, estar), defines location of section
};


// Section defined by I-D == drag(u) - dissipation(u) == 0
class DragDissipation : public PoincareCondition {
public:
  DragDissipation();
  Real operator()(const FlowField& u);
};


class DNSPoincare : public DNS {

public:
  DNSPoincare();

  DNSPoincare(const FlowField& u, PoincareCondition* h, Real nu, Real dt, const DNSFlags& flags, Real t0=0);

  DNSPoincare(const FlowField& u, const array<FlowField>& e, const array<FieldSymmetry>& sigma,
	      PoincareCondition* h, Real nu, Real dt, const DNSFlags& flags, Real t0=0);

  bool advanceToSection(FlowField& u, FlowField& q, int nSteps, int crosssign=0, Real Tmin=0,
			Real epsilon=1e-13);

  //Real f(const FlowField& u) const; // poincare condition is f(u) == 0


  const FlowField& ucrossing() const;
  const FlowField& pcrossing() const;
  Real hcrossing() const;      // value of poincare condition at crossing
  Real tcrossing() const;      // time of poincare crossing
  int  scrossing() const;      // -1 or 1 for dh/dt<0 or dh/dt>0

  Real hcurrent() const;       // return h(u) at current timestep
private:
  array<FlowField> e_;         // Defines fundamental domain. See comments above.
  array<FieldSymmetry> sigma_; // Maps u(t) back into fundamental domain. See above.
  PoincareCondition* h_;       // The poincare condition h(u) == 0.

  FlowField ucrossing_;        // velocity field at crossing
  FlowField pcrossing_;        // pressure field at crossing
  Real      tcrossing_;        // time at crossing
  int       scrossing_;        // sign of dh/dt at crossing
  Real      hcrossing_;        // value of (*h)(ucrossing_)
  Real      hcurrent_;         // value of (*h)(u)
  Real t0_;                    // starting time, used to check t-t0 >= Tmin
};
// END EXPERIMENTAL CODE

} //namespace channelflow
#endif


/* dns.h: time-integration class for spectral Navier-Stokes simulator
 *
 * channelflow-1.3, www.channelflow.org
 *
 * Copyright (C) 2001-2009  John F. Gibson
 *
 * gibson@cns.physics.gatech.edu
 * Center for Nonlinear Sciences, School of Physics
 * Georgia Institute of Technology
 * Atlanta, GA 30332-0430
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
