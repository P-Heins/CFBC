/* tausolver.cpp: solves vector "tau" eqns (Canuto & Hussaini eqn 7.3.18-20)
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

#include <iomanip>
#include "channelflow/tausolver.h"

using namespace std;

namespace channelflow {

inline int kprime_func(int k, int Nb) {
  return (k%2 == 0) ? Nb-1 : Nb;
}

int n_func(int k, int Nb) {
  if (k==0)
    return Nb-1;
  else if (k==Nb)
    return 0;
  else if (k%2 == 0)
    return 2*(Nb-1);
  else
    return 2*Nb;
}

Real divcheck(string& label, int kx, int kz, Real kxLx, Real kzLz,
	      const ComplexChebyCoeff& u, const ComplexChebyCoeff& v,
	      const ComplexChebyCoeff& w, bool verbose) {

  // Verify divergence
  int N = u.length();
  if (v.length() != N || w.length() != N) {
    cout << "divcheck length problem!" << endl;
    exit(1);
  }

  ComplexChebyCoeff tmp = diff(v);
  Real diverge=0.0;
  for (int n=N-1; n>=0; --n)
    diverge += abs2(tmp[n] + (pi*2*(kxLx*u[n] + kzLz*w[n]))*I);
  diverge = sqrt(diverge);

  if (diverge > 1e-13 || verbose) {
    cout << label << endl;
    cout << "kx, kz == " << kx << ", " << kz << endl;
    cout << "kxLx, kzLz == " << kxLx << ", " << kzLz << endl;
    cout << "divergence == " << diverge << endl;
  }
  return diverge;
}


const Real MINIMUM_DISCRIMINANT = 1e-4;
const Real EPSILON = 1e-7;
TauSolver::TauSolver()
  :
  N_(0),
  Nb_(0),
  kx_(0),
  kz_(0),
  two_pi_kxLx_(0),
  two_pi_kzLz_(0),
  kappa2_(0),
  a_(0),
  b_(0),
  lambda_(0),
  nu_(0),
  tauCorrection_(true),
  pressureHelmholtz_(),
  velocityHelmholtz_(),
  P_0_(),
  v_0_(),
  P_plus_(),
  v_plus_(),
  P_minus_(),
  v_minus_(),
  i00_(0),
  i01_(0),
  i10_(0),
  i11_(0)
{}

    
TauSolver::TauSolver(int kx, int kz, Real Lx, Real Lz, Real a, Real b,
		     Real lambda, Real nu, int nChebyModes, bool tauCorrection)
  :
  N_(nChebyModes),
  Nb_(nChebyModes-1),
  kx_(kx),
  kz_(kz),
  two_pi_kxLx_(2*pi*kx/Lx),
  two_pi_kzLz_(2*pi*kz/Lz),
  kappa2_(4*square(pi)*(square(kx/Lx)+square(kz/Lz))),
  a_(a),
  b_(b),
  lambda_(lambda),
  nu_(nu),
  tauCorrection_(tauCorrection),
  pressureHelmholtz_(N_, a_, b_, kappa2_),
  velocityHelmholtz_(N_, a_, b_, lambda_, nu_),
  P_0_(N_, a_, b_, Spectral),
  v_0_(N_, a_, b_, Spectral),
  P_plus_(N_, a_, b_, Spectral),
  v_plus_(N_, a_, b_, Spectral),
  P_minus_(N_, a_, b_, Spectral),
  v_minus_(N_, a_, b_, Spectral),
  i00_(0),
  i01_(0),
  i10_(0),
  i11_(0)
{

  // ======================================================================
  // Calculate the influence matrix.

  // Make some local aliases to temp storage, for readability
  // This is not an efficiency issue, since TauSolvers are constructed once.
  ChebyCoeff dvplus_dy(N_, a_, b_, Spectral);
  ChebyCoeff dvminus_dy(N_, a_, b_, Spectral);
  ChebyCoeff zero(N_, a_, b_, Spectral);
  ChebyCoeff dPdy(N_, a_, b_, Spectral);

  // Solve homogeneous Helmholtz with P(-1) = 0, P(1) = 1.
  pressureHelmholtz_.solve(P_plus_, zero, 0.0, 1.0);   // eqn 7.3.25 discrete
  diff(P_plus_, dPdy);
  velocityHelmholtz_.solve(v_plus_, dPdy, 0.0, 0.0);   // eqn 7.3.26discrete

  // Solve homogeneous Helmholtz with P(-1) = 1, P(1) = 0.
  pressureHelmholtz_.solve(P_minus_, zero, 1.0, 0.0); // eqn 7.3.25 discrete
  diff(P_minus_, dPdy);
  velocityHelmholtz_.solve(v_minus_, dPdy, 0.0, 0.0); // eqn 7.3.26 discrete

  // Calculate influence matrix elements.
  diff(v_plus_,  dvplus_dy);
  diff(v_minus_, dvminus_dy);

  Real A = dvplus_dy.eval_b();
  Real B = dvminus_dy.eval_b();
  Real C = dvplus_dy.eval_a();
  Real D = dvminus_dy.eval_a();
  Real discriminant = A*D - B*C;

  // We know influence matrix is rank-deficient for kx==kz==0. It shouldn't
  // be for other wave numbers.
  if (kx_ !=0 || kz_!=0)
    assert((abs(discriminant) / Greater(abs(A*D), abs(B*C))) > MINIMUM_DISCRIMINANT);

  // Go ahead and divide by zero for kx==kz==0 case. Influence matrix is
  // unused in that case. Pollute the entries NaN's to make sure.
  i00_ =  D/discriminant;
  i01_ =  -B/discriminant;
  i10_ =  -C/discriminant;
  i11_ =  A/discriminant;

  // ====================================================================
  // Solve the B0 problem for tau corrections in solve(P,v)
  // !!!!STEP 1 in book!!!!
  ChebyCoeff p0_rhs(N_, a_, b_, Spectral);
  Real c = 2/(b_-a_);
  for (int i=0; i<=Nb_; ++i)
    p0_rhs[i] = c*n_func(i, Nb_);

  pressureHelmholtz_.solve(P_0_, p0_rhs, 0.0, 0.0); // eqn 7.3.41a

  ChebyCoeff dP0dy = diff(P_0_);
  //dP0dy_Nb1_ = dP0dy[Nb_-1];
  //dP0dy_Nb_  = dP0dy[Nb_];
  velocityHelmholtz_.solve(v_0_, dP0dy, 0.0, 0.0);  // eqn 7.3.41b

  influenceCorrection(P_0_, v_0_);

  ChebyCoeff v0yy = diff2(v_0_);
  //sigma0_Nb1_ = nu_*v0yy[Nb_-1] -(lambda_*v_0_[Nb_-1] + dP0dy[Nb_-1]);
  //sigma0_Nb_  = nu_*v0yy[Nb_]   -(lambda_*v_0_[Nb_]   + dP0dy[Nb_]);
  sigma0_Nb1_ = lambda_*v_0_[Nb_-1] + dP0dy[Nb_-1] - nu_*v0yy[Nb_-1];
  sigma0_Nb_  = lambda_*v_0_[Nb_]   + dP0dy[Nb_]   - nu_*v0yy[Nb_];

}
    
    

   
    
    
/************************************************************************    
 ***********************************************************************
 
 Written by: Peter H Heins (Postgraduate Research Student) 
 
 Date: January 2012
 
 Institution: University of Sheffield (ACSE Department)
 
 Purpose: This function adds general non-zero BCs at y=±1
 
 Based on code written by Binh Lieu, University of Minnesota March 2009
 
 *********************************************************************** 
*************************************************************************/ 
 
TauSolver::TauSolver(int kx, int kz, Real Lx, Real Lz, Real a, Real b, Real lambda, Real nu, int nChebyModes, Real vr_lower, Real vi_lower, Real vr_upper, Real vi_upper, bool tauCorrection)
    
    :
    N_(nChebyModes),
    Nb_(nChebyModes-1),
    kx_(kx),
    kz_(kz),
    two_pi_kxLx_(2*pi*kx/Lx),
    two_pi_kzLz_(2*pi*kz/Lz),
    kappa2_(4*square(pi)*(square(kx/Lx)+square(kz/Lz))),
    a_(a),
    b_(b),
    lambda_(lambda),
    nu_(nu),
    tauCorrection_(tauCorrection),
    pressureHelmholtz_(N_, a_, b_, kappa2_),
    velocityHelmholtz_(N_, a_, b_, lambda_, nu_),
    P_0_(N_, a_, b_, Spectral),
    v_0_(N_, a_, b_, Spectral),
    P_plus_(N_, a_, b_, Spectral),
    v_plus_(N_, a_, b_, Spectral),
    P_minus_(N_, a_, b_, Spectral),
    v_minus_(N_, a_, b_, Spectral),
    i00_(0),
    i01_(0),
    i10_(0),
    i11_(0),
    P_0_r(N_, a_, b_, Spectral),
    v_0_r(N_, a_, b_, Spectral),
    P_0_i(N_, a_, b_, Spectral),
    v_0_i(N_, a_, b_, Spectral)
  {
        
// =====================================================================
// Calculate the influence matrix as normal.

// Make some local aliases to temp storage, for readability
// This is not an efficiency issue, since TauSolvers are constructed once.
  ChebyCoeff dvplus_dy(N_, a_, b_, Spectral);
  ChebyCoeff dvminus_dy(N_, a_, b_, Spectral);
  ChebyCoeff zero(N_, a_, b_, Spectral);
  ChebyCoeff dPdy(N_, a_, b_, Spectral);
  
  // Solve homogeneous Helmholtz with P(-1) = 0, P(1) = 1.
  pressureHelmholtz_.solve(P_plus_, zero, 0.0, 1.0);   // eqn 7.3.25 discrete
  diff(P_plus_, dPdy);
  velocityHelmholtz_.solve(v_plus_, dPdy, 0.0, 0.0);   // eqn 7.3.26discrete
  
  // Solve homogeneous Helmholtz with P(-1) = 1, P(1) = 0.
  pressureHelmholtz_.solve(P_minus_, zero, 1.0, 0.0); // eqn 7.3.25 discrete
  diff(P_minus_, dPdy);
  velocityHelmholtz_.solve(v_minus_, dPdy, 0.0, 0.0); // eqn 7.3.26 discrete
  
  // Calculate influence matrix elements.
  diff(v_plus_,  dvplus_dy);
  diff(v_minus_, dvminus_dy);
  
  Real A = dvplus_dy.eval_b();
  Real B = dvminus_dy.eval_b();
  Real C = dvplus_dy.eval_a();
  Real D = dvminus_dy.eval_a();
  Real discriminant = A*D - B*C;
  
  // We know influence matrix is rank-deficient for kx==kz==0. It shouldn't
  // be for other wave numbers.
  if (kx_ !=0 || kz_!=0)
      assert((abs(discriminant) / Greater(abs(A*D), abs(B*C))) > MINIMUM_DISCRIMINANT);
  
  // Go ahead and divide by zero for kx==kz==0 case. Influence matrix is
  // unused in that case. Pollute the entries NaN's to make sure.
  i00_ =  D/discriminant;
  i01_ =  -B/discriminant;
  i10_ =  -C/discriminant;
  i11_ =  A/discriminant;
  
  // ====================================================================
  // Solve the REAL B0 problem for tau corrections in solve(P,v)
  // !!!!STEP 1 in book!!!!
  ChebyCoeff p0_rhs(N_, a_, b_, Spectral);
  Real c = 2/(b_-a_);
  for (int i=0; i<=Nb_; ++i)
      p0_rhs[i] = c*n_func(i, Nb_);
  
  pressureHelmholtz_.solve(P_0_r, p0_rhs, 0.0, 0.0); // eqn 7.3.41a
  
  ChebyCoeff dP0dy = diff(P_0_r);
  //dP0dy_Nb1_ = dP0dy[Nb_-1];
  //dP0dy_Nb_  = dP0dy[Nb_];
  
  // BCs on velocity Helmholtz equation changed to V^ at y=±1    
  velocityHelmholtz_.solve(v_0_r, dP0dy, vr_lower, vr_upper);  // eqn 7.3.41b
  
  influenceCorrection(P_0_r, v_0_r);
  
  ChebyCoeff v0yy = diff2(v_0_r);
  //sigma0_Nb1_ = nu_*v0yy[Nb_-1] -(lambda_*v_0_[Nb_-1] + dP0dy[Nb_-1]);
  //sigma0_Nb_  = nu_*v0yy[Nb_]   -(lambda_*v_0_[Nb_]   + dP0dy[Nb_]);
  sigma0_Nb1_r = lambda_*v_0_r[Nb_-1] + dP0dy[Nb_-1] - nu_*v0yy[Nb_-1];
  sigma0_Nb_r  = lambda_*v_0_r[Nb_]   + dP0dy[Nb_]   - nu_*v0yy[Nb_];
        
  
  // ====================================================================
  // Solve the IMAGINARY B0 problem for tau corrections in solve(P,v)
  // !!!!STEP 1 in book!!!!
  //ChebyCoeff p0_rhs(N_, a_, b_, Spectral);
  //Real c = 2/(b_-a_);
  //for (int i=0; i<=Nb_; ++i)
     // p0_rhs[i] = c*n_func(i, Nb_);
  
  pressureHelmholtz_.solve(P_0_i, p0_rhs, 0.0, 0.0); // eqn 7.3.41a
  
  ChebyCoeff dP0dyi = diff(P_0_i);
  //dP0dy_Nb1_ = dP0dy[Nb_-1];
  //dP0dy_Nb_  = dP0dy[Nb_];
  velocityHelmholtz_.solve(v_0_i, dP0dyi, vi_lower, vi_upper);  // eqn 7.3.41b
  
  influenceCorrection(P_0_i, v_0_i);
  
  ChebyCoeff v0yyi = diff2(v_0_i);
  //sigma0_Nb1_ = nu_*v0yy[Nb_-1] -(lambda_*v_0_[Nb_-1] + dP0dy[Nb_-1]);
  //sigma0_Nb_  = nu_*v0yy[Nb_]   -(lambda_*v_0_[Nb_]   + dP0dy[Nb_]);
  sigma0_Nb1_i = lambda_*v_0_i[Nb_-1] + dP0dy[Nb_-1] - nu_*v0yyi[Nb_-1];
  sigma0_Nb_i  = lambda_*v_0_i[Nb_]   + dP0dy[Nb_]   - nu_*v0yyi[Nb_];
      
      //cout << "B0 - problem!!" << endl;
      
        
  }
     
/**************************************************************************
 
                       ---End: Peter H Heins---
 
 **************************************************************************/
   
    void TauSolver::influenceCorrection(ChebyCoeff& P, ChebyCoeff& v) const {
        
        ChebyCoeff tmp = diff(v);
        Real dvp_dy_plus  = tmp.eval_b();
        Real dvp_dy_minus = tmp.eval_a();
        Real delta_plus  = -i00_*dvp_dy_plus - i01_*dvp_dy_minus;
        Real delta_minus = -i10_*dvp_dy_plus - i11_*dvp_dy_minus;
        
        // Add the influence matrix corrections to the particular solutions
        // to get a solution that statisfies both v(+-1)==0 and v'(+-1)==0.
        for (int i=0; i<N_; ++i) {
            P[i] += delta_plus*P_plus_[i] + delta_minus*P_minus_[i];
            v[i] += delta_plus*v_plus_[i] + delta_minus*v_minus_[i];
        }
    }        
    
    
    /*void TauSolver::influenceCorrection_uw(ChebyCoeff& P, ChebyCoeff& v, ChebyCoeff& u, ChebyCoeff& w, int kx, int kz) const {
        
        // Only for imaginary part!
        
        Real iku_plus = kx*u.eval_b;
        Real iku_minus = kx*u.eval_a;
        Real ikw_plus = kz*w.eval_b;
        Real ikw_minus = kz*w.eval_a;
        
        ChebyCoeff tmp = diff(v);
        Real dvp_dy_plus  = tmp.eval_b();
        dvp_dy_plus += iku_plus + ikw_plus;
        Real dvp_dy_minus = tmp.eval_a();
        dvp_dy_minus += iku_minus + ikw_minus;
        Real delta_plus  = -i00_*dvp_dy_plus - i01_*dvp_dy_minus;
        Real delta_minus = -i10_*dvp_dy_plus - i11_*dvp_dy_minus;
        
        // Add the influence matrix corrections to the particular solutions
        // to get a solution that statisfies both v(+-1)==0 and v'(+-1)==0.
        for (int i=0; i<N_; ++i) {
            P[i] += delta_plus*P_plus_[i] + delta_minus*P_minus_[i];
            v[i] += delta_plus*v_plus_[i] + delta_minus*v_minus_[i];
        }
    } */     
    

  // !!!!STEP 2 from book!!!!    
    
  void TauSolver::solve_P_and_v(ChebyCoeff& P, ChebyCoeff& v, const ChebyCoeff& r,
  const ChebyCoeff& Ry, Real& sigmaNb1, Real& sigmaNb)

  const {

  // P is Canuto & Hussaini's Ppart particular solution after this solve
  pressureHelmholtz_.solve(P, r, 0.0, 0.0); // eqn 7.3.25 discrete HH1

  // kx==kz==0 is a degenerate case for which the influence matrix is
  // rank-deficient, the v solution is identically zero and P satisfies
  // P' == Ry (which is equivalent to P'' == div(R) == r, found above).
  if (kx_== 0 && kz_==0) {
    for (int i=0; i<N_; ++i)
      v[i] = 0.0;
    return;
  }

  // The rest of this method is for the case kx != 0 or kz != 0.
  ChebyCoeff tmp = diff(P);
  tmp -= Ry;

      
   /*   ORIGINAL!!!!!!!
  // v is Canuto & Hussaini's vpart particular solution after this solve
  velocityHelmholtz_.solve(v, tmp, 0.0, 0.0); // eqn 7.3.25 discrete HH2
  influenceCorrection(P,v);/ */
      
  // MY CODE!!!!!
  // v is Canuto & Hussaini's vpart particular solution after this solve
  velocityHelmholtz_.solve(v, tmp, 0.0, 0.0); // eqn 7.3.25 discrete HH2
  influenceCorrection(P,v);    

  // Jump ship if not doing tau correction
  if (!tauCorrection_)
    return;

  // Add up tau correction terms. Quotes are from Canuto & Hussaini pg 219
  // My Nb is Canuto's N.
  // My sigma1_Nb  is Canuto's sigma_{1m} for m=N
  // My sigma1_Nb1 is Canuto's sigma_{1m} for m=N-1

  // Set tmp = vyy;
  diff2(v,tmp);
      
  
  // "define sigma_{1m} and sigma_{0m} for m = N-1, N as the tau
  // terms that must be added to the v-momentum eqns for (P1,v1) and
  // (P0,v0) for them to hold"
  Real sigma1_Nb  = lambda_*v[Nb_]   - nu_*tmp[Nb_]   - Ry[Nb_];
  Real sigma1_Nb1 = lambda_*v[Nb_-1] - nu_*tmp[Nb_-1] - Ry[Nb_-1];
  diff(P, tmp);
  sigma1_Nb  += tmp[Nb_];
  sigma1_Nb1 += tmp[Nb_-1];

  // !!!!STEP 3 from book!!!!    
      
  // "One can show that sigma_m = sigma_{1m}/(1-sigma_{0m}, m=N-1,N"
  sigmaNb  = sigma1_Nb/(1.0 - sigma0_Nb_);
  sigmaNb1 = sigma1_Nb1/(1.0 - sigma0_Nb1_);
      
  // !!!!STEP 4 from book!!!!    

  // " and that ..."
  for (int i=0; i<=Nb_; ++i) {
    P[i]   += ((i%2 == 0) ? sigmaNb1 : sigmaNb) * P_0_[i];
    v[i]   += ((i%2 == 0) ? sigmaNb  : sigmaNb1) * v_0_[i];
  }

  //#ifdef DEBUG
  //verify_P_and_v(P,v,r,Ry, sigmaNb1, sigmaNb, true);
  //#endif

  return;
}

    
/************************************************************************    
 ***********************************************************************
 
 Written by: Peter H Heins (Postgraduate Research Student) 
 
 Date: January 2012
 
 Institution: University of Sheffield (ACSE Department)
 
 Purpose: This function solves the tau equations with non-zeros BCs at y=±1 for Pressure and v-component velocity
 
 Based on code written by Binh Lieu, University of Minnesota March 2009
 
 *********************************************************************** 
 *************************************************************************/  
    
    
void TauSolver::solve_P_and_v_inhom(ChebyCoeff& P, ChebyCoeff& v, ChebyCoeff& r, const ChebyCoeff& Ry, Real v_lower, Real v_upper, Real& sigmaNb1, Real& sigmaNb, const Real& sigma0_Nb_, const Real& sigma0_Nb1_)
    
const{
    
    // P is Canuto & Hussaini's Ppart particular solution after this solve
    pressureHelmholtz_.solve(P, r, 0.0, 0.0); // eqn 7.3.25 discrete HH1
    
    // kx==kz==0 is a degenerate case for which the influence matrix is
    // rank-deficient, the v solution is identically zero and P satisfies
    // P' == Ry (which is equivalent to P'' == div(R) == r, found above).
    if (kx_== 0 && kz_==0) {
        for (int i=0; i<N_; ++i)
            v[i] = 0.0;
        return;
    }
    
    // The rest of this method is for the case kx != 0 or kz != 0.
    ChebyCoeff tmp = diff(P);
    tmp -= Ry;
    
    // v is Canuto & Hussaini's vpart particular solution after this solve
    // with v^ at y=±1 as BCs on Helmholtz equation    
    velocityHelmholtz_.solve(v, tmp, v_lower, v_upper); // eqn 7.3.25 discrete HH2
    influenceCorrection(P,v);    
    
    // Jump ship if not doing tau correction
    if (!tauCorrection_)
        return;
    
    // Add up tau correction terms. Quotes are from Canuto & Hussaini pg 219
    // My Nb is Canuto's N.
    // My sigma1_Nb  is Canuto's sigma_{1m} for m=N
    // My sigma1_Nb1 is Canuto's sigma_{1m} for m=N-1
    
    // Set tmp = vyy;
    diff2(v,tmp);
    
    
    // "define sigma_{1m} and sigma_{0m} for m = N-1, N as the tau
    // terms that must be added to the v-momentum eqns for (P1,v1) and
    // (P0,v0) for them to hold"
    Real sigma1_Nb  = lambda_*v[Nb_]   - nu_*tmp[Nb_]   - Ry[Nb_];
    Real sigma1_Nb1 = lambda_*v[Nb_-1] - nu_*tmp[Nb_-1] - Ry[Nb_-1];
    diff(P, tmp);
    sigma1_Nb  += tmp[Nb_];
    sigma1_Nb1 += tmp[Nb_-1];
    
    // !!!!STEP 3 from book!!!!    
    
    // "One can show that sigma_m = sigma_{1m}/(1-sigma_{0m}, m=N-1,N"
    sigmaNb  = sigma1_Nb/(1.0 - sigma0_Nb_);
    sigmaNb1 = sigma1_Nb1/(1.0 - sigma0_Nb1_);
    
    // !!!!STEP 4 from book!!!!    
    
    // " and that ..."
    for (int i=0; i<=Nb_; ++i) {
        P[i]   += ((i%2 == 0) ? sigmaNb1 : sigmaNb) * P_0_[i];
        v[i]   += ((i%2 == 0) ? sigmaNb  : sigmaNb1) * v_0_[i];
    }
    
    //#ifdef DEBUG
    //verify_P_and_v(P,v,r,Ry, sigmaNb1, sigmaNb, true);
    //#endif
    
    //cout << "solve_P_V" << endl;
    
    return;
}

       
/**************************************************************************
 
 ---End: Peter H Heins---
 
 **************************************************************************/    
    
    
    

Real TauSolver::verify_P_and_v(const ChebyCoeff& P, const ChebyCoeff& v,
			       const ChebyCoeff& r, const ChebyCoeff& Ry,
			       Real sigmaNb1, Real sigmaNb, bool verbose) const {

  Real error = 0.0;

  if (verbose) {
    cout << "--------TauSolver::verify_P_and_v(P,v,r,Ry,sigmaNb1,sigmaNb) {" << endl;
    cout << "nu = " << nu_ << "; lambda = " << lambda_ << endl;
  }

  ChebyCoeff Py =  diff(P);
  ChebyCoeff Pyy = diff(Py);
  ChebyCoeff p_lhs(P);
  p_lhs *= -kappa2_;
  p_lhs += Pyy;

  Real l2err = L2Dist(p_lhs, r);
  Real t2err = tauDist(p_lhs, r);
  error += l2err;
  if (verbose) {
    cout << "Homog tauDist(P'' - k^2 P, r)   == " << t2err << endl;
    cout << "Homog  L2Dist(P'' - k^2 P, r)   == " << l2err << endl;
  }

  ChebyCoeff vy = diff(v);
  ChebyCoeff vyy = diff(vy);
  ChebyCoeff v_lhs = nu_*vyy - lambda_*v - Py;
  ChebyCoeff minusRy = Ry;
  minusRy *= -1.0;

  t2err = tauDist(v_lhs, minusRy);
  l2err = L2Dist(v_lhs, minusRy);
  error += l2err;
  if (verbose) {
    cout << "tauDist(v_lhs, -Ry) == " << t2err << endl;
    cout << " L2Dist(v_lhs, -Ry) == " << l2err << endl;
  }

  // Now compare P and v eqns with tau correction terms.
  ChebyCoeff sigma(N_,a_,b_,Spectral);
  sigma[Nb_-1] = sigmaNb1;
  sigma[Nb_]   = sigmaNb;
  ChebyCoeff p_rhs(N_,a_,b_,Spectral);
  diff(sigma, p_rhs);
  p_rhs += r;

  t2err = tauDist(p_lhs, p_rhs);
  l2err = L2Dist(p_lhs, p_rhs);
  error += l2err;
  if (verbose) {
    cout << "Homog tauDist(P'' - k^2 P, r + sigma')   == " << t2err << endl;
    cout << "Homog  L2Dist(P'' - k^2 P, r + sigma')   == " << l2err << endl;
  }

  ChebyCoeff v_rhs(minusRy);
  v_rhs -= sigma;

  t2err = tauDist(v_lhs, v_rhs);
  l2err = L2Dist(v_lhs, v_rhs);
  error += l2err;
  if (verbose) {
    cout << "tauDist(v_lhs, -Ry - sigma) == " << tauDist(v_lhs, v_rhs) << endl;
    cout << " L2Dist(v_lhs, -Ry - sigma) == " << L2Dist(v_lhs, v_rhs)  << endl;
  }

  Real vp = v.eval_b();
  Real vm = v.eval_a();
  Real vyp = vy.eval_b();
  Real vym = vy.eval_a();
  Real v_norm =
    (abs(nu_*L2Norm(vyy)) + abs(lambda_*L2Norm(v))) + L2Norm(Py);
  v_norm = (v_norm > EPSILON) ? v_norm : 1.0;

  Real p_norm = L2Norm(Pyy) + kappa2_*L2Norm(P);
  p_norm = (p_norm > EPSILON) ? p_norm : 1.0;


  error += abs(vp)/v_norm;
  error += abs(vm)/v_norm;
  error += abs(vyp)/v_norm;
  error += abs(vym)/v_norm;

  if (verbose) {
    cout << " v(+/- 1) == " << vp << ' ' << vm << endl;
    cout << "v'(+/- 1) == " << vyp << ' ' << vym << endl;
  }

  //Real v_err = tauDist(v_lhs,v_rhs)/v_norm;
  //Real p_err = tauDist(p_lhs,p_rhs)/p_norm;
  //assert(v_err < EPSILON);
  //assert(p_err < EPSILON);
  if (verbose)
    cout << "} TauSolver::verify_P_and_v" << endl;

  return error;

}


void TauSolver::solve(ComplexChebyCoeff& u, ComplexChebyCoeff& v,
		      ComplexChebyCoeff& w, ComplexChebyCoeff& P,
		      const ComplexChebyCoeff& Rx,
		      const ComplexChebyCoeff& Ry,
		      const ComplexChebyCoeff& Rz) const
{

  ComplexChebyCoeff r(N_,a_,b_,Spectral);
  Real sigmaNb1;
  Real sigmaNb;

  ChebyCoeff rr(N_,a_,b_,Spectral);

  // Always solve v(y) from momentum, and get P.
  // Re and Im parts of v and P eqns decouple. Solve them seperately.
  diff(Ry.re, rr);
  int n; // MSVC++ FOR-SCOPE BUG
  for (n=0; n<N_; ++n)
    rr[n] -= two_pi_kxLx_*Rx.im[n] + two_pi_kzLz_*Rz.im[n];
  solve_P_and_v(P.re, v.re, rr, Ry.re, sigmaNb1, sigmaNb);

  diff(Ry.im, rr);
  for (n=0; n<N_; ++n)
    rr[n] += two_pi_kxLx_*Rx.re[n] + two_pi_kzLz_*Rz.re[n];
  solve_P_and_v(P.im, v.im, rr, Ry.im, sigmaNb1, sigmaNb);

  // Re and Im parts of u and w eqns seperate.
  // Use r as temporary space to store RHS of eqns.
  for (n=0; n<N_; ++n)
    r.set(n, two_pi_kxLx_*I*P[n] - Rx[n]);
  //Complex c = pi2i*kxLx_*P[n] - Rx[n];
  //r.set(n,c);

  velocityHelmholtz_.solve(u.re, r.re, 0.0, 0.0);
  velocityHelmholtz_.solve(u.im, r.im, 0.0, 0.0);

  for (n=0; n<N_; ++n)
    r.set(n, two_pi_kzLz_*I*P[n] - Rz[n]);
  //Complex c = pi2i*kzLz_*P[n] - Rz[n];
  //r.set(n, c);

  velocityHelmholtz_.solve(w.re, r.re, 0.0, 0.0);
  velocityHelmholtz_.solve(w.im, r.im, 0.0, 0.0);

  // This is for debugging ONLY
  /**************************
  if ((two_pi_kxLx_ == 0.0 & kx_ != 0) ||
      (two_pi_kzLz_ == 0.0 & kz_ != 0)) {
    verify(u,v,w,P,Rx,Ry,Rz,true);
    u.setToZero();
    v.setToZero();
    w.setToZero();
    P.setToZero();
  }
  **************************/

#ifdef DEBUG
  //verify(u,v,w,P,Rx,Ry,Rz);
#endif
  return;
}

void TauSolver::solve(ComplexChebyCoeff& u, ComplexChebyCoeff& v,
		      ComplexChebyCoeff& w, ComplexChebyCoeff& P,
		      Real& dPdx,
		      const ComplexChebyCoeff& Rx,
		      const ComplexChebyCoeff& Ry,
		      const ComplexChebyCoeff& Rz, Real umean) const {

  // This function should only be called for kx==kz==0, since the enforcing
  // const velocity flux makes sense only for that case. Divergence is not a
  // problem here, so solve everything via momentum.
  assert(kx_==0 && kz_==0);

  ComplexChebyCoeff r(N_,a_,b_,Spectral);
  Real pi2 = 2*pi;
  Complex pi2i = pi2*I;
  Real sigmaNb1;            // tau correction term
  Real sigmaNb;             // tau correction term

  ChebyCoeff rr(N_,a_,b_,Spectral);

  // Re and Im parts of v and P eqns decouple. Solve them seperately.
  diff(Ry.re, rr);
  int n; // MSVC++ FOR-SCOPE BUG
  for (n=0; n<N_; ++n)
    rr[n] -= two_pi_kxLx_*Rx.im[n] + two_pi_kzLz_*Rz.im[n];
  solve_P_and_v(P.re, v.re, rr, Ry.re, sigmaNb1, sigmaNb);

  diff(Ry.im, rr);
  for (n=0; n<N_; ++n)
    rr[n] += two_pi_kxLx_*Rx.re[n] + two_pi_kzLz_*Rz.re[n];
  solve_P_and_v(P.im, v.im, rr, Ry.im, sigmaNb1, sigmaNb);

  // Re and Im parts of u and w eqns seperate.
  // Use r as temporary space to store RHS of eqns.
  for (n=0; n<N_; ++n)
    r.set(n, two_pi_kxLx_*I*P[n] - Rx[n]);

  velocityHelmholtz_.solve(u.re, dPdx, r.re, umean, 0.0, 0.0);
  velocityHelmholtz_.solve(u.im, r.im, 0.0, 0.0);

  for (n=0; n<N_; ++n)
    r.set(n, two_pi_kzLz_*I*P[n] - Rz[n]);

  velocityHelmholtz_.solve(w.re, r.re, 0.0, 0.0);
  velocityHelmholtz_.solve(w.im, r.im, 0.0, 0.0);

#ifdef DEBUG
  //verify(u,v,w,P, dPdx, Rx,Ry,Rz, umean);
#endif


  return;
}
    
    
/************************************************************************    
 ***********************************************************************
 
 Written by: Peter H Heins (Postgraduate Research Student) 
 
 Date: January 2012
 
 Institution: University of Sheffield (ACSE Department)
 
 Purpose: This function solves the tau equations with non-zeros BCs at y=±1 for the case when kx,kz != 0, finding the u and w components of velocity
 
 Based on code written by Binh Lieu, University of Minnesota March 2009
 
 *********************************************************************** 
 *************************************************************************/ 
    
void TauSolver::solve_Inhom(ComplexChebyCoeff& u, ComplexChebyCoeff& v,
                      ComplexChebyCoeff& w, ComplexChebyCoeff& P, const ComplexChebyCoeff& Rx,
                      const ComplexChebyCoeff& Ry, const ComplexChebyCoeff& Rz, Real vr_lower, Real vr_upper, Real vi_lower, Real vi_upper, Complex u_lower, Complex u_upper, Complex w_lower, Complex w_upper) const
{
    
    ComplexChebyCoeff r(N_,a_,b_,Spectral);
    Real sigmaNb1;
    Real sigmaNb;
    
    ChebyCoeff rr(N_,a_,b_,Spectral);
    
    // Always solve v(y) from momentum, and get P.
    // Re and Im parts of v and P eqns decouple. Solve them seperately.
    diff(Ry.re, rr);
    int n; // MSVC++ FOR-SCOPE BUG
    
    // Use new solve P_and_v_inhom function which includes v^ at y=±1.
    for (n=0; n<N_; ++n)
        rr[n] -= two_pi_kxLx_*Rx.im[n] + two_pi_kzLz_*Rz.im[n];
    solve_P_and_v_inhom(P.re, v.re, rr, Ry.re, vr_lower, vr_upper, sigmaNb1, sigmaNb, sigma0_Nb_r, sigma0_Nb1_r);
    
    diff(Ry.im, rr);
    for (n=0; n<N_; ++n)
        rr[n] += two_pi_kxLx_*Rx.re[n] + two_pi_kzLz_*Rz.re[n];
    solve_P_and_v_inhom(P.im, v.im, rr, Ry.im, vi_lower, vi_upper, sigmaNb1, sigmaNb,sigma0_Nb_i, sigma0_Nb1_i);
    
    // Re and Im parts of u and w eqns seperate.
    // Use r as temporary space to store RHS of eqns.
    for (n=0; n<N_; ++n)
        r.set(n, two_pi_kxLx_*I*P[n] - Rx[n]);
    //Complex c = pi2i*kxLx_*P[n] - Rx[n];
    //r.set(n,c);
    
    Real ur_lower = real(u_lower);
    Real ui_lower = imag(u_lower);
    Real ur_upper = real(u_upper);
    Real ui_upper = imag(u_upper); 
    
    // Uses U^ at y=±1 as BCs
    velocityHelmholtz_.solve(u.re, r.re, ur_lower, ur_upper);
    velocityHelmholtz_.solve(u.im, r.im, ui_lower, ui_upper);
    
       
    for (n=0; n<N_; ++n)
        r.set(n, two_pi_kzLz_*I*P[n] - Rz[n]);
    //Complex c = pi2i*kzLz_*P[n] - Rz[n];
    //r.set(n, c);
    
    Real wr_lower = real(w_lower);
    Real wi_lower = imag(w_lower);
    Real wr_upper = real(w_upper);
    Real wi_upper = imag(w_upper); 
    
    // Uses W^ at y=±1 as BCs
    velocityHelmholtz_.solve(w.re, r.re, wr_lower, wr_upper);
    velocityHelmholtz_.solve(w.im, r.im, wi_lower, wi_upper);
    
        
#ifdef DEBUG
    //verify(u,v,w,P,Rx,Ry,Rz);
#endif
    
    //cout << "Solve_Inhom" << endl;
    
    return;
}

void TauSolver::solve_Inhom(ComplexChebyCoeff& u, ComplexChebyCoeff& v,
                      ComplexChebyCoeff& w, ComplexChebyCoeff& P,
                      Real& dPdx,
                      const ComplexChebyCoeff& Rx,
                      const ComplexChebyCoeff& Ry,
                      const ComplexChebyCoeff& Rz, Real vr_lower, Real vr_upper, Real vi_lower, Real vi_upper, Real umean) const 
    {
    
    // This function should only be called for kx==kz==0, since the enforcing
    // const velocity flux makes sense only for that case. Divergence is not a
    // problem here, so solve everything via momentum.
    assert(kx_==0 && kz_==0);
    
    ComplexChebyCoeff r(N_,a_,b_,Spectral);
    Real pi2 = 2*pi;
    Complex pi2i = pi2*I;
    Real sigmaNb1;            // tau correction term
    Real sigmaNb;             // tau correction term
    
    ChebyCoeff rr(N_,a_,b_,Spectral);
    
    // Re and Im parts of v and P eqns decouple. Solve them seperately.
    diff(Ry.re, rr);
    int n; // MSVC++ FOR-SCOPE BUG
    for (n=0; n<N_; ++n)
        rr[n] -= two_pi_kxLx_*Rx.im[n] + two_pi_kzLz_*Rz.im[n];
        //solve_P_and_v(P.re, v.re, rr, Ry.re, sigmaNb1, sigmaNb);
        solve_P_and_v_inhom(P.re, v.re, rr, Ry.re, vr_lower, vr_upper, sigmaNb1, sigmaNb, sigma0_Nb_r, sigma0_Nb1_r);
    
    diff(Ry.im, rr);
    for (n=0; n<N_; ++n)
        rr[n] += two_pi_kxLx_*Rx.re[n] + two_pi_kzLz_*Rz.re[n];
        //solve_P_and_v(P.im, v.im, rr, Ry.im, sigmaNb1, sigmaNb);
        solve_P_and_v_inhom(P.im, v.im, rr, Ry.im, vi_lower, vi_upper, sigmaNb1, sigmaNb, sigma0_Nb_i, sigma0_Nb1_i);
    
    // Re and Im parts of u and w eqns seperate.
    // Use r as temporary space to store RHS of eqns.
    for (n=0; n<N_; ++n)
        r.set(n, two_pi_kxLx_*I*P[n] - Rx[n]);
    
    // Uses U^ at y=±1 as BCs
    //velocityHelmholtz_.solve(u.re, r.re, U_lower, U_upper);
    //velocityHelmholtz_.solve(u.im, r.im, 0.0, 0.0);
        
   /* Real ur_lower = real(u_lower);
    Real ui_lower = imag(u_lower);
    Real ur_upper = real(u_upper);
    Real ui_upper = imag(u_upper);   */ 
    
    velocityHelmholtz_.solve(u.re, dPdx, r.re, umean,0.0, 0.0);
    velocityHelmholtz_.solve(u.im, r.im, 0.0, 0.0);
    
    for (n=0; n<N_; ++n)
        r.set(n, two_pi_kzLz_*I*P[n] - Rz[n]);
    
    /* Real wr_lower = real(w_lower);
    Real wi_lower = imag(w_lower);
    Real wr_upper = real(w_upper);
    Real wi_upper = imag(w_upper);   */
        
    // Uses W^ at y=±1 as BCs
    velocityHelmholtz_.solve(w.re, r.re, 0.0, 0.0);
    velocityHelmholtz_.solve(w.im, r.im,0.0, 0.0);
    
#ifdef DEBUG
    //verify(u,v,w,P, dPdx, Rx,Ry,Rz, umean);
#endif
    
    
    return;
}
    
    
/**************************************************************************
 
 ---End: Peter H Heins---
 
 **************************************************************************/     
    
    
    
    
    

Real TauSolver::verify(const ComplexChebyCoeff& u, const ComplexChebyCoeff& v,
		       const ComplexChebyCoeff& w, const ComplexChebyCoeff& P,
		       const ComplexChebyCoeff& Rx,const ComplexChebyCoeff& Ry,
		       const ComplexChebyCoeff& Rz, bool verbose) const {

  Real umean = Re(u.mean());
  Real dPdx = 0.0;
  return verify(u,v,w,P,dPdx,Rx,Ry,Rz,umean,verbose);
}

Real TauSolver::verify(const ComplexChebyCoeff& u, const ComplexChebyCoeff& v,
		       const ComplexChebyCoeff& w, const ComplexChebyCoeff& P,
		       Real dPdx,
		       const ComplexChebyCoeff& Rx,const ComplexChebyCoeff& Ry,
		       const ComplexChebyCoeff& Rz,
		       Real umean, bool verbose) const {

  // verify   nu u''(y) - lambda u(y) - grad P = -R,
  //                                     div u = 0
  //                                    u(+-1) = 0

  if (verbose) {
    cout << "TauSolver::verify(u,v,w,P,dPdx,Rx,Ry,Rz,umean,verbose) {" << endl;
    cout << " kx kz == " << kx_ << ' ' << kz_ << endl;
  }
  ComplexChebyCoeff lhs(N_,a_,b_,Spectral);
  ComplexChebyCoeff tmp(N_,a_,b_,Spectral);
  Real error = 0.0;
  Real terr = 0.0;
  Real lerr = 0.0;

  // Verify u eqn, -nu u'' + lambda u + dP/dx == Rx
  lhs = u;
  lhs *= lambda_;
  diff2(u, tmp);
  tmp *= nu_;
  lhs -= tmp;
  tmp = P;
  //tmp *= pi2*I*kxLx_; gcc-2.95 bug makes the following necessary
  tmp *= Complex(0.0, two_pi_kxLx_);
  lhs += tmp;
  lhs.re[0] += dPdx;

  terr = tauDist(lhs, Rx);
  lerr = L2Dist(lhs, Rx);
  error += lerr;
  if (verbose) {
    cout << "L2Norm(Rx) == " << L2Norm(Rx) << endl;
    cout << "tauDist(nu u'' - lambda u - dP/dx, -Rx) == " << terr << endl;
    cout << " L2Dist(nu u'' - lambda u - dP/dx, -Rx) == " << lerr << endl;
  }

  // Verify v eqn.
  lhs = v;
  lhs *= lambda_;
  diff2(v, tmp);
  tmp *= nu_;
  lhs -= tmp;
  diff(P, tmp);
  lhs += tmp;
  terr = tauDist(lhs, Ry);
  lerr = L2Dist(lhs, Ry);
  error += lerr;
  if (verbose) {
    cout << "L2Norm(Ry) == " << L2Norm(Ry) << endl;
    cout << "tauDist(nu v'' - lambda v - dP/dy, -Ry) == " << terr << endl;
    cout << " L2Dist(nu v'' - lambda v - dP/dy, -Ry) == " << lerr << endl;
  }

  // Verify w eqn
  lhs = w;
  lhs *= lambda_;
  diff2(w, tmp);
  tmp *= nu_;
  lhs -= tmp;
  tmp = P;
  tmp *= Complex(0.0, two_pi_kzLz_);
  lhs += tmp;
  terr = tauDist(lhs, Rz);
  lerr = L2Dist(lhs, Rz);
  error += lerr;
  if (verbose) {
    cout << "L2Norm(Rz) == " << L2Norm(Rz) << endl;
    cout << "tauDist(nu w'' - lambda w - dP/dz, -Rz) == " << terr << endl;
    cout << " L2Dist(nu w'' - lambda w - dP/dz, -Rz) == " << lerr << endl;
  }

  // Verify P eqn. P'' - kappa^2 P = div R
  diff2(P, lhs);
  tmp = P;
  tmp *= -kappa2_;
  lhs += tmp;

  ComplexChebyCoeff r(N_,a_,b_,Spectral);
  ChebyCoeff r_re(N_,a_,b_,Spectral);
  ChebyCoeff r_im(N_,a_,b_,Spectral);

  // Re and Im parts of v and P eqns decouple. Solve them seperately.
  diff(Ry.re, r_re);
  int n; // MSVC++ FOR-SCOPE BUG
  for (n=0; n<N_; ++n)
    r_re[n] -= two_pi_kxLx_*Rx.im[n] + two_pi_kzLz_*Rz.im[n];

  diff(Ry.im, r_im);
  for (n=0; n<N_; ++n)
    r_im[n] += two_pi_kxLx_*Rx.re[n] + two_pi_kzLz_*Rz.re[n];

  r.re = r_re;
  r.im = r_im;
  terr = tauDist(lhs, r);
  lerr = L2Dist(lhs, r);
  error += lerr;
  if (verbose) {
    cout << "L2Norm(div R) == " << L2Norm(r) << endl;
    cout << "tauDist(P'' - k^2 P, div R)   == " << terr << endl;
    cout << " L2Dist(P'' - k^2 P, div R)   == " << lerr << endl;
  }
  // Don't assert on pressure eqn, since it's modified for tau correction.
  //else
  //assert(error<EPSILON);


  // Verify divergence
  diff(v, tmp);
  for (n=0; n<N_; ++n)
    tmp.add(n, I*(two_pi_kxLx_*u[n] + two_pi_kzLz_*w[n]));
    //tmp.add(n, I*pi2*(kxLx_*u[n] + kzLz_*w[n]));

  terr = tauNorm(tmp);
  lerr = L2Norm(tmp);
  error += lerr;
  if (verbose) {
    cout << "tauNorm(div) == " << terr << endl;
    cout << " L2Norm(div) == " << lerr << endl;
  }

  // BCs!!!!!!!!!!!!!!!!!!!
  Complex ua = u.eval_a();
  Complex ub = u.eval_b();
  error += abs(ua) + abs(ub);
  if (verbose)
    cout << "u(a),u(b) == " << ua << " " << ub << endl;


  Complex va = v.eval_a();
  Complex vb = v.eval_b();
  error += abs(va) + abs(vb);
  if (verbose)
    cout << "v(a),v(b) == " << va << " " << vb << endl;

  ComplexChebyCoeff vy = diff(v);
  Complex vya = vy.eval_a();
  Complex vyb = vy.eval_b();
  error += abs(vya) + abs(vyb);
  if (verbose)
    cout << "v' at a,b == " << vya << " " << vyb << endl;

  Complex wa = u.eval_a();
  Complex wb = u.eval_b();
  error += abs(wa) + abs(wb);
  if (verbose)
    cout << "w(a),w(b) == " << wa << " " << wb << endl;

  Real mean_error = abs2(Re(u.mean()) - umean);
  error += mean_error;
  if (verbose)
    cout << "abs2(u.mean() - umean) == " << mean_error << endl;
  else
    assert(mean_error < EPSILON);

  if (verbose) {
    cout << "total verification error == " << error << endl;
    cout << "} TauSolver::verify(...)" << endl;
  }
  return error;
}

Real TauSolver::tauNorm(const ChebyCoeff& u) const {
  ChebyCoeff tmp(u.numModes()-2, u);
  return L2Norm(u);
}
Real TauSolver::tauNorm(const ComplexChebyCoeff& u) const {
  ComplexChebyCoeff tmp(u.numModes()-2, u);
  return L2Norm(u);
}
Real TauSolver::tauDist(const ChebyCoeff& u, const ChebyCoeff& v) const {
  ChebyCoeff utmp(u.numModes()-2, u);
  ChebyCoeff vtmp(v.numModes()-2, v);
  return L2Dist(utmp,vtmp);
}
Real TauSolver::tauDist(const ComplexChebyCoeff& u, const ComplexChebyCoeff& v) const {
  ComplexChebyCoeff utmp(u.numModes()-2, u);
  ComplexChebyCoeff vtmp(v.numModes()-2, v);
  return L2Dist(utmp,vtmp);
}

} //namespace channelflow
