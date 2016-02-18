//
//  controller.h
//  
//
//  Created by Peter Heins on 12/06/2012.
//  Copyright 2012 University of Sheffield. All rights reserved.
//

#ifndef CHANNELFLOW_CONTROLLER_H
#define CHANNELFLOW_CONTROLLER_H

#include <iostream>
#include "channelflow/flowfield.h"

extern "C" {
#include <clapack.h>
#include <cblas.h>
}

namespace channelflow {

typedef double Real;
    
// General functions
double** MatIn(const char* FileName);           // Reads matrices in binary format
double** MatIn_asc(const char* FileName);       // Reads matrices in ascii format
int NumRow(const char* FileName);               // Outputs no. of rows of matrix
int NumCol(const char* FileName);               // Outputs no. of columns of matrix
double* MatVec(double** Mat, double Vec[], int Matrow, int Matcol);         // Performs matrix-vector multiplication
void gauss(double** A, double* B, int N);                                   // Performs Gauss-Seidel iteration
double** eyeF(int numCstates);                                              // Outputs identity matrix
    
    
class Controller {
public:
    // Creates controller object
    Controller();
    Controller(int Ny, int minKx, int maxKx, int minKz, int maxKz, int uvw, const char* StateInfoFile, int NumConInputs, double tau, bool Spectral_states);
    
    double** CStateInfo();      // 2D Array for no. of controller states for controlled wavenumber pairs
    double*** ConStates();      // 3D Arra containing controller states for controlled wavenumber pairs
    
    // Advances discrete-time controller forward in time with fixed time-step dT
    void advance_Con(FlowField& u, FlowField& q, FlowField& BC, double*** CStateMat, Real dT, Real t, double*** IO);
    
    // Advances continuous-time controller forward in time with variable time-step dT
    void advance_Con_CN(FlowField& u, FlowField& q, FlowField& BC, double*** CStateMat, Real dT, Real t, double*** IO);
    
    
    int kx_c(int mx_c);
private:  
    int Ny_;       // No. of simulation wall-normal gridpoints
    int minKx_;    // Minimum kx to control
    int maxKx_;    // Maximum kx to control
    int minKz_;    // Minimum kz to control
    int maxKz_;    // Maximum kz to control
    int uvw_;      // Actuation via u, v or w
    const char* StateInfoFile_;     // Filename for file containing information of no. of controller states for each controlled wavenumber pair
    int NumConInputs_;              // No. of inputs INTO the controller, i.e. flow measurements
    double tau_;                    // Actuator time-constant
    bool Spectral_states_;          // States physical or Chebyshev spectral state
};

    
}


#endif
