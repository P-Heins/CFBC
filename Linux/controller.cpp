//
//  controller.cpp
//  
//
//  Created by Peter Heins on 27/02/2012.
//  Copyright 2012 University of Sheffield. All rights reserved.
//
#include <iostream>

#include "channelflow/controller.h"


using namespace std;

namespace channelflow 
{
    // Reads in matrix from file
    double** MatIn(const char* FileName) {
        
        double* memblock;
        int size;
        
        ifstream Matfile(FileName, ios_base::in|ios_base::binary|ios::ate);
        if (Matfile.is_open()) {
            
            size = (int) Matfile.tellg();            
            memblock = new double[size];
            Matfile.seekg(0,ios::beg);
            Matfile.read((char*)memblock,size); 
            Matfile.close(); 
            
            // read in matrices' dimensions
            int num_row = iround(memblock[0]);
            int num_col = iround(memblock[1]);

            
            double **Mat;
            Mat = new double*[num_row];
            for (int i=0;i<num_row;++i){
                Mat[i] = new double[num_col];
            }        
            
            int Vecsize = num_row*num_col;
            double Vec[Vecsize];
            
            for (int i=0; i<Vecsize; ++i) { // vector
                 Vec[i] = memblock[i+2];
            }
            
            int count =0;
            for (int row=0; row<num_row; ++row) {
                for (int col=0; col<num_col; ++col) { // matrix
                    Mat[row][col] = Vec[count];
                    ++count;
                }
            }
            
            delete[] memblock;
            return Mat;
            //delete[] Mat;

            
        }
        else {
            cout << "Unable to open file" << endl;
            return 0;             
        }     
        
        
                
    }
    
    
    double** MatIn_asc(const char* FileName) {
        
        fstream Matfile(FileName, ios_base::in|ios_base::out);
        Matfile.seekg(0);
        
        // read in matrices' dimensions
        int num_row;
        int num_col;
        Matfile >> num_row >> num_col;
        
        double **Mat;
        Mat = new double*[num_row];
        for (int i=0;i<num_row;++i){
            Mat[i] = new double[num_col];
        }        
        
        int Vecsize = num_row*num_col;
        double Vec[Vecsize];
        
        for (int i=0; i<Vecsize; ++i) { // vector
            Matfile >> Vec[i];
        }
        
        int count =0;
        for (int row=0; row<num_row; ++row) {
            for (int col=0; col<num_col; ++col) { // matrix
                Mat[row][col] = Vec[count];
                ++count;
            }
        }
        
        return Mat;
        
        
        
    }

    
    // Reads in number of rows of matrix in file
    int NumRow(const char* FileName) {
        
        double* memblock;
        int size;
        
        ifstream Matfile(FileName, ios_base::in|ios_base::binary|ios::ate);
        size = (int) Matfile.tellg();            
        memblock = new double[size];
        Matfile.seekg(0,ios::beg);
        Matfile.read((char*)memblock,size); 
        Matfile.close(); 
        
        // read in matrices' dimensions
        int num_row = iround(memblock[0]);
       
        delete[] memblock;
        
        return num_row;    
        
        
    }
    // Reads in number of columnss of matrix in file
    int NumCol(const char* FileName) {
        
        double* memblock;
        int size;
        
        ifstream Matfile(FileName, ios_base::in|ios_base::binary|ios::ate);
        size = (int) Matfile.tellg();            
        memblock = new double[size];
        Matfile.seekg(0,ios::beg);
        Matfile.read((char*)memblock,size); 
        Matfile.close(); 
        
        // read in matrices' dimensions
        int num_col = iround(memblock[1]);       
        
        
        delete[] memblock;

        return num_col; 

    }
    // Multiplies matrix with vector
    double* MatVec(double** Mat, double Vec[], int Matrow, int Matcol) {
        
       
        double *Result;
        Result = new double[Matrow];
        
        for (int i=0;i<Matrow;++i){ // Initialise result vector
            Result[i] = 0.0;
        }
        
        double resultmat[Matrow][Matcol];
        for (int j=0;j<Matcol;++j) {
            for (int i=0;i<Matrow;++i) {
                resultmat[i][j] = Mat[i][j]*Vec[j];
            }
        }
        
        for (int i=0;i<Matrow;++i) {
            for (int j=0;j<Matcol;++j) {
                Result[i] += resultmat[i][j] ;
            }
        }
        return Result;  
        //delete[] Result;
        
    }
    
    
        
    
    // Solves Ax=B for x
    void gauss(double** A, double* B, int N){
        
                
        float* alpha;
	alpha = new float[N*N];
        float* beta;
        beta = new float[N];
        
        int cnt = 0;
        
        for (int i=0;i<N;++i){
	  beta[i] = B[i];
            for (int j=0;j<N;++j){
	        alpha[cnt] = A[j][i];
                ++cnt;
            }
        }
       
        int *IPIV = new int[N+1];
        
        int INFO;
	
        int NRHS = 1;
        
        
        INFO = clapack_sgetrf(CblasColMajor,N,N,alpha,N,IPIV); 

        if (INFO==0){
	  INFO = clapack_sgetrs(CblasColMajor,CblasNoTrans,N,NRHS,alpha,N,IPIV,beta,N);
      
	  if (INFO==0){

	    for (int i=0;i<N;++i){
	      B[i] = beta[i];
	    }                    
	  }
	  else {
	    cerr << "Error: Gauss LU" << endl;
	  }
	}
	else {
	  cerr << "Error: Gauss" << endl;
	}

        delete[] IPIV;
	delete[] alpha;
	delete[] beta;
    }



    
    // returns square identity matrix dimension numCstates
    double** eyeF(int numCstates) {
        
        // Identity matrix, size of Ak
        double** eyeA;
        eyeA = new double*[numCstates];
        for (int i=0;i<numCstates;++i){
            eyeA[i] = new double[numCstates];
        }
        
        
        for (int i=0;i<numCstates;++i){
            for (int j=0;j<numCstates;++j){
                if (i-j==0) {
                    eyeA[i][j] = 1.0;
                }
                else {
                    eyeA[i][j] = 0.0;
                }
            }
        }
        
        return eyeA;
        //delete[] eyeA;
        
    }
    
    
    
    
    
    
    // CONTROLLER CLASS
    
    
    Controller::Controller()
    :
    Ny_(0),
    minKx_(0),
    maxKx_(0),
    minKz_(0),
    maxKz_(0),
    uvw_(0),
    StateInfoFile_("Mult_Control_Mat/StateInfo.bin"),
    NumConInputs_(0),
    tau_(0.0),
    Spectral_states_(false)
    {}
    
  Controller::Controller(int Ny, int minKx, int maxKx, int minKz, int maxKz, int uvw, const char* StateInfoFile, int NumConInputs, double tau, bool Spectral_states)
    :
    Ny_(Ny),
    minKx_(minKx),
    maxKx_(maxKx),
    minKz_(minKz),
    maxKz_(maxKz),
    uvw_(uvw),
    StateInfoFile_(StateInfoFile),
    NumConInputs_(NumConInputs),
    tau_(tau),
    Spectral_states_(Spectral_states)
    {}
    
    // returns array containing number of controller states at each kx,kz pair
    double** Controller::CStateInfo() {
        
        
        double** InfoMat;
        int minKx = minKx_;
        int maxKx = maxKx_;
        int minKz = minKz_;
        int maxKz = maxKz_;
        
        int maxMx = (maxKx-minKx)+1;
        int maxMz = (maxKz-minKz)+1;
        
        
        InfoMat = new double*[maxMx];
        for (int mx=0;mx<maxMx;++mx){
            InfoMat[mx] = new double[maxMz];
        }
        
        double** InMat = MatIn(StateInfoFile_);
        int cnt=0;
        
        for (int mx=0;mx<maxMx;++mx){
            for (int mz=0;mz<maxMz;++mz){
                InfoMat[mx][mz] = InMat[cnt][2];
                cnt += 1;           
            }
        }
        
        
        for (int i=0;i<maxMx*maxMz;++i){
            delete[] InMat[i];
        }
        
        
        delete[] InMat;
        return InfoMat;
        //delete[] InfoMat;
        
    }
    
    
    // returns 3D array which is used to store controller states at each kx,kz pair
    double*** Controller::ConStates() {
        
        int minKx = minKx_;
        int maxKx = maxKx_;
        int minKz = minKz_;
        int maxKz = maxKz_;
        
        int maxMx = (maxKx-minKx)+1;
        int maxMz = (maxKz-minKz)+1;
        
        double** StateInfo = CStateInfo();
        int cnt=0;
        int tmpNS;
        double*** StateMat;
        StateMat = new double**[maxMx];
        for (int mx=0;mx<maxMx;++mx) {
            StateMat[mx] = new double*[maxMz];
            
            for (int mz=0;mz<maxMz;++mz){ 
                tmpNS = StateInfo[mx][mz];
                StateMat[mx][mz] = new double[tmpNS];
                cnt +=1;
            }    
        }
        
        for (int mx=0;mx<maxMx;++mx){
            for (int mz=0;mz<maxMz;++mz) {
                int tmp = StateInfo[mx][mz];
                for (int i=0;i<tmp;++i){
                    StateMat[mx][mz][i] = 0.0;            
                }
            }        
        }
        
        
        for (int i=0;i<maxMx;++i){
            delete[] StateInfo[i];
        }
        
        
        delete[] StateInfo;
        return StateMat;
        
        //delete[] StateMat;
    }
    
    // Returns kx wavenumber based on controlled mx_c mode number (which differs from sim mx mode number)
    int Controller::kx_c(int mx_c) {
        
        int minKx = minKx_;
        int maxKx = maxKx_;
        int maxMx_c = (maxKx-minKx)+1;
        int kx_c;
        
        if (mx_c >= 0 && mx_c <= maxKx) {
                kx_c = mx_c;               
        }
        else {
                kx_c = mx_c-maxMx_c;
        }

        return kx_c;
        
        
    }
    
    // integrates controller forward in time by dT
    void Controller::advance_Con(FlowField& u, FlowField& q, FlowField& BC, double*** CStateMat, Real dT, Real t, double*** IO) {
        
        int IOcnt = 0;
        int minKx = minKx_;
        int maxKx = maxKx_;
        int minKz = minKz_;
        int maxKz = maxKz_;
        int Ny = Ny_;
        int uvw = uvw_;
        
        int maxMx_c = (maxKx-minKx)+1;
        int maxMz_c = (maxKz-minKz)+1;
        
        
        
        
        double** StateInfo =  CStateInfo();
        
                
        for (int mx_c=0;mx_c<maxMx_c;++mx_c){
            
              
            int kx = kx_c(mx_c);           
                       
            
            
            for(int mz_c=0;mz_c<maxMz_c;++mz_c){ 
                                
                int kz = mz_c + minKz;
                
                                
                if (kx!=0 || kz!=0) {
                        
                    cout << "(" << kx << "," << kz << ") " ;
                    
                    
                    int mx = u.mx(kx); // Simulation mode numbers
                    int mz = u.mz(kz);
                    
                    // Get number of controller states, simulation states and control inputs
                    int numCStates = StateInfo[mx_c][mz_c]; // Controller states   
                    int numSStates = 4*Ny; // Sim states
                    int numInputs = NumConInputs_; // Controller inputs
                    int numOutputs = 4; // Controller outputs - cannot change from 4
                    
                    // Extract number of states from state array     
                    double CStates[numCStates];
                                       
                    for (int i=0;i<numCStates;++i){
                        CStates[i] = CStateMat[mx_c][mz_c][i]; // Controller state matrix
                    }                   

                                       
                    
                    
                    // Read in A,B,C,D and Cs matrices for correct kx,kz pair
                    string KaF = "Mult_Control_Mat/Ka_mat_kx"+i2s(kx)+"kz"+i2s(kz)+".bin";
                    string KbF = "Mult_Control_Mat/Kb_mat_kx"+i2s(kx)+"kz"+i2s(kz)+".bin";
                    string KcF = "Mult_Control_Mat/Kc_mat_kx"+i2s(kx)+"kz"+i2s(kz)+".bin";
                    string KdF = "Mult_Control_Mat/Kd_mat_kx"+i2s(kx)+"kz"+i2s(kz)+".bin";
                    string CsF = "Mult_Control_Mat/Cs_mat_kx"+i2s(kx)+"kz"+i2s(kz)+".bin";
                    const char* KaFile = KaF.c_str(); 
                    const char* KbFile = KbF.c_str();
                    const char* KcFile = KcF.c_str();
                    const char* KdFile = KdF.c_str();
                    const char* CsFile = CsF.c_str();
                    double** Ka = MatIn(KaFile);
                    double** Kb = MatIn(KbFile);
                    double** Kc = MatIn(KcFile);
                    double** Kd = MatIn(KdFile);
                    double** Cs = MatIn(CsFile);



                    
                    // Create sim state vector and find controller input
                    //int Nys = u.Ny(); // num of Cheb Ny points
                    double RSimStates[numSStates]; // real part of states
                    double ISimStates[numSStates]; // imag part of states
                    u.makeState(Spectral, Physical); // Make physical in y
                    q.makeState(Spectral, Physical);
                    
                    for (int i=0;i<Ny;++i){
                        RSimStates[i] =      real(u.cmplx(mx,i,mz,0));
                        RSimStates[i+Ny] =   real(u.cmplx(mx,i,mz,1));
                        RSimStates[i+2*Ny] = real(u.cmplx(mx,i,mz,2));
                        RSimStates[i+3*Ny] = real(q.cmplx(mx,i,mz,0));
                        ISimStates[i] =      imag(u.cmplx(mx,i,mz,0));
                        ISimStates[i+Ny] =   imag(u.cmplx(mx,i,mz,1));
                        ISimStates[i+2*Ny] = imag(u.cmplx(mx,i,mz,2));
                        ISimStates[i+3*Ny] = imag(q.cmplx(mx,i,mz,0));       
                    }      
                    
                    u.makeState(Spectral, Spectral);
                    q.makeState(Spectral, Spectral);
                    
                    

                    
                    // Collect real and imaginary sim states in one vector
                    double SimStates[2*numSStates];
                    for (int i=0;i<numSStates;++i){
                        SimStates[i] = RSimStates[i];
                    }
                    for (int i=numSStates;i<2*numSStates;++i){
                        SimStates[i] = ISimStates[i-numSStates];
                    }
                    
                    
                    // Multiply state vector by Cs to get sim output/controller input    
                    double* uk;   // Controller input vector
                    double* yk;   // Controller output vector
                    yk = new double[4];
                    uk = new double[numInputs];
                    
                    uk = MatVec(Cs,SimStates,numInputs,2*numSStates);                        
                    
                    
                    
                    
                                                 
                    // Integrate Discrete Controller forward in time
                    double* AkXk = MatVec(Ka,CStates,numCStates,numCStates);
                    double* BkUk = MatVec(Kb,uk,numCStates,numInputs);
                    
                    
                    // Update controller states
                    for (int i=0;i<numCStates;++i){
                        CStates[i] = AkXk[i] + BkUk[i];
                    }
                            
                                                                                   
                    
                    
                    
                    // Calculate controller output                    
                    double* CkXk;
                    double* DkUk;                            
                    
                    CkXk = MatVec(Kc,CStates,numOutputs,numCStates);
                    DkUk = MatVec(Kd,uk,numOutputs,numOutputs);
                    for (int i=0; i<numOutputs;++i){
                        yk[i] = CkXk[i] + DkUk[i]; // Calculate controller output
                    }
                                       
                    
                    
                    // Assign controller output to BCs flowfield    
                    double Ad = exp((-1/tau_)*dT);
                    double Bd = -( exp((-1/tau_)*dT) -1);
                    
                    // Simple low-pass filter
                    BC.cmplx(mx,1,mz,uvw) = (Ad*BC.cmplx(mx,1,mz,uvw)) + (Bd*(yk[0] + yk[2]*I));
                    BC.cmplx(mx,0,mz,uvw) = (Ad*BC.cmplx(mx,0,mz,uvw)) + (Bd*(yk[1] + yk[3]*I));
                    
                    
                    // Store updated states in state array    
                    for (int i=0;i<numCStates;++i){
                        CStateMat[mx_c][mz_c][i] = CStates[i];
                    }                
                    
                    
                    if (numInputs == 8){
                        // Sim Output matrix
                        IO[0][0][IOcnt]    = t;
                        IO[0][1][IOcnt]    = kx;
                        IO[0][2][IOcnt]    = kz;        
                        IO[0][3][IOcnt]    = uk[0]; // Real Str upper
                        IO[0][4][IOcnt]    = uk[4]; // Imag Str upper   
                        IO[0][5][IOcnt]    = uk[1]; // Real Str lower   
                        IO[0][6][IOcnt]    = uk[5]; // Imag Str lower  
                        IO[0][7][IOcnt]    = uk[2]; // Real Spa upper 
                        IO[0][8][IOcnt]    = uk[6]; // Imag Spa upper 
                        IO[0][9][IOcnt]    = uk[3]; // Real Spa lower 
                        IO[0][10][IOcnt]    = uk[7]; // Imag Spa lower 
                        
                        // Sim Input matrix
                        IO[1][0][IOcnt]    = t;
                        IO[1][1][IOcnt]    = kx;
                        IO[1][2][IOcnt]    = kz;        
                        IO[1][3][IOcnt]    = real(BC.cmplx(mx,1,mz,uvw));//yk[0]; // Real upper
                        IO[1][4][IOcnt]    = imag(BC.cmplx(mx,1,mz,uvw));//yk[2]; // Imag upper   
                        IO[1][5][IOcnt]    = real(BC.cmplx(mx,0,mz,uvw));//yk[1]; // Real lower   
                        IO[1][6][IOcnt]    = imag(BC.cmplx(mx,0,mz,uvw));//yk[3]; // Imag lower 
                        IO[1][7][IOcnt]    = 0.0;
                        IO[1][8][IOcnt]    = 0.0;
                        IO[1][9][IOcnt]    = 0.0;
                        IO[1][10][IOcnt]    = 0.0;

                        
                    }
                    else {
                        // Sim Output matrix
                        IO[0][0][IOcnt]    = t;
                        IO[0][1][IOcnt]    = kx;
                        IO[0][2][IOcnt]    = kz;        
                        IO[0][3][IOcnt]    = uk[0]; // Real Str upper
                        IO[0][4][IOcnt]    = uk[2]; // Imag Str upper   
                        IO[0][5][IOcnt]    = uk[1]; // Real Str lower   
                        IO[0][6][IOcnt]    = uk[3]; // Imag Str lower
                        
                        // Sim Input matrix
                        IO[1][0][IOcnt]    = t;
                        IO[1][1][IOcnt]    = kx;
                        IO[1][2][IOcnt]    = kz;        
                        IO[1][3][IOcnt]    = real(BC.cmplx(mx,1,mz,uvw));//yk[0]; // Real upper
                        IO[1][4][IOcnt]    = imag(BC.cmplx(mx,1,mz,uvw));//yk[2]; // Imag upper   
                        IO[1][5][IOcnt]    = real(BC.cmplx(mx,0,mz,uvw));//yk[1]; // Real lower   
                        IO[1][6][IOcnt]    = imag(BC.cmplx(mx,0,mz,uvw));//yk[3]; // Imag lower  
                    }
                    
                    
                    
                    
                    
                    IOcnt += 1;   
                    
                    for (int i=0;i<numCStates;++i){                       
                        
                        delete[] Ka[i];
                        //delete[] eyeA[i];
                        //delete[] alpha[i];
                        //delete[] beta[i];  
                        delete[] Kb[i];

                        
                    }
                    
                    for (int i=0;i<4;++i){
                        delete[] Kc[i];
                        delete[] Kd[i];                       
                        
                    }
                    
                    for (int i=0;i<numInputs;++i){
                        delete[] Cs[i];
                        
                    }
                    
                    
                    
                    
                    delete[] Ka;
                    delete[] Kb;
                    delete[] Kc;
                    delete[] Kd;
                    delete[] Cs;
                    //delete[] alpha;
                    //delete[] beta;
                    //delete[] bXk;
                    delete[] BkUk;
                    //delete[] Gam;
                    //delete[] eyeA;
                    delete[] CkXk;
                    delete[] DkUk;  
                    delete[] uk;
                    delete[] yk;
                    
                    
                    }  
                
                
                else {
                    IO[0][0][IOcnt]    = t;
                    IO[0][1][IOcnt]    = kx;
                    IO[0][2][IOcnt]    = kz;        
                    IO[0][3][IOcnt]    = 0; // Real upper
                    IO[0][4][IOcnt]    = 0; // Imag upper   
                    IO[0][5][IOcnt]    = 0; // Real lower   
                    IO[0][6][IOcnt]    = 0; // Imag lower  
                    
                    IO[1][0][IOcnt]    = t;
                    IO[1][1][IOcnt]    = kx;
                    IO[1][2][IOcnt]    = kz;        
                    IO[1][3][IOcnt]    = 0; // Real upper
                    IO[1][4][IOcnt]    = 0; // Imag upper   
                    IO[1][5][IOcnt]    = 0; // Real lower   
                    IO[1][6][IOcnt]    = 0; // Imag lower   
                    
                    IOcnt += 1;          
                    
                    
                }
                
            }
        
        } 
        cout << " " << endl;
        cout << "Control applied to WNs kx="+i2s(minKx)+":"+i2s(maxKx)+", kz="+i2s(minKz)+":"+i2s(maxKz) << endl;
        
        
        for (int i=0;i<maxMx_c;++i){                           
            delete[] StateInfo[i];           
            
        }
        
        
        delete[] StateInfo;

        
}
    
    
    

    
    
    
    // integrates controller forward in time by dT
  void Controller::advance_Con_CN(FlowField& u, FlowField& q, FlowField& BC, double*** CStateMat, Real dT, Real t, double*** IO) {
        
        int IOcnt = 0;
        int minKx = minKx_;
        int maxKx = maxKx_;
        int minKz = minKz_;
        int maxKz = maxKz_;
        int Ny = Ny_;
        int uvw = uvw_;
        
        int maxMx_c = (maxKx-minKx)+1;
        int maxMz_c = (maxKz-minKz)+1;
        
        bool Spectral_states = Spectral_states_;
        
        
        double** StateInfo =  CStateInfo();
        
        
        for (int mx_c=0;mx_c<maxMx_c;++mx_c){
            
            
            int kx = kx_c(mx_c);           
            
            
            
            for(int mz_c=0;mz_c<maxMz_c;++mz_c){ 
                
                int kz = mz_c + minKz;
                
                
                if (kx!=0 || kz!=0) {
                    
                    //cout << "(" << kx << "," << kz << ") " ;
                    
                    
                    int mx = u.mx(kx); // Simulation mode numbers
                    int mz = u.mz(kz);
                    
                    // Get number of controller states, simulation states and control inputs
                    int numCStates = StateInfo[mx_c][mz_c]; // Controller states   
                    int numSStates = 4*Ny; // Sim states
                    int numInputs = NumConInputs_; // Controller inputs
                    int numOutputs = 4; // Controller outputs - cannot change from 4
                    
                    // Extract number of states from state array     
                    double CStates[numCStates];
                    
                    for (int i=0;i<numCStates;++i){
                        CStates[i] = CStateMat[mx_c][mz_c][i]; // Controller state matrix
                    }                   
                    
                                        
                    
                    // Read in A,B,C,D and Cs matrices for correct kx,kz pair
                    string KaF = "Mult_Control_Mat/Ka_mat_kx"+i2s(kx)+"kz"+i2s(kz)+".bin";
                    string KbF = "Mult_Control_Mat/Kb_mat_kx"+i2s(kx)+"kz"+i2s(kz)+".bin";
                    string KcF = "Mult_Control_Mat/Kc_mat_kx"+i2s(kx)+"kz"+i2s(kz)+".bin";
                    string KdF = "Mult_Control_Mat/Kd_mat_kx"+i2s(kx)+"kz"+i2s(kz)+".bin";
                    string CsF = "Mult_Control_Mat/Cs_mat_kx"+i2s(kx)+"kz"+i2s(kz)+".bin";
                    const char* KaFile = KaF.c_str(); 
                    const char* KbFile = KbF.c_str();
                    const char* KcFile = KcF.c_str();
                    const char* KdFile = KdF.c_str();
                    const char* CsFile = CsF.c_str();
                    double** Ka = MatIn(KaFile);
                    double** Kb = MatIn(KbFile);
                    double** Kc = MatIn(KcFile);
                    double** Kd = MatIn(KdFile);
                    double** Cs = MatIn(CsFile);
                    
                    
                    
                    
                    // Create sim state vector and find controller input
                    //int Nys = u.Ny(); // num of Cheb Ny points
                    double RSimStates[numSStates]; // real part of states
                    double ISimStates[numSStates]; // imag part of states
                   
		     if (Spectral_states){
                        
                        u.makeState(Spectral, Spectral); // Make spectral in y
                        q.makeState(Spectral, Spectral);
                        
                        for (int i=0;i<Ny;++i){
                            RSimStates[i] =      real(u.cmplx(mx,i,mz,0));
                            RSimStates[i+Ny] =   real(u.cmplx(mx,i,mz,1));
                            RSimStates[i+2*Ny] = real(u.cmplx(mx,i,mz,2));
                            RSimStates[i+3*Ny] = real(q.cmplx(mx,i,mz,0));
                            ISimStates[i] =      imag(u.cmplx(mx,i,mz,0));
                            ISimStates[i+Ny] =   imag(u.cmplx(mx,i,mz,1));
                            ISimStates[i+2*Ny] = imag(u.cmplx(mx,i,mz,2));
                            ISimStates[i+3*Ny] = imag(q.cmplx(mx,i,mz,0));       
                        }                           
                        
                        
                    }
                    else {
                        
                        u.makeState(Spectral, Physical); // Make physical in y
                        q.makeState(Spectral, Physical);
                        
                        for (int i=0;i<Ny;++i){
                            RSimStates[i] =      real(u.cmplx(mx,i,mz,0));
                            RSimStates[i+Ny] =   real(u.cmplx(mx,i,mz,1));
                            RSimStates[i+2*Ny] = real(u.cmplx(mx,i,mz,2));
                            RSimStates[i+3*Ny] = real(q.cmplx(mx,i,mz,0));
                            ISimStates[i] =      imag(u.cmplx(mx,i,mz,0));
                            ISimStates[i+Ny] =   imag(u.cmplx(mx,i,mz,1));
                            ISimStates[i+2*Ny] = imag(u.cmplx(mx,i,mz,2));
                            ISimStates[i+3*Ny] = imag(q.cmplx(mx,i,mz,0));       
                        }      
                        
                        u.makeState(Spectral, Spectral);
                        q.makeState(Spectral, Spectral);                        
                        
                    }

                  
                    
                    
                    // Collect real and imaginary sim states in one vector
                    double SimStates[2*numSStates];
                    for (int i=0;i<numSStates;++i){
                        SimStates[i] = RSimStates[i];
                    }
                    for (int i=numSStates;i<2*numSStates;++i){
                        SimStates[i] = ISimStates[i-numSStates];
                    }
                    
                    
                    // Multiply state vector by Cs to get sim output/controller input    
                    double* uk;   // Controller input vector
                    double* yk;   // Controller output vector
                    yk = new double[4];
                    uk = new double[numInputs];
                    
                    uk = MatVec(Cs,SimStates,numInputs,2*numSStates);                        
                    
                    
                    
                                                
                    ////////////////////////////////////////////////////////////////
                    
                    // Integrate controller forward in time USING CRANK-NICOLSON
                    
                    // Calculate alpha and beta
                    // alpha = I-delT/2*A, beta = I+delT/2*A
                    double** alpha; //alpha matrix
                    alpha = new double*[numCStates];
                    for (int i=0;i<numCStates;++i){
                        alpha[i] = new double[numCStates];
                    }
                    
                    double** beta; // beta matrix
                    beta = new double*[numCStates];
                    for (int i=0;i<numCStates;++i){
                        beta[i] = new double[numCStates];
                    }
                    
                    double** eyeA = eyeF(numCStates); // Identity matrix                    
                    
                    for (int i=0;i<numCStates;++i){
                        for (int j=0;j<numCStates;++j){
                            alpha[i][j] = eyeA[i][j] - ((dT/2)*Ka[i][j]);
                            beta[i][j] = eyeA[i][j] + ((dT/2)*Ka[i][j]);
                        }
                    }     
                    
                    
                    // Calculate beta*Xk, Bk*uk and Gam    
                    double* bXk; // beta*Xk
                    double* Gam; // bXk + dt*Bkuk
                    double* BkUk; // Bk*uk

                    Gam = new double[numCStates];
                    
                    bXk = MatVec(beta,CStates,numCStates,numCStates);                   
                    BkUk = MatVec(Kb,uk,numCStates,numInputs);
                    
                    
                    // Calculate gamma
                    for (int i=0;i<numCStates;++i){
                        Gam[i] = bXk[i] + (dT*BkUk[i]);
                    }                  
                    
                                        
                    // Calculate new state vector Xk
                    gauss(alpha,Gam,numCStates); // Solves alpha*Xk=Gam
                    // Gam now equal to Xk
                    
                    // Update controller states
                    for (int i=0;i<numCStates;++i){
                        CStates[i] = Gam[i];
                    }
                    ////////////////////////////////////////////////////////////////
                            
                    
                    
                    // Calculate controller output                    
                    double* CkXk;
                    double* DkUk;                            
                    
                    CkXk = MatVec(Kc,CStates,numOutputs,numCStates);
                    DkUk = MatVec(Kd,uk,numOutputs,numOutputs);
                    for (int i=0; i<numOutputs;++i){
                        yk[i] = CkXk[i] + DkUk[i]; // Calculate controller output
                    }
                    
                    
                                        
                    // Crank-Nicolson low-pass filter integration
                    double alpha_lpf = 1 + (dT/(2*tau_));
                    double beta_lpf =  1 - (dT/(2*tau_));
                    
                    BC.cmplx(mx,1,mz,uvw) = (1/alpha_lpf) * (  (beta_lpf*BC.cmplx(mx,1,mz,uvw)) + ((dT/tau_)*(yk[0] + yk[2]*I))  );
                    BC.cmplx(mx,0,mz,uvw) = (1/alpha_lpf) * (  (beta_lpf*BC.cmplx(mx,0,mz,uvw)) + ((dT/tau_)*(yk[1] + yk[3]*I))  );
                    
                                       
                    
                    // Store updated states in state array    
                    for (int i=0;i<numCStates;++i){
                        CStateMat[mx_c][mz_c][i] = CStates[i];
                    }                
                    
		    if (numInputs == 12){
		        // Sim Output matrix
                        IO[0][0][IOcnt]    = t;
                        IO[0][1][IOcnt]    = kx;
                        IO[0][2][IOcnt]    = kz;        
                        IO[0][3][IOcnt]    = uk[0]; // Real Str upper
                        IO[0][4][IOcnt]    = uk[6]; // Imag Str upper   
                        IO[0][5][IOcnt]    = uk[1]; // Real Str lower   
                        IO[0][6][IOcnt]    = uk[7]; // Imag Str lower  
                        IO[0][7][IOcnt]    = uk[2]; // Real Spa upper 
                        IO[0][8][IOcnt]    = uk[8]; // Imag Spa upper 
                        IO[0][9][IOcnt]    = uk[3]; // Real Spa lower 
                        IO[0][10][IOcnt]   = uk[9]; // Imag Spa lower 
                        IO[0][11][IOcnt]   = uk[4]; // Real P upper
                        IO[0][12][IOcnt]   = uk[10]; // Imag P upper
                        IO[0][13][IOcnt]   = uk[5]; // Real P lower
                        IO[0][14][IOcnt]   = uk[11]; // Imag P lower
                     
			// Sim Input matrix
                        IO[1][0][IOcnt]    = t;
                        IO[1][1][IOcnt]    = kx;
                        IO[1][2][IOcnt]    = kz;        
                        IO[1][3][IOcnt]    = real(BC.cmplx(mx,1,mz,uvw));//yk[0]; // Real upper
                        IO[1][4][IOcnt]    = imag(BC.cmplx(mx,1,mz,uvw));//yk[2]; // Imag upper   
                        IO[1][5][IOcnt]    = real(BC.cmplx(mx,0,mz,uvw));//yk[1]; // Real lower   
                        IO[1][6][IOcnt]    = imag(BC.cmplx(mx,0,mz,uvw));//yk[3]; // Imag lower 
                        IO[1][7][IOcnt]    = 0.0;
                        IO[1][8][IOcnt]    = 0.0;
                        IO[1][9][IOcnt]    = 0.0;
                        IO[1][10][IOcnt]    = 0.0;
		    }  
		    else if (numInputs == 8){
                        // Sim Output matrix
                        IO[0][0][IOcnt]    = t;
                        IO[0][1][IOcnt]    = kx;
                        IO[0][2][IOcnt]    = kz;        
                        IO[0][3][IOcnt]    = uk[0]; // Real Str upper
                        IO[0][4][IOcnt]    = uk[4]; // Imag Str upper   
                        IO[0][5][IOcnt]    = uk[1]; // Real Str lower   
                        IO[0][6][IOcnt]    = uk[5]; // Imag Str lower  
                        IO[0][7][IOcnt]    = uk[2]; // Real Spa upper 
                        IO[0][8][IOcnt]    = uk[6]; // Imag Spa upper 
                        IO[0][9][IOcnt]    = uk[3]; // Real Spa lower 
                        IO[0][10][IOcnt]    = uk[7]; // Imag Spa lower 
                        
                        // Sim Input matrix
                        IO[1][0][IOcnt]    = t;
                        IO[1][1][IOcnt]    = kx;
                        IO[1][2][IOcnt]    = kz;        
                        IO[1][3][IOcnt]    = real(BC.cmplx(mx,1,mz,uvw));//yk[0]; // Real upper
                        IO[1][4][IOcnt]    = imag(BC.cmplx(mx,1,mz,uvw));//yk[2]; // Imag upper   
                        IO[1][5][IOcnt]    = real(BC.cmplx(mx,0,mz,uvw));//yk[1]; // Real lower   
                        IO[1][6][IOcnt]    = imag(BC.cmplx(mx,0,mz,uvw));//yk[3]; // Imag lower 
                        IO[1][7][IOcnt]    = 0.0;
                        IO[1][8][IOcnt]    = 0.0;
                        IO[1][9][IOcnt]    = 0.0;
                        IO[1][10][IOcnt]    = 0.0;
                        
                        
                    }
                    else if (numInputs == 4) {
                        // Sim Output matrix
                        IO[0][0][IOcnt]    = t;
                        IO[0][1][IOcnt]    = kx;
                        IO[0][2][IOcnt]    = kz;        
                        IO[0][3][IOcnt]    = uk[0]; // Real Str upper
                        IO[0][4][IOcnt]    = uk[2]; // Imag Str upper   
                        IO[0][5][IOcnt]    = uk[1]; // Real Str lower   
                        IO[0][6][IOcnt]    = uk[3]; // Imag Str lower
                        
                        // Sim Input matrix
                        IO[1][0][IOcnt]    = t;
                        IO[1][1][IOcnt]    = kx;
                        IO[1][2][IOcnt]    = kz;        
                        IO[1][3][IOcnt]    = real(BC.cmplx(mx,1,mz,uvw));//yk[0]; // Real upper
                        IO[1][4][IOcnt]    = imag(BC.cmplx(mx,1,mz,uvw));//yk[2]; // Imag upper   
                        IO[1][5][IOcnt]    = real(BC.cmplx(mx,0,mz,uvw));//yk[1]; // Real lower   
                        IO[1][6][IOcnt]    = imag(BC.cmplx(mx,0,mz,uvw));//yk[3]; // Imag lower  
                    }
		    else {
		        // Sim Input matrix
                        IO[1][0][IOcnt]    = t;
                        IO[1][1][IOcnt]    = kx;
                        IO[1][2][IOcnt]    = kz;        
                        IO[1][3][IOcnt]    = real(BC.cmplx(mx,1,mz,uvw));//yk[0]; // Real upper
                        IO[1][4][IOcnt]    = imag(BC.cmplx(mx,1,mz,uvw));//yk[2]; // Imag upper   
                        IO[1][5][IOcnt]    = real(BC.cmplx(mx,0,mz,uvw));//yk[1]; // Real lower   
                        IO[1][6][IOcnt]    = imag(BC.cmplx(mx,0,mz,uvw));//yk[3]; // Imag lower 
		    }
		       
                    
                    
                    
                    
                    
                    IOcnt += 1;   
                    
                    for (int i=0;i<numCStates;++i){                       
                        
                        delete[] Ka[i];
                        delete[] eyeA[i];
                        delete[] alpha[i];
                        delete[] beta[i];  
                        delete[] Kb[i];
                        
                        
                    }
                    
                    for (int i=0;i<4;++i){
                        delete[] Kc[i];
                        delete[] Kd[i];                       
                        
                    }
                    
                    for (int i=0;i<numInputs;++i){
                        delete[] Cs[i];
                        
                    }
                    
                    
                    
                    
                    delete[] Ka;
                    delete[] Kb;
                    delete[] Kc;
                    delete[] Kd;
                    delete[] Cs;
                    delete[] alpha;
                    delete[] beta;
                    delete[] bXk;
                    delete[] BkUk;
                    delete[] Gam;
                    delete[] eyeA;
                    delete[] CkXk;
                    delete[] DkUk;  
                    delete[] uk;
                    delete[] yk;
                    
                    
                }  
                
                
                else {
                    IO[0][0][IOcnt]    = t;
                    IO[0][1][IOcnt]    = kx;
                    IO[0][2][IOcnt]    = kz;        
                    IO[0][3][IOcnt]    = 0; // Real upper
                    IO[0][4][IOcnt]    = 0; // Imag upper   
                    IO[0][5][IOcnt]    = 0; // Real lower   
                    IO[0][6][IOcnt]    = 0; // Imag lower  
                    
                    IO[1][0][IOcnt]    = t;
                    IO[1][1][IOcnt]    = kx;
                    IO[1][2][IOcnt]    = kz;        
                    IO[1][3][IOcnt]    = 0; // Real upper
                    IO[1][4][IOcnt]    = 0; // Imag upper   
                    IO[1][5][IOcnt]    = 0; // Real lower   
                    IO[1][6][IOcnt]    = 0; // Imag lower   
                    
                    IOcnt += 1;          
                    
                    
                }
                
            }
            
        } 
        //cout << " " << endl;
        //cout << "Control applied to WNs kx="+i2s(minKx)+":"+i2s(maxKx)+", kz="+i2s(minKz)+":"+i2s(maxKz) << endl;
        cout << "C";
        
        
        for (int i=0;i<maxMx_c;++i){                           
            delete[] StateInfo[i];           
            
        }
        
        
        delete[] StateInfo;
        
        
    }


    
    
    
    
    
    
    
    
}
