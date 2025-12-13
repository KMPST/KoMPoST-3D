/*
 * Copyright (c) 2018, Aleksi Kurkela, Aleksas Mazeliauskas, Jean-Francois
 * Paquet, Soeren Schlichting and Derek Teaney
 * All rights reserved.
 *
 * KoMPoST is distributed under MIT license;
 * see the LICENSE file that should be present in the root
 * of the source distribution, or alternately available at:
 * https://github.com/KMPST/KoMPoST/
 */
#include<string> 
#include<iostream>
#include<fstream>
#include<sstream>
#include<omp.h>


// GSL INTEGRATION //
#include <gsl/gsl_integration.h>

#include "GreensFunctions.h"

struct GSLVariables {double r; double s;};

// GSL INTEGRATION MACRO //
#define EVALUATE_GSL_INTEGRAL(Integrand,Variables) \
double Value,Error;\
gsl_integration_workspace *Workspace=gsl_integration_workspace_alloc(NumberOfIntegrationPoints); \
gsl_function F; \
F.function=&Integrand; \
F.params=&Variables; \
gsl_integration_qag(&F,0.0,6.0/Variables.s,1E-3,1E-6,NumberOfIntegrationPoints,2,Workspace,&Value,&Error); \
gsl_integration_workspace_free (Workspace); \
return Value;


//gsl_integration_qag(&F,0.0,10.0/SigmaX,1E-3,1E-3,NumberOfIntegrationPoints,2,Workspace,&Value,&Error);
//gsl_integration_qagiu(&F,0.0,1E-3,1E-3,NumberOfIntegrationPoints,Workspace,&Value,&Error);

// GSL INTERPOLATION //
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

// EVALUATION MACRO FOR GSL INTERPOLATION FUNCTIONS //
#define EVALUATE_GSL_INTERPOLATOR(Interpolator,Value,Accelerator,MinValue,MaxValue)  \
if(Value<MinValue){return gsl_spline_eval(Interpolator,MinValue,Accelerator);} \
else if(Value>MaxValue){ return 0.0;} \
else{return gsl_spline_eval(Interpolator,Value,Accelerator);\
}\

#define EVALUATE_GSL_INTERPOLATOR_2D(Interpolator,xValue,yValue,xAccelerator,yAccelerator,xMinValue,xMaxValue,yMinValue,yMaxValue)  \
if((xValue)<(xMinValue) || (xValue)>(xMaxValue) || (yValue)<(yMinValue) || (yValue)>(yMaxValue)){\
    if((yValue)>=(yMinValue) && (yValue)<=(yMaxValue)){ \
        if((xValue)<(xMinValue)){return gsl_spline2d_eval(Interpolator,xMinValue,yValue,xAccelerator,yAccelerator);} \
        else{return 0.0;} \
    }\
    else if((yValue)>(yMaxValue)){\
        if((xValue)<(xMinValue)){return gsl_spline2d_eval(Interpolator,xMinValue,yMaxValue,xAccelerator,yAccelerator);} \
        else if((xValue)>(xMaxValue)){return 0.0;}\
	    else{return gsl_spline2d_eval(Interpolator,xValue,yMaxValue,xAccelerator,yAccelerator);} \
    } \
    else{\
	if((xValue)<(xMinValue)){return gsl_spline2d_eval(Interpolator,xMinValue,yMinValue,xAccelerator,yAccelerator);} \
	else if((xValue)>(xMaxValue)){return 0.0;} \
	else {return gsl_spline2d_eval(Interpolator,xValue,yMinValue,xAccelerator,yAccelerator);}\
    }\
}\
else{\
	return gsl_spline2d_eval(Interpolator,xValue,yValue,xAccelerator,yAccelerator);\
}\

// BESSEL FUNCTIONS //
#ifndef BesselJ0
#define BesselJ0(x) jn(0,x)
#endif

#ifndef BesselJ1
#define BesselJ1(x) jn(1,x)
#endif

#ifndef BesselJ2
#define BesselJ2(x) jn(2,x)
#endif

#ifndef BesselJ3
#define BesselJ3(x) jn(3,x)
#endif

namespace  GreensFunctions {

    // SMEARING //
    namespace Smearing{
        
        // SMEARING KERNEL IN MOMENTUM SPACE //
        double KernelK(double ks,double Sigma){
            return exp(-0.5*ks*ks*(Sigma*Sigma));
        }
 
    }
    
    
    ////////////////////////////////////////////////
    // GREENS FUNCTIONS FOR ENERGY PERTURBATIONS  //
    ////////////////////////////////////////////////
    
    namespace EnergyPerturbations{
        
        ////////////////////
        // FREE STREAMING //
        ////////////////////
        
        namespace FreeStreaming{
            
            // GREENS FUNCTIONS IN MOMENTUM SPACE //
            #include "ENERGYRESPONSE/MomentumSpaceFreeStreaming.inc"
            
            // PERFORM TRANSFORMATIONS TO COORDINATE SPACE //
            #include "ENERGYRESPONSE/BesselTransform.inc"
            
            // COORDINATE SPACE EXPRESSIONS FOR RESPONSE KERNELS //
            #include "ENERGYRESPONSE/CoordinateSpaceFreeStreaming.inc"
            
            
            // SETUP //
            void Setup(int NumberOfPoints_r, double rMin, double rMax, int NumberOfPoints_eta, double etaMin, double etaMax){
                
                
                CoordinateSpace::Setup(NumberOfPoints_r, rMin, rMax, NumberOfPoints_eta, etaMin, etaMax);
                
                // COMMANDLINE OUTPUT //
                std::cerr << "#FS setup DONE" << std::endl;
            }
            
        }
        
        
        ////////////////////
        // KINETIC THEORY //
        ////////////////////
        
        namespace KineticTheory{
            
            // GREENS FUNCTIONS IN MOMENTUM SPACE //
            #include "ENERGYRESPONSE/MomentumSpaceKineticTheory.inc"
            
            // PERFORM TRANSFORMATIONS TO COORDINATE SPACE //
            #include "ENERGYRESPONSE/BesselTransform.inc"
            
            // COORDINATE SPACE EXPRESSIONS FOR RESPONSE KERNELS //
            #include "ENERGYRESPONSE/CoordinateSpaceKineticTheory.inc"
            
            // SETUP //
            void Setup(int NumberOfPoints_r, double rMin, double rMax, int NumberOfPoints_eta, double etaMin, double etaMax){
                
                // COMMANDLINE OUTPUT //
                std::cerr << "#LOADING IN ENERGY-PERTURBATION PROPAGATOR IN KINETIC THEORY " << std::endl;
                
                // SETUP COORDINATE SPACE INTERPOLATION GRID //
                int NumberOfTimes=173;// int tStart=1; int tOffset=0;
                CoordinateSpace::SetupGrid(NumberOfPoints_r,rMin,rMax,NumberOfPoints_eta,etaMin,etaMax,NumberOfTimes);
                CoordinateSpace::SetValues(NumberOfPoints_r,rMin,rMax,NumberOfPoints_eta,etaMin,etaMax,NumberOfTimes);
              
                // SETUP INTERPOLATORS //
                CoordinateSpace::SetupInterpolators(NumberOfPoints_r,rMin,rMax,NumberOfPoints_eta,etaMin,etaMax,NumberOfTimes);
 
                // COMMANDLINE OUTPUT //
                std::cerr << "#Kinetic set up DONE" << std::endl;
            }
        }

        namespace LowKLimit { 
            // LowKLimit of Green functions
            #include "ENERGYRESPONSE/CoordinateSpaceLowKLimit.inc"
        }
    }
    
    
    //////////////////////////////////////////////////
    // GREENS FUNCTIONS FOR MOMENTUM PERTURBATIONS  //
    //////////////////////////////////////////////////
    
    namespace MomentumPerturbations{
        ////////////////////
        // FREE STREAMING //
        ////////////////////
        
        namespace FreeStreaming{
        
            // GREENS FUNCTIONS IN MOMENTUM SPACE //
            #include "MOMENTUMRESPONSE/MomentumSpaceFreeStreaming.inc"
        
            // PERFORM TRANSFORMATIONS TO COORDINATE SPACE //
            #include "MOMENTUMRESPONSE/BesselTransform.inc"
        
            // COORDINATE SPACE EXPRESSIONS FOR RESPONSE KERNELS //
            #include "MOMENTUMRESPONSE/CoordinateSpaceFreeStreaming.inc"
            
            // SETUP //
            void Setup(int NumberOfPoints_r, double rMin, double rMax, int NumberOfPoints_eta, double etaMin, double etaMax){
                
                CoordinateSpace::Setup(NumberOfPoints_r, rMin, rMax, NumberOfPoints_eta, etaMin, etaMax);
                
                // COMMANDLINE OUTPUT //
                std::cerr << "#DONE" << std::endl;
            }
            
        }
        
        
        namespace KineticTheory{
            
            // GREENS FUNCTIONS IN MOMENTUM SPACE //
            #include "MOMENTUMRESPONSE/MomentumSpaceKineticTheory.inc"
            
            // PERFORM TRANSFORMATIONS TO COORDINATE SPACE //
            #include "MOMENTUMRESPONSE/BesselTransform.inc"
            
            // COORDINATE SPACE EXPRESSIONS FOR RESPONSE KERNELS //
            #include "MOMENTUMRESPONSE/CoordinateSpaceKineticTheory.inc"
            
            // SETUP //
            void Setup(int NumberOfPoints_r, double rMin, double rMax, int NumberOfPoints_eta, double etaMin, double etaMax){
                
                // SETUP COORDINATE SPACE INTERPOLATION GRID //
                //int NumberOfTimes=173; int tStart=1; int tOffset=0;
                
		//CoordinateSpace::SetupGrid(NumberOfPoints_r,rMin,rMax,NumberOfPoints_eta,etaMin,etaMax,NumberOfTimes);
                
		//CoordinateSpace::SetValues(NumberOfPoints_r,rMin,rMax,NumberOfPoints_eta,etaMin,etaMax,NumberOfTimes);

                // SETUP INTERPOLATORS //
                //CoordinateSpace::SetupInterpolators(NumberOfPoints_r,rMin,rMax,NumberOfPoints_eta,etaMin,etaMax,NumberOfTimes);
                
                // COMMANDLINE OUTPUT //
                std::cerr << "#DONE" << std::endl;
            }
        }
    }
    
    ////////////////////////////
    // SETUP RESPONSE KERNELS //
    ////////////////////////////
    
    void Setup(int NumberOfPoints_r,double rMin,double rMax,int NumberOfPoints_eta,double etaMin,double etaMax,int ENERGY_PERTURBATIONS,int MOMENTUM_PERTURBATIONS){
        
        // SETUP ENERGY PERTURBATIONS //
        if(ENERGY_PERTURBATIONS){
            
            //SETUP FREE-STREAMING AND KINETIC THEORY EVOLUTION //
            EnergyPerturbations::FreeStreaming::Setup(NumberOfPoints_r,rMin,rMax,NumberOfPoints_eta,etaMin,etaMax);
            EnergyPerturbations::KineticTheory::Setup(NumberOfPoints_r,rMin,rMax,NumberOfPoints_eta,etaMin,etaMax);

        }
        
        // SETUP MOMENTUM PERTURBATIONS //
        if(MOMENTUM_PERTURBATIONS) {
            
            // SETUP FREE-STREAMING AND KINETIC THEORY EVOLUTION //
            MomentumPerturbations::FreeStreaming::Setup(NumberOfPoints_r,rMin,rMax,NumberOfPoints_eta,etaMin,etaMax);
            MomentumPerturbations::KineticTheory::Setup(NumberOfPoints_r,rMin,rMax,NumberOfPoints_eta,etaMin,etaMax);

        }
    }
    
    ///////////////////////////////
    // OUTPUT RESPONSE FUNCTIONS //
    ///////////////////////////////
    void Output(int ENERGY_PERTURBATIONS,int MOMENTUM_PERTURBATIONS){
        // Sets the final time to print the Greeen functions
        const double EtaOverS0 = 2.0/(4.0 * M_PI);
        double t1 = 100 * pow(EtaOverS0, 3.0/2.0);
        double KValue = M_HBARC ;
        
        // GREENS FUNCTIONS FOR ENERGY PERTURBATIONS //
        if(ENERGY_PERTURBATIONS){
            
            EnergyPerturbations::FreeStreaming::CoordinateSpace::Output("EnergyResponse_FS_X.txt",128,20,0.0,10.0);

            EnergyPerturbations::KineticTheory::CoordinateSpace::Output("EnergyResponse_KoMPoST_X.txt",128,100,0.0,t1,KValue,EtaOverS0) ;
            
        }
        
        // GREENS FUNCTIONS FOR MOMENTUM PERTURBATIONS //
        if(MOMENTUM_PERTURBATIONS) {
            
            MomentumPerturbations::FreeStreaming::CoordinateSpace::Output("MomentumResponse_FS_X.txt",128,20,0.0,10.0);
            MomentumPerturbations::KineticTheory::CoordinateSpace::Output("MomentumResponse_KoMPoST_X.txt",128,100,0.0,t1,KValue,EtaOverS0);
            
        }
    }
}
