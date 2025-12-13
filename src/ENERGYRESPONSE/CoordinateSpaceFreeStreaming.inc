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


/////////////////////////////////////////
// COORDINATE SPACE RESPONSE FUNCTIONS //
/////////////////////////////////////////

namespace CoordinateSpace{
    
    //////////////////////////////////////////////
    // EVALUTE RESPONSE KERNEL BY INTERPOLATION //
    //////////////////////////////////////////////
    
    // EVALUATE BETWENN rMin AND rMax //
    double rMin,rMax,etaMin,etaMax,deta;
    
    // GSL INTERPOLATION OBJECTS //
    gsl_interp_accel **GsXAcc,**GvXAcc,**GdXAcc,**GrXAcc;
    gsl_interp_accel **GsZAcc,**GvZAcc,**GdZAcc,**GrZAcc;
    gsl_spline2d *GsInt,*GvInt,*GdInt,*GrInt;
    
    // G00=Gs //
    double GsScalingCurve(double r, double eta){
        
        int tID=omp_get_thread_num();
        EVALUATE_GSL_INTERPOLATOR_2D(GsInt,r,eta,GsXAcc[tID],GsZAcc[tID],rMin,rMax,etaMin,etaMax);
    }
    
    double Gs(double dX,double dETA,double dT){
        return GsScalingCurve(dX/dT, dETA)/(dT*dT);
    }
    
    // G0i=r^{i}/r Gv //
    double GvScalingCurve(double r, double eta){

        int tID=omp_get_thread_num();
        EVALUATE_GSL_INTERPOLATOR_2D(GvInt,r,eta,GvXAcc[tID],GvZAcc[tID],rMin,rMax,etaMin,etaMax);
    }
    
    double Gv(double dX,double dETA,double dT){
        return GvScalingCurve(dX/dT, dETA)/(dT*dT);
    }
    
    // Gij=Gd \delta^{ij} + Gr r^{i}r^{j}/r^2 //
    double GdScalingCurve(double r, double eta){
        
        int tID=omp_get_thread_num();
        EVALUATE_GSL_INTERPOLATOR_2D(GdInt,r,eta,GdXAcc[tID],GdZAcc[tID],rMin,rMax,etaMin,etaMax);
    }
    
    double Gd(double dX,double dETA,double dT){
        return GdScalingCurve(dX/dT, dETA)/(dT*dT);
    }
    
    double GrScalingCurve(double r, double eta){
        
        int tID=omp_get_thread_num();
        EVALUATE_GSL_INTERPOLATOR_2D(GrInt,r,eta,GrXAcc[tID],GrZAcc[tID],rMin,rMax,etaMin,etaMax);
    }
    
    double Gr(double dX,double dETA,double dT){
        return GrScalingCurve(dX/dT, dETA)/(dT*dT);
    }
    
    /////////////////////////////////////////////////////////
    // COMPUTE SMEARED COORDINATE SPACE RESPONSE FUNCTIONS //
    /////////////////////////////////////////////////////////
    
    void Setup(int NumberOfPoints_r, double rMinVal, double rMaxVal, int NumberOfPoints_eta, double etaMinVal, double etaMaxVal){
        
        // SETUP DISCRETE r GRID //
        double *rValues =new double[NumberOfPoints_r];
        double *etaValues =new double[NumberOfPoints_eta];

        double *GsValues=new double[NumberOfPoints_r*NumberOfPoints_eta];
        double *GvValues=new double[NumberOfPoints_r*NumberOfPoints_eta];
        double *GdValues=new double[NumberOfPoints_r*NumberOfPoints_eta];
        double *GrValues=new double[NumberOfPoints_r*NumberOfPoints_eta];
        
        // SET r GRID //
        for(int rIndex=0;rIndex<NumberOfPoints_r;rIndex++){
            double r=rMinVal+(rMaxVal-rMinVal)*(rIndex+1.0)/(NumberOfPoints_r+2);
            rValues[rIndex]=r;
        }
        // SET eta GRID        
        for(int etaIndex=0;etaIndex<NumberOfPoints_eta;etaIndex++){
            double eta=etaMinVal+(etaMaxVal-etaMinVal)*(etaIndex)/(NumberOfPoints_eta-1);
            etaValues[etaIndex]=eta;
        }
       
 
        /////////////////////////
        // SETUP INTERPOLATORS //
        /////////////////////////
        
        // SET BOUNDARY //
        rMin=rValues[0]; rMax=rValues[NumberOfPoints_r-1];
        etaMin=etaValues[0]; etaMax=etaValues[NumberOfPoints_eta-1];deta=(etaMax-etaMin)/(NumberOfPoints_eta-1);
 
        // SETUP GSL INTERPOLATION //
        int NumberOfOpenMPThreads=omp_get_max_threads();
        
        GsXAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];
        GvXAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];
        GdXAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];
        GrXAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];
        
        GsZAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];
        GvZAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];
        GdZAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];
        GrZAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];
        
	#pragma omp parallel for
        for(int i=0;i<NumberOfOpenMPThreads;i++){
            
            GsXAcc[i] = gsl_interp_accel_alloc ();
            GvXAcc[i] = gsl_interp_accel_alloc ();
            GdXAcc[i] = gsl_interp_accel_alloc ();
            GrXAcc[i] = gsl_interp_accel_alloc ();
            
	    GsZAcc[i] = gsl_interp_accel_alloc ();
            GvZAcc[i] = gsl_interp_accel_alloc ();
            GdZAcc[i] = gsl_interp_accel_alloc ();
            GrZAcc[i] = gsl_interp_accel_alloc ();
        }
        
        GsInt=gsl_spline2d_alloc(gsl_interp2d_bilinear,NumberOfPoints_r,NumberOfPoints_eta);
        GvInt=gsl_spline2d_alloc(gsl_interp2d_bilinear,NumberOfPoints_r,NumberOfPoints_eta);
        GdInt=gsl_spline2d_alloc(gsl_interp2d_bilinear,NumberOfPoints_r,NumberOfPoints_eta);
        GrInt=gsl_spline2d_alloc(gsl_interp2d_bilinear,NumberOfPoints_r,NumberOfPoints_eta);

        std::ifstream InStream;
        InStream.open("./EKT/ENERGYRESPONSE/Gs_tot_FS_new.dat");
        std::string InputLine;
        if(InStream.good()){getline(InStream,InputLine);}

//      #pragma omp parallel for
        for(int rIndex=0;rIndex<NumberOfPoints_r;rIndex++){
            for(int etaIndex=0; etaIndex<NumberOfPoints_eta; etaIndex++){
                double Gss,Gsv,Gsd,Gsr,Gseta,Gsseta,Gsveta;
                getline(InStream,InputLine);

                if(!(InputLine.empty())){
                    std::stringstream InputValues(InputLine);
                    InputValues >> Gss; InputValues >> Gsv; InputValues >> Gsd; InputValues >> Gsr;
                    InputValues >> Gseta; InputValues >> Gsseta; InputValues >> Gsveta;

                    gsl_spline2d_set(GsInt, GsValues, rIndex, etaIndex, Gss);
                    gsl_spline2d_set(GvInt, GvValues, rIndex, etaIndex, Gsv);
                    gsl_spline2d_set(GdInt, GdValues, rIndex, etaIndex, Gsd);
                    gsl_spline2d_set(GrInt, GrValues, rIndex, etaIndex, Gsr);

                }
            }
        }
        InStream.close();

        gsl_spline2d_init(GsInt,rValues,etaValues,GsValues,NumberOfPoints_r,NumberOfPoints_eta);
        gsl_spline2d_init(GvInt,rValues,etaValues,GvValues,NumberOfPoints_r,NumberOfPoints_eta);
        gsl_spline2d_init(GdInt,rValues,etaValues,GdValues,NumberOfPoints_r,NumberOfPoints_eta);
        gsl_spline2d_init(GrInt,rValues,etaValues,GrValues,NumberOfPoints_r,NumberOfPoints_eta);
        
        // CLEAN-UP //
        delete[] rValues;
	delete[] etaValues;
        delete[] GsValues;
        delete[] GvValues;
        delete[] GdValues;
        delete[] GrValues;
        
    }
    
    
    
    ///////////////////
    // CREATE OUTPUT //
    ///////////////////
    
    void Output(std::string fname,int NdX,int NdT,double t0,double t1){
        
        std::ofstream Outstream;
        Outstream.open(fname.c_str());
        
        Outstream << "# 1:|x-x_0| 2:|tau-tau0| 3:Gss 4:Gsv 5:Gsd 6:Gsr" << std::endl;
        for(int tIndex=0;tIndex<NdT;tIndex++){
            
            double dT=t0+(tIndex+0.5)*(t1-t0)/(NdT);
            
            for(int xIndex=0;xIndex<NdX;xIndex++){
                
                //double dX=2.0*(xIndex+0.5)*dT/(NdX);
                
                //Outstream << dX << " " << dT  << " " << Gs(dX,dT) << " " << Gv(dX,dT) << " " << Gd(dX,dT) << " " << Gr(dX,dT) << std::endl;
                
            }
            
            Outstream << std::endl;
            Outstream << std::endl;
            
        }
        
        Outstream.close();
    }
    
    
    
}
