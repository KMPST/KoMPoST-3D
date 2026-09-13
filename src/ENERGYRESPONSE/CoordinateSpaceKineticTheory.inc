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

#ifndef M_HBARC
#define M_HBARC  0.197326979
#endif


/////////////////////////////////////////
// COORDINATE SPACE RESPONSE FUNCTIONS //
/////////////////////////////////////////

namespace CoordinateSpace{
    
    //////////////////////////////////////////////
    // EVALUTE RESPONSE KERNEL BY INTERPOLATION //
    //////////////////////////////////////////////
    
    // EVALUATE BETWENN rMin,rMax AND tMin,tMax //
    double rMin,rMax,etaMin,etaMax,deta,tMin,tMax; 
    // GSL INTERPOLATION OBJECTS //
    gsl_interp_accel ***GsTAcc,***GvTAcc,***GdTAcc,***GrTAcc,***GetaTAcc,***GsetaTAcc,***GvetaTAcc;
    gsl_interp_accel ***GsXAcc,***GvXAcc,***GdXAcc,***GrXAcc,***GetaXAcc,***GsetaXAcc,***GvetaXAcc;

    gsl_spline2d **GsInt,**GvInt,**GdInt,**GrInt,**GetaInt,**GsetaInt,**GvetaInt;
    
    // G00=Gs //
    double GsScalingCurve(double r,double t, int etaIndex){
        int tID=omp_get_thread_num();
        EVALUATE_GSL_INTERPOLATOR_2D(GsInt[etaIndex],r,t,GsXAcc[etaIndex][tID],GsTAcc[etaIndex][tID],rMin,rMax,tMin,tMax);
    }
    
    double Gs(double dX,double eta,double dT,double ScalingVarOut){ // ScaleingVarOut is x in 2D-KoMPoST paper.
        double tVal=ScalingVarOut;
        int iETA; double temp,tmp1,tmp2;
	int NumberOfPoints_eta=(etaMax-etaMin)/deta+1;
        if (eta<=etaMin){temp=GsScalingCurve(dX/dT,tVal,0)/(dT*dT);}
	else if (eta>=etaMax){temp=GsScalingCurve(dX/dT,tVal,NumberOfPoints_eta-1)/(dT*dT);}
	else{
            iETA = (int) ( (eta - etaMin )/ deta );
	    tmp1=GsScalingCurve(dX/dT,tVal,iETA  )/(dT*dT); 
	    tmp2=GsScalingCurve(dX/dT,tVal,iETA+1)/(dT*dT); 
            temp = tmp1 + (tmp2-tmp1) * (eta-iETA * deta-etaMin)/deta ;
	}
	return temp;
    }
 
    // G0i=r^{i}/r Gv //
    double GvScalingCurve(double r,double t, int etaIndex){
        int tID=omp_get_thread_num();
        EVALUATE_GSL_INTERPOLATOR_2D(GvInt[etaIndex],r,t,GvXAcc[etaIndex][tID],GvTAcc[etaIndex][tID],rMin,rMax,tMin,tMax);
    }

    double Gv(double dX,double eta,double dT,double ScalingVarOut){
        double tVal=ScalingVarOut;
        int iETA; double temp,tmp1,tmp2;
	int NumberOfPoints_eta=(etaMax-etaMin)/deta+1;
    
        if (eta<=etaMin){temp=GvScalingCurve(dX/dT,tVal,0)/(dT*dT);}
        else if (eta>=etaMax){temp=GvScalingCurve(dX/dT,tVal,NumberOfPoints_eta-1)/(dT*dT);}
        else{
            iETA = (int) ( (eta - etaMin )/ deta );
            tmp1=GvScalingCurve(dX/dT,tVal,iETA  )/(dT*dT);
            tmp2=GvScalingCurve(dX/dT,tVal,iETA+1)/(dT*dT);
            temp = tmp1 + (tmp2-tmp1) * (eta-iETA * deta-etaMin)/deta ;
        }
        return temp;
    }
    
    // Gij=Gd \delta^{ij} + Gr r^{i}r^{j}/r^2 //
    double GdScalingCurve(double r,double t,int etaIndex){
        int tID=omp_get_thread_num();
        EVALUATE_GSL_INTERPOLATOR_2D(GdInt[etaIndex],r,t,GdXAcc[etaIndex][tID],GdTAcc[etaIndex][tID],rMin,rMax,tMin,tMax);
    }

    double Gd(double dX,double eta, double dT,double ScalingVarOut){
        double tVal=ScalingVarOut;
        int iETA; double temp,tmp1,tmp2;
	int NumberOfPoints_eta=(etaMax-etaMin)/deta+1;
    
        if (eta<=etaMin){temp=GdScalingCurve(dX/dT,tVal,0)/(dT*dT);}
        else if (eta>=etaMax){temp=GdScalingCurve(dX/dT,tVal,NumberOfPoints_eta-1)/(dT*dT);}
        else{
            iETA = (int) ( (eta - etaMin )/ deta );
            tmp1=GdScalingCurve(dX/dT,tVal,iETA  )/(dT*dT);
            tmp2=GdScalingCurve(dX/dT,tVal,iETA+1)/(dT*dT);
            temp = tmp1 + (tmp2-tmp1) * (eta-iETA * deta-etaMin)/deta ;
        }
        return temp;
    }

    double GrScalingCurve(double r,double t,int etaIndex){
        int tID=omp_get_thread_num();
        EVALUATE_GSL_INTERPOLATOR_2D(GrInt[etaIndex],r,t,GrXAcc[etaIndex][tID],GrTAcc[etaIndex][tID],rMin,rMax,tMin,tMax);
    }

    double Gr(double dX,double eta,double dT,double ScalingVarOut){
        double tVal=ScalingVarOut;
        int iETA; double temp,tmp1,tmp2;
	int NumberOfPoints_eta=(etaMax-etaMin)/deta+1;
    
        if (eta<=etaMin){temp=GrScalingCurve(dX/dT,tVal,0)/(dT*dT);}
        else if (eta>=etaMax){temp=GrScalingCurve(dX/dT,tVal,NumberOfPoints_eta-1)/(dT*dT);}
        else{
            iETA = (int) ( (eta - etaMin )/ deta );
            tmp1=GrScalingCurve(dX/dT,tVal,iETA  )/(dT*dT);
            tmp2=GrScalingCurve(dX/dT,tVal,iETA+1)/(dT*dT);
            temp = tmp1 + (tmp2-tmp1) * (eta-iETA * deta-etaMin)/deta ;
        }
        return temp;
    }

    // G33=Geta //
    double GetaScalingCurve(double r,double t, int etaIndex){
        int tID=omp_get_thread_num();
        EVALUATE_GSL_INTERPOLATOR_2D(GetaInt[etaIndex],r,t,GetaXAcc[etaIndex][tID],GetaTAcc[etaIndex][tID],rMin,rMax,tMin,tMax);
    }    

    double Geta(double dX,double eta,double dT,double ScalingVarOut){
        double tVal=ScalingVarOut;
        int iETA; double temp,tmp1,tmp2;
	int NumberOfPoints_eta=(etaMax-etaMin)/deta+1;
    
        if (eta<=etaMin){temp=GetaScalingCurve(dX/dT,tVal,0)/(dT*dT);}
        else if (eta>=etaMax){temp=GetaScalingCurve(dX/dT,tVal,NumberOfPoints_eta-1)/(dT*dT);}
        else{
            iETA = (int) ( (eta - etaMin )/ deta );
            tmp1=GetaScalingCurve(dX/dT,tVal,iETA  )/(dT*dT);
            tmp2=GetaScalingCurve(dX/dT,tVal,iETA+1)/(dT*dT);
            temp = tmp1 + (tmp2-tmp1) * (eta-iETA * deta-etaMin)/deta ;
        }
        return temp;
    }
    
    // G03=Gseta //
    double GsetaScalingCurve(double r,double t, int etaIndex){
        int tID=omp_get_thread_num();
        EVALUATE_GSL_INTERPOLATOR_2D(GsetaInt[etaIndex],r,t,GsetaXAcc[etaIndex][tID],GsetaTAcc[etaIndex][tID],rMin,rMax,tMin,tMax);
    }    

    double Gseta(double dX,double eta, double dT,double ScalingVarOut){
        double tVal=ScalingVarOut;
        int iETA; double temp,tmp1,tmp2;
	int NumberOfPoints_eta=(etaMax-etaMin)/deta+1;
    
        if (eta<=etaMin){temp=GsetaScalingCurve(dX/dT,tVal,0)/(dT*dT);}
        else if (eta>=etaMax){temp=GsetaScalingCurve(dX/dT,tVal,NumberOfPoints_eta-1)/(dT*dT);}
        else{
            iETA = (int) ( (eta - etaMin )/ deta );
            tmp1=GsetaScalingCurve(dX/dT,tVal,iETA  )/(dT*dT);
            tmp2=GsetaScalingCurve(dX/dT,tVal,iETA+1)/(dT*dT);
            temp = tmp1 + (tmp2-tmp1) * (eta-iETA * deta-etaMin)/deta ;
        }
        return temp;
    }
    
    // Gi3=Gveta //
    double GvetaScalingCurve(double r,double t, int etaIndex){
        int tID=omp_get_thread_num();
        EVALUATE_GSL_INTERPOLATOR_2D(GvetaInt[etaIndex],r,t,GvetaXAcc[etaIndex][tID],GvetaTAcc[etaIndex][tID],rMin,rMax,tMin,tMax);
    }    

    double Gveta(double dX,double eta, double dT,double ScalingVarOut){
        double tVal=ScalingVarOut;
        int iETA; double temp,tmp1,tmp2;
	int NumberOfPoints_eta=(etaMax-etaMin)/deta+1;
    
        if (eta<=etaMin){temp=GvetaScalingCurve(dX/dT,tVal,0)/(dT*dT);}
        else if (eta>=etaMax){temp=GvetaScalingCurve(dX/dT,tVal,NumberOfPoints_eta-1)/(dT*dT);}
        else{
            iETA = (int) ( (eta - etaMin )/ deta );
            tmp1=GvetaScalingCurve(dX/dT,tVal,iETA  )/(dT*dT);
            tmp2=GvetaScalingCurve(dX/dT,tVal,iETA+1)/(dT*dT);
            temp = tmp1 + (tmp2-tmp1) * (eta-iETA * deta-etaMin)/deta ;
        }
        return temp;
    }

    /////////////////////////////////////////////////////////
    // COMPUTE SMEARED COORDINATE SPACE RESPONSE FUNCTIONS //
    /////////////////////////////////////////////////////////
    
    double *rValues,*etaValues,*tValues;

    double **GsValues;
    double **GvValues;
    double **GdValues;
    double **GrValues;
    double **GetaValues;
    double **GsetaValues;
    double **GvetaValues;
    
    void SetupGrid(int NumberOfPoints_r, double rMinVal, double rMaxVal, int NumberOfPoints_eta,double etaMinVal,double etaMaxVal,int NumberOfTimes){
        
        // ALLOCATE //
        rValues=new double[NumberOfPoints_r];
        etaValues=new double[NumberOfPoints_eta];
        tValues=new double[NumberOfTimes];
	// SET GValues   
    GsValues   =new double*[NumberOfPoints_eta];
    GvValues   =new double*[NumberOfPoints_eta];
    GdValues   =new double*[NumberOfPoints_eta];
    GrValues   =new double*[NumberOfPoints_eta];
    GetaValues =new double*[NumberOfPoints_eta];
    GsetaValues=new double*[NumberOfPoints_eta];
    GvetaValues=new double*[NumberOfPoints_eta];     
	for(int etaIndex=0; etaIndex<NumberOfPoints_eta; etaIndex++){
            GsValues[etaIndex]=new double[NumberOfTimes*NumberOfPoints_r];
            GvValues[etaIndex]=new double[NumberOfTimes*NumberOfPoints_r];
            GdValues[etaIndex]=new double[NumberOfTimes*NumberOfPoints_r];
            GrValues[etaIndex]=new double[NumberOfTimes*NumberOfPoints_r];
            
	    GetaValues[etaIndex] =new double[NumberOfTimes*NumberOfPoints_r];
            GsetaValues[etaIndex]=new double[NumberOfTimes*NumberOfPoints_r];
            GvetaValues[etaIndex]=new double[NumberOfTimes*NumberOfPoints_r];
	}        
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
      

            GsXAcc   =new gsl_interp_accel**[NumberOfPoints_eta];
            GvXAcc   =new gsl_interp_accel**[NumberOfPoints_eta];
            GdXAcc   =new gsl_interp_accel**[NumberOfPoints_eta];
            GrXAcc   =new gsl_interp_accel**[NumberOfPoints_eta];
            GetaXAcc =new gsl_interp_accel**[NumberOfPoints_eta];
            GsetaXAcc=new gsl_interp_accel**[NumberOfPoints_eta];
            GvetaXAcc=new gsl_interp_accel**[NumberOfPoints_eta];
            
            GsTAcc   =new gsl_interp_accel**[NumberOfPoints_eta];
            GvTAcc   =new gsl_interp_accel**[NumberOfPoints_eta];
            GdTAcc   =new gsl_interp_accel**[NumberOfPoints_eta];
            GrTAcc   =new gsl_interp_accel**[NumberOfPoints_eta];
            GetaTAcc =new gsl_interp_accel**[NumberOfPoints_eta];
            GsetaTAcc=new gsl_interp_accel**[NumberOfPoints_eta];
            GvetaTAcc=new gsl_interp_accel**[NumberOfPoints_eta];


	// SETUP GSL INTERPOLATION ACCELARATOR //
        int NumberOfOpenMPThreads=omp_get_max_threads();
        for(int etaIndex=0; etaIndex<NumberOfPoints_eta; etaIndex++){ 
            GsXAcc[etaIndex]=new gsl_interp_accel*[NumberOfOpenMPThreads];
            GvXAcc[etaIndex]=new gsl_interp_accel*[NumberOfOpenMPThreads];
            GdXAcc[etaIndex]=new gsl_interp_accel*[NumberOfOpenMPThreads];
            GrXAcc[etaIndex]=new gsl_interp_accel*[NumberOfOpenMPThreads];
            GetaXAcc[etaIndex] =new gsl_interp_accel*[NumberOfOpenMPThreads];
            GsetaXAcc[etaIndex]=new gsl_interp_accel*[NumberOfOpenMPThreads];
            GvetaXAcc[etaIndex]=new gsl_interp_accel*[NumberOfOpenMPThreads];
            
            GsTAcc[etaIndex]=new gsl_interp_accel*[NumberOfOpenMPThreads];
            GvTAcc[etaIndex]=new gsl_interp_accel*[NumberOfOpenMPThreads];
            GdTAcc[etaIndex]=new gsl_interp_accel*[NumberOfOpenMPThreads];
            GrTAcc[etaIndex]=new gsl_interp_accel*[NumberOfOpenMPThreads];
            GetaTAcc[etaIndex] =new gsl_interp_accel*[NumberOfOpenMPThreads];
            GsetaTAcc[etaIndex]=new gsl_interp_accel*[NumberOfOpenMPThreads];
            GvetaTAcc[etaIndex]=new gsl_interp_accel*[NumberOfOpenMPThreads];
        }
        for(int etaIndex=0; etaIndex<NumberOfPoints_eta; etaIndex++){ 
            #pragma omp parallel for
            for(int i=0;i<NumberOfOpenMPThreads;i++){
                GsXAcc[etaIndex][i] = gsl_interp_accel_alloc ();
                GvXAcc[etaIndex][i] = gsl_interp_accel_alloc ();
                GdXAcc[etaIndex][i] = gsl_interp_accel_alloc ();
                GrXAcc[etaIndex][i] = gsl_interp_accel_alloc ();
                GetaXAcc[etaIndex][i]  = gsl_interp_accel_alloc ();
                GsetaXAcc[etaIndex][i] = gsl_interp_accel_alloc ();
                GvetaXAcc[etaIndex][i] = gsl_interp_accel_alloc ();
                
                GsTAcc[etaIndex][i] = gsl_interp_accel_alloc ();
                GvTAcc[etaIndex][i] = gsl_interp_accel_alloc ();
                GdTAcc[etaIndex][i] = gsl_interp_accel_alloc ();
                GrTAcc[etaIndex][i] = gsl_interp_accel_alloc ();
                GetaTAcc[etaIndex][i]  = gsl_interp_accel_alloc ();
                GsetaTAcc[etaIndex][i] = gsl_interp_accel_alloc ();
                GvetaTAcc[etaIndex][i] = gsl_interp_accel_alloc ();
            }
	}
        // SET INTERPO FUNC
        GsInt   =new gsl_spline2d*[NumberOfPoints_eta];
        GvInt   =new gsl_spline2d*[NumberOfPoints_eta];
        GdInt   =new gsl_spline2d*[NumberOfPoints_eta];
        GrInt   =new gsl_spline2d*[NumberOfPoints_eta];
	GetaInt =new gsl_spline2d*[NumberOfPoints_eta];
        GsetaInt=new gsl_spline2d*[NumberOfPoints_eta];
        GvetaInt=new gsl_spline2d*[NumberOfPoints_eta];
        for(int etaIndex=0; etaIndex<NumberOfPoints_eta; etaIndex++){
            GsInt[etaIndex]=gsl_spline2d_alloc(gsl_interp2d_bilinear,NumberOfPoints_r,NumberOfTimes);
            GvInt[etaIndex]=gsl_spline2d_alloc(gsl_interp2d_bilinear,NumberOfPoints_r,NumberOfTimes);
            GdInt[etaIndex]=gsl_spline2d_alloc(gsl_interp2d_bilinear,NumberOfPoints_r,NumberOfTimes);
            GrInt[etaIndex]=gsl_spline2d_alloc(gsl_interp2d_bilinear,NumberOfPoints_r,NumberOfTimes);
            
	    GetaInt[etaIndex] =gsl_spline2d_alloc(gsl_interp2d_bilinear,NumberOfPoints_r,NumberOfTimes);
            GsetaInt[etaIndex]=gsl_spline2d_alloc(gsl_interp2d_bilinear,NumberOfPoints_r,NumberOfTimes);
            GvetaInt[etaIndex]=gsl_spline2d_alloc(gsl_interp2d_bilinear,NumberOfPoints_r,NumberOfTimes);
        }
    }
    
    void SetValues(int NumberOfPoints_r, double rMinVal, double rMaxVal, int NumberOfPoints_eta,double etaMinVal,double etaMaxVal,int NumberOfTimes){
	std::ifstream InStream_t;
	InStream_t.open("./RTA/ENERGYRESPONSE/tValues.dat");
	std::string InputLine_t;
	getline(InStream_t,InputLine_t);
	for (int tIndex=0;tIndex<NumberOfTimes;tIndex++){
		getline(InStream_t,InputLine_t);
		if(!(InputLine_t.empty())){
		std::stringstream InputValue(InputLine_t);
		double temp;
		InputValue >> temp;
		tValues[tIndex]=temp;
		}
	}	
        InStream_t.close();

        std::ifstream InStream;
        InStream.open("./RTA/ENERGYRESPONSE/Gs_tot_finer.dat");
        std::string InputLine;
        if(InStream.good()){getline(InStream,InputLine);}
    
//	#pragma omp parallel for
        for(int tIndex=0;tIndex<NumberOfTimes;tIndex++){
            for(int rIndex=0;rIndex<NumberOfPoints_r;rIndex++){
                for(int etaIndex=0; etaIndex<NumberOfPoints_eta; etaIndex++){
                        
                    double Gss,Gsv,Gsd,Gsr,Gseta,Gsseta,Gsveta;
                    getline(InStream,InputLine);
                    
                    if(!(InputLine.empty())){
                        std::stringstream InputValues(InputLine);
                        InputValues >> Gss; InputValues >> Gsv; InputValues >> Gsd; InputValues >> Gsr;
                        InputValues >> Gseta; InputValues >> Gsseta; InputValues >> Gsveta;

            		gsl_spline2d_set(GsInt[etaIndex], GsValues[etaIndex], rIndex, tIndex, Gss);
            		gsl_spline2d_set(GvInt[etaIndex], GvValues[etaIndex], rIndex, tIndex, Gsv);
            		gsl_spline2d_set(GdInt[etaIndex], GdValues[etaIndex], rIndex, tIndex, Gsd);
            		gsl_spline2d_set(GrInt[etaIndex], GrValues[etaIndex], rIndex, tIndex, Gsr);

            		gsl_spline2d_set(GetaInt[etaIndex],  GetaValues[etaIndex],  rIndex, tIndex, Gseta);
            		gsl_spline2d_set(GsetaInt[etaIndex], GsetaValues[etaIndex], rIndex, tIndex, Gsseta);
            		gsl_spline2d_set(GvetaInt[etaIndex], GvetaValues[etaIndex], rIndex, tIndex, Gsveta);
                    }
                }
            }
        }
        InStream.close();
    }
    
    /////////////////////////
    // SETUP INTERPOLATORS //
    /////////////////////////
    
    void SetupInterpolators(int NumberOfPoints_r, double rMinVal, double rMaxVal, int NumberOfPoints_eta,double etaMinVal,double etaMaxVal,int NumberOfTimes){
        // SET BOUNDARIES //
        rMin=rValues[0]; rMax=rValues[NumberOfPoints_r-1];
        etaMin=etaValues[0]; etaMax=etaValues[NumberOfPoints_eta-1];deta=(etaMax-etaMin)/(NumberOfPoints_eta-1);
        tMin=tValues[0]; tMax=tValues[NumberOfTimes-1];

	for(int etaIndex=0; etaIndex<NumberOfPoints_eta; etaIndex++){        
            gsl_spline2d_init(GsInt[etaIndex],rValues,tValues,GsValues[etaIndex],NumberOfPoints_r,NumberOfTimes);
            gsl_spline2d_init(GvInt[etaIndex],rValues,tValues,GvValues[etaIndex],NumberOfPoints_r,NumberOfTimes);
            gsl_spline2d_init(GdInt[etaIndex],rValues,tValues,GdValues[etaIndex],NumberOfPoints_r,NumberOfTimes);
            gsl_spline2d_init(GrInt[etaIndex],rValues,tValues,GrValues[etaIndex],NumberOfPoints_r,NumberOfTimes);
            
	    gsl_spline2d_init(GetaInt[etaIndex], rValues,tValues,GetaValues[etaIndex], NumberOfPoints_r,NumberOfTimes);
            gsl_spline2d_init(GsetaInt[etaIndex],rValues,tValues,GsetaValues[etaIndex],NumberOfPoints_r,NumberOfTimes);
            gsl_spline2d_init(GvetaInt[etaIndex],rValues,tValues,GvetaValues[etaIndex],NumberOfPoints_r,NumberOfTimes);
        } 
        // CLEAN-UP //
        delete[] rValues;
        delete[] etaValues;
        delete[] tValues;

        for(int etaIndex=0; etaIndex<NumberOfPoints_eta; etaIndex++){        
        	delete[] GsValues[etaIndex];
        	delete[] GvValues[etaIndex];
        	delete[] GdValues[etaIndex];
        	delete[] GrValues[etaIndex];
        	delete[] GetaValues[etaIndex];
        	delete[] GsetaValues[etaIndex];
        	delete[] GvetaValues[etaIndex];
	}
	delete[] GsValues;
        delete[] GvValues;
        delete[] GdValues;
        delete[] GrValues;
        delete[] GetaValues ;
        delete[] GsetaValues;
        delete[] GvetaValues;
    }
    
    
    ///////////////////
    // CREATE OUTPUT //
    ///////////////////
    
    void Output(std::string fname,int NdX,int NdT,double t0,double t1,double K,double EtaOverS0){
        
        std::ofstream Outstream;
        Outstream.open(fname.c_str());
        Outstream << "# tau_R(tau=1fm) = " << pow(K/M_HBARC,2.0/3.0)/EtaOverS0 <<std::endl;
        Outstream << "# 1:|x-x_0| 2:|tau-tau0| 3:Gss 4:Gsv 5:Gsd 6:Gsr" << std::endl;
        
        for(int tIndex=0;tIndex<NdT;tIndex++){
            
            double dT=t0+(tIndex+0.5)*(t1-t0)/(NdT);
            
            for(int xIndex=0;xIndex<NdX;xIndex++){
                
                double dX=2.0*(xIndex+0.5)*dT/(NdX);

                double ScalingVarOut =pow(K*dT/M_HBARC,2.0/3.0)/EtaOverS0 ;
                
                Outstream << dX << " " << dT  << " "
                   << Gs(dX,deta,dT,ScalingVarOut) << " "
                   << Gv(dX,deta,dT,ScalingVarOut) << " "
                   << Gd(dX,deta,dT,ScalingVarOut) << " "
                   << Gr(dX,deta,dT,ScalingVarOut) << std::endl;
                
            }
            
            Outstream << std::endl;
            Outstream << std::endl;
            
        }
        
        Outstream.close();
    }
}
