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
#ifndef EnergyMomentumTensorMap_h
#define EnergyMomentumTensorMap_h

#include "EventInput.h"

// ENERGYMOMENTUMTESORMAP -- data holder of input/output energy-momentum tensor
class EnergyMomentumTensorMap {
    
public:

    // Number of data cells allocated per cell
    static const int NCellData = 16;
    // Number of fluid velocity components
    static const int NUi=3 ;
    
    //////////////////////////////////////////
    // COMPONENTS OF ENERGY MOMENTUM TENSOR //
    //////////////////////////////////////////
    // The covariant and contravariant relation of EnergyMomentum Tensor on the right hand side are checked with gmunu=(1,-1,-1,-1)
    // tau is to make all components of Tmunu have same unit.
    // all response functions should be dimensionless.
    
    //* T00= T^{tau tau}
    //* TXX= T^{xx}
    //* TYY= T^{yy}
    //^ TZZ= 1/t^2 T_{zz} = t^2 T^{zz}
    
    //* T0X=     (F_{ty}F^{y}_{x}+F_{tz}F^{z}_{x})  =        T_{tx} = -T^{tx}
    //* T0Y=     (F_{tz}F^{z}_{y}+F_{tx}F^{x}_{y})  =        T_{ty} = -T^{ty}
    //^ T0Z=  1/t(F_{tx}F^{x}_{z}+F_{ty}F^{y}_{z})  =    1/t T_{tz} = -t T^{tz}
    
    //* TXY=    -(F_{xt}F^{t}_{y}+F_{xz}F^{z}_{y})  =       -T_{xy} = -T^{xy}
    //^ TYZ= -1/t(F_{xt}F^{t}_{z}+F_{xy}F^{y}_{z})  =   -1/t T_{yz} = -t T^{yz}
    //^ TZX= -1/t(F_{zt}F^{t}_{x}+F_{zy}F^{y}_{x})  =   -1/t T_{zx} = -t T^{zx}
    
public:
    
    // EVOLUTION TIME in fm//
    double tau;
    // ENERGY MOMENTUM TENSOR  -- T_{mu nu} in GeV**4 //
    double *T;
    // Additional data about the corresponding fluid state //
    double *CellData;
    // Energy density from reconstructed stress tensor
    double *Ed ;
    // Flow velocity  from reconstructed stress tensor
    double *Ui ;
    
    // INDEXING //
    int Index3D(int x,int y,int eta){
        using EventInput::Ns;
        using EventInput::Neta;
        return x+Ns*y+Ns*Ns*eta;
    }
    // Index of mu nu component of stress tensor
    int Index(int mu,int nu,int xS,int yS,int etaS){
        return mu+4*(nu+4*(Index3D(xS,yS,etaS)));
    }
    // Index of the id data slot of the CellData
    int IndexCellData(int id,int xS,int yS,int etaS){
        return id + NCellData*(Index3D(xS,yS,etaS));
    }
    // Index of the id data slot of the CellData
    int IndexUi(int i,int xS,int yS,int etaS){
        return i + NUi*(Index3D(xS,yS,etaS));
    }
    
    // GET VALUES //
    double Get(int mu,int nu,int xS,int yS,int etaS){
        return T[Index(mu,nu,xS,yS,etaS)];
    }
    // Get Data //
    double GetCellData(int id,int xS,int yS,int etaS){
        return CellData[IndexCellData(id,xS,yS,etaS)];
    }
    // Get Energy Density //
    double GetEd(int xS,int yS,int etaS){
        return Ed[Index3D(xS,yS,etaS)];
    }
    // Get Flow Velocity //
    double GetUi(int i,int xS,int yS,int etaS){
        return Ui[IndexUi(i,xS,yS,etaS)];
    }

    // SET COMPONENT VALUES //
    void SetComponent(int mu,int nu,int xS,int yS,int etaS,double Value){
        T[Index(mu,nu,xS,yS,etaS)]=Value;
    }
    // SET CELLDATA //
    void SetCellData(int id,int xS,int yS,int etaS,double Value){
        CellData[IndexCellData(id,xS,yS,etaS)]=Value;
    }
    // SET Ed //
    void SetEd(int xS,int yS,int etaS,double Value){
        Ed[Index3D(xS,yS,etaS)]=Value;
    }
    // SET Ui //
    void SetUi(int i,int xS,int yS,int etaS,double Value){
        Ui[IndexUi(i,xS,yS,etaS)]=Value;
    }
    // Set the entire stress tensor 
    void Set(int xS,int yS,int etaS,double T00,double TXX,double TYY,double TZZ,double T0X,double T0Y,double T0Z,double TXY,double TYZ,double TZX){
        T[Index(0,0,xS,yS,etaS)]=T00;    T[Index(0,1,xS,yS,etaS)]=T0X;    T[Index(0,2,xS,yS,etaS)]=T0Y;    T[Index(0,3,xS,yS,etaS)]=T0Z;
        T[Index(1,0,xS,yS,etaS)]=T0X;    T[Index(1,1,xS,yS,etaS)]=TXX;    T[Index(1,2,xS,yS,etaS)]=TXY;    T[Index(1,3,xS,yS,etaS)]=TZX;
        T[Index(2,0,xS,yS,etaS)]=T0Y;    T[Index(2,1,xS,yS,etaS)]=TXY;    T[Index(2,2,xS,yS,etaS)]=TYY;    T[Index(2,3,xS,yS,etaS)]=TYZ;
        T[Index(3,0,xS,yS,etaS)]=T0Z;    T[Index(3,1,xS,yS,etaS)]=TZX;    T[Index(3,2,xS,yS,etaS)]=TYZ;    T[Index(3,3,xS,yS,etaS)]=TZZ;
    }
    
    // The following two functions function oppositely
    // Set the stress tensor from a one dimensional array, consisting of the
    // raised components of T^{munu} stored in order
    void SetRaised(int xS, int yS, int etaS, double *tmunu_raised)  {
        double u = tau / M_HBARC ;
        double measure[16] = {1., -1., -1., -u ,
                             -1.,  1., -1., -u ,
                             -1., -1.,  1., -u ,
                             -u , -u , -u , u*u } ;
        for (int id = 0 ; id < 16 ; id++) {
            T[id + 16*Index3D(xS,yS,etaS) ] = tmunu_raised[id]*measure[id] ;
        }
    }

    // Get the components of the stress tensor as a one dimnensional array
    void GetRaised(int xS, int yS, int etaS, double *tmunu_raised) {
        double u = M_HBARC / tau ;
        double measure[16] = {1., -1., -1., -u ,
                             -1.,  1., -1., -u ,
                             -1., -1.,  1., -u ,
                             -u , -u , -u , u*u } ;
        for (int id = 0 ; id < 16 ; id++) {
            tmunu_raised[id] = T[id + 16*Index3D(xS,yS,etaS)] * measure[id] ;
        }
    }

    
    // RESET //
    void Reset(){
      using EventInput::Ns;
      using EventInput::Neta;

      // SET ALL ENTRIES TO ZERO //
      #pragma omp parallel for collapse(2) schedule(static)
      for(int etaS=0;etaS<Neta;etaS++){
        for(int yS=0;yS<Ns;yS++){
            for(int xS=0;xS<Ns;xS++){
                for(int nu=0;nu<4;nu++){
                    for(int mu=0;mu<4;mu++){
                        T[Index(mu,nu,xS,yS,etaS)]=0.0;
                    }
                }
                for(int id=0;id<NCellData;id++){
                    CellData[IndexCellData(id,xS,yS,etaS)]=0.0;
                }
                Ed[Index3D(xS,yS,etaS)]=0.0;
                for(int k=0;k<NUi;k++){
                    Ui[IndexUi(k,xS,yS,etaS)]=0.0;
                }
            }
        }
      }
    }
    // CONSTRUCTOR //
    EnergyMomentumTensorMap(double tau){
        
        using EventInput::Ns;
        using EventInput::Neta;

        // SET TIME //
        this->tau=tau;
        
        // ALLOCATE MEMORY for T//
        this->T=new double[4*4*Ns*Ns*Neta];

        // ALLOCATE MEMORY for Additional CellData //
        this->CellData = new double[NCellData*Ns*Ns*Neta];
        
        // ALLOCATE MEMORY for Energy Density //
        this->Ed = new double[Ns*Ns*Neta];

        // Allocate memory for velocity
        this->Ui = new double[NUi*Ns*Ns*Neta];
        
        // RESET //
        this->Reset();
    }
    
    // DE-STRUCTOR //
    ~EnergyMomentumTensorMap(){
        delete T;
        delete CellData;
        delete Ed;
        delete Ui;
    }
};
#endif
