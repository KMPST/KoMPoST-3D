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
#include <fstream>
#include <vector>
#include "EventInput.h"
#include "ChargeCurrent.h"
#include "EnergyMomentumTensorIO_hdf5.inc"
// routines for finding the hydro fields from TMUNU
#include "EnergyMomentumDecomposition.inc"
// routines for calculating Navier-Stokes estimate of PIMUNU
#include "ComputeSigma.inc"

// EnergyMomentumTensoIO_music.inc defines read in/out routines loading/writing out
// different components of EnergyMomentumTensorMap data.
// It should be replaced by your own favourite read in/out routine.

// SAVE TMUNU TO FILE //
// File format
//
// for ieta in Neta:
//   for iy in Ns:
//     for ix in Ns:
//         1:x 2:y 3:eta  
//         3:T^{tautau} 4:T^{xx} 5:T^{yy} 6:tau^2 T^{etaeta} 
//         7:-T^{tau x} 8:-T^{tauy} 9:-tau T^{taueta} 
//         10:-T^{xy} 11:-tau T^{yeta} 12:-tau T^{etax} 
//         13-28:CellData[0 -- 15] -- additonal cell data. See KineticEvolution.cpp for definitions
//
// where ix and iy and ieta are cell indecies with grid parameters given by EventInput
// with grid spacing afm in [fm] or deta and grid size [Ns x Ns x Neta]
// x, y given in [fm]
//       
// T^{munu} writen out in units [Gev/fm^3]
void EnergyMomentumTensorMapSave(EnergyMomentumTensorMap *TOut,std::string fname){
    using namespace EventInput ;
    std::ofstream OutStream;
    OutStream.open(fname.c_str());
    
    std::cerr << "#SAVING ENERGY-MOMENTUM TENSOR AT " << TOut->tau << " fm/c TO " << fname << std::endl;
    // Optionally CONVERT  E [GeV^4] TO E [GeV/fm^3] //
    double ScaleFactor= 1./(M_HBARC*M_HBARC*M_HBARC) ;
    
    for(int etaS=etaSTART;etaS<=etaEND;etaS++){
      for(int yS=ySTART;yS<=yEND;yS++){
        for(int xS=xSTART;xS<=xEND;xS++){
        
            // GET RELEVANT COMPONENTS //
            double T00=TOut->Get(0,0,xS,yS,etaS);
            double TXX=TOut->Get(1,1,xS,yS,etaS);
            double TYY=TOut->Get(2,2,xS,yS,etaS);
            double TZZ=TOut->Get(3,3,xS,yS,etaS);
            
            double T0X=TOut->Get(0,1,xS,yS,etaS);
            double T0Y=TOut->Get(0,2,xS,yS,etaS);
            double T0Z=TOut->Get(0,3,xS,yS,etaS);
            
            double TXY=TOut->Get(1,2,xS,yS,etaS);
            double TYZ=TOut->Get(2,3,xS,yS,etaS);
            double TZX=TOut->Get(3,1,xS,yS,etaS);

            double  x = (xS-Ns/2)*afm ;
            double  y = (yS-Ns/2)*afm ;
            double eta=(etaS-Neta/2)*deta ;

            // CREATE OUTPUT //
            OutStream << eta << " "  << x << " " << y << " ";
            OutStream << T00*ScaleFactor << " " << TXX*ScaleFactor << " " << TYY*ScaleFactor << " " << TZZ*ScaleFactor << " ";
            OutStream << T0X*ScaleFactor << " " << T0Y*ScaleFactor << " " << T0Z*ScaleFactor << " ";
            OutStream << TXY*ScaleFactor << " " << TYZ*ScaleFactor << " " << TZX*ScaleFactor << " ";

            // Additional CellData. See KineticEvolution.cpp for definitions 
            for (int id = 0 ; id < EnergyMomentumTensorMap::NCellData; id++) {
                OutStream << TOut->GetCellData(id,xS,yS,etaS) << " " ;
            }
            OutStream << std::endl;
        }//x loop
        OutStream << std::endl;
      }//y loop
      OutStream << std::endl;
    }//eta loop
    
    // CLOSE //
    OutStream.close();
    
    // COMMANDLINE OUTPUT //
    std::cerr << "#map save DONE" << std::endl;
    
}

// LOAD INITIAL SURFACE DENSITIES FROM FILE //
// File format
//
// for ieta in Neta:
//   for ix in Ns:
//     for iy in Ns:
//         1:ix 2:iy 3:ieta 4:tau*e
//         optional 5:tau*n_u 6:tau*n_d 7:tau*n_s
//
// where ix and iy and ieta are cell indecies with grid parameters given by EventInput
// with grid spacing afm in [fm] and grid size [Ns x Ns x Neta]
// tau*e and tau*n_a are surface densities supplied in physical units.
void EnergyMomentumTensorMapLoad(EnergyMomentumTensorMap *TIn,
                                 ChargeCurrentMap *JIn, std::string fname){
    using namespace EventInput ;
        
    // OPEN INPUT FILE //
    std::ifstream InStream;
    InStream.open(fname.c_str());
    if (!InStream.is_open()) {
        std::cerr << "Error: cannot open input file " << fname << std::endl;
        exit(EXIT_FAILURE);
    }
    
    // COMMANDLINE OUTPUT //
    std::cerr << "#GETTING INITIAL SURFACE DENSITIES FROM FILE "
              << fname << std::endl;
    
    // INPUT LINE STRING //
    std::string InputLine;
    
    // SPATIAL GRID INDICES //
    int xS,yS,etaS;
    // INITIAL ENERGY SURFACE DENSITY //
    double tauEd;
    double tauNUp,tauNDown,tauNStrange;
    // CONVERT E [GeV/fm^3] TO  E [GeV^4] //
    double ScaleFactor=(M_HBARC*M_HBARC*M_HBARC);
    
    int InputCount=0;
    bool printfirst=true; 
    // GET DATA LINE BY LINE //
    while(InStream.good()){
        // READ LINE //
        getline(InStream,InputLine);
        
        if (!InputLine.empty() && InputLine[0]!='#') {
            std::stringstream InputValues(InputLine);
            // GET GRID INDICES AND SURFACE DENSITIES //
            InputValues >> xS; InputValues >> yS; InputValues >> etaS;
            InputValues >> tauEd;
            if (InputValues.fail()) {
                std::cerr << "Error: invalid input on line "
                          << InputCount + 1 << std::endl;
                exit(EXIT_FAILURE);
            }
            if (JIn != NULL) {
                InputValues >> tauNUp;
                InputValues >> tauNDown;
                InputValues >> tauNStrange;
                if (InputValues.fail()) {
                    std::cerr << "Error: missing charge columns on input line "
                              << InputCount + 1 << std::endl;
                    exit(EXIT_FAILURE);
                }
            }

            if (xS < 0 || xS >= Ns || yS < 0 || yS >= Ns ||
                etaS < 0 || etaS >= Neta) {
                std::cerr << "Error: input grid index outside the grid on line "
                          << InputCount + 1 << std::endl;
                exit(EXIT_FAILURE);
            }
            const double ed=tauEd/TIn->tau;

            // SAVE ENERGY-MOMENTUM TENSOR TO MEMORY //
            const double small_num=1000*GSL_DBL_MIN;
            if (ed < small_num) {
		        if (printfirst) {
			        std::cerr << "Warning: encountered small energy density T^{\\tau\\tau}<" << small_num << " at (ix,iy)=(" << xS << "," << yS << ")!! Using T^{\\mu\\nu}=diag(" << 10*small_num << "," << 5*small_num << "," << 5*small_num << ",0) instead.\n"; 
			        std::cerr << "Further warnings of small T^{\\tau\\tau} will be suppressed." << std::endl;
			        printfirst=false;
		        }
                    TIn->Set(xS,yS,etaS,10*small_num,5*small_num,5*small_num,0,0,0,0,0,0,0);
            }
            else {
                    const double value=ScaleFactor*ed;
                    TIn->Set(xS,yS,etaS,value,value/2.0,value/2.0,0.0,
                             0.0,0.0,0.0,0.0,0.0,0.0);
            }
            if (JIn != NULL) {
                const double chargeScale = ScaleFactor/JIn->tau;
                JIn->Set(ChargeIndex::Up, xS,yS,etaS,
                         chargeScale*tauNUp, 0.0, 0.0, 0.0);
                JIn->Set(ChargeIndex::Down, xS,yS,etaS,
                         chargeScale*tauNDown, 0.0, 0.0, 0.0);
                JIn->Set(ChargeIndex::Strange, xS,yS,etaS,
                         chargeScale*tauNStrange, 0.0, 0.0, 0.0);
            }
            InputCount++;
        }
    }
    // CHECK //
    if(InputCount!=Ns*Ns*Neta){
        std::cerr << "#ERROR -- COULD NOT LOAD INPUT FILE CORRECTLY" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    // CLOSE //
    InStream.close();
    
    // COMMANDLINE OUTPUT //
    std::cerr << "#map load DONE" << std::endl;
}
    
    
//Take Tmunu from EnergyMomentumTensorMap class and return T^{munu} (hyperbolic coordinates)
void get_Tmunu_raised_hyperbolic(class EnergyMomentumTensorMap * Tmunu, int ix, int iy, int ieta, double Tmunu_raised[4][4]) {
	//tau
	const double tau_in_fm=Tmunu->tau;
	const double tau_in_GeV=tau_in_fm/(M_HBARC);
	const double tau_inv_in_GeV=1./tau_in_GeV;

	const double Ttautau=Tmunu->Get(0,0,ix,iy,ieta);
	const double Ttaux=-1*Tmunu->Get(0,1,ix,iy,ieta);
	const double Ttauy=-1*Tmunu->Get(0,2,ix,iy,ieta);
	const double Ttaueta=-1*tau_inv_in_GeV*Tmunu->Get(0,3,ix,iy,ieta);
	const double Txx=Tmunu->Get(1,1,ix,iy,ieta);
	const double Txy=-1*Tmunu->Get(1,2,ix,iy,ieta);
	const double Txeta=-1*tau_inv_in_GeV*Tmunu->Get(1,3,ix,iy,ieta);
	const double Tyy=Tmunu->Get(2,2,ix,iy,ieta);
	const double Tyeta=-1*tau_inv_in_GeV*Tmunu->Get(2,3,ix,iy,ieta);
	const double Tetaeta=tau_inv_in_GeV*tau_inv_in_GeV*Tmunu->Get(3,3,ix,iy,ieta);

	Tmunu_raised[0][0]=Ttautau;
	Tmunu_raised[0][1]=Ttaux;
	Tmunu_raised[1][0]=Ttaux;
	Tmunu_raised[0][2]=Ttauy;
	Tmunu_raised[2][0]=Ttauy;
	Tmunu_raised[0][3]=Ttaueta;
	Tmunu_raised[3][0]=Ttaueta;
	Tmunu_raised[1][1]=Txx;
	Tmunu_raised[1][2]=Txy;
	Tmunu_raised[2][1]=Txy;
	Tmunu_raised[1][3]=Txeta;
	Tmunu_raised[3][1]=Txeta;
	Tmunu_raised[2][2]=Tyy;
	Tmunu_raised[2][3]=Tyeta;
	Tmunu_raised[3][2]=Tyeta;
	Tmunu_raised[3][3]=Tetaeta;
}

// Get conformal entropy from energy density (given in [GeV^4
double get_conformal_entropy(double energy_density) {
	const double p=energy_density/3.0;
	const double T=pow(p/(16*M_PI*M_PI/90),1./4);	
	const double entropy=(energy_density+p)/T;
	return entropy;
};

// Write initial condition file that can be read by MUSIC
// for ieta in Neta:
//   for ix in Ns:
//     for iy in Ns:
//         1:ieta 2:ix 3:iy  
//         4:e [GeV/fm^3] 5:u^tau 6:u^x 7:u^y 8:u^eta [fm^{-1}] 
//         9: pit^{tautau} 10:pi^{taux} 11:pi^{tauy} 12: pi^{taueta} 
//         13:pi^{xx} 14:pi^{xy} 15:pi^{xeta} 16:pi^{yy} 17:pi^{yeta} 18:pi^{etaeta}
//         19:P_iso
// If charge evolution is enabled, append for u, d, s:
//         n_a V_a^tau V_a^x V_a^y V_a^eta
// giving 34 columns in total.
//
// where ieta, and ix and iy are cell indecies with grid parameters given by EventInput
// with grid spacing afm in [fm] and grid size [Ns x Ns]
//       
//\pi^{\mu\nu} is in units of fermi.
// Note that the loop order is different than in the input file
void write_initial_conditions_MUSIC(std::string outfile_name, bool use_sigmamunu_NavierStokes, double dx, double dy, double deta, int Nx, int Ny, int Neta, \
					EnergyMomentumTensorMap *TOut_Full, EnergyMomentumTensorMap *TOut_BG,
					ChargeCurrentMap *JOut) {
  //void write_initial_conditions_MUSIC(std::string outfile_name, double tau_in_fm, double dx, double dy, int Nx, int Ny, double ** energy_density, double *** flow, double **** shear_tensor, double ** bulk_pressure) {

  std::cerr << "#SAVING " << outfile_name.c_str() << std::endl;

  const bool write_hdf5=KoMPoSTInputs::OutputFormat=="h5";
  std::fstream outfile;
  if (!write_hdf5)
    outfile.open(outfile_name.c_str(),std::fstream::out);
  hid_t hdf5_file=-1;
  hid_t hydro_dataset=-1;
  hid_t charge_dataset=-1;

  //From GeV^4 to fm^-4
  const double shear_tensor_normalisation=1./(M_HBARC*M_HBARC*M_HBARC*M_HBARC);
  const double shear_tensor_normalisation_one_eta=1./(M_HBARC*M_HBARC*M_HBARC*M_HBARC*M_HBARC);
  const double shear_tensor_normalisation_two_eta=1./(M_HBARC*M_HBARC*M_HBARC*M_HBARC*M_HBARC*M_HBARC);
  const double charge_normalisation=1./(M_HBARC*M_HBARC*M_HBARC);
  const double charge_normalisation_one_eta=
      1./(M_HBARC*M_HBARC*M_HBARC*M_HBARC);

  const double tau_in_fm=TOut_Full->tau;
  const double tau_in_GeV=tau_in_fm/M_HBARC;
  if (write_hdf5) {
    int number_of_charges=0;
    if (JOut!=NULL)
      number_of_charges=ChargeIndex::NumberOfCharges;
    hdf5_file=create_hdf5_output_file(outfile_name,tau_in_fm,dx,dy,deta);
    create_hdf5_output_datasets(hdf5_file,Neta,Nx,Ny,number_of_charges,
                                hydro_dataset,charge_dataset);
  }

  double TtautauMax=0.0;
  for(int ieta=0;ieta<Neta;ieta++)
    for(int iy=0;iy<Ny;iy++)
      for(int ix=0;ix<Nx;ix++) {
        const double value=TOut_Full->Get(0,0,ix,iy,ieta);
        if (value>TtautauMax) TtautauMax=value;
      }
  const double TtautauWarn=0.001*TtautauMax;
  const double TtautauFail=0.1*TtautauMax;
  int nMatchingWarnings=0;

  //>> dummy >>dummy >>dummy >>dummy >>neta >> dummy >> nx >> dummy >> ny >> dummy >> deta >> dummy >> dx >> dummy >> dy;
  //# dummy 1 etamax= 1 xmax= 200 ymax= 200 deta= 0 dx= 0.17 dy= 0.17
  if (!write_hdf5) {
    std::stringstream header;
    header << "# tau_in_fm "<<tau_in_fm<<" etamax= "<<Neta<<" xmax= "<<Nx << " ymax= " << Ny << " deta= " << deta <<" dx= " << dx << " dy= " << dy << "\n";
    outfile << header.str().c_str();
  }

  std::vector<float> hydro_slice;
  if (write_hdf5) hydro_slice.resize(Nx*Ny*16);
  std::vector<float> charge_slice;
  if (write_hdf5 && JOut!=NULL)
    charge_slice.resize(Nx*Ny*5*ChargeIndex::NumberOfCharges);

  // Careful, the order of x and y is not the same as in input file
  for(int ieta=0;ieta<Neta;ieta++) {
   double etapos=deta*(ieta-Neta/2);
   for(int ix=0;ix<Nx;ix++) {
    double xpos=dx*(ix-Nx/2);
    for(int iy=0;iy<Ny;iy++) {
      double ypos=dy*(iy-Ny/2);

      // Get Tmunu in a convenient format
      void get_Tmunu_raised_hyperbolic(class EnergyMomentumTensorMap * Tmunu, int ix, int iy, int ieta, double Tmunu_raised[4][4]);
      double Tmunu_Full_raised[4][4];
      get_Tmunu_raised_hyperbolic(TOut_Full, ix, iy, ieta, Tmunu_Full_raised);
      double Tmunu_BG_raised[4][4];
      get_Tmunu_raised_hyperbolic(TOut_BG, ix, iy, ieta, Tmunu_BG_raised);

      // Get the energy density and flow
      //void get_energy_flow_linear_improved(double Tmunu[4][4], double Tmunu_bg[4][4], double tau_in_GeV, double &e, double flow[]);
      //double energy_density;
      //double flow[4];
      //get_energy_flow_linear_improved(Tmunu_Full_raised,Tmunu_BG_raised,tau_in_GeV,energy_density,flow);

      void get_energy_flow_from_Tmunu(double Tmunu_raised[4][4], double tau_in_GeV, double * energy_density, double flow[], enum error_codes * error_code);
      double energy_density;
      double flow[4];
      enum error_codes error;
      get_energy_flow_from_Tmunu(Tmunu_Full_raised,tau_in_GeV,&energy_density,flow,&error);
      if (error != decomp_ok) {
        const double Ttautau=Tmunu_Full_raised[0][0];
        if (Ttautau>=TtautauFail) {
          std::cerr << "Error: Landau matching failed at grid point ("
                    << ix << "," << iy << "," << ieta
                    << "), Ttautau/TtautauMax="
                    << Ttautau/TtautauMax << std::endl;
          exit(EXIT_FAILURE);
        }
        if (Ttautau>=TtautauWarn)
          nMatchingWarnings++;
        TOut_BG->SetCellData(15, ix, iy, ieta, 1);
      }
      else {
        TOut_BG->SetCellData(15, ix, iy, ieta, 0);
      }

      // Get the shear stress tensor and the bulk pressure
      void get_shear_and_bulk_from_Tmunu(double Tmunu_raised[4][4], double tau_in_GeV, double energy_density, double * flow, double pimunu_raised[4][4], double & bulk_pressure);
      double shear_tensor_raised[4][4];
      double bulk_pressure;
      get_shear_and_bulk_from_Tmunu(Tmunu_Full_raised, tau_in_GeV,energy_density, flow, shear_tensor_raised, bulk_pressure);
      const double pressure=energy_density/3.;
      const double P_iso=pressure+bulk_pressure;
      if (use_sigmamunu_NavierStokes) {
        void get_sigma_proper(class EnergyMomentumTensorMap *TOut_Full, class EnergyMomentumTensorMap *TOut_BG, const int ix, const int iy, int ieta, double sigma[4][4]);
        double get_conformal_entropy(double energy_density);
        const double entropy_density=get_conformal_entropy(energy_density);
        get_sigma_proper(TOut_Full, TOut_BG, ix, iy, ieta, shear_tensor_raised);
        for(int mu=0;mu<4;mu++) for(int nu=0;nu<4;nu++) shear_tensor_raised[mu][nu]*=2*entropy_density*KoMPoSTParameters::EtaOverS;
      }
      // Output information
      //profile >> dummy1 >> dummy2 >> dummy3 >> density >> utau >> ux >> uy  >> dummy  >> dummy  >> dummy  >> dummy;
      //energy density should be in GeV/fm^3
      //utau should be in fm^-1 (T^{\tau/x/y} in fm^-4, T^{\tau/x/y \eta} in fm^-5, T^{\eta\eta} in fm^-6)

      double hydro_fields[16];
      hydro_fields[0]=energy_density/(M_HBARC*M_HBARC*M_HBARC);
      hydro_fields[1]=flow[0];
      hydro_fields[2]=flow[1];
      hydro_fields[3]=flow[2];
      hydro_fields[4]=flow[3]/M_HBARC;
      hydro_fields[5]=shear_tensor_raised[0][0]*shear_tensor_normalisation;
      hydro_fields[6]=shear_tensor_raised[0][1]*shear_tensor_normalisation;
      hydro_fields[7]=shear_tensor_raised[0][2]*shear_tensor_normalisation;
      hydro_fields[8]=shear_tensor_raised[0][3]*shear_tensor_normalisation_one_eta;
      hydro_fields[9]=shear_tensor_raised[1][1]*shear_tensor_normalisation;
      hydro_fields[10]=shear_tensor_raised[1][2]*shear_tensor_normalisation;
      hydro_fields[11]=shear_tensor_raised[1][3]*shear_tensor_normalisation_one_eta;
      hydro_fields[12]=shear_tensor_raised[2][2]*shear_tensor_normalisation;
      hydro_fields[13]=shear_tensor_raised[2][3]*shear_tensor_normalisation_one_eta;
      hydro_fields[14]=shear_tensor_raised[3][3]*shear_tensor_normalisation_two_eta;
      hydro_fields[15]=P_iso/(M_HBARC*M_HBARC*M_HBARC);

      //\pi^{\mu\nu} should be in fermi
      //order: pitautau >> pitaux >> pitauy >> pitaueta >> pixx >> pixy >> pixeta >> piyy >> piyeta >> pietaeta
      // P_iso=P+Pi is the isotropic pressure of Tmunu.
      if (write_hdf5) {
        const int offset=(ix*Ny+iy)*16;
        for (int field=0;field<16;field++)
          hydro_slice[offset+field]=(float)hydro_fields[field];
      }
      else {
        outfile << etapos << " " << xpos << " " << ypos << " ";
        for (int field=0;field<16;field++)
          outfile << hydro_fields[field] << " ";
      }
      if (JOut != NULL) {
        double charge_fields[5*ChargeIndex::NumberOfCharges];
        for (int charge = 0; charge < ChargeIndex::NumberOfCharges; charge++) {
          const double density=JOut->GetDensity(charge,ix,iy,ieta,flow);
          charge_fields[5*charge]=density*charge_normalisation;
          for (int mu=0;mu<3;mu++) {
            const double diffusion=JOut->Get(charge,mu,ix,iy,ieta)-density*flow[mu];
            charge_fields[5*charge+1+mu]=diffusion*charge_normalisation;
          }
          const double diffusionEta=JOut->Get(charge,3,ix,iy,ieta)-density*flow[3];
          charge_fields[5*charge+4]=diffusionEta*charge_normalisation_one_eta;
        }
        const int number_of_fields=5*ChargeIndex::NumberOfCharges;
        if (write_hdf5) {
          const int offset=(ix*Ny+iy)*number_of_fields;
          for (int field=0;field<number_of_fields;field++)
            charge_slice[offset+field]=(float)charge_fields[field];
        }
        else
          for (int field=0;field<number_of_fields;field++)
            outfile << charge_fields[field] << " ";
      }
      if (!write_hdf5) outfile << "\n";
    }//y loop
   }// x loop
   if (write_hdf5) {
     hsize_t size[4]={1,(hsize_t)Nx,(hsize_t)Ny,16};
     hsize_t offset[4]={(hsize_t)ieta,0,0,0};
     write_hdf5_slice(hydro_dataset,hydro_slice.data(),4,size,offset);
     if (JOut!=NULL) {
       hsize_t charge_size[5]={1,(hsize_t)Nx,(hsize_t)Ny,
                               ChargeIndex::NumberOfCharges,5};
       hsize_t charge_offset[5]={(hsize_t)ieta,0,0,0,0};
       write_hdf5_slice(charge_dataset,charge_slice.data(),5,
                        charge_size,charge_offset);
     }
   }
  }//eta loop

  if (nMatchingWarnings>0)
    std::cerr << "Warning: Landau matching used the static fallback at "
              << nMatchingWarnings
              << " grid points with Ttautau/TtautauMax >= 0.001."
              << std::endl;

  if (write_hdf5)
    close_hdf5_output(hdf5_file,hydro_dataset,charge_dataset);
  else
    outfile.close();
}
