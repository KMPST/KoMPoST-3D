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
#include "KineticEvolution.h"
#include "ScalingVariable.h"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <sstream>
#include <vector>


// ENERGY CUT-OFF //
static double ENERGY_CUTOFF = 1E-5;

// EVENT INPUT PARAMETERS //
#include "EventInput.h"
#include "EnergyMomentumTensor.h"
// BACKGROUND EVOLUTION //
#include "BackgroundEvolution.h"
// GREENS FUNCTIONS FOR ENERGY-MOMENTUM PERTRUBATIONS //
#include "GreensFunctions.h"

namespace KoMPoST {
//! Store the perturbations from the first pass, so that the they can be used
//! for the second pass.  This is only relevant if KoMPoSTParameters::Regulator is
//! "TwoPass". The relevant information for each cite is stored in
//! TOutBG->CellData  CellData[9 -- 12].
//!
//! The program runs in two passes. In the first pass we simply record the
//! values of the unregulated output. This is used to modify the scaling
//! variable in the second pass
void PrepareKoMPoSTEstimate(EnergyMomentumTensorMap *TOutBG,
                        EnergyMomentumTensorMap *TOutFull) {
  using namespace EventInput;

  for (int etaS = etaSTART; etaS <= etaEND; etaS++) {
   for (int yS = ySTART; yS <= yEND; yS++) {
    for (int xS = xSTART; xS <= xEND; xS++) {
      double t00bg = TOutBG->Get(0, 0, xS, yS, etaS);
      double t00 = TOutFull->Get(0, 0, xS, yS, etaS);
      double dt00 = t00 - t00bg;
      double t0x = TOutFull->Get(0, 1, xS, yS, etaS);
      double t0y = TOutFull->Get(0, 2, xS, yS, etaS);

      // Store the initial estimate for the second pass.
      TOutBG->SetCellData(9,  xS, yS, etaS, t00bg);
      TOutBG->SetCellData(10, xS, yS, etaS, dt00);
      TOutBG->SetCellData(11, xS, yS, etaS, t0x);
      TOutBG->SetCellData(12, xS, yS, etaS, t0y);
    }
   }
  }
}

//! Regulate the size of the perturbations based on how large they are.  This
//! is relevant if KoMPoSTParameters::Regulator = "KoMPoSTAddition"
//!
//! The routine takes the unregulated result (*TOutFull) and the background
//! (*TOutBg), and modifies (*TOutFull) by modifying the size of the
//! perturbations according to the regulating procedure. In this case we simply
//! increase T00 (and also TXX, TYY, TYY) until a discriminant is satisfied.
void RegulateKoMPoSTAddition(EnergyMomentumTensorMap *TOutBG,
                         EnergyMomentumTensorMap *TOutFull) {
  using namespace EventInput;

  for (int etaS = etaSTART; etaS <= etaEND; etaS++) {
   for (int yS = ySTART; yS <= yEND; yS++) {
    for (int xS = xSTART; xS <= xEND; xS++) {
      // Extract the stress
      double t00bg = TOutBG->Get(0, 0, xS, yS, etaS);
      double t00 = TOutFull->Get(0, 0, xS, yS, etaS);
      double txx = TOutFull->Get(1, 1, xS, yS, etaS);
      double tyy = TOutFull->Get(2, 2, xS, yS, etaS);
      double tzz = TOutFull->Get(3, 3, xS, yS, etaS);

      double dt00 = t00 - t00bg;
      double gx = TOutFull->Get(0, 1, xS, yS, etaS);
      double gy = TOutFull->Get(0, 2, xS, yS, etaS);

      // Store the initial values for the stress tensor
      TOutBG->SetCellData(9,  xS, yS, etaS, t00bg);
      TOutBG->SetCellData(10, xS, yS, etaS, dt00);
      TOutBG->SetCellData(11, xS, yS, etaS, gx);
      TOutBG->SetCellData(12, xS, yS, etaS, gy);

      // Construct the discriminant
      double x1 = sqrt(gx * gx + gy * gy) / (t00bg * 4. / 3.);
      double x2 = dt00 / t00bg;
      double z = x1 - 2. / 3. * x2;

      // Implementing a "Breaking Function" so the regulated
      // z satisfies z < zstart + width
      const double zstart = 0.32;
      const double width = 0.15;
      double znew = z;
      if (z > zstart) {
        znew = zstart + width * std::tanh((z - zstart) / width);
      }
      double delta_t00 = 3. / 2. * (x1 - znew) * t00bg - dt00;

      TOutFull->SetComponent(0, 0, xS, yS, etaS, t00bg + dt00 + delta_t00);
      TOutFull->SetComponent(1, 1, xS, yS, etaS, txx + 1. / 3. * delta_t00);
      TOutFull->SetComponent(2, 2, xS, yS, etaS, tyy + 1. / 3. * delta_t00);
      TOutFull->SetComponent(3, 3, xS, yS, etaS, tzz + 1. / 3. * delta_t00);
    }
   }
  }
}

// COMPUTE EVOLUTION OF THE BACKGROUND ENERGY-MOMENTUM TENSOR.
//
// Given the initial condition *TIn we compute an average initial energy
// density in a causal circle set by SigmaBG.
//
// This is used to determine the scaling variable x at the initial time and the
// final time at each point in the  transverse plane. The scaling variable is
// computed by the Scaler class, which knows about eta/s etc.
//
// Finally the average vale of the stress tensor at each point is stored in
// TOutBG at a time tOut.
//
// On output the TOutBG structure contains the background stress tensor.
//
// It also contains additional information about each grid point which is
// stored in the CellData of the TOutBG structure, CellsData(0...8), which can
// be used for diagnostics, but is also essential for propagating the
// perturbations.   In particular the essential information is the average
// input energy density and TXX, and the scaling variable.
//
// The detailed information about what information is stored is documented
// in the code below -- see SetCellData.
void ComputeBackground(bool IsFirstPass, EnergyMomentumTensorMap *TIn,
                       EnergyMomentumTensorMap *TOutBG,
                       const ScalingVariable &ScalerIn, double SigmaBG, double SigmaBG_eta,
                       int EVOLUTION_MODE) {

  using namespace EventInput;

  // GET RELEVANT TIMES //
  double tIn = TIn->tau;
  double tOut = TOutBG->tau;

  // COMMANDLINE OUTPUT //
  std::cerr << "#COMPUTING EVOLUTION FROM " << tIn << " fm/c TO " << tOut
            << " fm/c" << std::endl;

  const double sigma_r = SigmaBG;
  const double sigma_eta = SigmaBG_eta;
  const double nsigma = 4.0;
  const int range_r = int(nsigma * sigma_r / afm) + 1;
  const int range_eta =
      EVOLUTION_MODE == 0 ? 0 : int(nsigma * sigma_eta / deta) + 1;
  const int kernel_r_size = 2 * range_r + 1;

  std::vector<double> weight_r_kernel(kernel_r_size);
  for (int d = -range_r; d <= range_r; d++) {
    weight_r_kernel[d + range_r] =
        afm * exp(-d * d * afm * afm / (sigma_r * sigma_r)) /
        (sqrt(M_PI) * sigma_r);
  }

  std::vector<double> normalization_r(Ns, 0.0);
  for (int i = 0; i < Ns; i++) {
    const int begin = std::max(i - range_r, 0);
    const int end = std::min(i + range_r, Ns);
    for (int j = begin; j < end; j++) {
      normalization_r[i] += weight_r_kernel[j - i + range_r];
    }
  }

  std::vector<double> weight_eta_kernel(2 * range_eta + 1);
  for (int de = -range_eta; de <= range_eta; de++) {
    const double etasquare = de * de * deta * deta;
    weight_eta_kernel[de + range_eta] =
        deta * exp(-etasquare / (sigma_eta * sigma_eta)) /
        sqrt(M_PI) / sigma_eta;
  }

  std::vector<double> normalization_eta(Neta, 1.0);
  if (EVOLUTION_MODE != 0) {
    for (int i = 0; i < Neta; i++) {
      const int begin = std::max(i - range_eta, 0);
      const int end = std::min(i + range_eta, Neta);
      normalization_eta[i] = 0.0;
      for (int j = begin; j < end; j++) {
        normalization_eta[i] +=
            weight_eta_kernel[j - i + range_eta];
      }
    }
  }

  const long NumberOfCells = long(Ns) * Ns * Neta;
  std::vector<double> T00AvgXY(NumberOfCells);
  std::vector<double> TXXAvgXY(NumberOfCells);
  std::vector<double> TYYAvgXY(NumberOfCells);
  std::vector<double> TZZAvgXY(NumberOfCells);

  {
    std::vector<double> T00AvgX(NumberOfCells);
    std::vector<double> TXXAvgX(NumberOfCells);
    std::vector<double> TYYAvgX(NumberOfCells);
    std::vector<double> TZZAvgX(NumberOfCells);

    #pragma omp parallel for collapse(2) schedule(static)
    for (int etaS = 0; etaS < Neta; etaS++) {
      for (int yS = 0; yS < Ns; yS++) {
        for (int xS = 0; xS < Ns; xS++) {
          double T00 = 0.0;
          double TXX = 0.0;
          double TYY = 0.0;
          double TZZ = 0.0;
          const int begin = std::max(xS - range_r, 0);
          const int end = std::min(xS + range_r, Ns);
          for (int xE = begin; xE < end; xE++) {
            const double weight =
                weight_r_kernel[xE - xS + range_r];
            T00 += weight * TIn->Get(0, 0, xE, yS, etaS);
            TXX += weight * TIn->Get(1, 1, xE, yS, etaS);
            TYY += weight * TIn->Get(2, 2, xE, yS, etaS);
            TZZ += weight * TIn->Get(3, 3, xE, yS, etaS);
          }
          const long index =
              xS + long(Ns) * (yS + long(Ns) * etaS);
          T00AvgX[index] = T00;
          TXXAvgX[index] = TXX;
          TYYAvgX[index] = TYY;
          TZZAvgX[index] = TZZ;
        }
      }
    }

    #pragma omp parallel for collapse(2) schedule(static)
    for (int etaS = 0; etaS < Neta; etaS++) {
      for (int yS = 0; yS < Ns; yS++) {
        for (int xS = 0; xS < Ns; xS++) {
          double T00 = 0.0;
          double TXX = 0.0;
          double TYY = 0.0;
          double TZZ = 0.0;
          const int begin = std::max(yS - range_r, 0);
          const int end = std::min(yS + range_r, Ns);
          for (int yE = begin; yE < end; yE++) {
            const double weight =
                weight_r_kernel[yE - yS + range_r];
            const long input_index =
                xS + long(Ns) * (yE + long(Ns) * etaS);
            T00 += weight * T00AvgX[input_index];
            TXX += weight * TXXAvgX[input_index];
            TYY += weight * TYYAvgX[input_index];
            TZZ += weight * TZZAvgX[input_index];
          }
          const long index =
              xS + long(Ns) * (yS + long(Ns) * etaS);
          T00AvgXY[index] = T00;
          TXXAvgXY[index] = TXX;
          TYYAvgXY[index] = TYY;
          TZZAvgXY[index] = TZZ;
        }
      }
    }
  }

  #pragma omp parallel for collapse(2) schedule(runtime)
  for (int etaS = etaSTART; etaS <= etaEND; etaS++) {
   for (int yS = ySTART; yS <= yEND; yS++) {
    // Class responsible for computing the scaling varaible.  We create a copy
    // for independently running  parallel process.
    ScalingVariable Scaler(ScalerIn);
    for (int xS = xSTART; xS <= xEND; xS++) {
      //////////////////////////
      // BACKGROUND EVOLUTION //
      //////////////////////////

      // ENERGY AVERAGED WITH GAUSSIAN PROFILE  //
      double T00InAvg = 0.0;
      double TXXInAvg = 0.0;
      double TYYInAvg = 0.0;
      double TZZInAvg = 0.0;

      const int etastart =
          EVOLUTION_MODE == 0 ? etaS : std::max(etaS - range_eta, 0);
      const int etaend =
          EVOLUTION_MODE == 0 ? etaS + 1
                              : std::min(etaS + range_eta, Neta);

      for (int etaE = etastart; etaE < etaend; etaE++) {
        const double weight = EVOLUTION_MODE == 0
            ? 1.0 : weight_eta_kernel[etaE - etaS + range_eta];
        const long index =
            xS + long(Ns) * (yS + long(Ns) * etaE);
        T00InAvg += weight * T00AvgXY[index];
        TXXInAvg += weight * TXXAvgXY[index];
        TYYInAvg += weight * TYYAvgXY[index];
        TZZInAvg += weight * TZZAvgXY[index];
      }

      const double Normalization =
          normalization_r[xS] * normalization_r[yS] *
          normalization_eta[etaS];
      // NORMALIZE //
      T00InAvg /= Normalization;
      TXXInAvg /= Normalization;
      TYYInAvg /= Normalization;
      TZZInAvg /= Normalization;

      // EVOLUTION OF BACKGROUND ENERGY DENSITY  //
      double T00BG = 0.0;
      double TXXBG = 0.0;
      double TYYBG = 0.0;
      double TZZBG = 0.0;
      double BackgroundKValue = 0.0;

      // CHECK CUT-OFF CRITERION //
      if (T00InAvg < ENERGY_CUTOFF) {
        T00InAvg = ENERGY_CUTOFF;
        TXXInAvg = 0.5 * ENERGY_CUTOFF;
        TYYInAvg = 0.5 * ENERGY_CUTOFF;
        TZZInAvg = 0.;
      }

      // KINETIC-THEORY //
      if (EVOLUTION_MODE == 1 || EVOLUTION_MODE == 2) {

        if (IsFirstPass) {
          Scaler.ClearEstimate();
        } else {
          // This is the second pass. We use the results from the first pass to
          // modify the second pass.  The relevant results from the first  pass
          // were stored in CellData(9...12) by PrepareKoMPoSTAddition.
          double t00bg = TOutBG->GetCellData(9,  xS, yS, etaS);
          double dt00  = TOutBG->GetCellData(10, xS, yS, etaS);
          double t0x   = TOutBG->GetCellData(11, xS, yS, etaS);
          double t0y   = TOutBG->GetCellData(12, xS, yS, etaS);
          Scaler.SetEstimate(t00bg, dt00, t0x, t0y);
        }

        double ScalingVarIn;
        BackgroundKValue =
            BackgroundEvolution::KineticTheory::DetermineScalingFactor(T00InAvg, tIn, Scaler, ScalingVarIn);

        // ScalingVarBG is an abosolution scaling variable with tIn = 0
        double ScalingVarBG = Scaler.ScalingVar(tOut, BackgroundKValue);
        double EtaByS0 = Scaler.GetEtaOverS0();

        BackgroundEvolution::KineticTheory::Propagate(tOut, BackgroundKValue, ScalingVarBG, EtaByS0, T00BG, TXXBG, TYYBG, TZZBG);

        // The perturbations are evolved with a difference in scaling variables
        // tOut**2/3 - tI**2/3 which is used for the perturbations, and not for
        // the background.
        double ScalingVarOut = Scaler.ScalingVar(tIn, tOut, BackgroundKValue);

        // Store information about the lattice cite in TOutBG's CellData.
        TOutBG->SetCellData(0, xS, yS, etaS, T00InAvg);

        // Store the K value and scaling var in CellData
        TOutBG->SetCellData(1, xS, yS, etaS, BackgroundKValue);

        // Store the scaling variable on output
        TOutBG->SetCellData(2, xS, yS, etaS, ScalingVarOut);

        // The average TXX pressure on in input
        TOutBG->SetCellData(3, xS, yS, etaS, TXXInAvg);

        TOutBG->SetCellData(4, xS, yS, etaS, EtaByS0);
        TOutBG->SetCellData(5, xS, yS, etaS, SigmaBG);

      }
      // FREE-STREAMING 
      else if (EVOLUTION_MODE == 0) {
        // COMPUTE EVOLUTION AND CHECK MATCHING EFFICICENCY //
        BackgroundEvolution::FreeStreaming::Propagate(
            T00InAvg, TXXInAvg, TYYInAvg, TZZInAvg, tIn, tOut, T00BG, TXXBG,
            TYYBG, TZZBG);

        // Store information about the lattice cite in TOutBG's CellData
        // structure -- see above. Some of this information is not relevant for
        // the free streaming case, in which case it is set to zero.
        TOutBG->SetCellData(0, xS, yS, etaS, T00InAvg);

        // Kvalue and scaling var make no sense for free sreaming
        TOutBG->SetCellData(1, xS, yS, etaS, 0.);
        TOutBG->SetCellData(2, xS, yS, etaS, 0.);
        TOutBG->SetCellData(3, xS, yS, etaS, TXXInAvg);
        TOutBG->SetCellData(4, xS, yS, etaS, 0.);
        TOutBG->SetCellData(5, xS, yS, etaS, SigmaBG);
      }
      else {
        std::cerr << "#ERROR -- EVOLUTION MODE NOT SPECIFICED CORRECTLY"
                  << std::endl;
        exit(1);
      }

      TOutBG->Set(xS, yS, etaS, T00BG, TXXBG, TYYBG, TZZBG, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    }// Loop over x
   }// Loop over y
  }// Loop over eta

  std::cerr<< "#BACKGROUND DONE" << std::endl;
}

//! Propagate the perturbations.
//!
//! For a given backgound TOutBG propagate the perturbations and add the
//! perturbations as the come out to TOutBG to produce TOutFull.
void ComputePerturbations(EnergyMomentumTensorMap *TIn,
                          EnergyMomentumTensorMap *TOutBG,
                          EnergyMomentumTensorMap *TOutFull,
                          int ENERGY_PERTURBATIONS, int MOMENTUM_PERTURBATIONS,
                          int EVOLUTION_MODE, double ResponseCutoff_r,
                          double ResponseCutoff_eta) {

  using namespace EventInput;

  // GET RELEVANT TIMES //
  double tIn = TIn->tau;
  double tOut = TOutFull->tau;

  // COMMANDLINE OUTPUT //
  std::cerr << "#COMPUTING EVOLUTION FROM " << tIn << " fm/c TO " << tOut
            << " fm/c" << std::endl;
  // EVOLUTION TIME AND RADIUS OF CAUSAL CIRCLE //
  double EvolutionTime = (tOut - tIn);
  const int range_r = static_cast<int>(std::ceil(ResponseCutoff_r / afm));
  const int range_eta = EVOLUTION_MODE == 0
      ? 0 : static_cast<int>(std::ceil(ResponseCutoff_eta / deta));
  const double response_measure =
      afm * afm * (EVOLUTION_MODE == 0 ? 1.0 : deta);

  #pragma omp parallel for collapse(2) schedule(runtime)
  for (int etaS = etaSTART; etaS <= etaEND; etaS++) {
   for (int yS = ySTART; yS <= yEND; yS++) {
    for (int xS = xSTART; xS <= xEND; xS++) {
      // Get information about the background
      double T00BG = TOutBG->Get(0, 0, xS, yS, etaS);
      double TXXBG = TOutBG->Get(1, 1, xS, yS, etaS);
      double TYYBG = TOutBG->Get(2, 2, xS, yS, etaS);
      double TZZBG = TOutBG->Get(3, 3, xS, yS, etaS);

      // Extract extra-information about the lattice point stored in the
      // background structure TOutBG.
      double T00InAvg = TOutBG->GetCellData(0, xS, yS, etaS);
      double TXXInAvg = TOutBG->GetCellData(3, xS, yS, etaS);
      double ScalingVarOut = TOutBG->GetCellData(2, xS, yS, etaS);

      // The values of the perturbations
      double T00Pert = 0.0;
      double TXXPert = 0.0;
      double TYYPert = 0.0;
      double TZZPert = 0.0;
      double T0XPert = 0.0;
      double T0YPert = 0.0;
      double TXYPert = 0.0;

      double T0ZPert = 0.0;
      double TXZPert = 0.0;
      double TYZPert = 0.0;

      const int xbegin = std::max(xS - range_r, 0);
      const int xend = std::min(xS + range_r, Ns - 1);
      const int ybegin = std::max(yS - range_r, 0);
      const int yend = std::min(yS + range_r, Ns - 1);
      const int etabegin = EVOLUTION_MODE == 0 ? etaS : std::max(etaS - range_eta, 0);
      const int etaend = EVOLUTION_MODE == 0 ? etaS : std::min(etaS + range_eta, Neta - 1);

      // CHECK CUT-OFF CRITERION //
      if (T00InAvg <= ENERGY_CUTOFF) {
        // If below the cutoff ignore the perturbations
        goto PERTURBATION_FINISHUP;
      }

      // CHECK WETHER ENERGY/MOMENTUM PERTURBATIONS ARE INCLUDED //
      if (!(ENERGY_PERTURBATIONS || MOMENTUM_PERTURBATIONS)) {
        // Nothing to do ignore the perturbations
        goto PERTURBATION_FINISHUP;
      }

      for (int etaE = etabegin; etaE <= etaend; etaE++) {
       double DeltaEta = (etaS - etaE) * deta;
       double Distance_eta = std::sqrt(DeltaEta * DeltaEta);
       if (EVOLUTION_MODE != 0 && Distance_eta >= ResponseCutoff_eta) {
         continue;
       }
       double Eta = 0.0;
       if (Distance_eta != 0.0) {
         Eta = DeltaEta / Distance_eta;
       }
       for (int yE = ybegin; yE <= yend; yE++) {
        double DeltaY = (yS - yE) * afm;
        for (int xE = xbegin; xE <= xend; xE++) {
          // GET COORDINATES RELATIVE TO POINT OF INTEREST //
          double DeltaX = (xS - xE) * afm;
          double Distance_r = std::sqrt(DeltaX * DeltaX + DeltaY * DeltaY);

          // CHECK THAT DISTANCES ARE RELEVANT FOR EVOLUTION //
          if (Distance_r >= ResponseCutoff_r) {
            // We are outside the causal circle
            continue;
          }

          double rX = 0.0;
          double rY = 0.0;
          if (Distance_r != 0.0) {
            rX = DeltaX / Distance_r;
            rY = DeltaY / Distance_r;
          }

          // COMPUTE AMPLITUDE OF INITIAL PERTURBATIONS //
          double dT00In = (TIn->Get(0, 0, xE, yE, etaE) - T00InAvg) / T00InAvg;
          double dT0XIn =  TIn->Get(0, 1, xE, yE, etaE) / T00InAvg;
          double dT0YIn =  TIn->Get(0, 2, xE, yE, etaE) / T00InAvg;

          //////////////////////////////////////////////////////
          // COMPUTE ENERGY-MOMENTUM PERTURBATION PROPAGATORS //
          //////////////////////////////////////////////////////

          // GREENS-FUNCTIONS FOR ENERGY PERTURBATIONS //
          double G00_00 = 0.0;
          double G0X_00 = 0.0;
          double G0Y_00 = 0.0;
          double GXX_00 = 0.0;
          double GYY_00 = 0.0;
          double GZZ_00 = 0.0;
          double GXY_00 = 0.0;
          
	  double G0Z_00 = 0.0;
	  double GXZ_00 = 0.0;
	  double GYZ_00 = 0.0;

          if (ENERGY_PERTURBATIONS) {

            // GREENS-FUNCTIONS IN TENSOR BASIS //
            double Gs=0., Gv=0., Gd=0., Gr=0., Geta=0., Gseta=0., Gveta=0.;

            // COMPUTE GREENS FUCNTIONS //

            // Kinetic Theory
            if (EVOLUTION_MODE == 1) {

              Gs = GreensFunctions::EnergyPerturbations::KineticTheory::
                  CoordinateSpace::Gs(Distance_r,  Distance_eta, EvolutionTime, ScalingVarOut);
              Gv = GreensFunctions::EnergyPerturbations::KineticTheory::
                  CoordinateSpace::Gv(Distance_r,  Distance_eta, EvolutionTime, ScalingVarOut);
              Gd = GreensFunctions::EnergyPerturbations::KineticTheory::
                  CoordinateSpace::Gd(Distance_r,  Distance_eta, EvolutionTime, ScalingVarOut);
              Gr = GreensFunctions::EnergyPerturbations::KineticTheory::
                  CoordinateSpace::Gr(Distance_r,  Distance_eta, EvolutionTime, ScalingVarOut);

              Geta  = GreensFunctions::EnergyPerturbations::KineticTheory::
                  CoordinateSpace::Geta(Distance_r,  Distance_eta,EvolutionTime, ScalingVarOut);
              Gseta = GreensFunctions::EnergyPerturbations::KineticTheory::
                  CoordinateSpace::Gseta(Distance_r, Distance_eta, EvolutionTime, ScalingVarOut);
              Gveta = GreensFunctions::EnergyPerturbations::KineticTheory::
                  CoordinateSpace::Gveta(Distance_r, Distance_eta, EvolutionTime, ScalingVarOut);

              // LowK limit of Kinetic Theory
            } else if (EVOLUTION_MODE == 2) {

              GreensFunctions::LowKArguments Lowk(
                  tOut, tIn, T00InAvg, T00BG, TXXInAvg, TXXBG,
                  TOutBG->GetCellData(4, xS, yS, etaS),
                  TOutBG->GetCellData(5, xS, yS, etaS));

              //Gs = GreensFunctions::EnergyPerturbations::LowKLimit::
              //    CoordinateSpace::Gs(Distance, EvolutionTime, Lowk);
              //Gv = GreensFunctions::EnergyPerturbations::LowKLimit::
              //    CoordinateSpace::Gv(Distance, EvolutionTime, Lowk);
              //Gd = GreensFunctions::EnergyPerturbations::LowKLimit::
              //    CoordinateSpace::Gd(Distance, EvolutionTime, Lowk);
              //Gr = GreensFunctions::EnergyPerturbations::LowKLimit::
              //    CoordinateSpace::Gr(Distance, EvolutionTime, Lowk);

              // Free Streaming
            } else if (EVOLUTION_MODE == 0) {

              Gs = GreensFunctions::EnergyPerturbations::FreeStreaming::
                  CoordinateSpace::Gs(Distance_r, EvolutionTime);
              Gv = GreensFunctions::EnergyPerturbations::FreeStreaming::
                  CoordinateSpace::Gv(Distance_r, EvolutionTime);
              Gd = GreensFunctions::EnergyPerturbations::FreeStreaming::
                  CoordinateSpace::Gd(Distance_r, EvolutionTime);
              Gr = GreensFunctions::EnergyPerturbations::FreeStreaming::
                  CoordinateSpace::Gr(Distance_r, EvolutionTime);
	      Geta=0.0;
	      Gseta=0.0;
	      Gveta=0.0;

            } else {
              std::cerr << "#ERROR -- EVOLUTION MODE NOT SPECIFICED CORRECTLY"
                        << std::endl;
              exit(1);
            }

            // ENERGY-MOMENTUM TENSOR RESPONSE //
            if (Distance_r == 0) {// take care of the zero distance
	      rX=0.0;
	      rY=0.0;
	    }
	    if (Distance_eta == 0.0) {
	      Eta=0.0;
	    }

            G00_00 = Gs;

            G0X_00 = rX * Gv;
            G0Y_00 = rY * Gv;

            GXX_00 = Gd + rX * rX * Gr;
            GYY_00 = Gd + rY * rY * Gr;
            GZZ_00 = (G00_00 - GXX_00 - GYY_00);

            GXY_00 = rX * rY * Gr;

	    G0Z_00 = Eta * Gseta;
	    GXZ_00 = rX * Eta * Gveta;
	    GYZ_00 = rY * Eta * Gveta;
          
	  } // If energy perturbations

          // GREENS-FUNCTIONS FOR MOMENTUM PERTRUBATIONS  //
          double G00_0X = 0.0;
          double G0X_0X = 0.0;
          double G0Y_0X = 0.0;
          double GXX_0X = 0.0;
          double GYY_0X = 0.0;
          double GZZ_0X = 0.0;
          double GXY_0X = 0.0;
          double G0Z_0X = 0.0;
          double GXZ_0X = 0.0;
          double GYZ_0X = 0.0;

          double G00_0Y = 0.0;
          double G0X_0Y = 0.0;
          double G0Y_0Y = 0.0;
          double GXX_0Y = 0.0;
          double GYY_0Y = 0.0;
          double GZZ_0Y = 0.0;
          double GXY_0Y = 0.0;
          double G0Z_0Y = 0.0;
          double GXZ_0Y = 0.0;
          double GYZ_0Y = 0.0;

          if (MOMENTUM_PERTURBATIONS) {

            // GREENS FUNCTIONS FOR MOMENTUM PERTURBATIONS //
            //double Hv, Hd, Hr, Htd, Htm, Htr;

            // COMPUTE GREENS FUCNTIONS //
            if (EVOLUTION_MODE == 1) {

              //Hv = GreensFunctions::MomentumPerturbations::KineticTheory::
              //    CoordinateSpace::Hv(Distance, EvolutionTime, ScalingVarOut);
              //Hd = GreensFunctions::MomentumPerturbations::KineticTheory::
              //    CoordinateSpace::Hd(Distance, EvolutionTime, ScalingVarOut);
              //Hr = GreensFunctions::MomentumPerturbations::KineticTheory::
              //    CoordinateSpace::Hr(Distance, EvolutionTime, ScalingVarOut);
              //Htd = GreensFunctions::MomentumPerturbations::KineticTheory::
              //    CoordinateSpace::Htd(Distance, EvolutionTime, ScalingVarOut);
              //Htm = GreensFunctions::MomentumPerturbations::KineticTheory::
              //    CoordinateSpace::Htm(Distance, EvolutionTime, ScalingVarOut);
              //Htr = GreensFunctions::MomentumPerturbations::KineticTheory::
              //    CoordinateSpace::Htr(Distance, EvolutionTime, ScalingVarOut);

            }

            else if (EVOLUTION_MODE == 0) {

              //Hv = GreensFunctions::MomentumPerturbations::FreeStreaming::
              //    CoordinateSpace::Hv(Distance, EvolutionTime);
              //Hd = GreensFunctions::MomentumPerturbations::FreeStreaming::
              //    CoordinateSpace::Hd(Distance, EvolutionTime);
              //Hr = GreensFunctions::MomentumPerturbations::FreeStreaming::
              //    CoordinateSpace::Hr(Distance, EvolutionTime);
              //Htd = GreensFunctions::MomentumPerturbations::FreeStreaming::
              //    CoordinateSpace::Htd(Distance, EvolutionTime);
              //Htm = GreensFunctions::MomentumPerturbations::FreeStreaming::
              //    CoordinateSpace::Htm(Distance, EvolutionTime);
              //Htr = GreensFunctions::MomentumPerturbations::FreeStreaming::
              //    CoordinateSpace::Htr(Distance, EvolutionTime);

            }

            else {
              std::cerr << "#ERROR -- EVOLUTION MODE NOT SPECIFICED CORRECTLY"
                        << std::endl;
              exit(1);
            }

            // ENERGY-MOMENTUM TENSOR RESPONSE //
            //if (Distance_r == 0) {

            //  G00_0X = 0.0;
            //  G00_0Y = 0.0;

            //  G0X_0X = Hd + 0.5 * Hr;
            //  G0Y_0X = 0.0;
            //  G0X_0Y = 0.0;
            //  G0Y_0Y = Hd + 0.5 * Hr;

            //  GXX_0X = 0.0;
            //  GYY_0X = 0.0;
            //  GXX_0Y = 0.0;
            //  GYY_0Y = 0.0;
            //  GZZ_0X = 0.0;
            //  GZZ_0Y = 0.0;

            //  GXY_0X = 0.0;
            //  GXY_0Y = 0.0;

            //}

            //else {

            //  G00_0X = rX * Hv;
            //  G00_0Y = rY * Hv;

            //  G0X_0X = Hd + rX * rX * Hr;
            //  G0Y_0X = rY * rX * Hr;
            //  G0X_0Y = rX * rY * Hr;
            //  G0Y_0Y = Hd + rY * rY * Hr;

            //  GXX_0X = rX * Htd + rX * Htm + rX * rX * rX * Htr;
            //  GYY_0X = rX * Htd + rY * rY * rX * Htr;
            //  GXX_0Y = rY * Htd + rX * rX * rY * Htr;
            //  GYY_0Y = rY * Htd + rY * Htm + rY * rY * rY * Htr;
            //  GZZ_0X = (G00_0X - GXX_0X - GYY_0X);
            //  GZZ_0Y = (G00_0Y - GXX_0Y - GYY_0Y);

            //  GXY_0X = 0.5 * (rY)*Htm + rX * rY * rX * Htr;
            //  GXY_0Y = 0.5 * (rX)*Htm + rX * rY * rY * Htr;
            //}
          } // If momentum perturbations

          // COMPUTE CONTRIBUTION TO ENERGY-MOMENTUM TENSOR AT POINT OF
          // INTEREST //
          T00Pert += response_measure *
                     (+G00_00 * dT00In - G00_0X * dT0XIn - G00_0Y * dT0YIn) *
                     T00BG;
          TXXPert += response_measure *
                     (+GXX_00 * dT00In - GXX_0X * dT0XIn - GXX_0Y * dT0YIn) *
                     T00BG;
          TYYPert += response_measure *
                     (+GYY_00 * dT00In - GYY_0X * dT0XIn - GYY_0Y * dT0YIn) *
                     T00BG;
          TZZPert += response_measure *
                     (+GZZ_00 * dT00In - GZZ_0X * dT0XIn - GZZ_0Y * dT0YIn) *
                     T00BG;

          T0XPert += response_measure *
                     (-G0X_00 * dT00In + G0X_0X * dT0XIn + G0X_0Y * dT0YIn) *
                     T00BG;
          T0YPert += response_measure *
                     (-G0Y_00 * dT00In + G0Y_0X * dT0XIn + G0Y_0Y * dT0YIn) *
                     T00BG;
          TXYPert += response_measure *
                     (-GXY_00 * dT00In + GXY_0X * dT0XIn + GXY_0Y * dT0YIn) *
                     T00BG;

          T0ZPert += response_measure *
                     (-G0Z_00 * dT00In + G0Z_0X * dT0XIn + G0Z_0Y * dT0YIn) *
                     T00BG;
          TXZPert += response_measure *
                     (-GXZ_00 * dT00In + GXZ_0X * dT0XIn + GXZ_0Y * dT0YIn) *
                     T00BG;
          TYZPert += response_measure *
                     (-GYZ_00 * dT00In + GYZ_0X * dT0XIn + GYZ_0Y * dT0YIn) *
                     T00BG;
        } // Loop over the  x-coordinate of causal patch
       }  // Loop over the y-coordiante of causal patch
      }   // Loop over the z of interest

    PERTURBATION_FINISHUP:

      TOutFull->Set(xS, yS, etaS, T00BG + T00Pert, TXXBG + TXXPert, TYYBG + TYYPert, 
                    TZZBG + TZZPert, T0XPert, T0YPert, T0ZPert, TXYPert, TYZPert, TXZPert);

    } // Loop over the x-coordinate of grid
   }  // Loop over the y-coordinate of grid
  }   // Loop over the eta-coordinate of grid

  // COMMANDLINE OUTPUT //
  std::cerr<< "#PERTURBATION DONE" << std::endl;
}

void Setup() {
  using namespace KoMPoSTParameters;

  if (EVOLUTION_MODE == 2) {
    std::cerr << "#ERROR: low-k mode is not implemented for 3D evolution."
              << std::endl;
    exit(1);
  }

  if (EVOLUTION_MODE == 3) {
    std::cerr << "#ExactFS does not use response functions." << std::endl;
    return;
  }

  if (MOMENTUM_PERTURBATIONS) {
    std::cerr << "#ERROR: momentum perturbations are not implemented."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // SETUP GREENS FUNCTIONS FOR ENERGY-MOMENTUM PERTURBATIONS //
  const int NumberOfPoints_r = 151;
  const double rMin=0;
  const double rMax=1.5;
  const int NumberOfPoints_eta = 151;
  const double etaMin=0;
  const double etaMax=3.;
  GreensFunctions::Setup(NumberOfPoints_r, rMin, rMax, NumberOfPoints_eta, etaMin, etaMax, ENERGY_PERTURBATIONS,
                         MOMENTUM_PERTURBATIONS, EVOLUTION_MODE);

  
  //GreensFunctions::Output(ENERGY_PERTURBATIONS, MOMENTUM_PERTURBATIONS);
}

void Run(EnergyMomentumTensorMap *TIn, EnergyMomentumTensorMap *TOutBG,
         EnergyMomentumTensorMap *TOutFull,
         ChargeCurrentMap *JIn, ChargeCurrentMap *JOut) {

  using namespace KoMPoSTParameters;
  if (EVOLUTION_MODE == 3) {
    // ExactFS does not use a background decomposition.
    // TOutBG remains zero-filled for the common output interface.
    ComputeExactFreeStreaming(TIn, TOutFull, JIn, JOut);
    return;
  }

  double EvolutionTime = TOutFull->tau - TIn->tau;
  double SigmaBG = 0.;     // Transverse background Gaussian width in fm
  double SigmaBG_eta = 0.; // Background Gaussian width in spacetime rapidity

  // The Scaler is responsible for computing the scaling variable for a given
  // eta/s. The parameters EtaOverS and EtaOverSTemperature scale are passed
  // from the common block KoMPoSTParameters.
  ScalingVariable Scaler(EtaOverS, EtaOverSTemperatureScale);

  // Some error checking to make sure that the system can evolve the low k limit
  if (KoMPoSTParameters::Regulator == "TwoPass" && EVOLUTION_MODE == 2) {
    std::cerr << "*** KoMPoST::Run *** Regulator TwoPass is not appropriate for "
                 "the low k limit, selected by EVOLUTION_MODE == 2. Aborting!"
              << std::endl;
    exit(1);
  }

  // If the two pass regulator is used, then we run through the program twice
  // The first time is just to estimate the size of the perturbations.  In the
  // second pass we regulate the perturbations
  const bool UseTwoPass =
      KoMPoSTParameters::Regulator == "TwoPass" && EVOLUTION_MODE == 1;
  int npass = UseTwoPass ? 2 : 1;
  bool IsFirstPass = true;

  // Start the run
  for (int ipass = 0; ipass < npass; ipass++) {

    if (EVOLUTION_MODE != 2) {
      SigmaBG = EvolutionTime / sqrt(2.);
      SigmaBG_eta = 3. / sqrt(2.);
    } else {
      SigmaBG = EvolutionTime / sqrt(2.);
      SigmaBG_eta = 3. / sqrt(2.);
    }

    // Evolve the backround, compute the scaling variable, K etc
    ComputeBackground(IsFirstPass, TIn, TOutBG, Scaler, SigmaBG, SigmaBG_eta,
                      EVOLUTION_MODE);

    // Maximum source-target separations used in the response convolution.
    // For EKT and legacy FS, the transverse cutoff includes the causal radius
    // and three response-regulator widths.
    double ResponseCutoff_r = 0.;
    const double ResponseCutoff_eta = 3.; // EKT table support in |Delta eta|
    if (EVOLUTION_MODE != 2) {
      ResponseCutoff_r = EvolutionTime * (1. + 3. * Sigma);
    } else {
      ResponseCutoff_r = 4 * SigmaBG; // Low-k cutoff set by background width
    }

    // Evolve the petrubations. The scaling variable and K are
    // stored in the CellData structure of TOutFull.

    ComputePerturbations(TIn, TOutBG, TOutFull, ENERGY_PERTURBATIONS, MOMENTUM_PERTURBATIONS, EVOLUTION_MODE, ResponseCutoff_r, ResponseCutoff_eta);

    // Regulate the perturbations so that the inversion problem is well posed.
    if (UseTwoPass) {
      if (IsFirstPass) {
        std::cerr << "#Preparing for second pass ... " << std::endl;
        PrepareKoMPoSTEstimate(TOutBG, TOutFull);
        IsFirstPass = false;
      }
    } else if (KoMPoSTParameters::Regulator == "TwoPass" && EVOLUTION_MODE == 0) {
      std::cerr << "#TwoPass is inactive for free streaming; using one pass." << std::endl;
    } else if (KoMPoSTParameters::Regulator == "KoMPoSTAddition") {
      std::cerr<< "###########Regulation with KoMPoSTAddition ... " << std::endl;
      RegulateKoMPoSTAddition(TOutBG, TOutFull);
    } else if (KoMPoSTParameters::Regulator == "NoRegulator") {
      std::cerr<< "#Returning unregulated output ... " << std::endl;
    } else {
      std::cerr << "#KoMPoSTParameters::Regulator string does not match any of the "
                   "expected choices TwoPass/KoMPoSTAddition/NoRegulator! Aborting"
                << std::endl;
      exit(1);
    }
  }
}
}// Namepsace KoMPoST
