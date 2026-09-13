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
#ifndef EVENTINPUT_CPP
#define EVENTINPUT_CPP
#include "EventInput.h"
#include <cstdlib>
#include <gsl/gsl_math.h>
#include <iostream>

#include "INIReader.h"

namespace KoMPoSTInputs {

std::string InputFile("");
std::string InputFormat("auto");
std::string OutputFileTag("");
std::string OutputFormat("dat");
double McDipperKFactor=1.0;
double tIn = 0.0;
double tOut = 0.0;

void Setup(INIReader &reader) {
  InputFile =
      reader.GetString("KoMPoSTInputs", "InputFile", InputFile);
  InputFormat =
      reader.GetString("KoMPoSTInputs", "InputFormat", InputFormat);
  McDipperKFactor =
      reader.GetReal("KoMPoSTInputs", "McDipperKFactor", McDipperKFactor);
  OutputFileTag =
      reader.GetString("KoMPoSTInputs", "OutputFileTag", OutputFileTag);
  OutputFormat =
      reader.GetString("KoMPoSTInputs", "OutputFormat", OutputFormat);
  tIn = reader.GetReal("KoMPoSTInputs", "tIn", tIn);
  tOut = reader.GetReal("KoMPoSTInputs", "tOut", tOut);
}

}

namespace EventInput {

double afm = 0.08;
double deta = 0.08;
int Ns = 512;
int Neta = 256;

int xSTART = 0;
int xEND = Ns - 1;

int ySTART = 0;
int yEND = Ns - 1;

int etaSTART = 0;
int etaEND = Neta - 1;

void Setup(INIReader &reader) {
  afm = reader.GetReal("EventInput", "afm", afm);
  deta= reader.GetReal("EventInput", "deta", deta);
  Ns = reader.GetInteger("EventInput", "Ns", Ns);
  Neta = reader.GetInteger("EventInput", "Neta", Neta);

  xSTART = reader.GetInteger("EventInput", "xSTART", xSTART);
  xEND = reader.GetInteger("EventInput", "xEND", Ns - 1);

  ySTART = reader.GetInteger("EventInput", "ySTART", ySTART);
  yEND = reader.GetInteger("EventInput", "yEND", Ns - 1);

  etaSTART = reader.GetInteger("EventInput", "etaSTART", etaSTART);
  etaEND = reader.GetInteger("EventInput", "etaEND", Neta - 1);
}

void Validate() {
  if (afm <= 0.0 || deta <= 0.0 || Ns <= 0 || Neta <= 0) {
    std::cerr << "Error: grid spacing and dimensions must be positive."
              << std::endl;
    exit(EXIT_FAILURE);
  }
  if (xSTART < 0 || xSTART > xEND || xEND >= Ns ||
      ySTART < 0 || ySTART > yEND || yEND >= Ns ||
      etaSTART < 0 || etaSTART > etaEND || etaEND >= Neta) {
    std::cerr << "Error: grid index range is outside the grid."
              << std::endl;
    exit(EXIT_FAILURE);
  }
}

void Print() {
  std::cerr << "** EventInput ** Initialized a grid layout:\n" 
            << "  afm    = " << afm << "\n"
            << "  deta   = " << deta<< "\n"
            << "  Ns     = " << Ns << "\n"
            << "  Neta   = " << Neta << "\n"
            << "  xSTART = " << xSTART << "\n"
            << "  xEND   = " << xEND << "\n"
            << "  ySTART = " << ySTART << "\n"
            << "  yEND   = " << yEND <<  std::endl
            << "  etaSTART = " << etaSTART << "\n"
            << "  etaEND   = " << etaEND <<  std::endl;
}
}

namespace KoMPoSTParameters {

double NuEff = 40.0;
double EtaOverS = 2. / (4 * M_PI);
double EtaOverSTemperatureScale = 0.1; // GeV cuttoff temperature scale

std::string Regulator("TwoPass");

double Sigma;                // Sigma_Res_Reg / (tOut - tIn), dimensionless
double Sigma_Res_Reg = 0.1; // Response regulator width in fm

int EVOLUTION_MODE = 1;
// Exact free streaming: transverse velocity and angular resolution
double V_FS = 1.0;
int NPHI = 128;
int EVOLVE_CHARGES = 0;
int ENERGY_PERTURBATIONS = 0;
int MOMENTUM_PERTURBATIONS = 0;

void Setup(INIReader &reader) {

  NuEff = reader.GetReal("KoMPoSTParameters", "NuEff", NuEff);

  // Determines the scaling variable
  EtaOverS = reader.GetReal("KoMPoSTParameters", "EtaOverS", EtaOverS);

  EtaOverSTemperatureScale = reader.GetReal(
      "KoMPoSTParameters", "EtaOverSTemperatureScale", EtaOverSTemperatureScale);

  Regulator = reader.GetString("KoMPoSTParameters", "Regulator", Regulator);
  Sigma_Res_Reg = reader.GetReal(
      "KoMPoSTParameters", "Sigma_Res_Reg", Sigma_Res_Reg);
  Sigma = Sigma_Res_Reg
        / (KoMPoSTInputs::tOut - KoMPoSTInputs::tIn);

  //Sigma = reader.GetReal("KoMPoSTParameters", "Sigma", Sigma);

  EVOLUTION_MODE =
      reader.GetInteger("KoMPoSTParameters", "EVOLUTION_MODE", EVOLUTION_MODE);

  // Exact free-streaming parameters
  V_FS =
      reader.GetReal("KoMPoSTParameters", "V_FS", V_FS);
  NPHI =
      reader.GetInteger("KoMPoSTParameters", "NPHI", NPHI);
  EVOLVE_CHARGES = reader.GetInteger(
      "KoMPoSTParameters", "EVOLVE_CHARGES", EVOLVE_CHARGES);

  if (EVOLVE_CHARGES && EVOLUTION_MODE != 3) {
    std::cerr << "Charge evolution is currently implemented only for ExactFS."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  if (V_FS < 0.0 || V_FS > 1.0 || NPHI < 4) {
    std::cerr << "ExactFS requires 0 <= V_FS <= 1 and NPHI >= 4."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  ENERGY_PERTURBATIONS = reader.GetInteger(
      "KoMPoSTParameters", "ENERGY_PERTURBATIONS", ENERGY_PERTURBATIONS);

  MOMENTUM_PERTURBATIONS = reader.GetInteger(
      "KoMPoSTParameters", "MOMENTUM_PERTURBATIONS", MOMENTUM_PERTURBATIONS);

  std::cerr << "** EventInput ** Initialized KoMPoST parameters:\n" 
            << "  NuEff                    = " << NuEff << "\n"
            << "  EtaOverS                 = " << EtaOverS << "\n"
            << "  EtaOverSTemperatureScale = " << EtaOverSTemperatureScale << "\n"
            << "  Regulator                = " << Regulator << "\n"
            << "  Sigma_Res_Reg            = " << Sigma_Res_Reg << "\n"
            << "  EVOLUTION_MODE           = " << EVOLUTION_MODE
            << " ; 0 -- legacy FS, 1 -- EKT, 2 -- low k (unsupported in 3D), 3 -- ExactFS\n"
            << "  V_FS                     = " << V_FS << "\n"
            << "  NPHI                     = " << NPHI << "\n"
            << "  EVOLVE_CHARGES           = " << EVOLVE_CHARGES << "\n"
            << "  ENERGY_PERTURBATIONS     = " << ENERGY_PERTURBATIONS << "\n"
            << "  MOMENTUM_PERTURBATIONS   = " << MOMENTUM_PERTURBATIONS << "\n"
            << "  Sigma                    = " << Sigma <<  std::endl;
}
}
#endif
