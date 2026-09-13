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
#include <climits>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <omp.h>
#include <string>
//input file parser
#include "INIReader.h"

// ENERGY-MOMENTUM TENSOR DATA HOLDER //
#include "EnergyMomentumTensor.h"
#include "ChargeCurrent.h"

// Code-for evolution
#include "KineticEvolution.h"

#include "GreensFunctions.h"

// ENERGY-MOMENTUM TENSOR IO //
#include "EnergyMomentumTensorIO_music.inc"

void WriteResolvedSetup() {
  std::ofstream output(
      (KoMPoSTInputs::OutputFileTag + "_setup.ini").c_str());
  output << std::setprecision(17);
  output << "[KoMPoSTInputs]\n"
         << "tIn = " << KoMPoSTInputs::tIn << "\n"
         << "tOut = " << KoMPoSTInputs::tOut << "\n"
         << "InputFile = " << KoMPoSTInputs::InputFile << "\n"
         << "InputFormat = " << KoMPoSTInputs::InputFormat << "\n"
         << "McDipperKFactor = " << KoMPoSTInputs::McDipperKFactor << "\n"
         << "OutputFileTag = " << KoMPoSTInputs::OutputFileTag << "\n"
         << "OutputFormat = " << KoMPoSTInputs::OutputFormat << "\n\n";
  output << "[KoMPoSTParameters]\n"
         << "NuEff = " << KoMPoSTParameters::NuEff << "\n"
         << "EtaOverS = " << KoMPoSTParameters::EtaOverS << "\n"
         << "EtaOverSTemperatureScale = "
         << KoMPoSTParameters::EtaOverSTemperatureScale << "\n"
         << "Sigma_Res_Reg = "
         << KoMPoSTParameters::Sigma_Res_Reg << "\n"
         << "EVOLUTION_MODE = " << KoMPoSTParameters::EVOLUTION_MODE << "\n"
         << "V_FS = " << KoMPoSTParameters::V_FS << "\n"
         << "NPHI = " << KoMPoSTParameters::NPHI << "\n"
         << "EVOLVE_CHARGES = " << KoMPoSTParameters::EVOLVE_CHARGES << "\n"
         << "ENERGY_PERTURBATIONS = "
         << KoMPoSTParameters::ENERGY_PERTURBATIONS << "\n"
         << "MOMENTUM_PERTURBATIONS = "
         << KoMPoSTParameters::MOMENTUM_PERTURBATIONS << "\n"
         << "Regulator = " << KoMPoSTParameters::Regulator << "\n\n";
  output << "[EventInput]\n"
         << "afm = " << EventInput::afm << "\n"
         << "deta = " << EventInput::deta << "\n"
         << "Ns = " << EventInput::Ns << "\n"
         << "Neta = " << EventInput::Neta << "\n"
         << "xSTART = " << EventInput::xSTART << "\n"
         << "xEND = " << EventInput::xEND << "\n"
         << "ySTART = " << EventInput::ySTART << "\n"
         << "yEND = " << EventInput::yEND << "\n"
         << "etaSTART = " << EventInput::etaSTART << "\n"
         << "etaEND = " << EventInput::etaEND << "\n";
}

double ParseDoubleArgument(const char *name, const char *value) {
  char *end = NULL;
  double result = std::strtod(value, &end);
  if (end == value || *end != '\0' || !std::isfinite(result)) {
    std::cerr << "Error: invalid value for " << name << ": " << value
              << std::endl;
    exit(EXIT_FAILURE);
  }
  return result;
}

int ParseIntegerArgument(const char *name, const char *value) {
  char *end = NULL;
  long result = std::strtol(value, &end, 10);
  if (end == value || *end != '\0' || result < INT_MIN || result > INT_MAX) {
    std::cerr << "Error: invalid value for " << name << ": " << value
              << std::endl;
    exit(EXIT_FAILURE);
  }
  return static_cast<int>(result);
}

int main(int argc, char **argv) {
  // Try to open input file
  if (argc < 2 || (argc - 2) % 2 != 0) {
    std::cerr << "Error: invalid command-line arguments" << std::endl;
    std::cerr << "USAGE:" << std::endl;
    std::cerr << "    " << argv[0]
              << " setup.ini [--input file] [--output tag]"
              << " [--input-format auto|dat|h5] [--k-factor value]"
              << " [--output-format dat|h5]"
              << " [--tau-in fm] [--tau-out fm]"
              << " [--afm fm] [--deta value]"
              << " [--ns points] [--neta points]"
              << " [--x-start index] [--x-end index]"
              << " [--y-start index] [--y-end index]"
              << " [--eta-start index] [--eta-end index]" << std::endl;
    exit(EXIT_FAILURE);
  }

  INIReader reader(argv[1]) ;
  if (reader.ParseError()){
    std::cerr << "Error: ** " << argv[0] << " ** failed to open " << argv[1] << std::endl;
    exit(EXIT_FAILURE);
  }

  KoMPoSTInputs::Setup(reader);
  EventInput::Setup(reader);

  for (int i = 2; i < argc; i += 2) {
    std::string option = argv[i];
    if (option == "--input") KoMPoSTInputs::InputFile = argv[i + 1];
    else if (option == "--input-format")
      KoMPoSTInputs::InputFormat = argv[i + 1];
    else if (option == "--k-factor")
      KoMPoSTInputs::McDipperKFactor =
          ParseDoubleArgument(argv[i],argv[i + 1]);
    else if (option == "--output")
      KoMPoSTInputs::OutputFileTag = argv[i + 1];
    else if (option == "--output-format")
      KoMPoSTInputs::OutputFormat = argv[i + 1];
    else if (option == "--tau-in")
      KoMPoSTInputs::tIn = ParseDoubleArgument(argv[i], argv[i + 1]);
    else if (option == "--tau-out")
      KoMPoSTInputs::tOut = ParseDoubleArgument(argv[i], argv[i + 1]);
    else if (option == "--afm")
      EventInput::afm = ParseDoubleArgument(argv[i], argv[i + 1]);
    else if (option == "--deta")
      EventInput::deta = ParseDoubleArgument(argv[i], argv[i + 1]);
    else if (option == "--ns")
      EventInput::Ns = ParseIntegerArgument(argv[i], argv[i + 1]);
    else if (option == "--neta")
      EventInput::Neta = ParseIntegerArgument(argv[i], argv[i + 1]);
    else if (option == "--x-start")
      EventInput::xSTART = ParseIntegerArgument(argv[i], argv[i + 1]);
    else if (option == "--x-end")
      EventInput::xEND = ParseIntegerArgument(argv[i], argv[i + 1]);
    else if (option == "--y-start")
      EventInput::ySTART = ParseIntegerArgument(argv[i], argv[i + 1]);
    else if (option == "--y-end")
      EventInput::yEND = ParseIntegerArgument(argv[i], argv[i + 1]);
    else if (option == "--eta-start")
      EventInput::etaSTART = ParseIntegerArgument(argv[i], argv[i + 1]);
    else if (option == "--eta-end")
      EventInput::etaEND = ParseIntegerArgument(argv[i], argv[i + 1]);
    else {
      std::cerr << "Error: unknown option: " << option << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  if (KoMPoSTInputs::OutputFormat != "dat" &&
      KoMPoSTInputs::OutputFormat != "h5") {
    std::cerr << "Error: OutputFormat must be dat or h5" << std::endl;
    exit(EXIT_FAILURE);
  }
  if (KoMPoSTInputs::InputFormat != "auto" &&
      KoMPoSTInputs::InputFormat != "dat" &&
      KoMPoSTInputs::InputFormat != "h5") {
    std::cerr << "Error: InputFormat must be auto, dat or h5" << std::endl;
    exit(EXIT_FAILURE);
  }
  if (KoMPoSTInputs::InputFormat == "auto") {
    const std::string &name=KoMPoSTInputs::InputFile;
    const bool h5=name.size()>=3 && name.substr(name.size()-3)==".h5";
    const bool hdf5=name.size()>=5 && name.substr(name.size()-5)==".hdf5";
    KoMPoSTInputs::InputFormat=(h5 || hdf5) ? "h5" : "dat";
  }

  EventInput::Validate();
  EventInput::Print();
  //  Get KoMPoST parameters from input file
  KoMPoSTParameters::Setup(reader);

  WriteResolvedSetup();

  if (KoMPoSTParameters::EVOLUTION_MODE == 3) {
    const double nphiEstimate = 4.0 * M_PI * KoMPoSTParameters::V_FS
                               * KoMPoSTInputs::tOut / EventInput::afm;
    int recommendedNPHI = 32 * int(std::ceil(nphiEstimate / 32.0));
    if (recommendedNPHI < 32) recommendedNPHI = 32;
    std::cerr << "#ExactFS NPHI CHECK -- configured = "
              << KoMPoSTParameters::NPHI << ", recommended = "
              << recommendedNPHI << " for V_FS = "
              << KoMPoSTParameters::V_FS << ", tau_out = "
              << KoMPoSTInputs::tOut
              << " fm/c, min(dx,dy) = " << EventInput::afm << " fm"
              << std::endl;
    if (KoMPoSTParameters::NPHI < recommendedNPHI)
      std::cerr << "#ExactFS NPHI WARNING -- configured NPHI is below "
                << "the recommended value" << std::endl;
  }

  // SETUP OpenMP
  int NumberOfOpenMPThreads = omp_get_max_threads();
  std::cerr << "#CALCULATING WITH " << NumberOfOpenMPThreads
            << " OPEN-MP THREADS" << std::endl;

  // ALLOCATE INITIAL ENERGY-MOMENTUM TENSOR //
  EnergyMomentumTensorMap *Tmunu_In =
      new EnergyMomentumTensorMap(KoMPoSTInputs::tIn);
  ChargeCurrentMap *Jmu_In = NULL;
  if (KoMPoSTParameters::EVOLVE_CHARGES)
    Jmu_In = new ChargeCurrentMap(KoMPoSTInputs::tIn);
  
  // LOAD INITIAL ENERGY MOMENTUM TENSOR //
  if (KoMPoSTInputs::InputFormat == "h5")
    EnergyMomentumTensorMapLoadMcDipperHDF5(
        Tmunu_In,Jmu_In,KoMPoSTInputs::InputFile);
  else
    EnergyMomentumTensorMapLoad(
        Tmunu_In,Jmu_In,KoMPoSTInputs::InputFile);

  // ALLOCATE FINAL AND BACKGROUND ENERGY-MOMENTUM TENSOR //
  EnergyMomentumTensorMap *Tmunu_OutFull =
      new EnergyMomentumTensorMap(KoMPoSTInputs::tOut);
  EnergyMomentumTensorMap *Tmunu_OutBG =
      new EnergyMomentumTensorMap(KoMPoSTInputs::tOut);
  ChargeCurrentMap *Jmu_Out = NULL;
  if (KoMPoSTParameters::EVOLVE_CHARGES)
    Jmu_Out = new ChargeCurrentMap(KoMPoSTInputs::tOut);

  // LOAD RESPONSE FUNCTIONS //
  KoMPoST::Setup() ;
  // COMPUTE EVOLVE ENERGY-MOMENTUM TENSOR //
  KoMPoST::Run(Tmunu_In, Tmunu_OutBG, Tmunu_OutFull, Jmu_In, Jmu_Out);
  // END OF KOMPOST EVOLUTION //

  // WRITE OUT EVOLVED ENERGY-MOMENTUM TENSOR.
  void write_initial_conditions_MUSIC(std::string outfile_name, bool use_sigmamunu_NavierStokes, double dx, double dy, double deta, int Nx, int Ny, int Neta,
					EnergyMomentumTensorMap *Tmunu_Out_Full, EnergyMomentumTensorMap *Tmunu_Out_BG,
					ChargeCurrentMap *Jmu_Out);

  //write initial condition for hydro
  std::string tmp = KoMPoSTInputs::OutputFileTag
                    + "_4hydro." + KoMPoSTInputs::OutputFormat;
  bool use_pimunu_NS=false;
  write_initial_conditions_MUSIC(tmp, use_pimunu_NS, EventInput::afm, EventInput::afm, EventInput::deta, EventInput::Ns, EventInput::Ns, EventInput::Neta,
				 Tmunu_OutFull, Tmunu_OutBG, Jmu_Out);
  
  //tmp = std::string(OutputFileTag) + "_pimunuNS.dat";
  //use_pimunu_NS=true;
  //write_initial_conditions_MUSIC(tmp, use_pimunu_NS, EventInput::afm, EventInput::afm, EventInput::deta, EventInput::Ns, EventInput::Ns, EventInput::Neta,
  //				 Tmunu_OutFull, Tmunu_OutBG);

  ////Backup for energy momentum tensor
  //tmp = std::string(OutputFileTag) + ".input.txt";
  //EnergyMomentumTensorMapSave(Tmunu_In, tmp);
  //tmp = std::string(OutputFileTag) + "_full.dat";
  //EnergyMomentumTensorMapSave(Tmunu_OutFull, tmp);
  //tmp = std::string(OutputFileTag) + "_background.dat";
  //EnergyMomentumTensorMapSave(Tmunu_OutBG, tmp);
}
