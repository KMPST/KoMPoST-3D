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
#ifndef EVENTINPUT_H
#define EVENTINPUT_H
#ifndef M_HBARC
#define M_HBARC 0.197326979
#endif
#include <string>

class INIReader;

namespace KoMPoSTInputs {
extern std::string InputFile;
extern std::string InputFormat;
extern std::string OutputFileTag;
extern std::string OutputFormat;
extern double McDipperKFactor;
extern double tIn;
extern double tOut;

void Setup(INIReader &reader);
}

namespace EventInput {

extern double afm;
extern double deta;
extern int Ns;
extern int Neta;

extern int xSTART;
extern int xEND;

extern int ySTART;
extern int yEND;

extern int etaSTART;
extern int etaEND;

void Setup(INIReader &reader);
void Validate();
void Print();
}

namespace KoMPoSTParameters {
extern double NuEff;
extern double EtaOverS;
extern double EtaOverSTemperatureScale;
extern double Sigma;
extern double Sigma_Res_Reg;

extern std::string Regulator;

extern int EVOLUTION_MODE;
// Exact free-streaming parameters
extern double V_FS;
extern int NPHI;
extern int EVOLVE_CHARGES;
extern int ENERGY_PERTURBATIONS;
extern int MOMENTUM_PERTURBATIONS;

void Setup(INIReader &reader);
}
#endif
