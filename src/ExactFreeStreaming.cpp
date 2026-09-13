#include "KineticEvolution.h"
#include "EnergyMomentumTensor.h"
#include "ChargeCurrent.h"
#include "EventInput.h"

#include <cmath>
#include <omp.h>

namespace KoMPoST {

// Interpolate between two neighboring grid values.
static double LinearInterpolate(double left, double right,
                                double fraction) {
  return left + fraction * (right - left);
}

// Sample T00 at a continuous transverse position.
static double InterpolateT00(EnergyMomentumTensorMap *TIn,
                             double x, double y, int etaS) {
  using namespace EventInput;
  if (x < 0.0 || x > Ns - 1 || y < 0.0 || y > Ns - 1)
    return 0.0;

  const int x0 = int(std::floor(x));
  const int y0 = int(std::floor(y));
  const int x1 = x0 < Ns - 1 ? x0 + 1 : x0;
  const int y1 = y0 < Ns - 1 ? y0 + 1 : y0;
  const double dx = x - x0;
  const double dy = y - y0;
  const double lower = LinearInterpolate(TIn->Get(0, 0, x0, y0, etaS),
                                         TIn->Get(0, 0, x1, y0, etaS), dx);
  const double upper = LinearInterpolate(TIn->Get(0, 0, x0, y1, etaS),
                                         TIn->Get(0, 0, x1, y1, etaS), dx);
  return LinearInterpolate(lower, upper, dy);
}

static void InterpolateCharges(ChargeCurrentMap *JIn, double x, double y,
                               int etaS, double density[3]) {
  using namespace EventInput;
  for (int charge = 0; charge < ChargeIndex::NumberOfCharges; charge++)
    density[charge] = 0.0;
  if (x < 0.0 || x > Ns - 1 || y < 0.0 || y > Ns - 1)
    return;

  const int x0 = int(std::floor(x));
  const int y0 = int(std::floor(y));
  const int x1 = x0 < Ns - 1 ? x0 + 1 : x0;
  const int y1 = y0 < Ns - 1 ? y0 + 1 : y0;
  const double dx = x - x0;
  const double dy = y - y0;
  for (int charge = 0; charge < ChargeIndex::NumberOfCharges; charge++) {
    const double lower = LinearInterpolate(
        JIn->Get(charge,0,x0,y0,etaS), JIn->Get(charge,0,x1,y0,etaS), dx);
    const double upper = LinearInterpolate(
        JIn->Get(charge,0,x0,y1,etaS), JIn->Get(charge,0,x1,y1,etaS), dx);
    density[charge] = LinearInterpolate(lower, upper, dy);
  }
}

void ComputeExactFreeStreaming(EnergyMomentumTensorMap *TIn,
                               EnergyMomentumTensorMap *TOutFull,
                               ChargeCurrentMap *JIn, ChargeCurrentMap *JOut) {
  using namespace EventInput;
  using namespace KoMPoSTParameters;

  // The input supplies a surface density defined at tau -> 0.
  const double tauIn = TIn->tau;
  const double tauOut = TOutFull->tau;
  const double distance = V_FS * tauOut;
  const double surfaceToDensity = tauIn / tauOut;
  const double angularWeight = 1.0 / NPHI;

  double *cosPhi = new double[NPHI];
  double *sinPhi = new double[NPHI];
  for (int iphi = 0; iphi < NPHI; iphi++) {
    const double phi = 2.0 * M_PI * (iphi + 0.5) / NPHI;
    cosPhi[iphi] = std::cos(phi);
    sinPhi[iphi] = std::sin(phi);
  }

  #pragma omp parallel for collapse(2) schedule(runtime)
  for (int etaS = etaSTART; etaS <= etaEND; etaS++) {
    for (int yS = ySTART; yS <= yEND; yS++) {
      for (int xS = xSTART; xS <= xEND; xS++) {
        // Components follow the EnergyMomentumTensorMap sign convention.
        double T00 = 0.0;
        double T0X = 0.0;
        double T0Y = 0.0;
        double TXX = 0.0;
        double TXY = 0.0;
        double TYY = 0.0;
        double JTau[ChargeIndex::NumberOfCharges] = {0.0};
        double JX[ChargeIndex::NumberOfCharges] = {0.0};
        double JY[ChargeIndex::NumberOfCharges] = {0.0};

        for (int iphi = 0; iphi < NPHI; iphi++) {
          const double c = cosPhi[iphi];
          const double s = sinPhi[iphi];
          const double xIn = xS - distance * c / afm;
          const double yIn = yS - distance * s / afm;
          const double energy = surfaceToDensity * angularWeight
                              * InterpolateT00(TIn, xIn, yIn, etaS);
          const double vx = V_FS * c;
          const double vy = V_FS * s;
          if (JIn != NULL) {
            double initialCharge[ChargeIndex::NumberOfCharges];
            InterpolateCharges(JIn, xIn, yIn, etaS, initialCharge);
            for (int charge = 0; charge < ChargeIndex::NumberOfCharges; charge++) {
              const double density = surfaceToDensity * angularWeight
                                   * initialCharge[charge];
              JTau[charge] += density;
              JX[charge] += vx * density;
              JY[charge] += vy * density;
            }
          }
          T00 += energy;
          T0X -= vx * energy;
          T0Y -= vy * energy;
          TXX += vx * vx * energy;
          TXY -= vx * vy * energy;
          TYY += vy * vy * energy;
        }

        // Components containing the eta direction vanish in ExactFS.
        TOutFull->Set(xS, yS, etaS, T00, TXX, TYY, 0.0,
                      T0X, T0Y, 0.0, TXY, 0.0, 0.0);
        if (JOut != NULL) {
          for (int charge = 0; charge < ChargeIndex::NumberOfCharges; charge++)
            JOut->Set(charge, xS,yS,etaS,
                      JTau[charge], JX[charge], JY[charge], 0.0);
        }
      }
    }
  }

  delete[] cosPhi;
  delete[] sinPhi;
}

}
