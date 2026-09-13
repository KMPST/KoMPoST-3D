#ifndef CHARGECURRENT_H
#define CHARGECURRENT_H

#include "EventInput.h"

namespace ChargeIndex {
const int Up = 0;
const int Down = 1;
const int Strange = 2;
const int NumberOfCharges = 3;
}

class ChargeCurrentMap {
public:
  double tau;
  double *J;

  int Index3D(int x, int y, int eta) {
    using EventInput::Ns;
    return x + Ns*y + Ns*Ns*eta;
  }

  int Index(int charge, int mu, int x, int y, int eta) {
    return mu + 4*(charge + ChargeIndex::NumberOfCharges*Index3D(x,y,eta));
  }

  double Get(int charge, int mu, int x, int y, int eta) {
    return J[Index(charge,mu,x,y,eta)];
  }

  void SetComponent(int charge, int mu, int x, int y, int eta, double value) {
    J[Index(charge,mu,x,y,eta)] = value;
  }

  void Set(int charge, int x, int y, int eta,
           double JTau, double JX, double JY, double JEta) {
    J[Index(charge,0,x,y,eta)] = JTau;
    J[Index(charge,1,x,y,eta)] = JX;
    J[Index(charge,2,x,y,eta)] = JY;
    J[Index(charge,3,x,y,eta)] = JEta;
  }

  double GetDensity(int charge, int x, int y, int eta,
                    const double flow[4]) {
    const double tauGeV = tau/M_HBARC;
    return flow[0]*Get(charge,0,x,y,eta)
         - flow[1]*Get(charge,1,x,y,eta)
         - flow[2]*Get(charge,2,x,y,eta)
         - tauGeV*tauGeV*flow[3]*Get(charge,3,x,y,eta);
  }

  void Reset() {
    using EventInput::Ns;
    using EventInput::Neta;
    const int size = 4*ChargeIndex::NumberOfCharges*Ns*Ns*Neta;
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < size; i++)
      J[i] = 0.0;
  }

  ChargeCurrentMap(double tau) {
    using EventInput::Ns;
    using EventInput::Neta;
    this->tau = tau;
    J = new double[4*ChargeIndex::NumberOfCharges*Ns*Ns*Neta];
    Reset();
  }

  ~ChargeCurrentMap() {
    delete[] J;
  }
};

#endif
