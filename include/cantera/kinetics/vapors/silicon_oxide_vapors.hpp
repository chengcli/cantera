#pragma once

inline double svp_sio_Visscher(double T) {
  double log10p = 8.203 - 25898.9 / T;
  return 1.E5 * pow(10., log10p);
}

inline double svp_sio_Visscher_logsvp_ddT(double T) {
  return 25898.9 / (T * T);
}
