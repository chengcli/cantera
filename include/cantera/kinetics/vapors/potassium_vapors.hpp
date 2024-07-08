#ifndef SRC_SNAP_THERMODYNAMICS_VAPORS_POTASSIUM_VAPORS_HPP_
#define SRC_SNAP_THERMODYNAMICS_VAPORS_POTASSIUM_VAPORS_HPP_

inline double sat_vapor_p_KCl_Lodders(double T) {
  double logp = 7.611 - 11382. / T;
  return 1.E5 * exp(logp);
}

#endif  // SRC_SNAP_THERMODYNAMICS_VAPORS_POTASSIUM_VAPORS_HPP_
