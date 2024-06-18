#include "cantera/thermo/IdealMistPhase.h"
#include "cantera/base/utilities.h"

namespace Cantera
{

double IdealMistPhase::entropy_mole() const
{
  return GasConstant * mean_X(entropy_R_ref());
}

void IdealMistPhase::getPartialMolarEntropies(double* sbar) const
{
  const vector<double>& _s = entropy_R_ref();
  scale(_s.begin(), _s.end(), sbar, GasConstant);
}

}
