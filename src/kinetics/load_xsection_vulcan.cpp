#include <cmath>
#include <stdio.h>

#include "cantera/kinetics/Reaction.h"
#include "cantera/kinetics/Photolysis.h"

namespace Cantera
{

pair<vector<double>, vector<double>> 
load_xsection_vulcan(vector<string> const& files, vector<Composition> const& branches)
{
  if (files.size() != 2) {
    throw CanteraError("load_xsection_Vulcan",
                       "Only two files can be loaded for Vulcan format.");
  }

  vector<double> wavelength;
  vector<double> xsection;

  return {std::move(wavelength), std::move(xsection)};
}

} // namespace Cantera
