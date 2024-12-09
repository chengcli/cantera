#ifndef SVP_FUNCS_INC_H
#define SVP_FUNCS_INC_H

#include "vapors/water_vapors.hpp"
#include "vapors/ammonia_vapors.hpp"
#include "vapors/ammonium_hydrosulfide_vapors.hpp"
#include "vapors/silicon_dioxide_vapors.hpp"

namespace Cantera
{

string concatenate(const Composition& comp, char sep);

std::function<double(double)> find_svp(const Composition& reactants,
                                       const string& svp_name);

std::function<double(double)> find_logsvp_ddT(const Composition& reactants,
                                              const string& svp_name);

std::function<double(double)> find_svp_function(const Composition& reactants, 
                                                const string& svp_name);
}

#endif // SVP_FUNCS_INC_H
