#include <cmath>
#include <cstdio>
#include <cstring>

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

  // read cross sections
  FILE* file1 = fopen(files[0].c_str(), "r");

  if (!file1) {
    throw CanteraError("load_xsection_vulcan",
                       "Could not open file: " + files[0]);
  }

  vector<double> wavelength;
  vector<double> xsection;
  vector<double> xdiss;

  // first branch is the photoabsorption cross section (no dissociation)
  int nbranch = branches.size();

  char *line = NULL;
  size_t len = 0;
  ssize_t read;

  // read one header line
  getline(&line, &len, file1);

  // read content
  while ((read = getline(&line, &len, file1)) != -1) {
    double wave, pabs, pdis, pion;
    int num = sscanf(line, "%lf, %lf, %lf, %lf", &wave, &pabs, &pdis, &pion);
    if (num != 4) {
      throw CanteraError("load_xsection_vulcan",
                         "Could not read line: " + string(line));
    }

    // nm -> m
    wavelength.push_back(wave * 1.e-9);

    // TODO(AB): check this and we ignore pion
    // cm^2 -> m^2
    xsection.push_back(std::max(pabs - pdis, 0.) * 1.e-4);
    // populate photodissociation cross sections for all branches
    for (int i = 1; i < nbranch; ++i) {
      xsection.push_back(pdis * 1.e-4);
    }
  }

  fclose(file1);

  // read branch ratios
  FILE* file2 = fopen(files[1].c_str(), "r");

  if (!file2) {
    throw CanteraError("load_xsection_vulcan",
                       "Could not open file: " + files[1]);
  }

  std::vector<double> bwave;
  std::vector<double> bratio;
  
  // read two header lines
  getline(&line, &len, file2);
  getline(&line, &len, file2);

  // read content
  while ((read = getline(&line, &len, file1)) != -1) {
    char *token = strtok(line, ",");
    // nm -> m
    bwave.push_back(atof(token) * 1.e-9);

    for (int i = 1; i < nbranch; ++i) {
      token = strtok(NULL, ",");
      if (!token) {
        throw CanteraError("load_xsection_vulcan",
                           "Error parsing line: " + string(line));
      }
      bratio.push_back(atof(token));
    }
  }

  fclose(file2);

  // revise branch cross sections
  len = bwave.size();
  std::vector<double> br(nbranch - 1);

  for (size_t i = 0; i < wavelength.size(); ++i) {
    interpn(br.data(), &wavelength[i], bratio.data(), bwave.data(), &len, 1, nbranch - 1);

    for (int j = 1; j < nbranch; ++j) {
      xsection[i * nbranch + j] *= br[j - 1];
    }
  }

  return {std::move(wavelength), std::move(xsection)};
}

} // namespace Cantera
