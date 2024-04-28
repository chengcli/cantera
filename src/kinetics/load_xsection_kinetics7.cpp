#include <cmath>
#include <cstdio>

#include "cantera/base/stringUtils.h"
#include "cantera/kinetics/Reaction.h"
#include "cantera/kinetics/Photolysis.h"

namespace Cantera
{

pair<vector<double>, vector<double>> 
load_xsection_kinetics7(vector<string> const& files, vector<Composition> const& branches)
{
  if (files.size() != 1) {
    throw CanteraError("load_xsection_kinetics7",
                       "Only one file can be loaded for Kinetics7 format.");
  }

  auto const& filename = files[0];

  FILE* file = fopen(filename.c_str(), "r");

  if (!file) {
    throw CanteraError("load_xsection_kinetics7",
                       "Could not open file '{}'", filename);
  }

  vector<double> wavelength;
  vector<double> xsection;

  // first cross section data is always the photoabsorption cross section (no dissociation)
  int nbranch = branches.size();
  int min_is = 9999, max_ie = 0;

  char *line = NULL;
  size_t len = 0;
  ssize_t read;

  // Read each line from the file
  while ((read = getline(&line, &len, file)) != -1) {
    // Skip empty lines or lines containing only whitespace
    if (line[0] == '\n') continue;

    char equation[61];
    int is, ie, nwave;
    float temp;

    // read header
    int num = sscanf(line, "%60c%4d%4d%4d%6f", equation, &is, &ie, &nwave, &temp);
    min_is = std::min(min_is, is);
    max_ie = std::max(max_ie, ie);

    if (num != 5) {
      throw CanteraError("PhotolysisBase::loadCrossSectionKinetics7",
                         "Header format from file '{}' is wrong.", filename);
    }

    // initialize wavelength and xsection for the first time
    if (wavelength.size() == 0) {
      wavelength.resize(nwave);
      xsection.resize(nwave * nbranch);
    }

    // read content
    int ncols = 7;
    int nrows = ceil(1. * nwave / ncols);
    
    equation[60] = '\0';
    auto product = parseCompString(equation);

    auto it = std::find(branches.begin(), branches.end(), product);

    if (it == branches.end()) {
      // skip this section
      for (int i = 0; i < nrows; i++)
        getline(&line, &len, file);
    } else {
      for (int i = 0; i < nrows; i++) {
        getline(&line, &len, file);

        for (int j = 0; j < ncols; j++) {
          float wave, cross;
          int num = sscanf(line + 17*j, "%7f%10f", &wave, &cross);
          if (num != 2) {
            throw CanteraError("PhotolysisBase::loadCrossSectionKinetics7",
                               "Cross-section format from file '{}' is wrong.", filename);
          }
          int b = it - branches.begin();
          int k = i * ncols + j;

          if (k >= nwave) break;
          // Angstrom -> m
          wavelength[k] = wave * 1.e-10;
          // cm^2 -> m^2
          xsection[k * nbranch + b] = cross * 1.e-4;
        }
      }
    }
  }

  // remove unused wavelength and xsection
  wavelength = vector<double>(wavelength.begin() + min_is - 1,
                              wavelength.begin() + max_ie);

  xsection = vector<double>(xsection.begin() + (min_is - 1) * nbranch,
                            xsection.begin() + max_ie * nbranch);

  // A -> A is the total cross section in kinetics7 format
  // need to subtract the other branches

  for (size_t i = 0; i < wavelength.size(); i++) {
    for (int j = 1; j < nbranch; j++) {
      xsection[i * nbranch] -= xsection[i * nbranch + j];
    }
    xsection[i * nbranch] = std::max(xsection[i * nbranch], 0.);
  }

  free(line);
  fclose(file);

  /* debug
  for (size_t i = 0; i < wavelength.size(); i++) {
    printf("%g ", wavelength[i]);
    for (int j = 0; j < nbranch; j++) {
      printf("%g ",xsection[i * nbranch + j]);
    }
    printf("\n");
  }*/

  return {std::move(wavelength), std::move(xsection)};
}

} // namespace Cantera
