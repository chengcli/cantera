#include <cmath>
#include <stdio.h>

#include "cantera/kinetics/Reaction.h"
#include "cantera/kinetics/Photolysis.h"

namespace Cantera
{

void load_xsection_kinetics7(vector<string> const& files, 
                             vector<Composition> const& branches,
                             vector<double>& wavelength,
                             vector<double>& xsection)
{
  if (files.size() != 1) {
    throw CanteraError("load_xsection_kinetics7",
                       "Only one file can be loaded for Kinetics7 format.");
  }

  auto const& file = files[0];

  FILE* f = fopen(file.c_str(), "r");

  if (!f) {
    throw CanteraError("load_xsection_kinetics7",
                       "Could not open file '{}'", file);
  }

  wavelength.clear();
  xsection.clear();

  int nbranch = branches.size();
  int min_is = 9999, max_ie = 0;

  char *line = NULL;
  char *dump = NULL;
  size_t len = 0;
  ssize_t read;

  // Read each line from the file
  while ((read = getline(&line, &len, f)) != -1) {
    // Skip empty lines or lines containing only whitespace
    if (line[0] == '\n' || (line[0] && isspace((unsigned char)line[0])))
      continue;

    char equation[61];
    int is, ie, nwave;
    float temp;

    // read header
    int num = sscanf(line, "%60c%4d%4d%4d%6f", equation, &is, &ie, &nwave, &temp);
    min_is = std::min(min_is, is);
    max_ie = std::max(max_ie, ie);

    // dump comment
    getline(&dump, &len, f);

    if (num != 5) {
      throw CanteraError("PhotolysisBase::loadCrossSectionKinetics7",
                         "Header format from file '{}' is wrong.", file);
    }

    // initialize wavelength and xsection for the first time
    if (wavelength.size() == 0) {
      wavelength.resize(nwave);
      xsection.resize(nwave * nbranch);
    }

    // read content
    int ncols = 7;
    int nrows = ceil(1. * nwave / ncols);

    Reaction reaction;
    parseReactionEquation(reaction, equation, {}, nullptr);

    auto it = std::find(branches.begin(), branches.end(), reaction.reactants);

    if (it == branches.end()) {
      // skip this section
      for (int i = 0; i < nrows; i++)
        getline(&dump, &len, f);
    } else {
      for (int i = 0; i < nrows; i++) {
        getline(&line, &len, f);

        for (int j = 0; j < ncols; j++) {
          float wave, cross;
          int num = sscanf(line, "%7f%10f", &wave, &cross);
          if (num != 2) {
            throw CanteraError("PhotolysisBase::loadCrossSectionKinetics7",
                               "Content format from file '{}' is wrong.", file);
          }
          int b = it - branches.begin();
          int k = i * ncols + j;
          wavelength[k] = wave / 10.; // Angstrom to nm
          xsection[k * nbranch + b] = cross * 10.;  // cm^2 / Angstrom to cm^2 / nm
        }
      }
    }
  }

  // remove unused wavelength and xsection
  wavelength = vector<double>(wavelength.begin() + min_is - 1,
                              wavelength.begin() + max_ie);

  xsection = vector<double>(xsection.begin() + (min_is - 1) * nbranch,
                            xsection.begin() + max_ie * nbranch);

  free(dump);
  free(line);
  fclose(f);
}

} // namespace Cantera
