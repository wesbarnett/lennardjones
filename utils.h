
#ifndef UTILS_H
#define UTILS_H

#include "coordinates.h"
#include "triclinicbox.h"

#include <fstream>
#include <math.h>
#include <string>
#include <vector>

using namespace std;

double distance(coordinates a, coordinates b, triclinicbox box);
double distance2(coordinates a, coordinates b, triclinicbox box);
double dot(coordinates a, coordinates b);
double magnitude(coordinates x);
coordinates pbc(coordinates a, triclinicbox box);
double volume(triclinicbox box);

#endif
