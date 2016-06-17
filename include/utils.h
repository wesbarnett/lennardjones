
#ifndef UTILS_H
#define UTILS_H

#include "coordinates.h"
#include "cubicbox.h"

#include <fstream>
#include <math.h>
#include <string>
#include <vector>

using namespace std;

double distance(coordinates a, coordinates b, cubicbox box);
double distance2(coordinates a, coordinates b, cubicbox box);
double dot(coordinates a, coordinates b);
double magnitude(coordinates x);
coordinates pbc(coordinates a, cubicbox box);
double volume(cubicbox box);

#endif
