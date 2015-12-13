
#ifndef THERMOSTAT_H
#define THERMOSTAT_H

#include "coordinates.h"
#include <random>
#include <math.h>

using namespace std;

class Thermostat {
    private:
        double coll_freq;
        double coll_freq_dt;
        double reft;
        double sigma;
    public:
        Thermostat();
        Thermostat(double reft, double coll_freq, double dt);
        void DoCollisions(vector <coordinates> &v);
};

#endif
