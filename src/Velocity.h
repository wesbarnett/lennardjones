
#ifndef VELOCITY_H
#define VELOCITY_H

#include "coordinates.h"

#include <fstream>
#include <iomanip>
#include <math.h>
#include <string>
#include <vector>

using namespace std;

class Velocity {
    private:
        double binwidth;
        double max;
        double min;
        double n;
        double shift;
        int freq;
        int nbins;
        string outfile;
        vector <coordinates> hist;
    public:
        Velocity(int nbins, double max, double min, string outfile);
        void sample(vector <coordinates> &v);
        void normalize(int natoms);
        void output();
};


#endif
