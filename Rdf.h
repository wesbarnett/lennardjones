
#ifndef RDF_H
#define RDF_H

#include "coordinates.h"
#include "triclinicbox.h"
#include "utils.h"

#include <fstream>
#include <iomanip>
#include <math.h>
#include <string>
#include <vector>

using namespace std;

class Rdf {
    private:
        double binwidth;
        double n;
        int nbins;
        string outfile;
        vector <double> g;
        int freq;
    public:
        Rdf(int nbins, triclinicbox &box, string outfile);
        void sample(vector <coordinates> &x, triclinicbox &box);
        void normalize(int natoms, triclinicbox &box);
        void output();
};


#endif
