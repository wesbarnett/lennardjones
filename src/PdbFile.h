
#ifndef PDBFILE_H
#define PDBFILE_H

#include "gmxcpp/Trajectory.h"
#include <string>
using namespace std;

class PdbFile {
    public:
        PdbFile(string filename);
        PdbFile();
        void open(string filename);
        void write_header(string compnd, string author, string remark);
        void write_line(int atom_no, string atom_name, string res, int res_no, vector <double> &xyz, double occupancy, double beta);
        void close();
    private:
        ofstream oFS;
};

#endif
