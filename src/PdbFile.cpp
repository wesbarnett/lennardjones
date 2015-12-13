
#include "PdbFile.h"

void PdbFile()
{

}

PdbFile::PdbFile(string filename)
{
    open(filename.c_str());
}

void PdbFile::open(string filename)
{
    oFS.open(filename.c_str());
    return;
}

void PdbFile::close()
{
    oFS.close();
}

void PdbFile::write_header(string compnd, string author, string remark)
{
    oFS << "COMPND    " << compnd << endl;
    oFS << "AUTHOR    " << author << endl;
    oFS << "REMARK    " << remark << endl;
    return;
}

void PdbFile::write_line(int atom_no, string atom_name, string res, int res_no, vector <double> &xyz, double occupancy, double beta)
{
    oFS << fixed;
    oFS << "ATOM";
    oFS << setw(7) << right << atom_no;
    oFS << setw(5) << right << atom_name;
    oFS << setw(4) << right << res;
    oFS << setw(6) << right << res_no;
    oFS << setprecision(3);
    oFS << setw(12) << right << xyz.at(X);
    oFS << setw(8) << right << xyz.at(Y);
    oFS << setw(8) << right << xyz.at(Z);
    oFS << setw(6) << right << occupancy;
    oFS << setprecision(2);
    oFS << setw(6) << beta;
    oFS << "            ";
    oFS << endl;
    return;
}
