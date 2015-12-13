
/*
 * Copyright (C) 2015 James W. Barnett <jbarnet4@tulane.edu>
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 51
 * Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * The full license is located in a text file titled "LICENSE" in the root
 * directory of the source.
 *
 */

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
