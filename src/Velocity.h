
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
        Velocity();
        Velocity(int nbins, double max, double min, string outfile);
        void sample(vector <coordinates> &v);
        void normalize(int natoms);
        void output();
};


#endif
