
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

#include "Velocity.h"

Velocity::Velocity() { }

Velocity::Velocity(int nbins, double max, double min, string outfile)
{
    this->max = max;
    this->min = min;
    this->shift = -min;
    this->nbins = nbins;
    this->hist.resize(this->nbins);
    for (int i = 0; i < this->nbins; i++)
    {
        this->hist.at(i) = 0.0;
    }
    this->binwidth = (max - min) / (nbins);
    this->n = 0;
    this->outfile = outfile;
}

void Velocity::sample(vector <coordinates> &v)
{
    this->n++;
    double iv;
    for (unsigned int i = 0; i < v.size(); i++)
    {
        iv = (v.at(i).at(X) + this->shift) / this->binwidth;
        this->hist.at(iv).at(X) += 1.0;
        iv = (v.at(i).at(Y) + this->shift) / this->binwidth;
        this->hist.at(iv).at(Y) += 1.0;
        iv = (v.at(i).at(Z) + this->shift) / this->binwidth;
        this->hist.at(iv).at(Z) += 1.0;
    }

    return;
}

void Velocity::normalize(int natoms)
{
    double norm_factor = (double) this->n * natoms;
    for (int i = 0; i < this->nbins; i++)
    {
        this->hist.at(i).at(X) /= norm_factor;
        this->hist.at(i).at(Y) /= norm_factor;
        this->hist.at(i).at(Z) /= norm_factor;
    }
    return;
}

void Velocity::output()
{
    ofstream oFS(this->outfile.c_str());
    oFS << setprecision(6) << fixed;
    for (unsigned int i = 0; i < this->hist.size(); i++)
    {
        oFS << setw(20) << i*this->binwidth - this->shift;
        oFS << setw(20) << this->hist.at(i).at(X);
        oFS << setw(20) << this->hist.at(i).at(Y);
        oFS << setw(20) << this->hist.at(i).at(Z);
        oFS << endl;
    }
    oFS.close();
    return;
}
