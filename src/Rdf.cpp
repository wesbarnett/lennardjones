

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

#include "Rdf.h"

Rdf::Rdf() { }

Rdf::Rdf(int nbins, triclinicbox &box, string outfile)
{
    this->nbins = nbins;
    this->g.resize(this->nbins, 0.0);
    this->binwidth = box.at(0).at(0) / (2.0 * this->nbins);
    this->n = 0;
    this->outfile = outfile;
}

void Rdf::sample(vector <coordinates> &x, triclinicbox &box)
{
    this->n++;
    #pragma omp parallel for schedule(guided, CHUNKSIZE)
    for (unsigned int i = 0; i < x.size()-1; i++)
    {
        for (unsigned int j = i+1; j < x.size(); j++)
        {
            double d = distance(x.at(i), x.at(j), box);
            if (d < box.at(0).at(0)/2.0)
            {
                int ig = d/this->binwidth;
                this->g.at(ig) += 2.0;
            }
        }
    }
    return;
}

void Rdf::normalize(int natoms, triclinicbox &box)
{
    double norm_factor = 4.0/3.0 * M_PI * natoms * (natoms-1.0) * this->n * pow(this->binwidth,3) / volume(box);

    #pragma omp parallel for schedule(guided, CHUNKSIZE)
    for (int i = 0; i < this->nbins; i++)
    {
        double r = (double) i;
        double binvol = pow(r+1.0,3) - pow(r,3);
        g.at(i) /= (binvol * norm_factor);
    }

    return;
}

void Rdf::output()
{
    ofstream oFS(this->outfile.c_str());
    oFS << setprecision(6) << fixed;
    for (int i = 0; i < this->nbins; i++)
    {
        oFS << setw(20) << i*this->binwidth;
        oFS << setw(20) << this->g.at(i);
        oFS << endl;
    }
    oFS.close();
    return;
}
