
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

#include "ThermodynamicVariable.h"

ThermodynamicVariable::ThermodynamicVariable() {

    this->avg = 0.0;
}

void ThermodynamicVariable::Sample(double value)
{
    this->all.push_back(value);
    this->avg += value;

    return;
}

void ThermodynamicVariable::Normalize()
{
    this->avg /= all.size();
    return;
}

void ThermodynamicVariable::ErrorAnalysis(int nblocks)
{
    double blockavg = 0.0;
    int nsample = all.size();
    vector <double> block(nblocks);

    #pragma omp for schedule(guided, 15)
    for (int i = 0; i < nblocks; i++)
    {
        int first = i * nsample / nblocks;
        int last;

        if (i == nblocks)
        {
            last = nsample;
        }
        else
        {
            last = (i + 1) * nsample / nblocks;
        }

        for (int j = first; j < last; j++)
        {
            block.at(i) += this->all.at(j);
        }

        block.at(i) /= (double) (last - first);
        blockavg += block.at(i);
    }

    blockavg /= nblocks;

    this->error = 0.0;

    #pragma omp for schedule(guided, 15)
    for (int i = 0; i < nblocks; i++)
    {
        this->error += pow(block.at(i),2) - pow(blockavg,2);
    }

    this->error /= (nblocks-1);
    this->error = sqrt(this->error);
    return;
}

double ThermodynamicVariable::GetAvg()
{
    return avg;
}

double ThermodynamicVariable::GetError()
{
    return error;
}
