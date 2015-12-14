
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

#include "Thermostat.h"

Thermostat::Thermostat() {} 

Thermostat::Thermostat(double reft, double coll_freq, double dt)
{
    this->reft = reft;
    this->sigma = sqrt(reft);
    this->coll_freq = coll_freq;
    this->coll_freq_dt = coll_freq*dt;
    return;
}

void Thermostat::DoCollisions(vector <coordinates> &v)
{
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dis(0.0, 1.0);
    normal_distribution<double> ndis(0.0, this->sigma);

    #pragma omp for schedule(guided, CHUNKSIZE)
    for (unsigned int i = 0; i < v.size(); i++)
    {

        if (dis(gen) < this->coll_freq_dt)
        {
            v.at(i).at(X) = ndis(gen);
            v.at(i).at(Y) = ndis(gen);
            v.at(i).at(Z) = ndis(gen);
        }

    }

    return;
}
