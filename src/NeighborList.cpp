
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

#include "NeighborList.h"

NeighborList::NeighborList()
{
}

NeighborList::NeighborList(int natoms, double rlist)
{
    list.resize(natoms);
    this->rlist2 = rlist*rlist;
    return;
}

void NeighborList::Update(vector <coordinates> &x, cubicbox &box)
{
    
    #pragma omp parallel for schedule(guided, CHUNKSIZE)
    for (unsigned int i = 0; i < this->list.size(); i++)
    {
        this->list.at(i).resize(0);
    }

    // Atoms are not double counted in the neighbor list. That is, when atom j
    // is on atom i's list, the opposite is not true.
    #pragma omp parallel for schedule(guided, CHUNKSIZE)
    for (unsigned int i = 0; i < x.size()-1; i++)
    {
        for (unsigned int j = i+1; j < x.size(); j++)
        {
            if (distance2(x.at(i), x.at(j), box) < rlist2)
            {
                this->list.at(i).push_back(j);
            }

        }
    }

    return;
}


int NeighborList::GetSize(int i)
{
    return list.at(i).size();
}

int NeighborList::GetNeighbor(int i, int j)
{
    return list.at(i).at(j);
}
