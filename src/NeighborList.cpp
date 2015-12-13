
#include "NeighborList.h"

NeighborList::NeighborList()
{
}

void NeighborList::Init(int natoms, double rlist)
{
    list.resize(natoms);
    this->rlist2 = rlist*rlist;
    return;
}

void NeighborList::Update(vector <coordinates> &x, triclinicbox &box)
{
    
    for (unsigned int i = 0; i < this->list.size(); i++)
    {
        this->list.at(i).resize(0);
    }

    // Atoms are not double counted in the neighbor list. That is, when atom j
    // is on atom i's list, the opposite is not true.
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
