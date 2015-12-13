
#ifndef NEIGHBORLIST_H
#define NEIGHBORLIST_H

#include "coordinates.h"
#include "triclinicbox.h"
#include "utils.h"

class NeighborList {
    private:
        vector <vector <int> > list;
        double rlist2;
    public:
        NeighborList();
        int GetNeighbor(int i, int j);
        int GetSize(int i);
        void Init(int natoms, double rlist);
        void Update(vector <coordinates> &x, triclinicbox &box);
};

#endif
