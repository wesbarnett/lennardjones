
#include "Rdf.h"

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
