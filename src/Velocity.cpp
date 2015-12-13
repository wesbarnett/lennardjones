
#include "Velocity.h"

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
