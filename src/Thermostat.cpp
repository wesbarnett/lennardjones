
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

    #pragma omp for schedule(guided, 15)
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
