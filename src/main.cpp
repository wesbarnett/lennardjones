
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


#include "chunksize.h"
#include "coordinates.h"
#include "NeighborList.h"
#include "PdbFile.h"
#include "Rdf.h"
#include "Thermostat.h"
#include "ThermodynamicVariable.h"
#include "triclinicbox.h"
#include "utils.h"
#include "Velocity.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

#include <xdrfile/xdrfile.h>
#include <xdrfile/xdrfile_xtc.h>

#include <iostream>
#include "omp.h"
#include <random>
#include <string>
#include <vector>

using namespace std;

const double kB = 1.3806485279; // J / K

class System {
    private:
        double dt;
        double ecut;
        double entot;
        double etail;
        double halfdt;
        double halfdt2;
        double i2natoms;
        double i3natoms;
        double ke;
        double pe;
        double press;
        double ptail;
        double rcut2;
        double rho;
        double temp;
        double vol;
        int natoms;
        int natoms2;
        int nsample;
        int nsteps;
        NeighborList nlist;
        Rdf rdf;
        triclinicbox box;
        ThermodynamicVariable KineticEnergy;
        ThermodynamicVariable PotentialEnergy;
        ThermodynamicVariable Pressure;
        ThermodynamicVariable Temperature;
        ThermodynamicVariable TotalEnergy;
        Thermostat tstat;
        vector <coordinates> f;
        vector <coordinates> v;
        vector <coordinates> x;
        Velocity vel;
        XDRFILE *xd;
    public:
        System(int natoms, int nsteps, double rho, double rcut, double rlist, double temp, double dt, double mindist, double maxtries, string pdbfile, double reft, double coll_freq, string xtcfile, int rdf_nbins, string rdf_outfile, int v_nbins, double v_max, double v_min, string v_outfile);
        void CalcForce(bool samplestep);
        void CloseXTC();
        void ErrorAnalysis(int nblocks);
        void Integrate(int a, bool tcoupl, bool samplestep);
        void NormalizeAverages();
        void NormalizeRdf();
        void NormalizeVel();
        void OutputVel();
        void OutputRdf();
        void Print(int step);
        void PrintAverages();
        void PrintHeader();
        void Sample();
        void SampleRdf();
        void SampleVel();
        void UpdateNeighborList();
        void WriteXTC(int step);
};

System::System(int natoms, int nsteps, double rho, double rcut, double rlist, double temp, double dt, double mindist, double maxtries, string pdbfile, double reft, double coll_freq, string xtcfile, int rdf_nbins, string rdf_outfile, int v_nbins, double v_max, double v_min, string v_outfile)
{

    this->x.resize(natoms);
    this->v.resize(natoms);
    this->f.resize(natoms);

    this->dt = dt;
    this->halfdt = 0.5*dt;
    this->halfdt2 = 0.5*dt*dt;

    this->natoms = natoms;
    this->natoms2 = natoms*natoms;
    this->i2natoms = 1.0/(2.0*(double)natoms);
    this->i3natoms = 1.0/(3.0*(double)natoms);

    this->rcut2 = rcut*rcut;
    this->ecut = 4.0 * (1.0/pow(rcut,12) - 1.0/pow(rcut,6));
    this->etail = 8.0/3.0 * M_PI * rho * (1.0/(3.0*pow(rcut,9)) - 1.0/pow(rcut,3));
    this->ptail = 16.0/3.0 * M_PI * rho*rho * (2.0/3.0*1.0/pow(rcut,9) - 1.0/pow(rcut,3));

    this->nsample = 0;
    this->rho = rho;
    this->nsteps = nsteps;
    this->xd = xdrfile_open(xtcfile.c_str(), "w");

    // Calculate box dimensions based on density and number of atoms.
    double box_side = pow(natoms/rho,1.0/3.0);
    this->box.at(X).at(X) = box_side;
    this->box.at(X).at(Y) = 0.0;
    this->box.at(X).at(Z) = 0.0;
    this->box.at(Y).at(Y) = box_side;
    this->box.at(Y).at(X) = 0.0;
    this->box.at(Y).at(Z) = 0.0;
    this->box.at(Z).at(Z) = box_side;
    this->box.at(Z).at(X) = 0.0;
    this->box.at(Z).at(Y) = 0.0;
    cout << "Box is " << box_side << " in each dimension." << endl << endl;

    this->vol = volume(box);
    this->nlist = NeighborList(natoms, rlist);
    this->rdf = Rdf(rdf_nbins, box, rdf_outfile);
    this->tstat = Thermostat(reft, coll_freq, dt);
    this->vel = Velocity(v_nbins, v_max, v_min, v_outfile);

    // Draw from a uniform distribution centered at the origin
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> disx(-box.at(X).at(X)/2.0,box.at(X).at(X)/2.0);
    uniform_real_distribution<double> disy(-box.at(Y).at(Y)/2.0,box.at(Y).at(Y)/2.0);
    uniform_real_distribution<double> disz(-box.at(Z).at(Z)/2.0,box.at(Z).at(Z)/2.0);
    normal_distribution<double> dis_vel(0.0, sqrt(temp));
    coordinates sumv(0.0, 0.0, 0.0);
    double sumv2 = 0.0;
    const double mindist2 = mindist*mindist;

    // Generate random locations and velocities drawn from Gaussian distribution
    cout << "Generating point initial configuration...";
    
    int i = 0;
    while (i < natoms)
    {

retrypoint:

        this->x.at(i).at(X) = disx(gen);
        this->x.at(i).at(Y) = disy(gen);
        this->x.at(i).at(Z) = disz(gen);

        for (int j = 0; j < i; j++)
        {

            // Too close to to other points?
            if (distance2(x.at(i), x.at(j), box) < mindist2)
            {
                if (i > maxtries)
                {
                    cout << "ERROR: Exceeded maximum number of tries in generating initial configuration." << endl;
                }
                goto retrypoint;
            }

        }

        // Point accepted if we're here
        
        this->v.at(i).at(X) = dis_vel(gen);
        this->v.at(i).at(Y) = dis_vel(gen);
        this->v.at(i).at(Z) = dis_vel(gen);

        sumv += this->v.at(i);
        sumv2 += dot(this->v.at(i), this->v.at(i));

        i++;

    }

    sumv /= this->natoms;
    sumv2 /= this->natoms;
    double fs = sqrt(3.0*temp/sumv2);

    sumv2 = 0.0;
    for (int i = 0; i < this->natoms; i++)
    {
        this->v.at(i) = (this->v.at(i) - sumv) * fs;
        sumv2 += dot(this->v.at(i), this->v.at(i));
    }
    this->temp = sumv2 / (3.0 * this->natoms);
    this->ke = 0.5 * sumv2 / this->natoms;

    PdbFile pdb(pdbfile.c_str());
    pdb.write_header(pdbfile, "LJ MD Simulator", "First frame");
    for (int i = 0; i < natoms; i++)
    {
        pdb.write_line(i+1, "Ar", "LIG", 1, x.at(i), 1.00, 0.00);
    }
    pdb.close();

    cout << "done." << endl << endl;
}

void System::CalcForce(bool samplestep)
{

    int ncut = 0;
    double pe = 0.0;
    for (int i = 0; i < this->natoms; i++)
    {
        this->f.at(i) = 0.0;
    }

    #pragma omp parallel
    {

        vector <coordinates> f_thread(natoms);
        for (int i = 0; i < this->natoms; i++)
        {
            f_thread.at(i) = 0.0;
        }

        // Uses neighbor lists to calculate forces and energies. We didn't
        // double count the atoms on the neighbor list, so we have to look at
        // each atom's list. The last atom never has it's own list since it will
        // always be on at least one other atom's list (or it is too far away to
        // interact with any other atom)
        #pragma omp for schedule(guided, CHUNKSIZE) reduction(+:ncut,pe)
        for (int i = 0; i < this->natoms-1; i++)
        {

            for (int neighb = 0; neighb < nlist.GetSize(i); neighb++)
            {

                int j = nlist.GetNeighbor(i, neighb);
                coordinates dr = pbc(x.at(i) - x.at(j), box);
                double r2 = dot(dr,dr);

                if (r2 <= this->rcut2)
                {

                    double r2i = 1.0/r2;
                    double r6i = pow(r2i,3);
                    coordinates fr = 48.0 * r2i * r6i * (r6i - 0.5) * dr;
                    
                    // We have to count the force both on atom i from j and on j
                    // from i, since we didn't double count on the neighbor
                    // lists
                    f_thread.at(i) += fr;
                    f_thread.at(j) -= fr;

                    if (samplestep)
                    {
                        pe += 4.0*r6i*(r6i-1.0) - this->ecut;
                        ncut++;
                    }

                }

            }

        }

        #pragma omp critical
        {

            for (int i = 0; i < natoms; i++)
            {
                this->f.at(i) += f_thread.at(i);
            }

        }

    }

    if (samplestep)
    {
        ncut /= (natoms-1);
        this->pe = pe/this->natoms2;
        this->pe += etail + 0.5*ecut*(double)ncut; 

        double vir = 0.0;
        #pragma omp parallel for schedule(guided, CHUNKSIZE) reduction(+:vir)
        for (int i = 0; i < this->natoms; i++)
        {
            vir += dot(f.at(i), x.at(i));
        }
        vir /= 3.0;

        this->press = this->rho*kB*this->temp + vir/this->vol + this->ptail;
    }

    return;

}

// Velocity Verlet integrator in two parts
void System::Integrate(int a, bool tcoupl, bool samplestep)
{

    if (a == 0) 
    {
        #pragma omp parallel for schedule(guided, CHUNKSIZE)
        for (int i = 0; i < this->natoms; i++)
        {
            this->x.at(i) += this->v.at(i)*this->dt + this->f.at(i)*this->halfdt2;
            this->v.at(i) += this->f.at(i)*this->halfdt;
        }
    }
    else if (a == 1)
    {

        double sumv2 = 0.0;

        #pragma omp parallel for schedule(guided, CHUNKSIZE) reduction(+:sumv2)
        for (int i = 0; i < natoms; i++)
        {
            this->v.at(i) += this->f.at(i)*this->halfdt;
            if (samplestep)
            {
                sumv2 += dot(this->v.at(i), this->v.at(i));
            }
        }

        if (tcoupl == true)
        {
            tstat.DoCollisions(v);
        }

        if (samplestep)
        {
            this->temp = sumv2 * this->i3natoms;
            this->ke = sumv2 * this->i2natoms;
        }

    }

    return;
}

void System::Print(int step)
{
    cout << setw(14) << step;
    cout << setw(14) << step*this->dt;
    cout << setw(14) << this->temp;
    cout << setw(14) << this->press;
    cout << setw(14) << this->ke;
    cout << setw(14) << this->pe;
    cout << setw(14) << this->pe+this->ke << endl;
    return;
}

void System::PrintHeader()
{
    cout << setw(14) << "Step";
    cout << setw(14) << "Time";
    cout << setw(14) << "Temp";
    cout << setw(14) << "Press";
    cout << setw(14) << "KE";
    cout << setw(14) << "PE";
    cout << setw(14) << "Tot. En." << endl;
    return;
}

void System::PrintAverages()
{
    cout << "AVERAGES & CONSTANTS (" << this->nsample << " steps sampled out of " << this->nsteps << " total steps)" << endl;
    cout << setw(20) << "Number: " << setw(14) << this->natoms << endl;
    cout << setw(20) << "Density: " << setw(14) << this->rho << endl;
    cout << setw(20) << "Volume: " << setw(14) << volume(this->box) << endl;
    cout << setw(20) << "Temperature: " << setw(14) << this->Temperature.GetAvg() << " +/- " << setw(14) << this->Temperature.GetError() << endl;
    cout << setw(20) << "Pressure: " << setw(14) << this->Pressure.GetAvg() << " +/- " << setw(14) << this->Pressure.GetError() << endl;
    cout << setw(20) << "Kinetic Energy: " << setw(14) << this->KineticEnergy.GetAvg() << " +/- " << setw(14) << this->KineticEnergy.GetError() << endl;
    cout << setw(20) << "Potential Energy: " << setw(14) << this->PotentialEnergy.GetAvg() << " +/- " << setw(14) << this->PotentialEnergy.GetError() << endl;
    cout << setw(20) << "Total Energy: " << setw(14) << this->TotalEnergy.GetAvg() << " +/- " << setw(14) << this->TotalEnergy.GetError() << endl;
    return;
}

void System::UpdateNeighborList()
{
    this->nlist.Update(this->x, this->box);
    return;
}

void System::SampleRdf()
{
    this->rdf.sample(this->x, this->box);
    return;
}

void System::NormalizeRdf()
{
    this->rdf.normalize(this->natoms, this->box);
    return;
}

void System::OutputRdf()
{
    this->rdf.output();
    return;
}

void System::SampleVel()
{
    this->vel.sample(this->v);
    return;
}

void System::NormalizeVel()
{
    this->vel.normalize(natoms);
    return;
}

void System::OutputVel()
{
    this->vel.output();
    return;
}

void System::Sample()
{
    this->nsample++;
    this->Temperature.Sample(this->temp);
    this->Pressure.Sample(this->press);
    this->KineticEnergy.Sample(this->ke);
    this->PotentialEnergy.Sample(this->pe);
    this->TotalEnergy.Sample(this->ke+this->pe);
    return;
}


void System::WriteXTC(int step)
{
    rvec *x_xtc;
    matrix box_xtc;

    // Convert to "nanometer" (even though we are in reduced units)
    box_xtc[0][0] = this->box.at(0).at(0)/10.0;
    box_xtc[0][1] = this->box.at(0).at(1)/10.0;
    box_xtc[0][2] = this->box.at(0).at(2)/10.0;
    box_xtc[1][0] = this->box.at(1).at(0)/10.0;
    box_xtc[1][1] = this->box.at(1).at(1)/10.0;
    box_xtc[1][2] = this->box.at(1).at(2)/10.0;
    box_xtc[2][0] = this->box.at(2).at(0)/10.0;
    box_xtc[2][1] = this->box.at(2).at(1)/10.0;
    box_xtc[2][2] = this->box.at(2).at(2)/10.0;

    x_xtc = new rvec[this->x.size()];
    #pragma omp for
    for (unsigned int i = 0; i < x.size(); i++)
    {

        // Shift all the points to the center of the box
        this->x.at(i) = pbc(this->x.at(i), this->box);
        this->x.at(i).at(X) += this->box.at(X).at(X)/2.0;
        this->x.at(i).at(Y) += this->box.at(Y).at(Y)/2.0;
        this->x.at(i).at(Z) += this->box.at(Z).at(Z)/2.0;

        // Convert to "nanometers"
        x_xtc[i][X] = this->x.at(i).at(X)/10.0;
        x_xtc[i][Y] = this->x.at(i).at(Y)/10.0;
        x_xtc[i][Z] = this->x.at(i).at(Z)/10.0;

    }

    write_xtc(this->xd, this->x.size(), step, this->dt*step, box_xtc, x_xtc, 1000);

    return;
}

void System::CloseXTC()
{
    xdrfile_close(this->xd);
}

void System::ErrorAnalysis(int nblocks)
{
    this->KineticEnergy.ErrorAnalysis(nblocks);
    this->PotentialEnergy.ErrorAnalysis(nblocks);
    this->Temperature.ErrorAnalysis(nblocks);
    this->Pressure.ErrorAnalysis(nblocks);
    this->TotalEnergy.ErrorAnalysis(nblocks);
    return;
}

void System::NormalizeAverages()
{
    this->KineticEnergy.Normalize();
    this->PotentialEnergy.Normalize();
    this->Temperature.Normalize();
    this->Pressure.Normalize();
    this->TotalEnergy.Normalize();
    return;
}

int main(int argc, char *argv[])
{

    if (argc != 2)
    {
        cout << "ERROR: Configuration file should be first command line argument." << endl;
        return -1;
    }

    cout << "Reading from " << argv[1] << "..." << endl;

    // BEGIN: Read ini file --------------------------
    boost::property_tree::ptree pt;
    boost::property_tree::ini_parser::read_ini(argv[1], pt);

    char *endptr;

    cout << endl << "[ setup ]" << endl;
    const double mindist = strtod(pt.get<std::string>("setup.mindist","1.0").c_str(), &endptr); 
    cout << "mindist = " << mindist << endl;
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'setup.mindist' needs to be a real number." << endl;
        return -1;
    }
    const double maxtries = strtod(pt.get<std::string>("setup.maxtries","10e6").c_str(), &endptr); 
    cout << "maxtries = " << maxtries << endl;
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'setup.maxtries' needs to be a real number." << endl;
        return -1;
    }

    cout << endl << "[ runcontrol ]" << endl;
    const double dt = strtod(pt.get<std::string>("runcontrol.dt","0.005").c_str(), &endptr); 
    cout << "dt = " << dt << endl;
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'runcontrol.dt' needs to be a real number." << endl;
        return -1;
    }
    const int nsteps = strtol(pt.get<std::string>("runcontrol.nsteps","5000000").c_str(), &endptr, 10);
    cout << "nsteps = " << nsteps << endl;
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'runcontrol.nsteps' needs to be an integer." << endl;
        return -1;
    }
    const int eql_steps = strtol(pt.get<std::string>("runcontrol.eql_steps","10000").c_str(), &endptr, 10);
    cout << "eql_steps = " << eql_steps << endl;
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'runcontrol.eql_steps' needs to be an integer." << endl;
        return -1;
    }
    const int step_sample = strtol(pt.get<std::string>("runcontrol.nsample","1000").c_str(), &endptr, 10);
    cout << "nsample = " << step_sample << endl;
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'runcontrol.nsample' needs to be an integer." << endl;
        return -1;
    }
    const int nblocks = strtol(pt.get<std::string>("runcontrol.nblocks","5").c_str(), &endptr, 10);
    cout << "nblocks = " << nblocks << endl;
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'runcontrol.nblocks' needs to be an integer." << endl;
        return -1;
    }


    cout << endl << "[ system ]" << endl;
    const int natoms = strtol(pt.get<std::string>("system.natoms","108").c_str(), &endptr, 10);
    cout << "natoms = " << natoms << endl;
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'system.natoms' needs to be an integer." << endl;
        return -1;
    }
    const double rho = strtod(pt.get<std::string>("system.rho","0.5").c_str(), &endptr); 
    cout << "rho = " << rho << endl;
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'system.rho' needs to be a real number." << endl;
        return -1;
    }
    // Note: not a constant
    double temp = strtod(pt.get<std::string>("system.inittemp","1.0").c_str(), &endptr); 
    cout << "inittemp = " << temp << endl;
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'system.inittemp' needs to be a real number." << endl;
        return -1;
    }
    const double rcut = strtod(pt.get<std::string>("system.rcut","2.5").c_str(), &endptr); 
    cout << "rcut = " << rcut << endl;
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'system.rcut' needs to be a real number." << endl;
        return -1;
    }
    const double rlist = strtod(pt.get<std::string>("runcontrol.rlist","3.5").c_str(), &endptr); 
    cout << "rlist = " << rlist << endl;
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'runcontrol.rlist' needs to be a real number." << endl;
        return -1;
    }
    const int nlist = strtol(pt.get<std::string>("runcontrol.nlist","10").c_str(), &endptr, 10); 
    cout << "nlist = " << nlist << endl;
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'runcontrol.nlist' needs to be a real number." << endl;
        return -1;
    }

    cout << endl << "[ output ]" << endl;
    const string pdbfile = pt.get<std::string>("output.pdbfile","init.pdb");
    cout << "pdbfile = " << pdbfile << endl;
    const string xtcfile = pt.get<std::string>("output.xtcfile","traj.xtc");
    cout << "xtcfile = " << xtcfile << endl;
    const int nxtc = strtol(pt.get<std::string>("output.nxtc","1000").c_str(), &endptr, 10);
    cout << "nxtc = " << nxtc << endl;
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'output.nxtc' needs to be an integer." << endl;
        return -1;
    }
    const int nlog = strtol(pt.get<std::string>("output.nlog","1000").c_str(), &endptr, 10);
    cout << "nlog = " << nlog << endl;
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'output.nlog' needs to be an integer." << endl;
        return -1;
    }

    cout << endl << "[ temperature ]" << endl;
    const string tcouplstr = pt.get<std::string>("temperature.coupl","no");
    bool tcoupl = false;
    cout << "coupl = " << tcouplstr << endl;
    if (tcouplstr == "yes")
    {
        tcoupl = true;
    }
    const double coll_freq = strtod(pt.get<std::string>("temperature.coll_freq","0.001").c_str(), &endptr); 
    cout << "coll_freq = " << coll_freq << endl;
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'temperature.coll_freq' needs to be a real number." << endl;
        return -1;
    }
    const double reft = strtod(pt.get<std::string>("temperature.reft","1.0").c_str(), &endptr); 
    cout << "reft = " << reft << endl;
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'temperature.reft' needs to be a real number." << endl;
        return -1;
    }

    cout << endl << "[ rdf ]" << endl;
    const string dordfstr = pt.get<std::string>("rdf.sample","no");
    bool dordf = false;
    cout << "sample = " << dordfstr << endl;
    if (dordfstr == "yes")
    {
        dordf = true;
    }
    const int rdf_nbins = strtol(pt.get<std::string>("rdf.nbins","100").c_str(), &endptr, 10);
    cout << "nbins = " << rdf_nbins << endl;
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'rdf.nbins' needs to be an integer." << endl;
        return -1;
    }
    const string rdf_outfile = pt.get<std::string>("rdf.outfile","rdf.dat");
    cout << "outfile = " << rdf_outfile << endl;
    const int rdf_freq = strtol(pt.get<std::string>("rdf.freq","1000").c_str(), &endptr, 10);
    cout << "freq = " << rdf_freq << endl;
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'rdf.freq' needs to be an integer." << endl;
        return -1;
    }

    cout << endl << "[ velocity ]" << endl;
    const string dovelstr = pt.get<std::string>("velocity.sample","no");
    bool dovel = false;
    cout << "sample = " << dovelstr << endl;
    if (dovelstr == "yes")
    {
        dovel = true;
    }
    const double v_max = strtod(pt.get<std::string>("velocity.max","10.0").c_str(), &endptr); 
    cout << "max = " << v_max << endl;
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'velocity.max' needs to be a real number." << endl;
        return -1;
    }
    const double v_min = strtod(pt.get<std::string>("velocity.min","-10.0").c_str(), &endptr); 
    cout << "min = " << v_min << endl;
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'velocity.min' needs to be a real number." << endl;
        return -1;
    }
    const int v_nbins = strtol(pt.get<std::string>("velocity.nbins","100").c_str(), &endptr, 10);
    cout << "nbins = " << v_nbins << endl;
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'velocity.nbins' needs to be an integer." << endl;
        return -1;
    }
    const string v_outfile = pt.get<std::string>("velocity.outfile","vel_dist.dat");
    cout << "outfile = " << v_outfile << endl;
    const int v_freq = strtol(pt.get<std::string>("velocity.freq","1000").c_str(), &endptr, 10);
    cout << "freq = " << v_freq << endl;
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'velocity.freq' needs to be an integer." << endl;
        return -1;
    }

    cout << endl;
    // END: Read ini file -------------------------
    
    #pragma omp parallel
    #pragma omp master
    cout << "Using " << omp_get_num_threads() << " OpenMP threads." << endl;

    cout << endl;

    bool samplestep = false;

    System sys(natoms, nsteps, rho, rcut, rlist, temp, dt, mindist, maxtries, pdbfile, reft, coll_freq, xtcfile, rdf_nbins, rdf_outfile, v_nbins, v_max, v_min, v_outfile);
    sys.UpdateNeighborList();
    sys.CalcForce(samplestep);
    sys.PrintHeader();
    sys.Print(0);

    for (int step = 1; step < nsteps; step++)
    {

        if ( (step % step_sample) == 0 && (step > eql_steps) )
        {
            samplestep = true;
        }
        else
        {
            samplestep = false;
        }

        // Main part of algorithm
        sys.Integrate(0, tcoupl, samplestep);
        sys.CalcForce(samplestep);
        sys.Integrate(1, tcoupl, samplestep);


        // Update the neighbor list this step?
        if (step % nlist == 0)
        {
            sys.UpdateNeighborList();
        }

        // Sample the RDF this step?
        if (( dordf == true) && (step % rdf_freq == 0) && (step > eql_steps))
        {
            sys.SampleRdf();
        }

        // Sample the velocity distribution this step?
        if (( dovel == true) && (step % v_freq == 0) && (step > eql_steps))
        {
            sys.SampleVel();
        }

        // Do other sampling this step?
        if (samplestep)
        {
            sys.Sample();
        }

        // Print to the log this step?
        if (step % nlog == 0)
        {
            sys.Print(step);
        }

        // Write to the xtc file this step?
        if (step % nxtc == 0)
        {
            sys.WriteXTC(step);
        }

    }

    sys.CloseXTC();

    if (dordf == true)
    {
        sys.NormalizeRdf();
        sys.OutputRdf();
    }

    if (dovel == true)
    {
        sys.NormalizeVel();
        sys.OutputVel();
    }

    sys.ErrorAnalysis(nblocks);
    sys.NormalizeAverages();
    sys.PrintAverages();

    return 0;
}
