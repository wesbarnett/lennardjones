
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
#include "cubicbox.h"
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

const double kB = 1.3806485279; // Boltzmann's Constant (J / K)
const double oneSixth = 1.0/6.0;

class System 
{
    private:
        double dt;              // time step
        double ecut;            // potential energy at cutoff
        double entot;           // instantaneous total energy (ke + pe)
        double etail;           // energy tail correction
        double halfdt;          // 0.5 * dt
        double halfdt2;         // 0.5 * dt*dt
        double halfecut;        // 0.5 * ecut
        double inatomsm1;       // 1.0/natoms - 1.0
        double i2natoms;        // 1.0(2.0*natoms)
        double i3natoms;        // 1.0(3.0*natoms)
        double ke;              // instantaneous kinetic energy
        double pe;              // instantaneous potential energy
        double press;           // instantaneous pressure
        double ptail;           // pressure tail correction
        double rcut2;           // rcut*rcut
        double rho;             // density (constant)
        double rhokB;           // density * boltzmann's constant
        double temp;            // instantaneous temperature
        double vol;             // volume (constant)
        int natoms;             // number of atoms in system
        int nsample;            // counter of number of samples
        int nsteps;             // number of steps for simulation to perform
        NeighborList nlist;
        Rdf rdf;
        cubicbox box;
        ThermodynamicVariable KineticEnergy;
        ThermodynamicVariable PotentialEnergy;
        ThermodynamicVariable Pressure;
        ThermodynamicVariable Temperature;
        ThermodynamicVariable TotalEnergy;
        Thermostat tstat;
        vector <coordinates> f; // forces
        vector <coordinates> v; // velocities
        vector <coordinates> x; // positions
        Velocity vel;
        XDRFILE *xd;
    public:
        System(int natoms, int nsteps, double rho, double rcut, double rlist, double temp, double dt, double mindist, double maxtries, string pdbfile, double reft, double coll_freq, string xtcfile, int rdf_nbins, string rdf_outfile, int v_nbins, double v_max, double v_min, string v_outfile);
        void CalcForce();
        void CloseXTC();
        void ErrorAnalysis(int nblocks);
        void Integrate(int a, bool tcoupl);
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
    this->inatomsm1 = 1.0/(natoms - 1);
    this->i2natoms = 1.0/(2.0*(double)natoms);
    this->i3natoms = 1.0/(3.0*(double)natoms);

    this->rcut2 = rcut*rcut;
    this->ecut = 4.0 * (1.0/pow(rcut,12) - 1.0/pow(rcut,6));
    this->halfecut = ecut/2.0;
    this->etail = 8.0/3.0 * M_PI * rho * (1.0/(3.0*pow(rcut,9)) - 1.0/pow(rcut,3));
    this->ptail = 16.0/3.0 * M_PI * rho*rho * (2.0/3.0*1.0/pow(rcut,9) - 1.0/pow(rcut,3));

    this->nsample = 0;
    this->rho = rho;
    this->nsteps = nsteps;
    this->xd = xdrfile_open(xtcfile.c_str(), "w");
    this->rhokB = rho*kB;

    // Calculate box dimensions based on density and number of atoms.
    double box_side = pow(natoms/rho,1.0/3.0);
    this->box[X] = box_side;
    this->box[Y] = box_side;
    this->box[Z] = box_side;
    cout << "Box is " << box_side << " in each dimension." << endl << endl;

    this->vol = volume(box);
    this->nlist = NeighborList(natoms, rlist);
    this->rdf = Rdf(rdf_nbins, box, rdf_outfile);
    this->tstat = Thermostat(reft, coll_freq, dt);
    this->vel = Velocity(v_nbins, v_max, v_min, v_outfile);

    // Draw from a uniform distribution centered at the origin
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> disx(-box[X]/2.0,box[X]/2.0);
    uniform_real_distribution<double> disy(-box[Y]/2.0,box[Y]/2.0);
    uniform_real_distribution<double> disz(-box[Z]/2.0,box[Z]/2.0);
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

        this->x[i][X] = disx(gen);
        this->x[i][Y] = disy(gen);
        this->x[i][Z] = disz(gen);

        for (int j = 0; j < i; j++)
        {

            // Too close to to other points?
            if (distance2(x[i], x[j], box) < mindist2)
            {
                if (i > maxtries)
                {
                    cout << "ERROR: Exceeded maximum number of tries in generating initial configuration." << endl;
                }
                goto retrypoint;
            }

        }

        // Point accepted if we're here
        
        this->v[i][X] = dis_vel(gen);
        this->v[i][Y] = dis_vel(gen);
        this->v[i][Z] = dis_vel(gen);

        sumv += this->v[i];
        sumv2 += dot(this->v[i], this->v[i]);

        i++;

    }

    sumv /= this->natoms;
    sumv2 /= this->natoms;
    double fs = sqrt(3.0*temp/sumv2);

    sumv2 = 0.0;
    for (int i = 0; i < this->natoms; i++)
    {
        this->v[i] = (this->v[i] - sumv) * fs;
        sumv2 += dot(this->v[i], this->v[i]);
    }
    this->temp = sumv2 / (3.0 * this->natoms);
    this->ke = 0.5 * sumv2 / this->natoms;

    PdbFile pdb(pdbfile.c_str());
    pdb.write_header(pdbfile, "LJ MD Simulator", "First frame");
    for (int i = 0; i < natoms; i++)
    {
        pdb.write_line(i+1, "Ar", "LIG", 1, x[i], 1.00, 0.00);
    }
    pdb.close();

    cout << "done." << endl << endl;
}

void System::CalcForce()
{

    int ncut = 0;
    double pe = 0.0;
    for (int i = 0; i < this->natoms; i++)
    {
        this->f[i] = 0.0;
    }

    #pragma omp parallel
    {

        vector <coordinates> f_thread(natoms);
        for (int i = 0; i < this->natoms; i++)
        {
            f_thread[i] = 0.0;
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
                coordinates dr = pbc(x[i] - x[j], box);
                double r2 = dot(dr,dr);

                if (r2 <= this->rcut2)
                {

                    double r2i = 1.0/r2;
                    double r6i = pow(r2i,3);
                    coordinates fr = 48.0 * r2i * r6i * (r6i - 0.5) * dr;
                    
                    // We have to count the force both on atom i from j and on j
                    // from i, since we didn't double count on the neighbor
                    // lists
                    f_thread[i] += fr;
                    f_thread[j] -= fr;

                    pe += 4.0*r6i*(r6i-1.0) - this->ecut;
                    ncut++;

                }

            }

        }

        #pragma omp critical
        {

            for (int i = 0; i < natoms; i++)
            {
                this->f[i] += f_thread[i];
            }

        }

    }

    double vir = 0.0;
    #pragma omp parallel for schedule(guided, CHUNKSIZE) reduction(+:vir)
    for (int i = 0; i < natoms; i++)
    {
        vir += dot(f[i], x[i]);
    }
    vir *= oneSixth;
    this->press = this->rhokB * this->temp + vir/this->vol + this->ptail;

    ncut *= inatomsm1;
    this->pe = pe/this->natoms + this->etail + halfecut*(double)ncut; 

    return;

}

// Velocity Verlet integrator in two parts
void System::Integrate(int a, bool tcoupl)
{

    if (a == 0) 
    {
        #pragma omp parallel for schedule(guided, CHUNKSIZE)
        for (int i = 0; i < this->natoms; i++)
        {
            this->x[i] += this->v[i]*this->dt + this->f[i]*this->halfdt2;
            this->v[i] += this->f[i]*this->halfdt;
        }
    }
    else if (a == 1)
    {

        double sumv2 = 0.0;

        #pragma omp parallel for schedule(guided, CHUNKSIZE) reduction(+:sumv2)
        for (int i = 0; i < natoms; i++)
        {
            this->v[i] += this->f[i]*this->halfdt;
            sumv2 += dot(this->v[i], this->v[i]);
        }

        if (tcoupl == true)
        {
            tstat.DoCollisions(v);
        }

        this->temp = sumv2 * this->i3natoms;
        this->ke = sumv2 * this->i2natoms;

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
    cout << setw(20) << "Volume: " << setw(14) << this->vol << endl;
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
    this->KineticEnergy.Sample(this->ke);
    this->PotentialEnergy.Sample(this->pe);
    this->Pressure.Sample(this->press);
    this->Temperature.Sample(this->temp);
    this->TotalEnergy.Sample(this->ke+this->pe);
    return;
}


void System::WriteXTC(int step)
{
    rvec *x_xtc;
    matrix box_xtc;

    // Convert to "nanometer" (even though we are in reduced units)
    box_xtc[X][X] = this->box[X]*0.1;
    box_xtc[X][Y] = 0.0;
    box_xtc[X][Z] = 0.0;
    box_xtc[Y][X] = 0.0;
    box_xtc[Y][Y] = this->box[Y]*0.1;
    box_xtc[Y][Z] = 0.0;
    box_xtc[Z][X] = 0.0;
    box_xtc[Z][Y] = 0.0;
    box_xtc[Z][Z] = this->box[Z]*0.1;

    x_xtc = new rvec[this->x.size()];
    #pragma omp for
    for (unsigned int i = 0; i < x.size(); i++)
    {

        // Shift all the points to the center of the box
        this->x[i] = pbc(this->x[i], this->box);
        this->x[i][X] += this->box[X]*0.5;
        this->x[i][Y] += this->box[Y]*0.5;
        this->x[i][Z] += this->box[Z]*0.5;

        // Convert to "nanometers"
        x_xtc[i][X] = this->x[i][X]*0.1;
        x_xtc[i][Y] = this->x[i][Y]*0.1;
        x_xtc[i][Z] = this->x[i][Z]*0.1;

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
    this->Pressure.ErrorAnalysis(nblocks);
    this->Temperature.ErrorAnalysis(nblocks);
    this->TotalEnergy.ErrorAnalysis(nblocks);
    return;
}

void System::NormalizeAverages()
{
    this->KineticEnergy.Normalize();
    this->PotentialEnergy.Normalize();
    this->Pressure.Normalize();
    this->Temperature.Normalize();
    this->TotalEnergy.Normalize();
    return;
}

int main(int argc, char *argv[])
{
    cout << endl;
    cout << "=============================================================================" << endl;
    cout << "Lennard Jones Molecular Dynamics Simulator" << endl;
    cout << "=============================================================================" << endl << endl;

    cout << "Copyright (C) 2015 James W. Barnett <jbarnet4@tulane.edu>" << endl << endl;
    cout << "https://github.com/wesbarnett/lennardjones" << endl << endl;

    cout << "This program is free software; you can redistribute it and/or modify it under" << endl;
    cout << "the terms of the GNU General Public License as published by the Free Software" << endl;
    cout << "Foundation; either version 2 of the License, or (at your option) any later" << endl;
    cout << "version. See 'LICENSE' in source code for full license." << endl << endl;

    cout << "-----------------------------------------------------------------------------" << endl << endl;

    const string defaultfilename = "md.ini";
    string inifile;

    if (argc != 2)
    {
        cerr << "No configuration file specified. Using default filename '" << defaultfilename << "'." << endl;
        inifile = defaultfilename;
    }
    else
    {
        inifile = argv[1];
    }

    cout << "Reading from " << inifile << "..." << endl;

    boost::property_tree::ptree pt;

    try 
    {
        boost::property_tree::ini_parser::read_ini(inifile, pt);
    }
    catch( std::exception &ex )
    {
        cerr << ex.what() << endl;
        cerr << "Using default settings for all options." << endl;
    }

    char *endptr;

    const double mindist = strtod(pt.get<std::string>("setup.mindist","1.0").c_str(), &endptr); 
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'setup.mindist' needs to be a real number." << endl;
        return -1;
    }
    const double maxtries = strtod(pt.get<std::string>("setup.maxtries","10e6").c_str(), &endptr); 
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'setup.maxtries' needs to be a real number." << endl;
        return -1;
    }
    const double dt = strtod(pt.get<std::string>("runcontrol.dt","0.005").c_str(), &endptr); 
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'runcontrol.dt' needs to be a real number." << endl;
        return -1;
    }
    const int nsteps = strtol(pt.get<std::string>("runcontrol.nsteps","5000000").c_str(), &endptr, 10);
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'runcontrol.nsteps' needs to be an integer." << endl;
        return -1;
    }
    const int eql_steps = strtol(pt.get<std::string>("runcontrol.eql_steps","10000").c_str(), &endptr, 10);
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'runcontrol.eql_steps' needs to be an integer." << endl;
        return -1;
    }
    const int step_sample = strtol(pt.get<std::string>("runcontrol.nsample","1000").c_str(), &endptr, 10);
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'runcontrol.nsample' needs to be an integer." << endl;
        return -1;
    }
    const int nblocks = strtol(pt.get<std::string>("runcontrol.nblocks","5").c_str(), &endptr, 10);
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'runcontrol.nblocks' needs to be an integer." << endl;
        return -1;
    }


    const int natoms = strtol(pt.get<std::string>("system.natoms","108").c_str(), &endptr, 10);
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'system.natoms' needs to be an integer." << endl;
        return -1;
    }
    const double rho = strtod(pt.get<std::string>("system.rho","0.5").c_str(), &endptr); 
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'system.rho' needs to be a real number." << endl;
        return -1;
    }
    // Note: not a constant
    double temp = strtod(pt.get<std::string>("system.inittemp","1.0").c_str(), &endptr); 
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'system.inittemp' needs to be a real number." << endl;
        return -1;
    }
    const double rcut = strtod(pt.get<std::string>("system.rcut","2.5").c_str(), &endptr); 
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'system.rcut' needs to be a real number." << endl;
        return -1;
    }
    const double rlist = strtod(pt.get<std::string>("runcontrol.rlist","3.5").c_str(), &endptr); 
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'runcontrol.rlist' needs to be a real number." << endl;
        return -1;
    }
    const int nlist = strtol(pt.get<std::string>("runcontrol.nlist","10").c_str(), &endptr, 10); 
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'runcontrol.nlist' needs to be a real number." << endl;
        return -1;
    }

    const string pdbfile = pt.get<std::string>("output.pdbfile","init.pdb");
    const string xtcfile = pt.get<std::string>("output.xtcfile","traj.xtc");
    const int nxtc = strtol(pt.get<std::string>("output.nxtc","1000").c_str(), &endptr, 10);
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'output.nxtc' needs to be an integer." << endl;
        return -1;
    }
    const int nlog = strtol(pt.get<std::string>("output.nlog","1000").c_str(), &endptr, 10);
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'output.nlog' needs to be an integer." << endl;
        return -1;
    }

    const string tcouplstr = pt.get<std::string>("temperature.coupl","no");
    bool tcoupl = false;
    if (tcouplstr == "yes")
    {
        tcoupl = true;
    }
    const double coll_freq = strtod(pt.get<std::string>("temperature.coll_freq","0.001").c_str(), &endptr); 
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'temperature.coll_freq' needs to be a real number." << endl;
        return -1;
    }
    const double reft = strtod(pt.get<std::string>("temperature.reft","1.0").c_str(), &endptr); 
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'temperature.reft' needs to be a real number." << endl;
        return -1;
    }

    const string dordfstr = pt.get<std::string>("rdf.sample","no");
    bool dordf = false;
    if (dordfstr == "yes")
    {
        dordf = true;
    }
    const int rdf_nbins = strtol(pt.get<std::string>("rdf.nbins","100").c_str(), &endptr, 10);
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'rdf.nbins' needs to be an integer." << endl;
        return -1;
    }
    const string rdf_outfile = pt.get<std::string>("rdf.outfile","rdf.dat");
    const int rdf_freq = strtol(pt.get<std::string>("rdf.freq","1000").c_str(), &endptr, 10);
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'rdf.freq' needs to be an integer." << endl;
        return -1;
    }

    const string dovelstr = pt.get<std::string>("velocity.sample","no");
    bool dovel = false;
    if (dovelstr == "yes")
    {
        dovel = true;
    }
    const double v_max = strtod(pt.get<std::string>("velocity.max","10.0").c_str(), &endptr); 
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'velocity.max' needs to be a real number." << endl;
        return -1;
    }
    const double v_min = strtod(pt.get<std::string>("velocity.min","-10.0").c_str(), &endptr); 
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'velocity.min' needs to be a real number." << endl;
        return -1;
    }
    const int v_nbins = strtol(pt.get<std::string>("velocity.nbins","100").c_str(), &endptr, 10);
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'velocity.nbins' needs to be an integer." << endl;
        return -1;
    }
    const string v_outfile = pt.get<std::string>("velocity.outfile","vel_dist.dat");
    const int v_freq = strtol(pt.get<std::string>("velocity.freq","1000").c_str(), &endptr, 10);
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'velocity.freq' needs to be an integer." << endl;
        return -1;
    }

    cout << endl;
    cout << setw(30) << left << "[ setup ]" << endl;
    cout << setw(30) << left << "maxtries = " << setw(30) << left << maxtries << endl;
    cout << setw(30) << left << "mindist = "  << setw(30) << left << mindist << endl;
    cout << setw(30) << left << "dt = " << dt << endl;
    cout << endl;
    cout << setw(30) << left << "[ runcontrol ]" << endl;
    cout << setw(30) << left << "nsteps = " << setw(30) << left << nsteps << endl;
    cout << setw(30) << left << "eql_steps = " << setw(30) << left << eql_steps << endl;
    cout << setw(30) << left << "nsample = " << setw(30) << left << step_sample << endl;
    cout << setw(30) << left << "nblocks = " << setw(30) << left << nblocks << endl;
    cout << endl;
    cout << setw(30) << left << "[ system ]" << endl;
    cout << setw(30) << left << "natoms = " << setw(30) << left << natoms << endl;
    cout << setw(30) << left << "rho = " << setw(30) << left << rho << endl;
    cout << setw(30) << left << "inittemp = " << setw(30) << left << temp << endl;
    cout << setw(30) << left << "rcut = " << setw(30) << left << rcut << endl;
    cout << setw(30) << left << "rlist = " << setw(30) << left << rlist << endl;
    cout << setw(30) << left << "nlist = " << setw(30) << left << nlist << endl;
    cout << endl;
    cout << setw(30) << left << "[ output ]" << endl;
    cout << setw(30) << left << "pdbfile = " << setw(30) << left << pdbfile << endl;
    cout << setw(30) << left << "xtcfile = " << setw(30) << left << xtcfile << endl;
    cout << setw(30) << left << "nxtc = " << setw(30) << left << nxtc << endl;
    cout << setw(30) << left << "nlog = " << setw(30) << left << nlog << endl;
    cout << endl;
    cout << setw(30) << left << "[ temperature ]" << endl;
    cout << setw(30) << left << "reft = " << setw(30) << left << reft << endl;
    cout << setw(30) << left << "coupl = " << setw(30) << left << tcouplstr << endl;
    cout << setw(30) << left << "coll_freq = " << setw(30) << left << coll_freq << endl;
    cout << endl;
    cout << setw(30) << left << "[ rdf ]" << endl;
    cout << setw(30) << left << "sample = " << setw(30) << left << dordfstr << endl;
    cout << setw(30) << left << "nbins = " << setw(30) << left << rdf_nbins << endl;
    cout << setw(30) << left << "outfile = " << setw(30) << left << rdf_outfile << endl;
    cout << setw(30) << left << "freq = " << setw(30) << left << rdf_freq << endl;
    cout << endl;
    cout << setw(30) << left << "[ velocity ]" << endl;
    cout << setw(30) << left << "sample = " << setw(30) << left << dovelstr << endl;
    cout << setw(30) << left << "min = " << setw(30) << left << v_min << endl;
    cout << setw(30) << left << "max = " << setw(30) << left << v_max << endl;
    cout << setw(30) << left << "nbins = " << setw(30) << left << v_nbins << endl;
    cout << setw(30) << left << "outfile = " << setw(30) << left << v_outfile << endl;
    cout << setw(30) << left << "freq = " << setw(30) << left << v_freq << endl;
    cout << endl;
    
    #pragma omp parallel
    #pragma omp master
    cout << "Using " << omp_get_num_threads() << " OpenMP threads." << endl;

    cout << endl;

    cout << setprecision(6) << fixed << right;

    System sys(natoms, nsteps, rho, rcut, rlist, temp, dt, mindist, maxtries, pdbfile, reft, coll_freq, xtcfile, rdf_nbins, rdf_outfile, v_nbins, v_max, v_min, v_outfile);
    sys.UpdateNeighborList();
    sys.CalcForce();
    sys.PrintHeader();
    sys.Print(0);

    for (int step = 1; step < nsteps; step++)
    {


        // Main part of algorithm
        sys.Integrate(0, tcoupl);
        sys.CalcForce();
        sys.Integrate(1, tcoupl);


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
        if ( (step % step_sample) == 0 && (step > eql_steps) )
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
