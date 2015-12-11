
/*
 * Author: James W. Barnett
 * Date: December 7, 2015
 *
 * Simple molecular dynamics program
 * Uses velocity Verlet integrator
 *
 * Currently can do NVT and NVE LJ particles. NVT uses Andersen thermostat. Box
 * is cubic. Reduced units are used.
 *
 */

#include "omp.h"

#include "coordinates.h"
#include "Rdf.h"
#include "triclinicbox.h"
#include "PdbFile.h"
#include "utils.h"
#include "Velocity.h"

#include <xdrfile/xdrfile.h>
#include <xdrfile/xdrfile_xtc.h>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

#include <iostream>
#include <random>
#include <vector>
#include <string>

using namespace std;

int init(vector <coordinates> &x, vector <coordinates> &v, vector <coordinates> &xm, int natoms, triclinicbox box, double &temp, double dt, double mindist, double maxtries, double &ke);
void integrate(int a, vector <coordinates> &x, vector <coordinates> &v, vector <coordinates> &f, int natoms, double dt, double &temp, double &ke, bool tcoupl, double reft, double coll_freq);
void force(vector <coordinates> &f, double &en, vector <coordinates> &x, triclinicbox box, int natoms, double rcut2, double ecut, vector < vector <int> > &neighb_list);

void update_neighb_list(vector < vector <int> > &neighb_list, vector <coordinates> &x, triclinicbox &box, double rlist2);

void write_frame(XDRFILE *xd, vector <coordinates> x, triclinicbox &box, int step, double dt);

int main(int argc, char *argv[])
{

    if (argc != 2)
    {
        cout << "ERROR: Configuration file should be first command line argument." << endl;
        return -1;
    }

    cout << "Reading from " << argv[1] << "..." << endl;

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
    const int nsteps = strtol(pt.get<std::string>("runcontrol.nsteps","10000").c_str(), &endptr, 10);
    cout << "nsteps = " << nsteps << endl;
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'runcontrol.nsteps' needs to be an integer." << endl;
        return -1;
    }
    const int eql_steps = strtol(pt.get<std::string>("runcontrol.eql_steps","1000").c_str(), &endptr, 10);
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
    const int threads_n = strtol(pt.get<std::string>("runcontrol.nthreads","-1").c_str(), &endptr, 10);
    cout << "nthreads = " << threads_n << endl;
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'runcontrol.nthreads' needs to be an integer." << endl;
        return -1;
    }
    if (threads_n < 1 && threads_n != -1)
    {
        cout << "ERROR: 'runcontrol.nthreads' needs to be a positive integer." << endl;
        return -1;
    }


    cout << endl << "[ system ]" << endl;
    const int natoms = strtol(pt.get<std::string>("system.natoms","100").c_str(), &endptr, 10);
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
    double temp = strtod(pt.get<std::string>("system.inittemp","1.0").c_str(), &endptr); 
    cout << "inittemp = " << temp << endl;
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'system.inittemp' needs to be a real number." << endl;
        return -1;
    }
    const double rcut = strtod(pt.get<std::string>("system.rcut","9.0").c_str(), &endptr); 
    cout << "rcut = " << rcut << endl;
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'system.rcut' needs to be a real number." << endl;
        return -1;
    }
    const double rlist = strtod(pt.get<std::string>("runcontrol.rlist","10.0").c_str(), &endptr); 
    const double rlist2 = rlist*rlist;
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
    const double reft = strtod(pt.get<std::string>("temperature.reft","1.00").c_str(), &endptr); 
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
    cout << "freq = " << rdf_nbins << endl;
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
    if (threads_n != -1)
    {
        cout << "Number of OpenMP threads set by user." << endl;
        omp_set_dynamic(0);
        omp_set_num_threads(threads_n);
    }
    #pragma omp parallel
    {
        if (omp_get_thread_num() == 0)
        {
            cout << "Using " << omp_get_num_threads() << " OpenMP threads." << endl;
        }
    }
    cout << endl;

    double box_side = pow(natoms/rho,1.0/3.0);
    triclinicbox box(box_side, box_side, box_side);
    cout << "Box is " << box_side << " in each dimension." << endl << endl;

    const double ecut = 4.0 * (1.0/pow(rcut,12) - 1.0/pow(rcut,6));
    const double rcut2 = rcut*rcut;

    double en = 0.0;
    double entot_avg = 0.0;
    double ke = 0.0;
    double ke_avg = 0.0;
    double pe_avg = 0.0;
    double temp_avg = 0.0;

    int nsample = 0;

    vector <coordinates> x(natoms);
    vector <coordinates> xm(natoms);
    vector <coordinates> v(natoms);
    vector <coordinates> f(natoms);
    vector <double> pe_all;
    vector <double> ke_all;
    vector <double> temp_all;
    vector < vector <int> > neighb_list(natoms);

    XDRFILE *xd;
    xd = xdrfile_open(xtcfile.c_str(), "w");

    if (init(x, v, xm, natoms, box, temp, dt, mindist, maxtries, ke) != 0)
    {
        return -1;
    }
    
    Rdf rdf(rdf_nbins, box, rdf_outfile);
    Velocity vel(v_nbins, v_max, v_min, v_outfile);

    update_neighb_list(neighb_list, x, box, rlist2);
    force(f, en, x, box, natoms, rcut2, ecut, neighb_list);

    cout << setw(14) << "Step";
    cout << setw(14) << "Time";
    cout << setw(14) << "Temp";
    cout << setw(14) << "KE";
    cout << setw(14) << "PE";
    cout << setw(14) << "Tot. En." << endl;

    cout << setw(14) << 0;
    cout << setw(14) << 0;
    cout << setw(14) << temp;
    cout << setw(14) << ke;
    cout << setw(14) << en;
    cout << setw(14) << en+ke << endl;

    for (int step = 1; step < nsteps; step++)
    {

        integrate(0, x, v, f, natoms, dt, temp, ke, tcoupl, reft, coll_freq);
        force(f, en, x, box, natoms, rcut2, ecut, neighb_list);
        integrate(1, x, v, f, natoms, dt, temp, ke, tcoupl, reft, coll_freq);

        if (step % nlist == 0)
        {
            update_neighb_list(neighb_list, x, box, rlist2);
        }

        if (( dordf == true) && (step % rdf_freq == 0) && (step > eql_steps))
        {
            rdf.sample(x, box);
        }

        if (( dovel == true) && (step % v_freq == 0) && (step > eql_steps))
        {
            vel.sample(v);
        }

        if ( (step % step_sample) == 0 && (step > eql_steps) )
        {
            nsample++;
            temp_all.push_back(temp);
            temp_avg += temp;
            ke_all.push_back(ke);
            ke_avg += ke;
            pe_all.push_back(en);
            pe_avg += en;
            entot_avg += ke + en;
        }

        if (step % nlog == 0)
        {
            cout << setw(14) << step;
            cout << setw(14) << dt*step;
            cout << setw(14) << temp;
            cout << setw(14) << ke;
            cout << setw(14) << en;
            cout << setw(14) << en+ke << endl;
        }
        if (step % nxtc == 0)
        {
            write_frame(xd, x, box, step, dt);
        }

    }
    xdrfile_close(xd);

    if (dordf == true)
    {
        rdf.normalize(natoms, box);
        rdf.output();
    }

    if (dovel == true)
    {
        vel.normalize(natoms);
        vel.output();
    }

    vector <double> temp_block(nblocks);
    vector <double> pe_block(nblocks);
    vector <double> ke_block(nblocks);

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
            temp_block.at(i) += temp_all.at(j);
            ke_block.at(i) += ke_all.at(j);
            pe_block.at(i) += pe_all.at(j);
        }

        temp_block.at(i) /= (double) (last - first);
        ke_block.at(i) /= (double) (last - first);
        pe_block.at(i) /= (double) (last - first);
    }

    double temp_blockavg = 0.0;
    double ke_blockavg = 0.0;
    double pe_blockavg = 0.0;

    for (int i = 0; i < nblocks; i++)
    {
        temp_blockavg += temp_block.at(i);
        ke_blockavg += ke_block.at(i);
        pe_blockavg += pe_block.at(i);
    }

    temp_blockavg /= nblocks;
    pe_blockavg /= nblocks;
    ke_blockavg /= nblocks;

    double temp_stdev = 0.0;
    double ke_stdev = 0.0;
    double pe_stdev = 0.0;
    for (int i = 0; i < nblocks; i++)
    {
        temp_stdev += pow(temp_block.at(i),2) - pow(temp_blockavg,2);
        ke_stdev += pow(ke_block.at(i),2) - pow(ke_blockavg,2);
        pe_stdev += pow(pe_block.at(i),2) - pow(pe_blockavg,2);
    }

    temp_stdev /= (nblocks-1);
    ke_stdev /= (nblocks-1);
    pe_stdev /= (nblocks-1);

    temp_stdev = sqrt(temp_stdev);
    pe_stdev = sqrt(pe_stdev);
    ke_stdev = sqrt(ke_stdev);

    temp_avg /= (double) nsample;
    ke_avg /= (double) nsample;
    pe_avg /= (double) nsample;
    entot_avg /= (double) nsample;


    cout << "AVERAGES & CONSTANTS (" << nsample << " steps sampled out of " << nsteps << " total steps)" << endl;
    cout << setw(20) << "Number: " << setw(14) << natoms << endl;
    cout << setw(20) << "Density: " << setw(14) << rho << endl;
    cout << setw(20) << "Volume: " << setw(14) << volume(box) << endl;
    cout << setw(20) << "Temperature: " << setw(14) << temp_avg << " +/- " << setw(14) << temp_stdev << endl;
    cout << setw(20) << "Kinetic Energy: " << setw(14) << ke_avg << " +/- " << setw(14) << ke_stdev << endl;
    cout << setw(20) << "Potential Energy: " << setw(14) << pe_avg << " +/- " << setw(14) << pe_stdev << endl;
    cout << setw(20) << "Total Energy: " << setw(14) << entot_avg << endl;

    return 0;
}

int init(vector <coordinates> &x, vector <coordinates> &v, vector <coordinates> &xm, int natoms, triclinicbox box, double &temp, double dt, double mindist, double maxtries, double &ke)
{
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

        x.at(i).at(X) = disx(gen);
        x.at(i).at(Y) = disy(gen);
        x.at(i).at(Z) = disz(gen);

        for (int j = 0; j < i; j++)
        {

            if (distance2(x.at(i), x.at(j), box) < mindist2)
            {
                if (i > maxtries)
                {
                    cout << "ERROR: Exceeded maximum number of tries in generating initial configuration." << endl;
                    return -1;
                }
                goto retrypoint;
            }

        }

        v.at(i).at(X) = dis_vel(gen);
        v.at(i).at(Y) = dis_vel(gen);
        v.at(i).at(Z) = dis_vel(gen);

        sumv += v.at(i);
        sumv2 += dot(v.at(i), v.at(i));

        i++;

    }

    sumv /= natoms;
    sumv2 /= natoms;
    double fs = sqrt(3.0*temp/sumv2);

    sumv2 = 0.0;
    for (int i = 0; i < natoms; i++)
    {
        v.at(i) = (v.at(i) - sumv) * fs;
        sumv2 += dot(v.at(i), v.at(i));
    }
    temp = sumv2 / (3.0 * natoms);
    ke = 0.5 * sumv2 / natoms;

    string file = "init.pdb";
    PdbFile pdb(file.c_str());
    pdb.write_header(file, "WES BARNETT", "First frame");
    for (int i = 0; i < natoms; i++)
    {
        pdb.write_line(i+1, "Ar", "LIG", 1, x.at(i), 1.00, 0.00);
    }
    pdb.close();

    cout << "done." << endl << endl;

    return 0;
}

void force(vector <coordinates> &f, double &en, vector <coordinates> &x, triclinicbox box, int natoms, double rcut2, double ecut, vector < vector <int> > &neighb_list)
{

    en = 0.0;
    for (int i = 0; i < natoms; i++)
    {
        f.at(i) = 0.0;
    }

    #pragma omp parallel
    {

        vector <coordinates> f_thread(natoms);
        double en_thread = 0.0;
        for (int i = 0; i < natoms; i++)
        {
            f_thread.at(i) = 0.0;
        }

        #pragma omp for schedule(guided, 15)
        for (int i = 0; i < natoms; i++)
        {

            if (neighb_list.at(i).size() > 0)
            {
                for (unsigned int j = 0; j < neighb_list.at(i).size(); j++)
                {

                    coordinates dr = pbc(x.at(i) - x.at(neighb_list.at(i).at(j)), box);
                    double r2 = dot(dr,dr);
                    if (r2 <= rcut2)
                    {
                        double r2i = 1.0/r2;
                        double r6i = pow(r2i,3);
                        double ff = 48.0 * r2i * r6i * (r6i - 0.5);
                        f_thread.at(i) += ff * dr;
                        f_thread.at(neighb_list.at(i).at(j)) -= ff * dr;
                        en_thread += 4.0*r6i*(r6i-1.0) - ecut;
                    }

                }
            }

        }

        #pragma omp critical
        for (int i = 0; i < natoms; i++)
        {
            f.at(i) += f_thread.at(i);
            en += en_thread;
        }
    }

    en /= (natoms * natoms);

    return;

}

// Velocity Verlet integrator in two parts
void integrate(int a, vector <coordinates> &x, vector <coordinates> &v, vector <coordinates> &f, int natoms, double dt, double &temp, double &ke, bool tcoupl, double reft, double coll_freq)
{

    if (a == 0) 
    {
        #pragma omp for schedule(guided, 15)
        for (int i = 0; i < natoms; i++)
        {
            x.at(i) = x.at(i) + v.at(i)*dt + 0.5*f.at(i)*dt*dt;
            v.at(i) = v.at(i) + 0.5*f.at(i)*dt;
        }
    }
    else if (a == 1)
    {

        double sumv2 = 0.0;

        // Andersen thermostat. If enabled, each particle is tested to see if it
        // collides with the temperature bath. If it does, it is assigned a
        // velocity drawn from a Gaussian distribution.
        if (tcoupl == true)
        {
            double sigma = sqrt(reft);

            #pragma omp parallel
            {
                random_device rd;
                mt19937 gen(rd());
                uniform_real_distribution<double> dis(0.0, 1.0);
                normal_distribution<double> ndis(0.0, sigma);
                double sumv2_thread = 0.0;

                #pragma omp for schedule(guided, 15)
                for (int i = 0; i < natoms; i++)
                {

                    v.at(i) = v.at(i) + 0.5*f.at(i)*dt;
                    sumv2_thread += dot(v.at(i), v.at(i));

                    if (dis(gen) < coll_freq*dt)
                    {
                        v.at(i).at(X) = ndis(gen);
                        v.at(i).at(Y) = ndis(gen);
                        v.at(i).at(Z) = ndis(gen);
                    }
                }
                
                #pragma omp critical
                {
                    sumv2 += sumv2_thread;
                }

            }

        }
        else
        {
            #pragma omp parallel
            {
                double sumv2_thread = 0.0;

                #pragma omp for schedule(guided, 15)
                for (int i = 0; i < natoms; i++)
                {
                    v.at(i) = v.at(i) + 0.5*f.at(i)*dt;
                    sumv2_thread += dot(v.at(i), v.at(i));
                }

                #pragma omp critical
                {
                    sumv2 += sumv2_thread;
                }

            }

        }

        temp = sumv2 / (3.0 * natoms);
        ke = 0.5 * sumv2 / natoms;

    }

    return;
}

void write_frame(XDRFILE *xd, vector <coordinates> x, triclinicbox &box, int step, double dt)
{
    rvec *x_xtc;
    matrix box_xtc;
    box_xtc[0][0] = box.at(0).at(0)/10.0;
    box_xtc[0][1] = box.at(0).at(1)/10.0;
    box_xtc[0][2] = box.at(0).at(2)/10.0;
    box_xtc[1][0] = box.at(1).at(0)/10.0;
    box_xtc[1][1] = box.at(1).at(1)/10.0;
    box_xtc[1][2] = box.at(1).at(2)/10.0;
    box_xtc[2][0] = box.at(2).at(0)/10.0;
    box_xtc[2][1] = box.at(2).at(1)/10.0;
    box_xtc[2][2] = box.at(2).at(2)/10.0;
    x_xtc = new rvec[x.size()];

    #pragma omp for
    for (unsigned int i = 0; i < x.size(); i++)
    {
        x.at(i) = pbc(x.at(i), box);
        x.at(i).at(X) += box.at(X).at(X)/2.0;
        x.at(i).at(Y) += box.at(Y).at(Y)/2.0;
        x.at(i).at(Z) += box.at(Z).at(Z)/2.0;
        x_xtc[i][X] = x.at(i).at(X)/10.0;
        x_xtc[i][Y] = x.at(i).at(Y)/10.0;
        x_xtc[i][Z] = x.at(i).at(Z)/10.0;
    }

    write_xtc(xd, x.size(), step, dt*step, box_xtc, x_xtc, 1000);

    return;
}

void update_neighb_list(vector < vector <int> > &neighb_list, vector <coordinates> &x, triclinicbox &box, double rlist2)
{
    
    for (unsigned int i = 0; i < neighb_list.size(); i++)
    {
        neighb_list.at(i).resize(0);
    }

    for (unsigned int i = 0; i < x.size()-1; i++)
    {
        for (unsigned int j = i+1; j < x.size(); j++)
        {
            if (distance2(x.at(i), x.at(j), box) < rlist2)
            {
                neighb_list.at(i).push_back(j);
            }

        }
    }

    return;
}
