/*
This code performs an energy relaxation of a network using the FIRE algorithm.
*/

#include <iostream>
#include <fstream>
#include "network_utils.h"
#include <map>
#include <unordered_set>
#include <cmath>
#include <vector>
#include <cfloat>
#include <sstream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>

using namespace std;

#define TENSILE 0
#define SHEAR 1

struct EdgeDatum{

    EdgeDatum(int ival, double sval, double len) : index(ival), stiffness(sval), l0(len){
    }

    EdgeDatum(double ival, double len) : index(ival), l0(len){
        stiffness = 1;
    }

    int index;
    double stiffness;
    double l0;
};

void read_edges(ifstream& datfile, vector<double>& point_list, map<int, vector<EdgeDatum>>& neighbor_map, double &miny, double &maxy){

    map<Point, int> pmap;
    string nextline;
    Point p1, p2;
    int pindex = 0, point_count = 0, index1, index2, mindex, maxdex;
    vector<double> edge_data;
    double length;

    miny = FLT_MAX;
    maxy = FLT_MIN;

    while(! datfile.eof()){
        getline(datfile, nextline);

        edge_data = parse_doubles(split(nextline, ' '));
        if(edge_data.size() >= 2){
            if(point_count == 0){
                p1 = Point(edge_data[0], edge_data[1]);
                point_count ++;
            }

            else{
                p2 = Point(edge_data[0], edge_data[1]);
                point_count ++;
            }

            if(point_count == 2){
                point_count = 0;

                if(pmap.find(p1) == pmap.end()){
                    pmap.insert(make_pair(p1, pindex));
                    point_list.push_back(p1.x);
                    point_list.push_back(p1.y);
                    neighbor_map.insert(make_pair(pindex, vector<EdgeDatum>()));
                    pindex++;
                    miny = p1.y < miny - FLOAT_TOL ? p1.y : miny;
                    maxy = p1.y > maxy + FLOAT_TOL ? p1.y : maxy;
                }
                if(pmap.find(p2) == pmap.end()){
                    pmap.insert(make_pair(p2, pindex));
                    point_list.push_back(p2.x);
                    point_list.push_back(p2.y);
                    neighbor_map.insert(make_pair(pindex, vector<EdgeDatum>()));
                    pindex++;
                    miny = p2.y < miny - FLOAT_TOL ? p2.y : miny;
                    maxy = p2.y > maxy + FLOAT_TOL ? p2.y : maxy;
                }

                index1 = pmap[p1];
                index2 = pmap[p2];
                length = sqrt((p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y));
                mindex = index1 < index2 ? index1 : index2;
                maxdex = index1 == mindex ? index2 : index1;

                neighbor_map[mindex].push_back(EdgeDatum(maxdex, length));
            }
        }
    }
}

void gauss_random_d(map<int, double> &d_map, double sdev){
    const gsl_rng_type *T;
    gsl_rng *r;
    FILE *ranfile;
    char num_string[11];
    unsigned seed;

    ranfile = fopen("/dev/urandom", "r");
    fread(&seed, sizeof(unsigned), 1, ranfile);
    fclose(ranfile);
    sprintf(num_string, "%u", seed);
    setenv("GSL_RNG_SEED", num_string, 1);

    gsl_rng_env_setup();
    T = gsl_rng_mt19937;
    r = gsl_rng_alloc(T);

    for(auto iter = d_map.begin(); iter != d_map.end(); iter++){
        d_map[iter->first] = gsl_ran_gaussian(r, sdev);
    }

    gsl_rng_free(r);
}

void gauss_random_f(map<int, double> &f_map, double sdev){
    const gsl_rng_type *T;
    gsl_rng *r;
    double rand_force, x_sum = 0, y_sum = 0;
    map<int, double>::iterator map_iter, last;
    int count = 0;

    map_iter = f_map.begin();
    last = f_map.end();
    if(f_map.size() > 2){
        advance(last, -2);
    }

    gsl_rng_env_setup();
    T = gsl_rng_mt19937;
    r = gsl_rng_alloc(T);

    for(; map_iter != last; map_iter++){
        count ++;
        rand_force = gsl_ran_gaussian(r, sdev);
        if((count % 2) == 1) x_sum -= rand_force;
        else y_sum -= rand_force;
        f_map[map_iter->first] = rand_force;
    }

    //Make the last force minus the sum of all others, so that forces balance
    if(f_map.size() > 2){
        f_map[last->first] = x_sum;
        last ++;
        f_map[last->first] = y_sum;
    }

    count = 0;
    x_sum = y_sum = 0;
    for(auto iter = f_map.begin(); iter != f_map.end(); iter ++){
        count ++;
        if((count % 2) == 1) x_sum += f_map[iter->second];
        else y_sum += f_map[iter->second];
    }

    cout << "X sum: " << x_sum << "\tY sum: " << y_sum << "\n";

    gsl_rng_free(r);
}

void gauss_random_v(vector<double>& vel, double sdev){
    const gsl_rng_type *T;
    gsl_rng *r;
    int iter;

    gsl_rng_env_setup();
    T = gsl_rng_mt19937;
    r = gsl_rng_alloc(T);

    for(iter = 0; iter < vel.size(); iter ++){
        vel[iter] = gsl_ran_gaussian(r, sdev);
    }

    gsl_rng_free(r);
}

void get_forces_plain(map<int, vector<EdgeDatum>> neighbor_map, vector<double> pos, vector<double>& forces){
    int index1, index2;
    double stiffness, x1, x2, y1, y2, l0, dist, fmult, fx, fy;

    forces.assign(pos.size(), 0);

    for(index1 = 0; index1 < pos.size() / 2; index1++){
        x1 = pos[index1*2];
        y1 = pos[index1*2 + 1];
        for(EdgeDatum edat : neighbor_map[index1]){
            index2 = edat.index;
            l0 = edat.l0;
            stiffness = edat.stiffness;
            x2 = pos[2*index2];
            y2 = pos[2*index2 + 1];
            dist = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2) * (y1 - y2));
            fmult = stiffness * (dist - l0)/dist;
            fx = fmult * (x2 - x1);
            fy = fmult * (y2 - y1);
            forces[index1*2] += fx;
            forces[index1*2 + 1] += fy;
            forces[index2*2] -= fx;
            forces[index2*2 + 1] -= fy;
        }
    }
}

void get_forces(map<int, vector<EdgeDatum>> neighbor_map, vector<double> pos, vector<double>& forces, unordered_set<int> bottom, unordered_set<int> top, bool topSlide){
    int index1, index2;
    double stiffness, x1, x2, y1, y2, l0, dist, fmult, fx, fy;

    forces.assign(pos.size(), 0);

    for(index1 = 0; index1 < pos.size() / 2; index1++){
        x1 = pos[index1*2];
        y1 = pos[index1*2 + 1];
        for(EdgeDatum edat : neighbor_map[index1]){
            index2 = edat.index;
            l0 = edat.l0;
            stiffness = edat.stiffness;
            x2 = pos[2*index2];
            y2 = pos[2*index2 + 1];
            dist = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2) * (y1 - y2));
            fmult = stiffness * (dist - l0)/dist;
            fx = fmult * (x2 - x1);
            fy = fmult * (y2 - y1);

            //If the node is not at the bottom or top of the network, relax
            //the x and y coordinates. Otherwise, relax just the x coordinate.
            if(bottom.find(index1) == bottom.end() && top.find(index1) == top.end()){
                forces[index1*2] += fx;
                forces[index1*2 + 1] += fy;
            }
            else if(topSlide){
                forces[index1*2] += fx;
            }
            if(bottom.find(index2) == bottom.end() && top.find(index2) == top.end()){
                forces[index2*2] -= fx;
                forces[index2*2 + 1] -= fy;
            }
            else if(topSlide){
                forces[index2*2] -= fx;
            }
        }
    }
}

void * calc_forces_mt(void *data){
    int index1, index2;
    double stiffness, x1, x2, y1, y2, l0, dist, fmult, fx, fy;
    map<int, vector<EdgeDatum>> *neighbor_map;
    vector<double> *pos;
    unordered_set<int> *bottom, *top;
    bool topSlide;
    int begin, end;
    double *forces;

    neighbor_map = ((map<int, vector<EdgeDatum>> **) data)[0];
    pos = ((vector<double> **) data)[1];
    forces = ((double **) data)[2];
    bottom = ((unordered_set<int> **) data)[3];
    top = ((unordered_set<int> **) data)[4];
    topSlide = *((bool **) data)[5];
    begin = *((int **) data)[6];
    end = *((int **) data)[7];

    for(index1 = begin; index1 <= end; index1++){
        x1 = (*pos)[index1*2];
        y1 = (*pos)[index1*2 + 1];
        for(EdgeDatum edat : (*neighbor_map)[index1]){
            index2 = edat.index;
            l0 = edat.l0;
            stiffness = edat.stiffness;
            x2 = (*pos)[2*index2];
            y2 = (*pos)[2*index2 + 1];
            dist = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2) * (y1 - y2));
            fmult = stiffness * (dist - l0)/dist;
            fx = fmult * (x2 - x1);
            fy = fmult * (y2 - y1);

            //If the node is not at the bottom or top of the network, relax
            //the x and y coordinates. Otherwise, relax just the x coordinate.
            if((*bottom).find(index1) == (*bottom).end() && (*top).find(index1) == (*top).end()){
                forces[index1*2] += fx;
                forces[index1*2 + 1] += fy;
            }
            else if(topSlide){
                forces[index1*2] += fx;
            }
            if((*bottom).find(index2) == (*bottom).end() && (*top).find(index2) == (*top).end()){
                forces[index2*2] -= fx;
                forces[index2*2 + 1] -= fy;
            }
            else if(topSlide){
                forces[index2*2] -= fx;
            }
        }
    }

    return NULL;
}

//Given a number of worker threads, determine assignments of pairwise force
//calculation
vector<int> assign_indices(map<int, vector<EdgeDatum>> neighbor_map, int nthreads){

    int total_interactions, num_per_thread, thread_num, threshold;
    map<int, int> count_map;
    map<int, int>::iterator count_iter;
    vector<int> index_list;

    total_interactions = 0;
    //First find the total number of interactions
    for(auto iter = neighbor_map.begin(); iter != neighbor_map.end(); iter++){
        total_interactions += iter->second.size();
        count_map.insert(make_pair(total_interactions,iter->first));
    }

    num_per_thread = total_interactions / nthreads;
    index_list.push_back(0);

    for(thread_num = 1; thread_num < nthreads; thread_num++){
        count_iter = count_map.upper_bound(thread_num * num_per_thread);
        index_list.push_back(count_iter->second);
        index_list.push_back(count_iter->second + 1);
    }

    index_list.push_back(neighbor_map.rbegin()->first);

    return index_list;
}

void get_forces_mt(vector<double> &forces, int nthreads, void **worker_data, double **force_lists){

    int indexIter, offset, forceIter;
    pthread_t *tids;
    pthread_t tid;
    pthread_attr_t attr;

    forces.assign(forces.size(), 0);
    for(indexIter = 0; indexIter < nthreads; indexIter ++){
        for(forceIter = 0; forceIter < forces.size(); forceIter ++){
            force_lists[indexIter][forceIter] = 0;
        }
    }

    //Prepare pthread structures
    pthread_attr_init(&attr);
    tids = (pthread_t *) calloc(nthreads, sizeof(pthread_t));

    //Assign calculations and wait for each thread to complete
    for(indexIter = 0; indexIter < nthreads; indexIter ++){
        pthread_create(&tids[indexIter], &attr, &calc_forces_mt, (void *) (worker_data + indexIter*8));
    }

    for(indexIter = 0; indexIter < nthreads; indexIter++){
        pthread_join(tids[indexIter], NULL);
    }

    for(indexIter = 0; indexIter < nthreads; indexIter++){
        for(forceIter = 0; forceIter < forces.size(); forceIter ++){
            forces[forceIter] += force_lists[indexIter][forceIter];
        }
    }

    free(tids);
    pthread_attr_destroy(&attr);
}

//Prepare data describing division of labor for calculating forces
void mt_force_workspace_init(map<int, vector<EdgeDatum>> *neighbor_map, vector<double> *pos, unordered_set<int> *bottom, unordered_set<int> *top, bool *topSlide, int nthreads, int *indices, void ***worker_data, double ***force_lists){

    int indexIter, offset;
    vector<int> index_list;

    //Determine division of labor among worker threads
    index_list = assign_indices(*neighbor_map, nthreads);
    indices = (int *) calloc(index_list.size(), sizeof(int));
    for(indexIter = 0; indexIter < index_list.size(); indexIter++){
        indices[indexIter] = index_list[indexIter];
    }

    //Allocate workspace and determine division of labor among threads
    *worker_data = (void **) calloc(nthreads * 8, sizeof(void *));

    *force_lists = (double **) calloc(nthreads, sizeof(double *));
    for(indexIter = 0; indexIter < nthreads; indexIter ++){
        (*force_lists)[indexIter] = (double *) calloc((*pos).size(), sizeof(double));
    }

    //Prepare force calculation assignments for worker threads
    for(indexIter = 0; indexIter < nthreads; indexIter++){
        offset = indexIter*8;
        (*worker_data)[offset] = (void *) neighbor_map;
        (*worker_data)[offset + 1] = (void *) pos;
        (*worker_data)[offset + 2] = (void *) (*force_lists)[indexIter];
        (*worker_data)[offset + 3] = (void *) bottom;
        (*worker_data)[offset + 4] = (void *) top;
        (*worker_data)[offset + 5] = (void *) topSlide;
        (*worker_data)[offset + 6] = (void *) (indices + indexIter*2);
        (*worker_data)[offset + 7] = (void *) (indices + indexIter*2 + 1);
    }
}

//Free space allocated for multi-threaded calculation of forces
void mt_force_workspace_destroy(int nthreads, int *indices, void **worker_data, double **force_lists){

    int indexIter;

    for(indexIter = 0; indexIter < nthreads; indexIter++){
        free(force_lists[indexIter]);
    }

    free(indices);
    free(force_lists);\
    free(worker_data);
}

double mag(vector<double> vec){
    double mag = 0;
    int index;

    for(index = 0; index < vec.size(); index++){
        mag += vec[index] * vec[index];
    }

    return sqrt(mag);
}

double get_pe(map<int, vector<EdgeDatum>> neighbor_map, vector<double> pos_vec){
    double pe = 0;
    int num_points = pos_vec.size() / 2, index1, index2;
    double distance, diff, x1, y1, x2, y2;

    for(index1 = 0; index1 < num_points; index1 ++){
        x1 = pos_vec[2*index1];
        y1 = pos_vec[2*index1 + 1];
        for(EdgeDatum edat : neighbor_map[index1]){
            index2 = edat.index;
            x2 = pos_vec[2*index2];
            y2 = pos_vec[2*index2 + 1];
            distance = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
            diff = distance - edat.l0;
            pe += .5 * edat.stiffness * diff * diff;
            //cout << "k, l0, distance, diff y1, y2: " << edat.stiffness << " " << edat.l0 << " " << distance << " " << diff << " " << y1 << " " << y2 << "\n";
        }
    }

    return pe;
}

int get_dec_digits(int value){
    int num_digits = 0;

    do{
        num_digits++;
        value /= 10;
    }while(value > 0);

    return num_digits;
}

void get_report_params(string& base, int& rep_freq){

    string response;
    int num_read;

    do{
        cout << "Enter base name: ";
        getline(cin, response);
    }while(response.compare("") == 0);

    base = response;

    do{
        cout << "Enter reporting frequency: ";
        getline(cin, response);
        num_read = sscanf(response.c_str(), "%d", &rep_freq);
    }while(num_read < 1 || rep_freq < 1);
}

string report_name(string base, int num_digits, int count){
    int cdigits, padding, iter;
    ostringstream oss;

    padding = num_digits - get_dec_digits(count);
    oss << base << "_";
    for(iter = 0; iter < padding; iter++) oss << 0;
    oss << count << ".dat";

    return oss.str();
}

void report_deformed(map<int, vector<EdgeDatum>> neighbor_map, vector<double> pos, string filename){
    ofstream datfile;
    int index1, index2;

    datfile.open(filename);

    for(index1 = 0; index1 < pos.size(); index1 ++){
        for(EdgeDatum edat : neighbor_map[index1]){
            index2 = edat.index;
            datfile << pos[2*index1] << "\t" << pos[2*index1 + 1] << "\n";
            datfile << pos[2*index2] << "\t" << pos[2*index2 + 1] << "\n\n";
        }
    }

    datfile.close();
}

void simple_md(vector<double>& pos, vector<double>& vel, map<int, vector<EdgeDatum>> neighbor_map, unordered_set<int> bottom, unordered_set<int> top, double mass, double dt, int num_steps){

    vector<double> forces;
    int step_count, index;
    double dt_sq = dt*dt;

    for(step_count = 0; step_count < num_steps; step_count ++){
        get_forces(neighbor_map, pos, forces, bottom, top, false);
        for(index = 0; index < pos.size(); index++){
            if(bottom.find(index/2) == bottom.end() && top.find(index/2) == top.end()){
                pos[index] += dt*vel[index] + .5 * forces[index]/mass * dt_sq;
                vel[index] += .5 * dt * forces[index] / mass;
            }
        }

        get_forces(neighbor_map, pos, forces, bottom, top, false);
        for(index = 0; index < pos.size(); index++){
            if(bottom.find(index/2) == bottom.end() && top.find(index/2) == top.end()){
                vel[index] += .5 * dt * forces[index] / mass;
            }
        }

    }
}

int run_fire_md(vector<double>& pos, vector<double>& vel, map<int, vector<EdgeDatum>> neighbor_map, unordered_set<int> bottom, unordered_set<int> top, double mass, double fcut, int max_steps, int nmin, double fire_params[5], bool& flag, bool topSlide, int nthreads, void **worker_data, double **force_lists){

    vector<double> forces;
    double alpha, finc, fdec, alpha_start, falpha, fire_dt, fire_dt_max;
    int num_points, num_dof, step_count = 0, since_leq_0 = 0, index, freq;
    int rep_freq, digits, rep_count = 1, indexIter;
    double power, dt, dt_sq, vmag, fmag, inv_fmag, sqrt_inner_dof;
    string base, full_name;
    bool report;

    /*if((report = yesno("Report intermediate state?"))){
        get_report_params(base, rep_freq);
        digits = get_dec_digits(max_steps / rep_freq);
    }*/

    //Find the number of points and initialize force vector
    num_dof = pos.size();
    num_points = num_dof / 2;
    sqrt_inner_dof = sqrt(num_dof - bottom.size() - top.size());

    //Unpack parameters of FIRE minimization scheme and initialize values
    alpha_start = fire_params[0];
    falpha = fire_params[1];
    fire_dt_max = fire_params[2];
    finc = fire_params[3];
    fdec = fire_params[4];

    alpha = alpha_start;
    dt = fire_dt_max;
    dt_sq = dt*dt;
   
    //Find the forces at the outset
    forces.assign(num_dof, 0);
    if(nthreads == 1) get_forces(neighbor_map, pos, forces, bottom, top, topSlide);
    else get_forces_mt(forces, nthreads, worker_data, force_lists);
    vmag = mag(vel);

    cout << "Starting RMS Force: " << mag(forces) / sqrt_inner_dof << "\n";
 
    //Perform molecular dynamics steps using velocity verlet method until the
    //kinetic energy cutoff is reached, or the maximum number of steps have
    //taken place
    while(step_count < max_steps){
        step_count ++;        

        //Update positions
        for(index = 0; index < num_dof; index++){
            if((bottom.find(index/2) == bottom.end() && top.find(index/2) == top.end()) || ((index % 2) == 0 && topSlide)){
                pos[index] += dt*vel[index] + .5 * forces[index]/mass * dt_sq;
                vel[index] += .5 * dt * forces[index] / mass;
            }
        }

        //Calculate forces
        if(nthreads == 1) get_forces(neighbor_map, pos, forces, bottom, top, topSlide);
        else get_forces_mt(forces, nthreads, worker_data, force_lists);

        //Update velocities and calculate power
        power = 0;
        for(index = 0; index < num_dof; index++){
            if((bottom.find(index/2) == bottom.end() && top.find(index/2) == top.end()) || ((index % 2) == 0 && topSlide)){
                vel[index] += .5 * dt * forces[index] / mass;
                power += vel[index] * forces[index];
            }
        }

        //Adjust velocities according to FIRE algorithm
        vmag = mag(vel);
        fmag = mag(forces);
        inv_fmag = 1 / fmag;

        for(index = 0; index < num_dof; index++){
            if((bottom.find(index) == bottom.end() && top.find(index) == top.end()) || ((index % 2) == 0 && topSlide)){
                vel[index] += alpha*(vmag*forces[index]*inv_fmag - vel[index]);
            }
        }
        
        //Adjust FIRE parameters according to current power
        if(power > 0){
            since_leq_0 ++;
            if(since_leq_0 > nmin){
                dt = min(dt*finc, fire_dt_max);
                dt_sq = dt*dt;
                alpha *= falpha;
            }
        }
        else{
            since_leq_0 = 0;
            dt *= fdec;
            dt_sq = dt * dt;
            alpha = alpha_start;
            vel.assign(num_dof, 0);
        }

        /*if(report && (step_count % rep_freq == 0)){
            full_name = report_name(base, digits, rep_count);
            report_deformed(neighbor_map, pos, full_name);
            rep_count ++;
        }*/

        //Check for kinetic energy convergence
        if(fmag / sqrt_inner_dof < fcut){
            flag = true;
            cout << "RMS Force: " << fmag / sqrt_inner_dof << "\n";
            return step_count;
        }
    }

    flag = false;
    cout << "RMS Force: " << fmag / sqrt_inner_dof << "\n";
    cout << "Ending energy: " << .5 * mass * vmag * vmag + get_pe(neighbor_map, pos) << "\n";

    return step_count;
}

void displace(vector<double>& pos_vec, unordered_set<int> top, double disp, int dmode){

    if(dmode == TENSILE){
        for(int index : top){
            pos_vec[index*2 + 1] += disp;
        }
    }
    else{
        for(int index : top){
            pos_vec[index*2] += disp;
        }
    }
}

void displace_2d(vector<double> &pos_vec, map<int, double> d_map){
    for(auto iter = d_map.begin(); iter != d_map.end(); iter++){
        pos_vec[iter->first] += iter->second;
    }
}

//Performs a series of relaxation runs with varying displacement to obtain
//a curve of energy vs. displacement to measure stiffness
void e_vs_d(){
}

//Searches for the strain at which a network begins to fracture
//Provides critical strain and energy storage at fracture
void fracture_search(){
}

//Given a set of candidate points and a target point, find the point within the
//candidate set closest to the target point
int get_closest(set<Point> candidates, Point target){
    set<Point>::iterator lbound, curr;
    double dist, min_dist = FLT_MAX, x, y;
    int best_index;

    //Get the first point in the candidate set not considered to be below the
    //target point according to the ordering scheme for Point structs
    lbound = candidates.lower_bound(target);

    //Start the search from the first point not considered less than the target
    //End the search when the vertical separation is equal to or greater than
    //the minimum separation found so far
    curr = lbound;
    while(curr != candidates.end()){
        x = (*curr).x;
        y = (*curr).y;
        dist = sqrt((x-target.x)*(x-target.x) + (y-target.y)*(y-target.y));
        if(dist < min_dist){
            min_dist = dist;
            best_index = (*curr).index;
        }
        if(abs(y - target.y) < min_dist) break;
        advance(curr, 1);
    }

    //Repeat the procedure above starting just below the "lower bound"
    curr = lbound;
    while(curr != candidates.begin()){
        advance(curr, -1);
        x = (*curr).x;
        y = (*curr).y;
        dist = sqrt((x-target.x)*(x-target.x) + (y-target.y)*(y-target.y));
        if(dist < min_dist){
            min_dist = dist;
            best_index = (*curr).index;
        }
        if(abs(y - target.y) < min_dist) break;
    }

    return best_index;
}

//Given a small-scale network and a large-scale, skeletal network, find the
//small-scale node closest to each skeleton node
void approx_skel_nodes(vector<double> point_list, vector<double> targets, map<int, double> &d_forces){
    set<Point> small_points;
    Point next_large;
    int iter, best_match;

    for(int i = 0; i < point_list.size() - 1; i += 2){
        small_points.insert(Point(point_list[i], point_list[i+1], i/2));
    }

    for(iter = 0; iter < targets.size() - 1; iter += 2){
        next_large = Point(targets[iter], targets[iter+1]);
        best_match = get_closest(small_points, next_large);
	d_forces.insert(make_pair(best_match*2, 0));
        d_forces.insert(make_pair(best_match*2 + 1, 0));
    }
}

int fire_w_forces(vector<double>& pos, vector<double>& vel, map<int, vector<EdgeDatum>> neighbor_map, map<int, double> f_map, double mass, double fcut, int max_steps, int nmin, double fire_params[5], bool& flag){

    vector<double> forces;
    double alpha, finc, fdec, alpha_start, falpha, fire_dt, fire_dt_max;
    int num_points, num_dof, step_count = 0, since_leq_0 = 0, index, freq;
    int rep_freq, digits, rep_count = 1;
    double power, dt, dt_sq, vmag, fmag, inv_fmag, sqrt_inner_dof;
    string base, full_name;
    bool report;

    if((report = yesno("Report intermediate state?"))){
        get_report_params(base, rep_freq);
        digits = get_dec_digits(max_steps / rep_freq);
    }

    //Find the number of points and initialize force vector
    num_dof = pos.size();
    num_points = num_dof / 2;
    sqrt_inner_dof = sqrt(num_dof);

    //Unpack parameters of FIRE minimization scheme and initialize values
    alpha_start = fire_params[0];
    falpha = fire_params[1];
    fire_dt_max = fire_params[2];
    finc = fire_params[3];
    fdec = fire_params[4];

    alpha = alpha_start;
    dt = fire_dt_max;
    dt_sq = dt*dt;
   
    //Find the forces at the outset
    forces.assign(num_dof, 0);
    get_forces_plain(neighbor_map, pos, forces);
    for(auto map_iter = f_map.begin(); map_iter != f_map.end(); map_iter ++){
        forces[map_iter->first] += map_iter->second;
    }
    vmag = mag(vel);

    cout << "Starting RMS Force: " << mag(forces) / sqrt_inner_dof << "\n";
 
    //Perform molecular dynamics steps using velocity verlet method until the
    //kinetic energy cutoff is reached, or the maximum number of steps have
    //taken place
    while(step_count < max_steps){
        step_count ++;        

        //Update positions
        for(index = 0; index < num_dof; index++){
            pos[index] += dt*vel[index] + .5 * forces[index]/mass * dt_sq;
            vel[index] += .5 * dt * forces[index] / mass;
        }

        //Calculate forces
        get_forces_plain(neighbor_map, pos, forces);
        for(auto map_iter = f_map.begin(); map_iter != f_map.end(); map_iter ++){
            forces[map_iter->first] += map_iter->second;
        }

        //Update velocities and calculate power
        power = 0;
        for(index = 0; index < num_dof; index++){
            vel[index] += .5 * dt * forces[index] / mass;
            power += vel[index] * forces[index];
        }

        //Adjust velocities according to FIRE algorithm
        vmag = mag(vel);
        fmag = mag(forces);
        inv_fmag = 1 / fmag;

        for(index = 0; index < num_dof; index++){
            vel[index] += alpha*(vmag*forces[index]*inv_fmag - vel[index]);
        }
        

        //Adjust FIRE parameters according to current power
        if(power > 0){
            since_leq_0 ++;
            if(since_leq_0 > nmin){
                dt = min(dt*finc, fire_dt_max);
                dt_sq = dt*dt;
                alpha *= falpha;
            }
        }
        else{
            since_leq_0 = 0;
            dt *= fdec;
            dt_sq = dt * dt;
            alpha = alpha_start;
            vel.assign(num_dof, 0);
        }

        if(report && step_count % rep_freq == 0){
            full_name = report_name(base, digits, rep_count);
            report_deformed(neighbor_map, pos, full_name);
            rep_count ++;
        }

        //Check for kinetic energy convergence
        if(fmag / sqrt_inner_dof < fcut){
            flag = true;
            cout << "RMS Force: " << fmag / sqrt_inner_dof << "\n";
            return step_count;
        }
    }

    flag = false;
    cout << "RMS Force: " << fmag / sqrt_inner_dof << "\n";
    cout << "Ending energy: " << .5 * mass * vmag * vmag + get_pe(neighbor_map, pos) << "\n";
    return step_count;
}

void do_disorder_run(bool abort, int nthreads){

    ifstream net_file, skel_file;
    vector<double> small_points, skel_points, run_params, vel;
    map<int, double> d_map;
    double mass, fcut, mag, sdev, miny, maxy;
    double fparams[5];
    string response, new_net, new_skel;
    int max_steps, scount, rep_count, numd, skel_iter, num_read, nmin;
    int num_incr, inc_count;
    map<int, vector<EdgeDatum>> small_scale_map, skel_map;
    unordered_set<int> fixed_points, dummy;
    bool flag;
    map<int, double>::iterator dmap_iter, last;
    vector<int> index_list;
    int *indices = NULL;
    void **worker_data = NULL;
    double **force_lists = NULL;
    bool topSlide = false;

    //Obtain small-scale and skeletal networks
    open_dat_file("Enter the small-scale network file: ", net_file);
    if(! net_file.is_open()){
        cerr << "No small-scale file was read. The attempt will be aborted.\n";
        return;
    }
    read_edges(net_file, small_points, small_scale_map, miny, maxy);
    net_file.close();

    if(nthreads > 1){
        mt_force_workspace_init(&small_scale_map, &small_points, &fixed_points, &dummy, &topSlide, nthreads, indices, &worker_data, &force_lists);
    }
 
    open_dat_file("Enter the skeleton file: ", skel_file);
    if(! skel_file.is_open()){
        cerr << "No large-scale file was read. The attempt will be aborted.\n";
        return;
    }
    read_edges(skel_file, skel_points, skel_map, miny, maxy);
    skel_file.close();


    //If either the list of skeletal points or small-scale points is empty,
    //abort.
    if(small_points.size() == 0 || skel_points.size() == 0){
        cerr << "Not enough information was obtained to build two networks.\n";
        return;
    }

    //Determine the small-scale nodes that best approximate each large-scale
    //node
    approx_skel_nodes(small_points, skel_points, d_map);

    last = d_map.end();
    advance(last, -1);
    dmap_iter = d_map.begin();
    do{
        fixed_points.insert(dmap_iter->first / 2);
        advance(dmap_iter, 2);
    }while(dmap_iter != last && dmap_iter != d_map.end());

    //Assign random forces to selected nodes
    sdev = 0;
    do{
        run_params = getdoubles("Enter displacement standard deviation: ");
        if(run_params.size() < 1) cerr << "Enter a number.\n";
        else if(run_params[0] <= 0) cerr << "Enter a positive number.\n";
        else sdev = run_params[0];
    }while(sdev < 0);

    do{
        cout << "Enter the number of displacement increments: ";
        getline(cin, response);
        num_read = sscanf(response.c_str(), "%d", &num_incr);
        if(num_read == 0) cerr << "Enter a number.\n";
        if(num_read == 1 && num_incr < 1){
            cerr << "Enter a positive number.\n";
            num_read = 0;
        }
    }
    while(num_read < 1);

    //Obtain parameters for a FIRE run
    do{
        cout << "Enter alpha, falpha, dt, finc, fdec, and nmin: ";
        getline(cin, response);
        num_read = sscanf(response.c_str(),"%lf %lf %lf %lf %lf %d", fparams, fparams+1, fparams+2, fparams+3, fparams+4, &nmin);
    }while(num_read < 6);

    do{
        cout << "Enter particle mass, cutoff force, and maximum md steps: ";
        getline(cin, response);
        num_read = sscanf(response.c_str(), "%lf %lf %d", &mass, &fcut, &max_steps);
    }while(num_read < 3);

    vel.assign(small_points.size(), 0);

    for(inc_count = 0; inc_count != num_incr; inc_count ++){
        scount = 0;
        flag = false;
        gauss_random_d(d_map, sdev / num_incr);
        displace_2d(small_points, d_map);
        while(! flag){
            scount += run_fire_md(small_points, vel, small_scale_map, fixed_points, dummy, mass, fcut, max_steps, nmin, fparams, flag, false, nthreads, worker_data, force_lists);
            //scount += fire_w_forces(small_points, vel, small_scale_map, d_forces, mass, fcut, max_steps, nmin, fparams, flag);
            if(! flag){
                cout << "Convergence not reached after "<< scount<< " steps.\n";
                if(abort) break;
                else if(!yesno("Continue?")) break;
            }
        }
    }

    if(yesno("Report final state?")){
        cout << "Enter the skeleton file: ";
        getline(cin, new_skel);
        skel_iter = 0;
        for(auto d_iter = d_map.begin(); d_iter != d_map.end(); d_iter++){
            skel_points[skel_iter] = small_points[d_iter->first];
            skel_iter ++;
        }
        report_deformed(skel_map, skel_points, new_skel);

        cout << "Enter the small-scale file: ";
        getline(cin, new_net);
        report_deformed(small_scale_map, small_points, new_net);
    }

    if(nthreads > 1){
        mt_force_workspace_destroy(nthreads, indices, worker_data, force_lists);
    }
}

void do_relaxation_run(bool abort, int dmode, bool topSlide, int nthreads){

    ifstream netfile;
    vector<double> point_list, pos_copy, vel, sim_results;
    map<int, vector<EdgeDatum>> neighbor_map;
    double miny, maxy, mass, fcut, dmin, dmax, dinc, disp, pe, pre_dt, sdev;
    unordered_set<int> bottom, top;
    int i, num_dof, point_count, index1, index2, num_read, max_steps, nmin;
    double fparams[5];
    string response, report_base, report_file;
    bool converged, precondition, randv, batch_report;
    int pre_steps, scount, rep_count, numd;
    ofstream ereport;
    ostringstream report_line;
    vector<int> index_list;
    int *indices = NULL;
    void **worker_data;
    double **force_lists;

    //Obtain file with information about the network
    open_dat_file("Enter the file containing network data: ", netfile);
    if(! netfile.is_open()){
        cerr << "The file stream was not open for reading.\n";
        return;
    }

    //Obtain information about vertex locations and connections in the
    //undeformed lattice
    read_edges(netfile, point_list, neighbor_map, miny, maxy);
    num_dof = point_list.size();
    point_count = num_dof / 2;
    netfile.close();

    //Find the points in the top and bottom sections, so they can be held fixed
    for(i = 0; i < point_count; i++){
        if(abs(point_list[2*i + 1] - miny) < FLOAT_TOL) bottom.insert(i);
        if(abs(point_list[2*i + 1] - maxy) < FLOAT_TOL) top.insert(i);
    }

    //Get FIRE parameters for adjusting velocity
    do{
        cout << "Enter alpha, falpha, dt, finc, fdec, and nmin: ";
        getline(cin, response);
        num_read = sscanf(response.c_str(),"%lf %lf %lf %lf %lf %d", fparams, fparams+1, fparams+2, fparams+3, fparams+4, &nmin);
    }while(num_read < 6);

    do{
        cout << "Enter particle mass, cutoff force, and maximum md steps: ";
        getline(cin, response);
        num_read = sscanf(response.c_str(), "%lf %lf %d", &mass, &fcut, &max_steps);
    }while(num_read < 3);

    //Get range of displacements
    do{
        cout << "Enter minimum and maximum displacement and increment: ";
        getline(cin, response);
        num_read = sscanf(response.c_str(), "%lf %lf %lf", &dmin, &dmax, &dinc);
    }while(num_read < 3);

    //Get preconditioning parameters
    if((precondition = yesno("Do a preconditioning run?"))){
        do{
            cout << "Enter preconditioning time step and number of steps: ";
            getline(cin, response);
            num_read = sscanf(response.c_str(), "%lf %d", &pre_dt, &pre_steps);
        }while(num_read < 2);
    }

    //Get information for optional initial random velocities
    if((randv = yesno("Add random velocities?"))){
        do{
            cout << "Enter the standard deviation: ";
            getline(cin, response);
            num_read = sscanf(response.c_str(), "%lf", &sdev);
        }while(num_read < 1 || sdev <= 0);
    }

    //Give the option to do batch reporting of final networks
    if((batch_report = yesno("Do batch reporting of final networks?"))){
        do{
            cout << "Enter a base name for reporting: ";
            getline(cin, report_base);
        }while(report_base.compare("") == 0);
        numd = get_dec_digits((int) (1 + (dmax - dmin)/dinc));
    }

    pos_copy.assign(point_list.begin(), point_list.end());
    if(nthreads > 1){
        mt_force_workspace_init(&neighbor_map, &pos_copy, &bottom, &top, &topSlide, nthreads, indices, &worker_data, &force_lists);
    }

    displace(pos_copy, top, dmin, dmode);
    rep_count = 1;
    for(disp = dmin; disp <= dmax; disp += dinc){
        vel.assign(num_dof, 0);
        scount = 0;

        do{

            //Create optional random velocities
            if(randv){
                gauss_random_v(vel, sdev);
            }

            //Do optional MD run without FIRE minimization
            if(precondition){
                simple_md(pos_copy, vel, neighbor_map, bottom, top, mass, pre_dt, pre_steps);
            }

            scount += run_fire_md(pos_copy, vel, neighbor_map, bottom, top, mass, fcut, max_steps, nmin, fparams, converged, topSlide, nthreads, worker_data, force_lists);
            if(! converged && !abort){
                cerr << "Convergence not reached in " << scount << " steps.\n";
                if(! yesno("Continue?")) break;
            }
            else{
                if(converged) cout << "Convergence reached after " << scount << " steps.\n";
                else cout << "Covergence not reached after " << scount << " steps.\n";
            }
        }while(! converged && ! abort);

        if(batch_report){// || yesno("Report final network?")){

            if(! batch_report){
                do{
                    cout << "Enter a file name: ";
                    getline(cin, report_file);
                }while(report_file.compare("") == 0);
            }
            else report_file = report_name(report_base, numd, rep_count);

            report_deformed(neighbor_map, pos_copy, report_file);
        }

        pe = get_pe(neighbor_map, pos_copy);
        sim_results.push_back(disp);
        sim_results.push_back(pe);
        displace(pos_copy, top, dinc, dmode);
        rep_count ++;
    }

    cout << "Enter a file name for reporting energies, or return for no file: ";
    getline(cin, response);
    if(! response.compare("") == 0){
        ereport.open(response, ofstream::trunc);
    }

    for(i = 0; i < sim_results.size() / 2; i++){
        report_line << sim_results[2*i] << "\t" << sim_results[2*i + 1] << "\n";
    }

    if(ereport.is_open()) ereport << report_line.str();
    else cout << report_line.str();
    if(ereport.is_open()) ereport.close();

    if(nthreads > 1){
        mt_force_workspace_destroy(nthreads, indices, worker_data, force_lists);
    }
}

int main(int argc, char **argv){

    string filename, nextline;
    ifstream netfile;
    double min_disp, max_disp, step;
    char c;
    bool abort = false, disorder = false, topSlide = false;
    int dmode = TENSILE, nthreads = 1, num_read;

    while((c = getopt(argc, argv, "adn:st")) != -1){
        switch(c) {
            case 'a':
                abort = true;
                break;
            case 'd':
                disorder = true;
                break;
            case 'n':
                num_read = sscanf(optarg, "%d", &nthreads);
                if(num_read < 1 || nthreads < 1){
                    cerr << "Option \"n\" requires a positive integer argument.\n";
                    nthreads = 1;
                }
                break;
            case 's':
                dmode = SHEAR;
                break;
            case 't':
                topSlide = true;
                break;
            case '?':
                if(optopt == 'n'){
                    cerr << "Option \"n\" requires a positive integer argument.\n";
                }
                else if(isprint(optopt)){
                    fprintf(stderr, "Unknown option: -%c.\n", optopt);
                }
                else{
                    fprintf(stderr, "Unknown option character.\n");
                }
            default:
                break;
         }
    }

    while(yesno("Perform a relaxation run?")){
        if(disorder) do_disorder_run(abort, nthreads);
        else do_relaxation_run(abort, dmode, topSlide, nthreads);
    }
}
