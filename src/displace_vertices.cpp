#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include "network_utils.h"
#include <string>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;

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

void gauss_random_f(map<int, double>& fmap, double sdev){
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


//Given a point of interest, find the closes point to that one in a set
int get_closest(Point target, map<int, Point> vertex_map){

}

int(main int argc, char **argv){
    string network_name, skeleton_name;
    ifstream network_file, skeleton_file;
    map<int, Point> vertex_map;
    double amp, svdev;

}
