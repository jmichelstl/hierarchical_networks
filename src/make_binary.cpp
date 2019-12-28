#include <iostream>
#include <stdio.h>
#include <fstream>
#include "network_utils.h"
#include <map>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

//Given a planar line graph specified as pairs of floating-point x and y
//values, create a list of unique points, assign each point an index, create
//a list of index pairs for points that share bonds, and list the relaxed 
//distances between pairs of points that share bonds.
void read_edges(ifstream& datfile, int *num_dof, int *npairs, double **pos, int **pair_list){

    map<int, vector<int>> neighbor_map;
    vector<double> point_list;
    map<Point, int> pmap;
    string nextline;
    Point p1, p2;
    int pindex = 0, index1, index2, mindex, maxdex, point_count, iter;
    vector<double> edge_data;
 
    *npairs = 0;

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
                *npairs = *npairs + 1;

                if(pmap.find(p1) == pmap.end()){
                    pmap.insert(make_pair(p1, pindex));
                    point_list.push_back(p1.x);
                    point_list.push_back(p1.y);
                    neighbor_map.insert(make_pair(pindex, vector<int>()));
                    pindex++;
                }
                if(pmap.find(p2) == pmap.end()){
                    pmap.insert(make_pair(p2, pindex));
                    point_list.push_back(p2.x);
                    point_list.push_back(p2.y);
                    neighbor_map.insert(make_pair(pindex, vector<int>()));
                    pindex++;
                }

                index1 = pmap[p1];
                index2 = pmap[p2];
                mindex = index1 < index2 ? index1 : index2;
                maxdex = index1 == mindex ? index2 : index1;

                neighbor_map[mindex].push_back(maxdex);
            }
        }
    }

    //Establish degrees of freedom and the number of bonds, and allocate
    //arrays used in simulating network mechanics
    *num_dof = point_list.size();
    *pos = (double *) calloc(*num_dof, sizeof(double));
    for(iter = 0; iter < *num_dof; iter++){
        (*pos)[iter] = point_list[iter];
    }

    //Create a list of index pairs. If the topSlide parameter is set, then
    //allow the x coordinates of points along the top and bottom to relax
    *pair_list = (int *) calloc(2 * *npairs, sizeof(int));
    iter = 0;
    for(auto map_iter = neighbor_map.begin(); map_iter != neighbor_map.end(); map_iter++){
        index1 = map_iter->first;

        for(int neighbor : map_iter->second){
            index2 = neighbor;
            (*pair_list)[iter * 2] = index1;
            (*pair_list)[iter * 2 + 1] = index2;
            iter ++;
        }
    }
}

int main(int argc, char **argv){

    ifstream datfile;
    FILE *binfile;
    double *pos;
    int *pair_list;
    int num_dof, num_pairs, count, iter, remaining, max_iter, nv, chunk = 1024;
    bool success;

    if(argc < 2){
        cerr << "Usage: Original data file, binary file name.\n";
        exit(1);
    }
    
    do{
        datfile.open(string(argv[1]), ifstream::in);
        read_edges(datfile, &num_dof, &num_pairs, &pos, &pair_list);
        datfile.close();
        if(num_dof == 0){
            free(pos);
            free(pair_list);
        }
        cout << "DOF: " << num_dof << " Pairs: " << num_pairs << "\n";
    }while(num_dof == 0);

    do{
        success = true;
        binfile = fopen(argv[2], "wb");
        fwrite(&num_dof, sizeof(int), 1, binfile);
        fwrite(&num_pairs, sizeof(int), 1, binfile);
        count = 0;
        remaining = num_dof;
        max_iter = remaining / chunk;
        if(remaining % chunk == 0) max_iter --;
        for(iter = 0; iter <= max_iter; iter ++){
            nv = remaining >= chunk ? chunk : remaining;
            count += fwrite(pos + iter*chunk, sizeof(double), nv, binfile);
            remaining -= nv;
        }
        if(count < num_dof){
            success = false;
            fclose(binfile);
            continue;
        }
        count = 0;
        remaining = 2 * num_pairs;
        max_iter = remaining / chunk;
        if(remaining % chunk == 0) max_iter --;
        for(iter = 0; iter <= max_iter; iter ++){
            nv = remaining >= chunk ? chunk : remaining;
            count += fwrite(pair_list + iter*chunk, sizeof(int), nv, binfile);
            remaining -= nv;
        }
        if(count < 2*num_pairs){
            success = false;
            fclose(binfile);
            continue;
        }
        fclose(binfile);
    }while(! success);

    free(pos);
    free(pair_list);

    return 0;
}
