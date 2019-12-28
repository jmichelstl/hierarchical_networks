#include <iostream>
#include "network_utils.h"
#include <fstream>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <unordered_map>
#include <vector>
#include <string>
#include <cmath>

using namespace std;

void build_angle_lists(ifstream &datfile, vector<vector<double>> &angle_lists){
    unordered_map<Point, vector<double>> point_map;
    string nextline;
    vector<double> edge_data;
    int point_count = 0;
    Point p1, p2;
    double dx, dy, angle;

    while(! datfile.eof()){
        getline(datfile, nextline);
        edge_data = parse_doubles(split(nextline, ' '));
        if(edge_data.size() >= 2){

            if(point_count == 0){
                p1 = Point(edge_data[0], edge_data[1]);
                if(point_map.find(p1) == point_map.end()){
                     point_map.insert(make_pair(p1, vector<double>()));
                }
            }

            else{
                p2 = Point(edge_data[0], edge_data[1]);
                if(point_map.find(p2) == point_map.end()){
                     point_map.insert(make_pair(p2, vector<double>()));
                }
            }

            point_count ++;

            if(point_count == 2){
                point_count = 0;
                dx = p2.x - p1.x;
                dy = p2.y - p1.y;
                angle = atan2(dy, dx);
                point_map[p1].push_back(angle);
                point_map[p2].push_back(angle - M_PI);
            }
        }
    }

    for(auto iter = point_map.begin(); iter != point_map.end(); iter ++){
        angle_lists.push_back(iter->second);
    }
}

double calc_psi_6_mag(vector<vector<double>> angle_lists){
    gsl_complex psi_6_total, local_psi_6, arg;
    double psi_6_mag;

    psi_6_total = gsl_complex_rect(0,0);

    for(vector<double> next_list : angle_lists){
        local_psi_6 = gsl_complex_rect(0,0);
        for(double angle : next_list){
           arg = gsl_complex_rect(0,6*angle);
           local_psi_6 = gsl_complex_add(local_psi_6, gsl_complex_exp(arg));
        }
        local_psi_6 = gsl_complex_div_real(local_psi_6, next_list.size());
        psi_6_total = gsl_complex_add(psi_6_total, local_psi_6);
    }

    psi_6_total = gsl_complex_div_real(psi_6_total, angle_lists.size());
    psi_6_mag = gsl_complex_abs(psi_6_total);

    return psi_6_mag;
}

void process_file(){
    ifstream netfile;
    vector<vector<double>> angle_lists;
    double psi_6_mag;
    int nan_check = 0;

    do{
        open_dat_file("Enter the file containing network data: ", netfile);
        if(! netfile.is_open()){
            cerr << "The file stream was not open for reading.\n";
            if(!yesno("Try again?")) return;
        }
    }while(! netfile.is_open());

    build_angle_lists(netfile, angle_lists);
    netfile.close();

    for(vector<double> next_list : angle_lists){
        for(double num : next_list){
            if(num != num) nan_check ++;
        }
    }
    /*
    cout << "There were " << angle_lists.size() << " entries in all.\n";
    cout << nan_check << " values were not a number.\n";
    */
    psi_6_mag = calc_psi_6_mag(angle_lists);

    cout << psi_6_mag << "\n"; 
}

int main(){
    string response;

    do{
        process_file();
    }while(yesno("Process another file?"));

    return 0;
}
