#include "network_utils.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <ctype.h>
#include <stdlib.h>
#include <unistd.h>

using namespace std;

#define tuple_2d tuple<double, double>

struct EdgeDatum{

    EdgeDatum(int ival, double sval, double len) : index(ival), stiffness(sval), l0(len){
    }

    EdgeDatum(double ival, double len) : index(ival), l0(len){}

    int index;
    double stiffness;
    double l0;
};

bool contains(vector<char> list, char target){
    for(char c : list){
        if(c == target) return true;
    }

    return false;
}

vector<string> split_string(string input, vector<char> delims){
    vector<string> result;
    size_t start, iter, len;
    bool reading_token = false;

    for(iter = 0; iter < input.size(); iter++){
        if(! contains(delims, input[iter])){
            if(!reading_token){
                reading_token = true;
                start = iter;
            }
        }

        else{
            if(reading_token){
                reading_token = false;
                result.push_back(input.substr(start, iter - start));
            }
        }
    }

    if(reading_token) result.push_back(input.substr(start, iter - start));


    return result;
}

vector<Edge> read_edges(ifstream &edge_file){
    Point p1, p2;
    int point_count = 0;
    vector<double> edge_data;
    string nextline;
    vector<Edge> edge_list;
    vector<char> delims = {' ', '\t', '\n'};

    while(!edge_file.eof()){
        getline(edge_file, nextline);
        edge_data = parse_doubles(split_string(nextline, delims));
        if(edge_data.size() >= 2){
            if(point_count == 0){
                p1 = Point(edge_data[0], edge_data[1]);
            }
            else{
                p2 = Point(edge_data[0], edge_data[1]);
            }

            point_count ++;
            if(point_count == 2){
                point_count = 0;
                edge_list.push_back(Edge(p1, p2));
            }
        }
    }

    edge_file.close();
    return edge_list;
}

void read_relax_style(ifstream& datfile, vector<double>& point_list, map<int, vector<EdgeDatum>>& neighbor_map){

    map<Point, int> pmap;
    string nextline;
    Point p1, p2;
    int pindex = 0, point_count = 0, index1, index2, mindex, maxdex;
    vector<double> edge_data;
    double length;

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
                }
                if(pmap.find(p2) == pmap.end()){
                    pmap.insert(make_pair(p2, pindex));
                    point_list.push_back(p2.x);
                    point_list.push_back(p2.y);
                    neighbor_map.insert(make_pair(pindex, vector<EdgeDatum>()));
                    pindex++;
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

    datfile.close();
}

map<Point, tuple_2d> build_displacement_map(vector<double> point_list, map<int, vector<EdgeDatum>> neighbor_map, vector<Edge> strained_edges){

    int ref_index, point_count = point_list.size() / 2;
    map<Point, tuple_2d> disp_map;
    Point p1, p2, pi;
    int edge_list_index = 0;

    for(auto iter = neighbor_map.begin(); iter != neighbor_map.end(); iter++){
        ref_index = iter->first;
        p1 = strained_edges[edge_list_index].p1;

        pi = Point(point_list[ref_index*2], point_list[ref_index*2+1]);
        if(disp_map.find(pi) == disp_map.end()){
            disp_map.insert(make_pair(pi, make_tuple(p1.x-pi.x,p1.y-pi.y)));
        }

        for(EdgeDatum edat : iter->second){
            p2 = strained_edges[edge_list_index].p2;

            pi = Point(point_list[edat.index*2],point_list[edat.index*2+1]);
            if(disp_map.find(pi) == disp_map.end()){
                disp_map.insert(make_pair(pi, make_tuple(p2.x-pi.x,p2.y-pi.y)));
            }
            edge_list_index ++;
        }

        if(disp_map.size() == point_count) break;
    }

    return disp_map;
}

int main(int argc, char **argv){
    vector<double> point_list;
    map<int, vector<EdgeDatum>> neighbor_map;
    ifstream unstrained_file, strained_file;
    FILE *disp_report;
    string init, strained, report_name;
    vector<Edge> strained_edges;
    map<Point, tuple_2d> ref_disp_map, hole_disp_map;
    tuple_2d ref_disp, hole_disp;
    double xdiff, ydiff;

    //First make a map from points in a complete network to the displacement
    //field for each of these points in the strained state
    cout << "Enter the name of the unstrained reference network: ";
    getline(cin, init);
    cout << "Enter the name of the strained reference network: ";
    getline(cin, strained);

    unstrained_file.open(init);
    strained_file.open(strained);

    read_relax_style(unstrained_file, point_list, neighbor_map);
    strained_edges = read_edges(strained_file);

    unstrained_file.close();
    strained_file.close();

    ref_disp_map = build_displacement_map(point_list, neighbor_map, strained_edges);
    point_list.clear();
    neighbor_map.clear();
    strained_edges.clear();

    do{
        //Now prompt for unstrained and strained networks from which holes have
        //been removed
        cout << "Enter the name of the unstrained network: with a hole: ";
        getline(cin, init);
        cout << "Enter the name of the strained network with a hole: ";
        getline(cin, strained);

        unstrained_file.open(init);
        strained_file.open(strained);

        read_relax_style(unstrained_file, point_list, neighbor_map);
        strained_edges = read_edges(strained_file);
        hole_disp_map = build_displacement_map(point_list, neighbor_map, strained_edges);
        point_list.clear();
        neighbor_map.clear();
        strained_edges.clear();

        //Prompt for a file name for reporting the discrepancy between the
        //reference displacement field and the displacement field with a hole
        cout << "Enter the name for the report file: ";
        getline(cin, report_name);
        disp_report = fopen(report_name.c_str(), "w");

        //Iterate over points in the network with a hole. Substract from the
        //corresponding displacement vector the displacement vector for the
        //same point in the network with no hole.
        for(auto iter=hole_disp_map.begin();iter!=hole_disp_map.end();iter++){
            hole_disp = iter->second;
            if(ref_disp_map.find(iter->first) != ref_disp_map.end()){
                ref_disp = ref_disp_map.find(iter->first)->second;
                xdiff = get<0>(hole_disp) - get<0>(ref_disp);
                ydiff = get<1>(hole_disp) - get<1>(ref_disp);
                fprintf(disp_report, "%2.8lf %2.8lf %2.8lf %2.8lf\n", iter->first.x, iter->first.y, xdiff, ydiff);
            }
            else cerr << "Point not found.\n";
        }
        fclose(disp_report);
        hole_disp_map.clear();
    }while(yesno("Process another network?"));

    return 0;
}
