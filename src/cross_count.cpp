#include <iostream>
#include "network_utils.h"
#include <vector>

using namespace std;

vector<Edge> read_edges(ifstream &edge_file){
    Point p1, p2;
    int point_count;
    vector<double> edge_data;
    string nextline;
    vector<Edge> edge_list;

    while(!edge_file.eof()){
        getline(edge_file, nextline);
        edge_data = parse_doubles(split(nextline, ' '));
        if(edge_data.size() >= 2){
            if(point_count == 0){
                point_count ++;
                p1 = Point(edge_data[0], edge_data[1]);
            }
            else{
                point_count ++;
                p2 = Point(edge_data[0], edge_data[1]);
            }

            if(point_count == 2){
                edge_list.push_back(Edge(p1, p2));
            }
        }
    }

    edge_file.close();
    return edge_list;
}

int main(int argc, char **argv){
    vector<Edge> edge_list;
    ifstream edge_file;
    int crossings = 0, iter1, iter2;
    size_t len;

    open_dat_file("Enter the dat file describing the nextwork: ", edge_file);
    if(! edge_file.is_open()) return 1;

    edge_list = read_edges(edge_file);
    len = edge_list.size();

    for(iter1 = 0; iter1 < len; iter1 ++){
        for(iter2 = iter1 + 1; iter2 < len; iter2 ++){
            if(intersection(edge_list[iter1], edge_list[iter2])){
                crossings ++;
            }
        }
    }

    cout << "There were " << crossings << " crossings.\n";

    return 0;
}
