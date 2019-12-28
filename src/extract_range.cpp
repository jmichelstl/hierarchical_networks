#include <iostream>
#include <fstream>
#include <vector>
#include "network_utils.h"

/*
Given a network specified as a set of x-y coordinates, extract those bonds
within a bounding box specified by the user.
*/

bool in_bounds(Point p, vector<double> bounds){
    return p.x >= bounds[0] && p.x <= bounds[2] && p.y >= bounds[1] && p.y <= bounds[3];
}

vector<Edge> get_edges_in_bounds(ifstream& datfile, vector<double> bounds){

    string nextline;
    Point p1, p2;
    int point_count;
    vector<double> edge_data;
    vector<Edge> toKeep;

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
                if(in_bounds(p1, bounds) || in_bounds(p2, bounds)){
                    p1.x -= bounds[0];
                    p1.y -= bounds[1];
                    p2.x -= bounds[0];
                    p2.y -= bounds[1];
                    toKeep.push_back(Edge(p1,p2));
                }
            }
        }
    }

    return toKeep;
}

int main(int argc, char **argv){
    vector<double> bounds;
    vector<Edge> toKeep;
    ifstream datfile;
    ofstream croppedFile;
    bool valid = false, tryAgain = true;

    //Obtain the boundaries in which points should lie
    do{
        bounds = getdoubles("Enter minimum and maximum x and y coordinates: ");
        if(bounds.size() != 4){
            cerr << "Four numbers are needed.";
            tryAgain = yesno("Try again?");
        }
        else valid = true;
    }while(!valid && tryAgain);

    if(! valid){
        cerr << "Bounds were not specified. The program will now exit.\n";
        return -1;
    }

    //Obtain the data file from which to read edges
    open_dat_file("Enter the network file: ", datfile);
    if(! datfile.is_open()){
        cerr << "No network file was specified. The program will now exit.\n";
        return -1;
    } 

    toKeep = get_edges_in_bounds(datfile, bounds);
    datfile.close();

    open_output_file("Enter the name for the cropped network: ", croppedFile);
    if(! croppedFile.is_open()){
        cerr << "No file was opened for writing. The program will now exit.\n";
        return -1;
    }

    for(Edge nextEdge : toKeep){
        croppedFile << nextEdge.p1.x << "\t" << nextEdge.p1.y << "\n";
        croppedFile << nextEdge.p2.x << "\t" << nextEdge.p2.y << "\n\n";
    }

    croppedFile.close();

    return 0;
}
