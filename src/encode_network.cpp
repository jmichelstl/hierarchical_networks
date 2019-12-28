#include "network_utils.h"
#include <iostream>
#include <fstream>
#include <cmath>

#define dx .5
#define dy 0.866025404
#define tol 1e-7
#define edge1 1
#define edge2 2
#define edge3 4

using namespace std;

/*Encode an edge as type 1, two or three, corresponding to the left and up,
  right and up, and right directions, respectively
*/
unsigned char encoding(Point p1, Point p2, double length){
    double dx_scaled = length * dx;

    if(abs(p2.x - p1.x + dx_scaled) <= tol) return (unsigned char) edge1;
    else if(abs(p2.x - p1.x - dx_scaled) <= tol) return (unsigned char) edge2;
    else return (unsigned char) edge3;
}

/* Read a planar line graph containing a dilute, triangular lattice, possibly
   with a hierarchical structure. Return a map from points to an encoding of
   edges beginning at those points.
*/
void read_edges(ifstream& datfile, map<Point, unsigned char>& code_map, double length, double &minx, double &maxx, double &miny, double &maxy){

    map<Point, vector<Point>> neighbor_map;
    string nextline;
    Point p1, p2;
    int point_count = 0;
    vector<double> edge_data;
    unsigned char edge_code; 

    minx = FLT_MAX;
    maxx = FLT_MIN;
    miny = FLT_MAX;
    maxy = FLT_MIN;

    //Produce a map from points to neighbors connected to those points
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

                if(neighbor_map.find(p1) == neighbor_map.end()){
                    neighbor_map.insert(make_pair(p1, vector<Point>()));
                    minx = p1.x < minx - FLOAT_TOL ? p1.x : minx;
                    maxx = p1.x > maxx + FLOAT_TOL ? p1.x : maxx;
                    miny = p1.y < miny - FLOAT_TOL ? p1.y : miny;
                    maxy = p1.y > maxy + FLOAT_TOL ? p1.y : maxy;
                }
                if(neighbor_map.find(p2) == neighbor_map.end()){
                    neighbor_map.insert(make_pair(p2, vector<Point>()));
                    minx = p2.x < minx - FLOAT_TOL ? p2.x : minx;
                    maxx = p2.x > maxx + FLOAT_TOL ? p2.x : maxx;
                    miny = p2.y < miny - FLOAT_TOL ? p2.y : miny;
                    maxy = p2.y > maxy + FLOAT_TOL ? p2.y : maxy;
                }
            }

            if(p1 < p2) neighbor_map[p1].push_back(p2);
            else neighbor_map[p2].push_back(p1);
        }
    }

    //Now convert neighbor map to a map from points to a coded list of edges
    //beginning at each point
    for(auto iter = neighbor_map.begin(); iter != neighbor_map.end(); iter++){
        edge_code = 0;
        for(Point neighbor : iter->second){
            edge_code |= encoding(iter->first, neighbor, length);
        }
        code_map.insert(make_pair(iter->first, edge_code));
    }
}

//Pad a stream with a specified number of zeros
void write_padding(FILE *stream, int count){
    int iter;

    for(iter = 0; iter < count; iter ++){
        fputc(0, stream);
    }
}

/*Prompt the user for a file containing a network specified in the format
  p1x    p1y
  p2x    p2y
    ...
  p(2n-1)x    p(2n-1)y
  p2nx    p2ny

  where n is the number of edges in the network. Use the helper function
  read_edges to convert this to a map from points to encodings of edges
  beginning at those points. Write these encodings in a grid, with filler for
  places where no points are present.
*/
int main(int argc, char **argv){
    double length = 0, miny, maxy, minx, maxx;
    double x_diff, curr_x, new_x, curr_y, new_y, remainder;
    int nrows, ncolumns, padding, leftInRow, rowIter = 0;
    ifstream datfile;
    FILE *coded_file = NULL;
    map<Point, unsigned char> code_map;
    string response;
    unsigned char offset = 0;
    int num_read;

    //Prompt for the initial file describing the network and encode bonds
    open_dat_file("Enter file describing network: ", datfile);
    do{
        cout << "Enter the bond length: ";
        getline(cin, response);
        num_read = sscanf(response.c_str(), "%lf", &length);
        if(num_read == 0 || length < 0) cerr << "Enter a positive number.\n";
    }while(! length > 0);
    read_edges(datfile, code_map, length, minx, maxx, miny, maxy);
    datfile.close();

    //Open the output file for writing
    do{
        cout << "Enter a file name for the output file: ";
        getline(cin, response);
        if(! response.empty()) coded_file = fopen(response.c_str(), "wb");
    }while(coded_file == NULL);

    //Find the number of cells per row needed to describe the lattice, and the
    //number of rows in the lattice
    ncolumns = (int) round((maxx - minx) / length) + 1;
    nrows = (int) round((maxy - miny) / (dy*length)) + 1;

    //Find the left padding needed for the first row, and determine whether a
    //horizontal offset is needed
    
    x_diff = (code_map.begin()->first).x - minx;
    if(x_diff > tol){
        padding = (int) floor(x_diff / length);
        if(fmod(x_diff, length) > tol) offset = 1;
    }

    //Report bond length, rows, columns and horizontal offset of the first row
    fwrite((void *) &length, sizeof(double), 1, coded_file);
    fwrite((void *) &nrows, sizeof(int), 1, coded_file);
    fwrite((void *) &ncolumns, sizeof(int), 1, coded_file);
    fwrite((void *) &offset, sizeof(char), 1, coded_file);

    //Add padding for the first row
    if(padding > 0) write_padding(coded_file, padding);

    leftInRow = ncolumns - padding;

    //Initialize current values of x and y to prepare to step along the lattice
    curr_x = (code_map.begin()->first).x - length;
    curr_y = (code_map.begin()->first).y;

    //Now write edge codes row by row, adding padding as needed
    for(auto mapiter = code_map.begin(); mapiter != code_map.end(); mapiter++){
        new_x = mapiter->first.x;
        new_y = mapiter->first.y;

        //If the new y value is higher than the old one, add any right padding
        //needed to complete the old row, then begin the new row, starting
        //with any left padding that may be needed
        if(new_y > curr_y + tol){
            if(leftInRow > 0) write_padding(coded_file, leftInRow);

            x_diff = new_x - minx;
            padding = (int) floor(x_diff / length);
            write_padding(coded_file, padding);
            leftInRow = ncolumns - padding;
            rowIter ++;
        }

        else{
            //If two successive x values in the same row differ by more than one
            //lattice spacing, add padding to fill the gap
            x_diff = new_x - curr_x;
            if(x_diff - length > tol){
                padding = (int) round(x_diff / length) - 1;
                write_padding(coded_file, padding);
                leftInRow -= padding;
            }
        }

        //Now add the edge code for the new point and update current x and
        //y values
        fputc(mapiter->second, coded_file);
        leftInRow --;
        if(leftInRow < 0) cerr << "Too many codes written in row " << rowIter << ".\n";

        curr_x = new_x;
        curr_y = new_y;
    }

    //Finally, add any needed final right padding
    if(leftInRow > 0) write_padding(coded_file, leftInRow);

    fclose(coded_file);

    return 0;
}
