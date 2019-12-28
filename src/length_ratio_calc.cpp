/*
This program takes as arguments the names of a large-scale "skeleton" and a
description of the small-scale structure of a hierarchical network. The program
then sums the lengths of all bonds on each scale, and gives the ratio of the
small-scale to the large-scale sum. The mean and standard deviation of 
small-scale bond length is determined, and bonds with lengths one standard
deviation or more greater than the mean are discarded.
*/

#include <fstream>
#include <vector>
#include <string>
#include "network_utils.h"
#include <cmath>

using namespace std;

//Sum a list of double-precision floating point numbers
double sum(vector<double> data){
    double sum = 0;

    for(double val : data){
        sum += val;
    }

    return sum;
}

//Find the mean of an array of double precision floating-point numbers
double mean(vector<double> data){
    return sum(data) / data.size();
}

//Given a list of double-precision numbers and their mean, find the standard
//deviation of the list
double sdev(vector<double> data, double mean){
    double var = 0;
    size_t size = data.size();

    for(double val : data){
        var += val*val;
    }

    var = var / size - mean*mean;
    return sqrt(var);
}

//Find the length of every edge in a list
vector<double> get_lengths(vector<Edge> edgeList){
    vector<double> length_list;
    double length, dx, dy;

    for(Edge e : edgeList){
        dx = e.p2.x - e.p1.x;
        dy = e.p2.y - e.p1.y;
        length = sqrt(dx*dx + dy*dy);
        length_list.push_back(length);
    }

    return length_list;
}

//Given an array of double-precision numbers, keep only those values that do not
//exceed the mean by one standard deviation or more 
vector<double> discard_outliers(vector<double> data){
    vector<double> replace;
    double list_mean, list_sdev;

    list_mean = mean(data);
    list_sdev = sdev(data, list_mean);

    for(double val : data){
        if(val - list_mean < list_sdev) replace.push_back(val);
    }

    return replace;
}

Point find_match(Point p1, Point p2, double y){
    if(p1.y == y) return p1;
    else return p2;
}

//Trim skeletal network so that edges extending beyond the y bounds for
//the small-scale network are truncated at the upper and lower bounds
void trim_top_bottom(vector<Edge> &edge_list, double ymin, double ymax){
    vector<double> ycuts;
    Point p1, p2, new_point;
    vector<Edge> replace;
    double intersect_x;
    Edge next;

    while(!edge_list.empty()){
        next = edge_list.front();
        p1 = next.p1;
        p2 = next.p2;
        edge_list.erase(edge_list.begin());
        if(next.bottom >= ymin && next.top <= ymax){
            replace.push_back(next);
        }

        else if(next.bottom < ymin && next.top >= ymin){
            intersect_x = x_intersect(next, ymin);
            new_point = Point(intersect_x, ymin);
            replace.push_back(Edge(new_point, find_match(p1, p2, next.top)));
        }
        
        else if(next.bottom <= ymax && next.top > ymax){
            intersect_x = x_intersect(next, ymax);
            new_point = Point(intersect_x, ymax);
            replace.push_back(Edge(find_match(p1, p2, next.bottom), new_point));
        }
    }

    edge_list = replace;
}

//Given a file stream pointing to a list of edges defining a small-scale
//network, extract the lengths of each edge and the minimum and maximum y
//coordinates for all points in the network
void extract_lengths(ifstream& datfile, vector<double>& lengths, double &miny, double &maxy){

    string nextline;
    double x1, y1, x2, y2;
    int point_count = 0;
    vector<double> edge_data;

    miny = FLT_MAX;
    maxy = FLT_MIN;

    while(! datfile.eof()){
        getline(datfile, nextline);

        edge_data = parse_doubles(split(nextline, ' '));
        if(edge_data.size() >= 2){
            if(point_count == 0){
                x1 = edge_data[0];
                y1 = edge_data[1];
                miny = y1 < miny ? y1 : miny;
                maxy = y1 > maxy ? y1 : maxy;
                point_count ++;
            }
            else{
                x2 = edge_data[0];
                y2 = edge_data[1];
                miny = y2 < miny ? y2 : miny;
                maxy = y2 > maxy ? y2 : maxy;
                point_count ++;
            }

            if(point_count == 2){
                point_count = 0;

                lengths.push_back(sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1)));
            }
        }
    }

    datfile.close();
}

//Given a file stream pointing to a description of a network, extract a
//list of edges
vector<Edge> extract_edges(ifstream& datfile){

    string nextline;
    Point p1, p2;
    int point_count = 0;
    vector<double> edge_data;
    vector<Edge> edge_list;

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

                edge_list.push_back(Edge(p1, p2));
            }
        }
    }

    datfile.close();
    return edge_list;
}

//Given two file streams pointing to data files with lists of edges, construct
//networks described by the data files, and find the sum of all edge lengths
//in each network, then compute the ratio of one of these sums to the other
void find_sums(ifstream& skel_stream, ifstream& small_stream, double &skel_sum, double &small_sum){
    double ymin, ymax, dummy1, dummy2;
    vector<double> skel_lengths, small_lengths;
    vector<Edge> skel_edges;
    bool trim, cull;

    trim = yesno("Trim skeleton?");
    cull = yesno("Remove outliers?");

    //Extract information about small-scale and skeletal networks
    extract_lengths(small_stream, small_lengths, ymin, ymax);
    if(cull) small_lengths = discard_outliers(small_lengths);

    if(trim){
        skel_edges = extract_edges(skel_stream);
        trim_top_bottom(skel_edges, ymin, ymax);
        skel_lengths = get_lengths(skel_edges); 
    }
    else{
        extract_lengths(skel_stream, skel_lengths, dummy1, dummy2);
    }

    cout << "Sizes: " << skel_lengths.size() << "\t" << small_lengths.size() << "\n";

    skel_sum = sum(skel_lengths);
    small_sum = sum(small_lengths);
}

int main(int argc, char **argv){
    ifstream skeleton, small_scale;
    double skel_sum, small_sum, ratio;

    do{
        //Prompt for files describing the large-scale "skeleton" and the
        //small-scale structure of a network
        open_dat_file("Enter the skeleton file: ", skeleton);
        open_dat_file("Enter the small-scale file: ", small_scale);

        //Calculate the ratio of total length of small-scale bonds to
        //total length of large scale bonds if two file streams are
        //successfully opened
        if(skeleton.is_open() && small_scale.is_open()){
            find_sums(skeleton, small_scale, skel_sum, small_sum);
            ratio = small_sum / skel_sum;
            cout << "Large: " << skel_sum << "\tSmall: " << small_sum;
            cout << "\tRatio: " << ratio << "\n";
        }

        //If there was not a skeleton and a small-scale file, close any
        //open file stream and print an error message
        else{
            if(skeleton.is_open()) skeleton.close();
            if(small_scale.is_open()) small_scale.close();
            cerr << "There were not enough open data files.\n";
        }
    }while(yesno("Process another network?"));
}
