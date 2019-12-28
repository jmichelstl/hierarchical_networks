/*
 * This code reads in a network, then cuts the network along a set of lines
 * passing through the origin. The user is prompted for the minimum and
 * maximum angle the cut line should make with the x axis, the increment
 * by which this angle should be varied. The user is also asked to choose
 * between two policies: one in which edges along the cut line are simply
 * discarded, and one in which edges that straddle the cut line are kept,
 * provided each end point is connected to at least two neighbors.
 */

#include <iostream>
#include "network_utils.h"
#include <vector>
#include <unordered_set>
#include <map>
#include <string>
#include <tuple>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <unistd.h>
#include <algorithm>

using namespace std;

//Theshold for considering a floating point number to be effectively zero
#define threshold 1e-8

//Structure for hashing tuples containing a pair of integers
namespace std{
template<> struct hash<tuple<int,int>>{

    size_t operator() (const tuple<int, int> &int_pair) const;
};

size_t hash<tuple<int,int>>::operator() (const tuple<int, int> &int_pair) const{
    return hash<int>()(get<0>(int_pair))*hash<int>()(get<1>(int_pair));
}
}

/*Structure for sorting a list of neighbors in an entry of an adjacency map,
which maps the integer index of one point to the indices of that point's
neighbors.
*/
struct adj_map_comparator{

    adj_map_comparator(vector<Point> plist,int base) : points(plist),key(base){}

    vector<Point> points;
    int key;

    bool operator() (int i, int j){
        double ang1=atan2(points[i].y-points[key].y,points[i].x-points[key].x);
        double ang2=atan2(points[j].y-points[key].y,points[j].x-points[key].x);
        return ang1 < ang2;
    }

};

//Import a network, and return a list of points, and a set of neighbor indices
//indicating the points to which a point is connected.
void read_network(ifstream &input, vector<Point> &pts, map<int, vector<int>> &nlists){
    map<Point, int> pmap;
    int pindex = 0, index1, index2, point_count = 0;
    string nextline;
    vector<double> edge_data;
    Point p1, p2;

    while(!input.eof()){
        getline(input, nextline);
        edge_data = parse_doubles(split(nextline, ' '));

        if(edge_data.size() >= 2){
            if(point_count  == 0){
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
                    nlists.insert(make_pair(pindex, vector<int>()));
                    pts.push_back(p1);
                    pindex++;
                }
                if(pmap.find(p2) == pmap.end()){
                    pmap.insert(make_pair(p2, pindex));
                    nlists.insert(make_pair(pindex, vector<int>()));
                    pts.push_back(p2);
                    pindex++;
                }

                index1 = pmap[p1];
                index2 = pmap[p2];
                nlists[index1].push_back(index2);
                nlists[index2].push_back(index1);
            }

        }
    }

    input.close();
}


//Find the next edge during a traversal of a graph to find its edges
void get_next(int &i, int &j, map<int, vector<int>> adj_map){

    int pos = 0;
    while(adj_map[j][pos] != i) pos++;

    pos = (pos + 1) % adj_map[j].size();
    i = j;
    j = adj_map[j][pos];
}

//Extract the border of a network. Start by choosing the lowest point in
//the network. Then enumerate the graph faces in which that point participates.
vector<int> get_border(vector<Point> points, map<int, vector<int>> nlists){

    int iter, initial, current, next, maxSize, maxIndex;
    unordered_set<tuple<int,int>> visited;
    vector<vector<int>> faces;
    vector<int> currFace;

    //Sort each point's list of neighbors according to the angle between the 
    //positive x axis and a line segment from that point to a neighbor.
    for(iter = 0; iter < points.size(); iter ++){
        sort(nlists[iter].begin(), nlists[iter].end(),adj_map_comparator(points,iter));
    }

    //Iterate over each edge in the graph. If an edge has not been visited yet,
    //use it as the beginning of a graph face, and traverse the face until
    //the search returns to the starting point
    for(auto miter = nlists.begin(); miter != nlists.end() && !(faces.size() > 0); miter ++){
        initial = miter->first;
        for(int neighbor : nlists[initial]){
            if(visited.find(make_tuple(initial,neighbor))==visited.end()){
                current = initial;
                next = neighbor;
                visited.emplace(make_tuple(initial, next));
                currFace.clear();
                currFace.push_back(initial);
                while(next != initial){
                    get_next(current, next, nlists);
                    visited.emplace(make_tuple(current, next));
                    currFace.push_back(current);
                }
                faces.push_back(currFace);
            }
        }
    }

    maxIndex = 0;
    maxSize = 0;
    for(iter = 0; iter < faces.size(); iter++){
        if(faces[iter].size() > maxSize){
            maxSize = faces[iter].size();
            maxIndex = iter;
        }
    }

    return faces[maxIndex];
}

//Update a set of neighbor lists to indicate that two vertices share a bond
void update_neighbor_lists(map<int,vector<int>> &nlists,int index1,int index2){
    if(nlists.find(index1) == nlists.end()){
        nlists.insert(make_pair(index1, vector<int>()));
    }
    if(nlists.find(index2) == nlists.end()){
        nlists.insert(make_pair(index2, vector<int>()));
    }
    nlists[index1].push_back(index2);
    nlists[index2].push_back(index1);
}

//Delete a given integer from a list of integers, if it is present.
void remove_from_vector(vector<int> &list, int target){
    int index = 0;
    while(index < list.size() && list[index] != target) index ++;
    if(index < list.size()) list.erase(list.begin() + index);
}

void remove_singletons(map<int, vector<int>> &nlists){
    vector<int> singletons;

    //Find singletons, and remove them from their neighbors lists of neighbors
    //while traversing the map.
    for(auto miter  = nlists.begin(); miter != nlists.end(); miter ++){
        if(miter->second.size() == 1){
            remove_from_vector(nlists[miter->second[0]], miter->first);
            singletons.push_back(miter->first);
        }
    }

    //Now purge all pairs from the map from indices to neighbors
    for(int next_index : singletons){
        nlists.erase(next_index);
    }
}

//Given the outer border of a network, and the border of a subset of that
//network, find the mean connectivity of the portion of the subset that does
//not belong to the outer border.
double calc_z(unordered_set<int> outer, vector<int> sub, map<int, vector<int>> nlists){

    int count = 0;
    double mean_z = 0;
    bool currIn = false, prevIn = false;

    for(int nextIndex : sub){
        if(outer.find(nextIndex) == outer.end()){
            count ++;
            mean_z += nlists[nextIndex].size();
            currIn = true;
        }
        else if(prevIn) break;

        prevIn = currIn;
    }

    return mean_z / count;
}

void report_net_and_border(vector<Point> points, unordered_set<int> outer, map<int, vector<int>> nlists, vector<int> border, char *net_name, char *border_name, char *hist_name){

    FILE *net_file, *border_file, *hist_file;
    Point p1, p2;

    //Report the edges of the sub-network
    net_file = fopen(net_name, "w");
    for(auto miter = nlists.begin(); miter != nlists.end(); miter ++){
        p1 = points[miter->first];
        for(int nextNeighbor : miter->second){
            //if(nextNeighbor > miter->first){
                p2 = points[nextNeighbor];
                fprintf(net_file, "%12.8f\t%12.8f\n", p1.x, p1.y);
                fprintf(net_file, "%12.8f\t%12.8f\n\n", p2.x, p2.y);
            //}
        }
    }
    fclose(net_file);

    //Report the border of the sub-network
    border_file = fopen(border_name, "w");
    for(int nextPoint : border){
        fprintf(border_file, "%12.8f\t%12.8f\n", points[nextPoint].x,points[nextPoint].y);
    }

    fclose(border_file);

    //Report a histogram of connectivity vales for nodes along the cut line
    if(hist_name != NULL){
        hist_file = fopen(hist_name, "w");
        for(int nextIndex : border){
            if(outer.find(nextIndex) == outer.end()){
                fprintf(hist_file, "%d\n", nlists[nextIndex].size());
            }
        }

        fclose(hist_file);
    }
}

void report_sub_networks(vector<Point> points, unordered_set<int> outer, map<int, vector<int>> upper_nlists, map<int, vector<int>> lower_nlists, vector<int> upper_border, vector<int> lower_border){

    string base = "";
    char *net_name, *border_name, *hist_name = NULL;

    while(base.compare("") == 0){
        cout << "Enter the base name: ";
        getline(cin, base);
    }

    net_name = (char *) malloc(sizeof(char) * (base.size() + 11));
    border_name = (char *) malloc(sizeof(char) * (base.size() + 18));

    if(yesno("Report connectivity histogram?")){
        hist_name = (char *) malloc(sizeof(char) * (base.size() + 16));
    }

    //Report the upper, and then the lower network
    sprintf(net_name, "%s_upper.dat", base.c_str());
    sprintf(border_name, "%s_upper_border.dat", base.c_str());
    if(hist_name != NULL){
        sprintf(hist_name, "%s_upper_hist.txt", base.c_str());
    }
    report_net_and_border(points, outer, upper_nlists, upper_border, net_name, border_name, hist_name);

    sprintf(net_name, "%s_lower.dat", base.c_str());
    sprintf(border_name, "%s_lower_border.dat", base.c_str());
    if(hist_name != NULL){
        sprintf(hist_name, "%s_lower_hist.txt", base.c_str());
    }
    report_net_and_border(points, outer, lower_nlists, lower_border, net_name, border_name, hist_name);

    free(net_name);
    free(border_name);
    if(hist_name != NULL) free(hist_name);
}

//Sort a network into a piece in which all nodes are above a cut line or
//below a cut line. Optionally, also keep bonds straddling the cut line.
double split_network(vector<Point> points, map<int,vector<int>> nlist, unordered_set<int> outer_border, Point center, double angle, char mode, bool report){

    map<int, vector<int>> upper_nlists, lower_nlists;
    vector<int> upper_border, lower_border;
    int upper_count, lower_count;
    double sine = sin(angle), cosine = cos(angle), upper_z, lower_z;
    double y1_rot, y2_rot, upperMax, lowerMin;
    Point p1, p2;

    upperMax = FLT_MIN;
    lowerMin = FLT_MAX;

    //Divide the network into the upper and lower groups, with edges straddling
    //the cut line retained
    for(auto miter = nlist.begin(); miter != nlist.end(); miter++){
        p1 = points[miter->first];
        y1_rot = (p1.y-center.y)*cosine+(p1.x-center.x)*sine;
        if(abs(y1_rot) < threshold) y1_rot = 0;

        for(int nextNeighbor : miter->second){
            if(nextNeighbor > miter->first) continue;

            p2 = points[nextNeighbor];
            y2_rot = (p2.y-center.y)*cosine+(p2.x-center.x)*sine;
            if(abs(y2_rot) < threshold) y2_rot = 0;

            if(y1_rot > 0 && y2_rot > 0){
                update_neighbor_lists(upper_nlists, miter->first,nextNeighbor);
            }
            else if(y1_rot < 0 && y2_rot < 0){
                update_neighbor_lists(lower_nlists, miter->first,nextNeighbor);
            }
            else if((y1_rot==0&&y2_rot>0)||(y1_rot>0&&y2_rot==0)){
                update_neighbor_lists(upper_nlists, miter->first,nextNeighbor);
            }
            else if((y1_rot==0&&y2_rot<0)||(y1_rot<0&&y2_rot==0)){
                update_neighbor_lists(lower_nlists, miter->first,nextNeighbor);
            }
            else if(y1_rot == 0 && y2_rot == 0){
                update_neighbor_lists(upper_nlists, miter->first,nextNeighbor);
                update_neighbor_lists(lower_nlists, miter->first,nextNeighbor);
            }
            else if(mode == 'k'){
                update_neighbor_lists(upper_nlists, miter->first,nextNeighbor);
                update_neighbor_lists(lower_nlists, miter->first,nextNeighbor);
            }
        }
    }

    //Remove nodes with just one neighbor.
    remove_singletons(upper_nlists);
    remove_singletons(lower_nlists);

    //Find the borders of the upper and lower networks, and determine which
    //vertices are in borders for the sub-networks, but not the overall border.
    upper_border = get_border(points, upper_nlists);
    upper_z = calc_z(outer_border, upper_border, upper_nlists);

    lower_border = get_border(points, lower_nlists);
    lower_z = calc_z(outer_border, lower_border, lower_nlists);

    //If the report option was selected, offer the opportunity to report the
    //upper and lower networks to a file.
    if(report && yesno("Report subnetworks?")){
        report_sub_networks(points, outer_border, upper_nlists, lower_nlists, upper_border, lower_border);
    }

    return (upper_z + lower_z) / 2;
}

int main(int argc, char **argv){

    double angle, minAngle, maxAngle, incr, meanx = 0, meany = 0;
    vector<double> parse_results;
    vector<Point> points;
    unordered_set<int> outer_border, subBorder1, subBorder2;
    map<int, vector<int>> nlists;
    ifstream input;
    ofstream report;
    string response;
    char choice, option;
    int iter;
    vector<double> results;
    bool make_report = false;

    //Check for a command line switch
    while((option = getopt(argc, argv, "r")) != -1){
        switch(option) {
            case 'r':
                make_report = true;
                break;
            case '?':
                if(isprint(optopt)){
                    fprintf(stderr, "Unknown option: -%c.\n", optopt);
                }
                else{
                    fprintf(stderr, "Unknown option character.\n");
                }
            default:
                break;
         }
    }

    //Prompt for the type of cut to perform
    while(true){
        cout << "Enter \"o\" to obliterate or \"k\" to keep: ";
        getline(cin, response);

        if(response.length() == 0) continue;

        switch((choice = response[0])){
            case 'o':
                break;
            case 'k':
                break;
            default:
                cerr << "Enter either \"o\" or \"k\".\n";
        }

        if(choice == 'o' || choice == 'k') break;
    }

    //Prompt for the minimum and maximum angles, and the increment in the
    //angle.
    while(true){
        parse_results = getdoubles("Enter the minimum and maximum angles, and increment: ");
        if(parse_results.size() < 3){
            cerr << "Enter three numbers.\n";
            continue;
        }
        if(parse_results[0] < parse_results[1] && parse_results[2] > 0){
            minAngle = parse_results[0];
            maxAngle = parse_results[1];
            incr = parse_results[2];
            break;
        }

        else{
            cerr << "The minimum angle must be less than the maximum angle, ";
            cerr << "and the increment must be positive.\n";
            continue;
        }
    }

    //Prompt for the input network, and extract the nodes and edges
    open_dat_file("Enter the input network: ", input);
    if(! input.is_open()){
        cerr << "No input file could be read.\n";
        exit(-1);
    }
    
    read_network(input, points, nlists);

    //Find the mean position of all nodes in the network
    for(Point p : points){
        meanx += p.x;
        meany += p.y;
    }
    meanx /= points.size();
    meany /= points.size();

    cout << "Center: " << meanx << ", " << meany << "\n";

    //Find the border of the outer network
    vector<int> bvec = get_border(points, nlists);
    outer_border = unordered_set<int>(bvec.begin(), bvec.end());

    //Step through each cut line angle, and find the mean connectivity of
    //nodes along the cut line.
    for(angle = minAngle; angle <= maxAngle; angle += incr){
        results.push_back(angle);
        results.push_back(split_network(points, nlists, outer_border, Point(meanx, meany), angle, choice, make_report));
    }

    //Report the angle and mean connectivity pairs to a file
    open_output_file("Enter the name of the report file: ", report);
    if(report.is_open()){
        for(iter = 0; iter < results.size() - 1; iter += 2){
            report << results[iter] << "\t" << results[iter + 1] << "\n";
        }
        report.close();
    }

    return 0;
}
