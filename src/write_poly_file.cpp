/*
This program reads data describing a planar line graph by specifying pairs
of end points for edges, and creates a .poly file. The .poly file lists all
points in the planar line graph, indexed by natural numbers starting from 1,
and specifies a set of polygons whose vertices are indicated by integer
indices. There will be one bounding polygon, and a set of interior polygons,
which constitute "holes". A point in each hole is located, so that the planar
line graph can be used to produce a triangulation suitable for finite element
analysis.
*/

#include <iostream>
#include <unordered_map>
#include "network_utils.h"
#include <vector>
#include <cmath>
#include <list>

using namespace std;

/*
If a point is not currently in an index map, add it and map it to the integer
index one greater than the current size of the map.
*/
void update_maps(Point p, map<Point, int> &pmap, map<int, list<int>> &emap){
    int index;

    if(pmap.find(p) == pmap.end()){
        index = pmap.size() + 1;
        pmap.insert(make_pair(p, index));
        emap.insert(make_pair(index, list<int>()));
    }
}

/*
This function processes a file describing a planar line graph by specifying
pairs of endpoints of line segments.

Arguments:
datfile - An open stream to a file describing the planar line graph
imap - A map from vertices in the planar line graph to integer indices
emap - A map in which each pair is a mapping from a point's integer index to
the indices of all other points to which that point is connected
*/
void read_edges(ifstream &datfile, map<int, Point> &imap, map<int, list<int>> &emap){

    string nextline;
    Point p1, p2;
    int point_count, index1, index2;
    vector<double> edge_data;
    map<Point, int> pmap;

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

                //Ensure that both points in the edge just read are accounted
                //for in both maps
                update_maps(p1, pmap, emap);
                update_maps(p2, pmap, emap);

                //Indicate that the two points are now known to be joined by
                //an edge
                index1 = pmap[p1];
                index2 = pmap[p2];
                emap[index1].push_back(index2);
                emap[index2].push_back(index1);
            }
        }
    }

    //Reverse the mapping of point to integers to one of integers to points
    for(auto iter = pmap.begin(); iter != pmap.end(); iter++){
        imap.insert(make_pair(iter->second, iter->first));
    }

    datfile.close();
}

vector<vector<int>> build_cycles(map<int, list<int>> emap){
    int curr_index, next_index;
    map<int, list<int>>::iterator emap_pos;
    vector<int> curr_cycle;
    vector<vector<int>> cycle_list;
    bool cycle_complete;

    //Build cycles by iterating over entries in the map from indices to
    //lists of neighbor indices
    while(! emap.empty()){
        //Start a new cycle by finding the first element of the map and
        //making the key the first element in the cycle
        vector<int> curr_cycle;
        curr_index = emap.begin()->first;
        cycle_complete = false;

        //Add successive elements to the current cycle by finding the first
        //neighbor listed for the current vertex in the cycle
        do{
            next_index = emap[curr_index].front();
            curr_cycle.push_back(next_index);

            //Remove the current index and the next index from each other's
            //lists of neighbors
            emap[curr_index].remove(next_index);
            if(emap[curr_index].empty()){
                emap_pos = emap.find(curr_index);
                emap.erase(emap_pos);
            }

            emap[next_index].remove(curr_index);
            if(emap[next_index].empty()){
                emap_pos = emap.find(next_index);
                emap.erase(emap_pos);
                //When the list to which next_index is mapped has no elements,
                //the current cycle is complete
                cycle_complete = true;
            }
            curr_index = next_index;
        }while(! cycle_complete);
        cycle_list.push_back(curr_cycle);
    }

    return cycle_list;
}

//Given a set of integer indices mapped to the vertices of a polygon, find a
//point known to be inside that polygon
Point get_interior_point(vector<int> cycle, map<int, Point> imap){
    int curr_index, next_index, length = cycle.size();
    Point prev_p, curr_p, next_p, trial1, trial2;
    double dx1, dx2, dy1, dy2, angle1, angle2, ang_diff;
    vector<Edge> poly_edges;

    //Map the indices back to points and construct a polygon
    curr_p = imap[cycle[0]];
    for(curr_index = 0; curr_index != length; curr_index ++){
        next_index = (curr_index + 1) % length;
        next_p = imap[cycle[next_index]];
        poly_edges.push_back(Edge(curr_p, next_p));
        curr_p = next_p;
    }

    PolyData pdatum = PolyData(poly_edges);
    prev_p = poly_edges[0].p1;
    curr_p = poly_edges[0].p2;
    next_p = poly_edges[1].p2;

    dx1 = prev_p.x - curr_p.x;
    dy1 = prev_p.y - curr_p.y;
    dx2 = next_p.x - curr_p.x;
    dy2 = next_p.y - curr_p.y;

    //Find points at the ends of two antiparallel line segments that bisect
    //the angle between two consecutive line segments in the polygon. If
    //neither point is in the polygon, halve the length of the segments until a
    //match is found
    angle1 = atan2(dy1, dx1);
    angle2 = atan2(dy2, dx2);
    ang_diff = (angle2 - angle1) / 2.;
    dx2 = (dx1 * cos(ang_diff) - dy1 * sin(ang_diff)) / 2.;
    dy2 = (dy1 * cos(ang_diff) + dx1 * sin(ang_diff)) / 2.;

    do{
        trial1 = Point(curr_p.x + dx2, curr_p.y + dy2);
        trial2 = Point(curr_p.x - dx2, curr_p.y - dy2);

        if(pdatum.contains(trial1)) return trial1;
        if(pdatum.contains(trial2)) return trial2;

        dx2 /= 2.;
        dy2 /= 2.;
    }while(dx2 >= FLOAT_TOL);

    cout << "Failed to find an interior point.\n";

    return trial1;
}

/*
*Locate a starting point in a cycle that is not part of a colinear group of
*three or more points
*/
int pick_start(vector<int> cycle, map<int, Point> imap){
    int start_index, curr_index, next_index_1, next_index_2;
    Point curr_p, next_p_1, next_p_2;
    double slope1, slope2;
    size_t length = cycle.size();

    for(curr_index = 0; curr_index < length; curr_index ++){
        next_index_1 = (curr_index + 1) % length;
        next_index_2 = (curr_index + 2) % length;
        curr_p = imap[cycle[curr_index]];
        next_p_1 = imap[cycle[next_index_1]];
        next_p_2 = imap[cycle[next_index_2]];

        //If the slopes of the next two edges are different, choose the
        //next index as the starting index
        slope1 = (next_p_1.y - curr_p.y) / (next_p_1.x - curr_p.x);
        slope2 = (next_p_2.y - next_p_1.y) / (next_p_2.x - next_p_1.x);

        if(abs(slope2 - slope1) > FLOAT_TOL){
            start_index = next_index_1;
            break;
        }
    }

    return start_index;
}

/*
*Produce a more compact list of points for the planar line graph that eliminates
groups of three or more colinear points, and generate the corresponding edges
and points inside interior polygons that indicate holes.
*/
void get_poly_file_info(map<int, Point> imap, vector<vector<int>> cycle_list, map<int, Point> &reduced_map, vector<int> &edges, vector<Point> &holes){
    int start_index, curr_index, next_index, cycle_index = 0, begindex;
    int min_cycle = 0;
    Point prev_point, curr_point, next_point;
    double slope1, slope2, dx1, dy1, dx2, dy2, angle1, angle2;
    double min = FLT_MAX, cycle_min;
    size_t length;

    for(vector<int> next_cycle : cycle_list){
        length = next_cycle.size();

        //Locate a starting index at a corner of a polygon, and not along
        //an edge
        start_index = pick_start(next_cycle, imap);

        //Iterate over indices in the current cycle, beginning with the
        //strategically chosen starting point, and coallesce edges, then add
        //end points to the minimal map and update edge list
        curr_index = start_index + 1;
        prev_point = imap[next_cycle[start_index]];
        reduced_map.insert(make_pair(reduced_map.size() + 1, prev_point));
        begindex = reduced_map.size();
        cycle_min = prev_point.y;

        do{
            next_index = (curr_index + 1) % length;
            curr_point = imap[next_cycle[curr_index]];
            next_point = imap[next_cycle[next_index]];
            if(curr_point.y < cycle_min) cycle_min = curr_point.y;

            slope1 = (curr_point.y-prev_point.y)/(curr_point.x-prev_point.x);
            slope2 = (next_point.y-curr_point.y)/(next_point.x-curr_point.x);
            if(abs(slope2 - slope1) > FLOAT_TOL){
                reduced_map.insert(make_pair(reduced_map.size()+1, curr_point));
                edges.push_back(reduced_map.size() - 1);
                edges.push_back(reduced_map.size());
            }

            prev_point = curr_point;
            curr_index = next_index;
        }while(curr_index != start_index);
        edges.push_back(reduced_map.size());
        edges.push_back(begindex);

        //Determine whether this cycle has the lowest y coordinate of any cycle
        //so far
        if(cycle_min < min){
            min = cycle_min;
            min_cycle = cycle_index;
        }

        cycle_index ++;
    }

    //The cycle with the lowest y coordinate is the border. Remove this from
    //the list of the cycles, then find a point within every other polygon in
    //order to indicate all holes in the planar line graph.
    cycle_list.erase(cycle_list.begin() + min_cycle, cycle_list.begin() + min_cycle + 1);

    for(vector<int> next_cycle : cycle_list){
        holes.push_back(get_interior_point(next_cycle, imap));
    }

}

void write_poly_file(ofstream &poly_file, map<int, Point> point_map, vector<int> edges, vector<Point> holes){
    int index;
    double x,y;
    Point p;

    //Write initial header
    poly_file << point_map.size() << " " << "2" <<" "<< "0" << " " << "0" << "\n";
    //Write table of points
    for(auto iter = point_map.begin(); iter != point_map.end(); iter++){
        index = iter->first;
        x = iter->second.x;
        y = iter->second.y;
        poly_file << index << " " << x << " " << y << "\n";
    }

    //Write description of segments followed by table of segments
    poly_file << edges.size() / 2 << " " << "0" << "\n";
    for(index = 0; index < edges.size(); index += 2){
        poly_file << index / 2 + 1 << " " << edges[index] << " ";
        poly_file << edges[index + 1] << "\n";
    }

    //Write hole count, followed by table of holes
    poly_file << holes.size() << "\n";
    for(index = 0; index < holes.size(); index++){
       p = holes[index];
       poly_file << index + 1 << " " << p.x << " " << p.y << "\n";
    }

    poly_file.close();
}

/*
*Prompt the user for the name of a data file and an output .poly file. Create
an open input stream, obtain a map of points to indices, and a map of indices
to lists of other indices. Use the two maps to build a table of cycles, then
convert the table of cycles to polygon data.
*/
int main(int argc, char **argv){
    ifstream datfile;
    ofstream poly_file;
    map<int, Point> imap;
    map<int, list<int>> emap;
    vector<vector<int>> cycle_list;
    map<int, Point> reduced_map;
    vector<int> edges;
    vector<Point> holes;

    //Obtain network data and obtain maps from points to integers and from
    //integers to lists of neighbor indices
    open_dat_file("Enter the data file describing the network: ", datfile);
    if(! datfile.is_open()) exit(1);

    read_edges(datfile, imap, emap);

    //Use maps to build cycles described by integer indices in cyclical order
    cycle_list = build_cycles(emap);

    //Get final data to construct a .poly file
    get_poly_file_info(imap, cycle_list, reduced_map, edges, holes); 

    open_output_file("Enter the name for the .poly file: ", poly_file);
    if(! poly_file.is_open()) exit(1);
    write_poly_file(poly_file, reduced_map, edges, holes);

    return 0;
}
