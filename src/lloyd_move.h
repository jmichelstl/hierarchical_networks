/*
This code generates a random point set, then evolves it using Lloyd moves to
reduce the sum of the integrals of the squared distance from each point in the
point set over that point's corresponding Voronoi tile. Iteration stops when
the cost function changes by less than some amount from one iteration to the
next, or a maximum number of iterations have been performed.
*/

#define REAL double
#define VOID int

#include <iostream>
#include <stdio.h>
#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <string>
#include "triangle.hpp"
#include <unistd.h>
#include <map>
#include <cmath>
#include <algorithm>
#include <unordered_set>
#include <memory>
#include <queue>
#include <thread>
#include <unistd.h>
#include <mutex>
#include <time.h>
#include "VoronoiDiagramGenerator.h"


#define min(a,b) a < b ? a : b

#define tri_list_size 513

/*This array defines a lookup table that describes how a convex polygon may
 *  * be broken up into a simple triangulation. There are entries in the table
 *   * for polygons with up to 20 sides.
 *    */
int triangle_lists[513] =
{0, 1, 2, 0, 1, 2, 2, 3, 0, 0, 1, 2, 2, 3, 4, 0, 2, 4, 0, 1, 2, 2, 3,
4, 4, 5, 0, 0, 2, 4, 0, 1, 2, 2, 3, 4, 4, 5, 6, 0, 2, 4, 4, 6, 0, 0,
1, 2, 2, 3, 4, 4, 5, 6, 6, 7, 0, 0, 2, 4, 4, 6, 0, 0, 1, 2, 2, 3, 4,
4, 5, 6, 6, 7, 8, 0, 2, 4, 4, 6, 8, 0, 4, 8, 0, 1, 2, 2, 3, 4, 4, 5,
6, 6, 7, 8, 8, 9, 0, 0, 2, 4, 4, 6, 8, 0, 4, 8, 0, 1, 2, 2, 3, 4, 4,
5, 6, 6, 7, 8, 8, 9, 10, 0, 2, 4, 4, 6, 8, 8, 10, 0, 0, 4, 8, 0, 1,
2, 2, 3, 4, 4, 5, 6, 6, 7, 8, 8, 9, 10, 10, 11, 0, 0, 2, 4, 4, 6, 8,
8, 10, 0, 0, 4, 8, 0, 1, 2, 2, 3, 4, 4, 5, 6, 6, 7, 8, 8, 9, 10, 10,
11, 12, 0, 2, 4, 4, 6, 8, 8, 10, 12, 0, 4, 8, 8, 12, 0, 0, 1, 2, 2,
3, 4, 4, 5, 6, 6, 7, 8, 8, 9, 10, 10, 11, 12, 12, 13, 0, 0, 2, 4, 4,
6, 8, 8, 10, 12, 0, 4, 8, 8, 12, 0, 0, 1, 2, 2, 3, 4, 4, 5, 6, 6, 7,
8, 8, 9, 10, 10, 11, 12, 12, 13, 14, 0, 2, 4, 4, 6, 8, 8, 10, 12, 12,
14, 0, 0, 4, 8, 8, 12, 0, 0, 1, 2, 2, 3, 4, 4, 5, 6, 6, 7, 8, 8, 9,
10, 10, 11, 12, 12, 13, 14, 14, 15, 0, 0, 2, 4, 4, 6, 8, 8, 10, 12,
12, 14, 0, 0, 4, 8, 8, 12, 0, 0, 1, 2, 2, 3, 4, 4, 5, 6, 6, 7, 8, 8,
9, 10, 10, 11, 12, 12, 13, 14, 14, 15, 16, 0, 2, 4, 4, 6, 8, 8, 10,
12, 12, 14, 16, 0, 4, 8, 8, 12, 16, 0, 8, 16, 0, 1, 2, 2, 3, 4, 4, 5,
6, 6, 7, 8, 8, 9, 10, 10, 11, 12, 12, 13, 14, 14, 15, 16, 16, 17, 0,
0, 2, 4, 4, 6, 8, 8, 10, 12, 12, 14, 16, 0, 4, 8, 8, 12, 16, 0, 8,
16, 0, 1, 2, 2, 3, 4, 4, 5, 6, 6, 7, 8, 8, 9, 10, 10, 11, 12, 12, 13,
14, 14, 15, 16, 16, 17, 18, 0, 2, 4, 4, 6, 8, 8, 10, 12, 12, 14, 16,
16, 18, 0, 0, 4, 8, 8, 12, 16, 0, 8, 16, 0, 1, 2, 2, 3, 4, 4, 5, 6,
6, 7, 8, 8, 9, 10, 10, 11, 12, 12, 13, 14, 14, 15, 16, 16, 17, 18,
18, 19, 0, 0, 2, 4, 4, 6, 8, 8, 10, 12, 12, 14, 16, 16, 18, 0, 0, 4,
8, 8, 12, 16, 0, 8, 16};

//Mutex for coordinating sharing of access to a map
mutex map_mutex;

/*Structure for sorting a list of neighbors in an entry of an adjacency map,
which maps the integer index of one point to the indices of that point's
neighbors.
*/
struct adj_map_comparator{

    adj_map_comparator(shared_ptr<vector<Point>> plist,int base) : points(plist),key(base){}

    shared_ptr<vector<Point>> points;
    int key;

    bool operator() (int i, int j){
        double ang1=atan2((*points)[i].y-(*points)[key].y,(*points)[i].x-(*points)[key].x);
        double ang2=atan2((*points)[j].y-(*points)[key].y,(*points)[j].x-(*points)[key].x);
        return ang1 < ang2;
    }

};

//Structure for hashing tuples containing a pair of integers
namespace std{
template<> struct hash<tuple<int,int>>{

    size_t operator() (const tuple<int, int> &int_pair) const;
};

size_t hash<tuple<int,int>>::operator() (const tuple<int, int> &int_pair) const{
    return hash<int>()(get<0>(int_pair))*hash<int>()(get<1>(int_pair));
}
}

/*Generate a random set of points within a prescribed bounding box*/
void create_rps(vector<double> bounds, int count, double *xcoords, double *ycoords){

    int iter;
    double xrange, yrange, randx, randy;
    const gsl_rng_type *T  = gsl_rng_mt19937;
    gsl_rng *rng = gsl_rng_alloc(T);

    xrange = bounds[2] - bounds[0];
    yrange = bounds[3] - bounds[1];

    for(iter = 0; iter < count; iter ++){
        randx = gsl_rng_uniform(rng) * xrange + bounds[0];
        randy = gsl_rng_uniform(rng) * yrange + bounds[1];
        xcoords[iter] = randx;
        ycoords[iter] = randy;
    }

    gsl_rng_free(rng);
}

//Given a point and a pair of values indicating the direction of a ray extending
//from that point, find the nearest point at which the ray intersects a bounding
//box.
Point get_closest_intersection(vector<double> bounds, Point p, double dx, double dy){

    double d_sq, min_d_sq = FLT_MAX, xval, yval;
    Point closest;

    //Find the point of intersection with each segment of the box. If a line
    //segment from the starting point to the intersection point is antiparallel
    //to the ray direction, reject this option. Otherwise, find the distance
    //between the starting point and the intersection point, and update the
    //information for the closest intersection if a new closest point is found.

    yval = y_intersect(Edge(p, Point(p.x+dx, p.y+dy)), bounds[0]);

    if(dx*(bounds[0] - p.x) + dy*(yval - p.y) > 0){
        d_sq = (bounds[0]-p.x)*(bounds[0]-p.x) + (yval-p.y)*(yval-p.y);
        if(d_sq < min_d_sq){
            min_d_sq = d_sq;
            closest = Point(bounds[0], yval);
        }
    }

    yval = y_intersect(Edge(p, Point(p.x+dx, p.y+dy)), bounds[2]);

    if(dx*(bounds[2] - p.x) + dy*(yval - p.y) > 0){
        d_sq = (bounds[2]-p.x)*(bounds[2]-p.x) + (yval-p.y)*(yval-p.y);
        if(d_sq < min_d_sq){
            min_d_sq = d_sq;
            closest = Point(bounds[2], yval);
        }
    }

    xval = x_intersect(Edge(p, Point(p.x+dx, p.y+dy)), bounds[1]);

    if(dx*(xval - p.x) + dy*(bounds[1] - p.y) > 0){
        d_sq = (xval - p.x)*(xval - p.x) + (bounds[1] - p.y)*(bounds[1] - p.y);
        if(d_sq < min_d_sq){
            min_d_sq = d_sq;
            closest = Point(xval, bounds[1]);
        }
    }

    xval = x_intersect(Edge(p, Point(p.x+dx, p.y+dy)), bounds[3]);

    if(dx*(xval - p.x) + dy*(bounds[3] - p.y) > 0){
        d_sq = (xval - p.x)*(xval - p.x) + (bounds[3] - p.y)*(bounds[3] - p.y);
        if(d_sq < min_d_sq){
            min_d_sq = d_sq;
            closest = Point(xval, bounds[3]);
        }
    }

    return closest;
}

bool inBounds(vector<double> bounds, Point p){
    return p.x>=bounds[0]&&p.x<=bounds[2]&&p.y>=bounds[1]&&p.y<=bounds[3];
}

//Given a mapping from a set of integers to itself, find all permutation cycles
vector<vector<int>> get_cycles(map<int, int> index_map){

    vector<vector<int>> cycles;
    unordered_set<int> visited;
    map<int, int>::iterator map_iter = index_map.begin();
    int initial, current, cycle_index = 0;

    while(map_iter != index_map.end() && visited.size() < index_map.size()){
        if(visited.find(map_iter->first) == visited.end()){
            initial = map_iter->first;
            current = initial;
            cycles.push_back(vector<int>());
            do{
                cycles[cycle_index].push_back(current);
                visited.emplace(current);
                current = index_map[current];
            }while(current != initial);
            cycle_index ++;
        }
        map_iter ++;
    }

    return cycles;
}

/*Given a mapping from a set of integers to itself, replace each integer in a
list with its image under the mapping*/
vector<int> reassign_elems(vector<int> original, map<int, int> replace_map){
    vector<int> new_list;

    for(int next_elem : original){
        new_list.push_back(replace_map[next_elem]);
    }

    return new_list;
}

/*Given an adjacency map for a graph, sort the vertices according to a breadth-
first traversal of the graph*/
void bfs_sort(shared_ptr<vector<Point>> points, shared_ptr<map<int, vector<int>>> adj_map){

    map<int, int> replace_map;
    queue<int> to_visit;
    vector<vector<int>> cycles;
    vector<int> curr_list, next_list;
    Point curr_point, next_point;
    int iter;

    to_visit.push(0);

    //Traverse the nodes in the graph. Whenever a node is encountered for the
    //first time, add it to the queue. When a node is popped from the queue,
    //add it as a key in the map, and map its position in the sequence in which
    //points are popped from the queue.
    while(! to_visit.empty()){
        replace_map.insert(make_pair(to_visit.front(), replace_map.size()));

        for(int neighbor : (*adj_map)[to_visit.front()]){
            if(replace_map.find(neighbor) == replace_map.end()){
                to_visit.push(neighbor);
            }
        }

        to_visit.pop();
    }

    //Find all cycles in the permutation defined by the map from indices to 
    //their ordering in the breadth-first sorting
    cycles = get_cycles(replace_map);

    //Iterate over cycles. Exchange points in cyclical order. Also resassign
    //neighbor lists, and replace the indices in each neigbor list with their
    //images under the mapping.
    for(vector<int> next_cycle : cycles){
        if(next_cycle.size() > 1){
            curr_point = (*points)[next_cycle[0]];
            next_point = (*points)[next_cycle[1]];
            curr_list = reassign_elems((*adj_map)[next_cycle[0]], replace_map);
            next_list = reassign_elems((*adj_map)[next_cycle[1]], replace_map);

            for(iter = 0; iter < next_cycle.size(); iter++){
                (*points)[next_cycle[(iter+1)%next_cycle.size()]] = curr_point;
                (*adj_map)[next_cycle[(iter+1)%next_cycle.size()]] = curr_list;
                curr_point = next_point;
                curr_list = next_list;
                next_point = (*points)[next_cycle[(iter+2)%next_cycle.size()]];
                next_list = reassign_elems((*adj_map)[next_cycle[(iter+2)%next_cycle.size()]], replace_map);
            }
        }

        else{
            (*adj_map)[next_cycle[0]] = reassign_elems((*adj_map)[next_cycle[0]],replace_map);
        }
    }

}

//Prepare a data structure for providing instructions for the calculation of a
//Deulaunay triangulation of a point set.
void prepare_input(triangulateio *in, double *xcoords, double *ycoords, size_t num_points){
    int iter;

    in->numberofpoints = num_points;
    in->pointlist = (double *) malloc(2 * sizeof(double) * num_points);
    for(iter = 0; iter < num_points; iter++){
        in->pointlist[iter*2] = xcoords[iter];
        in->pointlist[iter*2+1] = ycoords[iter];
    }

    in->numberofpointattributes = 0;
    in->pointmarkerlist = (int *) NULL;
    in->numberofsegments = 0;
    in->numberofholes = 0;
    in->numberofregions = 0;
}

//Prepare a data structure to hold the output after a Delaunay triangulation is
//calculated
void prepare_output(triangulateio *out){
    out->pointlist = (double *) NULL;
    out->numberofpointattributes = 0;
    out->pointmarkerlist = (int *) NULL;
    out->edgelist = (int *) NULL;
    out->edgemarkerlist = (int *) NULL;
}

void free_input(triangulateio *in){
    free(in->pointlist);
}

void free_output(triangulateio *out){
    free(out->pointlist);
    free(out->edgelist);
}

bool on_border(Point p, vector<double> bounds){
    bool result = false;
    result |= abs(p.x - bounds[0]) < FLOAT_TOL;
    result |= abs(p.y - bounds[1]) < FLOAT_TOL;
    result |= abs(p.x - bounds[2]) < FLOAT_TOL;
    result |= abs(p.y - bounds[3]) < FLOAT_TOL;
    return result;
}

/*Calculate a Voronoi tessellation for a random point set, constrained by a
bounding box.
*/
void get_voronoi_tessellation(double *xcoords, double *ycoords, int ip_count, vector<double> bounds,  shared_ptr<vector<Point>> outpoints, shared_ptr<map<int, vector<int>>> adj_map){

    triangulateio in, out;
    map<Point, int> pmap;
    double x1, y1, x2, y2, min_dist = FLT_MAX, curr_dist;
    map<double, int> border_map;
    int pcount = 0, index1, index2, iter;
    //vector<Point> boundary_points;
    Point p1, p2;
    VoronoiDiagramGenerator vdg;
    double midx = (bounds[0]+bounds[2])/2, midy = (bounds[1]+bounds[3])/2;
    double angle;
    time_t startTime, endTime;
    //bool p1in, p2in;

    //Indicate that Triangle should index from zero, output edges, refrain from
    //reporting triangules, and operate quietly
    char options[5];
    sprintf(options, "%s", "zeEQ");

    //Produce a Delaunay triangulation, then use its edges to determine the 
    //minimum separation between points.
    prepare_input(&in, xcoords, ycoords, ip_count);
    prepare_output(&out);

    triangulate(options, &in, &out, (triangulateio *) NULL);

    for(iter = 0; iter < out.numberofedges; iter++){
        x1 = out.pointlist[2*out.edgelist[2*iter]];
        y1 = out.pointlist[2*out.edgelist[2*iter] + 1];
        x2 = out.pointlist[2*out.edgelist[2*iter + 1]];
        y2 = out.pointlist[2*out.edgelist[2*iter + 1] + 1];

        curr_dist = (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1);
        if(curr_dist < min_dist) min_dist = curr_dist;
    }
    min_dist = sqrt(min_dist);

    free_input(&in);
    free_output(&out);

    //Produce a Voronoi tessellation of the point set, bounded by the user-
    //specified bounding box.
    //time(&startTime);
    vdg.generateVoronoi(xcoords, ycoords, ip_count, bounds[0], bounds[2], bounds[1], bounds[3], min_dist);
    //time(&endTime);
    //cout << difftime(endTime, startTime) << " seconds.\n";
    vdg.resetIterator();

    while(vdg.getNext(x1, y1, x2, y2)){
        p1 = Point(x1, y1);
        p2 = Point(x2, y2);

        if(y1 == y2 && x1 == x2) continue;

        //If either point in an edge has not been encountered yet, add it to a
        //map from points to integers.
        if(pmap.find(p1) == pmap.end()){
            pmap.insert(make_pair(p1, pcount));
            (*adj_map).insert(make_pair(pcount, vector<int>()));
            (*outpoints).push_back(p1);
            if(on_border(p1, bounds)){
                angle = atan2(p1.y - midy, p1.x - midx);
                if(angle < 0) angle += 2*M_PI;
                border_map.insert(make_pair(angle, pcount));
            }
            pcount ++;
        }
        if(pmap.find(p2) == pmap.end()){
            pmap.insert(make_pair(p2, pcount));
            (*adj_map).insert(make_pair(pcount, vector<int>()));
            (*outpoints).push_back(p2);
            if(on_border(p2, bounds)){
                angle = atan2(p2.y - midy, p2.x - midx);
                if(angle < 0) angle += 2*M_PI;
                border_map.insert(make_pair(angle, pcount));
            }
            pcount ++;
        }

        //Update the adjacency lists for the two end points of the current edge
        index1 = pmap[p1];
        index2 = pmap[p2];
        (*adj_map)[index1].push_back(index2);
        (*adj_map)[index2].push_back(index1);
    }

    //cout << "Made it here.\n";

    //Place each point in the tessellation in the Voronoi ponit set the first
    //time it is encountered, and assign the point an integer index. Prepare
    //a set of adjacency lists from integer indices of points to the indices of
    //their neighbors.

    /*
    //Set up data structures for computing the Voronoi tessellation
    in.numberofpoints = ip_count;
    in.pointlist = inpoints;
    in.pointmarkerlist = (int *) NULL;
    in.numberofpointattributes = 0;
    in.numberofsegments = 0;
    in.numberofholes = 0;
    in.numberofregions = 0;

    out.pointlist = (REAL *) NULL;
    out.numberofpointattributes = 0;
    out.numberoftriangleattributes = 0;
    out.edgelist = (int *) NULL;
    out.normlist = (REAL *) NULL;

    //Compute the tessellation

    triangulate(options, &in, &dummy, &out);

    //Retrieve all points and eges in the set. Reject duplicate points, which 
    //may be present due to the fact that Triangle does not check for these when
    //computing Voronoi tessellations. If an edge has two non-negative point
    //indices, just add it to the list of edges. If it is a ray, find its
    //intersection with the bounding box, and add just the segment of the ray
    //contained within the bounding box to the list. Also add the intersection
    //point with the boundary to a pool.
    for(iter = 0; iter < out.numberofedges; iter++){

        outIndex1 = out.edgelist[iter * 2];
        p1 = Point(out.pointlist[outIndex1*2],out.pointlist[outIndex1*2+1]);

        if((p1in = inBounds(bounds, p1))){
            if(pmap.find(p1) == pmap.end()){
                pmap.insert(make_pair(p1, pcount));
                listIndex1 = pcount;
                pcount ++;
                (*outpoints).push_back(p1);
            }
            else listIndex1 = pmap[p1];
        }
            
        outIndex2 = out.edgelist[iter * 2 + 1];

        if(outIndex2 >= 0){
            p2 = Point(out.pointlist[outIndex2*2],out.pointlist[outIndex2*2+1]);

            if((p2in = inBounds(bounds, p2))){
                if(pmap.find(p2) == pmap.end()){
                    pmap.insert(make_pair(p2, pcount));
                    listIndex2 = pcount;
                    (*outpoints).push_back(p2);
                    pcount ++;
                }
                else listIndex2 = pmap[p2];
            }

            if(! p1in && ! p2in) continue;

            else if(!p1in && p2in){

                p1 = get_closest_intersection(bounds, p2, p1.x-p2.x,p1.y-p2.y);

                angle = atan2(p1.y - midy, p1.x - midx);
                if(angle < 0) angle += 2*M_PI;
                listIndex1 = pcount;
                pcount ++;
                border_map.insert(make_pair(angle, listIndex1));
                (*outpoints).push_back(p1);
            }

            else if(p1in && !p2in){

                p2 = get_closest_intersection(bounds, p1, p2.x-p1.x,p2.y-p1.y);

                angle = atan2(p2.y - midy, p2.x - midx);
                if(angle < 0) angle += 2*M_PI;
                listIndex2 = pcount;
                pcount ++;
                border_map.insert(make_pair(angle, listIndex2));
                (*outpoints).push_back(p2);
            }

        }

        else if(p1in){
            p2 = get_closest_intersection(bounds, p1, out.normlist[iter*2], out.normlist[iter*2+1]);
            angle = atan2(p2.y - midy, p2.x - midx);
            if(angle < 0) angle += 2*M_PI;
            listIndex2 = pcount;
            pcount ++;
            border_map.insert(make_pair(angle, listIndex2));
            (*outpoints).push_back(p2);
        }

        else continue;

        //Update the adjacency map for the planar line graph
        if((*adj_map).find(listIndex1) == (*adj_map).end()){
            (*adj_map).insert(make_pair(listIndex1, vector<int>()));
        }

        if((*adj_map).find(listIndex2) == (*adj_map).end()){
            (*adj_map).insert(make_pair(listIndex2, vector<int>()));
        }

        (*adj_map)[listIndex1].push_back(listIndex2);
        (*adj_map)[listIndex2].push_back(listIndex1);
    }
    */

    //Once all boundary points have been found, augment this list with the
    //vertices of the bounding box, and sort the boundary points according to 
    //the angle between the positive x axis and a line segment from the bounding
    //box centroid to each point. Connect adjacent pairs of points with line
    //segments, and add these to the list of edges.
 
    //Add bottom left corner
    angle = atan2(bounds[1]-midy, bounds[0]-midx);
    if(angle < 0) angle += 2*M_PI;
    (*outpoints).push_back(Point(bounds[0], bounds[1]));
    border_map.insert(make_pair(angle, pcount));
    (*adj_map).insert(make_pair(pcount, vector<int>()));
    pcount ++;

    //Add top left corner
    angle = atan2(bounds[3]-midy, bounds[0]-midx);
    if(angle < 0) angle += 2*M_PI;
    (*outpoints).push_back(Point(bounds[0], bounds[3]));
    border_map.insert(make_pair(angle, pcount));
    (*adj_map).insert(make_pair(pcount, vector<int>()));
    pcount ++;

    //Add bottom right corner
    angle = atan2(bounds[1]-midy, bounds[2]-midx);
    if(angle < 0) angle += 2*M_PI;
    (*outpoints).push_back(Point(bounds[2], bounds[1]));
    border_map.insert(make_pair(angle, pcount));
    (*adj_map).insert(make_pair(pcount, vector<int>()));
    pcount ++;

    //Add top right corner
    angle = atan2(bounds[3]-midy, bounds[2]-midx);
    if(angle < 0) angle += 2*M_PI;
    (*outpoints).push_back(Point(bounds[2], bounds[3]));
    border_map.insert(make_pair(angle, pcount));
    (*adj_map).insert(make_pair(pcount, vector<int>()));

    //Now add edges around the perimeter of the bounding box
    for(auto biter=border_map.begin(); biter!=prev(border_map.end()); biter++){
        (*adj_map)[biter->second].push_back(next(biter)->second);
        (*adj_map)[next(biter)->second].push_back(biter->second);
    }

    (*adj_map)[border_map.begin()->second].push_back(border_map.crbegin()->second);
    (*adj_map)[border_map.crbegin()->second].push_back(border_map.begin()->second);
    
    //Add bottom left corner
    //Perform a breadth-first search sort of the points and adjacency matrix to
    //imporve locality.
    //bfs_sort(outpoints, adj_map);

    //free(out.pointlist);
    //free(out.edgelist);
    //free(out.normlist);
}

//Find the next edge during a traversal of a graph to find its edges
void get_next(int &i, int &j, shared_ptr<map<int, vector<int>>> adj_map){

    int pos = 0;
    while((*adj_map)[j][pos] != i) pos++;

    pos = (pos + 1) % (*adj_map)[j].size();
    i = j;
    j = (*adj_map)[j][pos];
}

bool face_sort(PolyData f1, PolyData f2){
    return f1.polyEdges.size() < f2.polyEdges.size();
}

/*Extract the faces of a planar graph*/
vector<PolyData> get_faces(shared_ptr<vector<Point>> points, shared_ptr<map<int, vector<int>>> adj_map){

    int iter, initial, current, next;//, face_count = 0;
    unordered_set<tuple<int,int>> visited;
    vector<PolyData> pdata;
    vector<Edge> polyEdges;

    //Sort each point's list of neighbors according to the angle between the 
    //positive x axis and a line segment that point to a neighbor.
    for(iter = 0; iter < (*points).size(); iter ++){
        sort((*adj_map)[iter].begin(),(*adj_map)[iter].end(),adj_map_comparator(points,iter));
    }

    //cout << "Ready to find faces.\n";

    //Iterate over each edge in the graph. If an edge has not been visited yet,
    //use it as the beginning of a graph face, and traverse the face until
    //the search returns to the starting point
    for(auto miter  = (*adj_map).begin(); miter != (*adj_map).end(); miter++){
        for(int neighbor : miter->second){
            if(visited.find(make_tuple(miter->first,neighbor))==visited.end()){
                initial = miter->first;
                current = miter->first;
                next = neighbor;
                visited.emplace(make_tuple(initial, next));
                polyEdges.clear();
                polyEdges.push_back(Edge((*points)[initial], (*points)[next]));
                while(next != initial){
                    get_next(current, next, adj_map);
                    visited.emplace(make_tuple(current, next));
                    polyEdges.push_back(Edge((*points)[current],(*points)[next]));
                }
                pdata.push_back(PolyData(polyEdges));
                //face_count ++;
                //if(face_count % 100 == 0) cout << face_count << " faces.\n";
            }
        }
    }

    sort(pdata.begin(), pdata.end(), face_sort);

    return pdata;
}

/*Create a grid data structure for efficiently locating polygons that may
contain a target point.*/
void make_poly_grid(shared_ptr<vector<PolyData>> pdata, vector<vector<int>> &grid_lists, double length, double &minx, double &miny, double &maxx, double &maxy, int &xdim){

    int data_iter, xStart, yStart, xEnd, yEnd, xIter, yIter, ydim;
    minx = miny = (double) FLT_MAX;
    maxx = maxy = (double) FLT_MIN;

    //Pass through once to find the minimum and maximum x and y coordinates
    for(PolyData data : *pdata){
        if(data.left < minx) minx = data.left;
        if(data.right > maxx) maxx = data.right;
        if(data.low < miny) miny = data.low;
        if(data.high > maxy) maxy = data.high;
    }

    //Determine the number of grid cells along the x and y directions and
    //allocate the according number of grid lists. Insert a terminating -1
    //in each list, after which any appropriate edge indices will be added later
    xdim = (int) ceil((maxx - minx) / length);
    ydim = (int) ceil((maxy - miny) / length);
    for(data_iter = 0; data_iter < xdim * ydim; data_iter++){
        grid_lists.push_back(vector<int>());
    }

    //Make a second pass through to find the grid boundaries of each polygon
    for(data_iter = 0; data_iter < (*pdata).size(); data_iter++){
        PolyData data = (*pdata)[data_iter];
        xStart = (int) floor((data.left - minx) / length);
        xEnd = (int) floor((data.right - minx) / length);
        yStart = (int) floor((data.low - miny) / length);
        yEnd = (int) floor((data.high - miny) / length);

        for(yIter = yStart; yIter <= yEnd; yIter ++){
            for(xIter = xStart; xIter <= xEnd; xIter ++){
                grid_lists[yIter * xdim + xIter].push_back(data_iter);
            }
        }
    }
}

/*Given a set of Voronoi tiles and points, match each point with its
corresponding tile.*/
void assign_partners(double *xcoords, double *ycoords, int num_points, shared_ptr<vector<PolyData>> cells, int *cell_map, vector<double> bounds, double spacing){

    vector<vector<int>> grid_lists;
    int xdim, ydim, xStart, xEnd, yStart, yEnd, xIter, yIter, pIter, matches = 0;
    vector<int> curr_list;
    bool match_found;
    Point currPoint;

    make_poly_grid(cells, grid_lists, spacing, bounds[0], bounds[1], bounds[2], bounds[3], xdim);

    for(pIter = 0; pIter < num_points; pIter ++){

        currPoint = Point(xcoords[pIter], ycoords[pIter]);
        match_found = false;
        xStart  = (int) floor((currPoint.x - bounds[0]) / spacing);
        xEnd = (int) ceil((currPoint.x - bounds[0]) / spacing);
        yStart  = (int) floor((currPoint.y - bounds[1]) / spacing);
        yEnd = (int) ceil((currPoint.y - bounds[1]) / spacing);

        for(yIter = yStart; yIter <= yEnd; yIter ++){
            for(xIter = xStart; xIter <= xEnd; xIter ++){
                curr_list = grid_lists[yIter * xdim + xIter];

                for(int pdat_index : curr_list){
                    if((*cells)[pdat_index].contains(currPoint)){
                        match_found = true;
                        matches ++;
                        cell_map[pIter] = pdat_index;
                        break;
                    }
                }
                if(match_found) break;
            }
            if(match_found) break;
        }

    }

    //cout << "Matches: " << matches << "\n";
}

void assign_partners_mt(double *xcoords, double *ycoords, int begin, int end, shared_ptr<vector<PolyData>> cells, int *cell_map, vector<double> bounds, double spacing, shared_ptr<vector<vector<int>>> grid_lists, int xdim){

    int xStart, xEnd, yStart, yEnd, xIter, yIter, pIter;
    vector<int> curr_list;
    bool match_found;
    Point currPoint;

    for(pIter = begin; pIter <= end; pIter ++){

        currPoint = Point(xcoords[pIter], ycoords[pIter]);
        match_found = false;
        xStart  = (int) floor((currPoint.x - bounds[0]) / spacing);
        xEnd = (int) ceil((currPoint.x - bounds[0]) / spacing);
        yStart  = (int) floor((currPoint.y - bounds[1]) / spacing);
        yEnd = (int) ceil((currPoint.y - bounds[1]) / spacing);

        for(yIter = yStart; yIter <= yEnd; yIter ++){
            for(xIter = xStart; xIter <= xEnd; xIter ++){
                curr_list = (*grid_lists)[yIter * xdim + xIter];

                for(int pdat_index : curr_list){
                    if((*cells)[pdat_index].contains(currPoint)){
                        match_found = true;
                        cell_map[pIter] = pdat_index;
                        break;
                    }
                }
                if(match_found) break;
            }
            if(match_found) break;
        }

    }

}

void partner_mt_wrapper(double *xcoords, double *ycoords, int num_points, shared_ptr<vector<PolyData>> cells, int *cell_map, vector<double> bounds, double spacing, int nthreads){

    vector<int> jobs;
    int jiter, step = num_points / nthreads, xdim;
    vector<thread> workers;
    shared_ptr<vector<vector<int>>> grid_lists = make_shared<vector<vector<int>>>();

    for(jiter = 0; jiter < num_points; jiter += step){
        jobs.push_back(jiter);
        jobs.push_back(min(jiter + step - 1, num_points - 1));
    }

    make_poly_grid(cells, *grid_lists, spacing, bounds[0], bounds[1], bounds[2], bounds[3], xdim);

    for(jiter = 0; jiter < nthreads; jiter ++){
        workers.push_back(thread(assign_partners_mt, xcoords, ycoords, jobs[2*jiter], jobs[2*jiter+1], cells, cell_map, bounds, spacing, grid_lists, xdim));
    }

    for(auto &th : workers) th.join();
}

/*Triangulate a Voronoi tile, find its ceter of mass, and compute the integral
of the square of the distance from the corresponding point in the random point
set over the area of the tile.*/
void calc_centroids_cost(double *xcoords, double *ycoords, int start, int end, shared_ptr<vector<PolyData>> faces, int *vorPairs, double &cost, shared_ptr<vector<Point>> centroids){

    //Structures for triangulation of polygons
    //triangulateio in, out;
    //Switches for triangulation: Indicate a set of segments will be provided,
    //indexing should start from zero, and execution should be quiet
    //char options[] = "pzBPQ";
    int iter, triIter, num_vertices, v1, v2, v3, offset;
    double area, curr_area, x1, y1, x2, y2, x3, y3, cx, cy;
    vector<Edge> edge_list;

    cost = 0;

    /*in.pointmarkerlist = (int *) NULL;
    in.numberofpointattributes = 0;
    in.segmentmarkerlist = (int *) NULL;
    in.numberofholes = 0;
    in.numberofregions = 0;*/

    //Triangulate faces one-by-one, find the centroid of each triangle, evaluate
    //the cost function over each triangle, then combine the results
    for(iter = start; iter <= end; iter++){
        cx = 0;
        cy = 0;
        area = 0;

        edge_list = (*faces)[vorPairs[iter]].polyEdges;
        num_vertices = edge_list.size();
        offset = 3*(num_vertices-3)*(num_vertices-2)/2;

        /*prepare_input(&in, (*faces)[vorPairs[iter]].polyEdges);
        prepare_output(&out);
        triangulate(options, &in, &out, (triangulateio *) NULL);*/

        for(triIter = 0; triIter < num_vertices-2; triIter++){
            v1 = triangle_lists[3*triIter+offset];
            v2 = triangle_lists[3*triIter+1+offset];
            v3 = triangle_lists[3*triIter+2+offset];
            x1 = edge_list[v1].p1.x;
            y1 = edge_list[v1].p1.y;
            x2 = edge_list[v2].p1.x;
            y2 = edge_list[v2].p1.y;
            x3 = edge_list[v3].p1.x;
            y3 = edge_list[v3].p1.y;

            curr_area = abs((x1-x3)*(y2-y3) - (x2-x3)*(y1-y3));
            area += curr_area;
            cx += curr_area*(x1 + x2 + x3) / 3;
            cy += curr_area*(y1 + y2 + y3) / 3;

            cost += curr_area*(x1*x1+x2*x2+x3*x3 + x1*x2+x1*x3+x2*x3)/12;
            cost += curr_area*xcoords[iter]*(-(x1+x2+x3)/3 + xcoords[iter]/2);
            cost += curr_area*(y1*y1+y2*y2+y3*y3 + y1*y2+y1*y3+y2*y3)/12;
            cost += curr_area*ycoords[iter]*(-(y1+y2+y3)/3+ycoords[iter]/2);
        }
        (*centroids)[iter].x = cx / area;
        (*centroids)[iter].y = cy / area;
        //free_input(&in);
        //free_output(&out);
    }
}

//Split the calculation of centroids and the cost function over two or more
//worker threads.
void calc_centroids_cost_mt(double *xcoords, double *ycoords, int num_points, shared_ptr<vector<PolyData>> faces, int *vorPairs, double &cost, shared_ptr<vector<Point>> centroids, int nthreads){

    vector<double> cost_vec(nthreads, 0);
    vector<int> jobs;
    int jiter, step = num_points / nthreads;
    vector<thread> workers;

    cost = 0;

    for(jiter = 0; jiter < num_points; jiter += step){
        jobs.push_back(jiter);
        jobs.push_back(min(jiter + step - 1, num_points - 1));
    }

    for(jiter = 0; jiter < nthreads; jiter ++){
        workers.push_back(thread(calc_centroids_cost, xcoords, ycoords, jobs[2*jiter], jobs[2*jiter+1], faces, vorPairs, ref(cost_vec[jiter]), centroids));
    }

    for(auto &th : workers) th.join();

    for(double next_cost : cost_vec) cost += next_cost;
}

void move_to_centroids(double *xcoords, double *ycoords, shared_ptr<vector<Point>> centroids){
    int iter;

    for(iter = 0; iter < (*centroids).size(); iter++){
        xcoords[iter] = (*centroids)[iter].x;
        ycoords[iter] = (*centroids)[iter].y;
    }
}

/*Given a random point set, compute a corresponding Voronoi tessellation of a box
containing that point set,  compute the Lloyd cost function for the tessellation,
replace each point by the center of mass of its Voronoi tile, and repeat until
either the cost function converges to a specified threshold, or a maximum
number of iterations have been performed.*/
void lloyd_optimize(double *xcoords, double *ycoords, int pcount, vector<double> bounds, double density, double max_diff, int max_iter, int nthreads){

    int iter = 0, piter;
    double cost, prevCost = FLT_MAX, diff = FLT_MAX, spacing = sqrt(1/density);
    shared_ptr<vector<Point>> centroids = make_shared<vector<Point>>();
    shared_ptr<vector<Point>> vorPoints = make_shared<vector<Point>>();
    shared_ptr<map<int, vector<int>>> adj_map = make_shared<map<int, vector<int>>>();
    int *vorPairs = (int *) malloc(sizeof(int) * pcount);
    //ofstream tesfile;
    //FILE *log_file;
    string log_name;
    shared_ptr<vector<PolyData>> faces = make_shared<vector<PolyData>>();
    vector<int> jobs;
    //bool use_log_file;
    time_t start, end;

    *centroids = vector<Point>(pcount, Point());

    /*if((use_log_file = yesno("Use a log file?"))){
        cout << "Enter the log name: ";
        getline(cin, log_name);
    }*/

    while(iter < max_iter && diff > max_diff){
        (*vorPoints).clear();
        (*adj_map).clear();


        //Compute edges of a Voronoi tessellation
        //time(&start);
        get_voronoi_tessellation(xcoords, ycoords, pcount, bounds, vorPoints, adj_map);
        //time(&end);
        //cout << difftime(end, start) << " seconds for finding tessellation.\n";
        //cout << "Number of Voronoi points: " << (*vorPoints).size() << "\n";

        /*for(auto amIter=((*adj_map).begin)();amIter!=(*adj_map).end();amIter++){
            for(int neighbor : amIter->second){
                if(neighbor > amIter->first){
                    cout << (*vorPoints)[amIter->first].x << "\t";
                    cout << (*vorPoints)[amIter->first].y << "\n";
                    cout << (*vorPoints)[neighbor].x << "\t";
                    cout << (*vorPoints)[neighbor].y << "\n\n";
                }
            }
        }*/

        //Find Voronoi tiles and match them with their corresponding points
        //time(&start);
        *faces = get_faces(vorPoints, adj_map);
        //time(&end);
        //cout << difftime(end, start) << " seconds for extracting faces.\n";

        /*if(yesno("Report faces?")){
            open_output_file("Enter the file name: ", tesfile);
            if(tesfile.is_open()){
                for(PolyData nextFace : faces){
                    for(Edge e : nextFace.polyEdges){
                        tesfile << e.p1.x << "\t" << e.p1.y << "\n";
                        tesfile << e.p2.x << "\t" << e.p2.y << "\n\n";
                    }
                }
                tesfile.close();
            }
        }*/

        //cout << "Got faces.\n";

        //time(&start);
        if(nthreads == 1){
            assign_partners(xcoords, ycoords, pcount, faces, vorPairs, bounds, spacing);
        }
        else{
            partner_mt_wrapper(xcoords, ycoords, pcount, faces, vorPairs, bounds, spacing, nthreads);
        }

        //cout << "Got partners.\n";

        //time(&end);
        //cout << difftime(end, start) << " seconds for assigning partners.\n";

        /*if(yesno("Report tessellation?")){
        open_output_file("Enter the name for reporting: ", tesfile);
        if(tesfile.is_open()){
            for(auto aiter = (*adj_map).begin(); aiter != (*adj_map).end(); aiter++){
                for(int nextNeighbor : aiter->second){
                    if(nextNeighbor > aiter->first){
                        tesfile << (*vorPoints)[aiter->first].x << "\t";
                        tesfile << (*vorPoints)[aiter->first].y << "\n";
                        tesfile << (*vorPoints)[nextNeighbor].x << "\t";
                        tesfile << (*vorPoints)[nextNeighbor].y << "\n\n";
                    }
                }
            }
        }
        tesfile.close();
        }*/


        //time(&start);
        //Find centroids of Voronoi tiles and evaluate the Lloyd cost function
        if(nthreads == 1){
            calc_centroids_cost(xcoords, ycoords, 0, pcount-1, faces, vorPairs, cost, centroids);
        }
        else{
            calc_centroids_cost_mt(xcoords, ycoords, pcount, faces, vorPairs, cost, centroids, nthreads);
        }
        //time(&end);
        //cout << difftime(end, start) << " seconds for centroids and cost.\n";

        //Log the decrease in the cost function and check for convergence
        diff = prevCost - cost;
        prevCost = cost;

        //Replace points with the centroids of their associated Voronoi tiles
        move_to_centroids(xcoords, ycoords, centroids);

        iter++;

        /*if(use_log_file){
            log_file = fopen(log_name.c_str(), "a");
            fprintf(log_file, "%d\t%12.8f\n", iter, cost);
            fclose(log_file);
        }
        else printf("%d\t%12.8f\n", iter, cost);*/
    }

    /*if(yesno("Report tessellation?")){
        open_output_file("Enter the name for reporting: ", tesfile);
        if(tesfile.is_open()){
            for(auto aiter = (*adj_map).begin(); aiter != (*adj_map).end(); aiter++){
                for(int nextNeighbor : aiter->second){
                    if(nextNeighbor > aiter->first){
                        tesfile << (*vorPoints)[aiter->first].x << "\t";
                        tesfile << (*vorPoints)[aiter->first].y << "\n";
                        tesfile << (*vorPoints)[nextNeighbor].x << "\t";
                        tesfile << (*vorPoints)[nextNeighbor].y << "\n\n";
                    }
                }
            }
        }
        tesfile.close();
    }*/

    free(vorPairs);
}

/*
 * Produce a Delaunay triangulation of a point set, and return the edges.
 */
vector<Edge> get_delaunay_triangulation(double *xcoords, double *ycoords, size_t num_points){

    triangulateio in, out;
    char options[5];
    int iter;
    vector<Edge> edges;
    double x1, y1, x2, y2;

    prepare_input(&in, xcoords, ycoords, num_points);
    prepare_output(&out);
    sprintf(options, "%s", "zeEQ");

    triangulate(options, &in, &out, (triangulateio *) NULL);

    for(iter = 0; iter < out.numberofedges; iter++){
        x1 = out.pointlist[2 * out.edgelist[2*iter]];
        y1 = out.pointlist[2 * out.edgelist[2*iter] + 1];
        x2 = out.pointlist[2 * out.edgelist[2*iter + 1]];
        y2 = out.pointlist[2 * out.edgelist[2*iter + 1] + 1];
        edges.push_back(Edge(Point(x1,y1),Point(x2,y2)));
    }

    free_input(&in);
    free_output(&out);

    return edges;
}

/*
vector<double> get_bounds(string message){
    vector<double> bounds;

    while(true){
        bounds = getdoubles(message);
        if(bounds.size() == 4){
            if((bounds[0] > bounds[2]) || (bounds[1] > bounds[3])){
                fprintf(stderr,"Bounds are improperly ordered.\n");
            }
            else break;
        }
        else{
            fprintf(stderr, "Enter four numbers.\n");
        }
    }

    return bounds;
}
*/

bool in_bounds(vector<double> bounds, Point p){
    return p.x>=bounds[0]&&p.x<=bounds[2]&&p.y>=bounds[1]&&p.y<=bounds[3];
}

/*
 * Trim a set of edges to a desired bounding box
 */
vector<Edge> trim_edges(vector<Edge> original, vector<double> new_bounds){

    vector<Edge> replace;

    for(Edge nextEdge : original){
        if(in_bounds(new_bounds, nextEdge.p1) || in_bounds(new_bounds, nextEdge.p2)) replace.push_back(nextEdge);
    }

    return replace;
}

/*Affinely re-scale the x and y values of a set of points, given an old and
 *a new boudning box.
*/
void rescale_points(double *xcoords, double *ycoords, size_t num_points, vector<double> bounds, vector<double> new_bounds){

    double x_range_ratio, y_range_ratio;
    int iter;

    x_range_ratio = (new_bounds[2]-new_bounds[0])/(bounds[2]-bounds[0]);
    y_range_ratio = (new_bounds[3]-new_bounds[1])/(bounds[3]-bounds[1]);

    for(iter = 0; iter < num_points; iter ++){
        xcoords[iter] = (xcoords[iter]-bounds[0])*x_range_ratio + new_bounds[0];
        ycoords[iter] = (ycoords[iter]-bounds[1])*y_range_ratio + new_bounds[1];
    }
}

/*
 * Given a bounding box, prompt for a density of points, prompt for a 
 * convergence threshold for Lloyd iterations and a maximum number of 
 * iterations. Create a random point set fitting within the bounding box, then
 * hand it off to be evolved until one of the stopping criteria has been met.
 */
vector<Edge> get_lloyd_edges(vector<double> input_bounds, double density, int nthreads){

    vector<double> wide_bounds, intermediate_bounds;
    double max_diff, xrange, yrange, large_xrange, large_yrange;
    double *xcoords, *ycoords;
    int max_iter, num_read, point_count, iter;
    char choice;
    string response;
    ofstream outfile;
    vector<Edge> edges;

    //Obtain the convergence threshold
    while(true){
        cout << "Enter the convergence threshold: ";
        getline(cin, response);
        num_read = sscanf(response.c_str(), "%lf", &max_diff);
        if(num_read < 1) cerr << "Enter a numbers.\n";
        else if(max_diff < 0) cerr << "Enter a positive value.\n";
        else break;
    }

    wide_bounds = input_bounds;
    xrange = input_bounds[2] - input_bounds[0];
    yrange = input_bounds[3] - input_bounds[1];
    wide_bounds[0] -= xrange / 10;
    wide_bounds[1] -= yrange / 10;
    wide_bounds[2] += xrange / 10;
    wide_bounds[3] += yrange / 10;
    point_count = (int) (density * 1.44 * xrange * yrange);
    xcoords = (double *) malloc(point_count * sizeof(double));
    ycoords = (double *) malloc(point_count * sizeof(double));

    //Obtain the iteration limit
    while(true){
        cout << "Enter the maximum number of Lloyd iterations: ";
        getline(cin, response);
        num_read = sscanf(response.c_str(), "%d", &max_iter);
        if(num_read < 0 || max_iter <= 0) cerr << "Enter a positive integer.\n";
        else break;
    }

    //Create the random point set

    //If the desired density may lead to numerical stability issues, use a lower
    //density to produce the point set, the re-scale the points at the end.
    if(density > .1){
        intermediate_bounds = wide_bounds;
        large_xrange = xrange * sqrt(density / .1);
        large_yrange = yrange * sqrt(density / .1);
        intermediate_bounds[2] = intermediate_bounds[0] + large_xrange;
        intermediate_bounds[3] = intermediate_bounds[1] + large_yrange;
        create_rps(intermediate_bounds, point_count, xcoords, ycoords);
    }

    else create_rps(wide_bounds, point_count, xcoords, ycoords);

    if(nthreads > point_count) nthreads = point_count;

    //Optimize the point set
    if(density <= .1){
        lloyd_optimize(xcoords, ycoords, point_count, wide_bounds, density, max_diff, max_iter, nthreads);
    }

    else{
        lloyd_optimize(xcoords, ycoords, point_count, intermediate_bounds, density, max_diff, max_iter, nthreads);
    }

    //Re-scale the points, if necessary
    if(density > .1){
        rescale_points(xcoords, ycoords, point_count, intermediate_bounds, wide_bounds);
    }

    //Triangulate the points
    edges = get_delaunay_triangulation(xcoords, ycoords, point_count);

    //Trim the edges to the original size
    edges = trim_edges(edges, input_bounds);

    free(xcoords);
    free(ycoords);

    return edges;
}
