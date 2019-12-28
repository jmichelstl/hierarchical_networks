/*
This program reads data describing a planar line graph by specifying pairs
of end points for edges, and creates an stl file describing a 3D volume produced
by extruding the 2D network.
*/

#define REAL double
#define VOID int

#include <iostream>
#include <unordered_map>
#include "network_utils.h"
#include <vector>
#include <cmath>
#include <list>
#include <sstream>
#include <algorithm>
#include "triangle.hpp"
#include <stdio.h>

using namespace std;

#define tuple_3d tuple<double, double, double>
#define tuple_2d tuple<double, double>
#define NEG1 -1
#define POS1 1

struct STLPoint{

    STLPoint(double myx, double myy) : x(myx), y(myy) {}

    STLPoint(){
        x = 0;
        y = 0;
    }

    double x, y, tol = 1e-7;

    bool operator == (const STLPoint& point) const{
        return(abs(point.x - x) < tol && abs(point.y - y) < tol);
    }

    friend bool operator < (const STLPoint &p1, const STLPoint &p2){
        if(p2.y > p1.y + p1.tol) return true;
        else if(p1.y > p2.y + p2.tol) return false;
        else if(p2.x > p1.x + p1.tol) return true;
        else return false;
    }
};

void write_header(FILE *stream, unsigned int num_triangles){
    const char *message = "imaginaryguitarnotesandimaginaryvocalsexistonlyintheimaginationoftheimaginerjoe\n";

    fwrite((const void *)message, sizeof(char), 80, stream);
    fwrite((const void *) &num_triangles, sizeof(int), 1, stream);
}

bool in_poly(vector<STLPoint> plist, STLPoint p){
    bool inpoly = false;
    STLPoint p1, p2;
    double x = p.x, y = p.y;
    int iter, pcount = plist.size();

    for(iter = 0; iter < pcount; iter++){
        p1 = plist[iter];
        p2 = plist[(iter+1)%pcount];
        if(p1.y < y && p2.y >= y || p2.y < y && p1.y >= y){
            if(p1.x + (y - p1.y)/(p2.y-p1.y)*(p2.x-p1.x) < x){
                inpoly = !inpoly;
            }
        }
    }

    return inpoly;
}

bool plist_comp(vector<STLPoint> list1, vector<STLPoint> list2){
    double meanx1 = 0, meanx2 = 0, meany1 = 0, meany2 = 0; 

    for(STLPoint p : list1){
        meanx1 += p.x;
        meany1 += p.y;
    }
    meanx1 /= list1.size();
    meany1 /= list1.size();

    for(STLPoint p : list2){
        meanx2 += p.x;
        meany2 += p.y;
    }
    meanx2 /= list2.size();
    meany2 /= list2.size();

    if(meany1 < meany2) return true;
    if(meany2 < meany1) return false;
    if(meanx1 < meanx2) return true;
    return false;
}

/*
if a point is not currently in an index map, add it and map it to the integer
index one greater than the current size of the map.
*/
void update_maps(STLPoint p, map<STLPoint, int> &pmap, map<int, list<int>> &emap){
    int index;

    if(pmap.find(p) == pmap.end()){
        index = pmap.size() + 1;
        pmap.insert(make_pair(p, index));
        emap.insert(make_pair(index, list<int>()));
    }
}

/*
this function processes a file describing a planar line graph by specifying
pairs of endpoints of line segments.

arguments:
datfile - an open stream to a file describing the planar line graph
imap - a map from vertices in the planar line graph to integer indices
emap - a map in which each pair is a mapping from a point's integer index to
the indices of all other points to which that point is connected
*/
void read_edges(ifstream &datfile, map<int, STLPoint> &imap, map<int, list<int>> &emap){

    string nextline;
    STLPoint p1, p2;
    int point_count, index1, index2;
    vector<double> edge_data;
    map<STLPoint, int> pmap;

    while(! datfile.eof()){
        getline(datfile, nextline);

        edge_data = parse_doubles(split(nextline, ' '));
        if(edge_data.size() >= 2){
            if(point_count == 0){
                p1 = STLPoint((double) edge_data[0], (double) edge_data[1]);
                point_count ++;
            }

            else{
                p2 = STLPoint(edge_data[0], edge_data[1]);
                point_count ++;
            }

            if(point_count == 2){
                point_count = 0;

                //ensure that both points in the edge just read are accounted
                //for in both maps
                update_maps(p1, pmap, emap);
                update_maps(p2, pmap, emap);

                //indicate that the two points are now known to be joined by
                //an edge
                index1 = pmap[p1];
                index2 = pmap[p2];
                emap[index1].push_back(index2);
                emap[index2].push_back(index1);
            }
        }
    }

    //reverse the mapping of point to integers to one of integers to points
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

    //build cycles by iterating over entries in the map from indices to
    //lists of neighbor indices
    while(! emap.empty()){
        //start a new cycle by finding the first element of the map and
        //making the key the first element in the cycle
        vector<int> curr_cycle;
        curr_index = emap.begin()->first;
        cycle_complete = false;

        //add successive elements to the current cycle by finding the first
        //neighbor listed for the current vertex in the cycle
        do{
            next_index = emap[curr_index].front();
            curr_cycle.push_back(next_index);

            //remove the current index and the next index from each other's
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
                //when the list to which next_index is mapped has no elements,
                //the current cycle is complete
                cycle_complete = true;
            }
            curr_index = next_index;
        }while(! cycle_complete);
        cycle_list.push_back(curr_cycle);
    }

    return cycle_list;
}

/*
*locate a starting point in a cycle that is not part of a colinear group of
*three or more points
*/
int pick_start(vector<int> cycle, map<int, STLPoint> imap){
    int curr_index, next_index, prev_index;
    STLPoint curr_p, next_p, prev_p;
    double slope1, slope2;
    size_t length = cycle.size();

    for(curr_index = 0; curr_index < length; curr_index ++){
        next_index = (curr_index + 1) % length;
        prev_index = (curr_index  + length - 1) % length;
        curr_p = imap[cycle[curr_index]];
        next_p = imap[cycle[next_index]];
        prev_p = imap[cycle[prev_index]];

        //if the slopes of the next two edges are different, choose the
        //next index as the starting index
        slope1 = (curr_p.y - prev_p.y) / (curr_p.x - prev_p.x);
        slope2 = (next_p.y - curr_p.y) / (next_p.x - curr_p.x);

        if(abs(slope2 - slope1) > FLOAT_TOL){
            return curr_index;
        }
    }

    return curr_index;
}

/*
*produce a more compact list of points for the planar line graph that eliminates
groups of three or more colinear points, and generate the corresponding edges.
*/
void get_polygons(map<int, STLPoint> imap, vector<vector<int>> cycle_list, vector<vector<STLPoint>> &polygons, vector<STLPoint> &border, int &pointcount){
    int start_index, curr_index, next_index, cycle_index = 0;
    int min_cycle = 0;
    STLPoint prev_point, curr_point, next_point;
    double slope1, slope2, dx1, dy1, dx2, dy2, angle1, angle2;
    double min = FLT_MAX, cycle_min;
    size_t length;

    cout << "Size: " << cycle_list.size() << "\n";

    for(vector<int> next_cycle : cycle_list){
        //list of points describing the next polygon in a set
        vector<STLPoint> next_poly;
        length = next_cycle.size();

        //locate a starting index at a corner of a polygon, and not along
        //an edge
        start_index = pick_start(next_cycle, imap);

        //iterate over indices in the current cycle, beginning with the
        //strategically chosen starting point, and coallesce edges, then add
        //end points to the minimal map and update edge list
        curr_index = start_index;
        prev_point = imap[next_cycle[(start_index + length - 1) % length]];
        cycle_min = prev_point.y;

        do{
            next_index = (curr_index + 1) % length;
            curr_point = imap[next_cycle[curr_index]];
            //proposed = next_cycle[curr_index];
            /*if(imap.find(proposed) == imap.end()){
                cerr << "the index " << proposed << " is not in the map.\n";
            }*/
            next_point = imap[next_cycle[next_index]];
            if(curr_point.y < cycle_min) cycle_min = curr_point.y;

            slope1 = (curr_point.y-prev_point.y)/(curr_point.x-prev_point.x);
            slope2 = (next_point.y-curr_point.y)/(next_point.x-curr_point.x);
            if(abs(slope2 - slope1) > FLOAT_TOL){
                next_poly.push_back(curr_point);
            }

            prev_point = curr_point;
            curr_index = next_index;
        }while(curr_index != start_index);

        //determine whether this cycle has the lowest y coordinate of any cycle
        //so far
        if(cycle_min < min){
            min = cycle_min;
            min_cycle = cycle_index;
        }

        polygons.push_back(next_poly);
        cycle_index ++;
        pointcount += next_poly.size();
    }

    //the border will contain the lowest point in the planar line graph. this
    //polygon is set aside so that grips may be added and so that it can be
    //color coded
    border = polygons[min_cycle];
    polygons.erase(polygons.begin() + min_cycle, polygons.begin() + min_cycle + 1);

    sort(polygons.begin(), polygons.end(), plist_comp);
}

//Don't bother coallescing points. Just translate indices into points, and
//find the border
void get_polygons_simple(map<int, STLPoint> imap, vector<vector<int>> cycle_list, vector<vector<STLPoint>> &polygons, vector<STLPoint> &border){
    double miny = FLT_MAX;
    int iter, min_cycle = 0;

    for(iter = 0; iter < cycle_list.size(); iter++){
        vector<STLPoint> next_poly;
        for(int next_vertex : cycle_list[iter]){
            next_poly.push_back(imap[next_vertex]);
            if(imap[next_vertex].y > miny){
                miny = imap[next_vertex].y;
                min_cycle = iter;
            }
        }
        polygons.push_back(next_poly);
    }

    border = polygons[min_cycle];
    polygons.erase(polygons.begin() + min_cycle, polygons.begin() + min_cycle + 1);
}

//Helper function for get_normals
void process_polygon(vector<STLPoint> nextpoly, vector<vector<tuple_3d>> &normals, int case1, int case2){

    double sign, midx, midy, dx, dy, inv_mag;
    int iter, pcount;
    STLPoint p1, p2, guess1, guess2;
    bool is1, is2;
    vector<tuple_3d> normlist;

    pcount = nextpoly.size();
    if(pcount < 3){
        cerr << "That's no polygon!\n";
    }

    p1 = nextpoly[0];
    p2 = nextpoly[1];
    dx = p2.x - p1.x;
    dy = p2.y - p1.y;
    midx = (p1.x + p2.x) / 2;
    midy = (p1.y + p2.y) / 2;
    //find the normal to the first edge pointing into the polygon
    guess1 = STLPoint(midx - dy, midy + dx);
    guess2 = STLPoint(midx + dy, midy - dx);
    do{
        guess1 = STLPoint(midx - dy, midy + dx);
        guess2 = STLPoint(midx + dy, midy - dx);
        is1 = in_poly(nextpoly, guess1);
        is2 = in_poly(nextpoly, guess2);
        if(! (is1 || is2)){
            dx /= 2;
            dy /= 2;
        }
    }while(! (is1 || is2));

    //set the procedure for finding normals for successive edges
    sign = is1 ? case1 : case2;
    inv_mag = 1 / sqrt(dx*dx + dy*dy);
    normlist.push_back(make_tuple((float)(-dy*sign*inv_mag),(float)(dx*sign*inv_mag),0));

    for(iter = 1; iter != pcount; iter ++){
        p1 = nextpoly[iter];
        p2 = nextpoly[(iter + 1)%pcount];
        dx = p2.x - p1.x;
        dy = p2.y - p1.y;
        normlist.push_back(make_tuple((float)(-dy*sign*inv_mag),(float)(dx*sign*inv_mag),0));
    }

    normals.push_back(normlist);
}

void get_normals(vector<vector<STLPoint>> polygons, vector<STLPoint> border, vector<vector<tuple_3d>> &normals){
double sign, midx, midy, dx, dy, inv_mag;

    for(vector<STLPoint> nextpoly : polygons){
        process_polygon(nextpoly, normals, POS1, NEG1);
    }

    //process the border separately. normals should point out of the border,
    //and no interior point is desired.
    process_polygon(border, normals, NEG1, POS1);
}

void get_holes(vector<vector<STLPoint>> polygons, vector<tuple_2d> &holes){
    double midx, midy, dx, dy;
    STLPoint p1, p2, guess1, guess2;
    bool is1, is2;

    for(vector<STLPoint> nextpoly : polygons){
        if(nextpoly.size() < 3){
            cerr << "That's no polygon!\n";
        }

        p1 = nextpoly[0];
        p2 = nextpoly[1];
        dx = p2.x - p1.x;
        dy = p2.y - p1.y;
        midx = (p1.x + p2.x) / 2;
        midy = (p1.y + p2.y) / 2;
        //find the normal to the first edge pointing into the polygon
        guess1 = STLPoint(midx - dy, midy + dx);
        guess2 = STLPoint(midx + dy, midy - dx);
        do{
            guess1 = STLPoint(midx - dy, midy + dx);
            guess2 = STLPoint(midx + dy, midy - dx);
            is1 = in_poly(nextpoly, guess1);
            is2 = in_poly(nextpoly, guess2);
            if(is1){
                holes.push_back(make_tuple(guess1.x, guess1.y));
            }
            else if(is2){
                holes.push_back(make_tuple(guess2.x, guess2.y));
            }
            else{
                dx /= 2;
                dy /= 2;
            }
        }while(! (is1 || is2));
    }
}


void add_segments(int *segmentlist, int pcount, int &plabel, int &seg_index){
    int iter, start_plabel = plabel;

    for(iter = 0; iter < pcount - 1; iter++){
        segmentlist[seg_index] = plabel;
        seg_index++;
        segmentlist[seg_index] = plabel + 1;
        seg_index++;
        plabel ++;
    }
    segmentlist[seg_index] = plabel;
    seg_index++;
    segmentlist[seg_index] = start_plabel;
    seg_index++;
    plabel++;
}

void create_input(triangulateio *in, int point_count, vector<vector<STLPoint>> polygons, vector<STLPoint> border, vector<tuple_2d> holes){
    int pindex = 0, plabel = 0, seg_index = 0, hole_index = 0;

    in->numberofpoints = point_count;
    in->numberofsegments = point_count;
    in->pointlist = (double *) malloc(point_count * 2 * sizeof(double));
    in->segmentlist = (int *) malloc(point_count * 2 * sizeof(int));
    in->numberofregions = 0;

    for(vector<STLPoint> nextpoly : polygons){
        for(STLPoint nextSTLPoint : nextpoly){
            in->pointlist[pindex] = nextSTLPoint.x;
            pindex++;
            in->pointlist[pindex] = nextSTLPoint.y;
            pindex++;
        }
        add_segments(in->segmentlist, nextpoly.size(), plabel, seg_index);
    }
    for(STLPoint nextSTLPoint : border){
        in->pointlist[pindex] = nextSTLPoint.x;
        pindex++;
        in->pointlist[pindex] = nextSTLPoint.y;
        pindex++;
    }
    add_segments(in->segmentlist, border.size(), plabel, seg_index);

    in->segmentmarkerlist = (int *) malloc(sizeof(int) * in->numberofsegments);
    in->numberofpointattributes = 0;
    in->pointmarkerlist = (int *) NULL;
    in->numberofedges = 0;
    in->edgelist = (int *) NULL;
    in->numberofholes = holes.size();
    in->holelist = (double *) malloc(in->numberofholes * 2 * sizeof(double));
    for(tuple_2d nexthole : holes){
        in->holelist[hole_index++] = get<0>(nexthole);
        in->holelist[hole_index++] = get<1>(nexthole);
    }
}

void setup_output(triangulateio *out){
    out->pointlist = (double *) NULL;
    out->numberofpointattributes = 0;
    out->numberoftriangleattributes = 0;
    out->trianglelist = (int *) NULL;
    out->edgelist = (int *) NULL;
    out->edgemarkerlist = (int *) NULL;
    out->segmentlist = (int *) NULL;
    out->pointmarkerlist = (int *) NULL;
    out->segmentmarkerlist = (int *) NULL;
}

void write_point(FILE *stlfile, tuple_3d point){
    float px = (float) get<0>(point), py = (float) get<1>(point), pz = (float) get<2>(point);
    fwrite((const void *) &px, sizeof(float), 1, stlfile);
    fwrite((const void *) &py, sizeof(float), 1, stlfile);
    fwrite((const void *) &pz, sizeof(float), 1, stlfile);
}

void to_farr(tuple_3d normal, float *farr){
    farr[0] = (float) get<0>(normal);
    farr[1] = (float) get<1>(normal);
    farr[2] = (float) get<2>(normal);
}

/*Ensure the order of the vertices of a triangle is consistent with its normal
*vector. This is accomplished by first taking the cross product of two vectors:
*the difference between point 2 and point 1, and the difference betwen point
*3 and point 2. This cross product is dotted with the normal. If the result is
*positive, the order of the points is consistent with the normal vector
*/
double check_orientation(tuple_3d normal, tuple_3d p1, tuple_3d p2, tuple_3d p3){
    double cross_x, cross_y, cross_z;

    cross_x = (get<1>(p2)-get<1>(p1))*(get<2>(p3)-get<2>(p2))*(get<2>(p2)-get<2>(p1))*(get<1>(p3)-get<1>(p2));
    cross_y = (get<2>(p2)-get<2>(p1))*(get<0>(p3)-get<0>(p2))*(get<0>(p2)-get<0>(p1))*(get<2>(p3)-get<2>(p2));
    cross_z = (get<0>(p2)-get<0>(p1))*(get<1>(p3)-get<1>(p2))*(get<1>(p2)-get<1>(p1))*(get<0>(p3)-get<0>(p2));
    return cross_x*get<0>(normal)+cross_y*get<1>(normal)+cross_z*get<2>(normal);
}

void write_side(FILE *stlfile, vector<STLPoint> polygon, vector<tuple_3d> normal_set, float height){

    int vertex_iter, poly_size = polygon.size();
    short unsigned int nada = 0;
    tuple_3d v1, v2, v3;
    double ocheck;
    float farr[3];
    STLPoint p1, p2;

    for(vertex_iter = 0; vertex_iter < poly_size; vertex_iter ++){
        p1 = polygon[vertex_iter];
        p2 = polygon[(vertex_iter + 1) % poly_size];

        //bottom triangle
        to_farr(normal_set[vertex_iter], farr);
        fwrite((const void*) farr,sizeof(float),3,stlfile);
        v1 = make_tuple(p1.x, p1.y, 0);
        v2  = make_tuple(p2.x, p2.y, 0);
        v3  = make_tuple(p2.x, p2.y, height);

        //Check orientation
        ocheck = check_orientation(normal_set[vertex_iter], v1, v2, v3);

        if(ocheck >0){
            write_point(stlfile, v1);
            write_point(stlfile, v2);
        }
        else{
            write_point(stlfile, v2);
            write_point(stlfile, v1);
        }
        write_point(stlfile, v3);
        fwrite((const void *)&nada, sizeof(short unsigned int), 1, stlfile);

        //top triangle
        fwrite((const void*) farr,sizeof(float),3,stlfile);
        v1 = make_tuple(p2.x, p2.y, height);
        v2 = make_tuple(p1.x, p1.y, height);
        v3 = make_tuple(p1.x, p1.y, 0);

        ocheck = check_orientation(normal_set[vertex_iter], v1, v2, v3);

        if(ocheck >0){
            write_point(stlfile, v1);
            write_point(stlfile, v2);
        }
        else{
            write_point(stlfile, v2);
            write_point(stlfile, v1);
        }
        write_point(stlfile, v3);
        fwrite((const void *)&nada, sizeof(short unsigned int), 1, stlfile);
    }
}

void write_as_float(FILE *stream, double *arr){
    float fvals[2];
    fvals[0] = arr[0];
    fvals[1] = arr[1];
    fwrite((const void *) fvals, sizeof(float), 2, stream);
}

void write_stl(FILE *stlfile, triangulateio out, vector<vector<STLPoint>> polygons, vector<STLPoint> border, vector<vector<tuple_3d>> normals,  float height){
    int tri_iter, vertex_iter, v_index, poly_iter;
    unsigned short int nada = 0;
    tuple_3d bnorm = make_tuple(0,0,-1), tnorm = make_tuple(0,0,1);
    float barr[3] = {0,0,-1}, tarr[3] = {0,0,1};
    STLPoint p1, p2;
    vector<tuple_3d> normal_set;
    vector<STLPoint> next_poly;
    double ocheck;
    tuple_3d verts[3];
    float farr[3];

    //report top and bottom triangulations
    for(tri_iter = 0; tri_iter < out.numberoftriangles; tri_iter++){
        //report the bottom
        fwrite((const void *) barr, sizeof(float), 3, stlfile);
        for(vertex_iter = 0; vertex_iter < 3; vertex_iter++){
            v_index = out.trianglelist[tri_iter * 3 + vertex_iter];
            verts[vertex_iter] = make_tuple(out.pointlist[2*v_index],out.pointlist[2*v_index+1],0);
        }
        ocheck = check_orientation(bnorm, verts[0], verts[1], verts[2]);
        if(ocheck >0){
            write_point(stlfile, verts[0]);
            write_point(stlfile, verts[1]);
        }
        else{
            write_point(stlfile, verts[1]);
            write_point(stlfile, verts[0]);
        }
        write_point(stlfile, verts[2]);
        fwrite((const void *) &nada, sizeof(short unsigned int), 1, stlfile);

        //now write the top
        fwrite((const void *) tarr, sizeof(float), 3, stlfile);
        for(vertex_iter = 0; vertex_iter < 3; vertex_iter++){
            v_index = out.trianglelist[tri_iter * 3 + vertex_iter];
            verts[vertex_iter] = make_tuple(out.pointlist[2*v_index],out.pointlist[2*v_index+1], (double) height);
        }
        ocheck = check_orientation(tnorm, verts[0], verts[1], verts[2]);
        if(ocheck >0){
            write_point(stlfile, verts[0]);
            write_point(stlfile, verts[1]);
        }
        else{
            write_point(stlfile, verts[1]);
            write_point(stlfile, verts[0]);
        }
        write_point(stlfile, verts[2]);
        fwrite((const void *) &nada, sizeof(short unsigned int), 1, stlfile);
    }

    //Report triangulation for the sides

    //Start with interior polygons
    for(poly_iter = 0; poly_iter < polygons.size(); poly_iter ++){
        next_poly = polygons[poly_iter];
        normal_set = normals[poly_iter];
        write_side(stlfile, next_poly, normal_set, height);
    }

    //Finally, make the border
    write_side(stlfile, border, *normals.rbegin(), height);

}

void report_triangles(FILE *stream, triangulateio out){
    int tri_iter, index1, index2, index3;
    double *p1, *p2, *p3;
    cout << "The number is " << out.numberoftriangles << "\n";

    for(tri_iter = 0; tri_iter < out.numberoftriangles; tri_iter++){
        index1 = out.trianglelist[tri_iter * 3];
        index2 = out.trianglelist[tri_iter * 3 + 1];
        index3 = out.trianglelist[tri_iter * 3 + 2];
        p1 = out.pointlist + 2*index1;
        p2 = out.pointlist + 2*index2;
        p3 = out.pointlist + 2*index3;
        fprintf(stream, "%f\t%f\n%f\t%f\n\n", p1[0], p1[1], p2[0], p2[1]);
        fprintf(stream, "%f\t%f\n%f\t%f\n\n", p2[0], p2[1], p3[0], p3[1]);
        fprintf(stream, "%f\t%f\n%f\t%f\n\n", p3[0], p3[1], p1[0], p1[1]);
    }
}

void report_segments(FILE *stream, triangulateio in){
    int seg_iter, p1, p2;

    for(seg_iter = 0; seg_iter < in.numberofsegments; seg_iter++){
        p1 = in.segmentlist[seg_iter * 2];
        p2 = in.segmentlist[seg_iter * 2 + 1];
        fprintf(stream, "%lf\t%lf\n%lf\t%lf\n\n", in.pointlist[2*p1], in.pointlist[2*p1+1], in.pointlist[2*p2], in.pointlist[2*p2+1]);
    }
}

void report_holes(FILE *stream, triangulateio in){
    int hole_iter;

    for(hole_iter = 0; hole_iter < in.numberofholes; hole_iter ++){
        fprintf(stream, "%f", in.holelist[hole_iter * 2]);
        fprintf(stream, "\t%f\n", in.holelist[hole_iter * 2 + 1]);
    }
}

void read_boundary_segments(triangulateio *out, map<int, STLPoint> &imap, map<int, list<int>> &emap){
    int iter, p1, p2;
    double x,y;

    for(iter = 0; iter < out->numberofsegments; iter++){
        p1 = out->segmentlist[2*iter];
        p2 = out->segmentlist[2*iter + 1];

        if(emap.find(p1) == emap.end()){
            emap.insert(make_pair(p1, list<int>()));
            x = out->pointlist[2*p1];
            y = out->pointlist[2*p1 + 1];
            imap.insert(make_pair(p1, STLPoint(x,y)));
        }
        if(emap.find(p2) == emap.end()){
            emap.insert(make_pair(p2, list<int>()));
            x = out->pointlist[2*p2];
            y = out->pointlist[2*p2 + 1];
            imap.insert(make_pair(p2, STLPoint(x,y)));
        }
        emap[p1].push_back(p2);
        emap[p2].push_back(p1);
    }
}

/*
*prompt the user for the name of a data file and an output stl file. create
an open input stream, obtain a map of points to indices, and a map of indices
to lists of other indices. use the two maps to build a table of cycles, then
convert the table of cycles to polygon data.
*/
int main(int argc, char **argv){
    ifstream datfile;
    FILE *stl_file, *trifile, *segment_file, *hole_file;
    map<int, STLPoint> imap;
    map<int, list<int>> emap;
    vector<vector<int>> cycle_list;
    vector<vector<STLPoint>> polygons;
    vector<STLPoint> border;
    //specify one point inside each interior polygon for the triangulation
    //process
    vector<tuple_2d> holes;
    //normals for the triangles perpendicular to the horizontal plane
    vector<vector<tuple_3d>> normals;
    int point_count = 0, numread, iter, p1, p2;
    unsigned int trianglecount;
    float height = -1;
    string response;
    char options[7] = "pqQze\0";

    //Structures to hold input and output data for triangulation of the network
    struct triangulateio in, out;

    //Get the thickness of the network
    do{
        cout << "Enter the thickness of the network: ";
        getline(cin, response);
        numread = sscanf(response.c_str(), "%f", &height);
        if(numread != 1 || height <= 0){
            cerr << "Enter a positive number.\n";
        }
    }while(height <= 0);

    //Obtain network data and obtain maps from points to integers and from
    //integers to lists of neighbor indices
    open_dat_file("Enter the data file describing the network: ", datfile);
    if(! datfile.is_open()) exit(1);

    read_edges(datfile, imap, emap);
    datfile.close();

    //Use maps to build cycles described by integer indices in cyclical order
    cycle_list = build_cycles(emap);
    cout << "Cycles built.\n";

    //Resolve the planar line graph describing the network into a set of 
    //polygons
    get_polygons(imap, cycle_list, polygons, border, point_count);
    cout << "Polygons made.\n";

    //Find holes for triangulation and normal vectors for stl file creation
    get_holes(polygons, holes);
    cout << "Holes found.\n";

    //Set up input to Triangle's triangulate function and perform triangulation
    create_input(&in, point_count, polygons, border, holes);
    setup_output(&out);
    if(yesno("Write segments?")){
        do{
            cout << "Enter the file name: ";
            getline(cin, response);
            segment_file = fopen(response.c_str(), "w");
        }while(segment_file == NULL);

        report_segments(segment_file, in);
        fclose(segment_file);
    }
    if(yesno("Write holes?")){
        do{
            cout << "Enter the file name: ";
            getline(cin, response);
            hole_file = fopen(response.c_str(), "w");
        }while(hole_file == NULL);

        report_holes(hole_file, in);
        fclose(hole_file);
    }

    cout << "Attempting to triangulate:\n";
    triangulate(options, &in, &out, (struct triangulateio *) NULL);
    cout << "Returned from Triangle.\n";

    //Write STL file
    do{
        cout << "Enter the name of the .stl file: ";
        getline(cin, response);
        stl_file = fopen(response.c_str(), "wb");
    }while(stl_file == NULL);

    //Now that new points have been added by Triangle, find the new boundary
    //segments, create a new map from boundary points to other boundary points,
    //and create new ordered cycles of boundary points.
    imap.clear();
    emap.clear();
    cycle_list.clear();
    read_boundary_segments(&out, imap, emap);
    cycle_list = build_cycles(emap);
    polygons.clear();
    border.clear();
    point_count = 0;
    get_polygons_simple(imap, cycle_list, polygons, border);
    get_normals(polygons, border, normals);

    //Find the total number of triangles, and report triangles to an stl file
    trianglecount = 2*out.numberofsegments + 2*out.numberoftriangles;
    write_header(stl_file, trianglecount);
    write_stl(stl_file, out, polygons, border, normals,  height);
    fclose(stl_file);

    if(yesno("Report triangles?")){
        do{
            cout << "Enter the file name: ";
            getline(cin, response);
            trifile = fopen(response.c_str(), "w");
        }while(trifile == NULL);

        report_triangles(trifile, out);
        fclose(trifile);
    }

    //Delete triangulateio structures and exit
    free(in.pointlist);
    free(in.segmentlist);
    free(in.segmentmarkerlist);
    free(out.pointlist);
    free(out.edgelist);
    free(out.edgemarkerlist);
    free(out.segmentlist);
    free(out.segmentmarkerlist);
    free(out.trianglelist);
    return 0;
}
