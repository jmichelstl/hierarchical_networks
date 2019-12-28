/*
This code provides a means of adding large-scale position disorder to networks.
The program accepts a large-scale skeleton and a description of the small-scale
structure of a network, then iteratatively distorts the large-scale skeleton.
For each configuration of the skeleton, a large-scale envelope of finite width
is produced. Those small-scale vertices that lie on an edge of the large-scale
envelope are identified, and they are affinely transformed to lie on the image
of the large-scale envelope after distortion. Small-scale vertices lying on
envelope edges are then pinned, and the rest of the small-scale vertices are
relaxed to produce the final network.
*/

#include <iostream>
#include "network_utils.h"
#include <map>
#include <vector>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <unordered_set>
#include <algorithm>

using namespace std;

struct EdgeDatum{

    EdgeDatum(int ival, double len) : index(ival), l0(len){
    }

    int index;
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

void create_triples(map<int, vector<int>> triple_map, vector<double> point_list, vector<int> &triples, vector<double> &ref_angles){

    //Map from an angle a bond makes with the x axis to the index of that
    //bond's terminal point
    map<double, int> ang_map;
    double x1, y1, dx, dy, curr_angle, prev_angle, first_angle;
    int index1;

    for(auto iter = triple_map.begin(); iter != triple_map.end(); iter++){
        //Only consider points with two or more neighbors
        if(iter->second.size() >= 2){
            index1 = iter->first;
            x1 = point_list[index1 * 2];
            y1 = point_list[index1 * 2 + 1];
            for(int index2 : iter->second){
                dx = point_list[index2 * 2] - x1;
                dy = point_list[index2 * 2 + 1] - y1;
                curr_angle = atan2(dx, dy);
                ang_map.insert(make_pair(curr_angle, index2));
            }
            for(auto aIter = ang_map.begin(); aIter != prev(ang_map.end());aIter++){
               triples.push_back(index1);
               triples.push_back(aIter->second); 
               triples.push_back(next(aIter)->second);
               ref_angles.push_back(next(aIter)->first - aIter->first);
            }
            if(iter->second.size() > 2){
                triples.push_back(index1);
                triples.push_back(ang_map.crbegin()->second);
                triples.push_back(ang_map.begin()->second);
                ref_angles.push_back(ang_map.begin()->first + 2*M_PI - ang_map.crbegin()->first);
            }
            ang_map.clear();
        }
    }
}

void read_network(ifstream& datfile, vector<double>& point_list, map<int, vector<EdgeDatum>>& neighbor_map, vector<int> &triples, vector<double> &ref_angles, bool get_triples){

    map<Point, int> pmap;
    string nextline;
    Point p1, p2;
    int pindex = 0, point_count = 0, index1, index2, mindex, maxdex;
    int curr_index, next_index;
    vector<double> edge_data;
    double length, curr_x, curr_y, temp_x, temp_y;
    vector<EdgeDatum> curr_list, next_list;
    map<int, vector<int>> triple_map;

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
                    if(get_triples){
                        triple_map.insert(make_pair(pindex, vector<int>()));
                    }
                    pindex++;
                    /*miny = p1.y < miny - FLOAT_TOL ? p1.y : miny;
                    maxy = p1.y > maxy + FLOAT_TOL ? p1.y : maxy;*/
                }
                if(pmap.find(p2) == pmap.end()){
                    pmap.insert(make_pair(p2, pindex));
                    point_list.push_back(p2.x);
                    point_list.push_back(p2.y);
                    neighbor_map.insert(make_pair(pindex, vector<EdgeDatum>()));
                    if(get_triples){
                        triple_map.insert(make_pair(pindex, vector<int>()));
                    }
                    pindex++;
                    /*miny = p2.y < miny - FLOAT_TOL ? p2.y : miny;
                    maxy = p2.y > maxy + FLOAT_TOL ? p2.y : maxy;*/
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

    //Reorder indices so points appear in order sorted order as defined in the
    //Point struct

    //Start by mapping each index to its new index
    vector<int> permute_map(pmap.size(), 0);
    int new_index = 0;
    for(auto iter = pmap.begin(); iter != pmap.end(); iter++){
        permute_map[iter->second] = new_index++;
    }

    //Pool to keep track of neighbor map indices that still need to be
    //reassigned
    set<int> index_pool(permute_map.begin(), permute_map.end());

    //Reassign neighbor lists
    curr_index = 0;
    curr_list = neighbor_map[0];
    curr_x = point_list[0];
    curr_y = point_list[1];

    while(! index_pool.empty()){
        //Fix neighbor indices
        for(int i = 0; i < curr_list.size(); i++){
            curr_list[i].index = permute_map[curr_list[i].index];
        }

        index_pool.erase(curr_index);

        //Find the destination for the current list
        next_index = permute_map[curr_index];
        if(index_pool.find(next_index) != index_pool.end()){
            //Make a record of displaced neighbor list and coordinates
            next_list = neighbor_map[next_index];
            temp_x = point_list[next_index * 2];
            temp_y = point_list[next_index * 2 + 1];

            //Move the neighbor list and coordinates
            neighbor_map[next_index] = curr_list;
            point_list[next_index * 2] = curr_x;
            point_list[next_index * 2 + 1] = curr_y;

            curr_index = next_index;
            curr_list = next_list;
            curr_x = temp_x;
            curr_y = temp_y;
        }

        else{
            neighbor_map[next_index] = curr_list;
            point_list[next_index * 2] = curr_x;
            point_list[next_index * 2 + 1] = curr_y;

            if(! index_pool.empty()){
                curr_index = *(index_pool.begin());
                curr_list = neighbor_map[curr_index];
                curr_x = point_list[2 * curr_index];
                curr_y = point_list[2 * curr_index + 1];
            }
        }
    }

    //If triples are to be found, in order to compute bending energies and
    //forces, make a map from each point to a list of that point's neighbors.
    //Pass this map to create_triples, which will populate the list of triples,
    //and the accompanying list of rest angles
    if(get_triples){
        for(auto iter = neighbor_map.begin();iter != neighbor_map.end();iter++){
            index1 = iter->first;
            if(triple_map.find(index1) == triple_map.end()){
                triple_map.insert(make_pair(index1, vector<int>()));
            }
            for(EdgeDatum edat : iter->second){
                if(triple_map.find(edat.index) == triple_map.end()){
                    triple_map.insert(make_pair(edat.index, vector<int>()));
                }

                triple_map[index1].push_back(edat.index);
                triple_map[edat.index].push_back(index1);
            }
        }
    }

    create_triples(triple_map, point_list,  triples, ref_angles);
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
                /*miny = p1.y < miny - FLOAT_TOL ? p1.y : miny;
                maxy = p1.y > maxy + FLOAT_TOL ? p1.y : maxy;*/
            }
            else{
                p2 = Point(edge_data[0], edge_data[1]);
                /*miny = p2.y < miny - FLOAT_TOL ? p2.y : miny;
                maxy = p2.y > maxy + FLOAT_TOL ? p2.y : maxy;*/
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

//Given an edge and a set of small-scale vertices, find those small-scale
//vertices that lie along the edge.
void on_edge(vector<Edge> border, vector<double> point_set, vector<vector<int>> &border_map){
    double slope, y_target;
    int p_index, e_index = 0;
    bool use_slope;
    Point min, max;
    double x, y, minx, maxx;
    set<int> matched;

    for(Edge next_edge : border){

        min = next_edge.p1 < next_edge.p2 ? next_edge.p1 : next_edge.p2;
        max = next_edge.p2 < next_edge.p1 ? next_edge.p1 : next_edge.p2;
        border_map.push_back(vector<int>());
        minx = min.x < max.x ? min.x : max.x;
        maxx = min.x > max.x ? min.x : max.x;

        if(min.x != max.x){
            use_slope = true;
            slope = (max.y - min.y) / (max.x - min.x);
        }
        else{
            use_slope = false;
        }

        for(p_index = 0; p_index < point_set.size() / 2; p_index ++){

            if(matched.find(p_index) != matched.end()) continue;

            x = point_set[p_index * 2];
            y = point_set[p_index * 2 + 1];

            //if(y < min.y) continue;
            if(y < min.y - FLOAT_TOL) continue;

            //if(y > max.y) break;
            if(y > max.y + FLOAT_TOL) break;

            if(use_slope){
                y_target = min.y + slope * (x - min.x);
                if(abs(y - y_target) < FLOAT_TOL){
                    if(y > min.y - FLOAT_TOL && y < max.y + FLOAT_TOL && x > minx - FLOAT_TOL && x < maxx + FLOAT_TOL){
                        border_map[e_index].push_back(p_index);
                        matched.insert(p_index);
                    }
                }
            }

            else{
                if(abs(min.x - x) < FLOAT_TOL && y >= min.y && y <= max.y){
                    border_map[e_index].push_back(p_index);
                    matched.insert(p_index);
                }
            }
        }

        e_index ++;
    }
}

/*
Given edges with end points p1 and p2, find the inverse of the matrix
|p1x	p2x|
|p1y	p2y|
for each edge for the purpose of computing an affine transform to map points 
along that edge to their image after that edge is transformed to a distorted 
state
*/
vector<vector<double>> get_inv_mats(vector<Edge> edge_list){
    double inv_det;
    double m11, m12, m21, m22;
    vector<vector<double>> imat_list;
    int index = 0;

    for(Edge next : edge_list){
        inv_det = 1 / (next.p1.x * next.p2.y - next.p1.y * next.p2.x);
        m11 = inv_det * next.p2.y;
        m12 = -inv_det * next.p2.x;
        m21 = -inv_det * next.p1.y;
        m22 = inv_det * next.p1.x;
        imat_list.push_back(vector<double>());
        imat_list[index].push_back(m11);
        imat_list[index].push_back(m12);
        imat_list[index].push_back(m21);
        imat_list[index].push_back(m22);
        index++;
    }

    return imat_list;
}        

/*
Given a list of deformed envelope edges, and reference coordinates for points
lying on those images, find the image of each point lying along an edge in
the reference state
*/
void get_deformed_vertices(vector<Edge> deformed, vector<vector<double>> imats,vector<vector<int>> edge_map, map<int, vector<double>> ref_coords, vector<double> &coords){
    //Elements of a transformation matrix from reference coordinates to
    //deformed coordinates
    double t11, t12, t21, t22;
    int edge_iter;
    vector<double> inv_mat;
    vector<double> ref_coord;
    Edge curr_edge;

    cout << "Size: " << deformed.size() << "\n";

    //Affinely transform points lying along each edge of the large-scale
    //envelope
    for(edge_iter = 0; edge_iter < deformed.size(); edge_iter++){
        inv_mat = imats[edge_iter];
        curr_edge = deformed[edge_iter];
        t11 = curr_edge.p1.x * inv_mat[0] + curr_edge.p2.x * inv_mat[2];
        t12 = curr_edge.p1.x * inv_mat[1] + curr_edge.p2.x * inv_mat[3];
        t21 = curr_edge.p1.y * inv_mat[0] + curr_edge.p2.y * inv_mat[2];
        t22 = curr_edge.p1.y * inv_mat[1] + curr_edge.p2.y * inv_mat[3];

        for(int p_index : edge_map[edge_iter]){
            ref_coord = ref_coords[p_index];
            coords[2 * p_index] = t11 * ref_coord[0] + t12 * ref_coord[1];
            coords[2 * p_index + 1] = t21 * ref_coord[0] + t22 * ref_coord[1];
        }
    }
}

double get_pe(map<int, vector<EdgeDatum>> neighbor_map, vector<double> pos, vector<int> triples, vector<double> ref_angles, double ks, double kb){
    double pe = 0;
    int num_points = pos.size() / 2, index1, index2, index3, titer;
    double distance, diff, x1, y1, x2, y2, x3, y3, dx2, dy2, dx3, dy3;
    double curr_angle, ref_angle;

    for(index1 = 0; index1 < num_points; index1 ++){
        x1 = pos[2*index1];
        y1 = pos[2*index1 + 1];
        for(EdgeDatum edat : neighbor_map[index1]){
            index2 = edat.index;
            x2 = pos[2*index2];
            y2 = pos[2*index2 + 1];
            distance = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
            diff = distance - edat.l0;
            pe += .5 * ks * diff * diff;
        }
    }

    if(kb > 0){
        for(titer = 0; titer < triples.size() / 3; titer++){
            index1 = triples[titer * 3];
            index2 = triples[titer * 3 + 1];
            index3 = triples[titer * 3 + 2];

            x1 = pos[2*index1];
            y1 = pos[2*index1 + 1];
            x2 = pos[2*index2];
            y2 = pos[2*index2 + 1];
            x3 = pos[2*index3];
            y3 = pos[2*index3 + 1];

            dx2 = x2 - x1;
            dy2 = y2 - y1;
            dx3 = x3 - x1;
            dy3 = y3 - y1;

            curr_angle = atan2((x3-x1)*(y1-y2)+(y3-y1)*(x2-x1),(x3-x1)*(x2-x1)+(y3-y1)*(y2-y1));
            ref_angle = ref_angles[titer];
            pe += .5 * kb * (curr_angle - ref_angle) * (curr_angle - ref_angle);
        }
    }

    return pe;
}

void get_forces(map<int, vector<EdgeDatum>> neighbor_map, vector<double> pos, vector<double>& forces, vector<int> triples, vector<double> ref_angles, unordered_set<int> pinned, double ks, double kb){
    int index1, index2, index3, titer;
    double x1, x2, y1, y2, x3, y3, l0, dist, fmult, fx, fy, dx2, dy2, dx3, dy3;
    double curr_angle, ref_angle;

    forces.assign(pos.size(), 0);

    for(index1 = 0; index1 < pos.size() / 2; index1++){
        x1 = pos[index1*2];
        y1 = pos[index1*2 + 1];
        for(EdgeDatum edat : neighbor_map[index1]){
            index2 = edat.index;
            l0 = edat.l0;
            x2 = pos[2*index2];
            y2 = pos[2*index2 + 1];
            dist = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2) * (y1 - y2));
            fmult = ks * (dist - l0)/dist;
            fx = fmult * (x2 - x1);
            fy = fmult * (y2 - y1);
            if(pinned.find(index1) == pinned.end()){
                forces[index1*2] += fx;
                forces[index1*2 + 1] += fy;
            }
            if(pinned.find(index2) == pinned.end()){
                forces[index2*2] -= fx;
                forces[index2*2 + 1] -= fy;
            }
        }
    }

    if(kb > 0){
        for(titer = 0; titer < triples.size() / 3; titer++){
            index1 = triples[titer * 3];
            index2 = triples[titer * 3 + 1];
            index3 = triples[titer * 3 + 2];

            x1 = pos[2*index1];
            y1 = pos[2*index1 + 1];
            x2 = pos[2*index2];
            y2 = pos[2*index2 + 1];
            x3 = pos[2*index3];
            y3 = pos[2*index3 + 1];

            dx2 = x2 - x1;
            dy2 = y2 - y1;
            dx3 = x3 - x1;
            dy3 = y3 - y1;

            curr_angle = atan2((x3-x1)*(y1-y2)+(y3-y1)*(x2-x1),(x3-x1)*(x2-x1)+(y3-y1)*(y2-y1));
            ref_angle = ref_angles[titer];

            fmult = kb * (curr_angle - ref_angle)/(dx2*dx2 + dy2*dy2);
            fx = fmult * (y1 - y2);
            fy = fmult * (x2 - x1);
            if(pinned.find(index2) == pinned.end()){
                forces[2*index2] += fx;
                forces[2*index2 + 1] += fy;
            }
            if(pinned.find(index1) == pinned.end()){
                forces[2*index1] -= fx;
                forces[2*index1 + 1] -= fy;
            }

            fmult = kb * (curr_angle - ref_angle)/(dx3*dx3 + dy3*dy3);
            fx = fmult * (y3 - y1);
            fy = fmult * (x1 - x3);
            if(pinned.find(index3) == pinned.end()){
                forces[2*index3] += fx;
                forces[2*index3 + 1] += fy;
            }
            if(pinned.find(index1) == pinned.end()){
                forces[2*index1] -= fx;
                forces[2*index1 + 1] -= fy;
            }
        }
    }
}

void report_deformed(map<int, vector<EdgeDatum>> neighbor_map, vector<double> pos, string filename){
    ofstream datfile;
    int index1, index2;

    datfile.open(filename);

    for(index1 = 0; index1 < pos.size(); index1 ++){
        for(EdgeDatum edat : neighbor_map[index1]){
            index2 = edat.index;
            datfile << pos[2*index1] << " " << pos[2*index1 + 1] << "\n";
            datfile << pos[2*index2] << " " << pos[2*index2 + 1] << "\n\n";
        }
    }

    datfile.close();
}

void getangles(vector<double> angvec,double angle, double& low,double& high){
    int index = 0, size = angvec.size();

    while(index < size-1 && angle-angvec[index] > FLOAT_TOL) index++;
    
    low = angvec[(index+size-1)%size];
    high = angvec[(index+size+1)%size];
}

void changes(double low, double high, double hwidth, double& dx, double& dy){

    double diff = (high-low)/2;
    if(diff < 0) diff += M_PI;

    if(diff == 0){
        fprintf(stderr, "Illegal argument to function changes.\n");
        cerr << "Low: " << low << " High: " << high << "\n";
        return;
    }

    dx = hwidth * (cos(low)/tan(diff) - sin(low));
    dy = hwidth * (sin(low)/tan(diff) + cos(low));
}

void flush_with_edge(Point &pnt, double slope, double y_ext){
    double x_old, y_old, x_new;

    x_old = pnt.x;
    y_old = pnt.y;
    x_new = slope != 0 ? (y_ext - y_old) / slope + x_old : x_old;

    pnt.x = x_new;
    pnt.y = y_ext;
}

vector<Edge> add_thickness(vector<Edge> current, double thickness, vector<PolyData>& pdata, bool makepoly, bool level){
    map<Point, vector<double>> angmap;
    vector<Edge> replace;
    Point p1, p2, key;
    Point p1f, p2f, p3f, p4f, p1fb, p2fb, p3fb, p4fb;
    double low, high, ang1, ang2, dx, dy, hwidth, slope;
    double ymin = FLT_MAX, ymax = FLT_MIN, ylow, yhigh, y_ext;
    bool p1fflag, p2fflag, p3fflag, p4fflag, p1_is_end, p2_is_end;
    vector<double> angvec;

    cout << "Current size: " << current.size() << "\n";

    hwidth = thickness/2;

    for(Edge e : current){
        p1 = e.p1;
        p2 = e.p2;

        ymin = p1.y < ymin ? p1.y : ymin;
        ymin = p2.y < ymin ? p2.y : ymin;
        ymax = p1.y > ymax ? p1.y : ymax;
        ymax = p2.y > ymax ? p2.y : ymax;

        ang1 = atan2(p2.y-p1.y,p2.x-p1.x);
        if(ang1 < 0) ang1 += 2*M_PI;
        ang2 = ang1 < M_PI ? ang1 + M_PI : ang1 - M_PI; 

        if(angmap.find(p1) == angmap.end()){
            angmap.insert(make_pair(p1,vector<double>()));
        }
        if(angmap.find(p2) == angmap.end()){
            angmap.insert(make_pair(p2,vector<double>()));
        }

        angmap[p1].push_back(ang1);
        angmap[p2].push_back(ang2);
    }

    ylow = ymin - hwidth;
    yhigh = ymax + hwidth;

    for(auto iter = angmap.begin(); iter != angmap.end(); iter++){
        key = iter->first;
        sort(angmap[key].begin(), angmap[key].end());
    }

    /*if(v){
        for(auto iter = angmap.begin(); iter != angmap.end(); iter++){
            key = iter->first;
            angvec = angmap[key];
            printf("Angles for point (%lf,%lf):\n", key.x, key.y);
            for(auto iter2 = angvec.begin(); iter2!= angvec.end(); iter2++){
                printf("%lf ", *iter2);
            }
            printf("\n");
        }
    }*/

    for(Edge e : current){
        p1 = e.p1;
        p2 = e.p2;
        ang1 = atan2(p2.y-p1.y,p2.x-p1.x);
        if(ang1 < 0) ang1 += 2*M_PI;
        ang2 = ang1 < M_PI ? ang1 + M_PI : ang1 - M_PI;
        slope = (p2.y - p1.y)/(p2.x - p1.x); 
        p1fflag = false;
        p2fflag = false;
        p3fflag = false;
        p4fflag = false;
        p1_is_end = false;
        p2_is_end = false;

        getangles(angmap[p1], ang1, low, high);
        if(ang1 == low){
            p1_is_end = true;
            p1f = Point(p1.x + hwidth*sin(ang1), p1.y - hwidth*cos(ang1));
            p3f = Point(p1.x - hwidth*sin(ang1), p1.y + hwidth*cos(ang1));
            
            if((p1.y == ymin || p1.y == ymax) && slope != 0){
                y_ext = p1.y == ymin ? ylow : yhigh;
                flush_with_edge(p1f, slope, y_ext);
                flush_with_edge(p3f, slope, y_ext);
            }
            
            replace.push_back(Edge(p1f,p3f));
        }
        else{
            changes(low, ang1, hwidth, dx, dy);
            p1f = Point(p1.x + dx, p1.y + dy);
            changes(ang1, high, hwidth, dx, dy);
            p3f = Point(p1.x + dx, p1.y + dy);

            if((p1f.y < ylow - FLOAT_TOL || p1f.y > yhigh + FLOAT_TOL) && level){
                y_ext = p1.y == ymin ? ylow : yhigh;
                flush_with_edge(p1f, slope, y_ext);
                p1fb = Point(p1.x, y_ext);
                p1fflag = true;
            }
            if((p3f.y < ylow - FLOAT_TOL || p3f.y > yhigh + FLOAT_TOL) && level){
                y_ext = p1.y == ymin ? ylow : yhigh;
                flush_with_edge(p3f, slope, y_ext);
                p3fb = Point(p1.x, y_ext);
                p3fflag = true;
            }
        }
 
        getangles(angmap[p2], ang2, low, high);
        if(ang2 == low){
            p2_is_end = true;
            p2f = Point(p2.x + hwidth*sin(ang1), p2.y - hwidth*cos(ang1));
            p4f = Point(p2.x - hwidth*sin(ang1), p2.y + hwidth*cos(ang1));

            if((p2.y == ymin || p2.y == ymax) && slope != 0){
                y_ext = p1.y == ymin ? ylow : yhigh;
                flush_with_edge(p2f, slope, y_ext);
                flush_with_edge(p4f, slope, y_ext);
            }
            
            replace.push_back(Edge(p2f,p4f));
        }
        else{
            changes(ang2, high, hwidth, dx, dy);
            p2f = Point(p2.x + dx, p2.y + dy);
            changes(low, ang2, hwidth, dx, dy);
            p4f = Point(p2.x + dx, p2.y + dy);

            if((p2f.y < ylow - FLOAT_TOL || p2f.y > yhigh + FLOAT_TOL) && level){
                y_ext = p2.y == ymin ? ylow : yhigh;
                flush_with_edge(p2f, slope, y_ext);
                p2fb = Point(p2.x, y_ext);
                p2fflag = true;
            }
            if((p4f.y < ylow - FLOAT_TOL || p4f.y > yhigh + FLOAT_TOL) && level){
                //cerr << "Current p4f: " << p4f.x << "\t" << p4f.y << "\n";
                y_ext = p2.y == ymin ? ylow : yhigh;
                flush_with_edge(p4f, slope, y_ext);
                p4fb = Point(p2.x, y_ext);
                p4fflag = true;
                //cerr << "New p4f: " << p4f.x << "\t" << p4f.y << "\n\n";
            }
        }

        if(p1fflag) replace.push_back(Edge(p1fb, p1f));
        replace.push_back(Edge(p1f,p2f));
        if(p2fflag) replace.push_back(Edge(p2f, p2fb));
        if(p3fflag) replace.push_back(Edge(p3fb, p3f));
        replace.push_back(Edge(p3f,p4f));
        if(p4fflag) replace.push_back(Edge(p4f, p4fb));

        if(makepoly){
            vector<Edge> pedges;
            pedges.push_back(Edge(p1f,p2f));

            if(!p2_is_end){
                if(p2fflag){
                    pedges.push_back(Edge(p2f, p2fb));
                    pedges.push_back(Edge(p2fb, p2));
                    printf("p2fb at %lf\t%lf\n", p2fb.x, p2fb.y);
                }
                else pedges.push_back(Edge(p2f,p2));
            
                if(p4fflag){
                    pedges.push_back(Edge(p2, p4fb));
                    pedges.push_back(Edge(p4fb, p4f));
                    printf("p4fb at %lf\t%lf\n", p4fb.x, p4fb.y);
                }
                else pedges.push_back(Edge(p2,p4f));
            }
            else pedges.push_back(Edge(p2f, p4f));

            pedges.push_back(Edge(p4f,p3f));

            if(! p1_is_end){
                if(p3fflag){
                    pedges.push_back(Edge(p3f, p3fb));
                    pedges.push_back(Edge(p3fb, p1));
                    printf("p3fb at %lf\t%lf\n", p3fb.x, p3fb.y);
                }
                else pedges.push_back(Edge(p3f,p1));
                if(p1fflag){
                    pedges.push_back(Edge(p1, p1fb));
                    pedges.push_back(Edge(p1fb, p1f));
                    printf("p1fb at %lf\t%lf\n", p1fb.x, p1fb.y);
                }
                else pedges.push_back(Edge(p1,p1f));
            }
            else pedges.push_back(Edge(p3f, p1f));

            pdata.push_back(PolyData(pedges));
        }
    }

    return replace;
}

vector<double> get_grn_displacements(int size, double sdev){
    vector<double> displacements;
    const gsl_rng_type *T;
    gsl_rng *r;
    
    T = gsl_rng_mt19937;
    r = gsl_rng_alloc(T);

    for(int i = 0; i < size; i++){
        displacements.push_back(gsl_ran_gaussian(r, sdev));
    }

    gsl_rng_free(r);

    return displacements;
}

/*
This function displaces nodes in a network according to specified displacements.
The function takes as arguments the original set of points and edges describing
the network, and a pair of displacements for each node in the network. Points 
are shifted, and edges are updated accordingly.
*/
void displace_points_grn(vector<Edge>& edges, vector<Point> &skel_points, vector<double> displacements){
    //Map from old to new point locations
    unordered_map<Point, Point> change_map;
    int index = 0, iter;
    double xdisp, ydisp;
    Point next, replace, newp1, newp2;
    Edge old_edge;

    cout << "Size here: " << edges.size() << "\n";

    for(iter = 0; iter < skel_points.size(); iter++){
        next = skel_points[iter];
        xdisp = displacements[index];
        ydisp = displacements[index + 1];
        replace = Point(next.x + xdisp, next.y + ydisp);
        change_map.insert(make_pair(next, replace));
        skel_points[iter] = replace;
        index += 2;
    }

    for(int iter = 0; iter < edges.size(); iter++){
        old_edge = edges[iter];
        newp1 = change_map[old_edge.p1];
        newp2 = change_map[old_edge.p2];
        edges[iter] = Edge(newp1, newp2);
    }

    cout << "Finished displacing.\n";
}

double mag(vector<double> vec){
    double mag = 0;
    int index;

    for(index = 0; index < vec.size(); index++){
        mag += vec[index] * vec[index];
    }

    return sqrt(mag);
}

int run_fire_md(vector<double>& pos, vector<double>& vel, map<int, vector<EdgeDatum>> neighbor_map, vector<int> triples, vector<double> ref_angles, unordered_set<int> pinned, double ks, double kb, double mass, double fcut, int max_steps, int nmin, double fire_params[5], bool& flag){

    vector<double> forces;
    double alpha, finc, fdec, alpha_start, falpha, fire_dt, fire_dt_max;
    int num_points, num_dof, step_count = 0, since_leq_0 = 0, index, freq;
    int rep_freq, digits, rep_count = 1;
    double power, dt, dt_sq, vmag, fmag, inv_fmag, sqrt_inner_dof;
    string base, full_name;
    bool report;

    //Find the number of points and initialize force vector
    num_dof = pos.size();
    num_points = num_dof / 2;
    sqrt_inner_dof = sqrt(num_dof - pinned.size());

    //Unpack parameters of FIRE minimization scheme and initialize values
    alpha_start = fire_params[0];
    falpha = fire_params[1];
    fire_dt_max = fire_params[2];
    finc = fire_params[3];
    fdec = fire_params[4];

    alpha = alpha_start;
    dt = fire_dt_max;
    dt_sq = dt*dt;
   
    //Find the forces at the outset
    forces.assign(num_dof, 0);
    get_forces(neighbor_map, pos, forces, triples, ref_angles, pinned, ks, kb);
    vmag = mag(vel);

    cout << "Starting RMS Force: " << mag(forces) / sqrt_inner_dof << "\n";
 
    //Perform molecular dynamics steps using velocity verlet method with FIRE
    //until the kinetic energy cutoff is reached, or the maximum number of 
    //steps have taken place
    while(step_count < max_steps){
        step_count ++;        

        //Update positions
        for(index = 0; index < num_dof; index++){
            if(pinned.find(index/2) == pinned.end()){
                pos[index] += dt*vel[index] + .5 * forces[index]/mass * dt_sq;
                vel[index] += .5 * dt * forces[index] / mass;
            }
        }

        //Calculate forces
        get_forces(neighbor_map, pos, forces, triples, ref_angles,pinned,ks,kb);

        //Update velocities and calculate power
        power = 0;
        for(index = 0; index < num_dof; index++){
            if(pinned.find(index/2) == pinned.end()){
                vel[index] += .5 * dt * forces[index] / mass;
                power += vel[index] * forces[index];
            }
        }

        //Adjust velocities according to FIRE algorithm
        vmag = mag(vel);
        fmag = mag(forces);
        inv_fmag = 1 / fmag;

        for(index = 0; index < num_dof; index++){
            if(pinned.find(index) == pinned.end()){
                vel[index] += alpha*(vmag*forces[index]*inv_fmag - vel[index]);
            }
        }
        

        //Adjust FIRE parameters according to current power
        if(power > 0){
            since_leq_0 ++;
            if(since_leq_0 > nmin){
                dt = min(dt*finc, fire_dt_max);
                dt_sq = dt*dt;
                alpha *= falpha;
            }
        }
        else{
            since_leq_0 = 0;
            dt *= fdec;
            dt_sq = dt * dt;
            alpha = alpha_start;
            vel.assign(num_dof, 0);
        }

        /*if(report && step_count % rep_freq == 0){
            full_name = report_name(base, digits, rep_count);
            report_deformed(neighbor_map, pos, full_name);
            rep_count ++;
        }*/

        //Check for kinetic energy convergence
        if(fmag / sqrt_inner_dof < fcut){
            flag = true;
            cout << "RMS Force: " << fmag / sqrt_inner_dof << "\n";
            return step_count;
        }
    }

    flag = false;
    cout << "RMS Force: " << fmag / sqrt_inner_dof << "\n";
    cout << "Ending energy: " << .5 * mass * vmag * vmag + get_pe(neighbor_map, pos, triples, ref_angles, ks, kb) << "\n";
    return step_count;
}

void affine_displace_run(){
    vector<Edge> skeleton, envelope;
    set<Point> pset;
    vector<Point> skel_points;
    vector<vector<int>> edge_map;
    ifstream skel_stream, small_stream;
    map<int, vector<EdgeDatum>> neighbor_map;
    vector<int> triples;
    vector<double> coords, displacements, vel, ref_angles;
    double sdev, thickness, miny, maxy, mass, fcut, ks, kb;
    map<int, vector<double>> ref_coords;
    int n_steps, num_read, step, max_md_steps, nmin, iter;
    string response;
    vector<PolyData> pdata;
    vector<vector<double>> inv_mats;
    unordered_set<int> pinned;
    bool converge_flag;
    double fparams[5];
    ofstream env_file;

    //Prompt for displacement standard deviation and displacement steps
    do{
        cout << "Enter displacement standard deviation: ";
        getline(cin, response);
        num_read = sscanf(response.c_str(), "%lf", &sdev);
    }while(num_read < 1 || sdev <= 0);

    do{
        cout << "Enter the number of displacement steps: ";
        getline(cin, response);
        num_read = sscanf(response.c_str(), "%d", &n_steps);
    }while(num_read < 1 || n_steps < 1);

    //Obtain envelope thickness and create envelope
    do{
        cout << "Enter envelope thickness: ";
        getline(cin, response);
        num_read = sscanf(response.c_str(), "%lf", &thickness);
    }while(num_read < 1 || thickness <= 0);

    do{
        cout << "Enter the stretching and bending stiffnesses: ";
        getline(cin, response);
        num_read = sscanf(response.c_str(), "%lf %lf", &ks, &kb);
        if(ks <= 0) cout << "Enter a positive stretching stiffness.\n";
        if(kb < 0) cout << "Enter a non-negative bending stiffness.\n";
    }while(num_read < 1 || ks <= 0 || kb < 0);

    //Obtain FIRE parameters for adjusting velocity
    do{
        cout << "Enter alpha, falpha, dt, finc, fdec, and nmin: ";
        getline(cin, response);
        num_read = sscanf(response.c_str(),"%lf %lf %lf %lf %lf %d", fparams, fparams+1, fparams+2, fparams+3, fparams+4, &nmin);
    }while(num_read < 6);

    do{
        cout << "Enter particle mass, cutoff force, and maximum md steps: ";
        getline(cin, response);
        num_read = sscanf(response.c_str(), "%lf %lf %d", &mass, &fcut, &max_md_steps);
    }while(num_read < 3);

    //Prompt for large-scale skeleton
    open_dat_file("Enter the skeleton file: ", skel_stream);
    skeleton = read_edges(skel_stream);
    pset = get_points(skeleton);
    skel_points.insert(skel_points.begin(), pset.begin(), pset.end());
    pset.clear();

    //Prompt for unstrained small-scale network
    open_dat_file("Enter the small-scale file: ", small_stream);
    if(kb > 0){
        read_network(small_stream, coords, neighbor_map, triples, ref_angles, true);
    }
    else{
        read_network(small_stream, coords, neighbor_map, triples, ref_angles, false);
    }

    envelope = add_thickness(skeleton, thickness, pdata, false, false);

    cout << "Original size: " << envelope.size() << "\n";

    //Make map from envelope edges to small-scale vertices
    on_edge(envelope, coords, edge_map);

    //Log starting coordinates of pinned small-scale vertices
    for(vector<int> index_list : edge_map){
        for(int i : index_list){
            double ref_coord[2] = {coords[2 * i], coords[2 * i + 1]};
            ref_coords.insert(make_pair(i, vector<double>(ref_coord, ref_coord + 2)));
            pinned.insert(i);
        }
    }

    //Find inverse matrices for calculating affine transformations
    inv_mats = get_inv_mats(envelope);

    //Obtain displacements for each skeleton node, then progressively deform
    //skeleton and large-scale envelope
    displacements = get_grn_displacements(2 * skel_points.size(), sdev);
    for(iter = 0; iter < displacements.size(); iter ++){
        displacements[iter] = displacements[iter] / n_steps;
    }

    vel.assign(coords.size(), 0);

    for(step = 1; step <= n_steps; step++){
        displace_points_grn(skeleton, skel_points, displacements);
        envelope = add_thickness(skeleton, thickness, pdata, false, false);

        cout << "Check 1.\n";
        get_deformed_vertices(envelope, inv_mats, edge_map, ref_coords, coords);
        cout << "Check 2.\n";
        run_fire_md(coords, vel, neighbor_map, triples, ref_angles, pinned, ks, kb, mass, fcut, max_md_steps, nmin, fparams, converge_flag);
        cout << "Check 3.\n";
        if(! converge_flag) cerr << "Convergence not reached.\n";
    }

    if(yesno("Add small-scale displacements?")){
        while(true){
            cout << "Enter the standard deviation: ";
            getline(cin, response);
            num_read = sscanf(response.c_str(), "%lf", &sdev);
            if(num_read == 1 && sdev > 0) break;
            else cerr << "Enter a positive number.\n";
        }

    }

    cout << "Enter a file name for the final network, or enter to decline: ";
    getline(cin, response);

    if(! response.compare("") == 0){
        report_deformed(neighbor_map, coords, response);
    }

    cout << "Enter a file name for the final skeleton, or enter to decline: ";
    getline(cin, response);
    if(! response.compare("") == 0){
        env_file.open(response);
        for(Edge env_e : envelope){
            env_file << env_e.p1.x << " " << env_e.p1.y << "\n";
            env_file << env_e.p2.x << " " << env_e.p2.y << "\n\n";
        }
        env_file.close();
    }
}

int main(int argc, char **argv){
    FILE *ranfile;
    unsigned seed;
    char num_string[11];

    //Initialize parameters for random number generation
    ranfile = fopen("/dev/urandom", "r");
    fread(&seed, sizeof(unsigned), 1, ranfile);
    fclose(ranfile);
    sprintf(num_string, "%u", seed);
    setenv("GSL_RNG_SEED", num_string, 1);
    gsl_rng_env_setup();

    //Until the user indicates otherwise, generate distorted networks
    while(yesno("Perform a displacement run?")){
        affine_displace_run();
        cout << "Still working.\n";
    }
}
