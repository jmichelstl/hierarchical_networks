/*This program produces a hierarchical network by first producing a large-scale
*envelope, and then triangulating that envelope using Triangle, by Jonathan
*Richard Shewchuck. The program acts as an intermediary between the user and
*Triangle, offering the option to set the minimum internal angle of any
*triangle, the maximum allowable area of a triangle, and whether or not the
*to enforce the condition that triangles obey the Delaunay condition. After the
*small-scale triangulation is made, the program will dilute the small scale,
*if desired, and report a large-scale tiling.
*/

#define REAL double
#define VOID int

#include <iostream>
#include <sstream>
#include "network_utils.h"
#include "triangle.hpp"
#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <map>
#include <iomanip>
#include <algorithm>
#include <list>
#include <string.h>
#include <unordered_set>
#include "lloyd_move.h"

using namespace std;

#define tuple_2d tuple<double, double>
#define cos60 0.5
#define sin60 0.866025404
#define get_gindex(x, y, nx) x + nx*y

enum class smallMode : char {
    grain = 'g',
    triangulate  = 't',
    random = 'r',
    import = 'i',
};

enum class largeMode : char {
    crystalline  = 'c',
    poisson = 'p',
    lloyd = 'l',
};

void scale_vector(vector<double>& in, double scale){
    int index;
    for(index = 0; index < in.size(); index++){
        in.at(index) = in.at(index)*scale;
    }
}

vector<Edge> import_edges(ifstream &input){
    int num_read, point_count = 0;
    string nextline;
    double x, y;
    vector<Edge> edges;
    Point p1, p2;

    while(! input.eof()){
        getline(input, nextline);

        num_read = sscanf(nextline.c_str(), "%lf %lf", &x, &y);
        if(num_read >= 2){
            if(point_count == 0){
                p1 = Point(x, y);
            }

            else{
                p2 = Point(x, y);
            }
            point_count ++;
        }

        if(point_count == 2){
            point_count = 0;
            edges.push_back(Edge(p1, p2));
        }
    }

    input.close();

    return edges;
}

bool import_lattice(vector<vector<double>>& rules, map<int,vector<vector<double>>>& nns, double scale){
    ifstream latfile;
    string name, nextline;
    bool again, blank, fileopen = false;
    vector<double> rule, nn;
    int lcount = 0, pointiter, nncount = 0, ruleiter = 0;

    //Prompt for file name
    do{
        printf("Enter the lattice file name: ");
        getline(cin, nextline);
        name = split(nextline, ' ')[0];

        if(!name.empty()) latfile.open(name);

        if(!latfile.is_open()){
            again = yesno("No file was read. Try again? ");
            if(!again) return false;
        }
    }while(!latfile.is_open());

    //Process rules until a blank line is reached
    blank = false;
    while(!latfile.eof()){
        lcount ++;
        getline(latfile, nextline);
        if(nextline.empty()) break;

        rule = parse_doubles(split(nextline, ' '));
        if(rule.size() < 3){
            fprintf(stderr, "Too few numbers in rule on line line %d\n",lcount);
            latfile.close();
            return false;
        }
        scale_vector(rule, scale);
        rules.push_back(rule);
    }
    
    //Read nearest neighbor rules
    for(vector<double> nextrule : rules){
       ruleiter ++;
       for(pointiter = 1; pointiter <= nextrule.size() - 2; pointiter ++){
           vector<vector<double>> nextset;
           while(!latfile.eof()){
               lcount ++;
               getline(latfile, nextline);
               if(nextline.empty()) break;
               nn = parse_doubles(split(nextline, ' '));
               if(nn.size() != 2) fprintf(stderr, "Insufficient information for nearest neighbor rule on line %d.\n", lcount);
               else{
                   scale_vector(nn, scale);
                   nextset.push_back(nn);
               }
           }
           nns.insert(pair<int, vector<vector<double>>>(nncount++, nextset));
       }
    }


    latfile.close();
    return true;
}

/*
This function adds Gaussian random noise to the location of each point.
The function takes as arguments the original set of points and edges describing
the network, and a standard deviation for Gaussian random noise. Points are 
shifted, and edges are updated accordingly.
*/
void displace_points_grn(vector<Edge>& edges, double sdev){
    //Map from old to new point locations
    unordered_map<Point, Point> change_map;

    //Structures to generate Gaussian random noise
    const gsl_rng_type *T;
    gsl_rng *r;

    double xdisp, ydisp;
    set<Point> old_points;
    Point replace, newp1, newp2;
    Edge old_edge, new_edge;

    //gsl_rng_env_setup();
    T = gsl_rng_mt19937;
    r = gsl_rng_alloc(T);

    old_points = get_points(edges);

    for(Point next : old_points){
        xdisp = gsl_ran_gaussian(r, sdev);
        ydisp = gsl_ran_gaussian(r, sdev);
        replace = Point(next.x + xdisp, next.y + ydisp);
        change_map.insert(make_pair(next, replace));
    }

    for(int iter = 0; iter < edges.size(); iter++){
        old_edge = edges[iter];
        newp1 = change_map[old_edge.p1];
        newp2 = change_map[old_edge.p2];
        edges[iter] = Edge(newp1, newp2);
    }

    gsl_rng_free(r);
}


vector<Edge> makeedges(vector<vector<double>> rules, map<int,vector<vector<double>>> nns, vector<double> bnds){
    double x, y, x2, y2;
    int ruleiter = 0, rulesize, xiter, nextindex = 0;
    vector<double> rule;
    int base = 0, index, numrules = 0;
    Point p1, p2;
    set<Point> points;
    set<Point>::iterator piter;
    vector<Edge> edges;

    for(vector<double> vec : rules){
        numrules += vec.size() - 2;
    }

    y = bnds[1];
    while(y <= bnds[3]){
        xiter = 0;
        rule = rules[ruleiter%rules.size()];
        rulesize = rule.size() - 2;
        x = rule[0] + bnds[0];

        while(x <= bnds[2]){
           p1 = Point(x,y);
           piter = points.find(p1);
           if(piter != points.end()){
               p1.x = piter->x;
               p1.y = piter->y;
           }
           else points.insert(Point(p1.x, p1.y, nextindex++));

           index = base + xiter % rulesize;
           for(vector<double> nextnn : nns[index]){
               x2 = x + nextnn[0];
               y2 = y + nextnn[1];
               if(x2>=bnds[0] && x2<=bnds[2] && y2>=bnds[1] && y2<=bnds[3]){

                   p2 = Point(x2, y2);
                   piter = points.find(p2);
                   if(piter != points.end()){
                       p2.x = piter->x;
                       p2.y = piter->y;
                   }
                   else points.insert(Point(x2,y2,nextindex++));
                   edges.push_back(Edge(p1, p2));
               }
           }
           x += rule[xiter++%rulesize + 2];
        }
        y += rule[1];
        base = (base + rulesize) % numrules;
        ruleiter++;
    }

    return edges;
}

vector<Point> poisson_packing(vector<double> bounds, double radius, int max_attempts){

    vector<Point> poisson_points;
    vector<int> active_list;
    double spacing = radius / sqrt(2);
    double xrange = bounds[2] - bounds[0], yrange = bounds[3] - bounds[1];
    double minx, miny, maxx, maxy, currX, currY, newX, newY, r, theta;
    double dist_sq, r_sq = radius * radius;
    int xsteps = (int) (xrange/spacing)+1, ysteps = (int) (yrange/spacing)+1;
    int xIndex, minXIndex, maxXIndex, yIndex, minYIndex, maxYIndex, gindex;
    int ngindex, xIter, yIter, list_index, point_index, attempts;
    vector<int> grid = vector<int>(xsteps * ysteps, -1);
    const gsl_rng_type *T = gsl_rng_mt19937;
    gsl_rng *rng = gsl_rng_alloc(T);
    unsigned long int randex;
    bool point_added, conflict;
    Point neighbor;

    //Make the first point, and put it in the list
    currX = gsl_rng_uniform(rng) * xrange + bounds[0];
    currY = gsl_rng_uniform(rng) * yrange + bounds[1];
    poisson_points.push_back(Point(currX, currY));
    active_list.push_back(0);
    xIndex = (int) ((currX - bounds[0]) / spacing);
    yIndex = (int) ((currY - bounds[1]) / spacing);
    gindex = get_gindex(xIndex, yIndex, xsteps);
    grid[gindex] = 0;

    do{
        attempts = 0;
        point_added = false;
        list_index = (int) gsl_rng_uniform_int(rng,active_list.size());
        point_index = active_list[list_index];
        currX = poisson_points[point_index].x;
        currY = poisson_points[point_index].y;

        do{
            do{
                r = radius * gsl_rng_uniform(rng) + radius;
                theta = gsl_rng_uniform(rng) * 2 * M_PI;
                newX = currX + r * cos(theta);
                newY = currY + r * sin(theta);
            }while(!(newX>=bounds[0] && newY>=bounds[1] && newX<=bounds[2] && newY<=bounds[3]));
            xIndex = (int) ((newX - bounds[0]) / spacing);
            yIndex = (int) ((newY - bounds[1]) / spacing);
            gindex = get_gindex(xIndex, yIndex, xsteps);
            if(grid[gindex] >= 0) continue;
            
            if(xIndex > 2) minXIndex = xIndex - 3;
            else if(xIndex > 1) minXIndex = xIndex - 2;
            else if(xIndex > 0) minXIndex = xIndex - 1;
	    else minXIndex = xIndex;
            if(yIndex > 2) minYIndex = yIndex - 3;
            else if(yIndex > 1) minYIndex = yIndex - 2;
            else if(yIndex > 0) minYIndex = yIndex - 1;
            else minYIndex = yIndex;

            if(xIndex < xsteps - 3) maxXIndex = xIndex + 3;
            else if(xIndex < xsteps - 2) maxXIndex = xIndex + 2;
            else if(xIndex < xsteps - 1) maxXIndex = xIndex + 1;
            else maxXIndex = xIndex;
            if(yIndex < ysteps - 3) maxYIndex = yIndex + 3;
            else if(yIndex < ysteps - 2) maxYIndex = yIndex + 2;
            else if(yIndex < ysteps - 1) maxYIndex = yIndex + 1;
            else maxYIndex = yIndex;

            conflict = false;
            yIter = minYIndex;
            do{
                xIter = minXIndex;
                do{
                    ngindex = get_gindex(xIter, yIter, xsteps);
                    if(grid[ngindex] >= 0){
                        neighbor = poisson_points[grid[ngindex]];
                        dist_sq = (neighbor.x - newX) * (neighbor.x - newX);
                        dist_sq += (neighbor.y - newY) * (neighbor.y - newY);
                        if(dist_sq < r_sq) conflict  = true;
                    }
                    xIter ++;
                }while(xIter <= maxXIndex && ! conflict);
                yIter ++;
            }while(yIter <= maxYIndex && ! conflict);

            if(! conflict){
                poisson_points.push_back(Point(newX, newY));
                active_list.push_back(poisson_points.size() - 1);
                grid[gindex] = poisson_points.size() - 1;
                point_added = true;
            }
            attempts ++;
        }while(! point_added && attempts < max_attempts);
        if(! point_added){
            active_list.erase(active_list.begin()+list_index);
        }
    }while(! active_list.empty());

    gsl_rng_free(rng);
    return poisson_points;
}

vector<Edge> trim_to_size(vector<Edge> old, vector<double> bounds){

    Point p1, p2;
    bool p1_in, p2_in;
    vector<Edge> replace;

    for(Edge next : old){
        p1 = next.p1;
        p2 = next.p2;

        p1_in = p1.x >= bounds[0] && p1.y >= bounds[1] && p1.x <= bounds[2] && p1.y <= bounds[3];
        p2_in = p2.x >= bounds[0] && p2.y >= bounds[1] && p2.x <= bounds[2] && p2.y <= bounds[3];

        if(p1_in || p2_in) replace.push_back(next);
    }

    return replace;
}

vector<Edge> make_poisson_edges(vector<double> bounds, double radius){

    vector<Point> poisson_points;
    triangulateio in, out;
    vector<Edge> poisson_edges;
    int iter, index1, index2;
    Point p1, p2;
    char options[5];
    vector<double> wide_bounds(bounds);

    wide_bounds[0] -= 2*radius;
    wide_bounds[1] -= 2*radius;
    wide_bounds[2] += 2*radius;
    wide_bounds[3] += 2*radius;

    sprintf(options, "%s", "zenQ");

    poisson_points = poisson_packing(wide_bounds, radius, 30);

    in.numberofpoints = poisson_points.size();
    in.pointlist = (REAL *) malloc(poisson_points.size() * 2 * sizeof(REAL));
    for(iter = 0; iter < poisson_points.size(); iter ++){
        in.pointlist[iter * 2] = poisson_points[iter].x;
        in.pointlist[iter * 2 + 1] = poisson_points[iter].y;
    }
    in.numberofpointattributes = 0;
    in.pointmarkerlist = (int *) NULL;
    in.numberofsegments = 0;
    in.numberofholes = 0;
    in.numberofregions = 0;

    out.pointlist = (REAL *) NULL;
    out.numberofpointattributes = 0;
    out.numberoftriangleattributes = 0;
    out.pointmarkerlist = (int *) NULL;
    out.trianglelist = (int *) NULL;
    out.neighborlist = (int *) NULL;
    out.edgelist = (int *) NULL;
    out.edgemarkerlist = (int *) NULL;

    triangulate(options, &in, &out, (struct triangulateio *) NULL);

    for(iter = 0; iter < out.numberofedges; iter++){
        index1 = out.edgelist[iter * 2];
        index2 = out.edgelist[iter * 2 + 1];
        p1 = Point(out.pointlist[2*index1], out.pointlist[2*index1+1]);
        p2 = Point(out.pointlist[2*index2], out.pointlist[2*index2+1]);
        poisson_edges.push_back(Edge(p1, p2));
    }

    free(in.pointlist);
    free(out.pointlist);
    free(out.pointmarkerlist);
    free(out.trianglelist);
    free(out.neighborlist);
    free(out.edgelist);
    free(out.edgemarkerlist);

    return trim_to_size(poisson_edges, bounds);
}

void update_plist(set<Point> &points, Point &p){
    auto iter = points.find(p);
    if(iter != points.end()){
        p.x = iter->x;
        p.y = iter->y;
    }
    else points.insert(p);
}

vector<Edge> generate_tri_lat(double minx, double miny, double maxx, double maxy, double scale){

    double offset = 0, x1, y1, x2, y2, dx = scale*cos60, dy = scale*sin60;
    vector<Edge> lat_edges;
    set<Point> lat_points;
    Point p1, p2;

    for(y1 = miny; y1 <= maxy + FLOAT_TOL; y1 += dy){
        for(x1 = minx + offset; x1 <= maxx + FLOAT_TOL; x1 += scale){
            p1 = Point(x1, y1);
            update_plist(lat_points, p1);
            x2 = x1 - dx;
            y2 = y1 + dy;
            if(x2 >= minx && y2 <= maxy + FLOAT_TOL){
                p2 = Point(x2,y2);
                update_plist(lat_points, p2);
                lat_edges.push_back(Edge(p1, p2));
            }
            x2 = x1 + dx;
            if(x2 <= maxx + FLOAT_TOL && y2 <= maxy + FLOAT_TOL){
                p2 = Point(x2,y2);
                update_plist(lat_points, p2);
                lat_edges.push_back(Edge(p1, p2));
            }
            x2 = x1 + scale;
            if(x2 <= maxx + FLOAT_TOL){
                p2 = Point(x2,y1);
                update_plist(lat_points, p2);
                lat_edges.push_back(Edge(p1, p2));
            }
        }
        offset = offset == 0 ? scale / 2 : 0;
    }

    return lat_edges;
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
        cerr<< "Low: " << low << " High: " << high << "\n";
        return;
    }

    dx = hwidth * (cos(low)/tan(diff) - sin(low));
    dy = hwidth * (sin(low)/tan(diff) + cos(low));
}

vector<Edge> add_thickness(vector<Edge> current, double thickness, vector<PolyData>& pdata, bool get_envelope){
    map<Point, vector<double>> angmap;
    vector<Edge> replace;
    Point p1, p2, key;
    Point p1f, p2f, p3f, p4f, p1fb, p2fb, p3fb, p4fb;
    double low, high, ang1, ang2, dx, dy, hwidth, slope;
    double ymin = FLT_MAX, ymax = FLT_MIN, ylow, yhigh, y_ext;
    bool p1_is_end, p2_is_end;
    vector<double> angvec;

    hwidth = thickness/2;

    for(Edge e : current){
        p1 = e.p1;
        p2 = e.p2;

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

    for(Edge e : current){
        p1 = e.p1;
        p2 = e.p2;
        ang1 = atan2(p2.y-p1.y,p2.x-p1.x);
        if(ang1 < 0) ang1 += 2*M_PI;
        ang2 = ang1 < M_PI ? ang1 + M_PI : ang1 - M_PI;
        slope = (p2.y - p1.y)/(p2.x - p1.x); 
        p1_is_end = false;
        p2_is_end = false;

        getangles(angmap[p1], ang1, low, high);
        if(ang1 == low){
            p1_is_end = true;
            p1f = Point(p1.x + hwidth*sin(ang1), p1.y - hwidth*cos(ang1));
            p3f = Point(p1.x - hwidth*sin(ang1), p1.y + hwidth*cos(ang1));
            
            if(get_envelope) replace.push_back(Edge(p1f,p3f));
        }
        else{
            changes(low, ang1, hwidth, dx, dy);
            p1f = Point(p1.x + dx, p1.y + dy);
            changes(ang1, high, hwidth, dx, dy);
            p3f = Point(p1.x + dx, p1.y + dy);
        }
 
        getangles(angmap[p2], ang2, low, high);
        if(ang2 == low){
            p2_is_end = true;
            p2f = Point(p2.x + hwidth*sin(ang1), p2.y - hwidth*cos(ang1));
            p4f = Point(p2.x - hwidth*sin(ang1), p2.y + hwidth*cos(ang1));
            
            if(get_envelope) replace.push_back(Edge(p2f,p4f));
        }
        else{
            changes(ang2, high, hwidth, dx, dy);
            p2f = Point(p2.x + dx, p2.y + dy);
            changes(low, ang2, hwidth, dx, dy);
            p4f = Point(p2.x + dx, p2.y + dy);
        }

        if(get_envelope) replace.push_back(Edge(p1f,p2f));
        if(get_envelope) replace.push_back(Edge(p3f,p4f));

        vector<Edge> pedges;
        pedges.push_back(Edge(p1f,p2f));

        if(!p2_is_end){
            pedges.push_back(Edge(p2f,p2));
            pedges.push_back(Edge(p2,p4f));
        }
        else pedges.push_back(Edge(p2f, p4f));

        pedges.push_back(Edge(p4f,p3f));

        if(! p1_is_end){
            pedges.push_back(Edge(p3f,p1));
            pedges.push_back(Edge(p1,p1f));
        }
        else pedges.push_back(Edge(p3f, p1f));

        pdata.push_back(PolyData(pedges));
    }

    return replace;
}

void shuffle_edges(vector<Edge>& edges){
    Edge temp;
    int i;
    unsigned long int randval;
    int len = edges.size();
    const gsl_rng_type *T;
    gsl_rng *r;

    T = gsl_rng_mt19937;
    r = gsl_rng_alloc(T);

    for(i = 0; i < len; i ++){
        randval = gsl_rng_uniform_int(r, len);
        temp = edges[i];
        edges[i] = edges[randval];
        edges[randval] = temp;
    }

    gsl_rng_free(r);
}

vector<Edge> randomMST(vector<Edge> in, vector<Edge>& rejects){
    unordered_map<Point,int> point_map;
    vector<Edge> keep;
    int index = 0, numkept = 0, numneeded;
    int root1, root2;
    Edge next;

    for(Edge next : in){
        if(point_map.find(next.p1) == point_map.end()){
            point_map.insert(make_pair(next.p1, index++));
        }
        if(point_map.find(next.p2) == point_map.end()){
            point_map.insert(make_pair(next.p2, index++));
        }
    }

    vector<int> union_find_table(point_map.size(), -1);
    numneeded = point_map.size() - 1;

    shuffle_edges(in);
    while(!in.empty()){
        next = in[0];
        in.erase(in.begin());
        root1 = find_root(union_find_table, point_map[next.p1]);
        root2 = find_root(union_find_table, point_map[next.p2]);
        if(root1 != root2){
            makeunion(union_find_table, root1, root2);
            keep.push_back(next);
            if(++numkept == numneeded) break;
        }
        else rejects.push_back(next);
    }

    for(Edge nextEdge : in){
        rejects.push_back(nextEdge);
    }

    return keep;
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

vector<Edge> random_connected(vector<Edge> in, double toKeep){
    int num_needed, reject_index;
    vector<Edge> kept, rejects;

    num_needed = (int) (toKeep * in.size());
    kept = randomMST(in, rejects);
    shuffle_edges(rejects);

    reject_index = 0;
    do{
        kept.push_back(rejects[reject_index++]);
    }while(reject_index < rejects.size() && kept.size() < num_needed);

    return kept;
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

void create_input(triangulateio *in, int point_count, vector<vector<Point>> polygons, vector<Point> border, vector<tuple_2d> holes){
    int pindex = 0, plabel = 0, seg_index = 0, hole_index = 0;

    in->numberofpoints = point_count;
    in->numberofsegments = point_count;
    in->pointlist = (double *) malloc(point_count * 2 * sizeof(double));
    in->segmentlist = (int *) malloc(point_count * 2 * sizeof(int));
    in->numberofregions = 0;

    for(vector<Point> nextpoly : polygons){
        for(Point nextPoint : nextpoly){
            in->pointlist[pindex] = nextPoint.x;
            pindex++;
            in->pointlist[pindex] = nextPoint.y;
            pindex++;
        }
        add_segments(in->segmentlist, nextpoly.size(), plabel, seg_index);
    }
    for(Point nextPoint : border){
        in->pointlist[pindex] = nextPoint.x;
        pindex++;
        in->pointlist[pindex] = nextPoint.y;
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

void create_grain_input(triangulateio *in, vector<Point> polygon){
    int pindex = 0, plabel = 0, seg_index = 0, point_count = polygon.size();

    in->numberofpoints = point_count;
    in->numberofsegments = point_count;
    in->pointlist = (double *) malloc(point_count * 2 * sizeof(double));
    in->segmentlist = (int *) malloc(point_count * 2 * sizeof(int));
    in->numberofregions = 0;

    for(Point nextPoint : polygon){
        in->pointlist[pindex] = nextPoint.x;
        pindex++;
        in->pointlist[pindex] = nextPoint.y;
        pindex++;
    }
    add_segments(in->segmentlist, polygon.size(), plabel, seg_index);

    in->segmentmarkerlist = (int *) malloc(sizeof(int) * in->numberofsegments);
    in->numberofpointattributes = 0;
    in->pointmarkerlist = (int *) NULL;
    in->numberofedges = 0;
    in->edgelist = (int *) NULL;
    in->numberofholes = 0;
    in->holelist = (double *) NULL;
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

//Don't bother coallescing points. Just translate indices into points, and
//find the border
void get_polygons_simple(map<int, Point> imap, vector<vector<int>> cycle_list, vector<vector<Point>> &polygons, vector<Point> &border){
    double miny = FLT_MAX;
    int iter, min_cycle = 0;

    for(iter = 0; iter < cycle_list.size(); iter++){
        vector<Point> next_poly;
        for(int next_vertex : cycle_list[iter]){
            next_poly.push_back(imap[next_vertex]);
            if(imap[next_vertex].y < miny){
                miny = imap[next_vertex].y;
                min_cycle = iter;
            }
        }
        polygons.push_back(next_poly);
    }

    border = polygons[min_cycle];
    polygons.erase(polygons.begin() + min_cycle, polygons.begin() + min_cycle + 1);
}

vector<vector<Point>> get_grain_polygons(map<int, Point> imap, vector<vector<int>> cycle_list){

    double miny = FLT_MAX;
    int iter;
    vector<vector<Point>> polygons;

    for(iter = 0; iter < cycle_list.size(); iter++){
        vector<Point> next_poly;
        for(int next_vertex : cycle_list[iter]){
            next_poly.push_back(imap[next_vertex]);
        }
        polygons.push_back(next_poly);
    }

    return polygons;
}

bool in_poly(vector<Point> plist, Point p){
    bool inpoly = false;
    Point p1, p2;
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

void get_holes(vector<vector<Point>> polygons, vector<tuple_2d> &holes){
    double midx, midy, dx, dy;
    Point p1, p2, guess1, guess2;
    bool is1, is2;

    for(vector<Point> nextpoly : polygons){
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
        guess1 = Point(midx - dy, midy + dx);
        guess2 = Point(midx + dy, midy - dx);
        do{
            guess1 = Point(midx - dy, midy + dx);
            guess2 = Point(midx + dy, midy - dx);
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

void build_network_maps(vector<Edge> edges, map<int, Point> &imap, map<int, list<int>> &emap){

    map<Point, int> pmap;
    int pindex = 0, index1, index2;
    Point p1, p2;

    for(Edge next : edges){
        p1 = next.p1;
        p2 = next.p2;
        if(pmap.find(p1) == pmap.end()){
            pmap.insert(make_pair(p1, pindex));
            imap.insert(make_pair(pindex, p1));
            emap.insert(make_pair(pindex, list<int>()));
            pindex ++;
        }
        if(pmap.find(p2) == pmap.end()){
            pmap.insert(make_pair(p2, pindex));
            imap.insert(make_pair(pindex, p2));
            emap.insert(make_pair(pindex, list<int>()));
            pindex ++;
        }

        index1 = pmap[p1];
        index2 = pmap[p2];
        emap[index1].insert(emap[index1].end(), index2);
        emap[index2].insert(emap[index2].end(), index1);
    }
}

//Find the left and right borders of a lattice of edges
void get_borders(vector<Edge> all, vector<Point> &left, vector<Point> &right){

    map<double, double> min_map, max_map;
    Point p1, p2;

    for(Edge next : all){
        p1 = next.p1;
        p2 = next.p2;
        if(min_map.find(p1.y) == min_map.end()){
            min_map.insert(make_pair(p1.y, p1.x));
        }
        if(min_map.find(p2.y) == min_map.end()){
            min_map.insert(make_pair(p2.y, p2.x));
        }
        if(max_map.find(p1.y) == max_map.end()){
            max_map.insert(make_pair(p1.y, p1.x));
        }
        if(max_map.find(p2.y) == max_map.end()){
            max_map.insert(make_pair(p2.y, p2.x));
        }
        if(min_map[p1.y] > p1.x) min_map[p1.y] = p1.x;
        if(min_map[p2.y] > p2.x) min_map[p2.y] = p2.x;
        if(max_map[p1.y] < p1.x) max_map[p1.y] = p1.x;
        if(max_map[p2.y] < p2.x) max_map[p2.y] = p2.x;
    }

    for(auto iter = min_map.begin(); iter != min_map.end(); iter ++){
        left.push_back(Point(iter->second, iter->first));
    }
    for(auto iter = max_map.begin(); iter != max_map.end(); iter ++){
        right.push_back(Point(iter->second, iter->first));
    }
}

//Given a reference edge and a description of a polygon, find those edges of
//the polygon that are parallel to the reference edge
void get_parallel_edges(Edge ref, PolyData pdat, vector<Edge> &parallel_edges, vector<Edge> &other_edges){

    double ref_angle, curr_angle;

    ref_angle = atan2(ref.p2.x - ref.p1.x, ref.p2.y - ref.p1.y);

    for(Edge curr : pdat.polyEdges){
        curr_angle = atan2(curr.p2.x - curr.p1.x, curr.p2.y - curr.p1.y);
        if(abs(curr_angle - ref_angle) < FLOAT_TOL || abs(abs(curr_angle-ref_angle)-M_PI)<FLOAT_TOL){
            parallel_edges.push_back(curr);
        }
        else other_edges.push_back(curr);
    }
}

//Given a top and bottom edge for a grain, find the top and bottom y values
//and middle x values to use as the bounds for a rectangular region in which
//to form a lattice
void get_grain_boundaries(vector<Edge> top_bottom, double &minx, double &miny, double &maxx, double &maxy){

    double y1, y2;
    vector<double> xvals;

    y1 = top_bottom[0].p1.y;
    y2 = top_bottom[1].p1.y;
    miny = y1 < y2 ? y1 : y2;
    maxy = y1 > y2 ? y1 : y2;

    xvals.push_back(top_bottom[0].p1.x);
    xvals.push_back(top_bottom[0].p2.x);
    xvals.push_back(top_bottom[1].p1.x);
    xvals.push_back(top_bottom[1].p2.x);
    sort(xvals.begin(), xvals.end());
    minx = xvals[1];
    maxx = xvals[2];
}

//Remove edges from a grain that lie on the grain boundary. These will possibly
//be broken by the addition of Steiner points.
vector<Edge> remove_border(vector<Edge> lat_edges, vector<Edge> border){
    unordered_set<Edge> lat_set(lat_edges.begin(), lat_edges.end());
    unordered_set<Edge>::iterator pos;

    for(Edge bEdge : border){
        if((pos = lat_set.find(bEdge)) != lat_set.end()) lat_set.erase(pos);
    }

    return vector<Edge>(lat_set.begin(), lat_set.end());
}

//Given a lattice placed within a grain, find overhanging segments of that
//grain's bounding edges that extend beyond the lattice. These will be used
//to form the envelopes of grain boundaries, which will need to be triangulated
//separately
void get_overhang(vector<Edge> top_bottom, vector<Point> &left, vector<Point> &right, vector<Edge> &overhang, double cutoff, map<Point, Point> &change_map){

    Point bottom, top;
    int iter;
    Edge low, high;
    double minx, maxx;

    if(top_bottom[0].p1.y < top_bottom[1].p1.y){
       low = top_bottom[0];
       high = top_bottom[1];
    }
    else{
       low = top_bottom[1];
       high = top_bottom[0];
    }

    bottom = *(left.begin());
    top = *(left.crbegin());
    minx = low.p1.x < low.p2.x ? low.p1.x : low.p2.x;
    if(minx < bottom.x){
        if(bottom.x - minx < cutoff){
            change_map.insert(make_pair(bottom, Point(minx,low.p1.y)));
            left[0] = Point(minx,low.p1.y);
        }
        else overhang.push_back(Edge(Point(minx, low.p1.y),Point(bottom.x, low.p1.y)));
    }
    minx = high.p1.x < high.p2.x ? high.p1.x : high.p2.x;
    if(minx < top.x){
        if(top.x - minx < cutoff){
            change_map.insert(make_pair(top, Point(minx,high.p1.y)));
            left[left.size() - 1] = Point(minx, high.p1.y);
        }
        else overhang.push_back(Edge(Point(minx, high.p1.y),Point(top.x, high.p1.y)));
    }

    bottom = *(right.begin());
    top = *(right.crbegin());
    maxx = low.p1.x > low.p2.x ? low.p1.x : low.p2.x;
    if(bottom.x < maxx){
        if(maxx - bottom.x < cutoff){
            right[0] = Point(maxx, low.p1.y);
            change_map.insert(make_pair(bottom, Point(maxx, low.p1.y)));
        }
        else overhang.push_back(Edge(Point(bottom.x, low.p1.y), Point(maxx, low.p1.y)));
    }
    maxx = high.p1.x > high.p2.x ? high.p1.x : high.p2.x;
    if(top.x < maxx){
        if(maxx - top.x < cutoff){
            change_map.insert(make_pair(top, Point(maxx, high.p1.y)));
            right[right.size() - 1] = Point(maxx, high.p1.y);
        }
        else overhang.push_back(Edge(Point(top.x, high.p1.y), Point(maxx, high.p1.y)));
    }
}

void fix_edges(vector<Edge> &lattice, map<Point, Point> change_map){

    int iter;
    bool fix1, fix2;
    Point p1, p2;

    for(iter = 0; iter < lattice.size(); iter++){
        fix1 = false;
        fix2 = false;
        if(change_map.find(lattice[iter].p1) != change_map.end()){
            fix1 = true;
            p1 = change_map[lattice[iter].p1];
        }
        if(change_map.find(lattice[iter].p2) != change_map.end()){
            fix2 = true;
            p2 = change_map[lattice[iter].p2];
        }

        if(fix1 || fix2){
            if(! fix1) p1 = lattice[iter].p1;
            if(! fix2) p2 = lattice[iter].p2;
            lattice[iter] = Edge(p1,p2);
        }
    }
}

//Determine whether a polygon has an end cap
bool has_end_point(vector<Point> polygon, unordered_set<Point> end_points){

    for(Point next_point : polygon){
        if(end_points.find(next_point) != end_points.end()) return true;
    }

    return false;
}

void make_grain_based_network(vector<Edge> skeleton, vector<PolyData> pdata, double width, int num_rows, double min_ang, double max_area, bool enforce_delaunay, vector<Edge> &tri_edges, set<Point> &points){

    vector<Edge> junction_pool, raw_pool, parallel_edges, other_edges, lat_edges;
    vector<Edge> overhang;
    vector<Edge> unfixed;
    unordered_set<Edge> end_caps;
    unordered_set<Point> end_points;
    int poly_iter, piter, tri_iter, len;
    double angle, minx, miny, maxx, maxy, length;
    vector<Point> left, right;
    map<int, Point> imap, tri_point_map;
    map<int, list<int>> emap;
    vector<vector<int>> cycle_list;
    vector<vector<Point>> polygons;
    string base = "pQze";
    ostringstream tri_opts;
    char *options;
    triangulateio in, out;
    FILE *jreport, *env_report;
    string response;
    Point pivot;
    map<Point, Point> change_map;
    bool report_tri_envs;
    //unordered_set<Edge> boundary_set;
    Edge toAdd;

    length = width / (num_rows - 1) / sin60;

    //Generate a lattice for each grain. Add the lattice edges to a pool.
    //Also create edges at the periphery of each grain, from which to form
    //the outlines of junctions between adjacent grains.
    for(poly_iter = 0; poly_iter < pdata.size(); poly_iter ++){

        parallel_edges.clear();
        other_edges.clear();

        get_parallel_edges(skeleton[poly_iter], pdata[poly_iter], parallel_edges, other_edges);

        //If an edge of a polygon occurs only once, it is an end cap, because
        //it is not shared with any adjacent polygon. This should be noted to
        //avoid triangulation of small regions around end caps.
        for(Edge next_edge : other_edges){
            if(end_caps.find(next_edge) == end_caps.end()){
                end_caps.emplace(next_edge);
            }
            else end_caps.erase(next_edge);
        }

        //If two edges in a polygon are found that are parallel to the
        //accompanying skeleton edege, generate a small-scale lattice within
        //this polygon. Note segments of the polygon extending beyond this
        //small-scale lattice, as well as left and right edges of the 
        //lattice, in order to create polygons to be triangulated to stitch
        //together grains at junctions.
        pivot = parallel_edges[0].p1;
        angle = atan2(parallel_edges[0].p2.y - parallel_edges[0].p1.y,parallel_edges[0].p2.x-parallel_edges[0].p1.x);
        rotate_edges(parallel_edges, -angle, pivot);
        //cout << "Number of edges: " << parallel_edges.size() << "\n";
        //cout << "Alive here.\n";
        get_grain_boundaries(parallel_edges, minx, miny, maxx, maxy);
        lat_edges = generate_tri_lat(minx, miny, maxx, maxy, length);
        vector<Edge> lat_copy = lat_edges;
        rotate_edges(lat_copy, angle, pivot);
        unfixed.insert(unfixed.end(), lat_copy.begin(), lat_copy.end());
        if(lat_edges.size() > 0){
            get_borders(lat_edges, left, right);
            get_overhang(parallel_edges, left, right, overhang, length/2, change_map);
            fix_edges(lat_edges, change_map);
            /*rotate_edges(overhang, angle, pivot);
            tri_edges.insert(tri_edges.end(), overhang.begin(), overhang.end());
            rotate_edges(overhang, -angle, pivot);
            */
            //Add the left and right edges of the lattice placed within the
            //current grain
            for(piter = 0; piter < left.size() - 1; piter++){
                overhang.push_back(Edge(left[piter], left[piter + 1]));
            }
            for(piter = 0; piter < right.size() - 1; piter++){
                overhang.push_back(Edge(right[piter], right[piter + 1]));
            }
            lat_edges = remove_border(lat_edges, overhang);
            rotate_edges(overhang, angle, pivot);
            raw_pool.insert(raw_pool.end(),overhang.begin(),overhang.end());
            junction_pool.insert(junction_pool.end(),overhang.begin(),overhang.end());
            rotate_edges(lat_edges, angle, pivot);
            raw_pool.insert(raw_pool.end(), lat_edges.begin(), lat_edges.end());
            tri_edges.insert(tri_edges.end(), lat_edges.begin(), lat_edges.end());
            overhang.clear();
            left.clear();
            right.clear();
            lat_edges.clear();
            change_map.clear();
        }
        parallel_edges.clear();
    }

    /*for(auto capiter = end_caps.begin(); capiter!=end_caps.end(); capiter++){
        junction_pool.push_back(*capiter);
    }*/
    if(yesno("Report all edges?")){
        cout << "Enter the name of the file: ";
        getline(cin, response);
        jreport = fopen(response.c_str(), "w");
        for(Edge jedge : raw_pool){
            fprintf(jreport, "%lf %lf\n", jedge.p1.x, jedge.p1.y);
            fprintf(jreport, "%lf %lf\n\n", jedge.p2.x, jedge.p2.y);
        }
        fclose(jreport);
    }
    points = get_points(tri_edges);

    //Find all points belonging to an edge forming an end cap
    for(Edge next_edge : end_caps){
        end_points.emplace(next_edge.p1);
        end_points.emplace(next_edge.p2);
    }

    //After all grains have been filled with lattices, group edges in the pool
    //used to define the envelopes of junctions into polygons. Triangulate
    //these polygons, and add the edges from the resulting triangulations to
    //the general set of edges
    build_network_maps(junction_pool, imap, emap);
    cycle_list = build_cycles(emap);
    polygons = get_grain_polygons(imap, cycle_list);

    /*
    cout << "\n";
    for(vector<Point> nextPoly : polygons){
       for(poly_iter = 0; poly_iter < nextPoly.size(); poly_iter ++){
           cout << nextPoly[poly_iter].x << "\t" << nextPoly[poly_iter].y << "\n";
           cout << nextPoly[(poly_iter+1)%nextPoly.size()].x << "\t" << nextPoly[(poly_iter+1)%nextPoly.size()].y << "\n\n";
        }
        cout << "\n";
    }*/

    tri_opts << base << "q" << setprecision(4) << min_ang << "a" << setprecision(10) << max_area;
    if(enforce_delaunay) tri_opts << "D";
    len = strlen(tri_opts.str().c_str());
    options = (char *) malloc(sizeof(char) * (len + 1));
    sprintf(options, "%s", tri_opts.str().c_str());

    if((report_tri_envs = yesno("Report triangulation envelopes?"))){
        cout << "Enter the name for reporting: ";
        getline(cin, response);
        env_report = fopen(response.c_str(), "w");
    }

    //Triangulate each polygon in turn
    for(vector<Point> polygon : polygons){
        if(polygon.size() < 3) continue;
        if(has_end_point(polygon, end_points)) continue;

        if(report_tri_envs){
            for(int poly_iter = 0; poly_iter < polygon.size(); poly_iter ++){
                fprintf(env_report, "%lf %lf\n", polygon[poly_iter].x, polygon[poly_iter].y);
                fprintf(env_report, "%lf %lf\n\n", polygon[(poly_iter+1)%polygon.size()].x, polygon[(poly_iter+1)%polygon.size()].y);
            }
        }

        /*
        boundary_set.clear();
        for(tri_iter = 0; tri_iter < polygon.size(); tri_iter++){
            boundary_set.emplace(Edge(polygon[tri_iter],polygon[(tri_iter+1)%polygon.size()]));
        }
        */

        create_grain_input(&in, polygon);
        setup_output(&out);
        triangulate(options, &in, &out,(struct triangulateio *)NULL);
        //Make a map from indices to points in the triangulation
        for(tri_iter = 0; tri_iter < out.numberofpoints; tri_iter++){
            tri_point_map.insert(make_pair(tri_iter, Point(out.pointlist[2*tri_iter],out.pointlist[2*tri_iter + 1])));
            points.emplace(Point(out.pointlist[2*tri_iter],out.pointlist[2*tri_iter + 1]));
        }

        //Gather all edges in the triangulation
        for(tri_iter = 0; tri_iter < out.numberofedges; tri_iter ++){
            tri_edges.push_back(Edge(tri_point_map[out.edgelist[tri_iter*2]], tri_point_map[out.edgelist[tri_iter*2 + 1]]));
            /*if(boundary_set.find(toAdd) == boundary_set.end()){
                tri_edges.push_back(toAdd);
            }*/
        }

        //Delete triangulateio structures
        free(in.pointlist);
        free(in.segmentlist);
        free(in.segmentmarkerlist);
        free(out.pointlist);
        free(out.edgelist);
        free(out.edgemarkerlist);
        free(out.segmentlist);
        free(out.segmentmarkerlist);
        free(out.trianglelist);
        tri_point_map.clear();
    }

    if(report_tri_envs) fclose(env_report);
    free(options);
}

//Given a set of edges making a planar line graph of a large-scale network,
//triangulate the large-scale network, according to specified triangulation
//parameters. Extract all unique edges from the triangulation.
void make_triangulated_network(vector<Edge> envelope, double min_ang, double max_area, bool enforce_delaunay, vector<Edge> &tri_edges, set<Point> &points){

    triangulateio in, out;
    map<int, Point> imap, tri_points;
    map<int, list<int>> emap;
    vector<vector<int>> cycle_list;
    vector<vector<Point>> polygons;
    vector<Point> border;
    string base = "pQze";
    ostringstream tri_opts;
    int iter, len;
    Point p1, p2;
    char *options;
    vector<tuple_2d> holes;
    FILE *hole_file;
    string response;

    //Build maps from integer indices to points in the large-scale network
    //envelope, and from points' indices to lists of indices of those points'
    //neighbors
    build_network_maps(envelope, imap, emap);

    //Build lists of cycles of points in the planar line graph
    cycle_list = build_cycles(emap);

    //Use the cycles to group the large-scale edges into polygons
    get_polygons_simple(imap, cycle_list, polygons, border);

    //Find holes in the planar line graph, so these are not filled with 
    //triangles
    get_holes(polygons, holes); 

    /*
    if(yesno("Report holes?")){
        cout << "Enter hole file name: ";
        getline(cin, response);
        hole_file = fopen(response.c_str(), "w");
        for(tuple_2d hole : holes){
            fprintf(hole_file, "%lf %lf\n", get<0>(hole), get<1>(hole));
        }
        fclose(hole_file);
    }
    */

    //Prepare input for the triangulate function of Triangle, and produce the
    //triangulation
    tri_opts << base << "q" << setprecision(4) << min_ang << "a" << setprecision(10) << max_area;
    if(enforce_delaunay) tri_opts << "D";
    create_input(&in, imap.size(), polygons, border, holes);
    setup_output(&out);
    len = strlen(tri_opts.str().c_str());
    options = (char *) malloc(sizeof(char) * (len + 1));
    sprintf(options, "%s", tri_opts.str().c_str());
    triangulate(options, &in, &out,(struct triangulateio *)NULL);
    free(options);

    //Make a map from indices to points in the triangulation
    for(iter = 0; iter < out.numberofpoints; iter++){
        tri_points.insert(make_pair(iter, Point(out.pointlist[2*iter],out.pointlist[2*iter + 1])));
        points.emplace(Point(out.pointlist[2*iter],out.pointlist[2*iter + 1]));
    }

    //Gather all edges in the triangulation
    for(iter = 0; iter < out.numberofedges; iter ++){
       tri_edges.push_back(Edge(tri_points[out.edgelist[iter*2]], tri_points[out.edgelist[iter*2 + 1]]));
    } 

    //Delete triangulateio structures
    free(in.pointlist);
    free(in.segmentlist);
    free(in.segmentmarkerlist);
    free(out.pointlist);
    free(out.edgelist);
    free(out.edgemarkerlist);
    free(out.segmentlist);
    free(out.segmentmarkerlist);
    free(out.trianglelist);
}

/*
Helper function for sieve_edges that calculates indices for use as keys
to map pairs of adjacent polygons from a tiling of a network to edges straddling
those polygons
*/
int adj_map_index(int index1, int index2, int total){
    int low = index1 <= index2 ? index1 : index2;
    int high = low == index1 ? index2 : index1;

    return low * (total + total + 1 - low) / 2 + high;
}

//Create a data structure mapping grid cells to lists of polygon indices
//partially contained in those grid cells, as well as lists of edge indices
void make_poly_grid(vector<PolyData> pdata, vector<vector<int>> &grid_lists, double length, double &minx, double &miny, double &maxx, double &maxy, int &xdim, int &ydim){

    int data_iter, xStart, yStart, xEnd, yEnd, xIter, yIter;
    minx = miny = (double) FLT_MAX;
    maxx = maxy = (double) FLT_MIN;

    //Pass through once to find the minimum and maximum x and y coordinates
    for(PolyData data : pdata){
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
    for(data_iter = 0; data_iter < pdata.size(); data_iter++){
        PolyData data = pdata[data_iter];
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

void make_edge_grid(vector<Edge> edges, vector<vector<int>> &grid_lists, double length, double minx, double miny, double maxx, double maxy, int xdim, int ydim){

    int edge_iter, xStart, yStart, xEnd, yEnd, xIter, yIter;

    //Determine the number of grid cells along the x and y directions and
    //allocate the according number of grid lists. Insert a terminating -1
    //in each list, after which any appropriate edge indices will be added later
    for(edge_iter = 0; edge_iter < xdim * ydim; edge_iter++){
        grid_lists.push_back(vector<int>());
    }

    //Make a second pass through to find the grid boundaries of each polygon
    for(edge_iter = 0; edge_iter < edges.size(); edge_iter++){
        Edge nextEdge = edges[edge_iter];
        xStart = (int) floor((nextEdge.left - minx) / length);
        xEnd = (int) floor((nextEdge.right - minx) / length);
        yStart = (int) floor((nextEdge.bottom - miny) / length);
        yEnd = (int) floor((nextEdge.top - miny) / length);

        for(yIter = yStart; yIter <= yEnd; yIter ++){
            for(xIter = xStart; xIter <= xEnd; xIter ++){
                grid_lists[yIter * xdim + xIter].push_back(edge_iter);
            }
        }
    }
}

//Group points according to the large-scale bonds in which they lie. Also find
//all edges within a given large-scale bond, and produce lists of edges
//on the small scale that span large-scale bonds.
void sort_edges(double length, vector<Edge> bottom, vector<PolyData> pdat, vector<vector<Edge>> &edge_collection, map<int, vector<Edge>> &adj_map, bool getAdjMap){

    bool p1here, p2here, inpoly;
    int  p1_index, p2_index, map_index, num_polys;
    int xdim, ydim, xStart, xEnd, yStart, yEnd, xIter, yIter;
    double minx, miny, maxx, maxy;
    vector<vector<int>> grid_lists;
    vector<int> currList;
    Edge large_edge;

    //Sort edges by vertical position to narrow search
    sort(bottom.begin(), bottom.end());
    //sort(pdat.begin(), pdat.end());

    for(int i = 0; i < pdat.size(); i ++){
        edge_collection.push_back(vector<Edge>());
    }

    //Set up grid data structure to determine in which polygons a point may lie
    make_poly_grid(pdat, grid_lists, length, minx, miny, maxx, maxy, xdim, ydim);

    for(Edge nextEdge : bottom){
        p1_index = -1;
        p2_index = -1;
        inpoly = false;

        if(nextEdge.left >= minx && nextEdge.right <= maxx && nextEdge.bottom >= miny && nextEdge.top <= maxy){

            xStart  = (int) floor((nextEdge.left - minx) / length);
            xEnd = (int) floor((nextEdge.right - minx) / length);
            yStart  = (int) floor((nextEdge.bottom - miny) / length);
            yEnd = (int) floor((nextEdge.top - miny) / length);

            for(yIter = yStart; yIter <= yEnd; yIter ++){
                if(inpoly) break;
                for(xIter = xStart; xIter <= xEnd; xIter ++){
                    if(inpoly) break;
                    currList = grid_lists[yIter * xdim + xIter];

                    for(int pdat_index : currList){
                    //while((pdat_index = currList[list_iter++]) >= 0){

                        p1here = false;
                        p2here = false;
                        PolyData datum = pdat[pdat_index];

                        if(datum.contains_permissive(nextEdge.p1)){
                            p1here = true;
                            p1_index = pdat_index;
                        }
                        if(datum.contains_permissive(nextEdge.p2)){
                            p2here = true;
                            p2_index = pdat_index;
                        }
                        if(p1here && p2here){
                            edge_collection[pdat_index].push_back(nextEdge);
                            inpoly = true;
                            break;
                        }
                    }
                }
            }

            if(! inpoly){
                if(p1_index >= 0 && p2_index < 0){
                    edge_collection[p1_index].push_back(nextEdge);
                }
                else if(p2_index >= 0 && p1_index < 0){
                    edge_collection[p2_index].push_back(nextEdge);
                }
                else if(p1_index >= 0 && p2_index >= 0){
                    if(getAdjMap){
                        map_index = adj_map_index(p1_index, p2_index, num_polys);
                        if(adj_map.find(map_index) == adj_map.end()){
                            adj_map.insert(make_pair(map_index, vector<Edge>()));
                        }
                        adj_map[map_index].push_back(nextEdge);
                    }
                    else{
                        edge_collection[p1_index].push_back(nextEdge);
                        edge_collection[p2_index].push_back(nextEdge);
                    }
                 }
            }
        }
    }
}

//Sort small edges made from a random packing into large-scale tiles
void sort_random_edges(double length, vector<Edge> bottom, vector<PolyData> pdat, vector<vector<Edge>> &edge_collection, map<int, vector<Edge>> &adj_map, bool getAdjMap){

    bool inpoly;
    int  p1_index, p2_index, map_index, num_polys;
    int xdim, ydim, xStart, xEnd, yStart, yEnd, xIter, yIter;
    double minx, miny, maxx, maxy;
    vector<vector<int>> grid_lists;
    vector<int> currList;
    Edge large_edge;

    //Sort edges by vertical position to narrow search
    sort(bottom.begin(), bottom.end());
    //sort(pdat.begin(), pdat.end());

    for(int i = 0; i < pdat.size(); i ++){
        edge_collection.push_back(vector<Edge>());
    }

    //Set up grid data structure to determine in which polygons a point may lie
    make_poly_grid(pdat, grid_lists, length, minx, miny, maxx, maxy, xdim, ydim);

    for(Edge nextEdge : bottom){

        if(nextEdge.left >= minx && nextEdge.right <= maxx && nextEdge.bottom >= miny && nextEdge.top <= maxy){

            xStart  = (int) floor((nextEdge.left - minx) / length);
            xEnd = (int) floor((nextEdge.right - minx) / length);
            yStart  = (int) floor((nextEdge.bottom - miny) / length);
            yEnd = (int) floor((nextEdge.top - miny) / length);

            p1_index = -1;
            p2_index = -1;
            inpoly = false;

            for(yIter = yStart; yIter <= yEnd; yIter ++){
                if(inpoly) break;
                for(xIter = xStart; xIter <= xEnd; xIter ++){
                    if(inpoly) break;
                    currList = grid_lists[yIter * xdim + xIter];

                    for(int pdat_index : currList){
                    //while((pdat_index = currList[list_iter++]) >= 0){

                        PolyData datum = pdat[pdat_index];

                        if(datum.contains(nextEdge.p1)){
                            p1_index = pdat_index;
                        }
                        if(datum.contains(nextEdge.p2)){
                            p2_index = pdat_index;
                        }
                        if(p1_index >= 0 && p2_index >= 0){
                            inpoly = true;
                            break;
                        }
                    }
                }
            }

            if(p1_index >= 0 && p2_index >= 0){
                if(p1_index == p2_index){
                    edge_collection[p1_index].push_back(nextEdge);
                }
                else{
                    if(getAdjMap){
                        map_index = adj_map_index(p1_index, p2_index, num_polys);
                        if(adj_map.find(map_index) == adj_map.end()){
                            adj_map.insert(make_pair(map_index, vector<Edge>()));
                        }
                        adj_map[map_index].push_back(nextEdge);
                    }
                    /*else{
                        edge_collection[p1_index].push_back(nextEdge);
                        edge_collection[p2_index].push_back(nextEdge);
                    }*/
                }
            }
            else if(p1_index >= 0) edge_collection[p1_index].push_back(nextEdge);
            else if(p2_index >= 0) edge_collection[p2_index].push_back(nextEdge);
        }
    }
}

//Find all points with just one neighbor, and remove them from the network, as
//well as any edges to which they are connected
void remove_flappers(vector<Edge> &edges){

    vector<Edge> newEdges;
    map<Point, int> neighbor_tallies;
    Point p1, p2;
    set<Point> flappers;

    //Visit each edge. Every time a point is encountered in an edge, increment
    //the count tracking the number of neighbors that point has
    for(Edge nextEdge : edges){
        p1 = nextEdge.p1;
        p2 = nextEdge.p2;

        if(neighbor_tallies.find(p1) == neighbor_tallies.end()){
            neighbor_tallies.insert(make_pair(p1, 1));
        }
        else neighbor_tallies[p1] += 1;
        if(neighbor_tallies.find(p2) == neighbor_tallies.end()){
            neighbor_tallies.insert(make_pair(p2, 1));
        }
        else neighbor_tallies[p2] += 1;
    }

    for(auto iter = neighbor_tallies.begin(); iter != neighbor_tallies.end(); iter ++){
        if(iter->second == 1) flappers.insert(iter->first);
    }

    for(Edge nextEdge : edges){
        p1 = nextEdge.p1;
        p2 = nextEdge.p2;
        if(flappers.find(p1) == flappers.end() && flappers.find(p2) == flappers.end()){
            newEdges.push_back(nextEdge);
        }
    }

    edges = newEdges;
}

/*After small-scale edges have been placed within the large-scale envelope,
this function stitches together the network at the smaller length scale, given a
knowledge of edges at the small length scale contained within a polygonal tiling
of edges at the large length scale.
*/
vector<Edge> stitch_network(vector<vector<Edge>> collection, map<int, vector<Edge>> adj_map, set<Point> points, double toKeep){

    int canonical_count, root1, root2, index = 0, cLen, nJobs;
    int edgeCount = 0, numNeeded, discardIter;
    unordered_map<Point, int> point_map;
    vector<Edge> pool, next_adj_list, retain, edgeMST, discard;
    Edge front;
    vector<Edge>::iterator pool_iter;

    canonical_count = points.size();

    for(Point nextPoint : points){
        point_map.insert(make_pair(nextPoint, index++));
    }

    vector<int> mst_table(point_map.size(), -1);

    for(vector<Edge> nextSet : collection){
        edgeCount += nextSet.size();
        edgeMST = randomMST(nextSet, pool);
        for(Edge nextEdge : edgeMST){
            retain.push_back(nextEdge);
            root1 = find_root(mst_table, point_map[nextEdge.p1]);
            root2 = find_root(mst_table, point_map[nextEdge.p2]);
            if(root1 != root2){
                makeunion(mst_table, root1, root2);
                canonical_count --;
            }
        }
    }

    for(auto it = adj_map.begin(); it != adj_map.end(); it++){
        next_adj_list = it->second;
        edgeCount += next_adj_list.size();

        shuffle_edges(next_adj_list);
        retain.push_back(next_adj_list[0]);

        if(next_adj_list.size() > 1){
            pool.insert(pool.end(), next_adj_list.begin() + 1, next_adj_list.end());
        }
    }

    numNeeded = (int) (edgeCount * toKeep);

    shuffle_edges(pool);
    for(pool_iter = pool.begin(); pool_iter != pool.end(); pool_iter ++){
        if(canonical_count == 1) break;
        front = *pool_iter;
        root1 = find_root(mst_table, point_map[front.p1]);
        root2 = find_root(mst_table, point_map[front.p2]);
        if(root1 != root2){
            makeunion(mst_table, root1, root2);
            canonical_count --;
            retain.push_back(front);
        }
        else discard.push_back(front);
    }

    discard.insert(discard.end(), pool_iter, pool.end());

    shuffle_edges(discard);

    discardIter = 0;
    while(retain.size() < numNeeded && discardIter < discard.size()){
        retain.push_back(discard[discardIter ++]);
    }

    return retain;
}

vector<Edge> sieve_edges(vector<Edge> top, double length, vector<Edge> bottom, vector<PolyData> pdat, double toKeep, bool accept_straddlers){

    vector<Edge> allEdges;
    vector<vector<Edge>> edge_collection;
    map<int, vector<Edge>> adj_map;
    set<Point> inPoints;
    bool p1in, p2in, p1here, p2here, cross, inpoly;
    int intersections;
    int  p1_index, p2_index, map_index, num_polys, iter;
    int xdim, ydim, xStart, xEnd, yStart, yEnd, xIter, yIter;
    double minx, miny, maxx, maxy;
    vector<vector<int>> grid_lists, edge_grid_lists;
    vector<int> currList;
    Edge large_edge;
    FILE *reject_file = NULL;
    string name;

    if(yesno("Report rejects?")){
        cout << "Enter the reject file name: ";
        getline(cin, name);
        reject_file = fopen(name.c_str(), "w");
    }


    //Sort edges and polygons by vertical position to narrow search
    sort(top.begin(), top.end());
    sort(bottom.begin(), bottom.end());

    /*
    for(int i = 0; i < pdat.size(); i ++){
        edge_collection.push_back(vector<Edge>());
    }
    */

    //Set up grid data structure to determine in which polygons a point may
    //lie and which large-scale edges a small-scale edge may cross
    make_poly_grid(pdat, grid_lists, length, minx, miny, maxx, maxy, xdim, ydim);
    make_edge_grid(top, edge_grid_lists, length, minx, miny, maxx, maxy, xdim, ydim);


    for(Edge nextEdge : bottom){
        p1in = false;
        p2in = false;
        cross = false;
        intersections = 0;
        inpoly = false;

        if(nextEdge.left >= minx && nextEdge.right <= maxx && nextEdge.bottom >= miny && nextEdge.top <= maxy){

            xStart  = (int) floor((nextEdge.left - minx) / length);
            xEnd = (int) floor((nextEdge.right - minx) / length);
            yStart  = (int) floor((nextEdge.bottom - miny) / length);
            yEnd = (int) floor((nextEdge.top - miny) / length);

            for(yIter = yStart; yIter <= yEnd; yIter ++){
                if(inpoly) break;
                for(xIter = xStart; xIter <= xEnd; xIter ++){
                    if(inpoly) break;
                    currList = grid_lists[yIter * xdim + xIter];

                    for(int pdat_index : currList){

                        p1here = false;
                        p2here = false;
                        PolyData datum = pdat[pdat_index];

                        if(datum.contains(nextEdge.p1)){
                            p1here = true;
                            p1in = true;
                            inPoints.emplace(nextEdge.p1);
                            p1_index = pdat_index;
                        }
                        if(datum.contains(nextEdge.p2)){
                            p2here = true;
                            p2in = true;
                            inPoints.emplace(nextEdge.p2);
                            p2_index = pdat_index;
                        }
                        if(p1here && p2here){
                            /*if(toKeep < 1) edge_collection[pdat_index].push_back(nextEdge);
                            else allEdges.push_back(nextEdge);*/
                            //edge_collection[pdat_index].push_back(nextEdge);
                            allEdges.push_back(nextEdge);
                            inpoly = true;
                            break;
                        }
                    }
                }
            }

        //If an edge on the smaller length scale intersects multiple external
        //edges on the larger length scale, it is not contained entirely within
        //one bond on the large length scale and should be rejected.
	    if(p1in && p2in && !inpoly){
                intersections = 0;

                for(yIter = yStart; yIter <= yEnd; yIter ++){
                    for(xIter = xStart; xIter <= xEnd; xIter ++){
                        currList = edge_grid_lists[yIter * xdim + xIter];

                        for(int eIndex : currList){
                            large_edge = top[eIndex];
                            if(intersection(nextEdge, large_edge)){
                                intersections ++;
                                if(intersections == 2) break;
                            }
                        }
                    }
                }

                if(intersections < 2){
                    /*if(toKeep < 1){
                        map_index = adj_map_index(p1_index, p2_index, num_polys);
                        if(adj_map.find(map_index) == adj_map.end()){
                            adj_map.insert(make_pair(map_index, vector<Edge>()));
                        }
                        adj_map[map_index].push_back(nextEdge);
                    }*/
                    allEdges.push_back(nextEdge);
                }
                else if(reject_file != NULL){
                    fprintf(reject_file, "%lf %lf\n", nextEdge.p1.x, nextEdge.p1.y);
                    fprintf(reject_file, "%lf %lf\n\n", nextEdge.p2.x, nextEdge.p2.y);
                }
            }
            else if((p1in || p2in) && !inpoly && accept_straddlers){
                allEdges.push_back(nextEdge);
                /*if(p1in){
                    if(toKeep < 1) edge_collection[p1_index].push_back(nextEdge);
                    else allEdges.push_back(nextEdge);
                }
                else{
                    if(toKeep < 1) edge_collection[p2_index].push_back(nextEdge);
                    else allEdges.push_back(nextEdge);
                }*/
            }
        }
    }

    if(reject_file != NULL) fclose(reject_file);

    /*
    for(iter = 0; iter < edge_collection.size(); iter++){
        remove_flappers(edge_collection[iter]);
    }
    */
    remove_flappers(allEdges);

    cout << "Size here: " << allEdges.size() << "\n";

    if(toKeep < 1){
        sort_random_edges(length, allEdges, pdat, edge_collection, adj_map, true);
        return stitch_network(edge_collection, adj_map, inPoints, toKeep);
    }

    else{
        return allEdges;
    }
}

//Find the largest connected component of a graph defined by an edge list
vector<Edge> get_biggest_component(vector<Edge> edges){

    //Map from points to their indices
    map<Point, int> pmap;
    //Point drawn from the list of edges
    set<Point> points = get_points(edges);
    //List of canonical vertices for the Union-Find algorithm
    vector<int> canonical(points.size(), -1);
    //Indices for Union-Find lookups
    int root1, root2, pindex = 0, max_size = 0, max_index;
    //Mapping from surviving canonical indices to lists of edges in the
    //respective connected subgraphs
    map<int, vector<Edge>> component_map;

    //Prepare the mapping from points to indices
    for(Point p : points){
        pmap.insert(make_pair(p, pindex++));
    }

    //Group vertices in the graph by mapping them to canonical vertices
    for(Edge e : edges){
        root1 = find_root(canonical, pmap[e.p1]);
        root2 = find_root(canonical, pmap[e.p2]);
        if(root1 != root2) makeunion(canonical, root1, root2);
    }

    //Make lists of edges within each connected component, and return the
    //largest list
    for(Edge e : edges){
        root1 = find_root(canonical, pmap[e.p1]);
        if(component_map.find(root1) == component_map.end()){
            component_map.insert(make_pair(root1, vector<Edge>()));
        }

        component_map[root1].push_back(e);
    }

    for(auto iter = component_map.begin(); iter != component_map.end(); iter++){
        if(iter->second.size() > max_size){
            max_index = iter->first;
            max_size = iter->second.size();
        }
    }

    return component_map[max_index];
}

vector<vector<double>> build_angle_lists(vector<Edge> edges){
    unordered_map<Point, vector<double>> point_map;
    Point p1, p2;
    double dx, dy, angle;
    vector<vector<double>> angle_lists;

    for(Edge nextEdge : edges){
        p1 = nextEdge.p1;
        if(point_map.find(p1) == point_map.end()){
            point_map.insert(make_pair(p1, vector<double>()));
        }

        p2 = nextEdge.p2;
        if(point_map.find(p2) == point_map.end()){
             point_map.insert(make_pair(p2, vector<double>()));
        }

        dx = p2.x - p1.x;
        dy = p2.y - p1.y;
        angle = atan2(dy, dx);
        point_map[p1].push_back(angle);
        point_map[p2].push_back(angle - M_PI);
    }

    for(auto iter = point_map.begin(); iter != point_map.end(); iter ++){
        angle_lists.push_back(iter->second);
    }

    return angle_lists;
}

double calc_psi_6_mag(vector<vector<double>> angle_lists){
    gsl_complex psi_6_total, local_psi_6, arg;
    double psi_6_mag;

    psi_6_total = gsl_complex_rect(0,0);

    for(vector<double> next_list : angle_lists){
        local_psi_6 = gsl_complex_rect(0,0);
        for(double angle : next_list){
           arg = gsl_complex_rect(0,6*angle);
           local_psi_6 = gsl_complex_add(local_psi_6, gsl_complex_exp(arg));
        }
        local_psi_6 = gsl_complex_div_real(local_psi_6, next_list.size());
        psi_6_total = gsl_complex_add(psi_6_total, local_psi_6);
    }

    psi_6_total = gsl_complex_div_real(psi_6_total, angle_lists.size());
    psi_6_mag = gsl_complex_abs(psi_6_total);

    return psi_6_mag;
}

double calc_alignment(vector<Edge> small_edges, Edge skel_edge){
    double alignment_sum = 0, mag, dx, dy, skel_nx, skel_ny;

    dx = skel_edge.p2.x - skel_edge.p1.x;
    dy = skel_edge.p2.y - skel_edge.p1.y;
    mag = sqrt(dx*dx + dy*dy);
    skel_nx = dx / mag;
    skel_ny = dy / mag;

    for(Edge nextEdge : small_edges){
        dx = nextEdge.p2.x - nextEdge.p1.x;
        dy = nextEdge.p2.y - nextEdge.p1.y;
        mag = sqrt(dx*dx + dy*dy);
        alignment_sum += abs((dx * skel_nx + dy * skel_ny) / mag);
    }

    return alignment_sum / small_edges.size();
}

Point find_x_match(Point p1, Point p2, double x){
    if(p1.x == x) return p1;
    else return p2;
}

Point find_y_match(Point p1, Point p2, double y){
    if(p1.y == y) return p1;
    else return p2;
}

vector<Edge> trim_network(vector<Edge> edge_list, bool connect_top_bottom){
    double low, high, left, right;
    vector<double> cuts;
    bool valid, curr_defined;
    set<Point> bps, tps;
    Point p1, p2, new_point, curr;
    vector<Edge> replace;
    double intersect_x, intersect_y;

    do{
        cuts = getdoubles("Enter the left and right bounds: ");

        if(cuts.size() < 2){
            cerr << "Enter two numbers.\n";
            valid = false;
            continue;
        }

        left = cuts[0];
        right = cuts[1];

        if(right < left){
            cerr << "The right bound must not be less than the left bound.\n";
            valid = false;
        }

        else valid = true;
    }while(! valid);

    do{
        cuts = getdoubles("Enter the  heights for the lower and upper grips: ");

        if(cuts.size() < 2){
            cerr << "Enter two numbers.\n";
            valid = false;
            continue;
        }

        low = cuts[0];
        high = cuts[1];

        if(high < low){
            cerr << "The upper bound must not be less than the lower bound.\n";
            valid = false;
        }

        else valid = true;
    }while(! valid);


    for(Edge next : edge_list){
        p1 = next.p1;
        p2 = next.p2;

        //Enforce the x limits
        if(next.left > right || next.right < left) continue;

        else if(next.left < left && next.right > left){
            intersect_y = y_intersect(next, left);
            new_point = Point(left, intersect_y);
            next.reset(new_point, find_x_match(p1, p2, next.right));
        }

        else if(next.left < right && next.right > right){
            intersect_y = y_intersect(next, right);
            new_point = Point(right, intersect_y);
            next.reset(new_point, find_x_match(p1, p2, next.left));
        }

        //Enforce the y limits
        if(next.bottom >= low && next.top <= high){
            replace.push_back(next);
            if(next.bottom == low) bps.insert(find_y_match(p1, p2, low));
            if(next.top == high) tps.insert(find_y_match(p1, p2, high));
        }

        else if(next.bottom < low && next.top > low){
            intersect_x = x_intersect(next, low);
            new_point = Point(intersect_x, low);
            replace.push_back(Edge(new_point, find_y_match(p1, p2, next.top)));
            if(connect_top_bottom) bps.insert(new_point);
        }
        
        else if(next.bottom < high && next.top > high){
            intersect_x = x_intersect(next, high);
            new_point = Point(intersect_x, high);
            replace.push_back(Edge(find_y_match(p1, p2, next.bottom), new_point));
            if(connect_top_bottom) tps.insert(new_point);
        }

        else continue;
    }

    if(connect_top_bottom){
        curr_defined = false;
        for(Point next : bps){
            if(! curr_defined){
                curr = next;
                curr_defined = true;
            }
            else{
                replace.push_back(Edge(curr, next));
                curr = next;
            }
        }

        curr_defined = false;
        for(Point next : tps){
            if(! curr_defined){
                curr = next;
                curr_defined = true;
            }
            else{
                replace.push_back(Edge(curr, next));
                curr = next;
            }
        }
    }

    return replace;
}

//Given a list of edges and a cutoff numSDEV, return only those edges with
//lengths that are no more than numSDEV times the standard deviation of length
//below the mean length.
vector<Edge> cull_outliers(vector<Edge> original, double numSDEV){

    double mean = 0, sdev = 0, length, cutoff;
    Point p1, p2;
    map<double, vector<Edge>> edgeMap;
    vector<Edge> keepers;

    for(Edge next : original){
        p1 = next.p1;
        p2 = next.p2;
        length = sqrt((p1.x-p2.x)*(p1.x - p2.x) + (p1.y-p2.y)*(p1.y-p2.y));
        mean += length;
        sdev += length * length;
        if(edgeMap.find(length) == edgeMap.end()){
            edgeMap.insert(make_pair(length, vector<Edge>()));
        }
        edgeMap[length].push_back(next);
    }

    mean /= original.size();
    sdev /= original.size();
    sdev = sqrt(sdev - mean * mean);
    cutoff = mean - numSDEV * sdev;

    for(auto iter = edgeMap.upper_bound(cutoff); iter != edgeMap.end(); iter++){
        for(Edge next : iter->second) keepers.push_back(next);
    }

    return keepers;
}

int main(int argc, char **argv){

    vector<double> bounds;
    vector<vector<double>> rules;
    map<int, vector<vector<double>>> nns;
    double scale = 0, width = 0, to_keep, min_ang, max_area, smallP;
    double min_dist, small_min_dist, numSDEV, density, sdev;
    bool success, get_alignment = false, get_psi_6 = false;
    int num_read, num_rows, iter, num_after = 0, nthreads = 1;
    string response, tri_args;
    vector<Edge> large_edges, small_edges, envelope;
    vector<vector<Edge>> edge_collection;
    map<int, vector<Edge>> adj_map;
    set<Point> small_points;
    bool enforce_delaunay, connect_top_bottom, accept_straddlers = true;
    vector<PolyData> pdata;
    FILE *netfile, *ranfile, *p6file = NULL, *alignment_file = NULL;
    unsigned seed;
    char num_string[11];
    vector<vector<double>> angle_lists;
    smallMode smode;
    largeMode lmode;
    char choice, flag;
    ifstream input;

    //Process optional command line switch
    while((flag = getopt(argc, argv, "r")) != -1){
        switch(flag){
            case 'r':
                accept_straddlers = false;
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

    //Prompt for bounds, the file describing the undisturbed large-scale
    //network, the length of undistorted large-scale bonds, the portion of
    //large-scale bonds to retain, the width of large-scale bonds,  and the 
    //standard deviation of random displacements of large-scale nodes.
    while(true){
        bounds = getdoubles("Enter bottom left and top right bounds: ");
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

    success = false;
    while(! success){
        cout << "Enter the approach for creating the large scale network.\n";
        cout <<"Enter c for crystalline, p for poisson, or l for lloyd: ";
        getline(cin, response);
        if(response.size() == 0){
            cerr << "No response read.\n";
            continue;
        }
        choice = response.at(0);
        switch(choice){
            case 'c':
                lmode = largeMode::crystalline;
                success = true;
                break;
            case 'p':
                lmode = largeMode::poisson;
                success = true;
                break;
            case 'l':
                lmode = largeMode::lloyd;
                success = true;
                break;
            default:
                cerr << "Enter either c, p, or l.\n";
                break;
        }
    }

    if(lmode == largeMode::crystalline){
        while(true){
            cout << "Enter the large-scale bond length: ";
            getline(cin, response);
            num_read = sscanf(response.c_str(), "%lf", &scale);
            if(num_read == 0 || scale <= 0) cerr<<"Enter a positive number.\n";
            else break;
        }
        /*
        while(true){
            cout << "Enter the standard deviation for displacement: ";
            getline(cin, response);
            num_read = sscanf(response.c_str(), "%lf", &sdev);
            if(num_read > 0 && sdev >= 0) break;
            else cerr << "Enter a non-negative number.\n";
        }
        */
    }

    else if(lmode == largeMode::poisson){
        while(true){
            cout << "Enter the minimum large-scale point separation: ";
            getline(cin, response);
            num_read = sscanf(response.c_str(), "%lf", &min_dist);
            if(num_read==0 || min_dist <= 0) cerr<<"Enter a positive number.\n";
            else break;
        }
        scale = min_dist;
    }

    else{
        while(true){
            cout << "Enter the density of large-scale points: ";
            getline(cin, response);
            num_read = sscanf(response.c_str(), "%lf", &density);
            if(num_read==0 || density <= 0) cerr<<"Enter a positive number.\n";
            else break;
        }
        scale = 1 / sqrt(density);
    }

    while(true){
        cout << "Enter the large-scale bond width: ";
        getline(cin, response);
        num_read = sscanf(response.c_str(), "%lf", &width);
        if(num_read == 0 || width <= 0) cerr << "Enter a positive number.\n";
        else break;
    }

    while(true){
        cout << "Enter the portion of large-scale bonds to keep: ";
        getline(cin, response);
        num_read = sscanf(response.c_str(), "%lf", &to_keep);
        if(num_read == 0) cerr << "Enter a number\n";
        else if(to_keep <= 0 || to_keep > 1){
            cerr << "Enter a number greater than 0 and less than or equal to 1.\n";
        }
        else break;
    }

    success = false;
    while(! success){
        cout << "Enter the approach for creating the small scale network.\n";
        cout <<"Enter g for grain-based, t for triangulated, r for random, or i to import: ";
        getline(cin, response);
        if(response.size() == 0){
            cerr << "No response read.\n";
            continue;
        }
        choice = response.at(0);
        switch(choice){
            case 'g':
                smode = smallMode::grain;
                success = true;
                break;
            case 't':
                smode = smallMode::triangulate;
                success = true;
                break;
            case 'r':
                smode = smallMode::random;
                success = true;
                break;
            case 'i':
                smode = smallMode::import;
                success = true;
            default:
                cerr << "Enter either g, t, r or i.\n";
                break;
        }
    }

    //If the user chooses a random small-scale network, ask for the minimum
    //distance between small-scale nodes
    if(smode == smallMode::random){
        while(true){
            cout << "Enter the minimum small-scale separation: ";
            getline(cin, response);
            num_read = sscanf(response.c_str(), "%lf", &small_min_dist);
            if(num_read > 0 && small_min_dist > 0) break;
            else cerr << "Enter a positive number.\n";
        }
    }
 
    //Ask if a "grain-based" approach should be used, in which large-scale
    //polygonal tiles are partially filled with lattices and stitched together
    //at junctions
    if(smode == smallMode::grain){
        while(true){
            cout << "Enter the number of small-scale rows per grain: ";
            getline(cin, response);
            num_read = sscanf(response.c_str(), "%d", &num_rows);
            if(num_read == 1 && num_rows > 0) break;
            else cerr << "Enter a positive integer.\n";
        }
    }
        

    //Prepare to generate pseudorandom numbers
    ranfile = fopen("/dev/urandom", "r");
    fread(&seed, sizeof(unsigned), 1, ranfile);
    fclose(ranfile);
    sprintf(num_string, "%u", seed);
    setenv("GSL_RNG_SEED", num_string, 1);
    gsl_rng_env_setup();

    if(lmode == largeMode::crystalline){
        do{
            success = import_lattice(rules, nns, scale);
        }while(! success);

        //Make the large-scale network, dilute as appropriate, and add thickness
        //to the large scale edges to create the envelope.
        large_edges = makeedges(rules, nns, bounds);

        if(yesno("Add disorder to large points?")){
            while(true){
                cout << "Enter the standard deviation: ";
                getline(cin, response);
                num_read = sscanf(response.c_str(), "%lf", &sdev);
                if(num_read == 1 && sdev > 0) break;
                else cerr << "Enter a positive number.\n";
            }
            displace_points_grn(large_edges, sdev);
        }
    
    }

    else if(lmode == largeMode::poisson){
        large_edges = make_poisson_edges(bounds, min_dist);
    }

    else large_edges = get_lloyd_edges(bounds, density, nthreads);

    if(to_keep < 1){
        large_edges = random_connected(large_edges, to_keep);
    }

    if(smode == smallMode::grain){
        add_thickness(large_edges, width, pdata, false);
    }
    else { 
        envelope = add_thickness(large_edges, width, pdata, true);
    }

    if(yesno("Report polygons?")){
        cout << "Enter the file name: ";
        getline(cin, response);
        netfile = fopen(response.c_str(), "w");
        for(PolyData pdatum : pdata){
            for(Edge nextEdge : pdatum.polyEdges){
                fprintf(netfile, "%lf %lf\n", nextEdge.p1.x, nextEdge.p1.y);
                fprintf(netfile, "%lf %lf\n\n", nextEdge.p2.x, nextEdge.p2.y);
            }
        }
        fclose(netfile);
    }

    if(smode != smallMode::grain && yesno("Report envelope?")){
        cout << "Enter the file name: ";
        getline(cin, response);
        netfile = fopen(response.c_str(), "w");
        for(Edge next_edge : envelope){
            fprintf(netfile, "%lf %lf\n", next_edge.p1.x, next_edge.p1.y);
            fprintf(netfile, "%lf %lf\n\n", next_edge.p2.x, next_edge.p2.y);
        }
        fclose(netfile);
    }
        
    //Next, prompt for the parameters of the small-scale triangulation, and the
    //portion of small-scale bonds to retain. Offer the opportunity to add
    //small-scale disorder.
    if(smode != smallMode::random && smode != smallMode::import){
        while(true){
            cout << "Enter the minimum internal angle and maximum area: ";
            getline(cin, response);
            num_read = sscanf(response.c_str(), "%lf %lf", &min_ang, &max_area);
            if(num_read == 2 && min_ang > 0 && max_area > 0) break;
            else cerr << "Enter two positive numbers.\n";
        }

        enforce_delaunay = yesno("Enforce Delaunay condition?");
    }

    while(true){
        cout << "Enter the small-scale bond portion: ";
        getline(cin, response);
        num_read = sscanf(response.c_str(), "%lf", &smallP);
        if(num_read < 1 || smallP <= 0 || smallP > 1){
            cerr << "Enter a number greater than 0 and less than or equal to 1\n";
        }
        else break;
    }

    if(smode == smallMode::grain){
        make_grain_based_network(large_edges, pdata, width, num_rows, min_ang, max_area, enforce_delaunay, small_edges, small_points);
    }
    else if(smode == smallMode::triangulate){
        make_triangulated_network(envelope, min_ang, max_area, enforce_delaunay, small_edges, small_points);
    }
    else if(smode == smallMode::random){
        /*bounds[0] -= scale;
        bounds[1] -= scale;
        bounds[2] += scale;
        bounds[3] += scale;*/
        small_edges = make_poisson_edges(bounds, small_min_dist);
        //cout << "Reached here.\n";
        small_edges = sieve_edges(envelope, scale, small_edges, pdata, smallP, accept_straddlers);
        cout << "New size: " << small_edges.size() << "\n";
        //remove_flappers(small_edges);
        //cout << "Size here: " << small_edges.size() << "\n";
        /*if(smallP < 1){
            sort_random_edges(scale, small_edges, pdata, edge_collection, adj_map, true);
            for(vector<Edge> evec : edge_collection) num_after += evec.size();
            for(auto mIter = adj_map.begin(); mIter != adj_map.end(); mIter++){
                num_after += mIter->second.size();
            }
            cout << "Number now: " << num_after << "\n";
            small_edges = stitch_network(edge_collection, adj_map, small_points, smallP);
        }*/
    }

    else{
        open_dat_file("Enter the name of the edge file: ", input);
        small_edges = import_edges(input);
        small_edges = sieve_edges(envelope, scale, small_edges, pdata, smallP, accept_straddlers);
    }

    if(smode != smallMode::random && smode != smallMode::import && smallP < 1){
        sort_edges(scale, small_edges, pdata, edge_collection, adj_map, true);
        /*if(yesno("Report edge collection edges?")){
            cout << "Enter the name of the file: ";
            getline(cin, response);
            netfile = fopen(response.c_str(), "w");
            for(vector<Edge> nextVec : edge_collection){
                for(Edge next_edge : nextVec){
                    fprintf(netfile, "%lf %lf\n", next_edge.p1.x, next_edge.p1.y);
                    fprintf(netfile, "%lf %lf\n\n", next_edge.p2.x, next_edge.p2.y);
                }
            }
            fclose(netfile);
        }*/
        small_edges = stitch_network(edge_collection, adj_map, small_points, smallP);
    }

    
    if(yesno("Add disorder to small points?")){
        while(true){
            cout << "Enter the standard deviation: ";
            getline(cin, response);
            num_read = sscanf(response.c_str(), "%lf", &sdev);
            if(num_read == 1 && sdev > 0) break;
            else cerr << "Enter a positive number.\n";
        }
        displace_points_grn(small_edges, sdev);
    }
    

    if(lmode != largeMode::crystalline){
        if(yesno("Remove short bonds?")){
            do{
                cout << "Enter the cutoff, in standard deviations: ";
                getline(cin, response);
                num_read = sscanf(response.c_str(), "%lf", &numSDEV);
                if(num_read == 0 || numSDEV < 0){
                    cerr << "Enter a positive number.\n";
                }
                else break;
            }while(true);
            small_edges = cull_outliers(small_edges, numSDEV);
        }
    }

    if((get_psi_6 = yesno("Report psi-6 values?"))||(get_alignment = yesno("Report alignment of small and large-scale bonds?"))){

        pdata.clear();
        add_thickness(large_edges, width*1.1, pdata, false);

        //Open the psi 6 file, if requested
        if(get_psi_6){
            cout << "Enter the file name for psi-6 report: ";
            getline(cin, response);
            p6file = fopen(response.c_str(), "w");
        }

        //Open the alignment file, if requested
        if(get_alignment){
            cout << "Enter the name for the alignment report: ";
            getline(cin, response);
            alignment_file = fopen(response.c_str(), "w");
        }

        //Group small-scale edges by large-scale tile
        edge_collection.clear();
        if(smode == smallMode::random){
            sort_random_edges(scale, small_edges, pdata, edge_collection, adj_map, false);
        }
        else{
            sort_edges(scale, small_edges, pdata, edge_collection, adj_map, false);
        }
        netfile = NULL;
        if(yesno("Report edge collection edges?")){
            cout << "Enter the name of the file: ";
            getline(cin, response);
            netfile = fopen(response.c_str(), "w");
            for(vector<Edge> nextVec : edge_collection){
                for(Edge next_edge : nextVec){
                    fprintf(netfile, "%lf %lf\n", next_edge.p1.x, next_edge.p1.y);
                    fprintf(netfile, "%lf %lf\n\n", next_edge.p2.x, next_edge.p2.y);
                }
            }
            fclose(netfile);
        }
        //Iterate over groups of small-scale edges, and calculate and report
        //psi-6 or aligment values, or both
        iter = 0;
        for(vector<Edge> nextSet : edge_collection){
            if(p6file != NULL){
                angle_lists = build_angle_lists(nextSet);
                fprintf(p6file, "%lf\n", calc_psi_6_mag(angle_lists));
            }
            if(alignment_file != NULL){
                fprintf(alignment_file, "%lf\t%ld\n", calc_alignment(nextSet, large_edges[iter]), nextSet.size());
            }
            iter ++;
        }
        if(p6file != NULL) fclose(p6file);
        if(alignment_file != NULL) fclose(alignment_file);
    }

    //Trim the final network to size, and write it to a file.
    if(yesno("Trim the network?")){
        connect_top_bottom = yesno("Prepare for grips?");
        small_edges = get_biggest_component(trim_network(small_edges, connect_top_bottom));
    }

    cout << "Enter the name for the network file: ";
    getline(cin, response);
    netfile = fopen(response.c_str(), "w");
    for(Edge next_edge : small_edges){
        fprintf(netfile, "%lf %lf\n", next_edge.p1.x, next_edge.p1.y);
        fprintf(netfile, "%lf %lf\n\n", next_edge.p2.x, next_edge.p2.y);
    }

    fclose(netfile);

    return 1;
}
