/*Author: Jonathan Michel
This program prompts a user for boundaries and rules for producing lattice
sites and nearest neighbor connections. From here, a set of line segments
joining nearest neighbor locations is produced. Connections can be removed
from the network in two ways; one method is completely random, while the other
uses randomness but also guarantees that the network remains fully connected.
*/

#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <string>
#include <iostream>
#include <map>
#include <set>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <fstream>
#include <stack>
#include "network_utils.h"
#include <unistd.h>
#include <ctype.h>
#include <float.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;

/*
#define PRIME1 (int) 1021297
#define PRIME2 (int) 1021651 //Primes for computing hash functions
#define FLOAT_TOL 1E-6 //Tolerance for deeming floating point numbers different
#define BIG_SLOPE 1E10
#define INF 1E12
*/

struct NetData{

    NetData(vector<vector<double>> rls, map<int,vector<vector<double>>> nn, double w, double tk, bool prtct) : rules(rls), nns(nn), width(w), to_keep(tk), protect(prtct) {
        displace = false;
        sdev = 0;
    }
    
    NetData(vector<vector<double>> rls, map<int,vector<vector<double>>> nn, double w, double tk, double dev, bool prtct, bool disp) : rules(rls), nns(nn), width(w), to_keep(tk), sdev(dev), protect(prtct), displace(disp) {}

    vector<vector<double>> rules;
    map<int,vector<vector<double>>> nns;
    double width;
    double to_keep;
    double sdev;
    bool protect;
    bool displace;
};

struct MST_JOB{

    MST_JOB(vector<vector<Edge>> *arg1, int arg2, int arg3, int *arg4, vector<Edge> *arg5, vector<Edge> *arg6, vector<int> *arg7, unordered_map<Point, int> *arg8, int *arg9) : collection(arg1), begin(arg2), end(arg3), edgeCount(arg4), pool(arg5), retain(arg6), mst_table(arg7), point_map(arg8), canonical_count(arg9){}

    vector<vector<Edge>> *collection; 
    int begin;
    int end;
    int *edgeCount;
    vector<Edge> *pool;
    vector<Edge> *retain; 
    vector<int> *mst_table;
    unordered_map<Point, int> *point_map;
    int *canonical_count;
};

void print_edges(vector<Edge> edges, string message){    
    string nextline, filename;
    FILE *out;
    vector<string> tokens;

    cout << message;
    getline(cin, nextline);
    tokens = split(nextline, ' ');
    if(tokens.size() > 0){
        filename = tokens[0];
        out = fopen(filename.c_str(), "w");

        for(Edge e : edges){
            fprintf(out, "%10.8lf %10.8lf \n%10.8lf %10.8lf \n\n", e.p1.x, e.p1.y, e.p2.x, e.p2.y);
        }

       fclose(out);
    }
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

double get_v_offset(vector<vector<double>> rules){
    double voffset = 0;

    for(vector<double> rule : rules){
        voffset += rule[1];
    }

    return voffset;
}

double get_h_offset(vector<vector<double>> rules){
    double hoffset = 0;
    int index;

    for(index = 2; index < rules[0].size(); index++){
        hoffset += rules[0][index];
    }

    return hoffset;
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

vector<Edge> top_bottom(vector<Edge> total, set<Point> points){
    Point p1, p2, bleft, bright, tleft, tright, curr;
    vector<Edge> topbottom;
    map<Point,Edge> neighbor_map;
    vector<Edge> perimeter;
    Edge next;
    double midx, midy, dx, dy, slope, xcent = 0, ycent = 0;
    double miny = INF, maxy = -INF;
    bool onborder;
    int iter;

    bleft = Point(INF,INF);
    bright = Point(-INF,INF);
    tleft = Point(INF,-INF);
    tright = Point(-INF,-INF);

    for(Point p : points){
        xcent += p.x;
        ycent += p.y;
    }
    xcent /= points.size();
    ycent /= points.size();

    for(Edge e : total){
        p1 = e.p1;
        p2 = e.p2;
        midx = (p1.x + p2.x)/2;
        midy = (p1.y + p2.y)/2;
        dx = midx - xcent;
        dy = midy - ycent;

        onborder = true;
        iter = 0;
        do{
            next = total[iter];
            if(!(next == e)){
                onborder = !intersection(Edge(Point(midx,midy), Point(midx+dx,midy+dy)),next);
            }
            iter ++;
        }while(onborder && iter < total.size());
        if(onborder) perimeter.push_back(e);
    }
    //for(Edge e : perimeter) printf("(%lf,%lf)->(%lf,%lf)\n",e.p1.x,e.p1.y,e.p2.x,e.p2.y);

    for(Edge e : perimeter){
        p1 = e.p1;
        p2 = e.p2;

        if(p1.y <= miny){
            if(p1.y < miny) miny = p1.y;
            if(p1.x < bleft.x) bleft = p1;
            if(p1.x > bright.x) bright = p1;
        }
        if(p1.y >= maxy){
            if(p1.y > maxy){
                maxy = p1.y;
                tleft = p1;
                tright = p2;
            }
            if(p1.x < tleft.x) tleft = p1;
            if(p1.x > tright.x) tright = p1;
        }
        
        if(p2.y <= miny){
            if(p2.y < miny) miny = p2.y;
            if(p2.x < bleft.x) bleft = p2;
            if(p2.x > bright.x) bright = p2;
        }
        if(p2.y >= maxy){
            if(p2.y > maxy){
                maxy = p2.y;
                tleft = p2;
                tright = p2;
            }
            if(p2.x < tleft.x) tleft = p2;
            if(p2.x > tright.x) tright = p2;
        }
        
        if(neighbor_map.find(p1) == neighbor_map.end()){
            neighbor_map.insert(make_pair(p1,e));
        }
        else if(p2.x > neighbor_map[p1].p2.x) neighbor_map[p1] = e;
    }

    curr = bleft;
    while(!(curr == bright)){
        next = neighbor_map[curr];
        topbottom.push_back(next);
        curr = next.p2;
    }
    curr = tleft;
    while(!(curr == tright)){
        next = neighbor_map[curr];
        topbottom.push_back(next);
        curr = next.p2;
    }
    for(Edge e : topbottom) printf("(%lf,%lf)->(%lf,%lf)\n",e.p1.x,e.p1.y,e.p2.x,e.p2.y);
    return topbottom;
}

vector<Edge> connected(vector<Edge> edges, double p, set<Point> points){
    int numkept = 0, size = points.size(), root1, root2;
    int numneeded = edges.size() * p;
    vector<int> canonical(size, -1);
    vector<Edge> kept, topbottom;
    unordered_set<Edge> tb;
    Edge e;
    int considered = 0, ignored = 0;

    if(yesno("Protect top and bottom edges?")){
        topbottom = top_bottom(edges, points);
        for(Edge next : topbottom){
            root1 = find_root(canonical, points.find(next.p1)->index);
            root2 = find_root(canonical, points.find(next.p2)->index);
            makeunion(canonical, root1, root2);
            kept.push_back(next);
            numkept ++;
        }
    }

    shuffle_edges(edges);
    while(numkept < size - 1){
        e = edges[0];
        edges.erase(edges.begin());
        if(tb.find(e) == tb.end()){
            considered++;
            root1 = find_root(canonical, points.find(e.p1)->index);
            root2 = find_root(canonical, points.find(e.p2)->index);
            if(root1 != root2){
                makeunion(canonical, root1, root2);
                kept.push_back(e);
                numkept ++;
            }
            else edges.push_back(e);
        }
    }

    shuffle_edges(edges);
    while(numkept < numneeded){
        e = edges[0];
        edges.erase(edges.begin());
        if(tb.find(e) == tb.end()){
            kept.push_back(e);
            numkept ++;
        }
    }

    return kept;
}

vector<Edge> true_random(vector<Edge> original, double toKeep){
    vector<Edge> remaining;
    int numNeeded = (int) (toKeep*original.size());
    int count = 0;
    
    shuffle_edges(original);
    for(Edge e : original){
        remaining.push_back(e);
        if(++count == numNeeded) return remaining;
    }
	
    return remaining;
}

vector<Edge> randremove(vector<Edge> original, set<Point> points){
    double cutoff, randval;
    string response;
    vector<Edge> remaining;
    int option = 0;

    while(true){
        vector<double> answer = getdoubles("Enter probability to keep: ");
        if(!answer.empty()){
		if(answer[0] >= 0 && answer[0] <= 1){
                    cutoff = answer[0];
                    break;
                }
                else fprintf(stderr, "Enter a number between 0 and 1.\n");
        }
        else fprintf(stderr, "Enter a number.");
    }

    if(original.size()*cutoff >= points.size() -1){
        while(true){
            printf("Enter 0 for true random or 1 for connected: ");
            getline(cin, response);
            if(sscanf(response.c_str(), "%d", &option)){
                if(option == 0 || option == 1) break;
                else fprintf(stderr, "Enter 0 or 1.\n");
            }
            else fprintf(stderr, "Enter an integer.\n");
        }
    }
    else fprintf(stderr, "Too few edges kept to for a connected graph.\n");

    switch(option){
        case 0:
            return true_random(original, cutoff);
        case 1:
            return connected(original, cutoff, points);
        default:
            break;
    }
	
    return remaining;
}

int numinset(vector<Edge> edges, set<Point> points){
    int count = 0;

    for(Edge e : edges){
        if(points.count(e.p1) > 0 && points.count(e.p2) > 0) count ++;
    }

    return count;
}

vector<int> network_groups(vector<Edge> edges, set<Point> points){
    vector<int> canonical(points.size(), -1);
    vector<int> results(points.size(), 0);
    int root1, root2, max_zero;

    for(Edge e : edges){
        root1 = find_root(canonical, points.find(e.p1)->index);
        root2 = find_root(canonical, points.find(e.p2)->index);
        if(root1 != root2) makeunion(canonical, root1, root2);
    }
    
    for(int i = 0; i < points.size(); i ++){
        results[find_root(canonical,i)]++;
    }

    sort(results.begin(), results.end());
    for(int i = 0; i < results.size(); i++){
        if(!results[i]) max_zero = i;
        else break;
    }
    results.erase(results.begin(), results.begin() + max_zero + 1);
    return results;
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

/*
void truncate_miter(Point original, double slope, double y_ext){
}
*/

vector<Edge> add_thickness(vector<Edge> current, double thickness, vector<PolyData>& pdata, bool makepoly, bool level){
    map<Point, vector<double>> angmap;
    vector<Edge> replace;
    Point p1, p2, key;
    Point p1f, p2f, p3f, p4f, p1fb, p2fb, p3fb, p4fb;
    double low, high, ang1, ang2, dx, dy, hwidth, slope;
    double ymin = FLT_MAX, ymax = FLT_MIN, ylow, yhigh, y_ext;
    bool p1fflag, p2fflag, p3fflag, p4fflag, p1_is_end, p2_is_end;
    vector<double> angvec;

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
                y_ext = p2.y == ymin ? ylow : yhigh;
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
                }
                else pedges.push_back(Edge(p2f,p2));
            
                if(p4fflag){
                    pedges.push_back(Edge(p2, p4fb));
                    pedges.push_back(Edge(p4fb, p4f));
                }
                else pedges.push_back(Edge(p2,p4f));
            }
            else pedges.push_back(Edge(p2f, p4f));

            pedges.push_back(Edge(p4f,p3f));

            if(! p1_is_end){
                if(p3fflag){
                    pedges.push_back(Edge(p3f, p3fb));
                    pedges.push_back(Edge(p3fb, p1));
                }
                else pedges.push_back(Edge(p3f,p1));
                if(p1fflag){
                    pedges.push_back(Edge(p1, p1fb));
                    pedges.push_back(Edge(p1fb, p1f));
                }
                else pedges.push_back(Edge(p1,p1f));
            }
            else pedges.push_back(Edge(p3f, p1f));

            pdata.push_back(PolyData(pedges));
        }
    }

    return replace;
}

vector<Edge> trim_top_bottom(vector<Edge> original, set<Point> points){
    vector<Edge> replacement, topbottom;
    unordered_set<Edge> cull;

    topbottom = top_bottom(original, points);
    cull.insert(topbottom.begin(), topbottom.end());

    for(Edge e : original){
        if(cull.find(e) == cull.end()) replacement.push_back(e);
    }

    return replacement;
}

vector<vector<double>> getrules(){
    vector<double> rule;
    vector<vector<double>> rules;
    bool newRule;

    while(true){
        newRule = yesno("Enter a new rule? ");
        if(!newRule){
            if(rules.size()){
                break;
            }
            else{
                fprintf(stderr, "Enter at least one rule.\n");
            }
        }
        else{
            rule = getdoubles("Enter rule: ");
            if(rule.size() < 3){
                fprintf(stderr, "Enter at least three numbers.\n");
            }
            else{
                rules.push_back(rule);
            }
        }
    }

    return rules;
}

map<int,vector<vector<double>>> getnns(vector<vector<double>> rules){
    map<int,vector<vector<double>>> nns;
    vector<double> nn;
    int ruleiter = 0, count = 0, pointiter;
    string prompt, base;

    base = "Add a nearest neighbor vector for rule ";

    for(vector<double> nextrule : rules){
       ruleiter ++;
       for(pointiter = 1; pointiter <= nextrule.size() - 2; pointiter ++){
           vector<vector<double>> nextset;
           while(true){
               prompt = base + to_string(ruleiter) + " point " + to_string(pointiter) + "?";
               if(!yesno(prompt)) break;
               nn = getdoubles("Enter the vector: ");
               if(nn.size() != 2) fprintf(stderr, "Enter two numbers.\n");
               else{
                   nextset.push_back(nn);
               }
           }
           nns.insert(pair<int, vector<vector<double>>>(count++, nextset));
       }
    }

    return nns;
}

void scale_vector(vector<double>& in, double scale){
    int index;
    for(index = 0; index < in.size(); index++){
        in.at(index) = in.at(index)*scale;
    }
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

void get_lattice_info(vector<vector<double>>& rules, map<int,vector<vector<double>>>& nns){
    bool success;
    double scale;

    if(yesno("Read lattice data from file?")){
        do{
            scale = getdoubles("Enter the scale factor: ").at(0);
            success = import_lattice(rules, nns, scale);
            if(!success) if(!yesno("Read failed. Try again?")) break;
        }while(!success);
    }

    if(!success){
        rules = getrules();
        nns = getnns(rules);
    }
}

void loadNetStack(stack<NetData>& netStack){
    //vector<vector<double>> rules;
    //map<int,vector<vector<double>>> nns;
    double width, toKeep;
    bool valid, protect, displace;
    vector<double> response;
    double sdev;

    do{
        vector<vector<double>> rules;
        map<int,vector<vector<double>>> nns;
        get_lattice_info(rules, nns);
        valid = true;
        do{
            response = getdoubles("Enter the bond width: ");
            if(!response.size() || response.at(0) < 0){
                fprintf(stderr,"Enter a non-negative number.\n");
                valid = false;
            }
            else if(response.at(0) == 0 && netStack.size()){
                fprintf(stderr, "The width must be greater than 0.\n");
                valid = false;
            }
            else{
                width = response.at(0);
                valid = true;
            }
        }while(!valid);
       
        valid = true;
        do{
            response = getdoubles("Enter portion of bonds to keep: ");
            if(!response.size() || response.at(0) < 0){
                fprintf(stderr, "Enter a non-negative number.\n");
                valid = false;
            }
            else{
                toKeep = response.at(0);
                valid = true;
            }
        }while(!valid); 
        protect = yesno("Protect top and bottom edges?");

        displace = yesno("Displace edges?");
        if(displace){
            do{
                response = getdoubles("Enter deviation standard deviation: ");
                if(! response.size() || response[0] < 0){
                    cerr << "Enter a non-negative number.\n";
                    valid = false;
                }
                else{
                    sdev = response[0];
                    valid = true;
                }
            }while(! valid);
        }

        if(!displace) netStack.push(NetData(rules, nns, width, toKeep, protect));
        else netStack.push(NetData(rules, nns, width, toKeep, sdev, protect, displace));

    }while(yesno("Add another layer?"));
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

vector<Edge> random_connected(vector<Edge> in, double toKeep, bool protect){
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

//Given a set of vectors to nearest neighbors, find the smallest length
//among these vectors
double get_min_dist(map<int, vector<vector<double>>> nns){

    double dist_sq, min_dist_sq;
    min_dist_sq = FLT_MAX;

    for(auto iter = nns.begin(); iter != nns.end(); iter ++){
        for(vector<double> disp : iter->second){
            dist_sq = disp[0]*disp[0] + disp[1]*disp[1];
            if(dist_sq < min_dist_sq) min_dist_sq = dist_sq;
        }
    }

    return sqrt(min_dist_sq);
}

/*After the edges at one length scale that fall within the edges of the next
larger length scale have been determined, this function stitches together the 
network at the smaller length scale, given a knowledge of edges at the small 
length scale contained within a polygonal tiling of edges at the large length 
scale.
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

//Add edge indices to the data structure intiated with a call to make_poly_grid
void make_edge_grid(vector<Edge> edgeList, vector<vector<int>> &grid_lists, double length, double minx, double miny, int xdim, int ydim){

    int edgeIndex, xStart, yStart, xEnd, yEnd, xIter, yIter;
    Edge curr;

    //Iterate over the list of edges, adding an edge's index to each cell in
    //a set of cells such that one of the edge's points is in either the bottom
    //left or bottom right cell, and the other point is in the diagonally
    //opposite cell
    for(edgeIndex = 0; edgeIndex < edgeList.size(); edgeIndex ++){
        curr = edgeList[edgeIndex];
        xStart  = (int) floor((curr.left - minx) / length);
        xEnd = (int) floor((curr.right - minx) / length);
        yStart  = (int) floor((curr.bottom - miny) / length);
        yEnd = (int) floor((curr.top - miny) / length);
 
        for(yIter = yStart; yIter <= yEnd; yIter ++){
            for(xIter = xStart; xIter <= xEnd; xIter ++){
                grid_lists[yIter * xdim + xIter].push_back(edgeIndex);
            }
        }
    }
}

//Create a data structure mapping grid cells to lists of polygon indices
//partially contained in those grid cells, as well as lists of edge indices
void make_poly_grid(vector<PolyData> pdata, vector<vector<int>> &grid_lists, double length, double &minx, double &miny, double &maxx, double &maxy, int &xdim, int &ydim){

    int data_iter, xStart, yStart, xEnd, yEnd, xIter, yIter;
    minx = miny = (double) FLT_MAX;
    maxx = maxy = (double) FLT_MIN;
    vector<int>::iterator front;

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
        grid_lists[data_iter].push_back(-1);
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
                front = grid_lists[yIter * xdim + xIter].begin();
                grid_lists[yIter * xdim + xIter].insert(front, data_iter);
            }
        }
    }
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

/*This function is called when the edges describing the network at a given
length scale have been given non-zero width and stored in the collection top,
and the resulting geometrical object has been broken into a tiling of polygons
represented by the members of the collection pdat. The edges describing the 
network at the next lowest length scale, stored in the collection bottom, are
examined to determine whether they are contained entirely in one of these
polygons. Those edges meeting this criterion are kept and returned in the
vector retain.
*/
vector<Edge> sieve_edges(vector<Edge> top, double length, vector<Edge> bottom, vector<PolyData> pdat, double toKeep){

    vector<Edge> allEdges;
    vector<vector<Edge>> edge_collection;
    map<int, vector<Edge>> adj_map;
    set<Point> inPoints;
    bool p1in, p2in, p1here, p2here, cross, inpoly;
    int intersections, pdat_index, pdat_start = 0;
    int  p1_index, p2_index, map_index, num_polys;
    int xdim, ydim, xStart, xEnd, yStart, yEnd, list_iter, xIter, yIter;
    double minx, miny, maxx, maxy;
    vector<vector<int>> grid_lists;
    vector<int> currList;
    Edge large_edge;


    //Sort edges and polygons by vertical position to narrow search
    sort(top.begin(), top.end());
    sort(bottom.begin(), bottom.end());
    sort(pdat.begin(), pdat.end());

    for(int i = 0; i < pdat.size(); i ++){
        edge_collection.push_back(vector<Edge>());
    }

    //Set up grid data structure to determine in which polygons a point may
    //lie and which large-scale edges a small-scale edge may cross
    make_poly_grid(pdat, grid_lists, length, minx, miny, maxx, maxy, xdim, ydim);
    make_edge_grid(top, grid_lists, length, minx, miny, xdim, ydim);


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
                    list_iter = 0;

                    while((pdat_index = currList[list_iter++]) >= 0){

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
                            if(toKeep < 1) edge_collection[pdat_index].push_back(nextEdge);
                            else allEdges.push_back(nextEdge);
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
                        currList = grid_lists[yIter * xdim + xIter];
                        list_iter = -1;
                        while(currList[list_iter++] >= 0);

                        while(list_iter < currList.size()){
                            large_edge = top[currList[list_iter]];
                            if(intersection(nextEdge, large_edge)){
                                intersections ++;
                                if(intersections == 2) break;
                            }
                            list_iter ++;
                        }
                    }
                }

                if(intersections < 2){
                    if(toKeep < 1){
                        map_index = adj_map_index(p1_index, p2_index, num_polys);
                        if(adj_map.find(map_index) == adj_map.end()){
                            adj_map.insert(make_pair(map_index, vector<Edge>()));
                        }
                        adj_map[map_index].push_back(nextEdge);
                    }
                    else allEdges.push_back(nextEdge);
            
                }
            }
        }
    }


    if(toKeep < 1){
        return stitch_network(edge_collection, adj_map, inPoints, toKeep);
    }

    else{
        return allEdges;
    }
}

//Find which edge in a given edge crosses. If the edge crosses none of the
//edges in the set, return the edge itself.
Edge get_cross_edge(Edge input, vector<Edge> target_set){

    for(Edge next : target_set){
        if(intersection(input, next)) return next;
    }

    return input;
}

//Given a set of polygonal tiles, find all edges common to two tiles and
//arrange them in a map for the purpose of stitching together networks
//contained in different tiles
unordered_map<Edge, vector<Edge>> get_border_map(vector<PolyData> pdat){
    unordered_map<Edge, int> count_map;
    unordered_map<Edge, vector<Edge>> border_edges;

    for(PolyData data : pdat){
        for(Edge e : data.polyEdges){
            if(count_map.find(e) == count_map.end()){
                count_map.insert(make_pair(e, 1));
            }

            else count_map[e] = count_map[e] + 1;
        }
    }

    for(auto iter = count_map.begin(); iter != count_map.end(); iter++){
        if(iter->second >= 2){
            border_edges.insert(make_pair(iter->first,vector<Edge>()));
        }
    }

    return border_edges;
}

/*
Given a set of tiles describing a network with disorder in node position,
fill tiles with small-scale networks, then stitch tiles together.
*/
vector<Edge> make_edges_deformed(vector<vector<double>> rules, map<int, vector<vector<double>>> nns, vector<Edge> tedges, vector<PolyData> pdat, double width, double toKeep){

    double hwidth = width / 2;
    //Map edges defining grain boundaries to points on either side of boundaries
    unordered_map<Edge, vector<Edge>> global_map, local_map;
    Edge cross_edge;
    //Map PolyData objects defining grain boundaries to the bonds within a grain
    //map<PolyData, vector<Edge>> grain_map;


    //Bounds for creating a set of points and edges within a grain
    double h_offset = get_h_offset(rules);
    double starting_bounds[] = {-h_offset,0,0,width};
    vector<double> bounds(starting_bounds, starting_bounds + 4);
    double length, angle, dx, dy, distance, min;
    Point midpoint, p1, p2, closest, pivot;
    int index, total_edges, num_needed;
    vector<Edge> large, small, curr, stitch_edges;
    vector<Edge> retain, final_network, pool, mst;
    double left, right, low, high; 

    //Find borders between neighboring tiles
    global_map = get_border_map(pdat);

    //cout << "The size: " << tedges.size() << "\n";

    //Create a set of nodes and edges for each tile, retaining interior points 
    //from edges that cross borders with neighboring tiles
    for(index = 0; index < tedges.size(); index++){
        /*length = tedges[index].length();
        midpoint = tedges[index].midpoint();
        bounds[2] = length;*/
        p1 = tedges[index].p1;
        p2 = tedges[index].p2;
        angle = atan2(p2.y - p1.y, p2.x - p1.x);
        PolyData grain = pdat[index];
        vector<Edge> copy;
        copy.insert(copy.begin(),grain.polyEdges.begin(),grain.polyEdges.end());
        /*printf("%4.8lf\t%4.8lf\n%4.8lf\t%4.8f\n\n", p1.x, p1.y, p2.x, p2.y);
        for(Edge pdEdge : grain.polyEdges){
            printf("%4.8lf\t%4.8lf\n%4.8lf\t%4.8f\n\n", pdEdge.p1.x, pdEdge.p1.y, pdEdge.p2.x, pdEdge.p2.y);
        }*/
        pivot = copy[0].p1;
        rotate_edges(copy, -angle, pivot);
        /*for(Edge pdEdge : copy){
            printf("%4.8lf\t%4.8lf\n%4.8lf\t%4.8f\n\n", pdEdge.p1.x, pdEdge.p1.y, pdEdge.p2.x, pdEdge.p2.y);
        }*/
        get_extremes(copy, left, low, right, high); 
        length = right - left;
        midpoint = Point((left + right)/2, (high+low)/2);
        /*cout << "Midpoint: " << left << "\t" << low << "\n";
        cout << "Midpoint: " << right << "\t" << high << "\n";*/
        midpoint = rotate_point(midpoint, angle, pivot);
        //cout << midpoint.x << "\t" << midpoint.y << "\n";
        bounds[2] = length + h_offset;
        vector<Edge> grain_edges = makeedges(rules, nns, bounds);
        get_extremes(grain_edges, left, low, right, high);
        displace(grain_edges, midpoint.x - (left+right)/2, midpoint.y - (high+low)/2);
        rotate_edges(grain_edges, angle, midpoint);
        

        //cout << "Still alive here.\n";
        //final_network.insert(final_network.end(), grain_edges.begin(), grain_edges.end());

        //Build local map from border edges to grain edges that cross border
        //edges
        for(Edge gEdge : grain.polyEdges){
            if(global_map.find(gEdge) != global_map.end()){
                local_map.insert(make_pair(gEdge, vector<Edge>()));
            }
        }

        for(Edge e : grain_edges){
            //Identify edges entirely within a grain
            if(grain.contains(e.p1) && grain.contains(e.p2)){
                retain.push_back(e);
            }

            //Identify edges with one point in the current grain and the other
            //point in a neighboring grain
            else if(grain.contains(e.p1) || grain.contains(e.p2)){
                //Find the grain edge the current edge intersects, and determine
                //whether that grain edge is shared by a neighboring edge
                cross_edge = get_cross_edge(e, grain.polyEdges);
                if(local_map.find(cross_edge) != local_map.end()){
                    if(! grain.contains(e.p1)) e.reset(e.p2, e.p1);
                    local_map[cross_edge].push_back(e);
                }
            }
        }

        for(auto iter = local_map.begin(); iter != local_map.end(); iter++){
            if(global_map[iter->first].size() == 0){
                curr = iter->second;
                global_map[iter->first].insert(global_map[iter->first].end(), curr.begin(), curr.end());
            }

            else{
                if(iter->second.size() > global_map[iter->first].size()){
                    large = iter->second;
                    small = global_map[iter->first];
                }
                else{
                    small = iter->second;
                    large = global_map[iter->first];
                }

                for(Edge small_edge : small){
                    min = FLT_MAX;
                    for(Edge large_edge : large){
                        dx = large_edge.p1.x - small_edge.p1.x;
                        dy = large_edge.p1.y - small_edge.p1.y;
                        distance = sqrt(dx*dx + dy*dy);
                        if(distance < min){
                            min = distance;
                            closest = large_edge.p1;
                        }
                    }
                    stitch_edges.push_back(Edge(small_edge.p1, closest));
                }
                global_map.erase(iter->first);
            }
            
            if(toKeep < 1 && stitch_edges.size() > 0){
                shuffle_edges(stitch_edges);
                final_network.push_back(stitch_edges[0]);
                pool.insert(pool.end(), stitch_edges.begin()+1,stitch_edges.end());
            }
            else final_network.insert(final_network.end(), stitch_edges.begin(),stitch_edges.end());

            stitch_edges.clear();
        }

        if(toKeep < 1){
            mst = randomMST(retain, pool);
            final_network.insert(final_network.end(), mst.begin(), mst.end());
            mst.clear();
        }
        else final_network.insert(final_network.end(), retain.begin(), retain.end());

        local_map.clear();
        retain.clear();

        //final_network.insert(final_network.end(), retain.begin(), retain.end());
    }

    if(toKeep < 1){
        total_edges = final_network.size() + pool.size();
        num_needed = (int) (toKeep * total_edges) - final_network.size();
        shuffle_edges(pool);
        if(num_needed > 0){
            final_network.insert(final_network.end(), pool.begin(), pool.begin() + num_needed);
        }
        else if(num_needed < 0){
            cerr << "Could not make a connected network with bond portion " << toKeep << ".\n";
        }
    }

    return final_network;
}

//Group points according to the large-scale bonds in which they lie. Also find
//all edges within a given large-scale bond, and produce lists of edges
//on the small scale that span large-scale bonds.
void sort_edges(double length, vector<Edge> bottom, vector<PolyData> pdat, vector<vector<Edge>> &edge_collection){

    bool p1here, p2here, inpoly;
    int  p1_index, p2_index, map_index, num_polys, list_iter;
    int xdim, ydim, xStart, xEnd, yStart, yEnd, xIter, yIter, pdat_index;
    double minx, miny, maxx, maxy;
    vector<vector<int>> grid_lists;
    vector<int> currList;

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

                    list_iter = 0;
                    while((pdat_index = currList[list_iter++]) >= 0){

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
                    edge_collection[p1_index].push_back(nextEdge);
                    edge_collection[p2_index].push_back(nextEdge);
                }
            }
        }
    }
}

//Sort small edges made from a random packing into large-scale tiles
void sort_random_edges(double length, vector<Edge> bottom, vector<PolyData> pdat, vector<vector<Edge>> &edge_collection){

    bool inpoly;
    int  p1_index, p2_index, map_index, num_polys, list_iter;
    int xdim, ydim, xStart, xEnd, yStart, yEnd, xIter, yIter, pdat_index;
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

                    list_iter = 0;
                    while((pdat_index = currList[list_iter++]) >= 0){

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
                    edge_collection[p1_index].push_back(nextEdge);
                    edge_collection[p2_index].push_back(nextEdge);
                }
            }
            else if(p1_index >= 0) edge_collection[p1_index].push_back(nextEdge);
            else if(p2_index >= 0) edge_collection[p2_index].push_back(nextEdge);
        }
    }
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

Point find_match(Point p1, Point p2, double y){
    if(p1.y == y) return p1;
    else return p2;
}

void prepare_grips_displaced(vector<Edge> &edge_list, bool connect_top_bottom){
    double low, high;
    vector<double> ycuts;
    bool valid, curr_defined;
    set<Point> bps, tps;
    Point p1, p2, new_point, curr;
    vector<Edge> replace;
    double intersect_x;
    Edge next;

    do{
        ycuts = getdoubles("Enter the heights for the lower and upper grips: ");

        if(ycuts.size() < 2){
            cerr << "Enter two numbers.\n";
            continue;
        }

        low = ycuts[0];
        high = ycuts[1];

        if(high < low){
            cerr << "The upper bound must not be less than the lower bound.\n";
        }high = ycuts[1];

        if(high < low){
            cerr << "The upper bound must not be less than the lower bound.\n";
        }

        else valid = true;
    }while(! valid);

    while(!edge_list.empty()){
        next = edge_list.front();
        p1 = next.p1;
        p2 = next.p2;
        edge_list.erase(edge_list.begin());
        if(next.bottom >= low && next.top <= high){
            replace.push_back(next);
            if(next.bottom == low) bps.insert(find_match(p1, p2, low));
            if(next.top == high) tps.insert(find_match(p1, p2, high));
        }

        else if(next.bottom < low && next.top >= low){
            intersect_x = x_intersect(next, low);
            new_point = Point(intersect_x, low);
            replace.push_back(Edge(new_point, find_match(p1, p2, next.top)));
            bps.insert(new_point);
        }
        
        else if(next.bottom <= high && next.top > high){
            intersect_x = x_intersect(next, high);
            new_point = Point(intersect_x, high);
            replace.push_back(Edge(find_match(p1, p2, next.bottom), new_point));
            tps.insert(new_point);
        }
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

    edge_list = replace;
}

void prepare_grips_simple(vector<Edge> &edge_list){
    unordered_set<Edge> low_edges, high_edges;
    set<Point> low_points, high_points;
    double miny, maxy;
    Edge proposed;

    if(edge_list.size() == 0) return;

    sort(edge_list.begin(), edge_list.end(), sort_by_bottom);
    miny = edge_list[0].bottom;
    for(Edge next_edge : edge_list){
        if(next_edge.bottom > miny) break;
        low_edges.insert(next_edge);
        if(next_edge.p1.y == miny) low_points.insert(next_edge.p1); 
        if(next_edge.p2.y == miny) low_points.insert(next_edge.p2); 
    }

    sort(edge_list.begin(), edge_list.end(), sort_by_top);
    maxy = (*edge_list.rbegin()).top;
    for(auto eiter = edge_list.rbegin(); eiter != edge_list.rend(); eiter ++){
        if((*eiter).top < miny) break;
        high_edges.insert(*eiter);
        if((*eiter).p1.y == maxy) high_points.insert((*eiter).p1); 
        if((*eiter).p2.y == maxy) high_points.insert((*eiter).p2); 
    }

    for(auto piter = low_points.begin(); piter != prev(low_points.end()); piter++){
        proposed = Edge(*piter, *(next(piter)));
        if(low_edges.find(proposed) == low_edges.end()){
            edge_list.push_back(proposed);
        }
    }

    for(auto piter = high_points.begin();piter != prev(high_points.end()); piter++){
        proposed = Edge(*piter, *(next(piter)));
        if(high_edges.find(proposed) == high_edges.end()){
            edge_list.push_back(proposed);
        }
    }
}

void prepare_grips(vector<Edge> &edge_list, bool displaced, bool connect){
    if(! displaced) prepare_grips_simple(edge_list);
    else prepare_grips_displaced(edge_list, connect);
}

void adjust_bounds(vector<double> &bounds, vector<vector<double>> rules, vector<Edge> edges){
    double hoffset, voffset, minx, miny, maxx, maxy;

    get_extremes(edges, minx, miny, maxx, maxy);
    hoffset = get_h_offset(rules);
    voffset = get_v_offset(rules);

    if(minx < bounds[0]){
        bounds[0] -= ceil((bounds[0] - minx) / hoffset) * hoffset;
    }
    if(miny < bounds[1]){
        bounds[1] -= ceil((bounds[1] - miny) / voffset) * voffset;
    }
    if(maxx > bounds[2]){
        bounds[2] += ceil((maxx - bounds[2]) / hoffset) * hoffset;
    }
    if(maxy > bounds[3]){
        bounds[3] += ceil((maxy - bounds[3]) / voffset) * voffset;
    }
}

vector<Edge> edge_hierarchy(vector<double> bounds, int polyflag, bool getAlign, bool connect){
    stack<NetData> netStack;
    vector<vector<Edge>> edge_collection;
    vector<Edge> tedges, bedges, backup;
    vector<PolyData> pdat;
    bool displacement = false;
    int iter;
    double length, alignment;
    FILE *align_report;
    string response;

    loadNetStack(netStack);
    NetData tdat = netStack.top();
    length = get_min_dist(tdat.nns);
    netStack.pop();
    tedges = makeedges(tdat.rules, tdat.nns, bounds);

    if(tdat.to_keep < 1){
        tedges = random_connected(tedges, tdat.to_keep, tdat.protect);
    }

    if(tdat.displace){
        displace_points_grn(tedges, tdat.sdev);
        displacement = true;
    }

    if(yesno("Report skeleton?")){
        print_edges(tedges, "Enter a file name or enter to decline: ");
    }

    if(netStack.empty() && yesno("Prepare for grips?")){
        prepare_grips(tedges, displacement, connect);
    }

    if(tdat.width > 0){
        if(tdat.displace || getAlign) backup = tedges;
        tedges = add_thickness(tedges, tdat.width, pdat, true, !tdat.displace);
    }

    if(polyflag == 1){
        for(iter = 0; iter < pdat.size(); iter++){
            if(displace){
                fprintf(stderr, "%lf\t%lf\n",backup[iter].p1.x,backup[iter].p1.y);
                fprintf(stderr, "%lf\t%lf\n\n",backup[iter].p2.x,backup[iter].p2.y);
            }
        //for(PolyData pdatum : pdat){
            PolyData pdatum = pdat[iter];
            for(Edge pedge : pdatum.polyEdges){
                fprintf(stderr, "%lf\t%lf\n", pedge.p1.x, pedge.p1.y);
                fprintf(stderr, "%lf\t%lf\n\n", pedge.p2.x, pedge.p2.y);
            }
        }
    }

    if(yesno("Report top network?")){
        print_edges(tedges, "Enter a file name or enter to decline: ");
    }
    
    while(!netStack.empty()){
        NetData bdat = netStack.top();
        netStack.pop();

        adjust_bounds(bounds, bdat.rules, tedges);
        if(displacement && yesno("Use grain-based approach?")){
            bedges = make_edges_deformed(bdat.rules, bdat.nns, backup, pdat, tdat.width, bdat.to_keep);
        }
        else{
            bedges = makeedges(bdat.rules, bdat.nns, bounds);
            //cout << "Check 1\n";
            bedges = sieve_edges(tedges, length, bedges, pdat, bdat.to_keep);
        }

        //cout << "Check2\n";

        if(bdat.displace){
            displace_points_grn(bedges, bdat.sdev);
            displacement = true;
        }

        //cout << "Check3\n";
        if(netStack.empty()){
            if(yesno("Prepare for grips?")) prepare_grips(bedges, displacement, connect);

            if(getAlign){

                cout << "Enter the name for the alignment report: ";
                getline(cin, response);
                align_report = fopen(response.c_str(), "w");

                if(align_report != NULL){
                    if(displacement){
                        sort_random_edges(length, bedges, pdat,edge_collection);
                    }
                    else{
                        sort_edges(length, bedges, pdat, edge_collection);
                    }

                    for(iter = 0; iter < edge_collection.size(); iter++){
                        alignment  = calc_alignment(edge_collection[iter], backup[iter]);
                        fprintf(align_report, "%lf\t%ld\n", alignment, edge_collection[iter].size());
                    }
                    fclose(align_report);
                }
            }
        }

        if(bdat.width > 0){
            pdat.clear();
            if(!netStack.empty() && (displace || getAlign)) backup = bedges;
            bedges = add_thickness(bedges, bdat.width, pdat, true, !bdat.displace);
        }

        //cout << "Check4\n";
        /*for(PolyData pdatum : pdat){
            for(Edge pedge : pdatum.polyEdges){
                fprintf(stderr, "%lf\t%lf\n", pedge.p1.x, pedge.p1.y);
                fprintf(stderr, "%lf\t%lf\n", pedge.p2.x, pedge.p2.y);
            }
        }*/

        tedges = bedges;
        length = get_min_dist(bdat.nns);
        //cout << "Check5\n";
    }

    cout << "There are " << tedges.size() << " edges.\n";
    return tedges;
}

int main(int argc, char **argv){
    vector<vector<double>> rules;
    map<int,vector<vector<double>>> nns;
    vector<double> rule, bounds, nn;
    set<Point> points;
    vector<Edge> edges;
    string nextline;
    int ruleiter, pointiter, count, numRead;
    string filename;
    FILE *out, *ranfile;
    double width;
    bool flag, success = false, align = false, connect = true;
    int c, polyflag = 0;
    unsigned seed;
    char num_string[11];

    opterr = 0;

    while((c = getopt(argc, argv, "adp")) != -1){
        switch(c) {
            case 'a':
                align = true;
                break;
            case 'd':
                connect = false;
                break;
            case 'p':
                polyflag = 1;
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

    //Read in network bounds
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

    ranfile = fopen("/dev/urandom", "r");
    fread(&seed, sizeof(unsigned), 1, ranfile);
    fclose(ranfile);
    sprintf(num_string, "%u", seed);
    setenv("GSL_RNG_SEED", num_string, 1);
    gsl_rng_env_setup();

    edges = edge_hierarchy(bounds, polyflag, align, connect);

    print_edges(edges, "Enter a file name for output or enter to decline: ");

    return 0;
}
