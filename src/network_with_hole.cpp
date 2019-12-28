/*This code produces a triangular lattice with a triangular hole of a given
size, centered within the lattice.
*/

#include <iostream>
#include <cmath>
#include "network_utils.h"
#include <string>
#include "stdio.h"
#include "stdlib.h"
#include "string.h"


#define sin60 0.866025404
#define cos60 0.5

void update_plist(set<Point> &points, Point &p){
    auto iter = points.find(p);
    if(iter != points.end()){
        p.x = iter->x;
        p.y = iter->y;
    }
    else points.insert(p);
}

//Make a triangular lattice with a polygon removed from it
vector<Edge> generate_tri_lat(double minx, double miny, double maxx, double maxy, double scale, PolyData hole){

    double offset = 0, x1, y1, x2, y2, dx = scale*cos60, dy = scale*sin60;
    vector<Edge> lat_edges;
    set<Point> lat_points;
    Point p1, p2, mid;

    for(y1 = miny; y1 <= maxy + FLOAT_TOL; y1 += dy){
        for(x1 = minx + offset; x1 <= maxx; x1 += scale){
            p1 = Point(x1, y1);
  
            if(hole.contains(p1)) continue;

            update_plist(lat_points, p1);
            x2 = x1 - dx;
            y2 = y1 + dy;
            if(x2 >= minx && y2 <= maxy + FLOAT_TOL){
                p2 = Point(x2,y2);
                mid = Point((x1 + x2)/2, (y1 + y2)/2);
                if(! (hole.contains(p2) || hole.contains(mid))){
                    update_plist(lat_points, p2);
                    lat_edges.push_back(Edge(p1, p2));
                }
            }
            x2 = x1 + dx;
            if(x2 <= maxx && y2 <= maxy + FLOAT_TOL){
                p2 = Point(x2,y2);
                mid = Point((x1 + x2)/2, (y1 + y2)/2);
                if(! (hole.contains(p2) || hole.contains(mid))){
                    update_plist(lat_points, p2);
                    lat_edges.push_back(Edge(p1, p2));
                }
            }
            x2 = x1 + scale;
            if(x2 <= maxx){
                p2 = Point(x2,y1);
                mid = Point((x1 + x2)/2, y1);
                if(! (hole.contains(p2) || hole.contains(mid))){
                    update_plist(lat_points, p2);
                    lat_edges.push_back(Edge(p1, p2));
                }
            }
        }
        offset = offset == 0 ? scale / 2 : 0;
    }

    return lat_edges;
}

int main(int argc, char **argv){

    double minx, miny, maxx, maxy, blen, xrange, yrange, xleft, xmid, xright;
    double ylower, yupper;
    int slen, num_read;
    Point p1, p2, p3;
    vector<Edge> poly_edges, final_network;
    string response, name;
    FILE *report_file;

    //Obtain the bounds for the network
    while(true){
        cout << "Enter the bounds for the network: ";
        getline(cin, response);
        num_read = sscanf(response.c_str(), "%lf %lf %lf %lf", &minx, &miny, &maxx, &maxy);
        if(num_read == 4 && minx < maxx && miny < maxy) break;
        else cerr << "The input was invalid.\n";
    }

    //Obtain the length of an individual bond
    while(true){
        cout << "Enter the bond length: ";
        getline(cin, response);
        num_read = sscanf(response.c_str(), "%lf", &blen);
        if(num_read == 1 && blen > 0) break;
        else cerr << "Enter a positive number.\n";
    }

    //Find the side length of the hole, in as a multiple of the bond length
    while(true){
        cout << "Enter the hole side length, as a multiple of bond length: ";
        getline(cin, response);
        num_read = sscanf(response.c_str(), "%d", &slen);
        if(num_read == 1 && slen > 0) break;
        else cerr << "Enter a positive integer.\n";
    }

    //Find the coordinates of the hole's vertices, and define the hole as a
    //PolyData object

    //First obtain x and y coordinates
    xrange = maxx - minx;
    yrange = maxy - miny;
    xleft = floor((xrange - slen*blen)/(2*blen)) * blen + minx + blen / 10;
    ylower = floor((yrange-slen*blen*sin60)/(4*blen*sin60))*2*blen*sin60 + miny + blen / 10;
    xright = xleft + blen * slen - blen / 5;
    xmid = xleft + blen * slen / 2 - blen / 20;
    yupper = ylower + blen * slen * sin60 - blen / 5;

    //Next, construct the edges and create the polygon data structure
    poly_edges.push_back(Edge(Point(xleft, ylower), Point(xright, ylower)));
    poly_edges.push_back(Edge(Point(xright, ylower), Point(xmid, yupper)));
    poly_edges.push_back(Edge(Point(xmid, yupper), Point(xleft, ylower)));
    PolyData hole = PolyData(poly_edges);

    if(yesno("Report hole?")){
        cout << "Enter the file name: ";
        getline(cin, name);
        report_file = fopen(name.c_str(), "w");
        for(Edge e : poly_edges){
            fprintf(report_file, "%2.8lf %2.8lf\n", e.p1.x, e.p1.y);
            fprintf(report_file, "%2.8lf %2.8lf\n\n", e.p2.x, e.p2.y);
        }
        fclose(report_file);
    }

    //Use information about the bonds, side length and hole to create the 
    //network
    final_network = generate_tri_lat(minx, miny, maxx, maxy, blen, hole);

    //Report the network to a specified file
    while(true){
        cout << "Enter the file name for the network: ";
        getline(cin, name);
        if(name.compare("") != 0) break;
        else cerr << "Enter a name.\n";
    }

    report_file = fopen(name.c_str(), "w");
    for(Edge e : final_network){
        fprintf(report_file, "%2.8lf %2.8lf\n", e.p1.x, e.p1.y);
        fprintf(report_file, "%2.8lf %2.8lf\n\n", e.p2.x, e.p2.y);
    }

    fclose(report_file);

    return 0;
}
