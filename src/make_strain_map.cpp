#include "network_utils.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <ctype.h>
#include <stdlib.h>
#include <unistd.h>

using namespace std;

struct EdgeDatum{

    EdgeDatum(int ival, double sval, double len) : index(ival), stiffness(sval), l0(len){
    }

    EdgeDatum(double ival, double len) : index(ival), l0(len){}

    int index;
    double stiffness;
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

vector<Edge> read_edges(ifstream &edge_file, double &miny, double &maxy){
    Point p1, p2;
    int point_count = 0;
    vector<double> edge_data;
    string nextline;
    vector<Edge> edge_list;
    vector<char> delims = {' ', '\t', '\n'};

    miny = FLT_MAX;
    maxy = -FLT_MAX;

    while(!edge_file.eof()){
        getline(edge_file, nextline);
        edge_data = parse_doubles(split_string(nextline, delims));
        if(edge_data.size() >= 2){
            if(point_count == 0){
                p1 = Point(edge_data[0], edge_data[1]);
                miny = p1.y < miny - FLOAT_TOL ? p1.y : miny;
                maxy = p1.y > maxy + FLOAT_TOL ? p1.y : maxy;
            }
            else{
                p2 = Point(edge_data[0], edge_data[1]);
                miny = p2.y < miny - FLOAT_TOL ? p2.y : miny;
                maxy = p2.y > maxy + FLOAT_TOL ? p2.y : maxy;
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

void read_relax_style(ifstream& datfile, vector<double>& point_list, map<int, vector<EdgeDatum>>& neighbor_map, double &minx, double &maxx, double &miny, double &maxy){

    map<Point, int> pmap;
    string nextline;
    Point p1, p2;
    int pindex = 0, point_count = 0, index1, index2, mindex, maxdex;
    vector<double> edge_data;
    double length;

    minx = miny = FLT_MAX;
    maxx = maxy = -FLT_MAX;

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
                    pindex++;
                    minx = p1.x < minx - FLOAT_TOL ? p1.x : minx;
                    maxx = p1.x > maxx + FLOAT_TOL ? p1.x : maxx;
                    miny = p1.y < miny - FLOAT_TOL ? p1.y : miny;
                    maxy = p1.y > maxy + FLOAT_TOL ? p1.y : maxy;
                }
                if(pmap.find(p2) == pmap.end()){
                    pmap.insert(make_pair(p2, pindex));
                    point_list.push_back(p2.x);
                    point_list.push_back(p2.y);
                    neighbor_map.insert(make_pair(pindex, vector<EdgeDatum>()));
                    pindex++;
                    minx = p2.x < minx - FLOAT_TOL ? p2.x : minx;
                    maxx = p2.x > maxx + FLOAT_TOL ? p2.x : maxx;
                    miny = p2.y < miny - FLOAT_TOL ? p2.y : miny;
                    maxy = p2.y > maxy + FLOAT_TOL ? p2.y : maxy;
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

    datfile.close();
}

void report_non_affine_par(ofstream &na_file, Point p_i, Point p_f, double strain, double xmid, double miny, bool contract){
    double y_affine_disp, x_affine_disp, x_disp, y_disp, non_affine_par;

    x_affine_disp = contract ? strain * (xmid - p_i.x) * strain : 0;
    y_affine_disp = (p_i.y - miny) * strain;
    x_disp = p_f.x - p_i.x - x_affine_disp;
    y_disp = p_f.y - p_i.y - y_affine_disp;
    //non_affine_par = sqrt(x_disp*x_disp + y_disp*y_disp);
    non_affine_par = abs(y_disp);

    na_file << p_i.x << " " << p_i.y << " " << non_affine_par << "\n";
} 

double min(double val1, double val2){
    return val1 < val2 ? val1 : val2;
}

double max(double val1, double val2){
    return val1 > val2 ? val1 : val2;
}

int main(int argc, char **argv){
    vector<double> point_list;
    map<int, vector<EdgeDatum>> neighbor_map;
    ifstream unstrained_file, strained_file;
    ofstream strain_file, na_file;
    string unstrained_name, strained_name, strain_report_name, non_affine_name;
    vector<Edge> strained_edges;
    int edge_list_index = 0, ref_index, num_read;
    Point p_i, p1, p2;
    double bottom_cutoff = -FLT_MAX, top_cutoff = FLT_MAX;
    double init_miny, init_maxy, final_miny, final_maxy, total_strain;
    double minx, maxx, xmid, x_disp, y_disp;
    double relaxed, strained, strain_mag, affine_disp, non_affine_par;
    int init_edge_count = 0, offset = 0;
    bool contract = false, affine_only = false, smap_only = false;
    char c;

    /*if(argc < 4){
        cerr << "Usage: original network, final network, strain file, non-affine file\n";
        return 1;
    }*/

    while((c = getopt(argc, argv, "ab:cst:")) != -1){
        switch(c) {
            //Only report the non-affine parameter
            case 'a':
                affine_only = true;
                break;
            case 'b':
                num_read = sscanf(optarg, "%lf", &bottom_cutoff);
                if(num_read < 1){
                    cerr << "Option \"b\" requires a numerical argument.\n";
                }
                break;
            case 'c':
                contract = true;
                break;
            case 's':
                smap_only = true;
                break;
            case 't':
                num_read = sscanf(optarg, "%lf", &top_cutoff);
                if(num_read < 1){
                    cerr << "Option \"t\" requires a numerical argument.\n";
                }
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

    /*if(argc > 4) offset = 1;

    unstrained_name = argv[1 + offset];
    strained_name = argv[2 + offset];
    strain_report_name = argv[3 + offset];
    non_affine_name = argv[4 + offset];
    */

    cout << "Enter the name of the unstrained file: ";
    getline(cin, unstrained_name);
    cout << "Enter the name of the strained file: ";
    getline(cin, strained_name);
    if(! affine_only){
        cout << "Enter the name of the strain map file: ";
        getline(cin, strain_report_name);
    }
    if(! smap_only){
        cout << "Enter the name of the non-affine file: ";
        getline(cin, non_affine_name);
    }

    unstrained_file.open(unstrained_name);
    strained_file.open(strained_name);
    if(! (unstrained_file.is_open() && strained_file.is_open())){
        cerr << "Not all input files could be opened.\n";
        if(unstrained_file.is_open()) unstrained_file.close();
        if(strained_file.is_open()) strained_file.close();
        return 1;
    }

    read_relax_style(unstrained_file, point_list, neighbor_map, minx, maxx, init_miny, init_maxy);
    xmid = (maxx + minx) / 2;
    strained_edges = read_edges(strained_file, final_miny, final_maxy);

    final_miny = max(final_miny, bottom_cutoff);
    final_maxy = min(final_maxy, top_cutoff);

    vector<bool> na_param_written(neighbor_map.size(), false);

    cout << "Initial min and max y: " << init_miny << ", " << init_maxy << "\n";
    cout << "Final min and max y: " << final_miny << ", " << final_maxy << "\n";

    total_strain = (final_maxy - final_miny) - (init_maxy - init_miny);
    total_strain /= (init_maxy - init_miny);
    cout << "Total strain: " << total_strain << "\n";
    cout << "Point list size: " << point_list.size() << "\n";
    cout << "Neighbor map size: " << neighbor_map.size() << "\n";
    cout << "Final edges: " << strained_edges.size() << "\n";
    for(auto iter = neighbor_map.begin(); iter != neighbor_map.end(); iter++){
        init_edge_count += iter->second.size();
    }
    cout << "Initial edge count here: " << init_edge_count << "\n";

    if(! affine_only){
        strain_file.open(strain_report_name);
    }
    if(! smap_only) na_file.open(non_affine_name);

    for(auto iter = neighbor_map.begin(); iter != neighbor_map.end(); iter++){
        ref_index = iter->first;
        p1 = strained_edges[edge_list_index].p1;

        if(! na_param_written[ref_index] && !smap_only){
            p_i = Point(point_list[ref_index*2], point_list[ref_index*2+1]);
            report_non_affine_par(na_file, p_i, p1, total_strain, xmid, init_miny, contract);
            na_param_written[ref_index] = true;
        }

        for(EdgeDatum edat : iter->second){
            p2 = strained_edges[edge_list_index].p2;

            if(! na_param_written[edat.index] && !smap_only){
                p_i = Point(point_list[edat.index*2],point_list[edat.index*2+1]);
                report_non_affine_par(na_file, p_i, p2, total_strain, xmid, init_miny, contract);
                na_param_written[edat.index] = true;
            }

            if(! affine_only){
                relaxed = edat.l0;
                strained = sqrt((p2.x-p1.x)*(p2.x-p1.x)+(p2.y-p1.y)*(p2.y-p1.y));
                strain_mag = abs((strained-relaxed)/relaxed);
                strain_file << p1.x << " " << p1.y << " " << p2.x << " " << p2.y;
                strain_file << " " << strain_mag << "\n";
            }
            edge_list_index ++;
        }
    }

    if(! affine_only) strain_file.close();
    if(! smap_only) na_file.close();

    return 0;
}
