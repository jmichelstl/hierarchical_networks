#include "network_utils.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <cmath>

using namespace std;

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
            }
            else{
                p2 = Point(edge_data[0], edge_data[1]);
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

void read_relaxed(ifstream& datfile, vector<double>& point_list, map<int, vector<int>>& neighbor_map){

    map<Point, int> pmap;
    string nextline;
    Point p1, p2;
    int pindex = 0, point_count = 0, index1, index2, mindex, maxdex;
    vector<double> edge_data;

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
                    neighbor_map.insert(make_pair(pindex, vector<int>()));
                    pindex++;
                }
                if(pmap.find(p2) == pmap.end()){
                    pmap.insert(make_pair(p2, pindex));
                    point_list.push_back(p2.x);
                    point_list.push_back(p2.y);
                    neighbor_map.insert(make_pair(pindex, vector<int>()));
                    pindex++;
                }

                index1 = pmap[p1];
                index2 = pmap[p2];
                mindex = index1 < index2 ? index1 : index2;
                maxdex = index1 == mindex ? index2 : index1;

                neighbor_map[mindex].push_back(maxdex);
            }
        }
    }

    datfile.close();
}

//Given a vector of displacements and a file giving the eigen vectors of a
//stiffness matrix, calculate the projection of the displacement vector along
//each normal mode, normalized by the overall magnitude of the displacement
//vector
int proj_calc(vector<double> &projections, vector<double> disps, string vec_name){
    double disp_mag = 0, dprod;
    int dim, size, target_size, num_read, vec_iter, disp_iter;
    double *vec;
    FILE *vec_file = NULL;

    vec_file = fopen(vec_name.c_str(), "rb");
    if(vec_file == NULL){
        cerr << "The eigen vector file could not be read.\n";
        return -1;
    }

    //Read the leading dimension of the stiffness matrix
    num_read = fread(&dim, sizeof(int), 1, vec_file);
    if(num_read == 0){
        cerr << "The dimension of the vectors could not be read.\n";
        fclose(vec_file);
        return -1;
    }

    //Confirm that the displacement vector has the same dimension as the eigen
    //vectors
    if(dim != disps.size()){
        cerr << "The dimensions of the displacement vector and eigen vectors do not match.\n";
        return -1;
        fclose(vec_file);
    }

    //Confirm that there is enough data to read all eigen vectors
    fseek(vec_file, 0, SEEK_END);
    size = ftell(vec_file) - sizeof(int);
    target_size = sizeof(double) * dim * dim;
    if(size < target_size){
        cerr << "There was too little data in the file.\n";
        fclose(vec_file);
        return -1;
    }

    vec = (double *) malloc(sizeof(double) * dim);

    //Find the magnitude of the displacement vector
    for(double dcomp : disps){
        disp_mag += dcomp * dcomp;
    }
    disp_mag = sqrt(disp_mag);


    //Set the position indicator for the data file just past the lead integer
    fseek(vec_file, sizeof(int), SEEK_SET);

    //Take the dot product of the displacement vector with each eigen vector
    for(vec_iter = 0; vec_iter < dim; vec_iter ++){
        num_read = fread(vec, sizeof(double), dim, vec_file);
        if(num_read < dim){
            cerr << "There was an error reading vector " << vec_iter+1 << ".\n";
            fclose(vec_file);
            return -1;
        }

        dprod = 0;

        for(disp_iter = 0; disp_iter < dim; disp_iter ++){
            dprod += disps[disp_iter] * vec[disp_iter];
        }
        dprod /= disp_mag;
        projections.push_back(abs(dprod));
    }

    free(vec);
    fclose(vec_file);
    return 0;
}

int main(int argc, char **argv){
    vector<double> point_list, projections, efreqs;
    map<int, vector<int>> neighbor_map;
    ifstream unstrained_file, strained_file, nm_file;
    ofstream proj_file;
    string unstrained_name, strained_name, vector_name, proj_name, nm_name;
    string nextline;
    vector<Edge> strained_edges;
    int edge_list_index = 0, ref_index, num_read, dindex;
    Point p1, p2;
    double efreq;
    //double init_miny, init_maxy, final_miny, final_maxy, total_strain;
    //double x_disp, y_disp;
    //double relaxed, strained, strain_mag, affine_disp, non_affine_par;
    //int init_edge_count = 0;

    if(argc < 4){
        cerr << "Usage: original network, final network, eigen vector file, projection file, optional eigen frequency file\n";
        return -1;
    }

    unstrained_name = argv[1];
    strained_name = argv[2];
    vector_name = argv[3];
    proj_name = argv[4];
    if(argc >= 5) nm_name = argv[5];

    unstrained_file.open(unstrained_name);
    strained_file.open(strained_name);
    if(! (unstrained_file.is_open() && strained_file.is_open())){
        cerr << "Not all input files could be opened.\n";
        if(unstrained_file.is_open()) unstrained_file.close();
        if(strained_file.is_open()) strained_file.close();
        return 1;
    }

    read_relaxed(unstrained_file, point_list, neighbor_map);
    strained_edges = read_edges(strained_file);

    vector<double> disp_vec(point_list.size(), 0);

    //Find the displacement of each point from its resting position
    for(auto iter = neighbor_map.begin(); iter != neighbor_map.end(); iter++){
        ref_index = iter->first;
        p1 = strained_edges[edge_list_index].p1;

        if(disp_vec[2 * ref_index] == 0){
            disp_vec[2 * ref_index] = p1.x - point_list[2 * ref_index];
            disp_vec[2 * ref_index + 1] = p1.y - point_list[2 * ref_index + 1];
        }

        for(int neighbor : iter->second){
            p2 = strained_edges[edge_list_index].p2;

            if(disp_vec[2*neighbor] == 0){
                disp_vec[2*neighbor] = p2.x - point_list[2*neighbor];
                disp_vec[2*neighbor + 1] = p2.y - point_list[2*neighbor + 1];
            }

            edge_list_index ++;
        }
    }

    //Once displacement vectors are found, project them along the stiffness
    //matrix eigen basis
    proj_calc(projections, disp_vec, vector_name); 

    //Load a list of eigen frequencies if an eigen frequency file has been
    //defined
    if(argc >= 5){
        nm_file.open(nm_name);
        if(nm_file.is_open()){
            while(! nm_file.eof()){
                getline(nm_file, nextline);
                num_read = sscanf(nextline.c_str(), "%lf", &efreq);
                if(num_read == 1) efreqs.push_back(efreq);
            }
            nm_file.close();
        }
    }

    //Report a list of frequency-projection pairs to the specified file
    proj_file.open(proj_name);
    proj_file.precision(10);

    if(efreqs.size() == projections.size()){
        for(dindex = 0; dindex < projections.size(); dindex++){
            proj_file << efreqs[dindex] << "\t" << projections[dindex] << "\n";
        }
    }

    else{
        for(double next_projection : projections){
            proj_file << next_projection << "\n";
        }
    }

    proj_file.close();    

    return 0;
}
