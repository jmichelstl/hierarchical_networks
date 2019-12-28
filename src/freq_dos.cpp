#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include "network_utils.h"
#include <stdio.h>
#include <iterator>

using namespace std;

int main(int argc, char **argv){

    int bin_size = 0, remainder, counter;
    vector<double> evals;
    map<double, int> eval_bins;
    double mean, next_eval, diff1, diff2, dos, prev_freq, correction;
    string filename, response, nextline;
    ifstream nmfile;
    ofstream outfile;
    map<double, int>::iterator it, prev_iter, next_iter;
    

    //Get the bin size
    do{
        cout << "Enter the bin size: ";
        getline(cin, response);
        if(sscanf(response.c_str(), "%i", &bin_size) == 1){
            if(bin_size > 0) break;
            else cerr << "The bin size must be greater than zero.\n";
        }
        else cerr << "No integer was read.\n";

        if(! yesno("Try again?")) exit(0);
    }while(bin_size <= 0);

    //Open the data file to be processed
    do{
        cout << "Enter the data file name: ";
        getline(cin, filename);
        nmfile.open(filename);
        if(! nmfile.is_open()){
            if(! yesno("The file could not be opened. Try again?")) exit(0);
        }
    }while(! nmfile.is_open());

    //Obtain eigen values
    while(! nmfile.eof()){
        getline(nmfile, nextline);
        if(sscanf(nextline.c_str(), "%lf", &next_eval) == 1){
            evals.push_back(next_eval);
        }
    }
    nmfile.close();

    //Create bins
    counter = 0;
    mean = next_eval;
    for(float eval : evals){
        mean += eval;
        counter ++;
        if(counter == bin_size){
            counter = 0;
            mean /= bin_size;
            if(eval_bins.find(mean) != eval_bins.end()){
                eval_bins[mean] += bin_size;
            }
            else eval_bins.insert(make_pair(mean, bin_size));
            mean = 0;
        }
    }

    if(counter != 0){
        mean /= counter;
        if(eval_bins.find(mean) != eval_bins.end()){
            eval_bins[mean] = eval_bins[mean] + counter;
        }
        else eval_bins.insert(make_pair(mean, counter));
    }

    if(eval_bins.size() < 2){
        cerr << "There are too few bins.\n";
        exit(0);
    }

    do{
        cout << "Enter a name for the output file: ";
        getline(cin, filename);
        if(filename.compare("") == 0) cerr << "Enter a name.\n";
        else{
            outfile.open(filename);
            if(!outfile.is_open()) cerr << "Opening " << filename <<"failed.\n";
        }
        if(! outfile.is_open()){
            if(! yesno("Try again?")) exit(0);
        }
    }while(! outfile.is_open());

    //Calculate DOS values and report to the output file
    it = eval_bins.begin();
    if(it != eval_bins.end()){
        next_iter = next(it, 1);
        dos = next_iter->second / (next_iter->first - it->first);
        outfile << it->first << "\t" << dos << "\n";
    }

    for(it = next(it,1); it != prev(eval_bins.end(), 1) ; it++){
        prev_iter = prev(it,1);
        next_iter = next(it, 1);
        diff1 = it->first - prev_iter->first;
        diff2 = next_iter->first - it->first;

        dos = next_iter->second * diff1 * diff1 + it->second * diff2 * diff2;
        dos /= diff1 * diff2*(diff1 + diff2);
        outfile << it->first << "\t" << dos << "\n";
    }

    prev_iter = prev(it, 1);
    dos = it->second / (it->first - prev_iter->first);
    outfile << it->first << "\t" << dos << "\n";
    outfile.close();

    return 0;
}
