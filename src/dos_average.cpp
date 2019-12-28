#include <stdio.h>
#include <iostream>
#include <fstream>
#include <dirent.h>
#include <string>
#include <iostream>
#include <regex>
#include "network_utils.h"
#include <map>
#include <algorithm>

using namespace std;

#define DEFAULT_EXT "nm"

map<double, double> get_dos(vector<double> freqs, int num_sets){
    int bin_size = 0, counter = 0;
    string response;
    map<double, double> dos_vs_freq;
    map<double, int> freq_bins;
    float mean, next_freq, diff1, diff2, dos, integral, correction;
    map<double, int>::iterator it, prev_iter, next_iter;
    map<double, double>::iterator next_fmap_iter;

    //Before bins are made, frequencies must be in ascending order
    sort(freqs.begin(), freqs.end());

    //Get the number of bins
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

    mean = 0;
    for(float freq : freqs){
        mean += freq;
        counter ++;
        if(counter == bin_size){
            counter = 0;
            mean /= bin_size;
            if(freq_bins.find(mean) != freq_bins.end()){
                freq_bins[mean] += bin_size;
            }
            else freq_bins.insert(make_pair(mean, bin_size));
            mean = 0;
        }
    }

    if(counter != 0){
        mean /= counter;
        if(freq_bins.find(mean) != freq_bins.end()){
            freq_bins[mean] = freq_bins[mean] + counter;
        }
        else freq_bins.insert(make_pair(mean, counter));
    }

    if(freq_bins.size() < 2){
        cerr << "There are too few bins.\n";
        return dos_vs_freq;
    }

    it = freq_bins.begin();
    if(it != freq_bins.end()){
        next_iter = next(it, 1);
        dos = next_iter->second / (next_iter->first - it->first) / num_sets;
        dos_vs_freq.insert(make_pair(it->first, dos));
    }

    for(it = next(it,1); it != prev(freq_bins.end(), 1) ; it++){
        prev_iter = prev(it,1);
        next_iter = next(it, 1);
        diff1 = it->first - prev_iter->first;
        diff2 = next_iter->first - it->first;

        dos = next_iter->second * diff1 * diff1 + it->second * diff2 * diff2;
        dos /= diff1 * diff2*(diff1 + diff2) * num_sets;
        dos_vs_freq.insert(make_pair(it->first, dos));
    }

    prev_iter = prev(it, 1);
    dos = it->second / (it->first - prev_iter->first) / num_sets;
    dos_vs_freq.insert(make_pair(it->first, dos));

    integral = 0;
    for(auto miter = dos_vs_freq.begin();miter != next(dos_vs_freq.end(),-1); miter++){
        next_fmap_iter = next(miter);
        integral += (next_fmap_iter->first - miter->first)*miter->second;
    }
    correction = (double) freqs.size() / integral;

    for(it = freq_bins.begin(); it != freq_bins.end(); it++){
        dos_vs_freq[it->first] = dos_vs_freq[it->first];// * correction;
    }

    return dos_vs_freq;
}

void run_dos_calc(vector<string> dos_files, string out_name){
    double next_freq;
    int num_sets = 0;
    vector<double> freqs;
    ifstream next_file;
    string nextline;
    map<double, double> dos_vs_freq;
    ofstream outfile;

    //Read in sets of frequncy, DOS coordinate pairs
    for(string filename : dos_files){
        next_file.open(filename);
        if(next_file.is_open()){

            while(! next_file.eof()){
                getline(next_file, nextline);
                if(sscanf(nextline.c_str(), "%lf", &next_freq) == 1){
                    freqs.push_back(next_freq);
                }
             }

             num_sets++;
             next_file.close();
         }
    }

    dos_vs_freq = get_dos(freqs, num_sets);

    if(out_name.compare("") == 0){
        for(auto iter = dos_vs_freq.begin(); iter != dos_vs_freq.end(); iter++){
            cout << iter->first << "\t" << iter->second << "\n";
        }
    }

    else{
        outfile.open(out_name);
        for(auto iter = dos_vs_freq.begin(); iter != dos_vs_freq.end(); iter++){
            outfile << iter->first << "\t" << iter->second << "\n";
        }

        outfile.close();
    }
}

int main(int argc, char **argv){

    DIR *dir = NULL;
    dirent *entry;
    string dir_name, extension, pattern, next_name, base_name, report_name = "";
    vector<string> matches;
    smatch results;
    vector<char> meta = {'[','\\','^','$','.','|','?','*','+','(',')','{','}'};

    cout << "Enter a directory name: ";
    getline(cin, dir_name);
    dir = opendir(dir_name.c_str());

    if(dir == NULL){
        cerr << "The directory " << dir_name << " was not found.\n";
        return -1;
    }

    cout << "Enter a base name: ";
    getline(cin, base_name);

    cout << "Enter a file extension: ";
    getline(cin, extension);
    if(extension == "") extension = DEFAULT_EXT;

    regex r("^" + base_name + ".*\\." + extension);

    while((entry = readdir(dir))) {
        next_name = entry->d_name;
        if(regex_search(next_name, results, r)){
            matches.push_back(results.str());
        }
    }

    closedir(dir);

    for(string next_match : matches){
        cout << next_match << "\n";
    }

    do{
        cout << "Enter a file name for reporting calculations: ";
        getline(cin, report_name);
        if(report_name.compare("") == 0){
            if(yesno("Proceed without a file name for reporting reults?")){
                break;
            }
        }
    }while(report_name.compare("") == 0);

    run_dos_calc(matches, report_name);

    return 0;
}
