#include <stdio.h>
#include "network_utils.h"
#include <iostream>
#include <stdlib.h>
#include <unistd.h>

using namespace std;

//Prompt the user for a binary file containing eigen vectors, and calculate the
//participation ratio of each eigen vector. Files should contain an integer
//giving the dimension of each vector, followed by the vectors.
void calc_part_ratios(bool normalized){
    //Information about file with eigen vector components
    string vec_file_name;
    FILE* vec_file = NULL;
    ofstream outfile;

    //Information about the dimension of eigen vectors and the data file size
    int dim, size, target_size, num_read, i, j;
    double dividend, inv_dim, sum, prat, ej, ej_sq;
    double *evec;

    //Prompt the user for a data file
    while(vec_file == NULL){
        cout << "Enter the eigen vector file: ";
        getline(cin, vec_file_name);
        if(vec_file_name.compare("") != 0){
            vec_file = fopen(vec_file_name.c_str(), "rb");
        }

        if(vec_file == NULL){
            if(! yesno("No file was read. Try again?")){
                return;
            }
        }
    }

    //Get the dimension of vectors
    num_read = fread(&dim, sizeof(int), 1, vec_file);
    if(num_read == 0){
        cerr << "The dimension of the vectors could not be read.\n";
        fclose(vec_file);
        return;
    }

    fseek(vec_file, 0, SEEK_END);
    size = ftell(vec_file) - sizeof(int);
    target_size = sizeof(double) * dim * dim;
    if(size < target_size){
        cerr << "There was too little data in the file.\n";
        fclose(vec_file);
        return;
    }

    open_output_file("Enter a file for reporting ratios: ", outfile);

    //Set the position indicator for the data file just past the lead integer
    fseek(vec_file, sizeof(int), SEEK_SET);

    inv_dim = 1 / (float) dim;

    //Read eigen vectors one by one and determine their participation ratios
    evec = (double *) malloc(dim * sizeof(double));

    for(i = 0; i < dim; i++){
        num_read = fread(evec, sizeof(double), dim, vec_file);

        if(num_read < dim){
            cerr << "A reading error occurred for vector " << i << ".\n";
            fclose(vec_file);
            if(outfile.is_open()) outfile.close();
            return;
        }

        sum = 0;
        if(! normalized){
            dividend = 0;
            for(j = 0; j < dim; j++){
                ej_sq = evec[j] * evec[j];
                dividend += ej_sq;
                sum += ej_sq * ej_sq;
            }
        }

        else{
            for(j = 0; j < dim; j++){
                ej = evec[j];
                sum += ej * ej * ej * ej;
            }
        }

        prat = inv_dim / sum;
        if(! normalized) prat *= dividend;

        if(outfile.is_open()) outfile << prat << "\n";
        else cout << prat << "\n";
    }

    fclose(vec_file);
    if(outfile.is_open()) outfile.close();

    return;
}

int main(int argc, char **argv){
    bool normalized = true;
    char c;

    while((c = getopt(argc, argv, "u")) != -1){
        switch(c) {
            case 'c':
                normalized = false;
                break;
            case '?':                if(isprint(optopt)){
                    fprintf(stderr, "Unknown option: -%c.\n", optopt);
                }
                else{
                    fprintf(stderr, "Unknown option character.\n");
                }
            default:
                break;
         }
    }

    //Until the user indicates otherwise, calculate the participation ratios
    //for the normal modes of a network, and report them to a user-specified
    //file
    do{
        calc_part_ratios(normalized);
    }while(yesno("Process another network?"));

    return 0;
}
