#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include "network_utils.h"
//#include <gsl/gsl_math.h>
//#include <gsl/gsl_eigen.h>
//#include <gsl/gsl_sort.h>
//#include <gsl/gsl_sort_vector.h>
#include <cmath>
#include <magma_v2.h>

using namespace std;

#define X 0
#define Y 1

enum calc_mode {lapack, magma_cpu, magma_gpu};

extern "C" void dsyevd_(char *jobz, char *uplo, int *n, double *
        a, int *lda, double *w, double *work, int *lwork,
        int *iwork, int *liwork, int *info); 

struct EdgeDatum{

    EdgeDatum(int ival, double sval, double ilc) : index(ival), stiffness(sval), inv_l_cubed(ilc){ 
    }

    EdgeDatum(double ival, double lval) : index(ival), inv_l_cubed(lval){
        stiffness = 1;
    }

    int index;
    double stiffness;
    double inv_l_cubed;
};


void read_edges(ifstream& datfile, vector<Point>& point_list, map<int, vector<EdgeDatum>>& neighbor_map){

    map<Point, int> pmap;
    string nextline;
    Point p1, p2;
    int pindex = 0, point_count = 0, index1, index2, mindex, maxdex;
    vector<double> edge_data;
    double length, inv_l_cube;

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
                    point_list.push_back(p1);
                    neighbor_map.insert(make_pair(pindex, vector<EdgeDatum>()));                    pindex++;
                }
                if(pmap.find(p2) == pmap.end()){
                    pmap.insert(make_pair(p2, pindex));
                    point_list.push_back(p2);
                    neighbor_map.insert(make_pair(pindex, vector<EdgeDatum>()));
                    pindex++;
                }

                index1 = pmap[p1];
                index2 = pmap[p2];
                length = sqrt((p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y));
                inv_l_cube = 1 / length / length / length;
                mindex = index1 < index2 ? index1 : index2;
                maxdex = index1 == mindex ? index2 : index1;

                neighbor_map[mindex].push_back(EdgeDatum(maxdex, inv_l_cube));
            }
        }
    }
}

int stiff_mat_index(int index1, int index2, int type1, int type2, int num_points, calc_mode mode){
    return 2 * (2 * index1 + type1) * num_points + 2 * index2 + type2;
}

void load_stiff_mat(double **stiff_mat, vector<Point> point_list, map<int, vector<EdgeDatum>> neighbor_map, calc_mode mode){
    int point_count = point_list.size();
    int index1, index2, mat_index, xx_index1, yx_index1, yy_index1;
    int xx_index2, yx_index2, yy_index2;
    Point p1, p2;
    double inv_l_cubed, stiffness, mat_elem;

    for(int i = 0; i < 4 * point_count * point_count; i++){
        stiff_mat[0][i] = 0;
    }

    for(auto it = neighbor_map.begin(); it != neighbor_map.end(); it++){
        index1 = it->first;
        p1 = point_list[index1];
        xx_index1 = stiff_mat_index(index1, index1, X, X, point_count, mode);
        yx_index1 = stiff_mat_index(index1, index1, Y, X, point_count, mode);
        yy_index1 = stiff_mat_index(index1, index1, Y, Y, point_count, mode);

        for(EdgeDatum edat : it->second){
            index2 = edat.index;
            stiffness = edat.stiffness;
            p2 = point_list[index2];
            inv_l_cubed = edat.inv_l_cubed;
        
            xx_index2 = stiff_mat_index(index2, index2, X, X, point_count,mode);
            yx_index2 = stiff_mat_index(index2, index2, Y, X, point_count,mode);
            yy_index2 = stiff_mat_index(index2, index2, Y, Y, point_count,mode);

            mat_index = stiff_mat_index(index2, index1, X, X, point_count,mode);
            mat_elem = stiffness * inv_l_cubed * (p1.x - p2.x)*(p2.x - p1.x);
            stiff_mat[0][mat_index] = mat_elem;
            stiff_mat[0][xx_index1] -= mat_elem;
            stiff_mat[0][xx_index2] -= mat_elem;

            mat_index = stiff_mat_index(index2, index1, X, Y, point_count,mode);
            mat_elem = stiffness * inv_l_cubed * (p1.x - p2.x)*(p2.y - p1.y);
            stiff_mat[0][mat_index] = mat_elem;
            mat_index = stiff_mat_index(index2, index1, Y, X, point_count,mode);
            stiff_mat[0][mat_index] = mat_elem;
            stiff_mat[0][yx_index1] -= mat_elem;
            stiff_mat[0][yx_index2] -= mat_elem;

            mat_index = stiff_mat_index(index2, index1, Y, Y, point_count,mode);
            mat_elem = stiffness * inv_l_cubed * (p1.y - p2.y)*(p2.y - p1.y);
            stiff_mat[0][mat_index] = mat_elem;
            stiff_mat[0][yy_index1] -= mat_elem;
            stiff_mat[0][yy_index2] -= mat_elem;
        }
    }
}

void calc_evals_magma_cpu(double *stiff_mat, double *evals, int size){
    magma_int_t lwork, liwork, info = 0, *iwork;
    double *work;
    string eval_filename, evec_filename;
    double aux_work[1];
    magma_int_t aux_iwork[1];

    //Obtain optimal workspace sizes
    magma_dsyevd(MagmaNoVec, MagmaUpper, size, stiff_mat, size, evals, aux_work, -1, aux_iwork, -1, &info);

    lwork = aux_work[0];
    liwork = aux_iwork[0];
    iwork = (magma_int_t *) malloc(sizeof(magma_int_t) * (int) liwork);
    work = (double *) malloc(sizeof(double) * (int) lwork);

    magma_dsyevd(MagmaNoVec, MagmaUpper, size, stiff_mat, size, evals, work, lwork, iwork, liwork, &info);

    if(info != 0){
        const char *message = magma_strerror(info);
        fprintf(stderr, "%s\n", message);
    }

    free(work);
    free(iwork);
}

void calc_evals_evecs_magma_cpu(double *stiff_mat, double *evals, int size){
    magma_int_t lwork, liwork, info = 0, *iwork;
    double *work;
    double aux_work[1];
    magma_int_t aux_iwork[1];

    //Obtain optimal workspace sizes
    magma_dsyevd(MagmaVec, MagmaUpper, size, stiff_mat, size, evals, aux_work, -1, aux_iwork, -1, &info);

    lwork = aux_work[0];
    liwork = aux_iwork[0];
    iwork = (magma_int_t *) malloc(sizeof(magma_int_t) * (int) liwork);
    work = (double *) malloc(sizeof(double) * (int) lwork);

    magma_dsyevd(MagmaVec, MagmaUpper, size, stiff_mat, size, evals, work, lwork, iwork, liwork, &info);

    if(info != 0){
        const char *message = magma_strerror(info);
        fprintf(stderr, "%s\n", message);
    }

    free(work);
    free(iwork);
}

void calc_evals_magma_gpu(double *stiff_mat, double *evals, int size){
    magma_int_t lwork, liwork, ldwa, info = 0, dev = 0, *iwork;
    magmaDouble_ptr gpuMat;
    double *wA, *work;
    magma_queue_t queue;
    char *message;
    double aux_work[1];
    magma_int_t aux_iwork[1];

    magma_queue_create(dev, &queue);
    magma_dmalloc(&gpuMat, size*size);
    magma_dsetmatrix(size, size, stiff_mat, size, gpuMat, size, queue);
    //free(stiff_mat);

    wA = (double *) malloc(sizeof(double) * size * size);

    //Obtain optimal workspace sizes
    magma_dsyevd_gpu(MagmaNoVec, MagmaUpper, size, gpuMat, size, evals, wA, size, aux_work, -1, aux_iwork, -1, &info);

    lwork = aux_work[0];
    liwork = aux_iwork[0];
    iwork = (magma_int_t *) malloc(sizeof(magma_int_t) * (int) liwork);
    work = (double *) malloc(sizeof(double) * (int) lwork);

    magma_dsyevd_gpu(MagmaNoVec, MagmaUpper, size, gpuMat, size, evals, wA, size, work, lwork, iwork, liwork, &info);

    if(info != 0){
        const char *message = magma_strerror(info);
        fprintf(stderr, "%s\n", message);
    }

    free(work);
    free(iwork);
    free(wA);
    magma_queue_destroy(queue);
    magma_free(gpuMat);
}

void calc_evals_evecs_magma_gpu(double **stiff_mat, double *evals, int size){
    magma_int_t lwork, liwork, ldwa, info = 0, dev = 0, *iwork;
    magmaDouble_ptr gpuMat;
    double *wA, *work;
    magma_queue_t queue;
    double aux_work[1];
    magma_int_t aux_iwork[1];

    magma_queue_create(dev, &queue);
    magma_dmalloc(&gpuMat, size*size);
    magma_dsetmatrix(size, size, *stiff_mat, size, gpuMat, size, queue);
    free(*stiff_mat);

    wA = (double *) malloc(sizeof(double) * size * size);

    //Obtain optimal workspace sizes
    magma_dsyevd_gpu(MagmaVec, MagmaUpper, size, gpuMat, size, evals, wA, size, aux_work, -1, aux_iwork, -1, &info);

    lwork = aux_work[0];
    liwork = aux_iwork[0];
    iwork = (magma_int_t *) malloc(sizeof(magma_int_t) * (int) liwork);
    work = (double *) malloc(sizeof(double) * (int) lwork);

    magma_dsyevd_gpu(MagmaVec, MagmaUpper, size, gpuMat, size, evals, wA, size, work, lwork, iwork, liwork, &info);

    if(info != 0){
        const char *message = magma_strerror(info);
        fprintf(stderr, "%s\n", message);
    }

    *stiff_mat = (double *) malloc(sizeof(double) * size * size);
    magma_dgetmatrix(size, size, gpuMat, size, *stiff_mat, size, queue);

    free(work);
    free(iwork);
    free(wA);
    magma_queue_destroy(queue);
    magma_free(gpuMat);
}
/*
void calc_evals(double *stiff_mat, gsl_vector *eval, int size){

    gsl_eigen_symm_workspace *wkspace = gsl_eigen_symm_alloc(size);
    gsl_matrix_view mview = gsl_matrix_view_array(stiff_mat, size, size);
    gsl_eigen_symm(&mview.matrix, eval, wkspace);
    gsl_sort_vector(eval);
    gsl_eigen_symm_free(wkspace);
}
*/

/*
void calc_evals_evecs(double * stiff_mat, gsl_vector *eval, gsl_matrix *evec, int dim){

    gsl_eigen_symmv_workspace * wkspace = gsl_eigen_symmv_alloc(dim);
    gsl_matrix_view mview = gsl_matrix_view_array(stiff_mat, dim, dim);
    gsl_eigen_symmv(&mview.matrix, eval, evec, wkspace);
    gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);
    gsl_eigen_symmv_free(wkspace);
}
*/

void calc_evals_lapack(double *stiff_mat, double *evals, int size){
    int lwork, liwork, info = 0, flag = -1, *iwork;
    double *work;
    string eval_filename, evec_filename;
    double aux_work[1];
    int aux_iwork[1];
    char job = 'N', upper = 'U';

    //Obtain optimal workspace sizes
    dsyevd_(&job, &upper, &size, stiff_mat, &size, evals, aux_work, &flag, aux_iwork, &flag, &info);

    lwork = aux_work[0];
    liwork = aux_iwork[0];
    iwork = (int *) malloc(sizeof(int) * liwork);
    work = (double *) malloc(sizeof(double) * lwork);

    dsyevd_(&job, &upper, &size, stiff_mat, &size, evals, work, &lwork, iwork, &liwork, &info);

    free(work);
    free(iwork);
}

void calc_evals_evecs_lapack(double *stiff_mat, double *evals, int size){
    int lwork, liwork, info = 0, flag = -1, *iwork;
    double *work;
    double aux_work[1];
    int aux_iwork[1];
    char job = 'V', upper = 'U';

    //Obtain optimal workspace sizes
    dsyevd_(&job, &upper, &size, stiff_mat, &size, evals, aux_work, &flag, aux_iwork, &flag, &info);

    lwork = aux_work[0];
    liwork = aux_iwork[0];
    iwork = (int *) malloc(sizeof(int) * liwork);
    work = (double *) malloc(sizeof(double) * lwork);

    dsyevd_(&job, &upper, &size, stiff_mat, &size, evals, work, &lwork, iwork, &liwork, &info);

    free(work);
    free(iwork);
}

void report_evals(double *evals, string filename, int dim){
    FILE *outfile;
    int i;

    outfile = fopen(filename.c_str(), "w");
    if(outfile != NULL){
        for(i = 0; i < dim; i++){
            fprintf(outfile, "%1.16le\n", sqrt(abs(evals[i])));
        }
        fclose(outfile);
    }

    else cerr << "Failed to open " << filename << " for writing.\n";
}

void report_evecs(double *evecs, string filename, int dim){
    FILE *outfile;
    int i, j;

    outfile = fopen(filename.c_str(), "w");
    if(outfile != NULL){
        for(int i = 0; i < dim; i++){
            for(int j = 0; j < dim; j++){
                fprintf(outfile, "%1.16le\t", evecs[dim*i + j]);
            }
            fprintf(outfile, "\n");
        }
        fclose(outfile);
    }

    else cerr << "Failed to open " << filename << " for writing.\n";
}

void report_evals_binary(double *evals, string filename, int dim){
    FILE *outfile;

    outfile = fopen(filename.c_str(), "wb");
    fwrite(&dim, sizeof(int), 1, outfile);
    fwrite(evals, sizeof(double), dim, outfile);
    fclose(outfile);
}

void report_evecs_binary(double *evecs, string filename, int dim){
    FILE *outfile;

    outfile = fopen(filename.c_str(), "wb");
    fwrite(&dim, sizeof(int) , 1, outfile);
    fwrite(evecs, sizeof(double), dim*dim, outfile);
    fclose(outfile);
}

/*
void report_evals(gsl_vector *evals, string filename, int dim){
    ofstream outfile;
    int i;

    outfile.open(filename, ofstream::out | ofstream::app);
    if(outfile.is_open()){
        for(i = 0; i < dim; i++){
            outfile << sqrt(abs(gsl_vector_get(evals, i))) << "\n";
        }
        outfile.close();
    }

    else cerr << "Failed to open  " << filename << " for writing.\n";
}
*/

/*
void report_evecs(gsl_matrix *evec, string filename, int dim){
    FILE * outfile;
    gsl_vector_view next_evec;
    gsl_vector vec;
    int i, j;

    outfile = fopen(filename.c_str(), "w");
    if(outfile != NULL){
        for(i = 0; i < dim; i++){
            next_evec = gsl_matrix_column(evec, i);
            vec = next_evec.vector;
            //gsl_vector_fprintf(outfile, &next_evec.vector, "%g");
            for(j = 0; j < dim; j++){
                fprintf(outfile, "%lf\t", gsl_vector_get(&vec, j)); 
            }
            fprintf(outfile, "\n");
        }
        fclose(outfile);
    }

    else cerr << "Failed to open " << filename << " for writing.\n";
}
*/

void process_dat_file(ifstream& datfile, calc_mode mode){

    vector<Point> point_list;
    map<int, vector<EdgeDatum>> neighbor_map;
    double *stiff_mat, *evals;
    int point_count, dim;
    /*gsl_matrix_view mview;
    gsl_vector *evals;
    gsl_eigen_symm_workspace * wkspace;
    gsl_matrix *evec;*/
    ofstream outfile;
    string eval_filename, evec_filename;

    //Read information about points and edges from data file specifying network
    read_edges(datfile, point_list, neighbor_map);
    point_count = point_list.size();
    datfile.close();

    if(point_count == 0){
        cout << "No edges were read.\n";
        return;
    }

    cout << "There are " << point_count << " points.\n";
    dim = 2*point_count;
    //Create the stiffness matrix to be diagonalized
    stiff_mat = (double *) malloc(4*point_count*point_count * sizeof(double));
    load_stiff_mat(&stiff_mat, point_list, neighbor_map, mode);

    /*
    for(int i = 0; i < 2*point_count; i++){
        for(int j = 0; j < 2*point_count; j++){
            cout << stiff_mat[i*2*point_count + j] << " ";
        }
        cout << "\n";
    }
    */

    //Calculate eignvalues of stiffness matrix to determine vibrational DOS
    //if(mode == calc_mode::gsl) evals = gsl_vector_alloc(2 * point_count);

    evals = (double *) malloc(sizeof(double) * dim);
    
    if(yesno("Calculate eigenvectors?")){
        /*if(mode == calc_mode::gsl){
            evec = gsl_matrix_alloc(dim, dim);
            calc_evals_evecs(stiff_mat, evals, evec, dim);
        }*/

        if(mode == calc_mode::lapack){
            calc_evals_evecs_lapack(stiff_mat, evals, dim);
        }

        else if(mode == calc_mode::magma_cpu){
            calc_evals_evecs_magma_cpu(stiff_mat, evals, dim);
        }
        else{
            calc_evals_evecs_magma_gpu(&stiff_mat, evals, dim);
        }
        evec_filename = enter_decline("Enter a name for eigen vector file");
        if(evec_filename.compare("") != 0){
            if(yesno("Use binary form?")){
                report_evecs_binary(stiff_mat, evec_filename, dim);
             }
            //if(mode == calc_mode::gsl) report_evecs(evec, evec_filename, dim);
            else report_evecs(stiff_mat, evec_filename, dim);
        }
        //if(mode == calc_mode::gsl) gsl_matrix_free(evec);
    }
    else{
        //if(mode == calc_mode::gsl) calc_evals(stiff_mat, evals, dim);
        if(mode == calc_mode::lapack) calc_evals_lapack(stiff_mat, evals, dim);
        else if(mode == calc_mode::magma_cpu){
            calc_evals_magma_cpu(stiff_mat, evals, dim);
        }
        else{
            calc_evals_magma_gpu(stiff_mat, evals, dim);
        }
    }

    eval_filename = enter_decline("Enter a name for eigen value file");
    if(eval_filename.compare("") != 0){
        if(yesno("Use binary form?")){
            report_evals_binary(evals, eval_filename, dim);
        }
        //if(mode == calc_mode::gsl) report_evals(evals, eval_filename, dim);
        else report_evals(evals, eval_filename, dim);
    }

    //if(mode == calc_mode::gsl) gsl_vector_free(evals);
    free(evals);
    free(stiff_mat);
    //datfile.close();
}

int main(int argc, char **argv){

    string filename, nextline;
    ifstream netfile;
    calc_mode mode = calc_mode::lapack;
    char c;

    do{

        if(yesno("Use magma?")){
            if(yesno("Use gpu?")) mode = calc_mode::magma_gpu;
            else mode = calc_mode::magma_cpu;
            magma_init();
        }

        do{
            cout << "Enter the next file name: ";
            getline(cin, filename);
            netfile.open(filename);
            if(! netfile.is_open()){
                if(! yesno("File could not be read. Try again?")) break;
            }
        }while(! netfile.is_open());

        if(netfile.is_open()) process_dat_file(netfile, mode);

        if(mode == calc_mode::magma_cpu || mode == calc_mode::magma_gpu){
            magma_finalize();
        }

    }while(yesno("Process another network?"));

    return 0;
}
