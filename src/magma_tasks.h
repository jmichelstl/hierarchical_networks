#include <magma_v2.h>

/*
Classes to provide a means of loading calculations into a magma_thread_queue so that reporting of calculation results will be carried out in a thread-safe 
manner
*/

//Carry out a calculation of just eigenvalues on a cpu
class evals_cpu: public magma_task {

public:
    evals_cpu(double *mat, double *vals, int my_size) : stiff_mat(mat), evals(vals), size(my_size) {}

    virtual void run(){
        magma_int_t lwork, liwork, info = 0, *iwork, nb, max_dim;
        double *work;
        string eval_filename, evec_filename;

        nb = magma_get_dsytrd_nb(size);
        max_dim = nb > size ? nb : size;
        lwork = 3*max_dim + 2*max_dim*max_dim;
        liwork = size;
        iwork = (magma_int_t *) malloc(sizeof(magma_int_t) * (int) liwork);
        work = (double *) malloc(sizeof(double) * (int) lwork);

        magma_dsyevd(MagmaNoVec, MagmaUpper, size, stiff_mat, size, evals, work, lwork, iwork, liwork, &info);

        free(work);
        free(iwork);
    }

private:
    double *stiff_mat, *evals;
    int size;
};

//Carry out a calculation of eigenvalues and eigenvectors on a cpu
class evals_evecs_cpu: public magma_task {

public:
    evals_evecs_cpu(double *mat, double *vals, int my_size) : stiff_mat(mat), evals(vals), size(my_size){}

    virtual void run(){
        magma_int_t lwork, liwork, info = 0, *iwork, nb, max_dim;
        double *work;
        magma_queue_t queue = NULL;

        nb = magma_get_dsytrd_nb(size);
        max_dim = nb > size ? nb : size;
        lwork = 7*max_dim + 3*max_dim*max_dim;
        liwork = 3 + 6*max_dim;
        iwork = (magma_int_t *) malloc(sizeof(magma_int_t) * (int) liwork);
        work = (double *) malloc(sizeof(double) * (int) lwork);

        magma_dsyevd(MagmaVec, MagmaUpper, size, stiff_mat, size, evals, work, lwork,iwork, liwork, &info);

        free(work);
        free(iwork);
    }

private:
    double *stiff_mat, *evals;
    int size;
};

//Carry out a calculation of just eigenvalues on a gpu
class evals_gpu: public magma_task {

public:
    evals_gpu(double *mat, double *vals, int my_size) : stiff_mat(mat), evals(vals), size(my_size){}

    virtual void run(){
        magma_int_t lwork, liwork, ldwa, info = 0, dev = 0, *iwork, nb, max_dim;
        magmaDouble_ptr gpuMat;
        double *wA, *work;
        magma_queue_t queue;

        magma_queue_create(dev, &queue);
        magma_dmalloc(&gpuMat, size*size);
        magma_setmatrix(size, size, sizeof(double), stiff_mat, size, gpuMat, size, queue);
        free(stiff_mat);

        nb = magma_get_dsytrd_nb(size);
        max_dim = nb > size ? nb : size;
        lwork = 3*max_dim + 2*max_dim*max_dim;
        liwork = size;
        iwork = (magma_int_t *) malloc(sizeof(magma_int_t) * (int) liwork);
        work = (double *) malloc(sizeof(double) * (int) lwork);
        wA = (double *) malloc(sizeof(double) *  2 * size * size);

        magma_dsyevd_gpu(MagmaNoVec, MagmaUpper, size, gpuMat, size, evals, wA, 2*size, work, lwork, iwork, liwork, &info);

        free(work);
        free(iwork);
        free(wA);
        magma_queue_destroy(queue);
        magma_free(gpuMat);
    }

private:
    double *stiff_mat, *evals;
    int size;
};

//Carry out a calculation of eigenvalues and eigenvectors on a gpu
class evals_evecs_gpu: public magma_task {

public:
    evals_evecs_gpu(double **mat, double *vals, int my_size) : stiff_mat(mat), evals(vals), size(my_size) {}

    virtual void run(){
        magma_int_t lwork, liwork, ldwa, info = 0, dev = 0, *iwork, nb, max_dim;
        magmaDouble_ptr gpuMat;
        double *wA, *work;
        magma_queue_t queue;

        magma_queue_create(dev, &queue);
        magma_dmalloc(&gpuMat, size*size);
        magma_setmatrix(size, size, sizeof(double), *stiff_mat, size, gpuMat, size, queue);
        free(*stiff_mat);

        nb = magma_get_dsytrd_nb(size);
        max_dim = nb > size ? nb : size;
        lwork = 7*max_dim + 3*max_dim*max_dim;
        liwork = 6 * size;
        iwork = (magma_int_t *) malloc(sizeof(magma_int_t) * (int) liwork);
        work = (double *) malloc(sizeof(double) * (int) lwork);
        wA = (double *) malloc(sizeof(double) *  2 * size * size);

        magma_dsyevd_gpu(MagmaVec, MagmaUpper, size, gpuMat, size, evals, wA, 2*size, work, lwork, iwork, liwork, &info);

        free(work);
        free(iwork);
        free(wA);

        *stiff_mat = (double *) malloc(sizeof(double) * size * size);
        magma_dgetmatrix(size, size, gpuMat, size, *stiff_mat, size, queue);

        magma_queue_destroy(queue);
        magma_free(gpuMat);
    }

private:
    double **stiff_mat, *evals;
    int size;
};
