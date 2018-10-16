/************************************************
 * 
 * Created: 10th Feb, 2017
 * Updated: 10th June, 2017
 *
 ************************************************/

#include <fstream>
#include <iostream>
#include <cstdio>
#include <fstream> 
#include "apply_Htarget.h"
#include <iostream>
#include <vector>
#include "Matrix.h"
#include <algorithm>
#include <string>
#include <tuple>
#include <cmath>
#include <dirent.h>
#include <sys/stat.h>
#include <errno.h>
#include <omp.h>

bool isDir(std::string dir)
{
    struct stat fileInfo;
    stat(dir.c_str(), &fileInfo);
    if (S_ISDIR(fileInfo.st_mode))
    {
        return true;
    }
    else
    {
        return false;
    }
}

void listFiles(std::string baseDir, bool recursive, std::vector<std::string> &all_C_files)
{
    DIR *dp;
    struct dirent *dirp;
    if ((dp = opendir(baseDir.c_str())) == NULL)
    {
        std::cout << "[ERROR: " << errno << " ] Couldn't open " << baseDir << "." << std::endl;
        return;
    }
    else
    {
        while ((dirp = readdir(dp)) != NULL)
        {
            if (dirp->d_name != std::string(".") && dirp->d_name != std::string(".."))
            {
                if (isDir(baseDir + dirp->d_name) == true && recursive == true)
                {
                    listFiles(baseDir + dirp->d_name + "/", true, all_C_files);
                }
                else
                {
                    all_C_files.push_back(baseDir + dirp->d_name);
                }
            }
        }
        closedir(dp);
    }
}

int main(int argc, char **argv)
{
    //Kokkos::initialize (argc, argv);
    char buf[256];
    int nrow, ncol, nnz;
    int ii, jj;
    double val;
    auto C = new Block_Matrix_t;

    if (argc != 7)
    {
        printf("Format: traverse C_DIR_NAME X_VEC Y_VEC L_PATCH_SIZE R_PATCH_SIZE DENSE_THRESHOLD\n");
        exit(1);
    }
    std::string c_dir_name(argv[1]);
    std::string x_file_name(argv[2]);
    std::string y_file_name(argv[3]);
    std::string l_patch_filename(argv[4]);
    std::string r_patch_filename(argv[5]);
    float dense_threshold = std::stof(std::string(argv[6]));

    ii = 0;
    jj = 0;
    int n_c_rows = 0;
    int n_c_cols = 0;
    std::string format = c_dir_name + "/C.%d.%d";
    DIR *dp;
    struct dirent *dirp;
    if ((dp = opendir(c_dir_name.c_str())) == NULL)
    {
        std::cout << "[ERROR: " << errno << " ] Couldn't open " << c_dir_name << "." << std::endl;
        exit(1);
    }
    while ((dirp = readdir(dp)) != NULL)
    {
        if (dirp->d_name == std::string(".") || dirp->d_name == std::string(".."))
            continue;
        std::string p = c_dir_name + "/" + dirp->d_name;
        if (isDir(p) == true)
        {
            sscanf(p.c_str(), format.c_str(), &ii, &jj);
            n_c_rows = std::max(n_c_rows, ii);
            n_c_cols = std::max(n_c_cols, jj);
        }
    }
    closedir(dp);

    C->cij.resize(n_c_rows);
    for (auto &a : C->cij)
    {
        a.resize(n_c_cols, nullptr);
    }
    std::string file_format = format + "/%c.%d";

    std::vector<std::string> all_C_files;
    if (c_dir_name.back() != '/')
    {
        c_dir_name += "/";
    }

    listFiles(c_dir_name, true, all_C_files);
    for (auto &p : all_C_files)
    {
        int ic, jc, imatrix;
        char AorB;
        sscanf(p.c_str(), file_format.c_str(), &ic, &jc, &AorB, &imatrix);
        ic -= 1;
        jc -= 1; 
        std::ifstream ifs(p, std::ifstream::in);
        ifs.getline(buf, 256); //Comment line
        ifs.getline(buf, 256);
        sscanf(buf, "%d %d %d", &nrow, &ncol, &nnz);
        auto newmat = new Matrix_t;
        newmat->nrow = nrow;
        newmat->ncol = ncol;
        if (((float)(nnz) / ((float)(nrow)*ncol)) <= dense_threshold)
        {
            //printf("C[%d][%d]->%c[%d] IS NOT DENSE\n", ic, jc, AorB, imatrix);
            newmat->is_dense = false;
            newmat->rowptr.resize(nrow + 1);
            newmat->col.resize(nnz);
            newmat->val.resize(nnz);
            newmat->rowptr[0] = 0;
        }
        else
        { 
            //Matrix treated as dense  - Store colmun major for BLAS reasons
            //printf("C[%d][%d]->%c[%d] IS DENSE\n", ic, jc, AorB, imatrix);
            newmat->is_dense = true;
            newmat->val.resize(nrow * ncol);
        }

        std::vector<std::tuple<int, int, double>> tmp_vec(nnz);
        for (int i = 0; i < nnz; i++)
        {
            ifs.getline(buf, 256);
            sscanf(buf, "%d  %d  %lf", &ii, &jj, &val);
            ii -= 1;
            jj -= 1; //Adjust for 0-based indexing
            auto t = std::make_tuple(ii, jj, val);
            tmp_vec[i] = t;
        }
        if (newmat->is_dense)
        {
            for (int i = 0; i < nnz; i++)
            {
                int ii = std::get<0>(tmp_vec[i]);
                int jj = std::get<1>(tmp_vec[i]);
                double val = std::get<2>(tmp_vec[i]);
                newmat->val[(ii) + (jj)*newmat->nrow] = val;
            }
        }
        else
        {
            std::sort(tmp_vec.begin(), tmp_vec.end());
            for (int i = 0; i < nnz; i++)
            {
                int ii = std::get<0>(tmp_vec[i]);
                int jj = std::get<1>(tmp_vec[i]);
                double val = std::get<2>(tmp_vec[i]);
                newmat->rowptr[ii + 1] += 1;
                newmat->col[i] = jj;
                newmat->val[i] = val;
            }
            for (int i = 1; i <= nrow; i++)
            {
                newmat->rowptr[i] += newmat->rowptr[i - 1];
            }
        }

        ifs.close();
        if (C->cij[ic][jc] == nullptr)
        {
            C->cij[ic][jc] = new CIJ_Elem_t;
        }
        CIJ_Elem celem = C->cij[ic][jc];
        if (AorB == 'A')
        {
            if (imatrix > celem->A.size())
            {
                celem->A.resize(imatrix);
            }
            celem->A[imatrix - 1] = newmat;
        }
        else
        {
            if (imatrix > celem->B.size())
            {
                celem->B.resize(imatrix);
            }
            celem->B[imatrix - 1] = newmat;
        }
    }

    int vec_size = 0;
    double element = 0.0;

    std::ifstream xfs(x_file_name, std::ifstream::in);
    xfs >> vec_size;

    std::vector<double> X;
    X.reserve(vec_size);
    while (xfs >> element)
    {
        X.push_back(element);
    }
    xfs.close();

    std::ifstream yfs(y_file_name, std::ifstream::in);
    yfs >> vec_size;

    std::vector<double> Y_target;
    Y_target.reserve(vec_size);
    while (yfs >> element)
    {
        Y_target.push_back(element);
    }
    yfs.close();


    printf("#######################################################\n");
    fflush(stdout);

    int isize;
    std::ifstream lfs(l_patch_filename, std::ifstream::in);
    lfs >> vec_size;
    std::vector<int> l_patch_size;
    l_patch_size.reserve(vec_size);
    while (lfs >> isize)
    {
        l_patch_size.push_back(isize);
    }
    lfs.close();
    

    std::ifstream rfs(r_patch_filename, std::ifstream::in);
    rfs >> vec_size;

    std::vector<int> r_patch_size;
    r_patch_size.reserve(vec_size);
    while (rfs >> isize)
    {
        r_patch_size.push_back(isize);
    } 
    rfs.close();

    std::vector<double> Y(X.size());

    double start,end;
    start = omp_get_wtime();
    apply_Htarget(*C, l_patch_size, r_patch_size, X, Y);
    end = omp_get_wtime();
    std::cout<<" Execution Time: "<< (double) ( (end-start)) << " : ";
    double diff_norm = 0.0;
    for (int i = 0; i < Y.size(); i++)
    {
        //printf("-- %d  %.8e   %.8e    %.8e\n", i, Y[i], Y_target[i], std::abs(Y[i] - Y_target[i]));
        diff_norm += std::abs(Y[i] - Y_target[i]);
    }
    printf("Correctness = %.6e\n", diff_norm);
    //Kokkos::finalize ();
}
