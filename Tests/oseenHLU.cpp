#include <medusa/Medusa_fwd.hpp>
#include <medusa/bits/domains/PolyhedronShape.hpp>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <Eigen/SparseCore>
#include <math.h>
#include <chrono>

#include "../h2libext.h"

using namespace mm;
using std::string;
using std::chrono::high_resolution_clock;
using std::chrono::duration;

void print_header_to_file(string path)
{
    std::ofstream file(path);

    file << "poly_grad;";
    file << "poly_div;";
    file << "poly_conv;";
    file << "poly_lap;";
    file << "nu;";
    file << "N_ui;";
    file << "N_p;";
    file << "N_total;";
    file << "iter;";
    file << "heps_vel;";
    file << "t_cluster_vel;";
    file << "t_cluster_p;";
    file << "t_lu_vel;";
    file << "t_grad_lower;";
    file << "t_grad_upper;";
    file << "t_schur_mul;";
    file << "t_schur_computation;";
    file << "t_lu_schur;";
    file << "t_prcd;";
    file << "t_solve;";
    file << "t_total;";
    file << "solution_error;";
    file << "solution_error_u;";
    file << "solution_error_p;";
    file << "solution_error_inf;";
    file << "solution_error_u_inf;";
    file << "solution_error_p_inf;";
    file << "residual_error;";
    file << "lu_vel_error;";
    file << "lu_schur_error;";
    file << "solver_tol;";

    file << std::endl;

    file.close();
}

// useful to work with matrix in Matlab
void writeMatrix2File(Eigen::SparseMatrix<double, Eigen::RowMajor> M) {
    int i = 0;
    Eigen::VectorXi rows(M.nonZeros());
    Eigen::VectorXi cols(M.nonZeros());
    Eigen::VectorXd values(M.nonZeros());
    for (int k=0; k<M.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(M,k); it; ++it) {
            rows(i) = it.row();
            cols(i) = it.col();
            values(i) = it.value();
            i++;
        }
    }

    std::ofstream rows_file("Daten/mat_rows.txt");
    for (int i = 0; i < rows.size(); i++) { rows_file << rows(i) << ";"; }
    rows_file.close();
    std::ofstream cols_file("Daten/mat_cols.txt");
    for (int i = 0; i < cols.size(); i++) { cols_file << cols(i) << ";"; }
    cols_file.close();
    std::ofstream values_file("Daten/mat_values.txt");
    values_file.precision(std::numeric_limits<double>::digits10 + 1);
    for (int i = 0; i < values.size(); i++) { values_file << values(i) << ";"; }
    values_file.close();
}

int write2vtk(Eigen::VectorXd vec, Range<Vec3d> positions) {
    // number of points
    int N = positions.size(); 

    // Create an ofstream object to manage file output
    std::ofstream outFile;

    // Open a file named "output.txt" in write mode
    outFile.open("output.vtk");

    // Check if the file was opened successfully
    if (!outFile) {
        std::cerr << "Error opening file" << std::endl;
        return 1;
    }

    // set precision: accurate to 9 decimal places
    outFile << std::setprecision(10);

    // Write formatted data to the file
    outFile << "# vtk DataFile Version 2.0\n";
    outFile << "Unstructured Grid Example\n";
    outFile << "ASCII\n\n";
    outFile << "DATASET UNSTRUCTURED_GRID\n";
    outFile << "POINTS " << N << " double\n";

    for (int i = 0; i < N; i++) {
        outFile << positions[i][0] << " " << positions[i][1] << " " << positions[i][2] << "\n";
    }

    outFile << "\nPOINT_DATA " << N << "\n";
    outFile << "VECTORS vectors double\n";
    for (int i = 0; i < N; i++) {
        outFile << vec[i] << " " << vec[i+N] << " " << vec[i+2*N] << "\n";
    }
    
    // Close the file
    outFile.close();

    // Inform the user that the data has been written
    std::cout << "Data has been written to output.vtk" << std::endl;

    return 0;
}

// create Matrix and RHS with Dirichlet boundary nodes eliminated
void createMatrixAndRHS(Eigen::SparseMatrix<double, Eigen::RowMajor>& mat, Eigen::VectorXd& rhs, OseenDiscretizationBetter& dc, int pressure_constraint = 0) {
    int N_all = 3*dc.N_u + dc.N_p + 1;
    int N_dofs = dc.N_dofs;
    // full system and righthand-side including dirichlet boundary conditions
    Eigen::SparseMatrix<double, Eigen::RowMajor> mat_temp(N_all,N_all);
    Eigen::VectorXd rhs_temp = setWeightsOseen_versatile_better(mat_temp, dc);

    // reduced system without dirichlet boundary conditions
    rhs = eliminateDirichlet(mat, mat_temp, dc, rhs_temp);

    if (pressure_constraint == 1)
    {
        // set last 2 rows and columns to zero and set 2x2 block to identity, removes pressure constraint
        // and sets pressure at node N_p to value of exact pressure solution at node N_p (N_p corresponds to N_dofs-2 overall)
        for (int k=0; k<mat.outerSize(); ++k)
            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(mat,k); it; ++it)
            {
                if (it.row()>N_dofs-3 || it.col()>N_dofs-2)
                    mat.coeffRef(it.row(),it.col()) = 0.0;
            }
        mat.coeffRef(N_dofs-2,N_dofs-2) = 1.0;
        mat.coeffRef(N_dofs-1,N_dofs-1) = 1.0;
        rhs(N_dofs-2) = dc.exact_solution(N_dofs-2); // set exact pressure solution at pressure node N_p (overall node N_dofs-2)
        mat.makeCompressed();
    }
}

psparsematrix symmetrizeVelocityBlockH2(Eigen::SparseMatrix<double, Eigen::RowMajor>& mat, int N_ui) {
    Eigen::SparseMatrix<double, Eigen::RowMajor> mat_vel = mat.block(0,0,N_ui,N_ui); // velocity block of main matrix
    // symmetrized A for nested dissection ordering
    Eigen::SparseMatrix<double, Eigen::RowMajor> symmat_vel = mat_vel + Eigen::SparseMatrix<double, Eigen::RowMajor>(mat_vel.transpose());
    psparsematrix symmat_vel_h2 = matMM2H(symmat_vel);

    return symmat_vel_h2;
}

void writeData2CSV(string file_basename,
    OseenDiscretizationBetter& dc, Eigen::VectorXd& rhs, Eigen::SparseMatrix<double, Eigen::RowMajor> mat, 
    pavector sol_h2, pspmatrix mat_h2, Block_HLU_Prcd* P, int iter, double heps,
    Block_HLU_Times *times, double time_prcd, double time_solve, double time_total, double tol_bicgstab) {
    // Compute distance to e
    Eigen::VectorXd sol_vec = vecH2MM(sol_h2);
    Eigen::VectorXd u_exact = dc.exact_solution;
    int N_ui = dc.N_ui;
    int N_p = dc.N_p;
    double norm_sol = (sol_vec-u_exact).norm()/u_exact.norm();
    double norm_sol_inf = (sol_vec-u_exact).lpNorm<Eigen::Infinity>();
    double norm_sol_u = (sol_vec.head(3*N_ui)-u_exact.head(3*N_ui)).norm()/u_exact.head(3*N_ui).norm();
    double norm_sol_u_inf = (sol_vec.head(3*N_ui)-u_exact.head(3*N_ui)).lpNorm<Eigen::Infinity>();
    double norm_sol_p = (sol_vec.segment(3*N_ui,N_p)-u_exact.segment(3*N_ui,N_p)).norm()/u_exact.segment(3*N_ui,N_p).norm();
    double norm_sol_p_inf = (sol_vec.segment(3*N_ui,N_p)-u_exact.segment(3*N_ui,N_p)).lpNorm<Eigen::Infinity>();
    double norm_res = (rhs-mat*sol_vec).norm()/rhs.norm();
    double norm_luvel = norm2diff_lr_sparsematrix_hmatrix(mat_h2->A, P->A);
    double norm_luschur = norm2diff_lr_schurcomplement_hmatrix_rbffd(P, P->S);

    // Print data to file
    std::stringstream outfile_name;
    outfile_name << file_basename << ".csv";
    // if file doesn't exist yet write header
    if (!std::ifstream(outfile_name.str()))
    {
      print_header_to_file(outfile_name.str());
    }

    std::ofstream outfile(outfile_name.str(), std::ofstream::app);
    auto default_precision = outfile.precision();

    outfile << dc.poly_grad << ";" << dc.poly_div << ";" << dc.poly_conv << ";" << dc.poly_lap << ";" << dc.nu
    << ";" << N_ui << ";" << N_p << ";" << dc.N_dofs << ";" << iter << ";";

    outfile.setf(std::ios::scientific);
    outfile.precision(2);
    outfile << heps << ";";
    outfile.unsetf(std::ios::scientific);
    outfile.precision(default_precision);

    outfile.setf(std::ios::fixed);
    outfile.precision(2);
    outfile << times->velocity_cluster << ";" << times->pressure_cluster << ";" << times->velocity_lu << ";" << times->grad_lower_solve << ";";
    outfile << times->grad_upper_solve << ";" << times->grad_schur_multiplication << ";" << times->schur_computation << ";" << times->schur_lu << ";";
    outfile << time_prcd << ";" << time_solve << ";" << time_total << ";";
    outfile.unsetf(std::ios::fixed);
    outfile.precision(default_precision);

    outfile.setf(std::ios::scientific);
    outfile << norm_sol << ";" << norm_sol_u << ";" << norm_sol_p << ";" 
    << norm_sol_inf << ";" << norm_sol_u_inf << ";" << norm_sol_p_inf << ";"
    << norm_res << ";" << norm_luvel << ";" << norm_luschur << ";" << tol_bicgstab << std::endl;
    outfile.unsetf(std::ios::scientific);

    outfile.close();
}

Block_HLU_Prcd* createBlockHLUPrcd(pspmatrix mat_h2, psparsematrix symmat_vel_h2 , double heps, Block_HLU_Times* times,
    OseenDiscretizationBetter& dc) {
    // Set up HLU_Options
    clustermode mode_vel = DD_ADAPTIVE;
    clustermode mode_p = H2_ADAPTIVE;

    ptruncmode tm = new_releucl_truncmode(); // truncation mode

    void *eta_vel, *eta_grad, *eta_schur; // etas for admissibility condtions
    real eta3 = 16.0;

    adm_sparse_data *sparse_value_vel, *sparse_value_grad;

    sparse_value_vel = new adm_sparse_data;
    sparse_value_vel->nmin = 1;
    sparse_value_vel->sp = mat_h2->A;
    eta_vel = (void *)&sparse_value_vel;

    sparse_value_grad = new adm_sparse_data;
    sparse_value_grad->nmin = 1;
    sparse_value_grad->sp = mat_h2->B[3];
    eta_grad = (void *)&sparse_value_grad;

    eta_schur = (void *)&eta3;

    uint clf_vel = 50; uint clf_p = 30;

    // admissibility conditions
    admissible adm_vel, adm_grad, adm_p;
    adm_vel = admissible_dd_sparse;
    adm_grad = admissible_sparse;
    adm_p = admissible_2_min_cluster;

    HLU_Options *opt = new HLU_Options(mat_h2, clf_vel, clf_p, mode_vel, mode_p,
                                        eta_vel, eta_grad, eta_schur,
                                        adm_vel, adm_grad, adm_p, tm,
                                        heps, heps, heps);

    
    Block_HLU_Prcd* P = new Block_HLU_Prcd(dc.d_u_int, dc.d_div, dc.idxs_p, 
                                            opt, times, symmat_vel_h2, BLOCK_TRIANGULAR);

    delete opt;
    del_truncmode(tm);

    return P;
}

// solve oseen with HLU
void test_strategy(double dx_u, real heps,
    string file_basename, double step_size_scale, string domain_name, double nu, int seed, 
    int poly_grad, int poly_lap, int poly_conv, int poly_div, int conv, int sol, int pressure_constraint)
{
    double dx_p = dx_u * step_size_scale; // seperation distance for pressure nodes

    // set input file for polyhedron
    std::stringstream OFF_File_base;
    OFF_File_base << "OFF_Files/" << domain_name << ".off";
    string OFF_file = OFF_File_base.str();
    PolyhedronShape<Vec3d> shape = PolyhedronShape<Vec3d>::fromOFF(OFF_file);

    // construct domains needed for discretization of Oseen equations and initialize parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    OseenDiscretizationBetter dc(shape, dx_u, dx_p, poly_grad, poly_conv, poly_lap, poly_div, nu, seed); 
    if (pressure_constraint == 0)
        dc.setConstraint(OFF_file); // set pressure constraint, needed to make computepressure at high accuracy
    dc.determineSupports();
    dc.setConvection(conv);
    dc.setSolution(sol, domain_name);

    std::cout << "N_ui = " << dc.N_ui << " and N_p = " << dc.N_p << std::endl;

    // reduced system without dirichlet boundary conditions
    Eigen::SparseMatrix<double, Eigen::RowMajor> mat(dc.N_dofs,dc.N_dofs);
    Eigen::VectorXd rhs(dc.N_dofs);
    createMatrixAndRHS(mat, rhs, dc, pressure_constraint);

    // H2Lib conversions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pspmatrix mat_h2 = matMM2HOseen(mat, dc.N_ui, dc.N_p); // convert to H2Lib-Sparsematrix format
    psparsematrix symmat_vel_h2 = symmetrizeVelocityBlockH2(mat, dc.N_ui); // symmetrized part of velocity block of mat in H2Lib-Sparsematrix format
    pavector rhs_h2 = vecMM2H(rhs); // rhs after eliminating boundary conditions in H2Lib-format
    pavector sol_h2 = new_zero_avector(dc.N_dofs); // vector for numeric solution in H2Lib-format

    // Set up preconditioner and store times %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    auto t4 = high_resolution_clock::now();
    Block_HLU_Times *times = new Block_HLU_Times();
    Block_HLU_Prcd *P = createBlockHLUPrcd(mat_h2, symmat_vel_h2, heps, times, dc);
    auto t5 = high_resolution_clock::now();
    
    // Solve problem %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    real tol_bicgstab = 1e-13; // bicgstab solver relative tolerance
    uint maxiter = 1024; // bicgstab maximal number of iterations
    BiCGStab solver(tol_bicgstab, maxiter);

    auto t6 = high_resolution_clock::now();
    uint iter = solver.psolve_avector_eigen(mat_h2, (addeval_t)addeval_spmatrix_avector_rbffd,
                                        P, rhs_h2, sol_h2);          
    auto t7 = high_resolution_clock::now();

    // post-processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    duration<double> duration_prcd = t5 - t4;    
    double time_prcd = duration_prcd.count();
    duration<double> duration_solve = t7 - t6;
    double time_solve = duration_solve.count();
    double time_total = time_prcd + time_solve;

    writeData2CSV(file_basename, dc, rhs, mat, sol_h2, mat_h2, P, iter, heps, 
        times, time_prcd, time_solve, time_total, tol_bicgstab);
    std::cout << "data written to file " << file_basename << std::endl;

    // clean up %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    delete P;
    delete times;
    del_avector(sol_h2);
    del_avector(rhs_h2);
    del_sparsematrix(symmat_vel_h2);
    del_spmatrix2(mat_h2);
}

int main(int argc, char **argv) {
    #ifdef USE_OPENMP
    int threads = 12; // number of threads
    if (argc > 14)
        threads = atoi(argv[14]);
    omp_set_num_threads(threads);
    std::cout << "Using OpenMP with " << threads << " threads" << std::endl;
    #else
    std::cout << "Using serial implementation" << std::endl;
    #endif

    init_h2lib(&argc, &argv);

    string domain_name[3] = {"cube", "bunny", "narrowing"};

    // accuracy parameter for HLU
    double heps = 0.1;
    if (argc > 1)
        heps = atof(argv[1]);

    // separation distance for velocity nodes
    double dx_u = 1.0/20.0;
    if (argc > 2)
        dx_u = atof(argv[2]);

    // polynomial augmentation degree for laplace, convection, gradient and divergence
    // only degrees 1 to 5 are supported but this can be extended in OseenDiscretizationBetter
    int poly_lap = 4;
    if (argc > 3)
        poly_lap = atoi(argv[3]);

    int poly_conv = 3;
    if (argc > 4)
        poly_conv = atoi(argv[4]);

    int poly_grad = 3;
    if (argc > 5)
        poly_grad = atoi(argv[5]);

    int poly_div = 3;
    if (argc > 6)
        poly_div = atoi(argv[6]);
    
    // choose domain: 0 - cube, 1 - bunny, 2 - narrowing
    int domain_num = 0;
    if (argc > 7)
        domain_num = atoi(argv[7]);

    // viscosity parameter
    double nu = 1e-2;
    if (argc > 8)
        nu = atof(argv[8]);

    // random seed for point generation
    int seed = 0;
    if (argc > 9)
        seed = atoi(argv[9]);

    // separation distance for pressure nodes = step_size_scale * separation distance for velocity nodes
    double step_size_scale = 1.4;
    if (argc > 10)
        step_size_scale = atof(argv[10]);

    // choose convection direction, see OseenDiscretizationBetter::setConvection for details
    int conv = 0;
    if (argc > 11)
        conv = atoi(argv[11]);

    // choose exact solution, see OseenDiscretizationBetter::setSolution for details
    int sol = 3;
    if (argc > 12)
        sol = atoi(argv[12]);

    // 0 - Glaubitz Integration, 1 - set pressure at one node, else - averaging
    int pressure_constraint = 0; 
    if (argc > 13)
        pressure_constraint = atoi(argv[13]);
    
    // set file to write data to
    std::stringstream file_base;
    file_base.setf(std::ios::scientific);
    auto p = file_base.precision();
    file_base.precision(0);
    file_base << "Daten/" << domain_name[domain_num] 
    << "_dc_lvl_lcgd_" << poly_lap <<  poly_conv << poly_grad  << poly_div
    << "_nu_" << nu << "_heps_" << heps  << "_seed_" << seed ;
    file_base.unsetf(std::ios::scientific);
    file_base.precision(p);
    file_base << "_step_size_scale_" << step_size_scale << "_conv_" << conv << "_sol_" << sol << "_pres_" << pressure_constraint;

    // compute for different separation distances
    for (; dx_u > 1.0/21.0; dx_u/=1.2)
        test_strategy(dx_u, heps, file_base.str(), step_size_scale, 
            domain_name[domain_num], nu, seed, poly_grad, poly_lap, poly_conv, poly_div, conv, sol, pressure_constraint);

    uninit_h2lib();

    return 0;
}
