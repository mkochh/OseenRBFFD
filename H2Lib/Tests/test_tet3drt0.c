
#include "tet3d.h"		/* 3-dimensional mesh */
#include "tet3drt0.h"		/* discretisation with Raviart-Thomas functions */
#include "sparsematrix.h"	/* Sparsematrices */
#include "krylov.h"		/* Iterative solvers of Krylov type */

#include <stdio.h>

static uint problems = 0;

#define IS_IN_RANGE(a, b, c) (((a) <= (b)) && ((b) <= (c)))

/*****************************************
 *Boundary conditions for domain
 ******************************************/
void set_boundary_unitcube_tet3drt0(ptet3drt0 rt0){
  
  pctet3d t3 = rt0->t3;
  uint faces = t3->faces;
  real (*x)[3] = t3->x;
  real xt[3][3];
  uint v[3];
  uint i, j;
  
  for(i=0; i<faces;i++){
    /*Get vertices of face i*/
    getvertices_face_tet3d(t3, i, v);
   /*Get coordinates of the vertices of face i*/
    for(j=0; j<3;j++){
      xt[j][0] = x[v[j]][0]; 
      xt[j][1] = x[v[j]][1];
      xt[j][2] = x[v[j]][2];
     }
    /*start and end vertex of edge i*/
    /*z=0*/
    if(xt[0][2]==0 && xt[1][2]==0 && xt[2][2]==0)
      rt0->is_dof[i] = 2;
    /*z=1*/
    else if(xt[0][2]==1 && xt[1][2]==1 && xt[2][2]==1)
      rt0->is_dof[i] = 2;
    /*y=0*/
    else if(xt[0][1]==0 && xt[1][1]==0 && xt[2][1]==0)
      rt0->is_dof[i] = 2;
    /*y=1*/
    else if(xt[0][1]==1 && xt[1][1]==1 && xt[2][1]==1)
      rt0->is_dof[i] = 2;
  }
}

/**********************************************
 *exact solution for simple model problem
 **********************************************/
field function_pressure(const real *x, void *data){
  
  return -1.0 * x[0] + 1.0;
}

void function_flux(const real *x, void *data, pavector v){

  v->v[0] = 1.0;
  v->v[1] = 0.0;
  v->v[2] = 0.0;
}

/*************************************************
 *Sets values for the permeability
 *************************************************/
void function_permeability(ptet3d t3, pavector k)
{
  uint i;
  
  for(i=0;i<k->dim;i++)
    k->v[i] = 1.0;
}

/********************************************
 *Functions for rhs and boundary conditions
 ********************************************/
field g_D(const real *x, void *fdata){
  
  uint v;

  v = 10;

  if(x[0] == 0.0) v = 1;
  else if(x[0] == 1.0) v = 0;
  else printf("error\n");

 return v; 
}

field g_N(const uint *e, void *data)
{
 return 0.0; 
}


field f(const real *e, void *data)
{
  
 field v;

  v =  0.0;
  
  return v;
}
/**************************************
 *Auxiliary structure for solver
 **************************************/
typedef struct _prcg prcg;
typedef prcg *pprcg;
struct _prcg {
  /** @brief Sparsematrix.*/
  psparsematrix A;
  /** @brief Accuracy or the cg-iteration*/
  real eps;
  /** @brief Number of maximal steps.*/
  uint steps;
};

pprcg new_prcg(psparsematrix A, real accuracy, uint steps)
{
 pprcg p; 
 p = (pprcg) allocmem(sizeof(prcg));
 p->A = A;
 p->eps = accuracy;
 p->steps = steps;
 return p; 
}

void del_prcg(pprcg p){

  freemem(p);
}

/***********************************
 *Solver
 ************************************/
/*CG-Solver for inner system*/
void solve_A(void *pdata, pavector b)
{
  addeval_t addeval_A;
  pavector x;
  pavector r, p, a;
  uint n, i;
  pprcg pr;
  real accuracy;
  uint steps;
  psparsematrix A;
  real norm;
  
  pr = (pprcg) pdata;
  A = pr->A;
  accuracy = pr->eps;
  steps = pr->steps;
  addeval_A = (addeval_t) addeval_sparsematrix_avector;
  n = b->dim;
  x = new_avector(n);
  random_avector(x);
  r = new_avector(n);
  p = new_avector(n);
  a = new_avector(n);
  
  init_cg(addeval_A, A, b, x, r, p, a);
  for(i=0;i<steps;i++){
      step_cg(addeval_A, A, b, x, r, p,a);
      norm = norm2_avector(r);
     if(norm < accuracy)break;
    
  }
  copy_avector(x,b);
  
  del_avector(x);
  del_avector(r);
  del_avector(p);
  del_avector(a);
  
}

/*Uzawa-Solver for Schur-System*/
uint solve_uzawa_sparsematrix(psparsematrix A, psparsematrix B, pavector b1, pavector b2,
                                      pavector x1, pavector x2, real accuracy, uint steps)
{
  prcd_t prcd_A;
  mvm_t mvm_B;
  pavector r2, p2, a1, s2;
  uint n_b1, n_b2;
  uint i;
  real norm;
  pprcg pr;
  
  norm = 0.0;
  prcd_A = (prcd_t) solve_A;
  mvm_B = (mvm_t) mvm_sparsematrix_avector;
  pr = new_prcg(A, accuracy, steps);
  n_b1 = b1->dim;
  n_b2= b2->dim;
  r2 = new_avector(n_b2);
  p2 = new_avector(n_b2);
  a1 = new_avector(n_b1);
  s2 = new_avector(n_b2);

  init_uzawa(prcd_A, pr, mvm_B, B, b1, b2, x1, x2, r2, p2, a1, s2);
  for(i=0;i<steps;i++){
    step_uzawa(prcd_A, pr, mvm_B, B, b1, b2, x1, x2, r2, p2, a1, s2);
    norm = norm2_avector(r2);
    if(norm < accuracy) break;
  }
  printf("\n");

  del_prcg(pr);
  del_avector(r2);
  del_avector(p2);
  del_avector(a1);
  del_avector(s2);
  
  return i;
}

int 
main (int argc, char **argv){
  
  ptet3d *gr;		/* 3d mesh hierarchy */
  ptet3drt0 *dc;	/* Raviart-Tomas basis function in 3d */
  psparsematrix sp_A, sp_Af, sp_B, sp_Bf;/* Sparsematrix objects */
  uint i;		/* Auxiliary variable for loops */
  uint L;		/* Numberof grid refinements */
  pavector k;		/* Vector for storing the permeabilities */
  pstopwatch sw;	/* Stopwatch for time measuring */
  real time;		/* Variable for time measuring */
  real accuracy;	/* Accuracy for solver */
  real error;		/* Auxiliary variable for error computation */
  uint dim_A, rows_B;	/* Auxiliary variables for matrix dimensions */
  uint steps;		/* Auxiliary variable for number of steps */
  uint max_steps;	/* Maximal number of steps */
  pavector b1, b2, g;	/* Vectors for the right-hand sides and the Neumann values */
  pavector x1, x2;	/* Vectors for the calculated solution */
  
  init_h2lib(&argc, &argv);
  
  sw = new_stopwatch();
  
  L = 4;
  max_steps = 1000;
  accuracy = 1e-12;
  
  (void) printf("----------------------------------------\n"
                "Testing tet3drt0\n");
  (void) printf("----------------------------------------\n"
		"Creating mesh hierarchy\n");
  gr = (ptet3d *) allocmem(sizeof(ptet3d) * (L + 1));

  gr[0] = new_unitcube_tet3d();
  (void) printf("  Level %2u: %u vertices, %u edges, %u faces, %u tetrahedra\n",
		0, gr[0]->vertices, gr[0]->edges, gr[0]->faces, gr[0]->tetrahedra);
  for (i = 0; i < L; i++) {
    gr[i + 1] = refine_tet3d(gr[i], 0);
    fixnormals_tet3d(gr[i]);
    (void) printf("  Level %2u: %u vertices, %u edges, %u faces, %u tetrahedra\n",
		  i + 1, gr[i + 1]->vertices, gr[i + 1]->edges,
		  gr[i + 1]->faces, gr[i + 1]->tetrahedra);
  }
  fixnormals_tet3d(gr[L]);
   
  (void) printf("Creating discretizations\n");
  dc = (ptet3drt0 *) allocmem(sizeof(ptet3drt0) * (L + 1));
  for (i = 0; i <= L; i++) {
    dc[i] = new_tet3drt0(gr[i]);
    (void) printf("  Level %2u: %u degrees of freedom, %u fixed\n",
		  i, dc[i]->ndof, dc[i]->nfix);
  }
  
  (void) printf("Setting boundary conditions\n");
  for(i = 0; i <= L; i++){
    set_boundary_unitcube_tet3drt0(dc[i]);  
    update_tet3drt0(dc[i]);
  }  
  
  for(i = 2; i <= L; i++){
    (void) printf("Testing level %u\n" "  Setting up matrix\n", i);
    start_stopwatch(sw);
    sp_A = build_tet3drt0_A_sparsematrix(dc[i]); //printf("  sp_A rows %u cols %u \n", sp_A->rows, sp_A->cols);
    sp_B = build_tet3drt0_B_sparsematrix(dc[i]); //printf("  sp_B rows %u cols %u \n", sp_B->rows, sp_B->cols);
    sp_Af = build_tet3drt0_A_interaction_sparsematrix(dc[i]);
    sp_Bf = build_tet3drt0_B_interaction_sparsematrix(dc[i]);
    k = new_avector(gr[i]->tetrahedra);
    function_permeability(gr[i], k);
    assemble_tet3drt0_darcy_A_sparsematrix(dc[i], sp_A, sp_Af, k);
    assemble_tet3drt0_darcy_B_sparsematrix(dc[i], sp_B, sp_Bf);
    time = stop_stopwatch(sw);
    (void) printf("  %.6f seconds\n"
		  "  sp_A:\n  %.1f MB (interaction %.1f MB)\n"
		  "  %.1f KB/DoF (interaction %.1f KB/DoF)\n"
		  "  %u non-zero entries (interaction %u)\n"
		  "  sp_B:\n  %.1f MB (interaction %.1f MB)\n"
		  "  %.1f KB/DoF (interaction %.1f KB/DoF)\n"
		  "  %u non-zero entries (interaction %u)\n",
		  time,
		  getsize_sparsematrix(sp_A) / 1048576.0,
		  getsize_sparsematrix(sp_Af) / 1048576.0,
		  getsize_sparsematrix(sp_A) / 1024.0 / dc[i]->ndof,
		  getsize_sparsematrix(sp_Af) / 1024.0 / dc[i]->ndof,
		  sp_A->nz, sp_Af->nz,
		  getsize_sparsematrix(sp_B) / 1048576.0,
		  getsize_sparsematrix(sp_Bf) / 1048576.0,
		  getsize_sparsematrix(sp_B) / 1024.0 / dc[i]->ndof,
		  getsize_sparsematrix(sp_Bf) / 1024.0 / dc[i]->ndof,
		  sp_B->nz, sp_Bf->nz);
    
    dim_A = sp_A->rows;
    rows_B = sp_B->rows;
  
    (void) printf("  Setting up b_D\n");
    b1 = new_zero_avector(dim_A);
    assemble_tet3drt0_b_D_avector(dc[i], g_D, dc[i] , b1); 
  
    (void) printf("  Setting up b_f\n");
    b2 = new_avector(rows_B);
    assemble_tet3drt0_b_f_avector(dc[i], f, 0, b2); /*constant zero*/
 
    (void) printf("  Setting up g (Neumann) \n");
    g = new_zero_avector(dc[i]->nfix); 
    assemble_tet3drt0_g_N_avector(dc[i], g_N, 0, g);
    addeval_sparsematrix_avector(-1.0, sp_Af, g, b1);
    addeval_sparsematrix_avector(-1.0, sp_Bf, g, b2);
    
    (void) printf("  Starting iteration");
    x1 = new_avector(dim_A);
    x2 = new_avector(rows_B);
    random_avector(x2);
    start_stopwatch(sw);
    steps = solve_uzawa_sparsematrix(sp_A, sp_B, b1, b2, x1, x2, accuracy, max_steps);
    time = stop_stopwatch(sw);
    (void) printf("  %u iterations\n", steps);

    (void) printf("\n"
		  "  %.1f seconds\n"
		  "  %.1f seconds per step\n", time, time / steps);

    error = norml2_pressure_centroid_tet3drt0(dc[i], function_pressure, 0, x2);
    (void) printf("  rel. L^2 error pressure %.4e %s\n", error,
		  (IS_IN_RANGE(1.0e-13, error, 1.0e-11) ? "    okay" :
		   "NOT okay"));
    if (!IS_IN_RANGE(1.0e-13, error, 1.0e-12))
      problems++;
      
    error = norml2_flux_centroid_tet3drt0(dc[i], function_flux, 0, x1, g);
    (void) printf("  rel. L^2 error flux %.4e     %s\n", error,
		  (IS_IN_RANGE(1.0e-13, error, 5.0e-10) ? "    okay" :
		   "NOT okay"));
    if (!IS_IN_RANGE(1.0e-13, error, 5.0e-10))
      problems++;
    
    del_sparsematrix(sp_A);
    del_sparsematrix(sp_B);
    del_sparsematrix(sp_Af);
    del_sparsematrix(sp_Bf);
    del_avector(k);
    del_avector(b1);
    del_avector(b2);
    del_avector(g);
    del_avector(x1);
    del_avector(x2);
  }
  
  (void) printf("----------------------------------------\n" "Cleaning up\n");
  for(i = 0; i <= L; i++)
    del_tet3drt0(dc[i]);
  freemem(dc);
  for(i = 0; i <=L; i++)
    del_tet3d(gr[i]);
  freemem(gr);
  del_stopwatch(sw);

  uninit_h2lib();
  
  return problems;
}