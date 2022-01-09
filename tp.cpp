/******************************************************************************
 * Copyright 1998-2019 Lawrence Livermore National Security, LLC and other
 * HYPRE Project Developers. See the top-level COPYRIGHT file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 ******************************************************************************/

/*
   Example 5

   Interface:    Linear-Algebraic (IJ)

   Compile with: make ex5

   Sample run:   mpirun -np 4 ex5

   Description:  This example solves the 2-D Laplacian problem with zero boundary
                 conditions on an n x n grid.  The number of unknowns is N=n^2.
                 The standard 5-point stencil is used, and we solve for the
                 interior nodes only.

                 This example solves the same problem as Example 3.  Available
                 solvers are AMG, PCG, and PCG with AMG or Parasails
                 preconditioners.  */

#include <math.h>
#include "_hypre_utilities.h"
#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
#include "eigen-3.3.7/unsupported/Eigen/SparseExtra"
#include "eigen-3.3.7/Eigen/Sparse"
#include "vis.c"
#include <iostream>
#include <vector>



int main (int argc, char *argv[])
{
   int i;
   int myid, num_procs;
   int N, n;

   int ilower, iupper;
   int local_size, extra;

   // Parameters
   int max_levels = 25;
   int gamma = 1;
   int nu = 1;
   int relax_type = 3;
   int interp_type = 0;
   int coarsen_type = 6;
   int agg_levels  = 0;
   int vis, print_system;
   double strong_threshold = 0.25;
   double tolerance = 1.e-8;

   double h, h2;

   HYPRE_IJMatrix A;
   HYPRE_ParCSRMatrix parcsr_A;
   HYPRE_IJVector b;
   HYPRE_ParVector par_b;
   HYPRE_IJVector x;
   HYPRE_ParVector par_x;

   HYPRE_Solver solver, precond;

   /* Used load and read a .mtx file */
   const char* filename = nullptr;
   typedef Eigen::SparseMatrix<double, Eigen::RowMajor>SMatrixXf;
   SMatrixXf my_matrix;

   /* Initialize MPI */
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &myid);
   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

   //printf("C > Hello, I am proc number %d/%d\n", myid+1, num_procs);
   //std::cout << "CXX > Hello, I am proc number " << myid+1 << "/" << num_procs << std::endl;

   /* Default problem parameters */
   n = 33;
   vis = 0;
   print_system = 0;

   /* Parse command line */
   {
      int arg_index = 0;
      int print_usage = 0;

      while (arg_index < argc)
      {
         if ( strcmp(argv[arg_index], "-n") == 0 )
         {
            arg_index++;
            n = atoi(argv[arg_index++]);
         }
         else if ( strcmp(argv[arg_index], "-max_levels") == 0 ) {
		 arg_index++; max_levels = atoi(argv[arg_index++]);
         }
         else if ( strcmp(argv[arg_index], "-gamma") == 0 ) {
		 arg_index++; gamma = atoi(argv[arg_index++]);
         }
         else if ( strcmp(argv[arg_index], "-nu") == 0 ) {
		 arg_index++; nu = atoi(argv[arg_index++]);
         }
         else if ( strcmp(argv[arg_index], "-relax_type") == 0 ) {
		 arg_index++; relax_type = atoi(argv[arg_index++]);
         }
         else if ( strcmp(argv[arg_index], "-interp_type") == 0 ) {
		 arg_index++; interp_type = atoi(argv[arg_index++]);
         }
         else if ( strcmp(argv[arg_index], "-coarsen_type") == 0 ) {
		 arg_index++; coarsen_type = atoi(argv[arg_index++]);
         }
         else if ( strcmp(argv[arg_index], "-agg_levels") == 0 ) {
		 arg_index++; agg_levels = atoi(argv[arg_index++]);
         }
         else if ( strcmp(argv[arg_index], "-strong_threshold") == 0 ) {
		 arg_index++; strong_threshold = atof(argv[arg_index++]);
         }
         else if ( strcmp(argv[arg_index], "-tolerance") == 0 ) {
		 arg_index++; tolerance = atof(argv[arg_index++]);
         }
         else if ( strcmp(argv[arg_index], "-vis") == 0 )
         {
            arg_index++;
            vis = 1;
         }
         else if ( strcmp(argv[arg_index], "-print_system") == 0 )
         {
            arg_index++;
            print_system = 1;
         }
         else if ( strcmp(argv[arg_index], "-file") == 0 ){
             arg_index++; filename = argv[arg_index++];
         }
         else if ( strcmp(argv[arg_index], "-help") == 0 )
         {
            print_usage = 1;
            break;
         }
         else
         {
            arg_index++;
         }
      }

      if ((print_usage) && (myid == 0))
      {
         printf("\n");
         printf("Usage: %s [<options>]\n", argv[0]);
         printf("\n");
         printf("  -file                 : Path to csr .mtx file to solve\n");
         printf("  -n <n>                : problem size in each direction (default: 33)\n");
         printf("  -max_levels <n>       : Sets maximum number of multigrid levels (default: 25)\n");
         printf("  -gamma <d>            : Sets AMG gamma (1=V-cycle, 2=W-cycle) (default: 1)\n");
         printf("  -nu <d>               : Sets the number of iteration sweep  (default: 1)\n");
         printf("  -relax_type <d>       : Defines the smoother used  (default: 3)\n");
         printf("  -interp_type <d>      : Defines the interpolation used (default: 0)\n");
         printf("  -coarsen_type <d>     : Defines the coarsening used (default: 6)\n");
         printf("  -agg_levels  <d>      : Defines the number of levels of agg. coarsening (default: 0)\n");
         printf("  -strong_threshold <f> : Sets AMG strength threshold (default: 0.25)\n");
         printf("  -tolerance <f>        : Convergence threshold (default: 1.e-8)\n");
         printf("  -vis                  : save the solution for GLVis visualization\n");
         printf("  -print_system         : print the matrix and rhs\n");
         printf("\n");
      }

      if (print_usage)
      {
         MPI_Finalize();
         return (0);
      }

      //Loading matrix from .mtx file
      //if(myid == 0 && filename != nullptr){
       if(filename != nullptr){
           std::cout << "Loading Matrix from " << filename << "...\n";
           Eigen::loadMarket(my_matrix, filename);

           my_matrix.makeCompressed();
           n = my_matrix.rows();

           std::cout << "Done ! n = " << n << std::endl;
       }
   }

   /* Preliminaries: want at least one processor per row */
   if (n*n < num_procs) n = sqrt(num_procs) + 1;

   /* global number of rows */
   if(filename == nullptr) N = n*n;
   else N = n;

   h = 1.0/(n+1); /* mesh size*/
   h2 = h*h;

   /* Each processor knows only of its own rows - the range is denoted by ilower
      and upper.  Here we partition the rows. We account for the fact that
      N may not divide evenly by the number of processors. */
   local_size = N/num_procs;
   extra = N - local_size*num_procs;

   ilower = local_size*myid;
   ilower += hypre_min(myid, extra);

   iupper = local_size*(myid+1);
   iupper += hypre_min(myid+1, extra);
   iupper = iupper - 1;

   /* How many rows do I have? */
   local_size = iupper - ilower + 1;

   /* Create the matrix.
      Note that this is a square matrix, so we indicate the row partition
      size twice (since number of rows = number of cols) */
   HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &A);

   /* Choose a parallel csr format storage (see the User's Manual) */
   HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);

   /* Initialize before setting coefficients */
   HYPRE_IJMatrixInitialize(A);

   /* Now go through my local rows and set the matrix entries.
      Each row has at most 5 entries. For example, if n=3:

      A = [M -I 0; -I M -I; 0 -I M]
      M = [4 -1 0; -1 4 -1; 0 -1 4]

      Note that here we are setting one row at a time, though
      one could set all the rows together (see the User's Manual).
   */
   if(filename == nullptr)
   {
      int nnz;
      double values[5];
      int cols[5];

      for (i = ilower; i <= iupper; i++)
      {
         nnz = 0;

         /* The left identity block:position i-n */
         if ((i-n)>=0)
         {
            cols[nnz] = i-n;
            values[nnz] = -1.0;
            nnz++;
         }

         /* The left -1: position i-1 */
         if (i%n)
         {
            cols[nnz] = i-1;
            values[nnz] = -1.0;
            nnz++;
         }

         /* Set the diagonal: position i */
         cols[nnz] = i;
         values[nnz] = 4.0;
         nnz++;

         /* The right -1: position i+1 */
         if ((i+1)%n)
         {
            cols[nnz] = i+1;
            values[nnz] = -1.0;
            nnz++;
         }

         /* The right identity block:position i+n */
         if ((i+n)< N)
         {
            cols[nnz] = i+n;
            values[nnz] = -1.0;
            nnz++;
         }

         /* Set the values for row i */
         HYPRE_IJMatrixSetValues(A, 1, &nnz, &i, cols, values);
      }
   }
   //We load matrix from file
   else{
       double* values = my_matrix.valuePtr();
       int* cols = my_matrix.innerIndexPtr();
       int* rows = my_matrix.outerIndexPtr();
       int nnz;

       std::vector<int> my_rows;
//       if(myid > 0 || num_procs == 1)
//            my_rows.assign(rows, rows + local_size + 1);
//       else
            my_rows.assign(rows+ilower, rows + ilower + local_size + 1);

       {
           if(myid==0){
               if(n < 10) std::cout << my_matrix << std::endl;

               std::cout << "Values : [";
               for(int i = 0; i < my_matrix.nonZeros(); ++i){
                   std::cout << values[i];
                   if(i < my_matrix.nonZeros() - 1){
                       std::cout << ", ";
                   }
                   else{
                       std::cout << "]" << std::endl;
                   }
               }

               std::cout << "Cols : [";
               for(int i = 0; i < my_matrix.nonZeros(); ++i){
                   std::cout << cols[i];
                   if(i < my_matrix.nonZeros() - 1){
                       std::cout << ", ";
                   }
                   else{
                       std::cout << "]" << std::endl;
                   }
               }

               std::cout << "Rows : [";
               for(int i = 0; i < n+1; ++i){
                   std::cout << rows[i];
                   if(i < n){
                       std::cout << ", ";
                   }
                   else{
                       std::cout << "]" << std::endl;
                   }
               }
           }

           std::cout << "+++++++++++++++++++++++++\nPROC " << myid << '\n';
           printf("\t>>> ilower : %d, iupper : %d, size : %d\n",    \
                                    ilower, iupper, local_size);
           //std::cout << "> First value : " << values[ilower] << '\n';
           //std::cout << "> First col   : " << cols[ilower] << '\n';
           //std::cout << "> First rows ptr : " << rows[ilower] << '\n';
           std::cout << "> My rows : [";
           int max = local_size+1;

           //if(myid > 0 || num_procs == 1) max++;

           for(int i = 0; i < max; ++i){
               std::cout << my_rows[i];
               if(i < max -1) std::cout << ", ";
               else std::cout << "]\n";
           }

           std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++" \
                        << std::endl;
       }
       return 0;


       HYPRE_IJMatrixSetValues(A, 1, &nnz, &i, cols, values);
   }

   /* Assemble after setting the coefficients */
   HYPRE_IJMatrixAssemble(A);

   /* Note: for the testing of small problems, one may wish to read
      in a matrix in IJ format (for the format, see the output files
      from the -print_system option).
      In this case, one would use the following routine:
      HYPRE_IJMatrixRead( <filename>, MPI_COMM_WORLD,
                          HYPRE_PARCSR, &A );
      <filename>  = IJ.A.out to read in what has been printed out
      by -print_system (processor numbers are omitted).
      A call to HYPRE_IJMatrixRead is an *alternative* to the
      following sequence of HYPRE_IJMatrix calls:
      Create, SetObjectType, Initialize, SetValues, and Assemble
   */


   /* Get the parcsr matrix object to use */
   HYPRE_IJMatrixGetObject(A, (void**) &parcsr_A);


   /* Create the rhs and solution */
   HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper,&b);
   HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR);
   HYPRE_IJVectorInitialize(b);

   HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper,&x);
   HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR);
   HYPRE_IJVectorInitialize(x);

   /* Set the rhs values to h^2 and the solution to zero */
   {
      double *rhs_values, *x_values;
      int    *rows;

      rhs_values =  (double*) calloc(local_size, sizeof(double));
      x_values =  (double*) calloc(local_size, sizeof(double));
      rows = (int*) calloc(local_size, sizeof(int));

      for (i=0; i<local_size; i++)
      {
         rhs_values[i] = h2;
         x_values[i] = 0.0;
         rows[i] = ilower + i;
      }

      HYPRE_IJVectorSetValues(b, local_size, rows, rhs_values);
      HYPRE_IJVectorSetValues(x, local_size, rows, x_values);

      free(x_values);
      free(rhs_values);
      free(rows);
   }


   HYPRE_IJVectorAssemble(b);
   /*  As with the matrix, for testing purposes, one may wish to read in a rhs:
       HYPRE_IJVectorRead( <filename>, MPI_COMM_WORLD,
                                 HYPRE_PARCSR, &b );
       as an alternative to the
       following sequence of HYPRE_IJVectors calls:
       Create, SetObjectType, Initialize, SetValues, and Assemble
   */
   HYPRE_IJVectorGetObject(b, (void **) &par_b);

   HYPRE_IJVectorAssemble(x);
   HYPRE_IJVectorGetObject(x, (void **) &par_x);


  /*  Print out the system  - files names will be IJ.out.A.XXXXX
       and IJ.out.b.XXXXX, where XXXXX = processor id */
   if (print_system)
   {
      HYPRE_IJMatrixPrint(A, "IJ.out.A");
      HYPRE_IJVectorPrint(b, "IJ.out.b");
   }


    double starttime, endtime;
    starttime = MPI_Wtime();

   /* AMG */
   int num_iterations;
   double final_res_norm;

   /* Create solver */
   HYPRE_BoomerAMGCreate(&solver);

   /* Set some parameters (See Reference Manual for more parameters) */
   HYPRE_BoomerAMGSetPrintLevel(solver, 3);  /* print solve info + parameters */
   // HYPRE_BoomerAMGSetOldDefault(solver); /* Falgout coarsening with modified classical interpolaiton */
   HYPRE_BoomerAMGSetRelaxOrder(solver, 1);   /* uses C/F relaxation */

   // Customization
   HYPRE_BoomerAMGSetMaxLevels(solver, max_levels);  /* maximum number of levels */
   HYPRE_BoomerAMGSetCycleType(solver, gamma);       /* Cycle type */
   HYPRE_BoomerAMGSetNumSweeps(solver, nu);          /* Sweeeps on each level */
   HYPRE_BoomerAMGSetRelaxType(solver, relax_type);  /* G-S/Jacobi hybrid relaxation */
   HYPRE_BoomerAMGSetInterpType(solver, interp_type);  /* Coarsening */
   HYPRE_BoomerAMGSetCoarsenType(solver, coarsen_type);  /* Coarsening */
   HYPRE_BoomerAMGSetAggNumLevels(solver, agg_levels);  /* Number of aggressive level */
   HYPRE_BoomerAMGSetStrongThreshold(solver, strong_threshold);
   HYPRE_BoomerAMGSetTol(solver, tolerance);      /* conv. tolerance */

   /* Now setup and solve! */
   HYPRE_BoomerAMGSetup(solver, parcsr_A, par_b, par_x);
   HYPRE_BoomerAMGSolve(solver, parcsr_A, par_b, par_x);

   /* Run info - needed logging turned on */
   HYPRE_BoomerAMGGetNumIterations(solver, &num_iterations);
   HYPRE_BoomerAMGGetFinalRelativeResidualNorm(solver, &final_res_norm);
   if (myid == 0)
   {
      printf("\n");
      printf("Problem size = %d\n", N);
      printf("Iterations = %d\n", num_iterations);
      printf("Final Relative Residual Norm = %e\n", final_res_norm);
      printf("\n");
   }

   endtime = MPI_Wtime();
   if(myid == 0) printf("That took %f seconds\n",endtime-starttime);

   /* Destroy solver */
   HYPRE_BoomerAMGDestroy(solver);

   /* Save the solution for GLVis visualization, see vis/glvis-ex5.sh */
   if (vis)
   {
      FILE *file;
      char filename[255];

      int nvalues = local_size;
      int *rows = (int*) calloc(nvalues, sizeof(int));
      double *values =  (double*) calloc(nvalues, sizeof(double));

      for (i = 0; i < nvalues; i++)
         rows[i] = ilower + i;

      /* get the local solution */
      HYPRE_IJVectorGetValues(x, nvalues, rows, values);

      sprintf(filename, "%s.%06d", "vis/ex5.sol", myid);
      if ((file = fopen(filename, "w")) == NULL)
      {
         printf("Error: can't open output file %s\n", filename);
         MPI_Finalize();
         exit(1);
      }

      /* save solution */
      for (i = 0; i < nvalues; i++)
         fprintf(file, "%.14e\n", values[i]);

      fflush(file);
      fclose(file);

      free(rows);
      free(values);

      /* save global finite element mesh */
      if (myid == 0)
         GLVis_PrintGlobalSquareMesh("vis/ex5.mesh", n-1);
   }

   /* Clean up */
   HYPRE_IJMatrixDestroy(A);
   HYPRE_IJVectorDestroy(b);
   HYPRE_IJVectorDestroy(x);

   /* Finalize MPI*/
   MPI_Finalize();

   return(0);
}
