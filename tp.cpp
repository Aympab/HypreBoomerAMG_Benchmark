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
   int vis, print_system, quiet;
   double strong_threshold = 0.25;
   double tolerance = 1.e-8;

   double h, h2;

   HYPRE_IJMatrix A;
   HYPRE_ParCSRMatrix parcsr_A;
   HYPRE_IJVector b;
   HYPRE_ParVector par_b;
   HYPRE_IJVector x;
   HYPRE_ParVector par_x;

   HYPRE_Solver solver;//, precond; removing this unused variable

   /* Used to load and read a .mtx file */
   const char* filename = nullptr;
   typedef Eigen::SparseMatrix<double, Eigen::RowMajor>SMatrixXf;
   SMatrixXf my_matrix;

   /* Initialize MPI */
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &myid);
   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

   /* Default problem parameters */
   n = 33;
   vis = 0;
   print_system = 0;
   quiet = 0;

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
         else if ( strcmp(argv[arg_index], "-quiet") == 0 )
         {
            arg_index++;
            quiet = 1;
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
         printf("  -quiet                : Shuts Boomer AMG parameters (default:0)\n");
         printf("\n");
      }

      if (print_usage)
      {
         MPI_Finalize();
         return (0);
      }

      //Loading matrix from .mtx file
       if(filename != nullptr){
           if(myid==0){
               std::cout << "Loading Matrix from " << filename << "...\n";
               Eigen::loadMarket(my_matrix, filename);

               n = my_matrix.rows();
               if(n == 0){
                 std::cout << "File doesn't exist or is corrupted" << std::endl;
                 return -1;
               }

               std::cout << "Matrix loaded !" << std::endl;
           }

            //Sending size of matrix to everyone
            MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
     }
   }

   /* Preliminaries: want at least one processor per row */
   if (n*n < num_procs) n = sqrt(num_procs) + 1;

   /* global number of rows
        If we load the matrix from file, N is the same as n,
        to simplify the work (teacher already did that so it's nicely coded)
   */
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


   //Settings arrays to use to split matrix with MPI
   int row_counts[num_procs];
   int row_displacements[num_procs];

   row_displacements[0] = 0;
   //Filling the counts and offset arrays for MPI_ScatterV
   if(myid == 0)
   {

       int plocal_size = N/num_procs;
       int pextra = N - plocal_size*num_procs;

       int plower;
       int pupper;
       int sum_count = 0;

       for(int p = 0; p < num_procs; p++){
           plower=0;
           pupper=0;

           plower = plocal_size*p;
           plower += hypre_min(p, pextra);

           pupper = plocal_size*(p+1);
           pupper += hypre_min(p+1, pextra);
           pupper = pupper - 1;

           //we add +1 because we take nbrows + 1 values for csr
           row_counts[p] = pupper - plower + 1 + 1;

           if(p > 0) row_displacements[p] = sum_count;
           sum_count += row_counts[p]-1;
       }
   }

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

       int my_rows[local_size + 1];

       int local_values_size;

       //Number of values for each proc, used in scatterV
       int v_counts[num_procs];
       int v_displacements[num_procs];

       //Send and receive rows ptr
       if(myid == 0){
           int* rows = my_matrix.outerIndexPtr();

           MPI_Scatterv(
                rows,
                row_counts,
                row_displacements,
                MPI_INT,
                &my_rows,
                local_size+1,
                MPI_INT,        //Datatype receive
                0,              //root MPI ID
                MPI_COMM_WORLD);


            //Displacement for root node is 0
            v_displacements[0] = 0;

            /*
            Here we are going to compute how many values will each proc have
            by using the values stored in rows[] for each proc (we use counts
            and displacement that we used MPI).
            We store this in another counts displacement for CSR values
            */
            int sum_count = 0;
            for(int i = 0; i < num_procs; i++){

                int start = row_displacements[i];
                int stop  = start + row_counts[i] - 1;

                int local_nnz = rows[stop] - rows[start];

                v_counts[i] = local_nnz;

                if(i>0){
                    v_displacements[i] = sum_count;

                    /*Sending local values size to the bros so they can allocate
                    we do this while computing the counts
                    */
                    MPI_Send(&local_nnz,
                         1,
                         MPI_INT,
                         i,
                         666,
                         MPI_COMM_WORLD);
                }

                 sum_count += v_counts[i];
            }

            /*settings size of values for root node so it can allocate
             with the other */
            local_values_size = v_counts[0];

       }
       else{
           //Receiving local rows
           MPI_Scatterv(NULL, NULL, NULL,
                MPI_INT,
                &my_rows,
                local_size+1,
                MPI_INT,        //Datatype receive
                0,              //root MPI ID
                MPI_COMM_WORLD);

           //Receiving size of our arrays of values and cols
           MPI_Recv(&local_values_size, 1, MPI_INT, 0, 666,
                    MPI_COMM_WORLD,
                    MPI_STATUS_IGNORE);
       }

       //Initialize local values and cols now that we have the sizes
       int my_cols[local_values_size];
       double my_values[local_values_size];

       if(myid==0){
           //Send local values
            MPI_Scatterv(my_matrix.valuePtr(),
                 v_counts,
                 v_displacements,
                 MPI_DOUBLE,
                 my_values,
                 local_values_size,
                 MPI_DOUBLE,
                 0,
                 MPI_COMM_WORLD);

             //Sending cols
             MPI_Scatterv(my_matrix.innerIndexPtr(),
                   v_counts,
                   v_displacements,
                   MPI_INT,
                   my_cols,
                   local_values_size,
                   MPI_INT,
                   0,
                   MPI_COMM_WORLD);
       }
       else{
           //Receive local values
            MPI_Scatterv(NULL, NULL, NULL, MPI_DOUBLE,
               my_values, local_values_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

            //Receive local cols
            MPI_Scatterv(NULL, NULL, NULL,
               MPI_INT, my_cols, local_values_size, MPI_INT, 0, MPI_COMM_WORLD);
       }

       /* We are iterating over the local rows ptr for each processes
       We take elements 2 by 2, compute the local number of nnz for this row,
       then we take the values from startindex + i to startindex + i + row_nnz
       and feed this to HYPRE_IJMatrixSetValues
       */
       for(int i = 0; i < local_size; ++i){
           int numrow = ilower + i;

           int startidx = my_rows[i] - my_rows[0];
           int stopidx = my_rows[i+1] - my_rows[0] - 1;
           int row_nnz = stopidx - startidx + 1;

           double* row_values = my_values + startidx;
           int* row_cols = my_cols + startidx;

           HYPRE_IJMatrixSetValues(A, 1, &row_nnz, &numrow, row_cols, row_values);
       }
   }

   /* Assemble after setting the coefficients */
   HYPRE_IJMatrixAssemble(A);

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
   if(quiet == 0)
   {
       HYPRE_BoomerAMGSetPrintLevel(solver, 3);  /* print solve info + parameters */
       HYPRE_BoomerAMGSetOldDefault(solver); /* Falgout coarsening with modified classical interpolaiton */
   }

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
