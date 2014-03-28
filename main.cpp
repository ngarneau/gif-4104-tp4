#include "Matrix.hpp"
#include "simpleCL.h"
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <stdexcept>

using namespace std;


// Inverser la matrice par la méthode de Gauss-Jordan; implantation séquentielle.
void invertSequential(Matrix& iA) {

    // vérifier que la matrice est carrée
    assert(iA.rows() == iA.cols());
    // construire la matrice [A I]
    MatrixConcatCols lAI(iA, MatrixIdentity(iA.rows()));

    // traiter chaque rangée
    for (size_t k=0; k<iA.rows(); ++k) {
        // trouver l'index p du plus grand pivot de la colonne k en valeur absolue
        // (pour une meilleure stabilité numérique).
        size_t p = k;
        double lMax = fabs(lAI(k,k));
        for(size_t i = k; i < lAI.rows(); ++i) {
            if(fabs(lAI(i,k)) > lMax) {
                lMax = fabs(lAI(i,k));
                p = i;
            }
        }
        // vérifier que la matrice n'est pas singulière
        if (lAI(p, k) == 0) throw runtime_error("Matrix not invertible");

        // échanger la ligne courante avec celle du pivot
        if (p != k) lAI.swapRows(p, k);

        double lValue = lAI(k, k);
        for (size_t j=0; j<lAI.cols(); ++j) {
            // On divise les éléments de la rangée k
            // par la valeur du pivot.
            // Ainsi, lAI(k,k) deviendra égal à 1.
            lAI(k, j) /= lValue;
        }

        // Pour chaque rangée...
        for (size_t i=0; i<lAI.rows(); ++i) {
            if (i != k) { // ...différente de k
                // On soustrait la rangée k
                // multipliée par l'élément k de la rangée courante
                double lValue = lAI(i, k);
                lAI.getRowSlice(i) -= lAI.getRowCopy(k)*lValue;
            }
        }
    }

    // On copie la partie droite de la matrice AI ainsi transformée
    // dans la matrice courante (this).
    for (unsigned int i=0; i<iA.rows(); ++i) {
        iA.getRowSlice(i) = lAI.getDataArray()[slice(i*lAI.cols()+iA.cols(), iA.cols(), 1)];
    }
}

// Inverser la matrice par la méthode de Gauss-Jordan; implantation MPI parallèle.
void invertParallel(Matrix& iA) {

   size_t global_size[2], local_size[2];
   int found;
   sclHard hardware;
   sclSoft software, software2;
   int n = iA.rows();

   MatrixConcatCols lAI(iA, MatrixIdentity(iA.rows()));

   int worksize = lAI.rows()*lAI.cols()*sizeof(double);
   int numberOfCols = lAI.cols();

   // Get the hardware
   hardware = sclGetCPUHardware( 0, &found );

   // Get the software
   software = sclGetCLSoftware( "example.cl", "findBiggest", hardware );
   software2 = sclGetCLSoftware("normalize.cl", "normalisation", hardware);

   // Setting NDRanges
   /*std::valarray<double> plainData2 = iA.getDataArray();
   global_size[0] = iA.rows(); global_size[1] = iA.cols();
   local_size[0] = 1; local_size[1] = 1;
   double test = 0;
   sclManageArgsLaunchKernel( hardware, software2, global_size, local_size,
        "%R %R",
        worksize, &plainData2[0], sizeof(double), &test, sizeof(double));
   cout << test << '\n';*/

  for (size_t k=0; k < lAI.rows(); ++k) {
    
      int maximumIndexValue = 0;
      std::valarray<double> plainData = lAI.getDataArray();
      global_size[0] = lAI.rows()*lAI.cols(); global_size[1] = 1;
      local_size[0] = global_size[0]; local_size[1] = 1;
      int max = 0;
      
      sclManageArgsLaunchKernel( hardware, software, global_size, local_size,
        "%r %R %r %r %R",
        worksize, &plainData[0], sizeof(double), &maximumIndexValue, sizeof(int), &k, sizeof(int), &numberOfCols, sizeof(int), &max);

      if (maximumIndexValue != k) lAI.swapRows(maximumIndexValue, k);

      double lValue = lAI(k, k);
      for (size_t j=0; j<lAI.cols(); ++j) {
          lAI(k, j) /= lValue;
      }

      //cout << "\nMatrice:\n" << lAI.str() << endl;

      // here we pass in the whole matrix, with its extension
      std::valarray<double> plainData2 = lAI.getDataArray();
      global_size[0] = lAI.rows(); global_size[1] = lAI.cols();
      local_size[0] = global_size[0]; local_size[1] = global_size[1];
      int worksize2 = lAI.rows()*lAI.cols()*sizeof(double);
      int nbrItems = lAI.rows()*lAI.cols();
      double test = 0;
      sclManageArgsLaunchKernel( hardware, software2, global_size, local_size,
        "%R %R %r %r %r",
        worksize2, &plainData2[0], sizeof(double), &test, sizeof(int), &k, sizeof(int), &n, sizeof(int), &nbrItems);

      /*for (size_t i=0; i<lAI.rows(); ++i) {
          if (i != k) { // ...différente de k
              // On soustrait la rangée k
              // multipliée par l'élément k de la rangée courante
              double lValue = lAI(i, k);
              lAI.getRowSlice(i) -= lAI.getRowCopy(k)*lValue;
          }
      }*/
  }



    /*int myrank, ranksize;
    ranksize = MPI::COMM_WORLD.Get_size(); // nombre de processus
    myrank = MPI::COMM_WORLD.Get_rank(); // numero du processus courant (me)
    int n = iA.rows();

    // vérifier que la matrice est carrée
    assert(iA.rows() == iA.cols());
    // construire la matrice [A I]
    MatrixConcatCols lAI(iA, MatrixIdentity(iA.rows()));

    //for (size_t k=0; k < 1; ++k) {
    for (size_t k=0; k < lAI.rows(); ++k) {

        // on trouve le q localement (par processus)
        double lMax = 0.0;
        size_t q = k;
        for(size_t i = k; i < lAI.rows(); ++i) {
            if( (i % ranksize) == myrank) {

                if(fabs(lAI(i,k)) > lMax) {
                    lMax = fabs(lAI(i,k));
                    q = i;
                }
            }
        }

        data in, out;
        in.index = q;
        in.value = lMax;

        MPI::COMM_WORLD.Allreduce(&in, &out, ranksize, MPI_FLOAT_INT, MPI_MAXLOC);
        q = out.index;
        int root = q%ranksize;

        // broadcast a tous les processus les elements de k a n-1 de l'indice q trouve
        MPI::COMM_WORLD.Bcast(&lAI(q,0), lAI.cols(), MPI::DOUBLE, root);

        // vérifier que la matrice n'est pas singulière
        if (lAI(q, k) == 0) throw runtime_error("Matrix not invertible");

        // on swap la ligne q avec la ligne k
        if (q != k) lAI.swapRows(q, k);

        // on normalise la ligne k afin que l'element (k,k) soit egale a 1
        double lValue = lAI(k, k);
        for (size_t j=0; j<lAI.cols(); ++j) {
            lAI(k, j) /= lValue;
        }

        //// Pour chaque rangée...
        for (size_t i=0; i<lAI.rows(); ++i) {
            if( (i % ranksize) == myrank) {
                if (i != k) { // ...différente de k
                    // On soustrait la rangée k
                    // multipliée par l'élément k de la rangée courante
                    double lValue = lAI(i, k);
                    lAI.getRowSlice(i) -= lAI.getRowCopy(k)*lValue;
                }
            }
        }

        for(size_t i = 0; i < lAI.rows(); ++i) {
            MPI::COMM_WORLD.Bcast(&lAI(i,0), lAI.cols(), MPI::DOUBLE, i%ranksize);
        }

    }


    //if(myrank == 1) {
        //cout << "\nMatrice parallele du processus " << myrank << "\n" << lAI.str() << endl;
    //}
    //MPI::COMM_WORLD.Barrier();
    */

    // On copie la partie droite de la matrice AI ainsi transformée
    // dans la matrice courante (this).
    for (unsigned int i=0; i<iA.rows(); ++i) {
        iA.getRowSlice(i) = lAI.getDataArray()[slice(i*lAI.cols()+iA.cols(), iA.cols(), 1)];
    }


}


// Multiplier deux matrices.
Matrix multiplyMatrix(const Matrix& iMat1, const Matrix& iMat2) {

    // vérifier la compatibilité des matrices
    assert(iMat1.cols() == iMat2.rows());
    // effectuer le produit matriciel
    Matrix lRes(iMat1.rows(), iMat2.cols());
    // traiter chaque rangée
    for(size_t i=0; i < lRes.rows(); ++i) {
        // traiter chaque colonne
        for(size_t j=0; j < lRes.cols(); ++j) {
            lRes(i,j) = (iMat1.getRowCopy(i)*iMat2.getColumnCopy(j)).sum();
        }
    }
    return lRes;
}

int main(int argc, char** argv) {


   // get input from user
   unsigned int lS = 5;
   if (argc >= 2) {
      lS = atoi(argv[1]);
   }

   int workUnits = 8;
   if (argc >= 3) {
      workUnits = atoi(argv[2]);
   }

   // on initialise une matrix sequential et parallel pour comparaison
   Matrix mainMatrixSequential = Matrix(lS, lS);
   Matrix mainMatrixParallel = Matrix(lS, lS);
   Matrix lA = Matrix(lS, lS);


   MatrixRandom lARandom(lS, lS);
   lA = lARandom;
   mainMatrixSequential = lA;

   // On traite la matrice en séquentiel
   invertSequential(mainMatrixSequential);
   //printf( "Elapsed time is %f\n", t2 - t1 );

   //cout << "Matrice inverse:\n" << mainMatrixSequential.str() << endl;
   Matrix lRes = multiplyMatrix(lA, mainMatrixSequential);
   //cout << "Produit des deux matrices:\n" << lRes.str() << endl;
   //cout << "Erreur: " << lRes.getDataArray().sum() - lS << endl;

   // on créer la matrice parallèle
   mainMatrixParallel = lA;
   cout << "Matrice:\n" << mainMatrixParallel.str() << endl;
   invertParallel(mainMatrixParallel);
   cout << "Matrice:\n" << mainMatrixParallel.str() << endl;

   //printf( "Elapsed time is %f\n", t2 - t1 );
   //cout << "\nMatrice inverse parallele:\n" << mainMatrixParallel.str() << endl;
   //Matrix lRes = multiplyMatrix(lA, mainMatrixParallel);
   //cout << "Produit des deux matrices:\n" << lRes.str() << endl;
   //cout << "Erreur: " << lRes.getDataArray().sum() - lS << endl;


   /*char buf[]="Hello, World!";
   size_t global_size[2], local_size[2];
   int found, worksize;
   sclHard hardware, *hardwares;
   sclSoft software;

   // Target buffer just so we show we got the data from OpenCL
   worksize = strlen(buf);
   char buf2[worksize];
   buf2[0]='?';
   buf2[worksize]=0;

   int foundHardwares;
   hardwares = sclGetAllHardware( &foundHardwares );
    
   // Get the hardware
   hardware = sclGetCPUHardware( 0, &found );
   // Get the software
   software = sclGetCLSoftware( "example.cl", "example", hardware );
   // Set NDRange dimensions
   global_size[0] = strlen(buf); global_size[1] = 1;
   local_size[0] = global_size[0]; local_size[1] = 1;
    
   sclManageArgsLaunchKernel( hardware, software, global_size, local_size,
                               " %r %w ",
                              worksize, buf, worksize, buf2 );
    
   // Finally, output out happy message.
   puts(buf2);*/

}
