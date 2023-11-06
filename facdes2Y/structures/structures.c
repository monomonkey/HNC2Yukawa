#include "structures.h"
#include "../math/math_aux.h"

void input(double fv, double xnu, species especie1, species especie2, double rmax, int potentialID) {

    int i;
    double dq;

    rho = (6.0 / M_PI) * fv;

    double sigma1 = especie1.diameter;
    double sigma2 = especie2.diameter;

    sigmaVec[0] = sigma1;
    sigmaVec[1] = (sigma1 + sigma2) / 2.0;
    sigmaVec[2] = sigma2;

    dq = M_PI / rmax;
    dr = rmax / ((double) nrows);

    for (i = 0; i < nrows; i++) {
        r[i] = i * dr;
        q[i] = i * dq;

//        r[i] = (i+1) * dr;
//        q[i] = (i+1) * dq;

    }

    POT(especie1, especie2, potentialID, xnu);
}

void POT(species especie1, species especie2, int potentialID, double xnu) {

    int i, k;
    int l_atr;
    double dmed, rlamb;
    double arg1, arg2, arg3, arg4;
    double *Ua, *Ur, *E, *E2, *z, *z2;

    // Allocate memory for an array of nrows doubles
    Ua = malloc(nrows*ncols * sizeof(double));
    Ur = malloc(nrows*ncols * sizeof(double));
    E = malloc(ncols * sizeof(double));
    E2 = malloc(ncols * sizeof(double));
    z = malloc(ncols * sizeof(double));
    z2 = malloc(ncols * sizeof(double));

    // Check if memory allocation was successful
    if (Ua == NULL || Ur == NULL || E == NULL || E2 == NULL || z == NULL || z2 == NULL ) {
        printf("Memory allocation failed 1.\n");
        return; // Exit with an error code
    }

    // dmed es la distancia media entre partículas.
    dmed = pow(rho, -1.0/3.0);
//    printf("=====================\n");
//    printf(" rho = %.10lf \n", rho);
//    printf(" dmed = %.10lf \n", dmed);
//    printf("=====================\n\n");

    // Se despliega el catálogo de potenciales que podemos utilizar
    
    switch(potentialID){
        
        case 1:

            //*************************************************************
            //--------------------- INVERSE POWER LAW ---------------------
            //
            //                      U = T* (σ/r)^(λ)
            //
            //*************************************************************

            E[0] = 1.0 / especie1.temperature;  //Temperature*
            E[2] = 1.0 / especie2.temperature;  //Temperature*
            E[1] = sqrt(E[0] * E[2]);

            z[0] = especie1.lambda;
            z[2] = especie2.lambda;
            z[1] = sqrt(z[0] * z[2]);

            for (k = 0; k < ncols; k++) {
                for (i = 0; i < nrows; i++) {
                    if (r[i] < (sigmaVec[k] / 2.0)) {
                        U[i*ncols + k] = 0.0;
                        Up[i*ncols + k] = 0.0;
                    } else {
                        U[i*ncols + k] = E[k] * pow(sigmaVec[k]/r[i], z[k]);
                        Up[i*ncols + k] = U[i*ncols + k] * (z[k]); // = -f(r)*r

                        // Up[i*ncols + k] = U[i*ncols + k] * (rlamb / r[i]); // = -f(r)
                    }
                }
            }

            break;
        
        case 2: 

            //*************************************************************
            //------- LENNARD JONES TRUNCADO (REPULSIVO SOLAMENTE) --------
            //
            //           U = T* 4 [(σ/r)^(12) - (σ/r)^(6) + 1/4]
            //
            //*************************************************************

            rlamb = 6.0;

            E[0] = 1.0 / especie1.temperature;  //Temperature*
            E[2] = 1.0 / especie2.temperature;  //Temperature*
            E[1] = sqrt(E[0] * E[2]);

            for (k = 0; k < ncols; k++) {
                for (i = 0; i < nrows; i++) {
                    arg4 = sigmaVec[k] * pow(2.0, 1.0/rlamb);
                    if ((r[i] < (sigmaVec[k] / 2.0)) || (r[i] > arg4)) {
                        U[i*ncols + k] = 0.0;
                        Up[i*ncols + k] = 0.0;
                    } else {
                        arg1 = pow(sigmaVec[k] / r[i], rlamb);
                        arg2 = arg1 * arg1;
                        arg3 = 1.0/4.0;
                        U[i*ncols + k] = 4.0 * E[k] * (arg2 - arg1 + arg3);
                        Up[i*ncols + k] = 4.0 * E[k] * rlamb * (2.0*arg2 - arg1);
                    }
                }
            }

            break;

        case 3: 

            //*************************************************************
            //------- LENNARD JONES TRUNCADO CON MÍNIMO EN r = σ ----------
            //
            //            U = T* [(σ/r)^(2λ) - (σ/r)^(λ) + 1]
            //
            //*************************************************************
                
            E[0] = 1.0 / especie1.temperature;  //Temperature*
            E[2] = 1.0 / especie2.temperature;  //Temperature*
            E[1] = sqrt(E[0] * E[2]);
            
            for (k = 0; k < ncols; k++) {
                for (i = 0; i < nrows; i++) {
                    if ((r[i] < (sigmaVec[k] / 2.0)) || (r[i] > sigmaVec[k])) {
                        U[i*ncols + k] = 0.0;
                        Up[i*ncols + k] = 0.0;
                    } else {
                        arg1 = pow(sigmaVec[k] / r[i], xnu);
                        arg2 = arg1 * arg1;
                        arg3 = 1.0;
                        U[i*ncols + k] = E[k] * (arg2 - 2.0*arg1 + arg3);
                        Up[i*ncols + k] = E[k] * 2.0*xnu * (arg2 - arg1);
                    }
                }
            }

            break;

        case 4: // ####---- Double Yukaka (atractive + repulsive) ----####

            //*************************************************************
            //----------- DOBLE YUKAWA (ATRACTIVO Y REPULSIVO) ------------
            //
            //     U = {-T₁* exp[-λ₁ (r-1)]/r} + {T₂* exp[-λ₂ (r-1)]/r}
            //
            //*************************************************************
                
            //rlamb = 0.149;
            //l_atr = 1.5;

            E[0] = 1.0 / especie1.temperature;  //Temperature₁*
            E[2] = 1.0 / especie2.temperature;  //Temperature₁*
            E[1] = sqrt(E[0] * E[2]);

            E2[0] = 1.0 / especie1.temperature2;  //Temperature₂*
            E2[2] = 1.0 / especie2.temperature2;  //Temperature₂*
            E2[1] = sqrt(E2[0] * E2[2]);

            z[0] = especie1.lambda;
            z[2] = especie2.lambda;
            z[1] = sqrt(z[0] * z[2]);

            z2[0] = especie1.lambda2;
            z2[2] = especie2.lambda2;
            z2[1] = sqrt(z2[0] * z2[2]);

            for (k = 0; k < ncols; k++) {
                for (i = 0; i < nrows; i++) {
                    if (r[i] < sigmaVec[k]) {
                        Ua[i*ncols + k] = 0.0;
                        Ur[i*ncols + k] = 0.0;
                        U[i*ncols + k] = 0.0;
                        Up[i*ncols + k] = 0.0;
                    } else {
                        Ua[i*ncols + k] = - E[k] * exp(- z[k] * (r[i] - 1.0)) / r[i];  // Atractive Potential
                        Ur[i*ncols + k] = E2[k] * exp(- z2[k] * (r[i] - 1.0)) / r[i];  // Repulsive Potential
                        U[i*ncols + k] = Ua[i*ncols + k] + Ur[i*ncols + k];            // Total Potential
                        Up[i*ncols + k] = (1.0 + z[k]*r[i]) * Ua[i*ncols + k] + \
                                          (1.0 + z2[k]*r[i]) * Ur[i*ncols + k];        // -[d U(r)/dr] * r
                    }
                }
            }
          
            break;

        case 5: // ####---- Atractive Yukaka ----####

            //*************************************************************
            //--------------------- ATRACTIVE YUKAWA ----------------------
            //
            //                   U = -T* exp[-λ (r-1)]/r                  
            //
            //*************************************************************

            E[0] = 1.0 / especie1.temperature;  //Temperature*
            E[2] = 1.0 / especie2.temperature;  //Temperature*
            E[1] = sqrt(E[0] * E[2]);

            z[0] = especie1.lambda;
            z[2] = especie2.lambda;
            z[1] = sqrt(z[0] * z[2]);

            for (k = 0; k < ncols; k++) {
                for (i = 0; i < nrows; i++) {
                    if (r[i] < sigmaVec[k]) {
                        U[i*ncols + k] = 0.0;
                        Up[i*ncols + k] = 0.0;
                    } else {
                        U[i*ncols + k] = - E[k] * exp(- z[k] * (r[i] - 1.0)) / r[i];
                        Up[i*ncols + k] = (1.0 + z[k]*r[i]) * U[i*ncols + k];        // -[d U(r)/dr] * r
                    }
                }
            }
          
            break;

        case 6: // ####---- Repulsive Yukaka ----####

            //*************************************************************
            //--------------------- REPULSIVE YUKAWA ----------------------
            //
            //                   U = T* exp[-λ (r-1)]/r                  
            //
            //*************************************************************

            E[0] = 1.0 / especie1.temperature;  //Temperature*
            E[2] = 1.0 / especie2.temperature;  //Temperature*
            E[1] = sqrt(E[0] * E[2]);

            z[0] = especie1.lambda;
            z[2] = especie2.lambda;
            z[1] = sqrt(z[0] * z[2]);


            for (k = 0; k < ncols; k++) {
                for (i = 0; i < nrows; i++) {
                    if (r[i] < sigmaVec[k]) {
                        U[i*ncols + k] = 0.0;
                        Up[i*ncols + k] = 0.0;
                    } else {
                        U[i*ncols + k] =  E[k] * exp( -z[k] * (r[i] - 1.0)) / r[i];
                        Up[i*ncols + k] = (1.0 + z[k]*r[i]) * U[i*ncols + k];        // -[d U(r)/dr] * r
                    }
                }
            }
          
            break;

        case 7: // ####---- Hard Sphere ----####

            //*************************************************************
            //----------------------- HARD SPHERE -------------------------
            //
            //                         U = 0.0                  
            //
            //*************************************************************

            for (k = 0; k < ncols; k++) {
                for (i = 0; i < nrows; i++) {
                    U[i*ncols + k] = 0.0;
                    Up[i*ncols + k] = 0.0;
                }
            }
          
            break;

        case 8: // ####---- Step Function ----####

            //*************************************************************
            //---------------------- STEP FUNCTION ------------------------
            //
            //                   U = T* λ   ∀ U ∈ [σ, T₂]                 
            //
            //*************************************************************

            E[0] = 1.0 / especie1.temperature;  //Temperature*
            E[2] = 1.0 / especie2.temperature;  //Temperature*
            E[1] = sqrt(E[0] * E[2]);

            z[0] = especie1.lambda;
            z[2] = especie2.lambda;
            z[1] = sqrt(z[0] * z[2]);

            E2[0] = especie1.temperature2;  //Width = T₂
            E2[2] = especie2.temperature2;  //Width = T₂
            E2[1] = sqrt(E2[0] * E2[2]);

            for (k = 0; k < ncols; k++) {
                for (i = 0; i < nrows; i++) {
                    if (r[i] < sigmaVec[k] || (r[i] > E2[k])) {
                        U[i*ncols + k] = 0.0;
                        Up[i*ncols + k] = 0.0;
                    } else {
                        U[i*ncols + k] =  E[k] * z[k];   // Step Potential with height = z[k]
                        Up[i*ncols + k] = 0.0;        // -[d U(r)/dr] * r
                    }
                }
            }
          
            break;

        case 9: // ####---- Down Hill Function ----####

            //*************************************************************
            //------------------- DOWN HILL FUNCTION ----------------------
            //
            //        U = T* (λ * (T₂ - r)/(T₂ - σ))   ∀ U ∈ [σ, T₂]                 
            //
            //*************************************************************

            // Please note that the inputs of this function are in kT/ε units (dimensionless)

            E[0] = 1.0 / especie1.temperature;  //Temperature*
            E[2] = 1.0 / especie2.temperature;  //Temperature*
            E[1] = sqrt(E[0] * E[2]);

            E2[0] = especie1.temperature2;  //Width = T₂
            E2[2] = especie2.temperature2;  //Width = T₂
            E2[1] = sqrt(E2[0] * E2[2]);

            z[0] = especie1.lambda / (E2[0] - sigmaVec[0]);  //slope
            z[2] = especie2.lambda / (E2[2] - sigmaVec[2]);  //slope
            z[1] = sqrt(z[0] * z[2]);
            
            for (k = 0; k < ncols; k++) {
                for (i = 0; i < nrows; i++) {
                    if (r[i] < sigmaVec[k] || (r[i] > E2[k])) {
                        U[i*ncols + k] = 0.0;
                        Up[i*ncols + k] = 0.0;
                    } else {
                        U[i*ncols + k] =  E[k] * (-z[k] * (r[i] - E2[k]));  // Down Hill Potential with height = λ
                        Up[i*ncols + k] = E[k] * z[k] * r[i];               // -[d U(r)/dr] * r
                    }
                }
            }
          
            break;

        case 10: // ####---- Gaussian Core Model ----####

            //*************************************************************
            //------------------- GAUSSIAN CORE MODEL ---------------------
            //
            //                   U = T* exp(- (r/σ)^2 )                 
            //
            //*************************************************************

            // Please note that the inputs of this function are in kT/ε units (dimensionless)

            E[0] = 1.0 / especie1.temperature;  //Temperature*
            E[2] = 1.0 / especie2.temperature;  //Temperature*
            E[1] = sqrt(E[0] * E[2]);

            for (k = 0; k < ncols; k++) {
                for (i = 0; i < nrows; i++) {
                    U[i*ncols + k] = E[k] * exp(- pow(r[i] / sigmaVec[k], 2.0));
                    Up[i*ncols + k] = 2.0 * pow(r[i] / sigmaVec[k], 2.0) * U[i*ncols + k];    // -[d U(r)/dr] * r
                }
            }
          
            break;
    }



    free(Ua);
    free(Ur);
    free(E);
    free(E2);
    free(z);
    free(z2);
}

void OZ2(double *Sk, double *Gr, int potentialID, int closureID, double alpha, double EZ, \
         double rmax, int nrho, char folderName[20], int *printFlag) {

    // Se modifica gammaOutput y cFuncMatrix
    // Ng(kj, gammaInput, gammaOutput, potentialID, closureID, cFuncMatrix, T, Tfin, alpha, EZ, rmax, nrho, printFlag)

    int i, k;
    int kj, IRY;
    double T, Tfin, pv, pv0, pv1, pv2;
    double chic, chic0, chic1, chic2;
    double ener, ener0, ener1, ener2;
    double rhoa, dT, drho, ddrho, dalpha, PexV;
    double *cFuncMatrix, *gammaInput1, *gammaInput2, *gammaOutput;
    
    // Allocate memory for an array of nrows doubles
    cFuncMatrix = malloc(nrows*ncols * sizeof(double));
    gammaInput1 = malloc(nrows*ncols * sizeof(double));
    gammaInput2 = malloc(nrows*ncols * sizeof(double));
    gammaOutput = malloc(nrows*ncols * sizeof(double));

    // Check if memory allocation was successful
    if (cFuncMatrix == NULL || gammaInput1 == NULL || gammaInput2 == NULL || gammaOutput == NULL) {
        printf("Memory allocation failed 2.\n");
        return; // Exit with an error code
    }

    Tfin = 0.0;

    for (k = 0; k < ncols; k++) {
        for (i = 0; i < nrows; i++) {
            gammaInput1[i*ncols + k] = 0.0;
//            cFuncMatrix[i*ncols + k] = 0.0;
        }
    }

    // ####################################### Anfang.
    rhoa = rho;
    dT = 1.0 / ((double) nrho);
    drho = rho / ((double) nrho);

    kj = 1;
    rho = kj * drho;
    T = dT * kj;

    Ng(kj, gammaInput1, gammaOutput, potentialID, closureID, cFuncMatrix, T, Tfin, alpha, EZ, rmax, nrho, printFlag);

    while (kj <= 1) {

        for (k = 0; k < ncols; k++) {
            for (i = 0; i < nrows; i++) {
                gammaInput1[i*ncols + k] = gammaOutput[i*ncols + k];
            }
        }

        kj++;
        rho = kj * drho;
        T = dT * kj;
        
        Ng(kj, gammaInput1, gammaOutput, potentialID, closureID, cFuncMatrix, T, Tfin, alpha, EZ, rmax, nrho, printFlag);

    }

    // Extrap modifica solo su primer argumento, en este caso gammaInput1
    Extrap(gammaInput1, gammaOutput, rho, drho);

    while(true) {

        for (k = 0; k < ncols; k++) {
            for (i = 0; i < nrows; i++) {
                gammaInput2[i*ncols + k] = gammaOutput[i*ncols + k];
            }
        }

        kj++;
        rho = kj * drho;
        T = dT * kj;

        Ng(kj, gammaInput1, gammaOutput, potentialID, closureID, cFuncMatrix, T, Tfin, alpha, EZ, rmax, nrho, printFlag);

        if (kj == nrho) {
            
            break;

        } else {
        
            // Extrap modifica solo su primer argumento, en este caso gammaInput2
            Extrap(gammaInput2, gammaOutput, rho, drho);

            for (k = 0; k < ncols; k++) {
                for (i = 0; i < nrows; i++) {
                    gammaInput1[i*ncols + k] = gammaInput2[i*ncols + k];
                }
            }
        }
    }

    // ----------------------------------- PY o HNC
    T = 1.0;
    rho = rhoa;

    Ng(kj, gammaInput1, gammaOutput, potentialID, closureID, cFuncMatrix, T, Tfin, alpha, EZ, rmax, nrho, printFlag);
    
    switch (closureID){

        case 1:

            printf("\n==============\n");
            printf("  Salida PY:\n");
            printf("==============\n\n");

            break;

        case 2:

            printf("\n==============\n");
            printf(" Salida HNC:\n");
            printf("==============\n\n");

            break;

        case 3:
            // ----------------------------------- RY

            printf("\n===========================\n");
            printf("  CALCULANDO VALOR ALPHA:\n");
            printf("===========================\n\n");

            ddrho = rhoa / 100.0;
            dalpha = -alpha / 50.0;

            do{
                T = 1.0;
                rho = rhoa;

                Ng(kj, gammaInput1, gammaOutput, potentialID, closureID, cFuncMatrix, T, Tfin, alpha, EZ, rmax, nrho, printFlag);
                
                for (k = 0; k < ncols; k++) {
                    for (i = 0; i < nrows; i++) {
                        gammaInput1[i*ncols + k] = gammaOutput[i*ncols + k];
                    }
                }

                Termo(gammaOutput, cFuncMatrix, &pv0, &chic0, &ener0);

                chic = chic0;

                // ------
                rho -= ddrho;

                Ng(kj, gammaInput1, gammaOutput, potentialID, closureID, cFuncMatrix, T, Tfin, alpha, EZ, rmax, nrho, printFlag);
                Termo(gammaOutput, cFuncMatrix, &pv1, &chic1, &ener1);

                // ------
                rho += 2.0*ddrho;

                Ng(kj, gammaInput1, gammaOutput, potentialID, closureID, cFuncMatrix, T, Tfin, alpha, EZ, rmax, nrho, printFlag);
                
                Termo(gammaOutput, cFuncMatrix, &pv2, &chic2, &ener2);

                RY(pv1, pv2, chic, ddrho, &alpha, dalpha, &IRY);

            } while(IRY == 1);
            
            // ------------------------
            printf("\n==============\n");
            printf("  Salida RY:\n");
            printf("==============\n\n");

            break;
    }

    // ----------------------------------- RY

    T = 1.0;
    Tfin = 1.0;
    
    printf(" Cálculo final: \n");

    rho = rhoa;

    printf("\n");
    Ng(kj, gammaInput1, gammaOutput, potentialID, closureID, cFuncMatrix, T, Tfin, alpha, EZ, rmax, nrho, printFlag);
    printf("\n");

    // --------------------------
    Escribe(gammaOutput, cFuncMatrix, Sk, Gr, potentialID, closureID, folderName);

    Termo(gammaOutput, cFuncMatrix, &pv, &chic, &ener);

    PexV = pv/rho - 1.0;

    // Write output
    printf("rho = %.10lf \n", (double) rho);
    printf("chic = %.10lf \n", (double) chic);
    printf("PexV = %.10lf \n", (double) PexV);
    printf("Eexe = %.10lf \n", (double) ener);

    free(cFuncMatrix);
    free(gammaInput1);
    free(gammaInput2);
    free(gammaOutput);

    return;
}


void Termo(double *gamma, double *cFuncMatrix, double *pv1, double *chic, double *ener) {
    
    int i, k;
    double ru1, ru2, ru3;
    double *r1;
    double *gMatrix;

    // Allocate memory for an array of nrows doubles
    r1 = malloc(nrows * sizeof(double));
    gMatrix = malloc(nrows*ncols * sizeof(double));

    // Check if memory allocation was successful
    if (r1 == NULL || gMatrix == NULL) {
        printf("Memory allocation failed 3.\n");
        return; // Exit with an error code
    }

    // dr = rmax / (1.0*nrows);

    for (i = 0; i < nrows; i++) {
        r1[i] = x[0]*x[0]*cFuncMatrix[i*ncols + 0] + 2.0*x[0]*x[1]*cFuncMatrix[i*ncols + 1] + x[1]*x[1]*cFuncMatrix[i*ncols + 2];
        r1[i] = r1[i] * r[i]*r[i];
    }

    *chic = 0.0;
    
    for (i = 1; i < nrows - 1; i++) {
        *chic += r1[i];
    }
    *chic = dr * (*chic + (r1[0] + r1[nrows-1]) / 2.0);
    *chic = 1.0 - 4.0 * M_PI * rho * (*chic);

    for (k = 0; k < ncols; k++) {
        for (i = 0; i < nrows; i++) {
            gMatrix[i + k*nrows] = gamma[i*ncols + k] + cFuncMatrix[i*ncols + k] + 1.0;
        }
    }

    for (i = 0; i < nrows; i++) {

        ru1 = x[0]*x[0] * gMatrix[i + 0*nrows] * Up[i*ncols + 0];
        ru2 = 2.0 * x[0]*x[1] * gMatrix[i + 1*nrows] * Up[i*ncols + 1];
        ru3 = x[1]*x[1] * gMatrix[i + 2*nrows] * Up[i*ncols + 2];
        r1[i] = (ru1 + ru2 + ru3) * r[i]*r[i];

    }

    *pv1 = 0.0;

    for (i = 1; i < (nrows-1); i++) {
        *pv1 += r1[i];
    }
    
    *pv1 = dr * ((*pv1) + (r1[0] + r1[nrows-1])/2.0);
    *pv1 = rho * (1.0 + 2.0*M_PI * rho*(*pv1)/3.0);

    for (i = 0; i < nrows; i++) {
        r1[i] = x[0]*x[0] * gMatrix[i + 0*nrows] * U[i*ncols + 0] + 2.0*x[0]*x[1] * gMatrix[i + 1*nrows]*U[i*ncols + 1];
        r1[i] = (r1[i] + x[1]*x[1] * gMatrix[i + 2*nrows] * U[i*ncols + 2]) * r[i]*r[i];
    }


    *ener = 0.0;

    for (i = 1; i < (nrows-1); i++) {
        *ener += r1[i];
    }

    *ener = dr * ((*ener) + (r1[0] + r1[nrows-1]) / 2.0);
    *ener *= 2.0 * M_PI * rho;

    free(r1);
    free(gMatrix);

    return;
}


void RY(double pv1, double pv2, double chic, double ddrho, double *alpha, double dalpha, int *IRY) {

    static double dif[2] = {0.0}; // Static array to maintain values across calls
    static int ix = 1; // Static index to alternate between 1 and 2
    double chiv, prod, A, B;

    *IRY = 0;
    chiv = (pv2 - pv1) / (2.0 * ddrho);
    dif[ix-1] = chic - chiv;

    prod = 0.0;

    if (ix < 2) {

        ix++;
        *alpha += dalpha;
        *IRY = 1;

//        printf("   CHIC = %.17g   CHIV = %.17g\n", chic, chiv);
//        printf("   DIFF = %.17g   \n", dif[ix - 1]);
//        printf("   ALPHA = %.17g  <------------------ \n\n", *alpha);
        printf("   ALPHA = %.10lf", (double) (*alpha));
        printf("   PROD = %.10lf \n", (double) prod);

        return;
    }

    prod = dif[0] * dif[1];

    if (prod < 0.0) {

        *alpha -= dalpha;
        A = (dif[1] - dif[0]) / dalpha;
        B = dif[0] - A * (*alpha);
        *alpha = -B / A;

        *IRY = 0;

//        printf("   CHIC = %.17g   CHIV = %.17g\n", chic, chiv);
//        printf("   DIFF = %.17g   \n", dif[ix - 1]);
//        printf("   ALPHA = %.17g  <------------------ \n\n", *alpha);
        printf("   ALPHA = %.10lf", (double) (*alpha));
        printf("   PROD = %.10lf \n", (double) prod);

        return;

    } else {

        dif[0] = dif[1];
        *alpha += dalpha;
        *IRY = 1;

//        printf("   CHIC = %.17g   CHIV = %.17g\n", chic, chiv);
//        printf("   DIFF = %.17g   \n", dif[ix - 1]);
//        printf("   ALPHA = %.17g  <------------------ \n\n", *alpha);
        printf("   ALPHA = %.10lf", (double) (*alpha));
        printf("   PROD = %.10lf \n", (double) prod);

        return;
    }
}


void Escribe(double *gamma, double *cFuncMatrix, double *Sk, double *Gr, int potentialID, int closureID, char folderName[20]) {

    int i, k;
    double dk, qmax, rk_max, sqmax, delta;

    double *rk, *c1, *gh, *gh2, *Ck, *S;

    // Allocate memory for an array of nrows doubles
    rk = malloc(nrows * sizeof(double));
    c1 = malloc(nrows*ncols * sizeof(double));
    gh = malloc(nrows*ncols * sizeof(double));
    gh2 = malloc(nrows*ncols * sizeof(double));
    Ck = malloc(nrows*ncols * sizeof(double));
    S = malloc(nrows*ncols * sizeof(double));

    // Check if memory allocation was successful
    if (rk == NULL || c1 == NULL || gh == NULL || gh2 == NULL || Ck == NULL || S == NULL) {
        printf("Memory allocation failed 4.\n");
        return; // Exit with an error code
    }

/*
    for (k = 0; k < ncols; k++) {
        
        g[0*ncols + k] = 0.0;
        
        for (i = 1; i < nrows; i++) {
            g[i*ncols + k] = gamma[i*ncols + k] + cFuncMatrix[i*ncols + k] + 1.0;
        }
    }
*/

    char fileGrName[120];
    char fileCrName[120];

    fileGrName[0] = '\0';
    fileCrName[0] = '\0';

    strcat(fileGrName, folderName);
    strcat(fileGrName, "/GdeR");
    appendclosureID(fileGrName, closureID);
    appendPotentialID(fileGrName, potentialID);
    strcat(fileGrName, ".dat");
    
    strcat(fileCrName, folderName);
    strcat(fileCrName, "/CdeR");
    appendclosureID(fileCrName, closureID);
    appendPotentialID(fileCrName, potentialID);
    strcat(fileCrName, ".dat");

    // salen las funciones de correlacion gij(r)
    FILE *file_Gr = fopen(fileGrName, "w");
    FILE *file_Cr = fopen(fileCrName, "w");

    for (i = 0; i < nrows; i++) {

        for (k = 0; k < ncols; k++) {
            gh[i*ncols + k] = gamma[i*ncols + k] + cFuncMatrix[i*ncols + k];
        }

        fprintf(file_Gr, "%.17lf\t%.17lf\t%.17lf\t%.17lf\n", r[i], gh[i*ncols + 0]+1.0, gh[i*ncols + 1]+1.0, gh[i*ncols + 2]+1.0);
        fprintf(file_Cr, "%.17lf\t%.17lf\t%.17lf\t%.17lf\n", r[i], \
                                                             cFuncMatrix[i*ncols + 0], cFuncMatrix[i*ncols + 1], cFuncMatrix[i*ncols + 2]);

        Gr[i*2 + 0] = r[i];
        Gr[i*2 + 1] = gh[i*ncols + 0] + 1.0;
    }

    fclose(file_Gr);
    fclose(file_Cr);

    // calculo de factores de estructura
    qmax = q[nrows - 1];
    rk_max = ((double) 0.5) * qmax;
    dk = rk_max / (1.0 * nrows);

    for (i = 0; i < nrows; i++) {
        rk[i] = (i+1) * dk;
    }

    // FT sólo modifica su segundo argumento
    FT(gh, gh2, rk, dr);
    FT(cFuncMatrix, c1, rk, dr);


    char file_hkij_Name[120];
    char file_ckij_Name[120];
    char file_ck_Name[120];

    file_hkij_Name[0] = '\0';
    file_ckij_Name[0] = '\0';
    file_ck_Name[0] = '\0';

    strcat(file_hkij_Name, folderName);
    strcat(file_hkij_Name, "/HdeKij");
    appendclosureID(file_hkij_Name, closureID);
    appendPotentialID(file_hkij_Name, potentialID);
    strcat(file_hkij_Name, ".dat");


    strcat(file_ckij_Name, folderName);
    strcat(file_ckij_Name, "/CdeKij");
    appendclosureID(file_ckij_Name, closureID);
    appendPotentialID(file_ckij_Name, potentialID);
    strcat(file_ckij_Name, ".dat");    
    
    
    strcat(file_ck_Name, folderName);
    strcat(file_ck_Name, "/CdeK11");
    appendclosureID(file_ck_Name, closureID);
    appendPotentialID(file_ck_Name, potentialID);
    strcat(file_ck_Name, ".dat");

    // salen las cij(k) y la Ck(k) para el cilindro
    FILE *file_hkij = fopen(file_hkij_Name, "w");
    FILE *file_ckij = fopen(file_ckij_Name, "w");
    FILE *file_ck = fopen(file_ck_Name, "w");

    for (i = 0; i < nrows; i++) {

        fprintf(file_hkij, "%.17lf\t%.17lf\t%.17lf\t%.17lf\n", rk[i], gh2[i*ncols + 0], gh2[i*ncols + 1], gh2[i*ncols + 2]);
        fprintf(file_ckij, "%.17lf\t%.17lf\t%.17lf\t%.17lf\n", rk[i], c1[i*ncols + 0], c1[i*ncols + 1], c1[i*ncols + 2]);
        fprintf(file_ck, "%.17lf\t%.17lf\n", rk[i], c1[i*ncols + 0]);

    }

    fclose(file_hkij);
    fclose(file_ckij);
    fclose(file_ck);

    for (i = 0; i < nrows; i++) {
        for (k = 0; k < ncols; k++) {
            Ck[i*ncols + k] = c1[i*ncols + k];
        }
    }

    for (i = 0; i < nrows; i++) {
        delta = (1.0 - rho*x[0] * Ck[i*ncols + 0]) * (1.0 - rho*x[1] * Ck[i*ncols + 2]);
        delta -= pow(rho, 2.0) * x[0]*x[1] * pow(Ck[i*ncols + 1], 2.0);
        
        S[i*ncols + 0] = (1.0 - rho*x[1] * Ck[i*ncols + 2]) * Ck[i*ncols + 0];
        S[i*ncols + 0] = (S[i*ncols + 0] + rho*x[1] * pow(Ck[i*ncols + 1], 2.0)) / delta;
        S[i*ncols + 1] = Ck[i*ncols + 1] / delta;
        S[i*ncols + 2] = (1.0 - rho*x[0] * Ck[i*ncols + 0]) * Ck[i*ncols + 2];
        S[i*ncols + 2] = (S[i*ncols + 2] + rho*x[0] * pow(Ck[i*ncols + 1], 2.0)) / delta;
        S[i*ncols + 0] = x[0] + rho*pow(x[0], 2.0) * S[i*ncols + 0];
        S[i*ncols + 1] = rho * x[0]*x[1] * S[i*ncols + 1];
        S[i*ncols + 2] = x[1] + rho*pow(x[1], 2.0) * S[i*ncols + 2];
    }

    printf("-----------------------------\n");
    printf("Salida para Sij(k): Sqij.dat \n");
    printf("--> orden: S11, S12, S22\n");
    printf("-----------------------------\n");

    char fileSqName[120];
    
    fileSqName[0] = '\0';
    
    strcat(fileSqName, folderName);
    strcat(fileSqName, "/Sq");
    appendclosureID(fileSqName, closureID);
    appendPotentialID(fileSqName, potentialID);
    strcat(fileSqName, ".dat");
    
    FILE *file_Sq = fopen(fileSqName, "w");

    sqmax = 0.0;

    for (i = 0; i < nrows; i++) {
        
        fprintf(file_Sq, "%.17lf\t%.17lf\t%.17lf\t%.17lf\n", rk[i], S[i*ncols + 0]/x[0], \
                                                             S[i*ncols + 1]/sqrt(x[0]*x[1]), S[i*ncols + 2]/x[1]);
        
        Sk[i*2 + 0] = rk[i];
        Sk[i*2 + 1] = S[i*ncols + 0];
        
        if (S[i*ncols + 0] > sqmax) {
            sqmax = S[i*ncols + 0];
        }

    }
    fclose(file_Sq);

    printf("\n");
    printf("sqmax = %.12lf\n", sqmax);
    printf("\n");

    free(rk);
    free(c1);
    free(gh);
    free(gh2);
    free(Ck);
    free(S);
}

void Ng(int kj, double *gammaInput, double *gammaOutput, int potentialID, int closureID, double *cFuncMatrix, \
        double T, double Tfin, double alpha, double EZ, double rmax, int nrho, int *printFlag) {

    // Se modifica gammaOutput y cFuncMatrix
    // ONg(gammaInput, gammaOutput, potentialID, closureID, cFuncMatrix, T, Tfin, alpha, rmax)

    int i, k, flag;

    double ETA, V, condition1;
    double *f;
    double *g1, *g2, *g3;
    double *d1, *d2, *d3, *d01, *d02;
    double *d01d01, *d01d02, *d02d02, *d3d01, *d3d02;
    double *const1, *const2;

    // Allocate memory for an array of nrows doubles
    f = malloc(nrows*ncols * sizeof(double));
    g1 = malloc(nrows*ncols * sizeof(double));
    g2 = malloc(nrows*ncols * sizeof(double));
    g3 = malloc(nrows*ncols * sizeof(double));
    d1 = malloc(nrows*ncols * sizeof(double));
    d2 = malloc(nrows*ncols * sizeof(double));
    d3 = malloc(nrows*ncols * sizeof(double));
    d01 = malloc(nrows*ncols * sizeof(double));
    d02 = malloc(nrows*ncols * sizeof(double));
    d01d01 = malloc(ncols * sizeof(double));
    d01d02 = malloc(ncols * sizeof(double));
    d02d02 = malloc(ncols * sizeof(double));
    d3d01 = malloc(ncols * sizeof(double));
    d3d02 = malloc(ncols * sizeof(double));
    const1 = malloc(ncols * sizeof(double));
    const2 = malloc(ncols * sizeof(double));

    // Check if memory allocation was successful
    if (f == NULL || g1 == NULL || g2 == NULL || g3 == NULL || d1 == NULL || d2 == NULL || d3 == NULL) {
        printf("Memory allocation failed 5.\n");
        return; // Exit with an error code
    }

    if (d01 == NULL || d02 == NULL || d01d01 == NULL || d01d02 == NULL || d02d02 == NULL) {
        printf("Memory allocation failed 6.\n");
        return; // Exit with an error code
    }

    if (d3d01 == NULL || d3d02 == NULL || const1 == NULL || const2 == NULL ) {
        printf("Memory allocation failed 7.\n");
        return; // Exit with an error code
    }


    for (k = 0; k < ncols; k++) {
        for (i = 0; i < nrows; i++) {
            f[i*ncols + k] = gammaInput[i*ncols + k];
        }
    }

    ONg(f, g1, potentialID, closureID, cFuncMatrix, T, Tfin, alpha, rmax);
    ONg(g1, g2, potentialID, closureID, cFuncMatrix, T, Tfin, alpha, rmax);
    ONg(g2, g3, potentialID, closureID, cFuncMatrix, T, Tfin, alpha, rmax);

    // Calculate d1 and d2 arrays
    for (k = 0; k < ncols; k++) {
        for (i = 0; i < nrows; i++) {
            d1[i*ncols + k] = (g1[i*ncols + k] - f[i*ncols + k]);
            d2[i*ncols + k] = (g2[i*ncols + k] - g1[i*ncols + k]);
        }
    }

/*
    while (kj >= 2) {
        
        // Calculate d3 array
        for (k = 0; k < ncols; k++) {
            for (i = 0; i < nrows; i++) {
                d3[i*ncols + k] = (g3[i*ncols + k] - g2[i*ncols + k]);
            }
        }

        // Calculate D01 and d02 arrays
        for (k = 0; k < ncols; k++) {
            for (i = 0; i < nrows; i++) {
                d01[i*ncols + k] = (d3[i*ncols + k] - d2[i*ncols + k]);
                d02[i*ncols + k] = (d3[i*ncols + k] - d1[i*ncols + k]);
            }
        }

        pp(dr, d01, d01, d01d01);
        pp(dr, d01, d02, d01d02);
        pp(dr, d3, d01, d3d01);
        pp(dr, d02, d02, d02d02);
        pp(dr, d3, d02, d3d02);

        V = (double) 1.0E-50;

        for (k = 0; k < ncols; k++) {

            if ((d01d01[k] <= V) || (d01d02[k] <= V)) {
                flag = 1;
                break;
            } else if((d02d02[k] - (d01d02[k]*d01d02[k]) / d01d01[k]) <= V) {
                flag = 1;
                break;
            } else {
                const2[k] = d3d02[k] - d01d02[k]*d3d01[k] / d01d01[k];
                const2[k] /= (d02d02[k] - d01d02[k]*d01d02[k] / d01d01[k]);
                const1[k] = (d3d02[k] - d02d02[k] * const2[k]) / d01d02[k];                
            
                for (i = 0; i < nrows; i++) {
                    f[i*ncols + k] = (1.0 - const1[k]-const2[k]) * g3[i*ncols + k];
                    f[i*ncols + k] += const1[k]*g2[i*ncols + k] + const2[k]*g1[i*ncols + k];
                }
                
                flag = 0; 
            }
        }


        if (flag == 0){
            ONg(f, g3, potentialID, closureID, cFuncMatrix, T, Tfin, alpha, rmax);

            // Calculate d3 array
            for (k = 0; k < ncols; k++) {
                for (i = 0; i < nrows; i++) {
                    d3[i*ncols + k] = (g3[i*ncols + k] - f[i*ncols + k]);
                }
            }
        }

        // La función Pres modifica la matriz d3 y el parámetro ETA
        Pres(d3, dr, &ETA);

        // Check ETA condition
        if (ETA <= EZ) {

            break;

        } else {

            for (k = 0; k < ncols; k++) {
                for (i = 0; i < nrows; i++) {
                    g1[i*ncols + k] = g2[i*ncols + k];
                    d1[i*ncols + k] = d2[i*ncols + k];
                    g2[i*ncols + k] = g3[i*ncols + k];
                    d2[i*ncols + k] = d3[i*ncols + k];
                }
            }

            ONg(g2, g3, potentialID, closureID, cFuncMatrix, T, Tfin, alpha, rmax);
        }
    }
*/

    while (kj >= 2) {
        
        // Calculate d3 array
        for (k = 0; k < ncols; k++) {
            for (i = 0; i < nrows; i++) {
                d3[i*ncols + k] = (g3[i*ncols + k] - g2[i*ncols + k]);
            }
        }

        // Calculate D01 and d02 arrays
        for (k = 0; k < ncols; k++) {
            for (i = 0; i < nrows; i++) {
                d01[i*ncols + k] = (d3[i*ncols + k] - d2[i*ncols + k]);
                d02[i*ncols + k] = (d3[i*ncols + k] - d1[i*ncols + k]);
            }
        }

        pp(dr, d01, d01, d01d01);
        pp(dr, d01, d02, d01d02);
        pp(dr, d3, d01, d3d01);
        pp(dr, d02, d02, d02d02);
        pp(dr, d3, d02, d3d02);

        V = (double) 1.0E-50;

        for (k = 0; k < ncols; k++) {

            condition1 = d02d02[k] - (d01d02[k]*d01d02[k]) / d01d01[k];

            if ( ((d01d01[k] > V) && (d01d02[k] > V)) && (condition1 > V) ) {
                    
                const2[k] = d3d02[k] - d01d02[k]*d3d01[k] / d01d01[k];
                const2[k] = const2[k] / (d02d02[k] - d01d02[k]*d01d02[k] / d01d01[k]);
                const1[k] = (d3d02[k] - d02d02[k] * const2[k]) / d01d02[k];
                
                for (i = 0; i < nrows; i++) {
                    f[i*ncols + k] = (1.0 - const1[k] - const2[k]) * g3[i*ncols + k];
                    f[i*ncols + k] = f[i*ncols + k] + const1[k]*g2[i*ncols + k] + const2[k]*g1[i*ncols + k];
                }

                flag = 0;

            } else {          
                    
                flag = 1;
                break;

            }
        }

        if (flag == 0){
            
            ONg(f, g3, potentialID, closureID, cFuncMatrix, T, Tfin, alpha, rmax);

            // Calculate d3 array
            for (k = 0; k < ncols; k++) {
                for (i = 0; i < nrows; i++) {
                    d3[i*ncols + k] = g3[i*ncols + k] - f[i*ncols + k];
                }
            }

        }

        // La función Pres modifica la matriz d3 y el parámetro ETA
        Pres(d3, dr, &ETA);

        // Check ETA condition
        if (ETA <= EZ) {
            break;
        }

        for (k = 0; k < ncols; k++) {
            for (i = 0; i < nrows; i++) {
                g1[i*ncols + k] = g2[i*ncols + k];
                d1[i*ncols + k] = d2[i*ncols + k];
                g2[i*ncols + k] = g3[i*ncols + k];
                d2[i*ncols + k] = d3[i*ncols + k];
            }
        }

        ONg(g2, g3, potentialID, closureID, cFuncMatrix, T, Tfin, alpha, rmax);
    }

    // Guardamos los resultados en gammaOutput
    for (k = 0; k < ncols; k++) {
        for (i = 0; i < nrows; i++) {
            gammaOutput[i*ncols + k] = g3[i*ncols + k];
        }
    }

    if ((*printFlag) == 0){
        printf("ITERACION %d \t %.17e \t %d\n", kj, (float) ETA, nrows);
    }

    if (kj >= nrho){
        *printFlag = 1;
    }
    
    free(f);
    free(g1);
    free(g2);
    free(g3);
    free(d1);
    free(d2);
    free(d3);
    free(d01);
    free(d02);
    free(d01d01);
    free(d01d02);
    free(d02d02);
    free(d3d01);
    free(d3d02);
    free(const1);
    free(const2);
}

// ##############################################################################################################

void closrel(double *gamma, int potentialID, int closureID, double *cFuncMatrix, double T, double alpha) {
    
    int i, k;
    double sigmaAux[ncols];
    double arg, F;

    for (k = 0; k < ncols; k++) {
        if (potentialID == 1 || potentialID == 2 || potentialID == 3) {
            sigmaAux[k] = (sigmaVec[k] / 2.0);
        } else if (potentialID == 10) {
            sigmaAux[k] = 0.0;
        } else {
            sigmaAux[k] = sigmaVec[k];
        }
    }


    for (k = 0; k < ncols; k++) {
        for (i = 0; i < nrows; i++) {
            cFuncMatrix[i*ncols + k] = gamma[i*ncols + k] + 1.0;
        }
    }

    switch(closureID){
        
        case 1: // ####---- PY ----####

            for (k = 0; k < ncols; k++) {

                for (i = 0; i < nrows; i++) {

                    if (r[i] < sigmaAux[k]) {
                        cFuncMatrix[i*ncols + k] = -cFuncMatrix[i*ncols + k];
                    
                    } else {
                        arg = U[i*ncols + k] * T;
                    
                        if (arg > 70.0) {
                            cFuncMatrix[i*ncols + k] = -cFuncMatrix[i*ncols + k];
                    
                        } else {
                            cFuncMatrix[i*ncols + k] = -(exp(-arg) - 1.0) * cFuncMatrix[i*ncols + k];
                    
                        }
                    }
                }
            }
            
            return;
        
        case 2: // ####---- HNC ----####

            for (k = 0; k < ncols; k++) {

                for (i = 0; i < nrows; i++) {

                    if (r[i] < sigmaAux[k]) {

                        cFuncMatrix[i*ncols + k] = -cFuncMatrix[i*ncols + k];

                    } else {

                        arg = U[i*ncols + k] * T - gamma[i*ncols + k];

                        if (arg > 70.0) {

                            cFuncMatrix[i*ncols + k] = -cFuncMatrix[i*ncols + k];

                        } else {

                            cFuncMatrix[i*ncols + k] = exp(-arg) - cFuncMatrix[i*ncols + k];

                        }
                    }

                }
            }

            return;
        
        case 3: // ####---- RY ----####

            for (k = 0; k < ncols; k++) {

                for (i = 0; i < nrows; i++) {

                    if (r[i] < sigmaAux[k]) {

                        cFuncMatrix[i*ncols + k] = -cFuncMatrix[i*ncols + k];

                    } else {

                        F = 1.0 - exp(-alpha * r[i]);

                        arg = (exp(gamma[i*ncols + k] * F) - 1.0) / F;

                        cFuncMatrix[i*ncols + k] = exp(-U[i*ncols + k] * T) * (1.0 + arg) - cFuncMatrix[i*ncols + k];

                    }
                }
            }

            return;
    }

}

void appendclosureID(char *inputString, int closureID) {
    
    switch (closureID){
    
        case 1:
            strcat(inputString, "_PY");
            break;
        
        case 2:
            strcat(inputString, "_HNC");
            break;
        
        case 3:
            strcat(inputString, "_RY");
            break;
    }
}

void appendPotentialID(char *inputString, int potentialID) {
    
    switch (potentialID){
    
        case 1:
            strcat(inputString, "_IPL");
            break;
        
        case 2:
            strcat(inputString, "_LJT");
            break;
        
        case 3:
            strcat(inputString, "_LJT2");
            break;

        case 4:
            strcat(inputString, "_DBLEYUK");
            break;

        case 5:
            strcat(inputString, "_ATRYUK");
            break;

        case 6:
            strcat(inputString, "_REPYUK");
            break;

        case 7:
            strcat(inputString, "_HS");
            break;

        case 8:
            strcat(inputString, "_STEPFUNC");
            break;

        case 9:
            strcat(inputString, "_DOWNHILL");
            break;
    }
}