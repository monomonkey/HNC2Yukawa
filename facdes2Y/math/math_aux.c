#include "math_aux.h"

void pp(double dr, double *matrix1, double *matrix2, double *prod) {

    int k, i;
    double average;

    for (k = 0; k < ncols; k++) {

        prod[k] = 0.0;
        
        for (i = 1; i < (nrows-1); i++) {
        
            prod[k] += matrix1[i*ncols + k] * matrix2[i*ncols + k];
        
        }
        
        average = (matrix1[0*ncols + k] * matrix2[0*ncols + k]);
        average += (matrix1[(nrows-1)*ncols + k] * matrix2[(nrows-1)*ncols + k]);
        average *= 0.5;

        prod[k] = (prod[k] + average) * dr;
    }

}

void ONg(double *gammaInput, double *gammaOutput, int potentialID, int closureID, double *cFuncMatrix, \
         double T, double Tfin, double alpha, double rmax) {
    
    int i, k;
    double sqmax, delta;
    double *S, *Ck;

    Ck = malloc(nrows*ncols * sizeof(double));
    S = malloc(nrows*ncols * sizeof(double));

    if (S == NULL || Ck == NULL) {
        printf("Memory allocation failed 6.\n");
        return; // Exit with an error code
    }

    // closrel modifica la matriz cFuncMatrix que se declara en OZ2
    closrel(gammaInput, potentialID, closureID, cFuncMatrix, T, alpha);
    
    // Se calcula la transformada seno de cada columna de cFuncMatrix
    // El 1 como último argumento en FFTM indica que es la transformada normal (no inversa)
    // NOTA: se sobreescriben los datos de la matriz cFuncMatrix
    FFTM(cFuncMatrix, rmax, 1);

    // Almacenamos los resultados de cFuncMatrix en Ck, ya que volveremos a sobreescribir las
    // entradas de cFuncMatrix más adelante.
    for (i = 0; i < nrows; i++) {
        for (k = 0; k < ncols; k++) {
            Ck[i*ncols + k] = cFuncMatrix[i*ncols + k];
        }
    }
/*
    for (i=0; i<nrows; i++){
        printf("%d\t%.15e\t%.15e\t%.15e\n", i, Ck[i*ncols + 0], Ck[i*ncols + 1], Ck[i*ncols + 2]);
    }
*/

    // NOTA: Aunque hallamos puesto la variable ncols como global
    // la siguiente estructura sólo contempla el valor de ncols = 3
    for (i = 0; i < nrows; i++) {
        delta = (1.0 - rho*x[0]*Ck[i*ncols + 0]) * \
                (1.0 - rho*x[1]*Ck[i*ncols + 2]);
        delta -= pow(rho, 2.0)*x[0]*x[1] * \
                 pow(Ck[i*ncols + 1], 2.0);

        gammaOutput[i*ncols + 0] = (1.0 - rho*x[1] * Ck[i*ncols + 2]) * Ck[i*ncols + 0];
        gammaOutput[i*ncols + 0] = (gammaOutput[i*ncols + 0] + \
                                    rho*x[1] * pow(Ck[i*ncols + 1], 2.0)) / \
                                    delta;
        gammaOutput[i*ncols + 0] = gammaOutput[i*ncols + 0] - Ck[i*ncols + 0];
        gammaOutput[i*ncols + 1] = Ck[i*ncols + 1] / delta - Ck[i*ncols + 1];
        gammaOutput[i*ncols + 2] = (1.0 - rho*x[0] * Ck[i*ncols + 0]) * Ck[i*ncols + 2];
        gammaOutput[i*ncols + 2] = (gammaOutput[i*ncols + 2] + rho*x[0] * pow(Ck[i*ncols + 1], 2.0)) / delta;
        gammaOutput[i*ncols + 2] = gammaOutput[i*ncols + 2] - Ck[i*ncols + 2];
    }

    // Se calcula la transformada seno INVERSA de cada columna de gammaOutput
    // El -1 como último argumento en FFTM indica que es la transformada inversa
    // NOTA: se sobreescriben los datos de la matriz cFuncMatrix

    FFTM(gammaOutput, rmax, -1);
/*
    for (i=0; i<nrows; i++){
        printf("%d\t%.15e\n", i, gammaOutput[i*ncols + 0]);
    }
*/
    closrel(gammaOutput, potentialID, closureID, cFuncMatrix, T, alpha);
/*
    for (i=0; i<nrows; i++){
        printf("%d\t%.15e\t%.15e\t%.15e\n", i, cFuncMatrix[i*ncols + 0], cFuncMatrix[i*ncols + 1], cFuncMatrix[i*ncols + 2]);
    }
*/
/*
    for (i=0; i<nrows; i++){
        printf("%d\t%.15e\t%.15e\t%.15e\n", i, gammaOutput[i*ncols + 0], gammaOutput[i*ncols + 1], gammaOutput[i*ncols + 2]);
    }
*/

    if (Tfin < 1.0) {
        free(Ck);
        free(S);

        return;
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

    sqmax = 0.0;

//    FILE *outFile = fopen("Sq2.out", "w");
    for (i = 0; i < nrows; i++) {

//        fprintf(outFile, "%.17e \t %.17e \t %.17e \t %.17e \n", q[i], S[i*ncols + 0], S[i*ncols + 1], S[i*ncols + 2]);
        
        if (S[i*ncols + 0] > sqmax) {
            sqmax = S[i*ncols + 0];
        }
    }
//    fclose(outFile);

    printf(" sqmax = %.17lf \n", sqmax);

    free(Ck);
    free(S);
}

// ##############################################################################################################

void Extrap(double *gammaOut, double *gammaIn, double rho, double drho) {

    int i, k;
    double a, b;

    for (k = 0; k < ncols; k++) {
        for (i = 0; i < nrows; i++) {
            a = (gammaIn[i*ncols + k] - gammaOut[i*ncols + k]) / drho;
            b = gammaIn[i*ncols + k] - a*rho;
            gammaOut[i*ncols + k] = a*(rho + drho) + b;
        }
    }
}

void interp(int m, int n, double *r1, double *gammaInput, double *gammaOutput) {
    
    int i, j, k;
    double a[nrows], b[nrows];

    for (k = 0; k < ncols; k++) {
        for (i = 0; i < (m-1); i++) {
            a[i] = (gammaInput[(i+1)*ncols + k] - gammaInput[i*ncols + k]) / (r1[i+1] - r1[i]);
            b[i] = gammaInput[(i+1)*ncols + k] - a[i] * r1[i+1];
        }
        for (j = 0; j < (m-1); j++) {
            for (i = ((j*n)/m); i < (((j+1)*n) / m); i++) {
                gammaOutput[i*ncols + k] = a[j] * r[i] + b[j];
            }
        }
        for (i = n-(n/m); i < n; i++) {
            gammaOutput[i*ncols + k] = gammaOutput[(n-(n/m))*ncols + k];
        }
    }
}


void Pres(double *f, double dr, double *eta) {

    int i, k;
    double sum, average;

    *eta = 0.0;

    for (k = 0; k < ncols; k++) {
        
        sum = 0.0;
        
        for (i = 1; i < (nrows-1); i++) {

            sum += pow(f[i*ncols + k], 2.0);

        }
        average = pow(f[0*ncols + k], 2.0) + \
                  pow(f[(nrows-1)*ncols + k], 2.0);
        average *= 0.5;

        sum += average;
        
        *eta += sum * dr;
    }

    *eta = sqrt(*eta);
}


void FFTM(double *inputDataMatrix, double rmax, int isDirect) {

    /* Este función se llama FFTM por Fast Fourier Transfor Matrix
       y calcula la transformada de fourier de una matriz, pero lo hace
       columna por columna.

       inputDataMatrix :> La matriz a la que se le calcula la FFT
       n  :> El número de filas de la matriz fa
       isDirect :> Parámetro que indica si se calcula la transformada directa o la inversa
         Si isDirect  = 1 se realiza la transformada seno de inputData
         Si isDirect != 1 se realiza la transformada inversa seno de inputData

       OJO: Esta función YA es aplicable a cualquier matriz aunque
            manualmente se introduce el número de columnas, que es ncols.
    */

    int i, k;
    double *tempVector;  //JJ:Este es un vector temporal para realizar la FFT

    // Allocate memory for an array of nrows doubles
    tempVector = malloc(nrows * sizeof(double));

    // Check if memory allocation was successful
    if (tempVector == NULL) {
        printf("Memory allocation failed 8.\n");
        return; // Exit with an error code
    }

    for (i = 0; i < ncols; i++) {
        //JJ:Guarda cada columna de inputDataMatrix en tempVector
        for (k = 0; k < nrows; k++) {
            tempVector[k] = inputDataMatrix[k*ncols + i];
        }

        //JJ:Le calcula la transformada de FFT a tempVector
        FFT(tempVector, rmax, isDirect);

/*
        for (k=0; k<nrows; k++){
            printf("%d\t%.15e\n", k, tempVector[k]);
        }
*/

        //JJ:El resultado lo guarda en inputDataMatrix de regreso
        for (k = 0; k < nrows; k++) {
            inputDataMatrix[k*ncols + i] = tempVector[k];
        }
    }

    free(tempVector);
}


double calint(double *f, double dr) {
    
    double rint = 0.0;
    
    for (int i = 1; i < (nrows-1); i++) {
        rint += f[i];
    }

    rint = dr * (rint + (f[0] + f[nrows-1]) / 2.0);
    
    return rint;
}


void intt(double *h, double dr, double *sft) {
    
    double *f;

    // Allocate memory for an array of nrows doubles
    f = malloc(nrows * sizeof(double));

    // Check if memory allocation was successful
    if (f == NULL) {
        printf("Memory allocation failed 9.\n");
        return; // Exit with an error code
    }
    
    for (int k = 0; k < ncols; k++) {
        for (int i = 0; i < nrows; i++) {

            f[i] = h[i*ncols + k];
        }

        sft[k] = calint(f, dr);
    }

    free(f);
}


void FT(double *c, double *c1, double *rk, double dr) {

    int i, j, k;
    double *ca, *sft;

    // Allocate memory for an array of nrows doubles
    ca = malloc(nrows*ncols * sizeof(double));
    sft = malloc(ncols * sizeof(double));

    // Check if memory allocation was successful
    if (ca == NULL || sft == NULL) {
        printf("Memory allocation failed 10.\n");
        return; // Exit with an error code
    }
    
    for (i = 0; i < nrows; i++) {
        
        for (k = 0; k < ncols; k++) {
            for (j = 0; j < nrows; j++) {

                ca[j*ncols + k] = r[j] * c[j*ncols + k] * sin(rk[i] * r[j]);
            
            }
        }

        intt(ca, dr, sft);

        for (int k = 0; k < ncols; k++) {
            c1[i*ncols + k] = 4.0 * M_PI * sft[k] / rk[i];
        }
    }

    free(ca);
    free(sft);
}


void FFT(double *inputData, double rmax, int isDirect) {
    
    /*
    Si isDirect =  1 se realiza la transformada seno de inputData
    Si isDirect != 1 se realiza la transformada inversa seno de inputData
    */

    int i;
    double a;

    for (i = 0; i < nrows; i++) {
        inputData[i] = (i * inputData[i]);
    }

    // Se ejecuta la transformada seno
    // Note que aquí no es necesario usar isDirect
    // porque la transformada seno es su propia inversa
    sinft(inputData, nrows); 

    if (isDirect == 1) {
//        printf("Se calcula la transformada seno directa");
        a = 4.0 * pow(rmax, 3.0) / (1.0 * nrows*nrows);
        inputData[0] = 0.0;
    }else{
//        printf("Se calcula la transformada seno inversa");
        a = nrows * (1.0 / (2.0 * pow(rmax, 3.0)));
    }
    
    for (i = 1; i < nrows; i++) {
        inputData[i] = a * inputData[i] / (1.0 * i);
//        printf("%d\t%.15e\n", i, inputData[i]);
    }
}


void sinft(double *y, int nmax) {
    
    int m, j, i;
    double wr, wi, wpr, wpi, wtemp, theta;
    double y1, y2, sum;

    theta = M_PI / nrows;
    wr = 1.0;
    wi = 0.0;
    wpr = -2.0 * pow(sin(theta/2.0), 2.0);
    wpi = sin(theta);
    y[0] = 0.0;
    m = nrows / 2;

    for (j = 1; j <= m; j++) {
        wtemp = wr;
        wr = wr*wpr - wi*wpi + wr;
        wi = wi * wpr + wtemp * wpi + wi;
        y1 = wi * (y[j] + y[nrows-j]);
        y2 = (y[j] - y[nrows-j]) / 2.0;
//        printf("%d\t%.15E\t%.15E\t%.15E\t%.15E\n", j, wi, (y[j] + y[nrows-j]), wi * (y[j] + y[nrows-j]), y1);
        y[j] = y1 + y2;
        y[nrows-j] = y1 - y2;
    }

    realft(y, m, 1, nmax);
/*
    for (i=0; i<nrows; i++){
        printf("%d\t%.15e\n", i, y[i]);
    }
*/
    sum = 0.0;
    y[0] = y[0] / 2.0;
    y[1] = 0.0;
    
    for (j = 0; j < (nrows-1); j+=2) {
        sum += y[j];
        y[j] = y[j+1];
        y[j+1] = sum;
//        printf("%.15E\n", sum);
    }

    return;
}


void realft(double *data, int n, int isign, int nmax) {
    
    int i, n2p3, i1, i2, i3, i4;
    double wr, wi, wpr, wpi, wtemp, theta;
    double c1, c2, wrs, wis, h1r, h1i, h2r, h2i;

    theta = M_PI / ((double) n);
    
    c1 = 0.5;

    if (isign == 1) {
        c2 = -0.5;

        four1(data, n, +1, nmax); 

/*
        for (i=0; i<n; i++){
            printf("%d\t%.15e\n", i, data[i]);
        }
*/


    } else {
        c2 = 0.5;
        theta = -theta;
    }

    wpr = -2.0 * pow(sin(theta/2.0), 2.0);
    wpi = sin(theta);
    wr = 1.0 + wpr;
    wi = wpi;
    n2p3 = 2*n + 3;

    for (i = 2; i <= n / 2; i++) {
        i1 = 2*i - 1;
        i2 = i1 + 1;
        i3 = n2p3 - i2;
        i4 = i3 + 1;
        wrs = ((float) wr);
        wis = ((float) wi);
        h1r = c1 * (data[i1-1] + data[i3-1]);
        h1i = c1 * (data[i2-1] - data[i4-1]);
        h2r = -c2 * (data[i2-1] + data[i4-1]);
        h2i = c2 * (data[i1-1] - data[i3-1]);
        data[i1-1] = h1r + wrs * h2r - wis * h2i;
        data[i2-1] = h1i + wrs * h2i + wis * h2r;
        data[i3-1] = h1r - wrs * h2r + wis * h2i;
        data[i4-1] = -h1i + wrs * h2i + wis * h2r;
        wtemp = wr;
        wr = wr * wpr - wi * wpi + wr;
        wi = wi * wpr + wtemp * wpi + wi;
    }

    if (isign == 1) {
        h1r = data[0];
        data[0] = h1r + data[1];
        data[1] = h1r - data[1];
    } else {
        h1r = data[0];
        data[0] = c1 * (h1r + data[1]);
        data[1] = c1 * (h1r - data[1]);
        four1(data, n, -1, nmax);
    }
}


void four1(double *data, int nn, int isign, int nmax) {

    int i, j, istep;
    int m, n, mmax;
    double wr, wi, wpr, wpi, wtemp, theta;
    float tempr, tempi;
/*
    for (i=0; i<nn; i++){
        printf("A %d\t%.15e\n", i, data[i]);
    }
*/

    n = 2 * nn;
    j = 1;

    for (i = 1; i <= n; i += 2) {

        if (j > i) {
            tempr = data[j-1];
            tempi = data[j];
            data[j-1] = data[i-1];
            data[j] = data[i];
            data[i-1] = tempr;
            data[i] = tempi;
        }

        m = n / 2;

        while ((m >= 2) && (j > m)) {
            j -= m;
            m /= 2;
        }

        j += m;
    }

/*    for (i=0; i<nn; i++){
        printf("A %d\t%.15e\n", i, data[i]);
    }
*/
    mmax = 2;

    while (n > mmax) {

        istep = 2 * mmax;
        theta = 2.0 * M_PI / (isign*mmax);
        wpr = - 2.0 * pow(sin(theta / 2.0), 2.0);
        wpi = sin(theta);
        wr = 1.0;
        wi = 0.0;

        for (m = 1; m <= mmax; m += 2) {
            for (i = m; i <= n; i += istep) {
                j = i + mmax;
                tempr = ((float) wr) * data[j-1] - ((float) wi) * data[j];
                tempi = ((float) wr) * data[j] + ((float) wi) * data[j-1];

//                printf("%d %d\t%.15lf\t%.15lf\n", m, i, (float) wr, (float) wi);
//                printf("%d %d\t%.15lf\t%.15lf\n", m, i, tempr, tempi);
//                printf("%d %d\t%.15lf\t%.15lf\n", m, i, data[i-1], data[i]);

                data[j-1] = data[i-1] - ((float) tempr);
                data[j] = data[i] - ((float) tempi);
                data[i-1] += ((float) tempr);
                data[i] += ((float) tempi);
//                printf("%d %d\t%.15lf\t%.15lf\n", m, i, data[i-1], data[i]);
            }
            wtemp = wr;
            wr = wr * wpr - wi * wpi + wr;
            wi = wi * wpr + wtemp * wpi + wi;
        }

        mmax = istep;

    }
/*
    for (i=0; i<nn; i++){
        printf("A %d\t%.15e\n", i, data[i]);
    }
*/
    return;
}
