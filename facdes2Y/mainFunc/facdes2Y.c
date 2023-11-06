#include "facdes2Y.h"

// We define some global parameters
int nrows, ncols;
double rho, dr;
double *r, *q, x[2];
double *U, *Up, *sigmaVec;

int facdes2YFunc(const int nodes, int nrho, double rmax, int potentialID, int closureID, species especie1, species especie2, \
                 double volumeFactor, double alpha, double xnu, double EZ) {

    nrows = nodes;
    ncols = 3;

    int i;
    int printFlag = 0;
    double *Sdek, *Gder;

    printf("NROWS  ----> %d \n", nrows);
    printf("NRHO  ----> %d \n", nrho);
    printf("RMAX  ----> %lf \n", rmax);
    printf("POTID ----> %d \n", potentialID);
    printf("Cr_ID  ----> %d \n", closureID);
    printf("phi    ----> %lf \n", volumeFactor);
    printf("ALPHA  ----> %lf \n", alpha);
    printf("XNU    ----> %lf \n", xnu);
    printf("EZ     ----> %lf \n", EZ);
    
    // Allocate memory for an array of nrows doubles
    r = malloc(nrows * sizeof(double));
    q = malloc(nrows * sizeof(double));
    U = malloc(nrows*ncols * sizeof(double));
    Up = malloc(nrows*ncols * sizeof(double));
    sigmaVec = malloc(ncols * sizeof(double));
    Sdek = malloc(nrows*2 * sizeof(double));
    Gder = malloc(nrows*2 * sizeof(double));

    // Check if memory allocation was successful
    if (r == NULL || q == NULL || U == NULL || Up == NULL || sigmaVec == NULL || Sdek == NULL || Gder == NULL ) {
        printf("Memory allocation failed 1.\n");
        return 1; // Exit with an error code
    }
    
    // x son las fracciones molrares de cada especie
    x[0] = 1.0;
    x[1] = 1.0 - x[0];

    //creates storage folder 
    char *folderName = getFolderID();

    if (directoryExists(folderName)) {
        printf("The directory exists :( .\n");
    } else {
        if (mkdir(folderName, S_IRWXU) != 0) {
            perror("Error creating directory");  /* This works on ubuntu*/
            return 1; // Exit with an error code
        }
    }

    char outputName1[120] = {0};
    char outputName2[120] = {0};
    
    strcat(outputName1, folderName);
    strcat(outputName1, "/GdeR");
    appendclosureID(outputName1, closureID);
    appendPotentialID(outputName1, potentialID);
    strcat(outputName1, ".dat");
    
    strcat(outputName2, folderName);
    strcat(outputName2, "/SdeK");
    appendclosureID(outputName2, closureID);
    appendPotentialID(outputName2, potentialID);
    strcat(outputName2, ".dat");


    FILE *output1 = fopen(outputName1, "w");
    FILE *output2 = fopen(outputName2, "w");

    // Read input data
    input(volumeFactor, xnu, especie1, especie2, rmax, potentialID);

    // Perform calculations
    OZ2(Sdek, Gder, potentialID, closureID, alpha, EZ, rmax, nrho, folderName, &printFlag);


    for (i = 0; i < nrows; i++) {
        fprintf(output1, "%lf %lf\n", Gder[i*2 + 0], Gder[i*2 + 1]);
        fprintf(output2, "%lf %lf\n", Sdek[i*2 + 0], Sdek[i*2 + 1]);
    }

    fclose(output1);
    fclose(output2);
    
    printf("\n\nRESULTS SAVED AT ---> %s\n", folderName);

    free(r);
    free(q);
    free(U);
    free(Up);
    free(sigmaVec);
    free(Sdek);
    free(Gder);
    free(folderName);

    return 0;
}

void initializeSpecies(species *especie){

    especie->diameter = 0.0;
    especie->temperature = 0.0;
    especie->lambda = 0.0;
    especie->temperature2 = 0.0;
    especie->lambda2 = 0.0;

    return;
}


char *getFolderID() {
    
    char *fullTime, *day, *month, *year, *hour, *min, *sec;

    time_t currentTime;
    currentTime = time(NULL);

    fullTime = ctime(&currentTime);

    month = strtok(fullTime, " "); // find the first double quote
    month = strtok(NULL, " ");    // find the second double quote
    day = strtok(NULL, " ");
    hour = strtok(NULL, ":");
    min = strtok(NULL, ":");
    sec = strtok(NULL, " ");
    year = strtok(NULL, "\n");

    // Calculate the length of the result string
    int resultLength = strlen(day) + strlen(month) + strlen(year) +
                       strlen(hour) + strlen(min) + strlen(sec) + 6; // 6 for underscores and null terminator

    // Allocate memory for fullTimeResult
    char *fullTimeResult = (char *)malloc(resultLength);

    if (fullTimeResult == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1); // Exit the program with an error code
    }

    fullTimeResult[0] = '\0'; // Initialize fullTimeResult as an empty string
    strcat(fullTimeResult, day);
    strcat(fullTimeResult, month);
    strcat(fullTimeResult, year);
    strcat(fullTimeResult, "_");
    strcat(fullTimeResult, hour);
    strcat(fullTimeResult, min);
    strcat(fullTimeResult, sec);

    return fullTimeResult;
}

// check if some directoriy exist
bool directoryExists(char *path) {

    struct stat statbuf;
  
    if (stat(path, &statbuf) != -1) {
        return S_ISDIR(statbuf.st_mode);
    }
  
    return false;
}