#include "facdes2Y/mainFunc/facdes2Y.h"
#include <stdio.h>

int main() {

    const int nodes = 1 << 12;  // Equivalent of 2**8
    int potentialID, closureID, nrho;
    double d, volumeFraction, xnu, alpha, EZ, rmax, sigma1, sigma2;
    species especie1, especie2;
    bool IsPolidispersed;

    initializeSpecies(&especie1);
    initializeSpecies(&especie2);

    volumeFraction = 0.2;//0.38;
    alpha = 1.092;
    xnu = 14.0;  // xnu is only used with potentialID = 3 and it's the power used for the Lennard-Jones Potential
    potentialID = 2;
    closureID = 3;
    nrho = 100;
    EZ = 1.0E-4;
    rmax = 160.0;

    IsPolidispersed = false;

    especie1.diameter = 1.0;
    especie1.temperature = 1.0/1.0;//56.1;
    especie1.lambda = 20.0;


    if (IsPolidispersed) {
        // PLEASE FILL THE VARIABLES FOR THE SECOND SPECIES

        especie2.diameter = 1.0;
        especie2.temperature = 1.0/56.1;
        especie2.lambda = 20.0;
        especie2.temperature2 = 1.0;
        especie2. lambda2 = 0.149;

    } else {

        especie2.diameter = especie1.diameter;
        especie2.temperature = especie1.temperature;
        especie2.lambda = especie1.lambda;
        especie2.temperature2 = especie1.temperature2;
        especie2. lambda2 = especie1. lambda2;

    }
    
    facdes2YFunc(nodes, nrho, rmax, potentialID, closureID, especie1, especie2, volumeFraction, alpha, xnu, EZ);

    return 0;
}
