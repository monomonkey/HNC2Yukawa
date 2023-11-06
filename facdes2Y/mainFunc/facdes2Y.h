#ifndef FACDES2Y_DOT_H    /* This is an "include guard" */
#define FACDES2Y_DOT_H    /* prevents the file from being included twice. */

#include "../structures/structures.h"
#include "../math/math_aux.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>

int facdes2YFunc(const int nodes, int nrho, double rmax, int potentialID, int closureID, species especie1, species especie2, \
                 double volumeFactor, double alpha, double xnu, double EZ);

void initializeSpecies(species *especie);

char *getFolderID();
bool directoryExists(char *path);

#endif /* FACDES2Y_DOT_H */
