#ifndef UTILITY_H
#define UTILITY_H

#include "nbody.h"

/* Calculates the total number of lines in a file */
extern int n_lines(char* filename);

/* Returns the smallest double */
extern double min(double a, double b);

/* Returns the largest double */
extern double max(double a, double b);

/* Extracts the specified amount of bodies from the given file */
extern struct body* extract_bodies_from_file(char* filename, int n_bodies);

/* Generates a random double within the min-max range */
extern double random_double(double min, double max);

/* Random body generator function */
extern struct body* generate_random_bodies(
    double coord_bound,
    double vel_bound,
    double mass_lower_bound,
    double mass_upper_bound,
    int n_bodies);

/* Determines the maximum mass of all bodies */
extern double max_mass(struct body* bodies, int n_bodies);

/* Determines the maximum absolute coordinate, either x, y or z, of all bodies */
extern double max_abs_coord(struct body* bodies, int n_bodies);

#endif