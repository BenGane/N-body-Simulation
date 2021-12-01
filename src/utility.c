#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <pthread.h>
#include "nbody.h"
#include "utility.h"

#define BUFFER_SIZE 2048

/* Calculates the total number of lines in a file */
int n_lines(char* filename) {
    FILE* f = fopen(filename, "r");

    if (NULL == f) {
        return -1;
    }

    char buffer[BUFFER_SIZE];

    int count = 0;

    while (fgets(buffer, BUFFER_SIZE, f)) {
        count++;
    }

    fclose(f);

    return count;
}

/* Returns the smallest double */
double min(double a, double b) {
    if (a < b) {
        return a;
    }
    return b;
}

/* Returns the largest double */
double max(double a, double b) {
    if (a > b) {
        return a;
    }
    return b;
}

/* Generates a random double within the min-max range */
double random_double(double min, double max) {
    double range = max - min;
    double div = RAND_MAX / range;
    return min + (rand() / div);
}

/* Extracts the specified amount of bodies from the given file */
struct body* extract_bodies_from_file(char* filename, int n_bodies) {
    FILE* f = fopen(filename, "r");

    if (NULL == f) {
        return NULL;
    }

    char buffer[BUFFER_SIZE];

    struct body* bodies = (struct body*) malloc(sizeof(struct body) * n_bodies);
    int i = 0;

    while (fgets(buffer, BUFFER_SIZE, f)) {
        struct body* body = bodies + i;
        sscanf(
            buffer, 
            "%lf,%lf,%lf,%lf,%lf,%lf,%lf",
            &body->x,
            &body->y,
            &body->z,
            &body->x_vel,
            &body->y_vel,
            &body->z_vel,
            &body->mass
        );
        i++;
    }

    return bodies;
}

/* Random body generator function */
struct body* generate_random_bodies(
    double coord_bound,
    double vel_bound,
    double mass_lower_bound,
    double mass_upper_bound,
    int n_bodies)    
{
    struct body* bodies = (struct body*) malloc(sizeof(struct body) * n_bodies);

    for (int i = 0; i < n_bodies; i++) {
		struct body* body = bodies + i; 
        body->x = random_double(-coord_bound, coord_bound);
        body->y = random_double(-coord_bound, coord_bound); 
        body->z = random_double(-coord_bound, coord_bound);
        body->x_vel= random_double(-vel_bound, vel_bound);
        body->y_vel = random_double(-vel_bound, vel_bound);
        body->z_vel = random_double(-vel_bound, vel_bound);
        body->mass = random_double(mass_lower_bound, mass_upper_bound);
	}   

    return bodies;
}

/* Determines the maximum mass of all bodies */
double max_mass(struct body* bodies, int n_bodies) {
	double max = bodies->mass;
	for (int i = 1; i < n_bodies; i++) {
		double current_mass = (bodies + i)->mass;
		if (current_mass > max) {
			max = current_mass;
		}
	}
	return max;
}

/* Determines the maximum absolute coordinate, x, y or z, of all bodies */
double max_abs_coord(struct body* bodies, int n_bodies) {
    double max_coord = max(fabs(bodies->x), max(fabs(bodies->y), fabs(bodies->z)));
    for (int i = 1; i < n_bodies; i++) {
        double abs_x = fabs((bodies + i)->x); 
        double abs_y = fabs((bodies + i)->y);
        double abs_z = fabs((bodies + i)->z);
        max_coord = max(max_coord, max(abs_x, max(abs_y, abs_z)));
    }
    return max_coord;
}