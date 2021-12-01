#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "nbody.h"
#include "utility.h"

int main(int argc, char** argv) {
	if (argc != 5 && argc != 6) {
		printf(
			"%s %s\n",
			"Arguments should be of the form",
			"<iterations> <dt> (-b <bodies> | -f <filename>) <number of threads>");
		return 1;
	}

	int iterations = atoi(argv[1]);
	double dt = atof(argv[2]);
	int n_threads = 1;

	if (argc == 6) {
		n_threads = atoi(argv[5]);
	}

	if (iterations <= 0) {
		printf("Iterations must be greater than zero\n");
		return 1;
	} else if (dt <= 0) {
		printf("dt must be greater than zero\n");
		return 1;
	} else if (n_threads <= 0){
		printf("number of threads must be greater than zero\n");
		return 1;
	}

    int n_bodies = 0;
    struct body* bodies = NULL;

	if (strcmp("-b", argv[3]) == 0) {
		n_bodies = atoi(argv[4]);

        if (0 >= n_bodies) {
            printf("Number of randomly generated bodies must be greater than zero\n");\
            return 1;
        }

		double coord_bound = 2.5e11;
		double vel_bound = 5.0e4;

		double mass_upper_bound = 2.0e30;
		double mass_lower_bound = 2.0e23;

		srand(time(NULL));

		bodies = generate_random_bodies(
			coord_bound, vel_bound, mass_lower_bound, mass_upper_bound, n_bodies
		);

	} else if (strcmp("-f", argv[3]) == 0) {
		n_bodies = n_lines(argv[4]);

        if (-1 == n_bodies) {
            printf("File '%s' not found\n", argv[4]);
            return 1;
        } else if (0 == n_bodies) {
			printf("File must contain at least one row of data\n");
			return 1;
		}

		bodies = extract_bodies_from_file(argv[4], n_bodies);

	} else {
        printf("Invalid flag provided. Should be (-f | -b)\n");
        return 1;
    }

	printf("Initial energy: %.6e\n", energy(bodies, n_bodies, n_threads));

    for (int i = 0; i < iterations; i++) {
		step(bodies, n_bodies, dt, n_threads);
	}

    printf("Final energy: %.6e\n", energy(bodies, n_bodies, n_threads));

    free(bodies);

	return 0;
}