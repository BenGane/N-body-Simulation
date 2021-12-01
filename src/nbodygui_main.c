#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <SDL2/SDL.h>
#include <SDL2/SDL2_gfxPrimitives.h>
#include "nbody.h"
#include "utility.h"

#define MAX_DIAMETER (20)

int main(int argc, char** argv) {
	if (argc != 8 && argc != 9) {
		printf(
			"%s %s\n",
			"Arguments should be of the form",
			"<resolution_width> <resolution_height> <iterations> <dt> (-b <bodies> | -f <filename>) <scale> <number of threads>");
		return 1;
	}

	double VIEW_WIDTH = atof(argv[1]);
	double VIEW_HEIGHT = atof(argv[2]);
	int iterations = atoi(argv[3]);
	double dt = atof(argv[4]);
	double scale = atof(argv[7]);
	int n_threads = 1;

	if (argc == 9) {
		n_threads = atoi(argv[8]);
	}

	int n_bodies = 0;
	struct body* bodies = NULL;

	if (VIEW_WIDTH <= 100 || VIEW_HEIGHT <= 100) {
		printf("Resolution width and height must be greater than 100");
		return 1;
	} else if (iterations <= 0) {
		printf("Iterations must be greater than zero\n");
		return 1;
	} else if (dt <= 0) {
		printf("dt must be greater than zero\n");
		return 1;
	} else if (scale <= 0){
		printf("Scale must be greater than zero\n");
		return 1;
	} else if (n_threads <= 0) {
		printf("Number of threads must be greater than zero\n");
		return 1;
	}

	if (strcmp("-b", argv[5]) == 0) {
		/* Run simulation with randomly generated bodies */

		n_bodies = atoi(argv[6]);

		if (0 >= n_bodies) {
            printf("Number of randomly generated bodies must be greater than zero.\n");
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

	} else if (strcmp("-f", argv[5]) == 0) {
		/* Run simulation with bodies extracted from file */
		n_bodies = n_lines(argv[6]);

		if (-1 == n_bodies) {
            printf("File '%s' not found.\n", argv[5]);
            return 1;
        } else if (0 == n_bodies) {
			printf("File must contain at least one row of data\n");
			return 1;
		}

		bodies = extract_bodies_from_file(argv[6], n_bodies);

	} else {
        printf("Invalid flag provided. Should be (-f | -b).\n");
        return 1;
    }

	SDL_Window* window = NULL;
	SDL_Renderer* renderer = NULL;
	SDL_Event event;
	
	if (SDL_Init(SDL_INIT_VIDEO) < 0) {
		printf("Error: SDL_Init returned a value < 0\n");
		return 1;
	}

	window = SDL_CreateWindow("SDL Template",
		SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
		VIEW_WIDTH, VIEW_HEIGHT, SDL_WINDOW_SHOWN);

	if (NULL == window) {
		printf("Error: SDL_CreateWindow() returned NULL\n");
		return 1;
	}

	renderer = SDL_CreateRenderer(window, -1,
		SDL_RENDERER_ACCELERATED 
		| SDL_RENDERER_PRESENTVSYNC);

	if (NULL == renderer) {
		printf("Error: SDL_CreateRenderer() returned NULL\n");
		return 1;
	}

	/* Scale used to make sure the bodies are displayed at a maximum of 10 pixels */
	double diameter_scale = MAX_DIAMETER / max_mass(bodies, n_bodies);
	
	/* Max absolute coordinate, x, y or z, in the system */
	double max_coord = max_abs_coord(bodies, n_bodies);

	/* Additional scaling to ensure the bodies are displayed nicely on the screen */
	double coord_scale = (VIEW_WIDTH / 2) / (max(VIEW_WIDTH / 2, max_coord));
	
	for (int i = 0; i < iterations; i++) {

		SDL_SetRenderDrawColor(renderer, 0, 0, 0, 0xFF);
		SDL_RenderClear(renderer);
		
		/* Display each body on the screen */
		for (int j = 0; j < n_bodies; j++) {
			struct body* body = bodies + j;

			/* Display coordinates are relative to the center of the screen */
			double x = VIEW_WIDTH / 2 + (scale * coord_scale * body->x);
			double y = VIEW_HEIGHT / 2 + -1 * (scale * coord_scale * body->y);

			double diameter = max(2, body->mass * diameter_scale);

			/* Display each body relative to the center of the view */
			filledCircleColor(renderer, x, y, diameter, 0xFF0000FF);
		}

		if (SDL_PollEvent(&event)) {
			if (SDL_QUIT == event.type) {
				break;
			}
		}
		
		SDL_RenderPresent(renderer);

		step(bodies, n_bodies, dt, n_threads);
	}

	SDL_DestroyRenderer(renderer);
	SDL_DestroyWindow(window);
	SDL_Quit();	
	free(bodies);

	return 0;	
}