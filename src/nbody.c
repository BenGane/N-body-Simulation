#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include "nbody.h"
#include "utility.h"

/* Calculates the distance between two bodies */
double distance(struct body* body_i, struct body* body_j) {
	double dx = body_i->x - body_j->x;
	double dy = body_i->y - body_j->y;
	double dz = body_i->z - body_j->z;

	double result = sqrt(dx * dx + dy * dy + dz * dz);

	if (result == 0) {
		result = 1e-10;
	} 

	return result;
}

/* Calculates the magnitude between two bodies */
double magnitude(struct body* body_i, struct body* body_j) {
	double dist = distance(body_i, body_j);
	return GCONST * body_i->mass * body_j->mass / (dist * dist);
}

/* Constructs a unit vector between two bodies */
void construct_unit_vector(struct body* body_i, struct body* body_j, double* dest) {
	double dx = body_i->x - body_j->x;
	double dy = body_i->y - body_j->y;
	double dz = body_i->z - body_j->z;
	double dist = distance(body_i, body_j);

	dest[0] = dx / dist;
	dest[1] = dy / dist;
	dest[2] = dz / dist;
}

/* Constructs the acceleration vector between two bodies */
void construct_acceleration_vector(struct body* body_i, struct body* body_j, double* dest) {
	double unit_vector[3];
	construct_unit_vector(body_i, body_j, unit_vector);

	double scale = magnitude(body_i, body_j) / body_j->mass;

	for (int i = 0; i < 3; i++) {
		dest[i] = unit_vector[i] * scale;
	}
}

/* Calculates the velocity of every body for the next time interval */
void compute_next_interval_velocities(struct step_thread_data* data) {
	for (int i = data->start; i < data->finish; i++) {
		struct body* body_i = data->bodies + i;

		register double next_x_vel = body_i->x_vel;
		register double next_y_val = body_i->y_vel;
		register double next_z_val = body_i->z_vel;

		for (int j = 0; j < data->n_bodies; j++) {
			if (i == j) {
				continue;
			}
			struct body* body_j = data->bodies + j;
			double acceleration_vector[3];
			construct_acceleration_vector(body_j, body_i, acceleration_vector);
			next_x_vel += acceleration_vector[0] * data->dt;
			next_y_val += acceleration_vector[1] * data->dt;
			next_z_val += acceleration_vector[2] * data->dt;
		}

		body_i->next_x_vel = next_x_vel;
		body_i->next_y_vel = next_y_val;
		body_i->next_z_vel = next_z_val;
	}
}

/* Updates the position's of all bodies w.r.t their current velocity. 
   Also updates the velocity's of all bodies to be the velocity computed in the next velocities function */
void update_positions_and_velocities(struct step_thread_data* data) {
	for (int i = data->start; i < data->finish; i++) {
		struct body* body = data->bodies + i;

		/* Update the body's position */
		body->x += body->x_vel * data->dt;
		body->y += body->y_vel * data->dt;
		body->z += body->z_vel * data->dt;

		/* Update the body's velocity */
		body->x_vel = body->next_x_vel;
		body->y_vel = body->next_y_vel;
		body->z_vel = body->next_z_vel;
	}
}

/* Advances the simulation with time difference dt using a specified number of threads */
void step(struct body* bodies, int n_bodies, double dt, int n_threads) {
	
	/* Make sure we don't have more threads than we do data */
	n_threads = min(n_threads, n_bodies);
	
	pthread_barrier_t barrier;
	pthread_barrier_init(&barrier, NULL, n_threads);

	pthread_t* threads = (pthread_t*) malloc(sizeof(pthread_t) * n_threads);
	struct step_thread_data* step_thread_data = 
				(struct step_thread_data*) malloc(sizeof(struct step_thread_data) * n_threads);	

	/* Number of bodies used in each individual thread */
	int n_bodies_per_thread = n_bodies / n_threads;

	for (int i = 0; i < n_threads; i++) {
		struct step_thread_data* data = step_thread_data + i;
		data->bodies = bodies;
		data->n_bodies = n_bodies;
		data->start = i * n_bodies_per_thread;
		data->finish = i == n_threads - 1 ? n_bodies : data->start + n_bodies_per_thread;
		data->dt = dt;
		data->barrier = &barrier;
	}

	for (int i = 0; i < n_threads; i++) {
		pthread_create(threads + i, NULL, thread_step, step_thread_data + i);
	}

	for (int i = 0; i < n_threads; i++) {
		pthread_join(threads[i], NULL);
	}

	free(threads);
	free(step_thread_data);

	pthread_barrier_destroy(&barrier);
}

/* Worker function executed by each individual thread */
void* thread_step(void* arg) {
	struct step_thread_data* data = (struct step_thread_data*) arg;
	compute_next_interval_velocities(data);
	pthread_barrier_wait(data->barrier);
	update_positions_and_velocities(data);
	return NULL;
}

/* Calculates the total energy in the system using a specified number of threads */
double energy(struct body* bodies, int n_bodies, int n_threads) {
	
	/* Make sure we don't have more threads than we do data */
	n_threads = min(n_threads, n_bodies);
	
	pthread_t* threads = (pthread_t*) malloc(sizeof(pthread_t) * n_threads);
	struct energy_thread_data* energy_thread_data = 
			(struct energy_thread_data*) malloc(sizeof(struct energy_thread_data) * n_threads);	

	/* Number of bodies used in each individual thread */
	int n_bodies_per_thread = n_bodies / n_threads;

	for (int i = 0; i < n_threads; i++) {
		struct energy_thread_data* data = energy_thread_data + i;
		data->bodies = bodies;
		data->n_bodies = n_bodies;
		data->start = i * n_bodies_per_thread;
		data->finish = i == n_threads - 1 ? n_bodies : data->start + n_bodies_per_thread;
		data->total_energy = 0;
	}

	for (int i = 0; i < n_threads; i++) {
		pthread_create(threads + i, NULL, thread_energy, energy_thread_data + i);
	}

	for (int i = 0; i < n_threads; i++) {
		pthread_join(threads[i], NULL);
	}

	double total_energy = 0;

	for (int i = 0; i < n_threads; i++) {
		struct energy_thread_data* data = energy_thread_data + i;
		total_energy += data->total_energy;
	}

	free(threads);
	free(energy_thread_data);
	
	return total_energy;
}

/* Calculates the kinetic energy of a body */
double kinetic_energy(struct body* body) {
	double velocity_squared = 
		body->x_vel * body->x_vel + 
		body->y_vel * body->y_vel +
		body->z_vel * body->z_vel;
	return body->mass * velocity_squared / 2;
}

/* Calculates the potential energy between body i and body j */
double potential_energy(struct body* body_i, struct body* body_j) {
	return -1 * GCONST * body_i->mass * body_j->mass / distance(body_i, body_j);
}

/* Worker function executed by each individual thread in the energy function */
void* thread_energy(void* arg) {
	struct energy_thread_data* data = (struct energy_thread_data*) arg;

	register double total_energy = 0;

	for (int i = data->start; i < data->finish; i++) {
		struct body* body_i = data->bodies + i;

		total_energy += kinetic_energy(body_i);

		for (int j = i + 1; j < data->n_bodies; j++) {
			struct body* body_j = data->bodies + j;
			total_energy += potential_energy(body_i, body_j);
		}
	}

	data->total_energy = total_energy;

	return NULL;
}