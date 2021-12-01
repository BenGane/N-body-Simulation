#ifndef NBODY_H
#define NBODY_H

#define PI (3.141592653589793)
#define SOLARMASS (4 * PI * PI) 
#define NDAYS (365.25) 
#define GCONST (6.67e-11)

/* Main struct */
struct body {
	double x;
	double y;
	double z;
	double x_vel;
	double y_vel;
	double z_vel;
	double next_x_vel;
	double next_y_vel;
	double next_z_vel;
	double mass;
};

/* Thread data for parallel execution of the step method */
struct step_thread_data {
	struct body* bodies;
	int n_bodies;
	int start;
	int finish;
	double dt;
	pthread_barrier_t* barrier;
};

/* Thread data for parallel execution of the energy method */
struct energy_thread_data {
	struct body* bodies;
	int n_bodies;
	int start;
	int finish;
	double total_energy;
};

/* Calculates the distance between two bodies */
extern double distance(struct body* body_i, struct body* body_j);

/* Calculates the magnitude between two bodies */
extern double magnitude(struct body* body_i, struct body* body_j);

/* Constructs the unit vector between two bodies */
extern void construct_unit_vector(struct body* body_i, struct body* body_j, double* dest);

/* Constructs the acceleration vector between two bodies */
extern void construct_acceleration_vector(struct body* body_i, struct body* body_j, double* dest);

/* Calculates the velocity of every body for the next time interval */
extern void compute_next_interval_velocities(struct step_thread_data* data);

/* Updates the position's of all bodies w.r.t their current velocity. 
   Also updates the velocity's of all bodies to be the velocity computed in the next velocities function */
extern void update_positions_and_velocities(struct step_thread_data* data);

/* Advances the simulation with time difference dt using a specified number of threads */
extern void step(struct body* bodies, int n_bodies, double dt, int n_threads);

/* Worker function executed by each individual thread in the step function */
extern void* thread_step(void* arg);

/* Calculates the kinetic energy of a body */
extern double kinetic_energy(struct body* body);

/* Calculates the potential energy between body i and body j */
extern double potential_energy(struct body* body_i, struct body* body_j);

/* Calculates the total energy in the system using a specified number of threads */
extern double energy(struct body* bodies, int n_bodies, int n_threads);

/* Worker function executed by each individual thread in the energy function */
extern void* thread_energy(void* arg);

#endif