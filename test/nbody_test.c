#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <cmocka.h>
#include <pthread.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "../src/nbody.h"
#include "../src/utility.h"

int verification(double* a, double* b, int num_elements) {
    for (int i = 0; i < num_elements; i++) {
        if (a[i] != b[i]) {
            return 0;
        }
    }
    return 1;
}

int body_equals(struct body* body_i, struct body* body_j) {
    return body_i->x == body_j->x &&
           body_i->y == body_j->y &&
           body_i->z == body_j->z &&
           body_i->x_vel == body_j->x_vel &&
           body_i->y_vel == body_j->y_vel &&
           body_i->z_vel == body_j->z_vel;
}

void distance_test_one() {
    struct body body_i = {
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    }; 

    struct body body_j = {
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    };

    double expected = 1e-10;
    double actual = distance(&body_i, &body_j);

    assert_true(fabs(expected - actual) == 0);
}

void distance_test_two() {
    struct body body_i = {
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    }; 

    struct body body_j = {
        1, 1, 0, 0, 0, 0, 0, 0, 0, 0
    };

    double expected = sqrt(2);
    double actual = distance(&body_i, &body_j);

    assert_true(fabs(expected - actual) == 0);
}

void distance_test_three() {
    struct body body_i = {
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    }; 

    struct body body_j = {
        1, 1, 1, 0, 0, 0, 0, 0, 0, 0
    };

    double expected = sqrt(3);
    double actual = distance(&body_i, &body_j);

    assert_true(fabs(expected - actual) == 0);
}

void distance_test_four() {
    struct body body_i = {
        -1000, -1000, 1000, 0, 0, 0, 0, 0, 0, 0
    }; 

    struct body body_j = {
        100, 100, 100, 0, 0, 0, 0, 0, 0, 0
    };

    double expected = sqrt(pow(1100, 2) + pow(1100, 2) + pow(900, 2));
    double actual = distance(&body_i, &body_j);

    assert_true(fabs(expected - actual) == 0);
}

void distance_test_five() {
    struct body body_i = {
        -1.25e20, -2.5e21, 4.5e15, 0, 0, 0, 0, 0, 0, 0
    }; 

    struct body body_j = {
        5.6e19, 2.4e24, 2.1e16, 0, 0, 0, 0, 0, 0, 0
    };

    double expected = sqrt(
        pow(-1.25e20 - 5.6e19, 2) + 
        pow(-2.5e21 - 2.4e24, 2) + 
        pow(4.5e15 - 2.1e16, 2)
    );

    double actual = distance(&body_i, &body_j);

    assert_true(fabs(expected - actual) < 1e-10);
}

void distance_test_six() {
    struct body body_i = {
        -1.25e20, -2.5e21, 4.5e15, 0, 0, 0, 0, 0, 0, 0
    }; 

    double expected = 1e-10;

    double actual = distance(&body_i, &body_i);

    assert_true(fabs(expected - actual) == 0);
}

void magnitude_test_one() {
    struct body body_i = {
        -9e10, -7e15, -1e20, 0, 0, 0, 0, 0, 0, 0
    }; 

    struct body body_j = {
        6e20, 2e24, -1e12, 0, 0, 0, 0, 0, 0, 0
    };

    double expected = 0;
    double actual = magnitude(&body_i, &body_j);

    assert_true(fabs(expected - actual) == 0);
}

void magnitude_test_two() {
    struct body body_i = {
        -9e10, -7e15, -1e20, 0, 0, 0, 0, 0, 0, 4e15
    }; 

    struct body body_j = {
        6e20, 2e24, -1e12, 0, 0, 0, 0, 0, 0, 2e17
    };

    double expected = GCONST * 1 * 4e15 * 2e17 / (pow(distance(&body_i, &body_j), 2));
    double actual = magnitude(&body_i, &body_j);

    assert_true(fabs(expected - actual) == 0);
}

void magnitude_test_three() {
    struct body body_i = {
        -100, 0, 0, 0, 0, 0, 0, 0, 0, 1e10
    }; 

    struct body body_j = {
        100, 0, 0, 0, 0, 0, 0, 0, 0, 1e10
    };

    double expected = GCONST * 1 * 1e10 * 1e10 / (pow(distance(&body_i, &body_j), 2));
    double actual = magnitude(&body_i, &body_j);

    assert_true(fabs(expected - actual) < 1e-10);
}

void magnitude_test_four() {
    struct body body_i = {
        -100, 0, 0, 0, 0, 0, 0, 0, 0, 0
    }; 

    struct body body_j = {
        100, 0, 0, 0, 0, 0, 0, 0, 0, 9e36
    };

    double expected = 0;
    double actual = magnitude(&body_i, &body_j);

    assert_true(fabs(expected - actual) == 0);
}

void magnitude_test_five() {
    struct body body_i = {
        -100, 0, 0, 0, 0, 0, 0, 0, 0, 9e36
    }; 

    struct body body_j = {
        100, 0, 0, 0, 0, 0, 0, 0, 0, 9e36
    };

    double expected = GCONST * 9e36 * 9e36 / pow(distance(&body_i, &body_j), 2);
    double actual = magnitude(&body_i, &body_j);

    assert_true(fabs(expected - actual) == 0);
}

void unit_vector_test_one() {
    struct body body_i = {
        -1000, 0, 0, 0, 0, 0, 0, 0, 0, 0
    }; 

    struct body body_j = {
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    };

    double actual[3];
    construct_unit_vector(&body_i, &body_j, actual);

    double expected[] = { -1, 0, 0 };

    assert_true(verification(actual, expected, 3));
}

void unit_vector_test_two() {
    struct body body_i = {
        -2.3e24, 1.2e10, 3.4e9, 0, 0, 0, 0, 0, 0, 0
    }; 

    struct body body_j = {
        -9e10, -9e23, 8.7e12, 0, 0, 0, 0, 0, 0, 0
    };

    double actual[3];
    construct_unit_vector(&body_i, &body_j, actual);

    double dist = distance(&body_i, &body_j);

    double expected[] = { 
        (-2.3e24 - -9e10) / dist, 
        (1.2e10 - -9e23) / dist, 
        (3.4e9 - 8.7e12) / dist 
    };

    assert_true(verification(actual, expected, 3));
}

void unit_vector_test_three() {
    struct body body_i = {
        0, 0, 1e10, 0, 0, 0, 0, 0, 0, 0
    }; 

    struct body body_j = {
        0, 0, -1e10, 0, 0, 0, 0, 0, 0, 0
    };

    double actual[3];
    construct_unit_vector(&body_i, &body_j, actual);

    double expected[] = { 0, 0, 1 };

    assert_true(verification(actual, expected, 3));
}

void unit_vector_test_four() {
    struct body body_i = {
        1e10, -1e10, 1e10, 0, 0, 0, 0, 0, 0, 0
    }; 

    struct body body_j = {
        -1e10, 1e10, -1e10, 0, 0, 0, 0, 0, 0, 0
    };

    double actual[3];
    construct_unit_vector(&body_i, &body_j, actual);

    double dist = distance(&body_i, &body_j);

    double expected[] = {
        2e10 / dist,
        -2e10 / dist,
        2e10 / dist
     };

    assert_true(verification(actual, expected, 3));
}

void unit_vector_test_five() {
    struct body body_i = {
        0, 1, 0, 0, 0, 0, 0, 0, 0, 0
    }; 

    struct body body_j = {
        0, -1, 0, 0, 0, 0, 0, 0, 0, 0
    };

    double actual[3];
    construct_unit_vector(&body_i, &body_j, actual);

    double dist = distance(&body_i, &body_j);

    double expected[] = { 0, 1, 0 };

    assert_true(verification(actual, expected, 3));
}

void acceleration_vector_test_one() {
    struct body body_i = {
        -1000, 0, 0, 0, 0, 0, 0, 0, 0, 100
    }; 

    struct body body_j = {
        0, 0, 0, 0, 0, 0, 0, 0, 0, 1000
    };

    double actual[3];
    construct_acceleration_vector(&body_i, &body_j, actual);

    double unit_vect[3];
    construct_unit_vector(&body_i, &body_j, unit_vect);

    double scale = magnitude(&body_i, &body_j) / body_j.mass;

    double expected[] = { 
        unit_vect[0] * scale, 
        unit_vect[1] * scale, 
        unit_vect[2] * scale 
    };

    assert_true(verification(actual, expected, 3));
}

void acceleration_vector_test_two() {
    struct body body_i = {
        -1e21, 2e23, 1e42, 5e4, 4e4, 3.3e2, 0, 0, 0, 100
    }; 

    struct body body_j = {
        1e23, -1e24, -4.2e10, -4.3e3, 6e4, 2.23e5, 0, 0, 0, 1000
    };

    double actual[3];
    construct_acceleration_vector(&body_i, &body_j, actual);

    double unit_vect[3];
    construct_unit_vector(&body_i, &body_j, unit_vect);

    double scale = magnitude(&body_i, &body_j) / body_j.mass;

    double expected[] = { 
        unit_vect[0] * scale, 
        unit_vect[1] * scale, 
        unit_vect[2] * scale 
    };

    assert_true(verification(actual, expected, 3));
}

void acceleration_vector_test_three() {
    struct body body_i = {
        -1, 0, 0, 0, 0, 0, 0, 0, 0, 1e25
    }; 

    struct body body_j = {
         1, 0, 0, 0, 0, 0, 0, 0, 0, 1e10
    };

    double actual[3];
    construct_acceleration_vector(&body_i, &body_j, actual);

    double unit_vect[3];
    construct_unit_vector(&body_i, &body_j, unit_vect);

    double scale = magnitude(&body_i, &body_j) / body_j.mass;

    double expected[] = { 
        unit_vect[0] * scale, 
        unit_vect[1] * scale, 
        unit_vect[2] * scale 
    };

    assert_true(verification(actual, expected, 3));
}

void acceleration_vector_test_four() {
    struct body body_i = {
        1e20, 1e20, 1e20, 1e5, 1e5, 1e5, 0, 0, 0, 1e30
    }; 

    struct body body_j = {
         -1e20, -1e20, -1e20, -1e5, -1e5, -1e5, 0, 0, 0, 1e30
    };

    double actual[3];
    construct_acceleration_vector(&body_i, &body_j, actual);

    double unit_vect[3];
    construct_unit_vector(&body_i, &body_j, unit_vect);

    double scale = magnitude(&body_i, &body_j) / body_j.mass;

    double expected[] = { 
        unit_vect[0] * scale, 
        unit_vect[1] * scale, 
        unit_vect[2] * scale 
    };

    assert_true(verification(actual, expected, 3));
}

void acceleration_vector_test_five() {
    struct body body_i = {
        -1, 0, 0, 0, 0, 0, 0, 0, 0, 1
    }; 

    struct body body_j = {
        1, 0, 0, 0, 0, 0, 0, 0, 0, -1
    };

    double actual[3];
    construct_acceleration_vector(&body_i, &body_j, actual);

    double unit_vect[3];
    construct_unit_vector(&body_i, &body_j, unit_vect);

    double scale = magnitude(&body_i, &body_j) / body_j.mass;

    double expected[] = { 
        unit_vect[0] * scale, 
        unit_vect[1] * scale, 
        unit_vect[2] * scale 
    };

    assert_true(verification(actual, expected, 3));
}

void extract_bodies_from_file_test() {
    char* filename = "planets.csv";

    int n_bodies = n_lines(filename);

    struct body* actual = extract_bodies_from_file(filename, n_bodies);

    struct body body_a = {
        0, 0, 0, 0, 0, 0, 0, 0, 0, 1988999999999999901909255192576.0
    };

    struct body body_b = {
        150000000000.0, 0, 0, 0, 29800.0, 0, 0, 0, 0, 5974000000000000373293056.0
    };

    struct body body_c = {
        230000000000.0, 0, 0, 0, 24100.0, 0, 0, 0, 0, 641899999999999963299840.0
    };

    struct body body_d = {
        55000000000.0, 0, 0, 0, 47900.0, 0, 0, 0, 0, 330199999999999993708544.0
    };

    struct body body_e = {
        100000000000.0, 0, 0, 0, 35000.0, 0, 0, 0, 0, 4869000000000000115343360.0
    };


    struct body expected[] = {
        body_a, body_b, body_c, body_d, body_e
    };

    for (int i = 0; i < n_bodies; i++) {
        assert_true(body_equals(expected + i, actual + i));
    }

    free(actual);
}

void step_test_one() {
    char* filename = "planets.csv";

    int dt = 0.01;
    int iterations = 100;
    int n_threads = 1;

    int n_bodies = n_lines(filename);
    struct body* bodies = extract_bodies_from_file(filename, n_bodies);

    double initial_energy = energy(bodies, n_bodies, n_threads);

    for (int i = 0; i < iterations; i++) {
        step(bodies, n_bodies, dt, n_threads);
    }

    double final_energy = energy(bodies, n_bodies, n_threads);
    double delta = fabs(0.01 * initial_energy);

    assert_true(fabs(initial_energy - final_energy) < delta);

    free(bodies);
}

void step_test_two() {
    char* filename = "planets.csv";

    int dt = 0.01;
    int iterations = 1000;
    int n_threads = 1;

    int n_bodies = n_lines(filename);
    struct body* bodies = extract_bodies_from_file(filename, n_bodies);

    double initial_energy = energy(bodies, n_bodies, n_threads);

    for (int i = 0; i < iterations; i++) {
        step(bodies, n_bodies, dt, n_threads);
    }

    double final_energy = energy(bodies, n_bodies, n_threads);
    double delta = fabs(0.01 * initial_energy);

    assert_true(fabs(initial_energy - final_energy) < delta);

    free(bodies);
}

void step_test_three() {
    char* filename = "planets.csv";

    int dt = 0.01;
    int iterations = 10000;
    int n_threads = 1;

    int n_bodies = n_lines(filename);
    struct body* bodies = extract_bodies_from_file(filename, n_bodies);

    double initial_energy = energy(bodies, n_bodies, n_threads);

    for (int i = 0; i < iterations; i++) {
        step(bodies, n_bodies, dt, n_threads);
    }

    double final_energy = energy(bodies, n_bodies, n_threads);
    double delta = fabs(0.01 * initial_energy);

    assert_true(fabs(initial_energy - final_energy) < delta);

    free(bodies);
}

void step_test_four() {
    char* filename = "planets.csv";

    int dt = 1;
    int iterations = 100;
    int n_threads = 1;

    int n_bodies = n_lines(filename);
    struct body* bodies = extract_bodies_from_file(filename, n_bodies);

    double initial_energy = energy(bodies, n_bodies, n_threads);

    for (int i = 0; i < iterations; i++) {
        step(bodies, n_bodies, dt, n_threads);
    }

    double final_energy = energy(bodies, n_bodies, n_threads);
    double delta = fabs(0.01 * initial_energy);

    assert_true(fabs(initial_energy - final_energy) < delta);

    free(bodies);
}

void step_test_five() {
    char* filename = "planets.csv";

    int dt = 1;
    int iterations = 1000;
    int n_threads = 1;

    int n_bodies = n_lines(filename);
    struct body* bodies = extract_bodies_from_file(filename, n_bodies);

    double initial_energy = energy(bodies, n_bodies, n_threads);

    for (int i = 0; i < iterations; i++) {
        step(bodies, n_bodies, dt, n_threads);
    }

    double final_energy = energy(bodies, n_bodies, n_threads);
    double delta = fabs(0.01 * initial_energy);

    assert_true(fabs(initial_energy - final_energy) < delta);

    free(bodies);
}

void step_test_six() {
    char* filename = "planets.csv";

    int dt = 1;
    int iterations = 10000;
    int n_threads = 1;

    int n_bodies = n_lines(filename);
    struct body* bodies = extract_bodies_from_file(filename, n_bodies);

    double initial_energy = energy(bodies, n_bodies, n_threads);

    for (int i = 0; i < iterations; i++) {
        step(bodies, n_bodies, dt, n_threads);
    }

    double final_energy = energy(bodies, n_bodies, n_threads);
    double delta = fabs(0.01 * initial_energy);

    assert_true(fabs(initial_energy - final_energy) < delta);

    free(bodies);
}

void step_test_seven() {
    char* filename = "planets.csv";

    int dt = 0.01;
    int iterations = 100;
    int n_threads = 4;

    int n_bodies = n_lines(filename);
    struct body* bodies = extract_bodies_from_file(filename, n_bodies);

    double initial_energy = energy(bodies, n_bodies, n_threads);

    for (int i = 0; i < iterations; i++) {
        step(bodies, n_bodies, dt, n_threads);
    }

    double final_energy = energy(bodies, n_bodies, n_threads);
    double delta = fabs(0.01 * initial_energy);

    assert_true(fabs(initial_energy - final_energy) < delta);

    free(bodies);
}

void step_test_eight() {
    char* filename = "planets.csv";

    int dt = 0.01;
    int iterations = 1000;
    int n_threads = 4;

    int n_bodies = n_lines(filename);
    struct body* bodies = extract_bodies_from_file(filename, n_bodies);

    double initial_energy = energy(bodies, n_bodies, n_threads);

    for (int i = 0; i < iterations; i++) {
        step(bodies, n_bodies, dt, n_threads);
    }

    double final_energy = energy(bodies, n_bodies, n_threads);
    double delta = fabs(0.01 * initial_energy);

    assert_true(fabs(initial_energy - final_energy) < delta);

    free(bodies);
}

void step_test_nine() {
    char* filename = "planets.csv";

    int dt = 0.01;
    int iterations = 10000;
    int n_threads = 4;

    int n_bodies = n_lines(filename);
    struct body* bodies = extract_bodies_from_file(filename, n_bodies);

    double initial_energy = energy(bodies, n_bodies, n_threads);

    for (int i = 0; i < iterations; i++) {
        step(bodies, n_bodies, dt, n_threads);
    }

    double final_energy = energy(bodies, n_bodies, n_threads);
    double delta = fabs(0.01 * initial_energy);

    assert_true(fabs(initial_energy - final_energy) < delta);

    free(bodies);
}

void step_test_ten() {
    char* filename = "planets.csv";

    int dt = 1;
    int iterations = 100;
    int n_threads = 4;

    int n_bodies = n_lines(filename);
    struct body* bodies = extract_bodies_from_file(filename, n_bodies);

    double initial_energy = energy(bodies, n_bodies, n_threads);

    for (int i = 0; i < iterations; i++) {
        step(bodies, n_bodies, dt, n_threads);
    }

    double final_energy = energy(bodies, n_bodies, n_threads);
    double delta = fabs(0.01 * initial_energy);

    assert_true(fabs(initial_energy - final_energy) < delta);

    free(bodies);
}

void step_test_eleven() {
    char* filename = "planets.csv";

    int dt = 1;
    int iterations = 1000;
    int n_threads = 4;

    int n_bodies = n_lines(filename);
    struct body* bodies = extract_bodies_from_file(filename, n_bodies);

    double initial_energy = energy(bodies, n_bodies, n_threads);

    for (int i = 0; i < iterations; i++) {
        step(bodies, n_bodies, dt, n_threads);
    }

    double final_energy = energy(bodies, n_bodies, n_threads);
    double delta = fabs(0.01 * initial_energy);

    assert_true(fabs(initial_energy - final_energy) < delta);

    free(bodies);
}

void step_test_twelve() {
    char* filename = "planets.csv";

    int dt = 1;
    int iterations = 10000;
    int n_threads = 4;

    int n_bodies = n_lines(filename);
    struct body* bodies = extract_bodies_from_file(filename, n_bodies);

    double initial_energy = energy(bodies, n_bodies, n_threads);

    for (int i = 0; i < iterations; i++) {
        step(bodies, n_bodies, dt, n_threads);
    }

    double final_energy = energy(bodies, n_bodies, n_threads);
    double delta = fabs(0.01 * initial_energy);

    assert_true(fabs(initial_energy - final_energy) < delta);

    free(bodies);
}

void step_test_thirteen() {
    int dt = 0.01;
    int iterations = 10;
    int n_threads = 1;
    int n_bodies = 100;

    struct body* bodies = 
        generate_random_bodies(
			2.5e11, 5.0e4, 2.0e30, 2.0e23, n_bodies
		);

    double initial_energy = energy(bodies, n_bodies, n_threads);

    for (int i = 0; i < iterations; i++) {
        step(bodies, n_bodies, dt, n_threads);
    }

    double final_energy = energy(bodies, n_bodies, n_threads);
    double delta = fabs(0.01 * initial_energy);

    assert_true(fabs(initial_energy - final_energy) < delta);

    free(bodies);
}

void step_test_fourteen() {
    int dt = 0.01;
    int iterations = 100;
    int n_threads = 1;
    int n_bodies = 100;

    struct body* bodies = 
        generate_random_bodies(
			2.5e11, 5.0e4, 2.0e30, 2.0e23, n_bodies
		);

    double initial_energy = energy(bodies, n_bodies, n_threads);

    for (int i = 0; i < iterations; i++) {
        step(bodies, n_bodies, dt, n_threads);
    }

    double final_energy = energy(bodies, n_bodies, n_threads);
    double delta = fabs(0.01 * initial_energy);

    assert_true(fabs(initial_energy - final_energy) < delta);

    free(bodies);
}

void step_test_fifteen() {
    int dt = 0.01;
    int iterations = 1000;
    int n_threads = 1;
    int n_bodies = 100;

    struct body* bodies = 
        generate_random_bodies(
			2.5e11, 5.0e4, 2.0e30, 2.0e23, n_bodies
		);

    double initial_energy = energy(bodies, n_bodies, n_threads);

    for (int i = 0; i < iterations; i++) {
        step(bodies, n_bodies, dt, n_threads);
    }

    double final_energy = energy(bodies, n_bodies, n_threads);
    double delta = fabs(0.01 * initial_energy);

    assert_true(fabs(initial_energy - final_energy) < delta);

    free(bodies);
}

void step_test_sixteen() {
    int dt = 1;
    int iterations = 10;
    int n_threads = 1;
    int n_bodies = 100;

    struct body* bodies = 
        generate_random_bodies(
			2.5e11, 5.0e4, 2.0e30, 2.0e23, n_bodies
		);

    double initial_energy = energy(bodies, n_bodies, n_threads);

    for (int i = 0; i < iterations; i++) {
        step(bodies, n_bodies, dt, n_threads);
    }

    double final_energy = energy(bodies, n_bodies, n_threads);
    double delta = fabs(0.01 * initial_energy);

    assert_true(fabs(initial_energy - final_energy) < delta);

    free(bodies);
}

void step_test_seventeen() {
    int dt = 1;
    int iterations = 100;
    int n_threads = 1;
    int n_bodies = 100;

    struct body* bodies = 
        generate_random_bodies(
			2.5e11, 5.0e4, 2.0e30, 2.0e23, n_bodies
		);

    double initial_energy = energy(bodies, n_bodies, n_threads);

    for (int i = 0; i < iterations; i++) {
        step(bodies, n_bodies, dt, n_threads);
    }

    double final_energy = energy(bodies, n_bodies, n_threads);
    double delta = fabs(0.01 * initial_energy);

    assert_true(fabs(initial_energy - final_energy) < delta);

    free(bodies);
}

void step_test_eighteen() {
    int dt = 1;
    int iterations = 1000;
    int n_threads = 1;
    int n_bodies = 100;

	struct body* bodies = 
        generate_random_bodies(
			2.5e11, 5.0e4, 2.0e30, 2.0e23, n_bodies
		);

    double initial_energy = energy(bodies, n_bodies, n_threads);

    for (int i = 0; i < iterations; i++) {
        step(bodies, n_bodies, dt, n_threads);
    }

    double final_energy = energy(bodies, n_bodies, n_threads);
    double delta = fabs(0.01 * initial_energy);

    assert_true(fabs(initial_energy - final_energy) < delta);

    free(bodies);
}

void step_test_nineteen() {
    int dt = 0.01;
    int iterations = 10;
    int n_threads = 4;
    int n_bodies = 100;

    struct body* bodies = 
        generate_random_bodies(
			2.5e11, 5.0e4, 2.0e30, 2.0e23, n_bodies
		);

    double initial_energy = energy(bodies, n_bodies, n_threads);

    for (int i = 0; i < iterations; i++) {
        step(bodies, n_bodies, dt, n_threads);
    }

    double final_energy = energy(bodies, n_bodies, n_threads);
    double delta = fabs(0.01 * initial_energy);

    assert_true(fabs(initial_energy - final_energy) < delta);

    free(bodies);
}

void step_test_twenty() {
    int dt = 0.01;
    int iterations = 100;
    int n_threads = 4;
    int n_bodies = 100;

    struct body* bodies = 
        generate_random_bodies(
			2.5e11, 5.0e4, 2.0e30, 2.0e23, n_bodies
		);

    double initial_energy = energy(bodies, n_bodies, n_threads);

    for (int i = 0; i < iterations; i++) {
        step(bodies, n_bodies, dt, n_threads);
    }

    double final_energy = energy(bodies, n_bodies, n_threads);
    double delta = fabs(0.01 * initial_energy);

    assert_true(fabs(initial_energy - final_energy) < delta);

    free(bodies);
}

void step_test_twenty_one() {
    int dt = 0.01;
    int iterations = 1000;
    int n_threads = 4;
    int n_bodies = 100;

    struct body* bodies = 
        generate_random_bodies(
			2.5e11, 5.0e4, 2.0e30, 2.0e23, n_bodies
		);

    double initial_energy = energy(bodies, n_bodies, n_threads);

    for (int i = 0; i < iterations; i++) {
        step(bodies, n_bodies, dt, n_threads);
    }

    double final_energy = energy(bodies, n_bodies, n_threads);
    double delta = fabs(0.01 * initial_energy);

    assert_true(fabs(initial_energy - final_energy) < delta);

    free(bodies);
}

void step_test_twenty_two() {
    int dt = 1;
    int iterations = 10;
    int n_threads = 4;
    int n_bodies = 100;

    struct body* bodies = 
        generate_random_bodies(
			2.5e11, 5.0e4, 2.0e30, 2.0e23, n_bodies
		);

    double initial_energy = energy(bodies, n_bodies, n_threads);

    for (int i = 0; i < iterations; i++) {
        step(bodies, n_bodies, dt, n_threads);
    }

    double final_energy = energy(bodies, n_bodies, n_threads);
    double delta = fabs(0.01 * initial_energy);

    assert_true(fabs(initial_energy - final_energy) < delta);

    free(bodies);
}

void step_test_twenty_three() {
    int dt = 1;
    int iterations = 100;
    int n_threads = 4;
    int n_bodies = 100;

    struct body* bodies = 
        generate_random_bodies(
			2.5e11, 5.0e4, 2.0e30, 2.0e23, n_bodies
		);

    double initial_energy = energy(bodies, n_bodies, n_threads);

    for (int i = 0; i < iterations; i++) {
        step(bodies, n_bodies, dt, n_threads);
    }

    double final_energy = energy(bodies, n_bodies, n_threads);
    double delta = fabs(0.01 * initial_energy);

    assert_true(fabs(initial_energy - final_energy) < delta);

    free(bodies);
}

void step_test_twenty_four() {
    int dt = 1;
    int iterations = 1000;
    int n_threads = 4;
    int n_bodies = 100;

	struct body* bodies = 
        generate_random_bodies(
			2.5e11, 5.0e4, 2.0e30, 2.0e23, n_bodies
		);

    double initial_energy = energy(bodies, n_bodies, n_threads);

    for (int i = 0; i < iterations; i++) {
        step(bodies, n_bodies, dt, n_threads);
    }

    double final_energy = energy(bodies, n_bodies, n_threads);
    double delta = fabs(0.01 * initial_energy);

    assert_true(fabs(initial_energy - final_energy) < delta);

    free(bodies);
}

void step_test_twenty_five() {
    int dt = 0.01;
    int iterations = 10;
    int n_threads = 16;
    int n_bodies = 100;

    struct body* bodies = 
        generate_random_bodies(
			2.5e11, 5.0e4, 2.0e30, 2.0e23, n_bodies
		);

    double initial_energy = energy(bodies, n_bodies, n_threads);

    for (int i = 0; i < iterations; i++) {
        step(bodies, n_bodies, dt, n_threads);
    }

    double final_energy = energy(bodies, n_bodies, n_threads);
    double delta = fabs(0.01 * initial_energy);

    assert_true(fabs(initial_energy - final_energy) < delta);

    free(bodies);
}

void step_test_twenty_six() {
    int dt = 0.01;
    int iterations = 100;
    int n_threads = 16;
    int n_bodies = 100;

    struct body* bodies = 
        generate_random_bodies(
			2.5e11, 5.0e4, 2.0e30, 2.0e23, n_bodies
		);

    double initial_energy = energy(bodies, n_bodies, n_threads);

    for (int i = 0; i < iterations; i++) {
        step(bodies, n_bodies, dt, n_threads);
    }

    double final_energy = energy(bodies, n_bodies, n_threads);
    double delta = fabs(0.01 * initial_energy);

    assert_true(fabs(initial_energy - final_energy) < delta);

    free(bodies);
}

void step_test_twenty_seven() {
    int dt = 0.01;
    int iterations = 1000;
    int n_threads = 16;
    int n_bodies = 100;

    struct body* bodies = 
        generate_random_bodies(
			2.5e11, 5.0e4, 2.0e30, 2.0e23, n_bodies
		);

    double initial_energy = energy(bodies, n_bodies, n_threads);

    for (int i = 0; i < iterations; i++) {
        step(bodies, n_bodies, dt, n_threads);
    }

    double final_energy = energy(bodies, n_bodies, n_threads);
    double delta = fabs(0.01 * initial_energy);

    assert_true(fabs(initial_energy - final_energy) < delta);

    free(bodies);
}

void step_test_twenty_eight() {
    int dt = 1;
    int iterations = 10;
    int n_threads = 16;
    int n_bodies = 100;

    struct body* bodies = 
        generate_random_bodies(
			2.5e11, 5.0e4, 2.0e30, 2.0e23, n_bodies
		);

    double initial_energy = energy(bodies, n_bodies, n_threads);

    for (int i = 0; i < iterations; i++) {
        step(bodies, n_bodies, dt, n_threads);
    }

    double final_energy = energy(bodies, n_bodies, n_threads);
    double delta = fabs(0.01 * initial_energy);

    assert_true(fabs(initial_energy - final_energy) < delta);

    free(bodies);
}

void step_test_twenty_nine() {
    int dt = 1;
    int iterations = 100;
    int n_threads = 16;
    int n_bodies = 100;

    struct body* bodies = 
        generate_random_bodies(
			2.5e11, 5.0e4, 2.0e30, 2.0e23, n_bodies
		);

    double initial_energy = energy(bodies, n_bodies, n_threads);

    for (int i = 0; i < iterations; i++) {
        step(bodies, n_bodies, dt, n_threads);
    }

    double final_energy = energy(bodies, n_bodies, n_threads);
    double delta = fabs(0.01 * initial_energy);

    assert_true(fabs(initial_energy - final_energy) < delta);

    free(bodies);
}

void step_test_thirty() {
    int dt = 1;
    int iterations = 1000;
    int n_threads = 16;
    int n_bodies = 100;

	struct body* bodies = 
        generate_random_bodies(
			2.5e11, 5.0e4, 2.0e30, 2.0e23, n_bodies
		);

    double initial_energy = energy(bodies, n_bodies, n_threads);

    for (int i = 0; i < iterations; i++) {
        step(bodies, n_bodies, dt, n_threads);
    }

    double final_energy = energy(bodies, n_bodies, n_threads);
    double delta = fabs(0.01 * initial_energy);

    assert_true(fabs(initial_energy - final_energy) < delta);

    free(bodies);
}

void step_test_thirty_one() {
    int dt = 0.01;
    int iterations = 10;
    int n_threads = 32;
    int n_bodies = 100;

    struct body* bodies = 
        generate_random_bodies(
			2.5e11, 5.0e4, 2.0e30, 2.0e23, n_bodies
		);

    double initial_energy = energy(bodies, n_bodies, n_threads);

    for (int i = 0; i < iterations; i++) {
        step(bodies, n_bodies, dt, n_threads);
    }

    double final_energy = energy(bodies, n_bodies, n_threads);
    double delta = fabs(0.01 * initial_energy);

    assert_true(fabs(initial_energy - final_energy) < delta);

    free(bodies);
}

void step_test_thirty_two() {
    int dt = 0.01;
    int iterations = 100;
    int n_threads = 32;
    int n_bodies = 100;

    struct body* bodies = 
        generate_random_bodies(
			2.5e11, 5.0e4, 2.0e30, 2.0e23, n_bodies
		);

    double initial_energy = energy(bodies, n_bodies, n_threads);

    for (int i = 0; i < iterations; i++) {
        step(bodies, n_bodies, dt, n_threads);
    }

    double final_energy = energy(bodies, n_bodies, n_threads);
    double delta = fabs(0.01 * initial_energy);

    assert_true(fabs(initial_energy - final_energy) < delta);

    free(bodies);
}

void step_test_thirty_three() {
    int dt = 0.01;
    int iterations = 1000;
    int n_threads = 32;
    int n_bodies = 100;

    struct body* bodies = 
        generate_random_bodies(
			2.5e11, 5.0e4, 2.0e30, 2.0e23, n_bodies
		);

    double initial_energy = energy(bodies, n_bodies, n_threads);

    for (int i = 0; i < iterations; i++) {
        step(bodies, n_bodies, dt, n_threads);
    }

    double final_energy = energy(bodies, n_bodies, n_threads);
    double delta = fabs(0.01 * initial_energy);

    assert_true(fabs(initial_energy - final_energy) < delta);

    free(bodies);
}

void step_test_thirty_four() {
    int dt = 1;
    int iterations = 10;
    int n_threads = 32;
    int n_bodies = 100;

    struct body* bodies = 
        generate_random_bodies(
			2.5e11, 5.0e4, 2.0e30, 2.0e23, n_bodies
		);

    double initial_energy = energy(bodies, n_bodies, n_threads);

    for (int i = 0; i < iterations; i++) {
        step(bodies, n_bodies, dt, n_threads);
    }

    double final_energy = energy(bodies, n_bodies, n_threads);
    double delta = fabs(0.01 * initial_energy);

    assert_true(fabs(initial_energy - final_energy) < delta);

    free(bodies);
}

void step_test_thirty_five() {
    int dt = 1;
    int iterations = 100;
    int n_threads = 32;
    int n_bodies = 100;

    struct body* bodies = 
        generate_random_bodies(
			2.5e11, 5.0e4, 2.0e30, 2.0e23, n_bodies
		);

    double initial_energy = energy(bodies, n_bodies, n_threads);

    for (int i = 0; i < iterations; i++) {
        step(bodies, n_bodies, dt, n_threads);
    }

    double final_energy = energy(bodies, n_bodies, n_threads);
    double delta = fabs(0.01 * initial_energy);

    assert_true(fabs(initial_energy - final_energy) < delta);

    free(bodies);
}

void step_test_thirty_six() {

    int dt = 1;
    int iterations = 1000;
    int n_threads = 32;
    int n_bodies = 100;

	struct body* bodies = 
        generate_random_bodies(
			2.5e11, 5.0e4, 2.0e30, 2.0e23, n_bodies
		);

    double initial_energy = energy(bodies, n_bodies, n_threads);

    for (int i = 0; i < iterations; i++) {
        step(bodies, n_bodies, dt, n_threads);
    }

    double final_energy = energy(bodies, n_bodies, n_threads);
    double delta = fabs(0.01 * initial_energy);

    assert_true(fabs(initial_energy - final_energy) < delta);

    free(bodies);
}

void step_test_thirty_seven() {

    int dt = 100;
    int iterations = 10;
    int n_threads = 32;
    int n_bodies = 1000;

	struct body* bodies = 
        generate_random_bodies(
			2.5e11, 5.0e4, 2.0e30, 2.0e23, n_bodies
		);

    double initial_energy = energy(bodies, n_bodies, n_threads);

    for (int i = 0; i < iterations; i++) {
        step(bodies, n_bodies, dt, n_threads);
    }

    double final_energy = energy(bodies, n_bodies, n_threads);
    double delta = fabs(0.01 * initial_energy);

    assert_true(fabs(initial_energy - final_energy) < delta);

    free(bodies);
}

void step_test_thirty_eight() {
    int dt = 1000;
    int iterations = 10;
    int n_threads = 32;
    int n_bodies = 1000;

	struct body* bodies = 
        generate_random_bodies(
			2.5e11, 5.0e4, 2.0e30, 2.0e23, n_bodies
		);

    double initial_energy = energy(bodies, n_bodies, n_threads);

    for (int i = 0; i < iterations; i++) {
        step(bodies, n_bodies, dt, n_threads);
    }

    double final_energy = energy(bodies, n_bodies, n_threads);
    double delta = fabs(0.01 * initial_energy);

    assert_true(fabs(initial_energy - final_energy) < delta);

    free(bodies);
}

int main(void) {
    const struct CMUnitTest tests[] = {
        cmocka_unit_test(distance_test_one),
        cmocka_unit_test(distance_test_two),
        cmocka_unit_test(distance_test_three),
        cmocka_unit_test(distance_test_four),
        cmocka_unit_test(distance_test_five),
        cmocka_unit_test(distance_test_six),
        cmocka_unit_test(magnitude_test_one),
        cmocka_unit_test(magnitude_test_two),
        cmocka_unit_test(magnitude_test_three),
        cmocka_unit_test(magnitude_test_four),
        cmocka_unit_test(magnitude_test_five),
        cmocka_unit_test(unit_vector_test_one),
        cmocka_unit_test(unit_vector_test_two),
        cmocka_unit_test(unit_vector_test_three),
        cmocka_unit_test(unit_vector_test_four),
        cmocka_unit_test(unit_vector_test_five),
        cmocka_unit_test(acceleration_vector_test_one),
        cmocka_unit_test(acceleration_vector_test_two),
        cmocka_unit_test(acceleration_vector_test_three),
        cmocka_unit_test(acceleration_vector_test_four),
        cmocka_unit_test(acceleration_vector_test_five),
        cmocka_unit_test(extract_bodies_from_file_test),
        cmocka_unit_test(step_test_one),
        cmocka_unit_test(step_test_two),
        cmocka_unit_test(step_test_three),
        cmocka_unit_test(step_test_four),
        cmocka_unit_test(step_test_five),
        cmocka_unit_test(step_test_six),
        cmocka_unit_test(step_test_seven),
        cmocka_unit_test(step_test_eight),
        cmocka_unit_test(step_test_nine),
        cmocka_unit_test(step_test_ten),
        cmocka_unit_test(step_test_eleven),
        cmocka_unit_test(step_test_twelve),
        cmocka_unit_test(step_test_thirteen),
        cmocka_unit_test(step_test_fourteen),
        cmocka_unit_test(step_test_fifteen),
        cmocka_unit_test(step_test_sixteen),
        cmocka_unit_test(step_test_seventeen),
        cmocka_unit_test(step_test_eighteen),
        cmocka_unit_test(step_test_nineteen),
        cmocka_unit_test(step_test_twenty),
        cmocka_unit_test(step_test_twenty_one),
        cmocka_unit_test(step_test_twenty_two),
        cmocka_unit_test(step_test_twenty_three),
        cmocka_unit_test(step_test_twenty_four),
        cmocka_unit_test(step_test_twenty_five),
        cmocka_unit_test(step_test_twenty_six),
        cmocka_unit_test(step_test_twenty_seven),
        cmocka_unit_test(step_test_twenty_eight),
        cmocka_unit_test(step_test_twenty_nine),
        cmocka_unit_test(step_test_thirty),
        cmocka_unit_test(step_test_thirty_one),
        cmocka_unit_test(step_test_thirty_two),
        cmocka_unit_test(step_test_thirty_three),
        cmocka_unit_test(step_test_thirty_four),
        cmocka_unit_test(step_test_thirty_five),
        cmocka_unit_test(step_test_thirty_six),
        cmocka_unit_test(step_test_thirty_seven),
        cmocka_unit_test(step_test_thirty_eight),
          
    };
    return cmocka_run_group_tests(tests, NULL, NULL);
}