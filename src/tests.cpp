#include "Inc/tests.h"

using namespace std;

tests::tests() {
    init_x_test = 4983;
    init_y_test = 5029;
    init_theta_test = 1.201;
    std_test[0] = 2.0;
    std_test[1] = 2.0;
    std_test[2] = 0.5;
}

bool tests::tester() {
    ParticleFilter pf;

    // Test initalization function
    pf.init(init_x_test, init_y_test, init_theta_test, std_test);
    double sum_x = 0;
    double sum_y = 0;
    double sum_theta = 0;
    for (size_t i = 0; i < pf.particles.size(); ++i) {
        sum_x += pf.particles[i].x;
        sum_y += pf.particles[i].y;
        sum_theta += pf.particles[i].theta;
    }
    bool passed = true;
    if (init_x_test - std_test[0] < sum_x / 100 &&
        sum_x / 100 < init_x_test + std_test[0]) {
        passed = passed && true;
    } else {
        passed = false;
        cout << "Particle Filter init out of bounds on X init" << endl;
    }
    if (init_y_test - std_test[1] < sum_y / 100 &&
        sum_y / 100 < init_y_test + std_test[1]) {
        passed = passed && true;
    } else {
        passed = false;
        cout << "Particle Filter init out of bounds on Y init" << endl;
    }
    if (init_theta_test - std_test[2] < sum_theta / 100 &&
        sum_theta / 100 < init_theta_test + std_test[2]) {
        passed = passed && true;
    } else {
        passed = false;
        cout << "Particle Filter init out of bounds on Theta init" << endl;
    }
    // End of Testing initalization function

    // Test Prediction function
    for (vector<Particle>::iterator it = pf.particles.begin();
         it != pf.particles.end(); ++it) {
        it->x = 102;
        it->y = 65;
        it->theta = 5 * M_PI / 8;
    }
    pf.prediction(0.1, std_test, 110, M_PI / 8);
    sum_x = 0;
    sum_y = 0;
    sum_theta = 0;
    for (size_t i = 0; i < pf.particles.size(); ++i) {
        sum_x += pf.particles[i].x;
        sum_y += pf.particles[i].y;
        sum_theta += pf.particles[i].theta;
    }
    if (97.59 - std_test[0] < sum_x / 100 &&
        sum_x / 100 < 97.59 + std_test[0]) {
        passed = passed && true;
    } else {
        passed = false;
        cout << "Particle Filter prediction out of bounds on X prediction"
             << endl;
    }
    if (75.08 - std_test[1] < sum_y / 100 &&
        sum_y / 100 < 75.08 + std_test[1]) {
        passed = passed && true;
    } else {
        passed = false;
        cout << "Particle Filter prediction out of bounds on Y prediction"
             << endl;
    }
    if (51 * M_PI / 80 - std_test[2] < sum_theta / 100 &&
        sum_theta / 100 < 51 * M_PI / 80 + std_test[2]) {
        passed = passed && true;
    } else {
        passed = false;
        cout << "Particle Filter prediction out of bounds on Theta prediction"
             << endl;
    }

    return passed;
}