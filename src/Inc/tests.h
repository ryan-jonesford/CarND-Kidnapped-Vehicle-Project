#include "helper_functions.h"
#include "particle_filter.h"
#include <iostream>

#ifndef TESTS_H
#define TESTS_H

class tests {
   public:
    tests();
    bool tester(void);

   private:
    double init_x_test;
    double init_y_test;
    double init_theta_test;
    double std_test[3];
};
#endif /* TESTS_H */