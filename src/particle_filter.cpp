/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <sstream>
#include <string>

#include "Inc/particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
    // Set the number of particles
    num_particles = 100;
    // Need this for sampling from distributions
    default_random_engine gen;
    // Standard deviations for x, y, and theta
    double std_x = std[0];
    double std_y = std[1];
    double std_theta = std[2];

    // Set standard deviations for x, y, and theta.
    std_x = 2;
    std_y = std_x;
    std_theta = 0.5;

    // Create a normal (Gaussian) distribution for x,y,theta.
    normal_distribution<double> dist_x(x, std_x);
    normal_distribution<double> dist_y(y, std_y);
    normal_distribution<double> dist_theta(theta, std_theta);

    Particle p;
    p.weight = 1;
    for (int inx = 0; inx < num_particles; ++inx) {
        p.id = inx;
        p.x = dist_x(gen);
        p.y = dist_y(gen);
        p.theta = dist_theta(gen);
        particles.push_back(p);
        // where "gen" is the random engine initialized earlier.
    }
    max_weight_ = 1;
    is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[],
                                double velocity, double yaw_rate) {
    // check for divide by zero error
    if (!yaw_rate) {
        // make it something really small
        yaw_rate = 1.0e-50;
    }
    double xf, yf, tf;
    // Need this for sampling from distributions
    default_random_engine gen;
    for (int inx = 0; inx < num_particles; ++inx) {
        xf = particles[inx].x;
        yf = particles[inx].y;
        tf = particles[inx].theta;
        xf =
            xf + velocity / yaw_rate * (sin(tf + yaw_rate * delta_t) - sin(tf));
        yf =
            yf + velocity / yaw_rate * (cos(tf) - cos(tf + yaw_rate * delta_t));
        tf = tf + yaw_rate * delta_t;
        // Add in gaussian Noise
        normal_distribution<double> dist_x(xf, std_pos[0]);
        normal_distribution<double> dist_y(yf, std_pos[1]);
        normal_distribution<double> dist_theta(tf, std_pos[2]);
        particles[inx].x = dist_x(gen);
        particles[inx].y = dist_y(gen);
        particles[inx].theta = dist_theta(gen);
    }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted,
                                     std::vector<LandmarkObs>& observations) {
    // TODO: Find the predicted measurement that is closest to each observed
    // measurement and assign the
    //   observed measurement to this particular landmark.
    // NOTE: this method will NOT be called by the grading code. But you will
    // probably find it useful to
    //   implement this method and use it as a helper during the updateWeights
    //   phase.
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   const std::vector<LandmarkObs>& observations,
                                   const Map& map_landmarks) {
    // TODO: Update the weights of each particle using a mult-variate Gaussian
    // distribution. You can read
    //   more about this distribution here:
    //   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
    // NOTE: The observations are given in the VEHICLE'S coordinate system. Your
    // particles are located
    //   according to the MAP'S coordinate system. You will need to transform
    //   between the two systems. Keep in mind that this transformation requires
    //   both rotation AND translation (but no scaling). The following is a good
    //   resource for the theory:
    //   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
    //   and the following is a good resource for the actual equation to
    //   implement (look at equation 3.33
    //   http://planning.cs.uiuc.edu/node99.html

    // transform each observation marker from the vehicle's coordinates to the
    // map's coordinates, with respect to our particle.
    vector<LandmarkObs> transformed_observations;
    LandmarkObs lo;
    for (int i = 0; i < num_particles; ++i) {
        double weight = 1.0;
        double prob = 1.0;
        for (size_t inx = 0; inx < observations.size(); ++inx) {
            array<double, 2> o_in_m_c = homogenous_transform(
                particles[i].x, particles[i].y, particles[i].theta,
                observations[inx].x, observations[inx].y);
            lo.id = inx;
            lo.x = o_in_m_c[0];
            lo.y = o_in_m_c[1];
            transformed_observations.push_back(lo);
            // find the associated landmark with the transformed observation
            double nearest;
            Map::single_landmark_s lm;
            bool nearest_initd = false;
            for (size_t iny = 0; iny < map_landmarks.landmark_list.size(); ++iny) {
                double mag =
                    dist(lo.x, lo.y, map_landmarks.landmark_list[iny].x_f,
                              map_landmarks.landmark_list[iny].y_f);
                if (mag > sensor_range) {
                    continue;
                }
                if (!nearest_initd) {
                    nearest = mag;
                    lm = map_landmarks.landmark_list[iny];
                }
                if (nearest != min(mag, nearest)) {
                    nearest = mag;
                    lm = map_landmarks.landmark_list[iny];
                }
                nearest_initd = true;
            }

            //update weight
            double obsx_m_near_2 = pow(
                (observations[inx].x - lm.x_f),
                2);
            double obsy_m_near_2 = pow(
                (observations[inx].y - lm.y_f),
                2);
            double stdx2 = std_landmark[0] * std_landmark[0];
            double stdy2 = std_landmark[1] * std_landmark[1];
            double gauss_norm =
                1.0 / (2.0 * M_PI * std_landmark[0] * std_landmark[1]);
            cout<<"this weight = "<<gauss_norm * exp(-1 * ((obsx_m_near_2 / (2.0 * stdx2)) +
                                             (obsy_m_near_2 / (2.0 * stdy2))))<<endl;
            weight *= gauss_norm * exp(-1 * ((obsx_m_near_2 / (2.0 * stdx2)) +
                                             (obsy_m_near_2 / (2.0 * stdy2))));
        }
        cout<<"weight = "<<weight<<endl;
        particles[i].weight = weight;
        max_weight_ = max(weight, max_weight_);
    }
}

void ParticleFilter::resample() {
    // TODO: Resample particles with replacement with probability proportional
    // to their weight. NOTE: You may find std::discrete_distribution helpful
    // here.
    //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    default_random_engine gen;
    vector<Particle> p;
    default_random_engine generator;
    for (int i = 0; i < num_particles; ++i) {
        discrete_distribution<int> weight_dist(0.0,2*max_weight_);
        double rand_weight = weight_dist(gen);
        uniform_int_distribution<int> index_dist(0,num_particles);
        int index = index_dist(gen);
        while(particles[index].weight < rand_weight ){
            rand_weight = rand_weight - particles[index].weight;
            index += 1;
            if(index > num_particles){
                index = 0;
            }
            p.push_back(particles[index]);
        }
    }
    particles = p;
}

Particle ParticleFilter::SetAssociations(Particle& particle,
                                         const std::vector<int>& associations,
                                         const std::vector<double>& sense_x,
                                         const std::vector<double>& sense_y) {
    // particle: the particle to assign each listed association, and
    // association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed
    // association sense_x: the associations x mapping already converted to
    // world coordinates sense_y: the associations y mapping already converted
    // to world coordinates

    particle.associations = associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
    return particle;
}

string ParticleFilter::getAssociations(Particle best) {
    vector<int> v = best.associations;
    stringstream ss;
    copy(v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length() - 1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best) {
    vector<double> v = best.sense_x;
    stringstream ss;
    copy(v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length() - 1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best) {
    vector<double> v = best.sense_y;
    stringstream ss;
    copy(v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length() - 1);  // get rid of the trailing space
    return s;
}

array<double, 2> ParticleFilter::homogenous_transform(double particle_x,
                                                      double particle_y,
                                                      double heading,
                                                      double map_x,
                                                      double map_y) {
    array<double, 2> transform;
    transform[0] = particle_x + (cos(heading) * map_x) - (sin(heading) * map_y);
    transform[1] = particle_y + (sin(heading) * map_x) + (cos(heading) * map_y);
    return transform;
}
