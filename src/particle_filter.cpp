/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <cmath> 
#include <vector>
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>
#include <random>

#include "particle_filter.h"

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	std::default_random_engine gen;
	std::normal_distribution<double> dist_x(x, std[0]);
	std::normal_distribution<double> dist_y(y, std[1]);
	std::normal_distribution<double> dist_theta(theta, std[2]);

	//set the number of particles
	num_particles = 50; 

	//create a temp vector of particles
	std::vector<Particle> n_particles;
	//create the particles
	for (unsigned int i=0; i < num_particles; i++) {
		Particle p = {};
		p.id       = i;
		//add noise
		p.x        = dist_x(gen);
		p.y        = dist_y(gen);
		p.theta    = dist_theta(gen);
		p.weight   = 1;
		n_particles.push_back(p);
	}
	particles = n_particles;
	is_initialized = true;
	std::cout << "ParticleFilter::init()" << std::endl;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	//create normal distributions (Gaussian) for x, y and theta
	std::default_random_engine gen;
	std::normal_distribution<double> dist_x(0, std_pos[0]);
	std::normal_distribution<double> dist_y(0, std_pos[1]);
	std::normal_distribution<double> dist_theta(0, std_pos[2]);

	//loop particles and access by reference
	for (Particle& p: particles) {
		//compute new position and theta
		double new_theta = 0.0;
		double new_x = 0.0;
		double new_y = 0.0;
		//check if yaw_rate is equal to zero
		if (std::abs(yaw_rate) == 0.0) {
			// moving straight
			new_x = p.x + (velocity*delta_t*cos(p.theta));
			new_y = p.y + (velocity*delta_t*sin(p.theta));
			new_theta = p.theta;
		}
		else {
			// turning
			new_theta = p.theta + yaw_rate*delta_t;
			new_x = p.x + (velocity/yaw_rate)*(sin(new_theta) - sin(p.theta));
			new_y = p.y + (velocity/yaw_rate)*(cos(p.theta) - cos(new_theta));      
		}
		//update particle
		p.x = new_x + dist_x(gen);
		p.y = new_y + dist_y(gen);
		p.theta = new_theta + dist_theta(gen); 
	}

	std::cout << "ParticleFilter::prediction()" << std::endl;
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y) {
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	particle.associations = associations;
	particle.sense_x = sense_x;
	particle.sense_y = sense_y;
}

std::string ParticleFilter::getAssociations(Particle best) {
	std::vector<int> v = best.associations;
	std::stringstream ss;
	std::copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
	std::string s = ss.str();
	s = s.substr(0, s.length()-1);  // get rid of the trailing space
	return s;
}

std::string ParticleFilter::getSenseX(Particle best) {
	std::vector<double> v = best.sense_x;
	std::stringstream ss;
	std::copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
	std::string s = ss.str();
	s = s.substr(0, s.length()-1);  // get rid of the trailing space
	return s;
}

std::string ParticleFilter::getSenseY(Particle best) {
	std::vector<double> v = best.sense_y;
	std::stringstream ss;
	std::copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
	std::string s = ss.str();
	s = s.substr(0, s.length()-1);  // get rid of the trailing space
	return s;
}
