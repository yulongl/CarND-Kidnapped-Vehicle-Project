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
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>
#include <time.h> 

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	particles.clear();
	weights.clear();

	num_particles = 30;

	default_random_engine gen(time(0));

	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

	
	
	for (int i = 0; i < num_particles; i++) {

		Particle ptc;
		ptc.id = i;
		ptc.x = dist_x(gen);
		ptc.y = dist_y(gen);
		ptc.theta = dist_theta(gen);
		ptc.weight = 1.0;

		//weights.push_back(ptc.weight);
		particles.push_back(ptc);
		weights.push_back(ptc.weight);

	}
	

	is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/


	default_random_engine gen(time(0));

	normal_distribution<double> dist_x_noise(0, std_pos[0]);
	normal_distribution<double> dist_y_noise(0, std_pos[1]);
	normal_distribution<double> dist_theta_noise(0, std_pos[2]);

	// predict location of each particle at current timestamp using previous timestamp location
	// and control parameters
	if (fabs(yaw_rate) > 0.000001) {
		double v_d_yr = velocity / yaw_rate;
		double yr_m_dt = yaw_rate * delta_t;
	
		for (int i = 0; i < num_particles; i++) {
			// Not sure why we are provided the GPS measurement uncertainty here.
			// I think it's the control noise that should be added here, not GPS noise.
			particles[i].x = particles[i].x + v_d_yr * (sin(particles[i].theta + yr_m_dt) - sin(particles[i].theta)) + dist_x_noise(gen);
			particles[i].y = particles[i].y + v_d_yr * (cos(particles[i].theta) - cos(particles[i].theta + yr_m_dt)) + dist_y_noise(gen);
			particles[i].theta = particles[i].theta + yr_m_dt + dist_theta_noise(gen);
	
		}
	}
	else { 
		// condition of yaw_rate == 0
		double v_m_dt = velocity * delta_t;

		for (int i = 0; i < num_particles; i++) {
			
			particles[i].x = particles[i].x + v_m_dt * cos(particles[i].theta) + dist_x_noise(gen);
			particles[i].y = particles[i].y + v_m_dt * sin(particles[i].theta) + dist_y_noise(gen);
			//particles[i].theta = particles[i].theta + yr_m_dt;
	
		}
	}

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

	vector<double> weights_temp;

	// add land mark noise
	default_random_engine gen(time(0));
	vector<LandmarkObs> observations_temp;
	for (int i = 0; i < observations.size(); i++) {
		normal_distribution<double> dist_x(observations[i].x, std_landmark[0]);
		normal_distribution<double> dist_y(observations[i].y, std_landmark[1]);
		LandmarkObs obs_temp;
		obs_temp.x = dist_x(gen);
		obs_temp.y = dist_y(gen);

		observations_temp.push_back(obs_temp);

	}
	

	// calculate new weight for each particle
	for (int i = 0; i < num_particles; i++) {
		
		particles[i].weight = 1;

		// iterate through each observation
		for (int j = 0; j < observations_temp.size(); j++) {


			// convert to map coordinates
			double m_x = particles[i].x + (cos(particles[i].theta) * observations_temp[j].x) - (sin(particles[i].theta) * observations_temp[j].y);
			double m_y = particles[i].y + (sin(particles[i].theta) * observations_temp[j].x) + (cos(particles[i].theta) * observations_temp[j].y);
			
			if (map_landmarks.landmark_list.size() > 0) {

				// dataAssociation - find the closest map object
				double dist_min = dist(m_x, m_y, map_landmarks.landmark_list[0].x_f, map_landmarks.landmark_list[0].y_f);
				int dist_min_ind = 0;

				for (int j = 0; j < map_landmarks.landmark_list.size(); j++) {

					double temp = dist(m_x, m_y, map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f);
					if (temp <= dist_min) {
						dist_min = temp;
						dist_min_ind = j;
					}

				}

				// update weights
				double gauss_norm = (1.0 / (2.0 * M_PI * std_landmark[0] * std_landmark[1]));
				double exponent = pow((m_x - map_landmarks.landmark_list[dist_min_ind].x_f), 2.0) / (2.0 * pow(std_landmark[0], 2.0)) + pow((m_y - map_landmarks.landmark_list[dist_min_ind].y_f), 2.0) / (2.0 * pow(std_landmark[1], 2.0));
				particles[i].weight = particles[i].weight * gauss_norm * exp(-exponent);
			}
			
		}
		
		weights_temp.push_back(particles[i].weight);

	}

	weights = weights_temp;

	// check if the weights are all zero
	// if weights are all zero, an exception will be thrown out by discrete_distribution
	int eq_flag = 1;
	double sigma_pos[3] = { 0.3, 0.3, 0.01 };

	for (auto w = weights.begin(); w != (weights.end() - 1); ) {
		if (*w != *(++w))
			eq_flag = 0;
	}
	// if weights are all zero, which means non of the particle is reasonable,
	// then particles will be reinitiated.
	if ((eq_flag == 1) && (weights[0] < 0.00000001))
		init(particles[0].x, particles[0].y, particles[0].theta, sigma_pos);

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	vector<Particle> particles_temp;
	vector<double> weights_temp;

	default_random_engine gen(time(0));
	
	discrete_distribution<int> d(begin(weights), end(weights));

	// discrete_distribution will resample the particle vector by weights
	for (int n = 0; n < num_particles; n++) {
		particles_temp.push_back(particles[d(gen)]);
		//weights_temp.push_back(weights(d(gen)));
		
	}

	particles = particles_temp;



}

/*
Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}
*/

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
