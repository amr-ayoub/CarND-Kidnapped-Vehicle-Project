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

#include "particle_filter.h"
#include "helper_functions.h"

using namespace std;

//random engine
static default_random_engine gen;




void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	
	
	//Set the number of particles
	num_particles = 301;
	
	//random Gaussian noise
	normal_distribution<double> dist_x(0, std[0]);
	normal_distribution<double> dist_y(0, std[1]);
	normal_distribution<double> dist_theta(0, std[2]);
	
	//initializing particles and weights

		
		for (int particle = 0; particle < num_particles; particle++)
		{
			Particle p;
			
			p.id = particle;
			
			p.x = x ;
			p.y = y ;
			p.theta = theta;
			
			p.weight = 1.0;
			
			//add noise
			p.x = p.x + + dist_x(gen);
			p.y = p.y + + dist_y(gen);
			p.theta = p.theta + + dist_theta(gen);
			
			particles.push_back(p);
			
		}
		
		is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	
	//random Gaussian noise
	normal_distribution<double> dist_x(0, std_pos[0]);
	normal_distribution<double> dist_y(0, std_pos[1]);
	normal_distribution<double> dist_theta(0, std_pos[2]);

	for (int particle = 0; particle < num_particles; particle++)
	{
		//Calculate Prediction
		if (fabs(yaw_rate) < 0.0001)
		{
			particles[particle].x = particles[particle].x + velocity * delta_t + cos(particles[particle].theta);
			particles[particle].y = particles[particle].y + velocity * delta_t + sin(particles[particle].theta);
		}
		else
		{
			//Recall the equations for updating x, y and the yaw angle when the yaw rate is not equal to zero
			particles[particle].x = particles[particle].x + velocity / yaw_rate * (sin(particles[particle].theta + yaw_rate*delta_t) - sin(particles[particle].theta));
			particles[particle].y = particles[particle].y + velocity / yaw_rate * (cos(particles[particle].theta) - cos(particles[particle].theta + yaw_rate*delta_t));
			particles[particle].theta = particles[particle].theta + yaw_rate * delta_t;
		}
		
		// Add noise
		particles[particle].x = particles[particle].x + dist_x(gen);
		particles[particle].y = particles[particle].y + dist_y(gen);
		particles[particle].theta = particles[particle].theta + dist_theta(gen);
	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	
	for (unsigned int i = 0; i < observations.size(); i++)
	{
		//current observation
		LandmarkObs cur_obsv = observations[i];
		
		//landmark id from map 
		int map_id = 0;
		
		double min_dist = numeric_limits<double>::max();
		
		for (unsigned int j = 0; j < predicted.size(); j++)
		{
			//current prediction
			LandmarkObs cur_pred = predicted[j];
			
			// calculate current distance between current observation and predicted observation
			double cur_dist = dist(cur_obsv.x, cur_obsv.y, cur_pred.x, cur_pred.y);
			
			// find the prediction nearest the current observation
			if (cur_dist < min_dist)
			{
				min_dist = cur_dist;
				map_id = cur_pred.id;
			}
			
			
		}
		
		// set the observation to the nearest prediction id
		observations[i].id = map_id;
	}

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account 
	//   for the fact that the map's y-axis actually points downwards.)
	//   http://planning.cs.uiuc.edu/node99.html
	
	//loop over all particles
	for (int particle = 0; particle < num_particles; particle++)
	{
		//read the particle coordinates
		double p_x = particles[particle].x;
		double p_y = particles[particle].y;
    	double p_theta = particles[particle].theta;
		
		// landmarks within sensor range
		vector<LandmarkObs> predictions;
		
		//loop over landmarks
		for (unsigned int landmark = 0; landmark < map_landmarks.landmark_list.size(); landmark++)
		{
			//read the landmark coordinates
			double lm_x = map_landmarks.landmark_list[landmark].x_f;
      		double lm_y = map_landmarks.landmark_list[landmark].y_f;
      		int lm_id = map_landmarks.landmark_list[landmark].id_i;
			
			//find the landmarks which is near the sensor
			if ( fabs(lm_y - p_y) <= sensor_range  && fabs(lm_x - p_x) <= sensor_range)
				predictions.push_back(LandmarkObs{ lm_id, lm_x, lm_y });
			
			
		}
		
		
		//transform observation from vehicle coordinates to map coordinates
		vector<LandmarkObs> transformed_obsrv;
		for (unsigned int i = 0; i < observations.size(); i++)
		{
			double x = cos(p_theta)*observations[i].x - sin(p_theta)*observations[i].y + p_x;
			double y = sin(p_theta)*observations[i].x + cos(p_theta)*observations[i].y + p_y;
			transformed_obsrv.push_back(LandmarkObs{ observations[i].id, x, y });
		}
		
		
		//use dataAssociation for the predictions and transformed observations
		dataAssociation(predictions, transformed_obsrv);
		
		//init weight
    	particles[particle].weight = 1.0;
		
		for (unsigned int j = 0; j < transformed_obsrv.size(); j++)
		{
			//observation and prediction coordinates
			double obs_x = transformed_obsrv[j].x;
			double obs_y = transformed_obsrv[j].y;
			double pred_x;
			double pred_y;
			
			int prediction_ = transformed_obsrv[j].id;

			for (unsigned int m = 0; m < predictions.size(); m++)
			{
				if (predictions[m].id == prediction_)
				{
					pred_x = predictions[m].x;
					pred_y = predictions[m].y;
				}
			}
			
			//calculate the weight using multivariate Gaussian
			double std_x = std_landmark[0];
			double std_y = std_landmark[1];
			double weight_ = ( 1/(2*M_PI*std_x*std_y)) * exp( -( pow(pred_x-obs_x,2)/(2*pow(std_x, 2)) + (pow(pred_y-obs_y,2)/(2*pow(std_y, 2))) ) );
			
			particles[particle].weight = particles[particle].weight * weight_;
		}
		
		
	}
	
	
}


void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution.
	
	//creat new particles
	vector<Particle> new_particles;
	
	//read the current weights for the resampling wheel
	vector<double> weights;
	for (unsigned int particle = 0; particle < num_particles; particle++)
		weights.push_back(particles[particle].weight);
	
	//random start index for the resampling wheel
	uniform_int_distribution<int> uniintdist(0, num_particles-1);
	auto index = uniintdist(gen);
	
	//find the maximum weight
	double max_weight = *max_element(weights.begin(), weights.end());
	
	// get a uniform random distribution
	uniform_real_distribution<double> uni_dist(0.0, max_weight);
	
	//beta
	double beta = 0.0;
	
	//resampling wheel
	for (int i = 0; i < num_particles; i++)
	{
		beta = beta + uni_dist(gen) * 2.0;
		
		while(beta > weights[index])
		{
			beta = beta - weights[index];
			index = (index + 1) % num_particles;
		}
		
		new_particles.push_back(particles[index]);
	}
	
	//resample
	particles = new_particles;
	
}

void ParticleFilter::write(std::string filename) {
	// You don't need to modify this file.
	std::ofstream dataFile;
	dataFile.open(filename, std::ios::app);
	for (int i = 0; i < num_particles; ++i) {
		dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
	}
	dataFile.close();
}
