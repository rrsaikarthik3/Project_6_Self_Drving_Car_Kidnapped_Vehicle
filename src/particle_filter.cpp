/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"
#ifndef M_PI
const double M_PI = 3.14159265358979323846;
#endif

using std::string;
using std::vector;

void ParticleFilter::init(double x, double y, double theta, double std[]) 
{
  /**
 
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   
   */
   std::default_random_engine gen;
   std::normal_distribution<double> dist_x(x, std[0]);
   std::normal_distribution<double> dist_y(y, std[1]);
   std::normal_distribution<double> dist_theta(theta, std[2]);
  
  num_particles = 10;  // TODO: Set the number of particles
  for (int i = 0; i < num_particles; i++)
  {
    Particle temp ; 
    temp.x = dist_x(gen);
    temp.y = dist_y(gen);
    temp.theta = dist_theta(gen);
    temp.weight = 1;
    particles.push_back(temp);
    
    
  }
  for (int i = 0; i < num_particles; i++)
  {
    //std::cout<<"\n"<<particles[i].x<<" "<<particles[i].y;
  }
  //std::cout<<" "<<particles.size();
  is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  //std::cout<<"\n Velocity "<<velocity<<" "<<yaw_rate;
  double temp = 0;
  if ( yaw_rate < 0.0000001 and  yaw_rate > -0.000001)
  {
    temp = 0.0;
  }
  else
  {
    temp = velocity/yaw_rate;
  }

  std::default_random_engine gen;
  for(int i = 0; i < num_particles; i++)
  {
    
   float temp_x,temp_y,temp_th;
    temp_x = particles[i].x;
    temp_y = particles[i].y;
    
    temp_th = particles[i].theta;
      
   temp_x =  temp_x + temp*(sin(temp_th + yaw_rate*delta_t) - sin(temp_th));
   temp_y =  temp_y + temp*(-cos(temp_th + yaw_rate*delta_t) + cos(temp_th));
   temp_th = temp_th + yaw_rate*delta_t;
      
   std::normal_distribution<double> dist_x(temp_x, std_pos[0]);
   std::normal_distribution<double> dist_y(temp_y, std_pos[1]);
   std::normal_distribution<double> dist_theta(temp_th, std_pos[2]);
    
   particles[i].x = dist_x(gen) ;
   particles[i].y = dist_y(gen) ;
   particles[i].theta = dist_theta(gen);
   
  }
  
  

}

vector<int> ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, 
                                     std::vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  int min_id;
  double min_dis,temp;
  vector<int> assoc;
  
  for(int i = 0; i < observations.size(); i++)
  {
    min_dis = 9999;
    for(int j = 0; j < predicted.size(); j++)
    {
       temp = dist(observations[i].x, observations[i].y , predicted[j].x ,predicted[j].y);
       if( temp < min_dis )
       {
         min_dis = temp;
         min_id = predicted[j].id;
         
       }
    }
    assoc.push_back(min_id);
    
  }
  return assoc;
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  

  
  double sum_weight = 0;
  std::vector<double> wt;
  
  for (int k = 0; k < num_particles; k++ )
  {
    std::vector<double> sense_x, sense_y;
    vector<LandmarkObs> predicted;
    LandmarkObs temp;
    vector<int> associations;
    wt.push_back(1);
     for( int j = 0; j < map_landmarks.landmark_list.size(); j++ )
     {
       
                    
        temp.x = (map_landmarks.landmark_list[j].x_f - particles[k].x)*cos(particles[k].theta) + (map_landmarks.landmark_list[j].y_f - particles[k].y)*sin(particles[k].theta);
        temp.y = (-map_landmarks.landmark_list[j].x_f + particles[k].x)*sin(particles[k].theta) + (map_landmarks.landmark_list[j].y_f - particles[k].y)*cos(particles[k].theta);
        temp.id = map_landmarks.landmark_list[j].id_i;
       
       if((dist(temp.x,temp.y,particles[k].x, particles[k].y) <= sensor_range) or true)
          {
            predicted.push_back(temp);
          }
        
      
     }
       
   associations = ParticleFilter::dataAssociation(predicted,observations);
   
     for( int i = 0; i< observations.size(); i++ )
      {
        sense_x.push_back(particles[k].x + observations[i].x*cos(particles[k].theta) - observations[i].y*sin(particles[k].theta));
        sense_y.push_back(particles[k].y + observations[i].x*sin(particles[k].theta) + observations[i].y*cos(particles[k].theta));                          
     }
    
    ParticleFilter::SetAssociations(particles[k],associations,sense_x,sense_y);
    
    
    for(int i = 0; i < particles[k].associations.size(); i++ )
    {
      //std::cout<<"\n Assoc "<<i<<" "<<particles[k].associations[i];
      for(int j = 0; j < map_landmarks.landmark_list.size(); j++)
      {
        if( (map_landmarks.landmark_list[j].id_i == particles[k].associations[i]))
        {
         
         long double temp_w = multiv_prob(std_landmark[0],std_landmark[1] , particles[k].sense_x[i], particles[k].sense_y[i],map_landmarks.landmark_list[j].x_f,map_landmarks.landmark_list[j].y_f);

      
            wt[k] = wt[k]*temp_w;

        } 
      }
    } 
    
    sum_weight = sum_weight + wt[k];
    
 }
  for (int k = 0; k < num_particles; k++ )
  {
    particles[k].weight = wt[k]/sum_weight;
    }
   
  
  
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  vector<double> v;
  std::discrete_distribution<int> distribution;
  for(int i = 0; i< num_particles; i++)
  {
    v.push_back(particles[i].weight);

  }
  
    std::random_device rd;
    std::mt19937 gen(rd());
    std::discrete_distribution<> d(v.begin(),v.end());
  
	std::vector<Particle> new_particles;
    for(int i = 0; i< num_particles; i++)
  {
    new_particles.push_back(particles[d(gen)]);

  }
  particles = new_particles;
  
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}
long double ParticleFilter::multiv_prob(double sig_x, double sig_y, double x_obs, double y_obs,
                   double mu_x, double mu_y) {
  // calculate normalization term
  double gauss_norm;

  gauss_norm = 1 / (2 * M_PI);

  
  gauss_norm = gauss_norm/pow(sig_x * sig_y,0.5);

  // calculate exponent
  double exponent = 0;

  
  exponent = pow((x_obs - mu_x), 2)/sig_x;
  
  exponent = exponent+pow((y_obs - mu_y), 2)/sig_y;    
  
  exponent = exponent + pow(y_obs - mu_y, 2)/sig_y;
  exponent = exponent/(2);

    
  // calculate weight using normalization terms and exponent
  long double weight;
  weight = gauss_norm * exp(-exponent);
    
  return weight;
}
