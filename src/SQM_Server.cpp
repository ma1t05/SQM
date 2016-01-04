
#include "SQM_Server.h"

server::server () {
  location = UNASIGNED_LOCATION;
  past_location = UNASIGNED_LOCATION;
  v = 1.0;
  beta = 2.0;
}

server::server (int i) {
  location = i;
  past_location = i;
  v = 1.0;
  beta = 2.0;
}

server::~server () {
  
}

void server::set_speed(double speed,double b) {
  v = speed;
  beta = b;
}

double server::get_speed () const {
  return v;
}

double server::get_beta () const {
  return beta;
}

double server::get_rate () const {
  return beta / v;
}

int server::get_location () const {
  return location;
}

int server::get_past_location () const {
  return past_location;
}

void server::set_location (int i) {
  past_location = location;
  location = i;
}

void server::test_location (int i) {
  if (location == UNASIGNED_LOCATION) past_location = i;
  location = i;
}
