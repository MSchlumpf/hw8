/* Homework 8
 * Solving a Kepler problem
 * 
 * Written by Michael Stumpf
 */

#include<cmath>
#include<fstream>
using namespace std;

//-----------------
//declaration of sub-functions
void step(double* const p, double* const q, const double& dt, double& H, const double& dim);
//-----------------
//main function
int main(){
  const int dim = 2;
  double q[dim], p[dim];
  const double pi = 3.14159265359, e = 0.6;
  const double tStart = 0, tEnd = 20*pi;
  const double dt = 0.0005;
  const int N = (tEnd-tStart)/dt;
  double t = tStart;
  //-----------------
  //starting values
  q[0] = 1-e;
  q[1] = 0;
  p[0] = 0;
  p[1] = sqrt((1+e)/(1-e));
  double H = 0.5*(p[0]*p[0]+p[1]*p[1])-1.0/sqrt(q[0]*q[0]+q[1]*q[1]);
  //-----------------
  ofstream out("data.txt");
  out << t << "\t" << q[0] << "\t" << q[1] << "\t" << p[0] << "\t" << p[1] << "\t" << H << endl;
  //-----------------
  //start numerical solving
  for(int i=0; i<N; i++){
    step(p, q, dt, H, dim);
    t += dt;
    out << t << "\t" << q[0] << "\t" << q[1] << "\t" << p[0] << "\t" << p[1] << "\t" << H << endl;
  }
  out.close();
  return 0;
}

//-----------------
//implementation of step-function
void step(double* const p, double* const q, const double& dt, double& H, const double& dim){
  double r = q[0]*q[0]+q[1]*q[1];
  for(int i=0; i<dim; i++){
    p[i] -= dt*q[i]/pow(r, 3.0/2.0);
    q[i] += dt*p[i];
  }
  H = 0.5*(p[0]*p[0]+p[1]*p[1])-1.0/sqrt(q[0]*q[0]+q[1]*q[1]);
}