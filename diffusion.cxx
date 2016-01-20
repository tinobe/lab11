#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
//---------------------------------------
using namespace std;
//---------------------------------------
void writeToFile(const double* const u, const string s, const double dx,
                 const double xmin, const int N, const double t);
void initialize(double* const u0, const double dx,const double dt, const double xmin,
                const int N);
void step(double* const u1,  double* const u0,  const double dt,
          const double dx, const double D,const int N);

//---------------------------------------
int main(){

  const double tEnd = 20 ;
  const double D = 1;

  const int N  = 200;
  const double xmin = -20;
  const double xmax = 20;
  const double dx = (xmax-xmin)/(N-1) ;

  double dt = dx;
  double t = 0;
  const int Na = 10;
  const int Nk = int(tEnd/Na/dt);


  double* u0 = new double[N];
  double* u1 = new double[N];
  double* h;
  stringstream strm;

  initialize(u0,dx,dt, xmin,N);

  writeToFile(u0, "u_0", dx, xmin, N,t);

//   cout << "Nk = " << Nk << endl;

  for(int i=1; i<=Na; i++)
  {
   for(int j=0; j<Nk; j++){

      step(u1,u0,dt,dx,D,N);

      h = u0;
      u0 = u1;
      u1 = h;
      t +=dt;
   }
   strm.str("");
   strm << "u_" << i;
   writeToFile(u0, strm.str(), dx, xmin, N,t);
  }

//   cout << "t = " << t << endl;

  delete[] u0;
  delete[] u1;
  return 0;
}
//-----------------------------------------------
void step(double* const f1, double* const f0,
          const double dt, const double dx,
          const double D, const int N)
{

  double* d=new double[N];
  double* u=new double[N];
  double* l=new double[N];

  for(int i=0;i<N;i++) d[i] = 1.0 + 2.0*D*dt/(dx*dx);
  for(int i=0;i<N;i++) u[i] = - D*dt/(dx*dx);
  for(int i=0;i<N;i++) l[i] = - D*dt/(dx*dx);
  
  for(int i=1; i<N; i++){
    d[i]-=l[i]/d[i-1]*u[i-1];
    f0[i]-=l[i]/d[i-1]*f0[i-1];
    l[i]=0;
  }
  
  f1[N-1]=f0[N-1]/d[N-1];
  for(int i=N-2; i>0; i--) f1[i]=(f0[i]-u[i]*f1[i+1])/d[i];

  delete[] d;
  delete[] u;
  delete[] l;
}
//-----------------------------------------------
void initialize(double* const u0, const double dx,
                const double dt, const double xmin,  const int N)
{
   double u,ux, uxx;
   for(int i=0; i<N; i++)
   {
     double x = xmin + i*dx;
     u0[i] = 1.0/sqrt(4*M_PI)*exp(-x*x/4.0);

   }
}
//-----------------------------------------------
void writeToFile(const double* const u, const string s, const double dx,
                 const double xmin, const int N, const double t)
{
   ofstream out(s.c_str());
   for(int i=0; i<N; i++){
     double x = xmin + i * dx;
     double ana = 1.0/sqrt(4*M_PI*(t+1.0)) * exp(-x*x/(4*(t+1)));
     out << x << "\t" << u[i] << "\t" << ana << endl;
   }
   out.close();
   double xs = xmin + N/2*dx;
   cout << t << "\t" << u[N/2] * sqrt(4*M_PI*(t+1)) / exp(-xs*xs/(4*(t+1))) << "\t" << u[N/2] << "\t" << 1./ sqrt(4*M_PI*(t+1)) * exp(-xs*xs/(4*(t+1))) << endl;
}
