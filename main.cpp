//SEMI IMPLICIT FINITE VOLUME METHOD
//using Navier stokes solver + temperature equation
// for differentially heated cavity
//Both upwind and quick schemes
//---------------------------------------------------------------------------


#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;
const int N=102;   // Grid Size
//const double Re=100.0;
const double Gr = 1e4;  // Grashof's Number
const float pr = 1;     // Prandtl Number

#include"pressure.h"
#include"uvupwind.h"
#include"vupwind.h"
#include"uvnew.h"
#include"vvnew.h"
#include"uquick.h"
#include"vquick.h"
#include"unewquick.h"
#include"vnewquick.h"
#include"Energy.h"
#include"Energyquick.h"


int main()
{
    int i,j, n =0, Nreal, choice;
    Nreal=(N-2)*(N-2);
    float dt=0.001, dx=0.01, dy=0.01, volp, r, rms;   // dt time step dx, dy grid size
    volp=dx*dy;
    double  up[N][N] ={0}, vp[N][N]={0}, unew[N][N]={0}, vnew[N][N] ={0}, u[N][N]={0}, v[N][N]={0};
    double p[N][N]={0}, fe[N][N]={0}, fw[N][N]={0}, fn[N][N]={0}, fs[N][N]={0}, theta[N][N]={0}, thetanew[N][N]={0};
    double U[N],V[N];


    cout<<"Enter 1 for upwind scheme or 2 for quick scheme "<<endl;
    cin>>choice;

     do
    {
     rms=0;
 //predictor step for u and v---------------
     if (choice == 1)
       {
           uvupwind( u, up, fe, fw, fn, fs, dx, dy, volp, dt, Nreal);
           vupwind( v, vp, fe, fw, fn, fs, theta, dx, dy, volp, dt, Nreal);
       }
     else if (choice ==2)
        {
           uquick( u, up, fe, fw, fn, fs, dx, dy, volp, dt, Nreal);
           vquick( v, vp, fe, fw, fn, fs, theta, dx, dy, volp, dt, Nreal);
        }

 //pressure Poisson equation ----------

pressure(p, up, vp, dx, dy, Nreal, dt);

 //u and v corrector step----------------
    if (choice == 1)
    {
        uvnew(unew, u, p, fe, fw, fn, fs, dx, dy, volp, dt, Nreal);
        vvnew(vnew, v, p, fe, fw, fn, fs, theta, dx, dy, volp, dt, Nreal);
    }
    else if (choice == 2)
    {
        unewquick(unew, u, p, fe, fw, fn, fs, dx, dy, volp, dt, Nreal);
        vnewquick(vnew, v, p, fe, fw, fn, fs, theta, dx, dy, volp, dt, Nreal);
    }
//Flux calculation for each face --------
  for( i=1;i<N-1;++i)
    {
    for(j=1;j<N-1;++j)
        {
        fe[i][j] = (((up[i][j]+up[i+1][j])/2.0)-dt*(p[i+1][j]-p[i][j])/dx)*dy ;
        fw[i][j] = -(((up[i][j]+up[i-1][j])/2.0)-dt*(p[i][j]-p[i-1][j])/dx)*dy ;
        fn[i][j] = (((vp[i][j]+vp[i][j+1])/2.0)-dt*(p[i][j+1]-p[i][j])/dy)*dx ;
        fs[i][j] = -(((vp[i][j]+vp[i][j-1])/2.0)-dt*(p[i][j]-p[i][j-1])/dy)*dx ;
        }
        }

//Energy equation
    if (choice == 1)
   {
     Energy(thetanew, theta, unew, vnew, fe, fw, fn, fs, dx, dy, volp, dt, Nreal);
   }
    else if (choice == 2)
    {
       Energyquick(thetanew, theta, unew, vnew, fe, fw, fn, fs, dx, dy, volp, dt, Nreal);
    }

    for( i=0;i<N;++i)
    {for(j=0;j<N;++j)
    {
         r=pow((unew[i][j]-u[i][j]),2) + pow((vnew[i][j]-v[i][j]),2) +  pow((thetanew[i][j]-theta[i][j]),2);
       //r = pow((thetanew[i][j]-theta[i][j]),2);
        rms=rms+r;
        }}

rms = sqrt(rms/Nreal);

    for( i=0;i<N;++i)
    for(j=0;j<N;++j)
        {
            u[i][j] = unew[i][j];
            v[i][j] = vnew[i][j];
            theta[i][j] = thetanew[i][j];
        }

    n=n+1;
    cout<<"n: "<<n<<" rms: "<< rms<<endl;
    if(n>100000 || rms >10000){
	    cout <<" convergence Error, Maximum iterations" << n <<" reached"<<endl;
	    return 1;
    }
}
while ((rms/dt)>1e-3);         //steady state solution criteria
cout << "Solution Converged\n"; 

//plotting center line velocities

    for(j=0;j<N;++j)
    {
       U[j] = (u[50][j]+u[51][j])/2.0;
       V[j] = (v[j][50]+v[j][51])/2.0;
    }
    U[0]=(U[0]+U[1])/2.0;
    U[N-1]=(U[N-2]+U[N-1])/2.0;
    V[0]=(V[0]+V[1])/2.0;
    V[N-1]=(V[N-2]+V[N-1])/2.0;

//boundary points
        for(j=0;j<N;++j)
    {
       theta[0][j] = (theta[0][j]+theta[1][j])/2.0;
       theta[101][j] = (theta[100][j]+theta[101][j])/2.0;
       theta[j][0] = (theta[j][0]+theta[j][1])/2.0;
       theta[j][101] = (theta[j][100]+theta[j][101])/2.0;
    }

            for(j=0;j<N;++j)
    {
       u[0][j] = (u[0][j]+u[1][j])/2.0;
      u[101][j] = (u[100][j]+u[101][j])/2.0;
      u[j][0] = (u[j][0]+u[j][1])/2.0;
     u[j][101] = (u[j][100]+u[j][101])/2.0;
    }

           for(j=0;j<N;++j)
    {
       v[0][j] = (v[0][j]+v[1][j])/2.0;
      v[101][j] = (v[100][j]+v[101][j])/2.0;
      v[j][0] = (v[j][0]+v[j][1])/2.0;
     v[j][101] = (v[j][100]+v[j][101])/2.0;
    }

 //writing data to text file

  ofstream myfile;
  myfile.open ("u_v_data.txt");
//  myfile << "Writing U center vertical line velocity to a file.\n";
  for (j=0;j<N;++j)
        myfile << U[j] <<",";
  myfile << "\n";
  for(j=0;j<N;++j)
        myfile << V[j]<< ",";
    myfile.close();


   myfile.open("u_data.txt");
//    myfile << "\n Writing u velocity to file.\n";
    for (j=0;j<N;++j)
    {
      for(i=0;i<N;++i)
        myfile <<u[i][j] <<",";
        myfile <<"\n";
    }
   myfile.close();

       myfile.open("v_data.txt");
//    myfile << "\n Writing v velocity to file.\n";
    for (j=0;j<N;++j)
    {
      for(i=0;i<N;++i)
        myfile <<v[i][j] <<",";
         myfile <<"\n";
    }
   myfile.close();

          myfile.open("temperature.txt");
 //   myfile << "\n Writing temperature to file.\n";
    for (j=0;j<N;++j)
    {
    for(i=0;i<N;++i)
        {myfile <<theta[i][j] <<",";}
        myfile <<"\n";
    }
   myfile.close();

return 0;
    }


























