#include <iostream>
# include <cmath>
using namespace std;
const int N=102;

void pressure(double p[N][N], double up[N][N], double vp[N][N], float dx, float dy, int Nreal, float dt)
{
    int gs=0,i,j;
    float rms,r;
   do
 {
     rms=0;
   for( i=1;i<N-1;++i)
    {
    for(j=1;j<N-1;++j)
        {
  r = -dt*((p[i+1][j]-2.0*p[i][j]+p[i-1][j])/(dx*dx)+(p[i][j+1]-2.0*p[i][j]+p[i][j-1])/(dy*dy))+((up[i+1][j]-up[i-1][j])/(2.0*dx)+(vp[i][j+1]-vp[i][j-1])/(2.0*dy));
 p[i][j]=p[i][j] - (r)/((2.0*dt/(dx*dx))+(2.0*dt/(dy*dy)));    //real cells
 rms=rms+r*r;
        }
    }

    //fictitious cells boundary conditions

    for (j=1;j<N-1;++j)       //left side
    p[0][j] = p[1][j];


    for (j=1;j<N-1;++j)       // right side
    p[N-1][j] = p[N-2][j];


    for (j=1;j<N-1;++j)       //bottom side
    p[j][0] = p[j][1];



    for (j=1;j<N-1;++j)       // top side
    p[j][N-1] = p[j][N-2];


rms=sqrt(rms/Nreal);
gs=gs+1;
 }
while (rms > 1e-6);
}
