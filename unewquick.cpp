#include <iostream>
# include <cmath>

using namespace std;
const int N=102;
const double Gr = 1e4;

//corrector step

void unewquick( double unew[N][N],  double u[N][N], double p[N][N], double fe[N][N],double fw[N][N],double fn[N][N],double fs[N][N], float dx, float dy, float volp, float dt, int Nreal)
{
   int i,j,k, gs=0;
    int arr[]={1,N-2};
     float Rp, ap, rms;
     double ffe,ffw,ffn,ffs,ff,tf,ue,uw,un,us;

        do
        {
        rms=0;
        //interior cells quick scheme
   for( i=2;i<N-2;++i)
    {
    for(j=2;j<N-2;++j)
        {
Rp = (volp/dt)*(u[i][j]-unew[i][j])+(volp/sqrt(Gr))*((unew[i+1][j]-2*unew[i][j]+unew[i-1][j])/(dx*dx)+(unew[i][j+1]-2*unew[i][j]+unew[i][j-1])/(dy*dy)) - (p[i+1][j]*dy-p[i-1][j]*dy)/2.0 ;

if (fe[i][j]>=0 )
   {
     ffe=0.75*fe[i][j]; ue = (0.75)*unew[i][j]+(0.375)*unew[i+1][j]-(0.125)*unew[i-1][j];
   }
 else
     {
    ffe=0.375*fe[i][j];  ue = (0.75)*unew[i+1][j]+(0.375)*unew[i][j]-(0.125)*unew[i+2][j];
     }
 if (fw[i][j]>=0 )
    {
         ffw=0.75*fw[i][j]; uw = (0.75)*unew[i][j]+(0.375)*unew[i-1][j]-(0.125)*unew[i+1][j];
    }
 else
    {
         ffw=0.375*fw[i][j];  uw = (0.75)*unew[i-1][j]+(0.375)*unew[i][j]-(0.125)*unew[i-2][j];
    }
 if (fn[i][j]>=0 )
   {
      ffn=0.75*fn[i][j]; un = (0.75)*unew[i][j]+(0.375)*unew[i][j+1]-(0.125)*unew[i][j-1];
   }
 else
    {
       ffn=0.375*fn[i][j]; un = (0.75)*unew[i][j+1]+(0.375)*unew[i][j]-(0.125)*unew[i][j+2];
    }
 if (fs[i][j]>=0 )
    {
       ffs=0.75*fs[i][j];  us = (0.75)*unew[i][j]+(0.375)*unew[i][j-1]-(0.125)*unew[i][j+1];
    }
 else
  {
      ffs=0.375*fs[i][j];  us = (0.75)*unew[i][j-1]+(0.375)*unew[i][j]-(0.125)*unew[i][j-2];
  }
     tf=fe[i][j]*ue + fw[i][j]*uw + fn[i][j]*un + fs[i][j]*us;
     ff=ffe+ffw+ffn+ffs;
    Rp = Rp - tf;
ap = (volp/dt) + (2.0*volp/sqrt(Gr))*(1.0/(dx*dx)+1.0/(dy*dy)) + ff;
unew[i][j]=unew[i][j] + Rp/ap;
rms=rms+Rp*Rp;
}
    }
//top and bottom cells// upwind scheme
    for (k=0;k<2;++k){
    j=arr[k];
for (i=1;i<N-1;++i)
        {
Rp = (volp/dt)*(u[i][j]-unew[i][j])+(volp/sqrt(Gr))*((unew[i+1][j]-2*unew[i][j]+unew[i-1][j])/(dx*dx)+(unew[i][j+1]-2*unew[i][j]+unew[i][j-1])/(dy*dy)) - (p[i+1][j]*dy-p[i-1][j]*dy)/2.0 ;

if (fe[i][j]>=0 )
   {
     ffe=fe[i][j]; ue=unew[i][j];
   }
 else
     {
    ffe=0; ue=unew[i+1][j];
     }
 if (fw[i][j]>=0 )
    {
       ffw=fw[i][j]; uw=unew[i][j];
    }
 else
    {
        ffw=0;  uw=unew[i-1][j];
    }
 if (fn[i][j]>=0 )
   {
      ffn=fn[i][j]; un=unew[i][j];
   }
 else
    {
        ffn=0; un=unew[i][j+1];
    }
 if (fs[i][j]>=0 )
    {
      ffs=fs[i][j]; us=unew[i][j];
    }
 else
  {
      ffs=0; us=unew[i][j-1];
  }
     tf=fe[i][j]*ue + fw[i][j]*uw + fn[i][j]*un + fs[i][j]*us;
     ff=ffe+ffw+ffn+ffs;
    Rp = Rp - tf;
ap = (volp/dt) + (2.0*volp/sqrt(Gr))*(1.0/(dx*dx)+1.0/(dy*dy)) + ff;
unew[i][j]=unew[i][j] + Rp/ap;
rms=rms+Rp*Rp;
}
    }



//left and right cells  // upwind scheme
    for (k=0;k<2;++k){
    i=arr[k];
for (j=2;j<N-2;++j){
Rp = (volp/dt)*(u[i][j]-unew[i][j])+(volp/sqrt(Gr))*((unew[i+1][j]-2*unew[i][j]+unew[i-1][j])/(dx*dx)+(unew[i][j+1]-2*unew[i][j]+unew[i][j-1])/(dy*dy)) - (p[i+1][j]*dy-p[i-1][j]*dy)/2.0 ;

if (fe[i][j]>=0 )
   {
     ffe=fe[i][j]; ue=unew[i][j];
   }
 else
     {
    ffe=0; ue=unew[i+1][j];
     }
 if (fw[i][j]>=0 )
    {
       ffw=fw[i][j]; uw=unew[i][j];
    }
 else
    {
        ffw=0;  uw=unew[i-1][j];
    }
 if (fn[i][j]>=0 )
   {
      ffn=fn[i][j]; un=unew[i][j];
   }
 else
    {
        ffn=0; un=unew[i][j+1];
    }
 if (fs[i][j]>=0 )
    {
      ffs=fs[i][j]; us=unew[i][j];
    }
 else
  {
      ffs=0; us=unew[i][j-1];
  }
     tf=fe[i][j]*ue + fw[i][j]*uw + fn[i][j]*un + fs[i][j]*us;
     ff=ffe+ffw+ffn+ffs;
    Rp = Rp - tf;
ap = (volp/dt) + (2.0*volp/sqrt(Gr))*(1.0/(dx*dx)+1.0/(dy*dy)) + ff;
unew[i][j]=unew[i][j] + Rp/ap;
rms=rms+Rp*Rp;
}
    }

//boundary conditions ---- fictitious cells

    for (j=1;j<N-1;++j)    //left wall
        unew[0][j] = -unew[1][j];


    for (j=1;j<N-1;++j)      // right wall
        unew[N-1][j] = -unew[N-2][j];


    for (j=1;j<N-1;++j)      // bottom wall
        unew[j][0] = -unew[j][1];


    for (j=1;j<N-1;++j)    // top wall
        unew[j][N-1] = -unew[j][N-2];


        rms=sqrt(rms/Nreal);
   gs=gs+1;
        }
    while (rms >1e-6);
}




