#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <inttypes.h>
#include <libgen.h>
#include <math.h>

#include <omp.h>

#define pow2(x) ((x)*(x))
#define pow3(x) ((x)*(x)*(x))
#define piG4  8.38176e-10


#define DEBUG

typedef struct{
  double mass;
  double x;
  double y;
  double z;
  double vx;
  double vy;
  double vz;
  double forcex;
  double forcey;
  double forcez;
} Particle;

typedef struct{
  double density;
  double potential;
  double accelerationx;
  double accelerationy;
  double accelerationz;
} Grid_Point;





//----------------------------------------
// PROTOTYPES
//----------------------------------------
Particle        *particles_creator(int, int, double);
Grid_Point    ***grid_creator( int, double);
void             Mass_GridPoint(Particle *, Grid_Point ***, int, int, double);
void             Potential_GridPoint(Grid_Point ***, int, double);
void             Acceleration_GridPoint(Grid_Point ***, int, double);
void             Force_particles(Particle *, Grid_Point ***, int, int, double);
void             LeapFrog_particles(Particle *, int, double, double, int);

//----------------------------------------
// MAIN
//----------------------------------------
void main()
{
  Grid_Point ***grid;
  Particle     *particles;
  double        H, delta_t;
  int           N_part, N_gridpoints;
  FILE         *fp_tmp;
  int           m, j, i;
  double        deltax, deltay, deltaz, deltavx, deltavy, deltavz, delta_aux;
  
  //N_part       = 50;
  N_gridpoints = 50;
  H            = 3.0e10;
  delta_t      = 4064.0;
  
////////////////////////////////////////////////////
  //for a sphere of particles
  //Count number of particles
N_part =0;
  for (m=0;m<N_gridpoints;m++){
    for (i=0;i<N_gridpoints;i++){
      for (j=0;j<N_gridpoints;j++){
        if ((pow2(i-25) + pow2(j-25) + pow2(m-25))< 12.1*12.1 && (pow2(i-25) + pow2(j-25) + pow2(m-25))> 11.9*11.9)
          N_part += 1;
}
}
}

fprintf(stderr,"%d",N_part);
///////////////////////////////////////////////////
  
  particles = particles_creator(N_part, N_gridpoints, H);
  grid      = grid_creator( N_gridpoints, H);
  Mass_GridPoint(particles, grid, N_part, N_gridpoints, H);
  Potential_GridPoint(grid, N_gridpoints, H);
  Acceleration_GridPoint(grid, N_gridpoints, H);
  Force_particles(particles, grid, N_part, N_gridpoints, H);

  
   //LEAP-FROG
   //Initial conditions
  fp_tmp = fopen("datos2.dat", "w");

  for(m=0; m<N_part; m++)
  {
    particles[m].vx = particles[m].vx + delta_t/2.0*(particles[m].forcex/particles[m].mass); //vx_1/2
    particles[m].vy = particles[m].vy + delta_t/2.0*(particles[m].forcey/particles[m].mass); //vy_1/2
    particles[m].vz = particles[m].vz + delta_t/2.0*(particles[m].forcez/particles[m].mass); //vz_1/2
  }


  for(m=0; m<N_part; m++){fprintf(fp_tmp, "%f   %f  %f \n", particles[m].x, particles[m].y, particles[m].z);}

  //Leap-Frog iterations
  for (i=0; i<500; i++)
  {
    fprintf(stderr, "step: %d \n" ,i);
/*
    //For time step refinement!
    delta_aux=1e18;
    for(j=0; j<N_part; j++)
    {

      deltax  = delta_t*particles[j].vx;
      deltay  = delta_t*particles[j].vy;
      deltaz  = delta_t*particles[j].vz;
      deltavx = delta_t*(particles[j].forcex/particles[j].mass);
      deltavy = delta_t*(particles[j].forcey/particles[j].mass);
      deltavz = delta_t*(particles[j].forcez/particles[j].mass);
      if ((sqrt(pow2(deltax)+pow2(deltay)+pow2(deltaz))/sqrt(pow2(deltavx)+pow2(deltavy)+pow2(deltavz)))<delta_aux)
      {
        delta_aux=(sqrt(pow2(deltax)+pow2(deltay)+pow2(deltaz))/sqrt(pow2(deltavx)+pow2(deltavy)+pow2(deltavz)));
      }
 

    }
    delta_t = delta_aux;
    fprintf(stderr, "deltat %e \n",delta_aux);
*/
  
    LeapFrog_particles(particles, N_part, delta_t, H, N_gridpoints);
    Mass_GridPoint(particles, grid, N_part, N_gridpoints, H);
    Potential_GridPoint(grid, N_gridpoints, H);
    Acceleration_GridPoint(grid, N_gridpoints, H);
    Force_particles(particles, grid, N_part, N_gridpoints, H);

    for(m=0; m<N_part; m++){fprintf(fp_tmp, "%f  %f  %f \n", particles[m].x, particles[m].y, particles[m].z);}


  }
  fclose(fp_tmp);


  //Free memory 
  //careful here: *every* calloc() has to have a corresponding free()!
  for(i=0; i<N_gridpoints; i++) 
  {
    for(j=0; j<N_gridpoints; j++)
    {
      free(grid[i][j]);
    }
    free(grid[i]);
  }
  free(grid);
  free(particles);
   
   
  
  
}




//********************************************************************************************************************
//----------------------------------------
//grid_creator
//----------------------------------------
Grid_Point    ***grid_creator( int N_gridpoints, double H)
{
  int i, j;
  Grid_Point ***grid;
  grid = (Grid_Point ***)calloc(N_gridpoints, sizeof(Grid_Point **)); // note: this is a memory overhead!
  for(i=0; i<N_gridpoints; i++)
  {
    grid[i] = (Grid_Point **)calloc(N_gridpoints, sizeof(Grid_Point *));
    for(j=0; j<N_gridpoints; j++)
    {
      grid[i][j] = (Grid_Point *)calloc(N_gridpoints, sizeof(Grid_Point));
    }
  }
  
  return(grid);
}




//----------------------------------------
//particles_creator
//----------------------------------------
Particle        *particles_creator(int N_part, int N_gridpoints, double H)
{
  int i,m,j,counter;
  double max_dist; //maximum distance for each dimension. It will let a space un the countours for the periodic boundaries
  Particle *particles;
  
  srand(8974);     //seed for the random generator (rand)
  max_dist  = (N_gridpoints)*H;
  particles = (Particle *)calloc(N_part, sizeof(Particle));

//Initial condition for a sphere of particles
  counter =0;
  for (m=0;m<N_gridpoints;m++){
    for (i=0;i<N_gridpoints;i++){
      for (j=0;j<N_gridpoints;j++){
        if ((pow2(i-25) + pow2(j-25) + pow2(m-25))< 12.1*12.1 && (pow2(i-25) + pow2(j-25) + pow2(m-25))> 11.9*11.9){
          particles[counter].x = H*m;
          particles[counter].y = H*j;
          particles[counter].z = H*i;
          particles[counter].mass = 1.0e29;
          counter +=1;

}
          
}
}
}
fprintf(stderr,"%d \n",counter);
  return(particles);
}





//----------------------------------------
//Mass_GridPoint calculator ---> rho(g_klm)
//----------------------------------------
void Mass_GridPoint(Particle *particles, Grid_Point ***grid, int N_part, int N_gridpoints, double H)
{
  int    i,j,k;
  int    count_x, count_y, count_z;//it says if the particle is in the first or in the second half of the grid poitn interval
  double density, Wx1, Wy1, Wz1, Wx2, Wy2, Wz2;
  double nx, ny, nz; // (double) number of grid points
  int    mx, my, mz; // (int) number of grid points
  int    mxm1, mxp1, mym1, myp1, mzm1, mzp1;
  double dx, dy, dz; // distance between center of the grid point and the particle

  //make grid density equal to ZERO AGAIN
    for(i=0; i<=(N_gridpoints-1); i++) //for the "red" grid points
    {
      for(j=0; j<=(N_gridpoints-1); j++)
      {
        for(k=0; k<=(N_gridpoints-1); k++)
        {
          grid[i][j][k].density = 0.0;
        }
      }
    }

 //#pragma omp parallel for default(none) private(i) shared(N_part,N_gridpoints, particles, grid,H,nx,ny,nz,mx,my,mz,mxm1,mym1,mzm1,mxp1,myp1,mzp1,Wx1,Wx2,Wy1,Wy2,Wz1,Wz2,dx,dy,dz,count_x,count_y,count_z) 
  for(i=0; i<N_part; i++)
  {
    nx = particles[i].x/H;
    ny = particles[i].y/H;
    nz = particles[i].z/H;
    
    mx = (int) floor(nx);
    my = (int) floor(ny);
    mz = (int) floor(nz);
    mx = (int) floor(fmod((double)(mx),(double)(N_gridpoints)));
    my = (int) floor(fmod((double)(my),(double)(N_gridpoints)));
    mz = (int) floor(fmod((double)(mz),(double)(N_gridpoints)));
    
    // get indices of neighbouring nodes taking into account periodic boundaries
    mxm1 = (int) floor(fmod((double)(mx-1+N_gridpoints),(double)(N_gridpoints)));
    mym1 = (int) floor(fmod((double)(my-1+N_gridpoints),(double)(N_gridpoints)));
    mzm1 = (int) floor(fmod((double)(mz-1+N_gridpoints),(double)(N_gridpoints)));
    mxp1 = (int) floor(fmod((double)(mx+1),(double)(N_gridpoints)));
    myp1 = (int) floor(fmod((double)(my+1),(double)(N_gridpoints)));
    mzp1 = (int) floor(fmod((double)(mz+1),(double)(N_gridpoints)));
    
    //look for the half part of the grid point in which the particle is:
    dx = (nx - mx)*H;
    if( (nx - mx)<0.5 )
    {
      count_x = 0;
      Wx1 = 1.0 - ((dx + H/2.0)/H);
      Wx2 = 1.0 - ((H/2.0 - dx)/H);
    }
    else
    {
      count_x = 1;
      Wx1 = 1.0 - ((dx - H/2.0)/H);
      Wx2 = 1.0 - (((H-dx) + H/2.0)/H);
    }
    
    
    dy = (ny - my)*H;
    if( (ny - my)<0.5 )
    {
      count_y = 0;
      Wy1 = 1.0 - ((dy + H/2.0)/H);
      Wy2 = 1.0 - ((H/2.0 - dy)/H);
    }
    else
    {
      count_y = 1;
      Wy1 = 1.0 - ((dy - H/2.0)/H);
      Wy2 = 1.0 - (((H-dy) + H/2.0)/H);
    }
    
    
    dz = (nz - mz)*H;
    if( (nz - mz)<0.5 )
    {
      count_z = 0;
      Wz1 = 1.0 - ((dz + H/2.0)/H);
      Wz2 = 1.0 - ((H/2.0 - dz)/H);
    }
    else
    {
      count_z = 1;
      Wz1 = 1.0 - ((dz - H/2.0)/H);
      Wz2 = 1.0 - (((H-dz) + H/2.0)/H);
    }
    
    // assing the mass to the corresponnding grind points
    if (count_x==0 && count_y==0 && count_z==0)
    {
      grid[mx][my][mz].density       += (particles[i].mass*Wx2*Wy2*Wz2)/pow3(H);
      grid[mx][my][mzm1].density     += (particles[i].mass*Wx2*Wy2*Wz1)/pow3(H);
      grid[mx][mym1][mz].density     += (particles[i].mass*Wx2*Wy1*Wz2)/pow3(H);
      grid[mx][mym1][mzm1].density   += (particles[i].mass*Wx2*Wy1*Wz1)/pow3(H);
      grid[mxm1][my][mz].density     += (particles[i].mass*Wx1*Wy2*Wz2)/pow3(H);
      grid[mxm1][my][mzm1].density   += (particles[i].mass*Wx1*Wy2*Wz1)/pow3(H);
      grid[mxm1][mym1][mz].density   += (particles[i].mass*Wx1*Wy1*Wz2)/pow3(H);
      grid[mxm1][mym1][mzm1].density += (particles[i].mass*Wx1*Wy1*Wz1)/pow3(H);
    }
    else if (count_x==0 && count_y==0 && count_z==1)
    {
      grid[mx][my][mz].density       += (particles[i].mass*Wx2*Wy2*Wz1)/pow3(H);
      grid[mx][my][mzp1].density     += (particles[i].mass*Wx2*Wy2*Wz2)/pow3(H);
      grid[mx][mym1][mz].density     += (particles[i].mass*Wx2*Wy1*Wz1)/pow3(H);
      grid[mx][mym1][mzp1].density   += (particles[i].mass*Wx2*Wy1*Wz2)/pow3(H);
      grid[mxm1][my][mz].density     += (particles[i].mass*Wx1*Wy2*Wz1)/pow3(H);
      grid[mxm1][my][mzp1].density   += (particles[i].mass*Wx1*Wy2*Wz2)/pow3(H);
      grid[mxm1][mym1][mz].density   += (particles[i].mass*Wx1*Wy1*Wz1)/pow3(H);
      grid[mxm1][mym1][mzp1].density += (particles[i].mass*Wx1*Wy1*Wz2)/pow3(H);
    }
    else if (count_x==0 && count_y==1 && count_z==0)
    {
      grid[mx][my][mz].density       += (particles[i].mass*Wx2*Wy1*Wz2)/pow3(H);
      grid[mx][my][mzm1].density     += (particles[i].mass*Wx2*Wy1*Wz1)/pow3(H);
      grid[mx][myp1][mz].density     += (particles[i].mass*Wx2*Wy2*Wz2)/pow3(H);
      grid[mx][myp1][mzm1].density   += (particles[i].mass*Wx2*Wy2*Wz1)/pow3(H);
      grid[mxm1][my][mz].density     += (particles[i].mass*Wx1*Wy1*Wz2)/pow3(H);
      grid[mxm1][my][mzm1].density   += (particles[i].mass*Wx1*Wy1*Wz1)/pow3(H);
      grid[mxm1][myp1][mz].density   += (particles[i].mass*Wx1*Wy2*Wz2)/pow3(H);
      grid[mxm1][myp1][mzm1].density += (particles[i].mass*Wx1*Wy2*Wz1)/pow3(H);
    }
    else if (count_x==0 && count_y==1 && count_z==1)
    {
      grid[mx][my][mz].density       += (particles[i].mass*Wx2*Wy1*Wz1)/pow3(H);
      grid[mx][my][mzp1].density     += (particles[i].mass*Wx2*Wy1*Wz2)/pow3(H);
      grid[mx][myp1][mz].density     += (particles[i].mass*Wx2*Wy2*Wz1)/pow3(H);
      grid[mx][myp1][mzp1].density   += (particles[i].mass*Wx2*Wy2*Wz2)/pow3(H);
      grid[mxm1][my][mz].density     += (particles[i].mass*Wx1*Wy1*Wz1)/pow3(H);
      grid[mxm1][my][mzp1].density   += (particles[i].mass*Wx1*Wy1*Wz2)/pow3(H);
      grid[mxm1][myp1][mz].density   += (particles[i].mass*Wx1*Wy2*Wz1)/pow3(H);
      grid[mxm1][myp1][mzp1].density += (particles[i].mass*Wx1*Wy2*Wz2)/pow3(H);
    }
    
    else if (count_x==1 && count_y==0 && count_z==0)
    {
      grid[mx][my][mz].density       += (particles[i].mass*Wx1*Wy2*Wz2)/pow3(H);
      grid[mx][my][mzm1].density     += (particles[i].mass*Wx1*Wy2*Wz1)/pow3(H);
      grid[mx][mym1][mz].density     += (particles[i].mass*Wx1*Wy1*Wz2)/pow3(H);
      grid[mx][mym1][mzm1].density   += (particles[i].mass*Wx1*Wy1*Wz1)/pow3(H);
      grid[mxp1][my][mz].density     += (particles[i].mass*Wx2*Wy2*Wz2)/pow3(H);
      grid[mxp1][my][mzm1].density   += (particles[i].mass*Wx2*Wy2*Wz1)/pow3(H);
      grid[mxp1][mym1][mz].density   += (particles[i].mass*Wx2*Wy1*Wz2)/pow3(H);
      grid[mxp1][mym1][mzm1].density += (particles[i].mass*Wx2*Wy1*Wz1)/pow3(H);
    }
    else if (count_x==1 && count_y==0 && count_z==1)
    {
      grid[mx][my][mz].density       += (particles[i].mass*Wx1*Wy2*Wz1)/pow3(H);
      grid[mx][my][mzp1].density     += (particles[i].mass*Wx1*Wy2*Wz2)/pow3(H);
      grid[mx][mym1][mz].density     += (particles[i].mass*Wx1*Wy1*Wz1)/pow3(H);
      grid[mx][mym1][mzp1].density   += (particles[i].mass*Wx1*Wy1*Wz2)/pow3(H);
      grid[mxp1][my][mz].density     += (particles[i].mass*Wx2*Wy2*Wz1)/pow3(H);
      grid[mxp1][my][mzp1].density   += (particles[i].mass*Wx2*Wy2*Wz2)/pow3(H);
      grid[mxp1][mym1][mz].density   += (particles[i].mass*Wx2*Wy1*Wz1)/pow3(H);
      grid[mxp1][mym1][mzp1].density += (particles[i].mass*Wx2*Wy1*Wz2)/pow3(H);
    }
    else if (count_x==1 && count_y==1 && count_z==0)
    {
      grid[mx][my][mz].density       += (particles[i].mass*Wx1*Wy1*Wz2)/pow3(H);
      grid[mx][my][mzm1].density     += (particles[i].mass*Wx1*Wy1*Wz1)/pow3(H);
      grid[mx][myp1][mz].density     += (particles[i].mass*Wx1*Wy2*Wz2)/pow3(H);
      grid[mx][myp1][mzm1].density   += (particles[i].mass*Wx1*Wy2*Wz1)/pow3(H);
      grid[mxp1][my][mz].density     += (particles[i].mass*Wx2*Wy1*Wz2)/pow3(H);
      grid[mxp1][my][mzm1].density   += (particles[i].mass*Wx2*Wy1*Wz1)/pow3(H);
      grid[mxp1][myp1][mz].density   += (particles[i].mass*Wx2*Wy2*Wz2)/pow3(H);
      grid[mxp1][myp1][mzm1].density += (particles[i].mass*Wx2*Wy2*Wz1)/pow3(H);
    }
    else if (count_x==1 && count_y==1 && count_z==1)
    {
      grid[mx][my][mz].density       += (particles[i].mass*Wx1*Wy1*Wz1)/pow3(H);
      grid[mx][my][mzp1].density     += (particles[i].mass*Wx1*Wy1*Wz2)/pow3(H);
      grid[mx][myp1][mz].density     += (particles[i].mass*Wx1*Wy2*Wz1)/pow3(H);
      grid[mx][myp1][mzp1].density   += (particles[i].mass*Wx1*Wy2*Wz2)/pow3(H);
      grid[mxp1][my][mz].density     += (particles[i].mass*Wx2*Wy1*Wz1)/pow3(H);
      grid[mxp1][my][mzp1].density   += (particles[i].mass*Wx2*Wy1*Wz2)/pow3(H);
      grid[mxp1][myp1][mz].density   += (particles[i].mass*Wx2*Wy2*Wz1)/pow3(H);
      grid[mxp1][myp1][mzp1].density += (particles[i].mass*Wx2*Wy2*Wz2)/pow3(H);
    }
  }
  
}


//----------------------------------------
//potential_gridpoint------> phi(g_klm)
//----------------------------------------
void Potential_GridPoint(Grid_Point ***grid, int N_gridpoints, double H)
{
  int i, j, k, iterations;
  int mxm1, mym1, mzm1, mxp1, myp1, mzp1;

  //make grid potential equal to ZERO AGAIN
    for(i=0; i<=(N_gridpoints-1); i++) //for the "red" grid points
    {
      for(j=0; j<=(N_gridpoints-1); j++)
      {
        for(k=0; k<=(N_gridpoints-1); k++)
        {
          grid[i][j][k].potential = 0.0;
        }
      }
    }

  //start iterations
  for(iterations=0; iterations<60; iterations++)
  {
    for(i=0; i<=(N_gridpoints-1); i++) //for the "red" grid points
    {
      for(j=0; j<=(N_gridpoints-1); j++)
      {
        for(k=0; k<=(N_gridpoints-1); k++)
        {
          if((i+j+k)%2 == 0)
          {
            // get indices of neighbouring nodes taking into account periodic boundaries
            mxm1 = (int) floor(fmod((double)(i-1+N_gridpoints),(double)(N_gridpoints)));
            mym1 = (int) floor(fmod((double)(j-1+N_gridpoints),(double)(N_gridpoints)));
            mzm1 = (int) floor(fmod((double)(k-1+N_gridpoints),(double)(N_gridpoints)));
            mxp1 = (int) floor(fmod((double)(i+1),(double)(N_gridpoints)));
            myp1 = (int) floor(fmod((double)(j+1),(double)(N_gridpoints)));
            mzp1 = (int) floor(fmod((double)(k+1),(double)(N_gridpoints)));
            
            //printf("PAR %d  %d  %d \n", i ,j ,k);
            grid[i][j][k].potential = grid[mxp1][j][k].potential + grid[mxm1][j][k].potential + grid[i][myp1][k].potential + grid[i][mym1][k].potential;
            grid[i][j][k].potential = grid[i][j][k].potential + grid[i][j][mzp1].potential + grid[i][j][mzm1].potential - piG4*H*H*grid[i][j][k].density;
            grid[i][j][k].potential = grid[i][j][k].potential/6.0;
          }
        }
      }
    }
    for(i=0; i<=(N_gridpoints-1); i++) //for the "black" grid points
    {
      for(j=0; j<=(N_gridpoints-1); j++)
      {
        for(k=0; k<=(N_gridpoints-1); k++)
        {
          if((i+j+k)%2 == 1)
          {
            // get indices of neighbouring nodes taking into account periodic boundaries
            mxm1 = (int) floor(fmod((double)(i-1+N_gridpoints),(double)(N_gridpoints)));
            mym1 = (int) floor(fmod((double)(j-1+N_gridpoints),(double)(N_gridpoints)));
            mzm1 = (int) floor(fmod((double)(k-1+N_gridpoints),(double)(N_gridpoints)));
            mxp1 = (int) floor(fmod((double)(i+1),(double)(N_gridpoints)));
            myp1 = (int) floor(fmod((double)(j+1),(double)(N_gridpoints)));
            mzp1 = (int) floor(fmod((double)(k+1),(double)(N_gridpoints)));
            
            //printf("IMPAR %d  %d  %d \n", i ,j ,k);
            grid[i][j][k].potential = grid[mxp1][j][k].potential + grid[mxm1][j][k].potential + grid[i][myp1][k].potential + grid[i][mym1][k].potential;
            grid[i][j][k].potential = grid[i][j][k].potential + grid[i][j][mzp1].potential + grid[i][j][mzm1].potential - piG4*H*H*grid[i][j][k].density;
            grid[i][j][k].potential = grid[i][j][k].potential/6.0;
          }
        }
      }
    }
  }

}


//----------------------------------------
//Force_GridPoint calculator ---> F(g_klm)
//----------------------------------------
void Acceleration_GridPoint(Grid_Point ***grid, int N_gridpoints, double H)
{
  int i, j, k;
  int mxm1, mym1, mzm1, mxp1, myp1, mzp1;
  
  
  for(i=0; i<(N_gridpoints); i++)
  {
    for(j=0; j<(N_gridpoints); j++)
    {
      for(k=0; k<(N_gridpoints); k++)
      {
        // get indices of neighbouring nodes taking into account periodic boundaries
        mxm1 = (int) floor(fmod((double)(i-1+N_gridpoints),(double)(N_gridpoints)));
        mym1 = (int) floor(fmod((double)(j-1+N_gridpoints),(double)(N_gridpoints)));
        mzm1 = (int) floor(fmod((double)(k-1+N_gridpoints),(double)(N_gridpoints)));
        mxp1 = (int) floor(fmod((double)(i+1),(double)(N_gridpoints)));
        myp1 = (int) floor(fmod((double)(j+1),(double)(N_gridpoints)));
        mzp1 = (int) floor(fmod((double)(k+1),(double)(N_gridpoints)));
 
        grid[i][j][k].accelerationx = -(grid[mxp1][j][k].potential - grid[mxm1][j][k].potential)/2.0/H;
        grid[i][j][k].accelerationy = -(grid[i][myp1][k].potential - grid[i][mym1][k].potential)/2.0/H;
        grid[i][j][k].accelerationz = -(grid[i][j][mzp1].potential - grid[i][j][mzm1].potential)/2.0/H;
      }
    }
  }

}

//----------------------------------------
//Interpolate_Force_particles calculator ---> F(r)
//----------------------------------------
void Force_particles(Particle *particles, Grid_Point ***grid, int N_part, int N_gridpoints, double H)
{
  
  
  int    i,j,k;
  int    count_x, count_y, count_z;//it says if the particle is in the first or in the second half of the grid poitn interval
  double Wx1, Wy1, Wz1, Wx2, Wy2, Wz2;
  double nx, ny, nz; // (double) number of grid points
  int    mx, my, mz; // (int) number of grid points
  int    mxm1, mxp1, mym1, myp1, mzm1, mzp1;
  double dx, dy, dz; // distance between center of the grid point and the particle
//#pragma omp parallel for default(none) private(i) shared(N_part,N_gridpoints, particles, grid,H,nx,ny,nz,mx,my,mz,mxm1,mym1,mzm1,mxp1,myp1,mzp1,Wx1,Wx2,Wy1,Wy2,Wz1,Wz2,dx,dy,dz,count_x,count_y,count_z) 
  for(i=0; i<N_part; i++)
  {
    nx = particles[i].x/H;
    ny = particles[i].y/H;
    nz = particles[i].z/H;
    
    mx = (int)floor(nx);
    my = (int)floor(ny);
    mz = (int)floor(nz);
    
    // get indices of neighbouring nodes taking into account periodic boundaries
    mxm1 = (int) floor(fmod((double)(mx-1+N_gridpoints),(double)(N_gridpoints)));
    mym1 = (int) floor(fmod((double)(my-1+N_gridpoints),(double)(N_gridpoints)));
    mzm1 = (int) floor(fmod((double)(mz-1+N_gridpoints),(double)(N_gridpoints)));
    mxp1 = (int) floor(fmod((double)(mx+1),(double)(N_gridpoints)));
    myp1 = (int) floor(fmod((double)(my+1),(double)(N_gridpoints)));
    mzp1 = (int) floor(fmod((double)(mz+1),(double)(N_gridpoints)));

    //look for the half part of the grid point in which the particle is:
    dx = (nx - mx)*H;
    if( (nx - mx)<0.5 )
    {
      count_x = 0;
      Wx1 = 1.0 - ((dx + H/2.0)/H);
      Wx2 = 1.0 - ((H/2.0 - dx)/H);
    }
    else
    {
      count_x = 1;
      Wx1 = 1.0 - ((dx - H/2.0)/H);
      Wx2 = 1.0 - (((H-dx) + H/2.0)/H);
    }
    
    
    dy = (ny - my)*H;
    if( (ny - my)<0.5 )
    {
      count_y = 0;
      Wy1 = 1.0 - ((dy + H/2.0)/H);
      Wy2 = 1.0 - ((H/2.0 - dy)/H);
    }
    else
    {
      count_y = 1;
      Wy1 = 1.0 - ((dy - H/2.0)/H);
      Wy2 = 1.0 - (((H-dy) + H/2.0)/H);
    }
    
    
    dz = (nz - mz)*H;
    if( (nz - mz)<0.5 )
    {
      count_z = 0;
      Wz1 = 1.0 - ((dz + H/2.0)/H);
      Wz2 = 1.0 - ((H/2.0 - dz)/H);
    }
    else
    {
      count_z = 1;
      Wz1 = 1.0 - ((dz - H/2.0)/H);
      Wz2 = 1.0 - (((H-dz) + H/2.0)/H);
    }

    // assing the x-force to the corresponding particle 
    if (count_x==0 && count_y==0 && count_z==0)
    {
      particles[i].forcex = grid[mx][my][mz].accelerationx*Wx2*Wy2*Wz2     + grid[mx][my][mzm1].accelerationx*Wx2*Wy2*Wz1 
                          + grid[mx][mym1][mz].accelerationx*Wx2*Wy1*Wz2   + grid[mx][mym1][mzm1].accelerationx*Wx2*Wy1*Wz1 
                          + grid[mxm1][my][mz].accelerationx*Wx1*Wy2*Wz2   + grid[mxm1][my][mzm1].accelerationx*Wx1*Wy2*Wz1 
                          + grid[mxm1][mym1][mz].accelerationx*Wx1*Wy1*Wz2 + grid[mxm1][mym1][mzm1].accelerationx*Wx1*Wy1*Wz1;
      particles[i].forcex = particles[i].forcex*particles[i].mass;
    }
    else if (count_x==0 && count_y==0 && count_z==1)
    {
      particles[i].forcex = grid[mx][my][mz].accelerationx*Wx2*Wy2*Wz1     + grid[mx][my][mzp1].accelerationx*Wx2*Wy2*Wz2 
                          + grid[mx][mym1][mz].accelerationx*Wx2*Wy1*Wz1   + grid[mx][mym1][mzp1].accelerationx*Wx2*Wy1*Wz2 
                          + grid[mxm1][my][mz].accelerationx*Wx1*Wy2*Wz1   + grid[mxm1][my][mzp1].accelerationx*Wx1*Wy2*Wz2 
                          + grid[mxm1][mym1][mz].accelerationx*Wx1*Wy1*Wz1 + grid[mxm1][mym1][mzp1].accelerationx*Wx1*Wy1*Wz2;
      particles[i].forcex = particles[i].forcex*particles[i].mass;
    }
    else if (count_x==0 && count_y==1 && count_z==0)
    {
      particles[i].forcex = grid[mx][my][mz].accelerationx*Wx2*Wy1*Wz2     + grid[mx][my][mzm1].accelerationx*Wx2*Wy1*Wz1
                          + grid[mx][myp1][mz].accelerationx*Wx2*Wy2*Wz2   + grid[mx][myp1][mzm1].accelerationx*Wx2*Wy2*Wz1
                          + grid[mxm1][my][mz].accelerationx*Wx1*Wy1*Wz2   + grid[mxm1][my][mzm1].accelerationx*Wx1*Wy1*Wz1 
                          + grid[mxm1][myp1][mz].accelerationx*Wx1*Wy2*Wz2 + grid[mxm1][myp1][mzm1].accelerationx*Wx1*Wy2*Wz1;
      particles[i].forcex = particles[i].forcex*particles[i].mass;
    }
    else if (count_x==0 && count_y==1 && count_z==1)
    {
      particles[i].forcex = grid[mx][my][mz].accelerationx*Wx2*Wy1*Wz1     + grid[mx][my][mzp1].accelerationx*Wx2*Wy1*Wz2 
                          + grid[mx][myp1][mz].accelerationx*Wx2*Wy2*Wz1   + grid[mx][myp1][mzp1].accelerationx*Wx2*Wy2*Wz2 
                          + grid[mxm1][my][mz].accelerationx*Wx1*Wy1*Wz1   + grid[mxm1][my][mzp1].accelerationx*Wx1*Wy1*Wz2 
                          + grid[mxm1][myp1][mz].accelerationx*Wx1*Wy2*Wz1 + grid[mxm1][myp1][mzp1].accelerationx*Wx1*Wy2*Wz2;
      particles[i].forcex = particles[i].forcex*particles[i].mass;
    }
    
    else if (count_x==1 && count_y==0 && count_z==0)
    {
      particles[i].forcex = grid[mx][my][mz].accelerationx*Wx1*Wy2*Wz2     + grid[mx][my][mzm1].accelerationx*Wx1*Wy2*Wz1 
                          + grid[mx][mym1][mz].accelerationx*Wx1*Wy1*Wz2   + grid[mx][mym1][mzm1].accelerationx*Wx1*Wy1*Wz1 
                          + grid[mxp1][my][mz].accelerationx*Wx2*Wy2*Wz2   + grid[mxp1][my][mzm1].accelerationx*Wx2*Wy2*Wz1
                          + grid[mxp1][mym1][mz].accelerationx*Wx2*Wy1*Wz2 + grid[mxp1][mym1][mzm1].accelerationx*Wx2*Wy1*Wz1;
      particles[i].forcex = particles[i].forcex*particles[i].mass;
    }
    else if (count_x==1 && count_y==0 && count_z==1)
    {
      particles[i].forcex = grid[mx][my][mz].accelerationx*Wx1*Wy2*Wz1     + grid[mx][my][mzp1].accelerationx*Wx1*Wy2*Wz2 
                          + grid[mx][mym1][mz].accelerationx*Wx1*Wy1*Wz1   + grid[mx][mym1][mzp1].accelerationx*Wx1*Wy1*Wz2 
                          + grid[mxp1][my][mz].accelerationx*Wx2*Wy2*Wz1   + grid[mxp1][my][mzp1].accelerationx*Wx2*Wy2*Wz2 
                          + grid[mxp1][mym1][mz].accelerationx*Wx2*Wy1*Wz1 + grid[mxp1][mym1][mzp1].accelerationx*Wx2*Wy1*Wz2;
      particles[i].forcex = particles[i].forcex*particles[i].mass;
    }
    else if (count_x==1 && count_y==1 && count_z==0)
    {
      particles[i].forcex = grid[mx][my][mz].accelerationx*Wx1*Wy1*Wz2     + grid[mx][my][mzm1].accelerationx*Wx1*Wy1*Wz1 
                          + grid[mx][myp1][mz].accelerationx*Wx1*Wy2*Wz2   + grid[mx][myp1][mzm1].accelerationx*Wx1*Wy2*Wz1 
                          + grid[mxp1][my][mz].accelerationx*Wx2*Wy1*Wz2   + grid[mxp1][my][mzm1].accelerationx*Wx2*Wy1*Wz1 
                          + grid[mxp1][myp1][mz].accelerationx*Wx2*Wy2*Wz2 + grid[mxp1][myp1][mzm1].accelerationx*Wx2*Wy2*Wz1;
      particles[i].forcex = particles[i].forcex*particles[i].mass;
    }
    else if (count_x==1 && count_y==1 && count_z==1)
    {
      particles[i].forcex = grid[mx][my][mz].accelerationx*Wx1*Wy1*Wz1     + grid[mx][my][mzp1].accelerationx*Wx1*Wy1*Wz2 
                          + grid[mx][myp1][mz].accelerationx*Wx1*Wy2*Wz1   + grid[mx][myp1][mzp1].accelerationx*Wx1*Wy2*Wz2 
                          + grid[mxp1][my][mz].accelerationx*Wx2*Wy1*Wz1   + grid[mxp1][my][mzp1].accelerationx*Wx2*Wy1*Wz2
                          + grid[mxp1][myp1][mz].accelerationx*Wx2*Wy2*Wz1 + grid[mxp1][myp1][mzp1].accelerationx*Wx2*Wy2*Wz2;
      particles[i].forcex = particles[i].forcex*particles[i].mass;
    }
    
    // assing the y-force to the corresponding particle 
    if (count_x==0 && count_y==0 && count_z==0)
    {
      particles[i].forcey = grid[mx][my][mz].accelerationy*Wx2*Wy2*Wz2     + grid[mx][my][mzm1].accelerationy*Wx2*Wy2*Wz1 
                          + grid[mx][mym1][mz].accelerationy*Wx2*Wy1*Wz2   + grid[mx][mym1][mzm1].accelerationy*Wx2*Wy1*Wz1 
                          + grid[mxm1][my][mz].accelerationy*Wx1*Wy2*Wz2   + grid[mxm1][my][mzm1].accelerationy*Wx1*Wy2*Wz1 
                          + grid[mxm1][mym1][mz].accelerationy*Wx1*Wy1*Wz2 + grid[mxm1][mym1][mzm1].accelerationy*Wx1*Wy1*Wz1;
      particles[i].forcey = particles[i].forcey*particles[i].mass;
    }
    else if (count_x==0 && count_y==0 && count_z==1)
    {
      particles[i].forcey = grid[mx][my][mz].accelerationy*Wx2*Wy2*Wz1     + grid[mx][my][mzp1].accelerationy*Wx2*Wy2*Wz2 
                          + grid[mx][mym1][mz].accelerationy*Wx2*Wy1*Wz1   + grid[mx][mym1][mzp1].accelerationy*Wx2*Wy1*Wz2 
                          + grid[mxm1][my][mz].accelerationy*Wx1*Wy2*Wz1   + grid[mxm1][my][mzp1].accelerationy*Wx1*Wy2*Wz2 
                          + grid[mxm1][mym1][mz].accelerationy*Wx1*Wy1*Wz1 + grid[mxm1][mym1][mzp1].accelerationy*Wx1*Wy1*Wz2;
      particles[i].forcey = particles[i].forcey*particles[i].mass;
    }
    else if (count_x==0 && count_y==1 && count_z==0)
    {
      particles[i].forcey = grid[mx][my][mz].accelerationy*Wx2*Wy1*Wz2     + grid[mx][my][mzm1].accelerationy*Wx2*Wy1*Wz1
                          + grid[mx][myp1][mz].accelerationy*Wx2*Wy2*Wz2   + grid[mx][myp1][mzm1].accelerationy*Wx2*Wy2*Wz1
                          + grid[mxm1][my][mz].accelerationy*Wx1*Wy1*Wz2   + grid[mxm1][my][mzm1].accelerationy*Wx1*Wy1*Wz1 
                          + grid[mxm1][myp1][mz].accelerationy*Wx1*Wy2*Wz2 + grid[mxm1][myp1][mzm1].accelerationy*Wx1*Wy2*Wz1;
      particles[i].forcey = particles[i].forcey*particles[i].mass;
    }
    else if (count_x==0 && count_y==1 && count_z==1)
    {
      particles[i].forcey = grid[mx][my][mz].accelerationy*Wx2*Wy1*Wz1     + grid[mx][my][mzp1].accelerationy*Wx2*Wy1*Wz2 
                          + grid[mx][myp1][mz].accelerationy*Wx2*Wy2*Wz1   + grid[mx][myp1][mzp1].accelerationy*Wx2*Wy2*Wz2 
                          + grid[mxm1][my][mz].accelerationy*Wx1*Wy1*Wz1   + grid[mxm1][my][mzp1].accelerationy*Wx1*Wy1*Wz2 
                          + grid[mxm1][myp1][mz].accelerationy*Wx1*Wy2*Wz1 + grid[mxm1][myp1][mzp1].accelerationy*Wx1*Wy2*Wz2;
      particles[i].forcey = particles[i].forcey*particles[i].mass;
    }
    
    else if (count_x==1 && count_y==0 && count_z==0)
    {
      particles[i].forcey = grid[mx][my][mz].accelerationy*Wx1*Wy2*Wz2     + grid[mx][my][mzm1].accelerationy*Wx1*Wy2*Wz1 
                          + grid[mx][mym1][mz].accelerationy*Wx1*Wy1*Wz2   + grid[mx][mym1][mzm1].accelerationy*Wx1*Wy1*Wz1 
                          + grid[mxp1][my][mz].accelerationy*Wx2*Wy2*Wz2   + grid[mxp1][my][mzm1].accelerationy*Wx2*Wy2*Wz1
                          + grid[mxp1][mym1][mz].accelerationy*Wx2*Wy1*Wz2 + grid[mxp1][mym1][mzm1].accelerationy*Wx2*Wy1*Wz1;
      particles[i].forcey = particles[i].forcey*particles[i].mass;
    }
    else if (count_x==1 && count_y==0 && count_z==1)
    {
      particles[i].forcey = grid[mx][my][mz].accelerationy*Wx1*Wy2*Wz1     + grid[mx][my][mzp1].accelerationy*Wx1*Wy2*Wz2 
                          + grid[mx][mym1][mz].accelerationy*Wx1*Wy1*Wz1   + grid[mx][mym1][mzp1].accelerationy*Wx1*Wy1*Wz2 
                          + grid[mxp1][my][mz].accelerationy*Wx2*Wy2*Wz1   + grid[mxp1][my][mzp1].accelerationy*Wx2*Wy2*Wz2 
                          + grid[mxp1][mym1][mz].accelerationy*Wx2*Wy1*Wz1 + grid[mxp1][mym1][mzp1].accelerationy*Wx2*Wy1*Wz2;
      particles[i].forcey = particles[i].forcey*particles[i].mass;
    }
    else if (count_x==1 && count_y==1 && count_z==0)
    {
      particles[i].forcey = grid[mx][my][mz].accelerationy*Wx1*Wy1*Wz2     + grid[mx][my][mzm1].accelerationy*Wx1*Wy1*Wz1 
                          + grid[mx][myp1][mz].accelerationy*Wx1*Wy2*Wz2   + grid[mx][myp1][mzm1].accelerationy*Wx1*Wy2*Wz1 
                          + grid[mxp1][my][mz].accelerationy*Wx2*Wy1*Wz2   + grid[mxp1][my][mzm1].accelerationy*Wx2*Wy1*Wz1 
                          + grid[mxp1][myp1][mz].accelerationy*Wx2*Wy2*Wz2 + grid[mxp1][myp1][mzm1].accelerationy*Wx2*Wy2*Wz1;
      particles[i].forcey = particles[i].forcey*particles[i].mass;
    }
    else if (count_x==1 && count_y==1 && count_z==1)
    {
      particles[i].forcey = grid[mx][my][mz].accelerationy*Wx1*Wy1*Wz1     + grid[mx][my][mzp1].accelerationy*Wx1*Wy1*Wz2 
                          + grid[mx][myp1][mz].accelerationy*Wx1*Wy2*Wz1   + grid[mx][myp1][mzp1].accelerationy*Wx1*Wy2*Wz2 
                          + grid[mxp1][my][mz].accelerationy*Wx2*Wy1*Wz1   + grid[mxp1][my][mzp1].accelerationy*Wx2*Wy1*Wz2
                          + grid[mxp1][myp1][mz].accelerationy*Wx2*Wy2*Wz1 + grid[mxp1][myp1][mzp1].accelerationy*Wx2*Wy2*Wz2;
      particles[i].forcey = particles[i].forcey*particles[i].mass;
    }
    

    
    // assing the z-force to the corresponding particle 
    if (count_x==0 && count_y==0 && count_z==0)
    {
      particles[i].forcez = grid[mx][my][mz].accelerationz*Wx2*Wy2*Wz2     + grid[mx][my][mzm1].accelerationz*Wx2*Wy2*Wz1 
                          + grid[mx][mym1][mz].accelerationz*Wx2*Wy1*Wz2   + grid[mx][mym1][mzm1].accelerationz*Wx2*Wy1*Wz1 
                          + grid[mxm1][my][mz].accelerationz*Wx1*Wy2*Wz2   + grid[mxm1][my][mzm1].accelerationz*Wx1*Wy2*Wz1 
                          + grid[mxm1][mym1][mz].accelerationz*Wx1*Wy1*Wz2 + grid[mxm1][mym1][mzm1].accelerationz*Wx1*Wy1*Wz1;
      particles[i].forcez = particles[i].forcez*particles[i].mass;
    }
    else if (count_x==0 && count_y==0 && count_z==1)
    {
      particles[i].forcez = grid[mx][my][mz].accelerationz*Wx2*Wy2*Wz1     + grid[mx][my][mzp1].accelerationz*Wx2*Wy2*Wz2 
                          + grid[mx][mym1][mz].accelerationz*Wx2*Wy1*Wz1   + grid[mx][mym1][mzp1].accelerationz*Wx2*Wy1*Wz2 
                          + grid[mxm1][my][mz].accelerationz*Wx1*Wy2*Wz1   + grid[mxm1][my][mzp1].accelerationz*Wx1*Wy2*Wz2 
                          + grid[mxm1][mym1][mz].accelerationz*Wx1*Wy1*Wz1 + grid[mxm1][mym1][mzp1].accelerationz*Wx1*Wy1*Wz2;
      particles[i].forcez = particles[i].forcez*particles[i].mass;
    }
    else if (count_x==0 && count_y==1 && count_z==0)
    {
      particles[i].forcez = grid[mx][my][mz].accelerationz*Wx2*Wy1*Wz2     + grid[mx][my][mzm1].accelerationz*Wx2*Wy1*Wz1
                          + grid[mx][myp1][mz].accelerationz*Wx2*Wy2*Wz2   + grid[mx][myp1][mzm1].accelerationz*Wx2*Wy2*Wz1
                          + grid[mxm1][my][mz].accelerationz*Wx1*Wy1*Wz2   + grid[mxm1][my][mzm1].accelerationz*Wx1*Wy1*Wz1 
                          + grid[mxm1][myp1][mz].accelerationz*Wx1*Wy2*Wz2 + grid[mxm1][myp1][mzm1].accelerationz*Wx1*Wy2*Wz1;
      particles[i].forcez = particles[i].forcez*particles[i].mass;
    }
    else if (count_x==0 && count_y==1 && count_z==1)
    {
      particles[i].forcez = grid[mx][my][mz].accelerationz*Wx2*Wy1*Wz1     + grid[mx][my][mzp1].accelerationz*Wx2*Wy1*Wz2 
                          + grid[mx][myp1][mz].accelerationz*Wx2*Wy2*Wz1   + grid[mx][myp1][mzp1].accelerationz*Wx2*Wy2*Wz2 
                          + grid[mxm1][my][mz].accelerationz*Wx1*Wy1*Wz1   + grid[mxm1][my][mzp1].accelerationz*Wx1*Wy1*Wz2 
                          + grid[mxm1][myp1][mz].accelerationz*Wx1*Wy2*Wz1 + grid[mxm1][myp1][mzp1].accelerationz*Wx1*Wy2*Wz2;
      particles[i].forcez = particles[i].forcez*particles[i].mass;
    }
    
    else if (count_x==1 && count_y==0 && count_z==0)
    {
      particles[i].forcez = grid[mx][my][mz].accelerationz*Wx1*Wy2*Wz2     + grid[mx][my][mzm1].accelerationz*Wx1*Wy2*Wz1 
                          + grid[mx][mym1][mz].accelerationz*Wx1*Wy1*Wz2   + grid[mx][mym1][mzm1].accelerationz*Wx1*Wy1*Wz1 
                          + grid[mxp1][my][mz].accelerationz*Wx2*Wy2*Wz2   + grid[mxp1][my][mzm1].accelerationz*Wx2*Wy2*Wz1
                          + grid[mxp1][mym1][mz].accelerationz*Wx2*Wy1*Wz2 + grid[mxp1][mym1][mzm1].accelerationz*Wx2*Wy1*Wz1;
      particles[i].forcez = particles[i].forcez*particles[i].mass;
    }
    else if (count_x==1 && count_y==0 && count_z==1)
    {
      particles[i].forcez = grid[mx][my][mz].accelerationz*Wx1*Wy2*Wz1     + grid[mx][my][mzp1].accelerationz*Wx1*Wy2*Wz2 
                          + grid[mx][mym1][mz].accelerationz*Wx1*Wy1*Wz1   + grid[mx][mym1][mzp1].accelerationz*Wx1*Wy1*Wz2 
                          + grid[mxp1][my][mz].accelerationz*Wx2*Wy2*Wz1   + grid[mxp1][my][mzp1].accelerationz*Wx2*Wy2*Wz2 
                          + grid[mxp1][mym1][mz].accelerationz*Wx2*Wy1*Wz1 + grid[mxp1][mym1][mzp1].accelerationz*Wx2*Wy1*Wz2;
      particles[i].forcez = particles[i].forcez*particles[i].mass;
    }
    else if (count_x==1 && count_y==1 && count_z==0)
    {
      particles[i].forcez = grid[mx][my][mz].accelerationz*Wx1*Wy1*Wz2     + grid[mx][my][mzm1].accelerationz*Wx1*Wy1*Wz1 
                          + grid[mx][myp1][mz].accelerationz*Wx1*Wy2*Wz2   + grid[mx][myp1][mzm1].accelerationz*Wx1*Wy2*Wz1 
                          + grid[mxp1][my][mz].accelerationz*Wx2*Wy1*Wz2   + grid[mxp1][my][mzm1].accelerationz*Wx2*Wy1*Wz1 
                          + grid[mxp1][myp1][mz].accelerationz*Wx2*Wy2*Wz2 + grid[mxp1][myp1][mzm1].accelerationz*Wx2*Wy2*Wz1;
      particles[i].forcez = particles[i].forcez*particles[i].mass;
    }
    else if (count_x==1 && count_y==1 && count_z==1)
    {
      particles[i].forcez = grid[mx][my][mz].accelerationz*Wx1*Wy1*Wz1     + grid[mx][my][mzp1].accelerationz*Wx1*Wy1*Wz2 
                          + grid[mx][myp1][mz].accelerationz*Wx1*Wy2*Wz1   + grid[mx][myp1][mzp1].accelerationz*Wx1*Wy2*Wz2 
                          + grid[mxp1][my][mz].accelerationz*Wx2*Wy1*Wz1   + grid[mxp1][my][mzp1].accelerationz*Wx2*Wy1*Wz2
                          + grid[mxp1][myp1][mz].accelerationz*Wx2*Wy2*Wz1 + grid[mxp1][myp1][mzp1].accelerationz*Wx2*Wy2*Wz2;
      particles[i].forcez = particles[i].forcez*particles[i].mass;
    }

  }

}


//----------------------------------------
// LeapFrog_particles()
//----------------------------------------
void LeapFrog_particles(Particle *particles, int N_part, double delta_t, double H, int N_gridpoints)
{
  int    i;
  double nx, ny, nz;
  double max_dist;
  int    mx, my, mz;

  max_dist = N_gridpoints*H;
  
  for(i=0; i<N_part; i++)
  {
    particles[i].x  = particles[i].x + delta_t*particles[i].vx;
    particles[i].y  = particles[i].y + delta_t*particles[i].vy;
    particles[i].z  = particles[i].z + delta_t*particles[i].vz;
    particles[i].vx = particles[i].vx + delta_t*(particles[i].forcex/particles[i].mass);
    particles[i].vy = particles[i].vy + delta_t*(particles[i].forcey/particles[i].mass);
    particles[i].vz = particles[i].vz + delta_t*(particles[i].forcez/particles[i].mass);
    

   //Periodic boundary conditions
    particles[i].x =  (fmod((double)(particles[i].x + max_dist),(double)(max_dist)));
    particles[i].y =  (fmod((double)(particles[i].y + max_dist),(double)(max_dist)));
    particles[i].z =  (fmod((double)(particles[i].z + max_dist),(double)(max_dist)));
  }

}







