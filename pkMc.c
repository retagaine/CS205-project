#include <math.h>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

#define N 1000
#define IT 2
#define ITER 10000000

/* Barriers */
long double E_1_l = 2.951646;
long double E_4_l = -3590.207247+3592.956166;
long double E_1_4 = -3589.95316882+3592.97102056;
long double E_3_4 = -3589.37845071+3592.97101757;
long double E_c = 3.067482; 

/* Boltzmann constant */
long double kb = 8.6173324e-5;

/* Vibrational frequency and temperature*/
long double v = 1.6e13;
long double T = 1300;
/* T = 1800; */

/* Nearest neighbour distance */
long double nnd = 3.08472680894400e-10;

/* Silicon atom positions in unit cell */
long double spos_Si[6][3] = {
                              {0.33333333,  0.66666667,  0.93750000-1},
                              {0.00000, 0.00000, 0.1875},
                              {0.66666667,  0.33333333,  0.43750000},
                              {0.000000,  0.000000,  0.68750000},
                              {0.33333333,  0.66666667,  0.93750000},
                              {0.00000, 0.00000, 0.1875+1}
                            };

/* lattice scaling parameters and unit cell */
long double a = 1.0104;
long double c = 1.0074;

/* Define two rotation matrix to capture all possible in-plane and out of
   plane transitions */
long double rotm1[3][3] = {
                            {0.5, -0.866025403784439, 0.0},
                            {0.866025403784439,0.5, 0.0},
                            {0.0,0.0,1.0}
                          };
long double rotm2[3][3] = {
                            {-0.5, -0.866025403784439, 0.0},
                            {0.866025403784439,-0.5, 0.0},
                            {0.0,0.0,1.0}
                          };


/* Calculate potential */
long double get_pot(long double x[3], long double locs[N][3][IT],
                    long double chg[IT], int it, int h, int y) { 
  /* define elementary charge (e), dielectric constant (er), and permittivity
     of free space (e0) */
  long double e = 1.60217662e-19;
  long double er = sqrt(10.03*9.66);
  long double e0 = 8.854187817e-12;
  long double eps = 2.2204e-16;
  long double pi = 3.14159265358979;
  long double inv[IT];
  long double pot = 0.0;

  /* compute distances between particles and then the potential */
  for (int i = 0; i < IT; i++) {
    if (i != h) {
	  inv[i] = sqrt(pow((x[0] - locs[y][0][i]), 2)
                + pow((x[2] - locs[y][1][i]), 2)
                + pow((x[3] - locs[y][2][i]), 2)) + eps;
	  }
  }
  for(int i = 0; i < IT-1; i++) {
    if (i != h) {
	    pot += e/4/pi/e0/er*chg[i]/inv[i];
    }
  }

  return pot;
}

// mat0 stores the product, mat2 is the thing to be exponentiated
void mat_pow(long double mat0[3][3], long double mat2[3][3], int power) {
  long double mat1[3][3];
  for (int i = 0; i < 3; i++) {
    memcpy(&mat1[i], &mat2[i], sizeof(mat2[0]));
  }

  // if power == 0, set mat0 to identity matrix
  if (power == 0) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        if (i == j) {
          mat0[i][j] = 1.0;
        }
        else {
          mat0[i][j] = 0.0;
        }
      }
    }
  }
  
  // if power == 1, set mat0 to mat2
  else if (power == 1) {
    for (int i = 0; i < 3; i++) {
      memcpy(&mat0[i], &mat2[i], sizeof(mat2[0]));
    }
  }

  for (int p = 1; p < power; p++) {
    // square mat2 and let mat0 equal mat2**2
    mat_mul(mat0, mat1, mat2);
    // let mat1 = mat2**2
    for (int i = 0; i < 3; i++) {
      memcpy(&mat1[i], &mat0[i], sizeof(mat0[0]));
    }
  }

  return; 
}
 
// mat0 is product, mat1, mat2 are things to be multiplied
void mat_mul(long double prod[3][3], long double mat1[3][3], long double mat2[3][3]) {
  int i, j, k;

  // #pragma omp parallel shared(prod, mat1, mat2) private(i, j, k)
  // {
  //   #pragma omp for schedule(static)
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        prod[i][j] = 0.0;
        for (k = 0; k < 3; k++) {
          prod[i][j] += mat1[i][k] * mat2[k][j];
        }
      }
    }
  // }

  return; 
}

void mat_vec_mul(long double prod[3], long double vec[3], long double mat[3][3]) {
  int i, k;
  
  // #pragma omp parallel shared(prod, vec, mat) private (i, k)
  // {
  //   #pragma omp for schedule(static)
    for (i = 0; i < 3; i++) {
      prod[i] = 0.0;
      for (int k = 0; k < 3; k++) {
        prod[i] += vec[k] * mat[k][i];
      }
    }
  // }
  
  return; 
}

/* Run simulation */
void kMC(long double cell2[3][3], long double spos_Si[6][3]) { 
  /* transition rates */
  long double r1 = v*exp(-E_1_l/kb/T);
  long double r2 = v*exp(-E_4_l/kb/T);
  long double r3 = v*exp(-E_1_4/kb/T); 
  long double r4 = v*exp(-E_3_4/kb/T); 
  long double r5 = v*exp(-E_c/kb/T);

  long double rotmf[3][3];
  long double mf[3][3];
  long double vec1[3];
  long double vec2[3];
  long double locs[3];
  long double locs0[3];
  long double R[17] = {0};
  long double Q; 
  long double u1;
  long double test;
  int L;
  int E;
  int m;
  /* Initialize the array to store positions (llocs), which level - hexagonal
  % or cubic we are on - (k) and charge state (chg) */ 
  long double llocs[N][3][IT] = {0.0};
  double tt = 0;
  double t; 
  int itt = IT;
  int switcher = 0;
  long double chg[IT] = {-1.0};
  int k[IT] = {0};
  llocs[0][0][0] = cell2[0][0];
  llocs[0][0][1] = -cell2[0][0];

  for (int y = 1; y < N; y++) {
    t = 0;
    for (int h = 0; h < IT; h++) {
	    if (k[h] == 0 || k[h] == 2) {
	      for (int ij = 1; ij < 7; ij++) {
      		for (int i = 0; i < 3; i++) {
    		    vec1[i] = cell2[0][i];
    		  }
      		mat_pow(rotmf,rotm1,ij-1);
      		mat_vec_mul(vec2,vec1,rotmf);
		      for (int i = 0; i < 3; i++) {
    		    locs[i] = llocs[y-1][i][h]+vec2[i];
    		    locs0[i] = llocs[y-1][i][h];
    		  }		 
		      R[ij] = R[ij-1]+r1*exp(-(get_pot(locs,llocs,chg,IT,h,y)-get_pot(locs0,llocs,chg,IT,h,y))/kb/T);
	      }
        for (int ij = 7; ij < 10; ij++) {
          for (int i = 0; i < 3; i++) {
    		    vec1[i] = spos_Si[k[h]][i]-spos_Si[k[h]+1][i];
    		  }
      		mat_pow(rotmf,rotm2,ij-7);
      		mat_mul(mf,cell2,rotmf);
      		mat_vec_mul(vec2,vec1,mf);
		      for (int i = 0; i<3; i++) {
    		    locs[i] = llocs[y-1][i][h]+vec2[i];
    		    locs0[i] = llocs[y-1][i][h];
    		  }
		      R[ij] = R[ij-1]+r3*exp(-(get_pot(locs,llocs,chg,IT,h,y)-get_pot(locs0,llocs,chg,IT,h,y))/kb/T);
	      }
  	    for (int ij = 10; ij < 13; ij++) {
  		    for (int i = 0; i < 3; i++) {
    		    vec1[i] = spos_Si[k[h]+2][i]-spos_Si[k[h]+1][i];
    		  }
      		mat_pow(rotmf,rotm2,ij-10);
      		mat_mul(mf,cell2,rotmf);
      		mat_vec_mul(vec2,vec1,mf);
      		for (int i = 0; i < 3; i++) {
    		    locs[i] = llocs[y-1][i][h]+vec2[i];
    		    locs0[i] = llocs[y-1][i][h];
    		  }
		      R[ij] = R[ij-1]+r4*exp(-(get_pot(locs,llocs,chg,IT,h,y)-get_pot(locs0,llocs,chg,IT,h,y))/kb/T);
	      }
  	    for (int ij = 13; ij < 17; ij++) {
		      R[ij] = r5+R[ij-1];
	      }
	    }
    	else if (k[h] == 1 || k[h] == 3) {
	      for (int ij = 1; ij < 7; ij++) {
		      for (int i = 0; i < 3; i++) {
    		    vec1[i] = cell2[0][i];
    		  }
      		mat_pow(rotmf,rotm1,ij-1);
      		mat_vec_mul(vec2,vec1,rotmf);
      		for (int i = 0; i<3; i++) {
    		    locs[i] = llocs[y-1][i][h]+vec2[i];
    		    locs0[i] = llocs[y-1][i][h];
    		  }
		      R[ij] = R[ij-1]+r2*exp(-(get_pot(locs,llocs,chg,IT,h,y)-get_pot(locs0,llocs,chg,IT,h,y))/kb/T);
	      }
	      for (int ij = 7; ij < 10; ij++) {
      		for (int i = 0; i < 3; i++) {
  		      vec1[i] = spos_Si[k[h]][i]-spos_Si[k[h]+1][i];
		      }
      		mat_pow(rotmf,rotm2,ij-7);
      		mat_mul(mf,cell2,rotmf);
      		mat_vec_mul(vec2,vec1,mf);
      		for (int i = 0; i < 3; i++) {
    		    locs[i] = llocs[y-1][i][h]+vec2[i];
    		    locs0[i] = llocs[y-1][i][h];
    		  }
        	R[ij] = R[ij-1]+r4*exp(-(get_pot(locs,llocs,chg,IT,h,y)-get_pot(locs0,llocs,chg,IT,h,y))/kb/T);
	      }
  	    for (int ij = 10; ij < 13; ij++) {
      		for (int i = 0; i < 3; i++) {
    		    vec1[i] = spos_Si[k[h]+2][i]-spos_Si[k[h]+1][i];
    		  }
          mat_pow(rotmf,rotm2,ij-10);
          mat_mul(mf,cell2,rotmf);
          mat_vec_mul(vec2,vec1,mf);
          for (int i = 0; i < 3; i++) {
    		    locs[i] = llocs[y-1][i][h]+vec2[i];
    		    locs0[i] = llocs[y-1][i][h];
    		  }
      		R[ij] = R[ij-1]+r3*exp(-(get_pot(locs,llocs,chg,IT,h,y)-get_pot(locs0,llocs,chg,IT,h,y))/kb/T);
	      }
	      for (int ij = 13; ij < 17; ij++) {
      		R[ij] = r5+R[ij-1];
	      }
	    }
	    else if (k[h] == 4) {
        continue;
	    }
    	Q = R[16];
    	u1 = (1.0*((rand()%(RAND_MAX-1))+1))/(1.0*RAND_MAX);
    	test = Q*u1;
    	L = 0;
    	E = 15;
    	for (int j = 0; j < 17; j++) {
        if (L > E) {
	        printf("\nunsuccessful\n");
        }
        m = (int) floor((E+L)/2);
        if (R[m] < test && R[m+1] < test) {
	        L = m+1;
	      }
        else if (R[m]>test && R[m+1] > test) {
		     E = m-1;
	      }
        else {
  	      break;
	      }
	    }
      if (m < 6) {
  	    for (int i = 0; i < 3; i++) {
  		    vec1[i] = cell2[0][i];
  		  }
        mat_pow(rotmf,rotm1,m);
        mat_vec_mul(vec2,vec1,rotmf);
        for (int i = 0; i < 3; i++) {
		      llocs[y][i][h] = llocs[y-1][i][h]+vec2[i];
	      }
	    }
      else if (m < 9) {
	      for (int i = 0; i < 3; i++) {
		      vec1[i] = spos_Si[k[h]][i]-spos_Si[k[h]+1][i];
	      }
  	    mat_pow(rotmf,rotm2,m-6);
  	    mat_mul(mf,cell2,rotmf);
  	    mat_vec_mul(vec2,vec1,mf);
  	    for (int i = 0; i < 3; i++) {
		      llocs[y][i][h] = llocs[y-1][i][h]+vec2[i];
	      }
	      k[h] = (k[h]+3)%4;
	    }
    	else if (m < 12) {
  	    for (int i = 0; i < 3; i++) {
		      vec1[i] = spos_Si[k[h]+2][i]-spos_Si[k[h]+1][i];
	      }
  	    mat_pow(rotmf,rotm2,m-9);
  	    mat_mul(mf,cell2,rotmf);
  	    mat_vec_mul(vec2,vec1,mf);
  	    for (int i = 0; i < 3; i++) {
		      llocs[y][i][h] = llocs[y-1][i][h]+vec2[i];
	      }
	      k[h] = (k[h]+1)%4;
	    }
	    else if (m > 11) {
	    for (int i = y; i < N; i++) {
    		llocs[i][0][h] = llocs[y-1][0][h];
    		llocs[i][1][h] = llocs[y-1][1][h];
    		llocs[i][2][h] = llocs[y-1][2][h];
	    }
	    k[h] = 4;
	    chg[h] = 2;
	    itt = (itt-1) < 1? (switcher = 1) : itt;  
	    itt = (itt-1) > 1? (itt-1): 1;
	  }
	  else {
	    printf("\nunsuccessful\n");
	  }
  	u1 = (1.0*((rand()%(RAND_MAX-1))+1))/(1.0*RAND_MAX);
  	t = t +1/Q*log(1/u1)/itt;
  }
  if (switcher == 1) {
	  break;
  }
  tt = tt + t;
}
  /*  for (int i = 0; i < IT; i++)
    {
      if (llocs[5][0][i] != llocs[6][0][i])
	{
	  printf("%Lf %Lf\n",llocs[5][0][i]*1e9,llocs[5][1][i]*1e9);
	}
    } 
}
*/
 
int main(int argc, char** argv)
{
  time_t t;
      
  long double cell2[3][3] = {
                              {nnd, 0*a, 0*a},
                              {-nnd/2,nnd/2*sqrt(3), 0*a},
                              {0*1.0, 0*1.0, 10.086*c*nnd/3.078*a/2.57218587467527*2.51866888630220}
                            };
  /* Intializes random number generator */
  srand((unsigned) time(&t));

  for (int i = 0; i < ITER; i++) {
    kMC(cell2,spos_Si);
  }
  
  return 0;
}
