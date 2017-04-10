/*****************************************************************************/
/*
  controller.c: control strategy.
*/
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "matrix3.h"
#include "main.h"
#include "main2.h"
#include "sdfast/alien.h"
#include "sdfast/lander.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* call this once to do one-time initialize: allocate memeory etc. */

void init_controller( SIM *s )
{

}

/*****************************************************************************/
/* call this many times to restart a controller */

void reinit_controller( SIM *s )
{

}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
// static FILE *fid=NULL;
int controller( SIM *s )
{
	int i;
	static int count = 0;
	double k_x = 2.0;//1.0;
	double b_x = 2.5;//2.0;

	double k_r = 10;//1.0;
	double b_r = 3; //2.0;

	double q_minus[N_Q];
	double q_diff[N_Q];
  
	 for ( i = 0; i < 3; i++ )
    {
      s->lander_thrust_world[i] = 
		k_x*( s->lander_x_d[i] - s->lander_x[i] ) +
		b_x*( - s->lander_xd[i] );
    }
  multiply_transpose_m3_v3( s->lander_r, s->lander_thrust_world, 
			    s->lander_thrust );

  // lander orientation control
  // Attempt to do PD control, but orientations not vectors, so complicated
  // "subtract quaternions"
  invert_q( s->lander_q, q_minus );
  compose_q( q_minus, s->lander_q_d, q_diff );
  q_to_rotvec( q_diff, s->rotvec ); 
  printf( "%f %f %f %f\n", s->lander_q_d[0], s->lander_q_d[1], s->lander_q_d[2], s->lander_q_d[3] );
  // printf( "%g %g %g\n", s->rotvec[0], s->rotvec[1], s->rotvec[2] );

  // PD servo for orientation
  // w and rotvec are in body coordinates.
  for ( i = 0; i < N_XYZ; i++ )
    {
      s->lander_torque[i] = (k_r*s->rotvec[i] - b_r*s->lander_w[i]);
    }
  multiply_m3_v3( s->lander_r, s->lander_torque, s->lander_torque_world );

  // To get better control should compute desired acceleration and use
  // inverse dynamics to compute torque.

  return 0;
}

/*****************************************************************************/

