/* DYNAMIC.H - Header file for interace of user interface to simulation.
   Copyright 2020 by Radford M. Neal
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <SFML/Graphics.h>
#include <SFML/System.h>


/* STATE OF DYNAMIC SIMULATION. */

struct dynamic_state
{ 
  sfClock *clock;		/* Clock used to control speed */

  sfTime real_time;		/* Real elapsed time after last update */
  double sim_time;		/* Simulation time after last update */
  double sim_speed;		/* Speed of simulation */

  void *i;			/* Additional information */
};


/* PROCEDURES IMPLEMENTING THE DYNAMIC SIMULATION. */

extern void dynui_start (struct dynamic_state *ds, int argc, char **argv);
