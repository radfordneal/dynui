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
  double sim_time;		/* Current simulation time */
  double sim_speed;		/* Speed of simulation */

  void *i;			/* Additional information */
};


/* PROCEDURES IMPLEMENTING THE DYNAMIC SIMULATION. */

extern void dynui_start (struct dynamic_state *ds, int argc, char **argv);
extern void dynui_destroy (struct dynamic_state *ds);
extern void dynui_advance (struct dynamic_state *ds);
extern void dynui_view (struct dynamic_state *ds, sfRenderWindow *window,
                        int width, int height);
