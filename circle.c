/* CIRCLE.C - Test dynamic program that moves dot in a circle.
   Copyright 2020 by Radford M. Neal
 */

#include "dynui.h"


/* INITIALIZE THE SIMULATION. */

void dynui_start (struct dynamic_state *ds, int argc, char **argv)
{
  ds->i = calloc (1, sizeof (double));
  if (ds->i == NULL) exit(3);

  *(double*)ds->i = 0;
}


/* RECOVER MEMORY FOR SIMULATION. */

void dynui_destroy (struct dynamic_state *ds)
{
  free (ds->i);
}


/* ADVANCE SIMULATION BY A SMALL TIME AMOUNT. */

static double delta_t = 0.02;

void dynui_advance (struct dynamic_state *ds)
{
  *(double*)ds->i = fmod (*(double*)ds->i + delta_t, 2*M_PI);
  ds->sim_time += delta_t;
  // sfSleep (sfSeconds (0.006));  /* for testing running behind indicator */
}


/* DRAW VIEW OF SIMULATION. */

void dynui_view (struct dynamic_state *ds, struct window_state *ws)
{ 
  int width = ws->width;
  int height = ws->height - c_height;

  sfCircleShape *dot = sfCircleShape_create();

  sfVector2f origin = { ws->offset.x + width/2 - 2, 
                        ws->offset.y + height/2 - 2 };
  sfCircleShape_setRadius (dot, 4);
  sfCircleShape_setFillColor (dot, sfWhite);
  sfCircleShape_setPosition (dot, origin);

  sfRenderWindow_drawCircleShape (ws->window, dot, NULL);

  double theta = *(double*)ds->i;
  sfVector2f pos = { ws->offset.x + width/2 + 100*cos(theta) - 2, 
                     ws->offset.y + height/2 + 100*sin(theta) - 2 };

  sfCircleShape_setRadius (dot, 4);
  sfCircleShape_setFillColor (dot, sfRed);
  sfCircleShape_setPosition (dot, pos);

  sfRenderWindow_drawCircleShape (ws->window, dot, NULL);

  sfCircleShape_destroy(dot);
}
