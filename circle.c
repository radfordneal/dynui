/* CIRCLE.C - Test dynamic program that moves dot in a circle.
   Copyright 2020 by Radford M. Neal
 */

#include "dynui.h"


/* MAIN PROGRAM. */

int main (int argc, char **argv)
{
  /* Allocate and initialize dynamic state. */

  struct dynamic_state ds;

  ds.i = calloc (1, sizeof (double));
  if (ds.i == NULL) exit(3);
  *(double*)ds.i = 0;

  ds.sim_time = 0;
  ds.sim_info = "";

  /* Allocate window state and initialize application-specified fields. */

  struct window_state ws;

  ws.width = 600;
  ws.height = 500;
  ws.title = "circle";
  ws.running = 0;
  ws.full_screen = 0;

  /* Create window and let user interact with it, with simulation running
     when requested. */

  dynui_window (&ds, &ws);

  dynui_terminate();
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

static int dot_radius = 4;

void dynui_view (struct dynamic_state *ds, struct window_state *ws)
{ 
  sfVector2f origin = { ws->scale*ws->offset.x + ws->width/2, 
                        ws->scale*ws->offset.y + ws->height/2 };

  sfVector2f dot_origin = { ws->scale*dot_radius, ws->scale*dot_radius };

  sfCircleShape *dot = sfCircleShape_create();
  sfCircleShape_setOrigin (dot, dot_origin);
  sfCircleShape_setRadius (dot, ws->scale*dot_radius);

  sfCircleShape_setFillColor (dot, sfWhite);
  sfCircleShape_setPosition (dot, origin);

  sfRenderWindow_drawCircleShape (ws->window, dot, NULL);

  double theta = *(double*)ds->i;
  sfVector2f pos = { origin.x + ws->scale*100*cos(theta), 
                     origin.y + ws->scale*100*sin(theta) };

  sfCircleShape_setFillColor (dot, sfRed);
  sfCircleShape_setPosition (dot, pos);

  sfRenderWindow_drawCircleShape (ws->window, dot, NULL);

  sfCircleShape_destroy(dot);
}


/* SAVE STATE.  Returns 1 if succesful, 0 if not. */

int dynui_save (struct dynamic_state *ds)
{ return 0;
}


/* TERMINATE PROGRAM. */

void dynui_terminate()
{ exit(0);
}
