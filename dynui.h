/* DYNUI.H - Header file for interace of user interface to simulation.
   Copyright 2020 by Radford M. Neal
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <SFML/Graphics.h>
#include <SFML/System.h>


/* CONTROL AREA PARAMETERS. */

static int c_height = 22;	/* c_height-4 must be divisible by 2, and
                                   c_height-10 must be divisible by 3 */

#define N_SPEEDS 4		/* Number speed settings */
#define N_SCALES 4		/* Number of zoom settings */


/* WINDOW STATE. */

struct window_state
{
  /* The following fields the must be set before calling dynui_window 
     (except width and height need not be set if full_screen is 1). */

  char *title;			/* Title for window */
  int running;			/* Should simulation be running? */
  int full_screen;		/* Is this a full-screen window? */
  int width;			/* Width of window */
  int height;			/* Height of window */

  /* The following fields should be used (but not changed) in the 
     dynui_view function provided by the application. */

  sfVector2f offset;		/* Offset of centre of view area */
  double scale;			/* Scaling factor for view area */

  /* The remaining fields should not be referenced by the application. */

  sfRenderWindow *window;	/* SFML window */
  sfFont *font;			/* Font used for text items */

  int view_area_pressed;	/* Mouse was pressed in view area */
  sfVector2i view_pressed;	/* Coordinates of press in view area */
  int control_pressed;		/* Mouse was pressed in control region */

  sfRectangleShape *boundary;	/* Line (rect) separating controls from view */
  sfRectangleShape *controls;	/* Rectangle in which controls reside */

  sfVertexArray *run_button;	/* Button to let simulation run */
  sfVertexArray *pause_button;	/* Button to pause simulation (replaces run) */
  sfCircleShape *speeds[N_SPEEDS];    /* Speed control buttons */
  sfRectangleShape *scales[N_SCALES]; /* Zoom buttons */
  sfVertexArray *save_button;	/* Button to ask application to save state */
  sfVertexArray *full_button;	/* Button to enter or exit full screen mode */

  sfText *sim_time_display;	/* Textual display of simulation time */

  sfClock *clock;		/* Clock used to control speed */
  double start_real_time;	/* Real elapsed time from start of run */
  double start_sim_time;	/* Simulation time from start of run */
  double sim_speed;		/* Speed of simulation */
  int running_behind;		/* Was simulation too slow for desired speed? */

  int exit;			/* Should this window be closed? */
};


/* STATE OF DYNAMIC SIMULATION. */

struct dynamic_state
{ 
  double sim_time;		/* Current simulation time */

  void *i;			/* Additional information */
};


/* PROCEDURE USED BY THE APPLICATION TO CREATE AND HANDLE A WINDOW. */

extern void dynui_window (struct dynamic_state *ds, struct window_state *ws);


/* PROCEDURES PROVIDED BY THE APPLICATION. */

extern void dynui_advance (struct dynamic_state *ds);
extern void dynui_view (struct dynamic_state *ds, struct window_state *window);
extern int dynui_save (struct dynamic_state *ds);
extern void dynui_terminate (void);
