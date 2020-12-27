/* DYNUI - User interface for dynamic simulations.
   Copyright 2020 by Radford M. Neal
 */

#include "dynamic.h"


/* WINDOW STATE. */

struct window_state
{
  int width;			/* Width of window */
  int height;			/* Height of window */
  char *title;			/* Title for window */

  sfRenderWindow *window;	/* SFML window */
  sfFont *font;			/* Font used for text items */

  int control_pressed;		/* Control where mouse was pressed, 0 = none */
  int running;			/* Is simulation running? */

  sfRectangleShape *boundary;	/* Line (rect) separating controls from view */
  sfRectangleShape *controls;	/* Rectangle in which controls reside */

  sfVertexArray *run_button;	/* Button to let simulation run */
  sfVertexArray *pause_button;	/* Button to pause simulation (replaces run) */
};


/* MAIN PROGRAM. */

static void dynui_window (struct dynamic_state *ds, struct window_state *ws);

int main (int argc, char **argv)
{
  /* Create dynamic state and start simulation (initially paused). */

  struct dynamic_state ds;
  ds.sim_time = 0;
  ds.sim_speed = 1;
  ds.clock = sfClock_create();
  ds.real_time = sfClock_getElapsedTime (ds.clock);

  dynui_start (&ds, argc, argv);

  /* Create and interact with window for display and control of simulation. */

  struct window_state ws;
  ws.width = 400;
  ws.height = 500;
  ws.title = "DYNUI TEST";

  dynui_window (&ds, &ws);

  exit(0);
}


/* HANDLE INTERACTION WITH A WINDOW. */

static void create_controls (struct window_state *ws);
static void draw_controls (struct window_state *ws);
static void destroy_controls (struct window_state *ws);

static void mouse_press (struct window_state *ws, int x, int y);
static void mouse_release (struct window_state *ws, int x, int y);

static void dynui_window (struct dynamic_state *ds, struct window_state *ws)
{
  /* Create the window, with specified size and title. */

  sfVideoMode mode = {ws->width, ws->height, 32};
  ws->window = sfRenderWindow_create (mode, ws->title, sfClose, NULL);
  if (!ws->window) exit(1);

  sfVector2i window_pos = { 10, 10 };
  sfRenderWindow_setPosition (ws->window, window_pos);

  /* Load a font to use for textual displays. */

  ws->font = sfFont_createFromFile("ubuntu-fonts/Ubuntu-R.ttf");
  if (!ws->font) exit(2);

  /* Create user control items. */

  create_controls (ws);

  /* The main user interaction loop. */

  ws->control_pressed = 0;
  ws->running = 0;

  while (sfRenderWindow_isOpen(ws->window))
  {
    /* Process events */

    sfEvent event;
    while (sfRenderWindow_pollEvent (ws->window, &event))
    {
      if (event.type == sfEvtClosed)
      { sfRenderWindow_close (ws->window);
      }
      else if (event.type == sfEvtMouseButtonPressed)
      { mouse_press (ws, event.mouseButton.x, event.mouseButton.y);
      }
      else if (event.type == sfEvtMouseButtonReleased)
      { mouse_release (ws, event.mouseButton.x, event.mouseButton.y);
      }
    }

    /* Redraw the window. */

    sfRenderWindow_clear (ws->window, sfBlack);
    draw_controls (ws);

    /* Render the window onto the display. */

    sfRenderWindow_display (ws->window);
  }

  /* Clean up resources used. */

  destroy_controls (ws);
  sfFont_destroy (ws->font);
  sfRenderWindow_destroy (ws->window);
}


/* CREATE USER CONTROLS. */

static int c_height = 22;	/* c_height-4 must be divisible by 2, and
                                   c_height-10 must be divisible by 3 */

static sfVector2f zero_vector = { 0, 0 };

static void create_controls (struct window_state *ws)
{
  /* Control area and boundary. */

  ws->boundary = sfRectangleShape_create();

  sfVector2f boundary_pos = { 0, ws->height - c_height - 2 };
  sfRectangleShape_setPosition (ws->boundary, boundary_pos);

  sfVector2f boundary_size = { ws->width, 2 };
  sfRectangleShape_setSize (ws->boundary, boundary_size);

  sfRectangleShape_setFillColor (ws->boundary, sfWhite);

  ws->controls = sfRectangleShape_create();

  sfVector2f controls_pos = { 0, ws->height - c_height };
  sfRectangleShape_setPosition (ws->controls, controls_pos);

  sfVector2f controls_size = { ws->width, c_height };
  sfRectangleShape_setSize (ws->controls, controls_size);

  sfRectangleShape_setFillColor (ws->controls, sfBlack);

  /* Run and pause buttons. */

  int x = 5;
  int y = ws->height - c_height + 2;
  int h = c_height - 4;
  int w = c_height - 10;

  sfVertex v;
  v.position = zero_vector;
  v.color = sfWhite;
  v.texCoords = zero_vector;

  ws->run_button = sfVertexArray_create();
  v.position.x = x;
  v.position.y = y;
  sfVertexArray_append (ws->run_button, v);
  v.position.y = y+h;
  sfVertexArray_append (ws->run_button, v);
  v.position.x = x+h/2;
  v.position.y = y+h/2;
  sfVertexArray_append (ws->run_button, v);
  sfVertexArray_setPrimitiveType (ws->run_button, sfTriangles);

  ws->pause_button = sfVertexArray_create();
  v.position.x = x;
  v.position.y = y;
  sfVertexArray_append (ws->pause_button, v);
  v.position.y = y+h;
  sfVertexArray_append (ws->pause_button, v);
  v.position.x = x+w/3;
  sfVertexArray_append (ws->pause_button, v);
  v.position.y = y;
  sfVertexArray_append (ws->pause_button, v);
  v.position.x = x+2*w/3;
  v.position.y = y;
  sfVertexArray_append (ws->pause_button, v);
  v.position.y = y+h;
  sfVertexArray_append (ws->pause_button, v);
  v.position.x = x+w;
  sfVertexArray_append (ws->pause_button, v);
  v.position.y = y;
  sfVertexArray_append (ws->pause_button, v);
  sfVertexArray_setPrimitiveType (ws->pause_button, sfQuads);
}


/* DRAW USER CONTROLS. */

static void draw_controls (struct window_state *ws) 
{
  sfRenderWindow_drawRectangleShape (ws->window, ws->boundary, NULL);
  sfRenderWindow_drawRectangleShape (ws->window, ws->controls, NULL);

  if (ws->running)
  { sfRenderWindow_drawVertexArray (ws->window, ws->pause_button, NULL);
  }
  else
  { sfRenderWindow_drawVertexArray (ws->window, ws->run_button, NULL);
  }
}


/* DESTROY USER CONTROLS. */

static void destroy_controls (struct window_state *ws)
{ sfRectangleShape_destroy (ws->boundary);
  sfRectangleShape_destroy (ws->controls);
}


static void mouse_press (struct window_state *ws, int x, int y) {}
static void mouse_release (struct window_state *ws, int x, int y) {}

