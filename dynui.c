/* DYNUI.C - User interface for dynamic simulations.
   Copyright 2020 by Radford M. Neal
 */

#include "dynui.h"


/* ZERO VECTORS.  For convenience. */

static sfVector2f zero_vector_f = { 0, 0 };
static sfVector2i zero_vector_i = { 0, 0 };


/* MAIN PROGRAM. */

static void dynui_window (struct dynamic_state *ds, struct window_state *ws);

int main (int argc, char **argv)
{
  /* Create dynamic state and start simulation (initially paused). */

  struct dynamic_state ds;
  ds.sim_time = 0;

  dynui_start (&ds, argc, argv);

  /* Create and interact with window for display and control of simulation. */

  struct window_state ws;
  ws.width = 600;
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
static void mouse_release (struct dynamic_state *ds, struct window_state *ws,
                           int x, int y);

static void dynui_window (struct dynamic_state *ds, struct window_state *ws)
{
  int i;

  /* Create the window, with specified size and title. */

  sfVideoMode mode = {ws->width, ws->height, 32};
  ws->window = sfRenderWindow_create(mode, ws->title, sfClose|sfTitlebar, NULL);
  if (ws->window == NULL) exit(1);

  sfVector2i window_pos = { 10, 10 };
  sfRenderWindow_setPosition (ws->window, window_pos);

  sfRenderWindow_setVerticalSyncEnabled (ws->window, 1);

  /* Load a font to use for textual displays. */

  ws->font = sfFont_createFromFile("ubuntu-fonts/UbuntuMono-R.ttf");
  if (ws->font == NULL) exit(2);

  /* Create a clock used in controlling speed. */

  ws->clock = sfClock_create();

  /* Set initial values for state. */

  ws->offset = zero_vector_f;
  ws->scale = 1;
  ws->sim_speed = 1;
  ws->view_area_pressed = 0;
  ws->control_pressed = 0;
  ws->running = 0;
  ws->running_behind = 0;

  /* Create user control items. */

  create_controls (ws);

  /* The main user interaction loop. */

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
      { mouse_release (ds, ws, event.mouseButton.x, event.mouseButton.y);
      }
    }

    /* Shift view during press of the mouse in the view area. */

    if (ws->view_area_pressed)
    { sfVector2i p = sfMouse_getPositionRenderWindow (ws->window);
      ws->offset.x += (p.x - ws->view_pressed.x) / ws->scale;
      ws->offset.y += (p.y - ws->view_pressed.y) / ws->scale;
      ws->view_pressed = p;
    }

    /* Run some more of the simulation (unless we're shifting the view).
       We run only as much as needed to keep up with the desired
       speed.  But if the simulation is slow, we give up if we're
       making negligible progress at matching the desired speed (since
       otherwise we'd never finish here). */

    if (!ws->view_area_pressed)
    {
      if (ws->running_behind)
      { for (i = 0; i < N_SPEEDS; i++)
        { sfCircleShape_setFillColor (ws->speeds[i], ws->sim_speed == 1<<i ?
                                      sfWhite : sfColor_fromRGB(150, 150, 150));
        }
      }

      if (ws->running)
      { double current_time;
        double gap, oldgap;
        oldgap = HUGE_VAL;
        ws->running_behind = 0;
        for (;;) 
        { current_time 
            = 1e-6 * sfTime_asMicroseconds (sfClock_getElapsedTime (ws->clock));
          gap = ws->sim_speed * (current_time - ws->start_real_time) 
                 - (ds->sim_time - ws->start_sim_time);
          if (gap <= 0)
          { break;
          }
          if (gap > 0.9*oldgap) 
          { ws->running_behind = 1;
            break;
          }
          dynui_advance (ds);
          oldgap = gap;
        }
      }

      if (ws->running_behind)
      { for (i = 0; i < N_SPEEDS; i++)
        { sfCircleShape_setFillColor (ws->speeds[i], ws->sim_speed == 1<<i ?
                                      sfRed : sfColor_fromRGB (150, 150, 150));
        }
      }

      char s[12];
      sprintf (s, "%10.1f", ds->sim_time);
      sfText_setString (ws->sim_time_display, s);
    }
          
    /* Redraw the window. */

    sfRenderWindow_clear (ws->window, sfBlack);

    dynui_view (ds, ws);
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

static void create_controls (struct window_state *ws)
{
  int x, y, h, w, i;
  sfVector2f p;

  sfVertex v;
  v.position = zero_vector_f;
  v.color = sfWhite;
  v.texCoords = zero_vector_f;

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

  x = 5;
  y = ws->height - c_height + 2;
  h = c_height - 4;
  w = c_height - 10;

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

  /* Speed controls. */

  x = 30;
  y = ws->height - c_height + 6;
  h = c_height - 12;

  for (i = 0; i < N_SPEEDS; i++)
  { ws->speeds[i] = sfCircleShape_create();
    p.x = x + (h+3)*i;
    p.y = y;
    sfCircleShape_setPosition (ws->speeds[i], p);
    sfCircleShape_setRadius (ws->speeds[i], h/2);
    sfCircleShape_setFillColor (ws->speeds[i], ws->sim_speed == 1<<i ? sfWhite :
                                               sfColor_fromRGB (150, 150, 150));
  }

  /* Zoom controls. */

  sfVector2f s = { h, h };
  x = ws->width/2 - (h+3)*N_SCALES/2;
  y = ws->height - c_height + 6;
  h = c_height - 12;

  for (i = 0; i < N_SPEEDS; i++)
  { ws->scales[i] = sfRectangleShape_create();
    p.x = x + (h+3)*i;
    p.y = y;
    sfRectangleShape_setPosition (ws->scales[i], p);
    sfRectangleShape_setSize (ws->scales[i], s);
    sfRectangleShape_setFillColor (ws->scales[i], ws->scale == 1<<i ? sfWhite :
                                   sfColor_fromRGB (150, 150, 150));
  }

  /* Text of simulation time. */

  ws->sim_time_display = sfText_create();
  sfVector2f t = { ws->width-5.2*(c_height-4), ws->height-c_height-1 }; 
  sfText_setPosition (ws->sim_time_display, t);
  sfText_setFont (ws->sim_time_display, ws->font);
  sfText_setCharacterSize (ws->sim_time_display, c_height-4);
  sfText_setString (ws->sim_time_display, "");
}


/* DRAW USER CONTROLS. */

static void draw_controls (struct window_state *ws) 
{
  int i;

  sfRenderWindow_drawRectangleShape (ws->window, ws->boundary, NULL);
  sfRenderWindow_drawRectangleShape (ws->window, ws->controls, NULL);

  if (ws->running)
  { sfRenderWindow_drawVertexArray (ws->window, ws->pause_button, NULL);
  }
  else
  { sfRenderWindow_drawVertexArray (ws->window, ws->run_button, NULL);
  }

  for (i = 0; i < N_SPEEDS; i++)
  { sfRenderWindow_drawCircleShape (ws->window, ws->speeds[i], NULL);
  }

  for (i = 0; i < N_SCALES; i++)
  { sfRenderWindow_drawRectangleShape (ws->window, ws->scales[i], NULL);
  }

  sfRenderWindow_drawText (ws->window, ws->sim_time_display, NULL);
}


/* DESTROY USER CONTROLS. */

static void destroy_controls (struct window_state *ws)
{ int i;
  sfRectangleShape_destroy (ws->boundary);
  sfRectangleShape_destroy (ws->controls);
  for (i = 0; i < N_SPEEDS; i++)
  { sfCircleShape_destroy (ws->speeds[i]);
  }
}


/* HANDLE A MOUSE PRESS EVENT.  

   For a press in the view area, we record that the mouse was pressed
   there, so shifting of the view will commence.

   For a press in the control area, we just record that it happened.  
   Nothing is done until/unless the release event happens, also in
   the control area. */

static void mouse_press (struct window_state *ws, int x, int y) 
{ if (y > ws->height - c_height)
  { ws->control_pressed = 1;
  }
  else
  { ws->view_area_pressed = 1;
    ws->view_pressed.x = x;
    ws->view_pressed.y = y;
  }
}


/* HANDLE A MOUSE RELEASE EVENT.  

   For a release after a press in the view area, we record that it is
   no longer pressed, so shifting of the view will stop.

   For a release in the control area, something may actually be done now, 
   if the press was also in the control area. */

static void mouse_release (struct dynamic_state *ds, struct window_state *ws,
                           int x, int y)
{ 
  int in_control_area = ws->control_pressed && y > ws->height - c_height;

  if (ws->view_area_pressed)
  { ws->view_area_pressed = 0;
    goto reset_running;
  }

  if (in_control_area)
  {
    sfFloatRect bounds;
    int i, j;

    ws->control_pressed = 0;

    bounds = sfVertexArray_getBounds(ws->pause_button);
    if (sfFloatRect_contains(&bounds,x,y))
    { ws->running = !ws->running;
      goto reset_running;
    }

    for (i = 0; i < N_SPEEDS; i++)
    { bounds = sfCircleShape_getGlobalBounds(ws->speeds[i]);
      if (sfFloatRect_contains(&bounds,x,y))
      { ws->sim_speed = 1<<i;
        for (j = 0; j < N_SPEEDS; j++)
        { sfCircleShape_setFillColor (ws->speeds[j], j==i ? sfWhite : 
                                       sfColor_fromRGB (150, 150, 150));
        }
        goto reset_running;
      }
    }

    for (i = 0; i < N_SCALES; i++)
    { bounds = sfRectangleShape_getGlobalBounds(ws->scales[i]);
      if (sfFloatRect_contains(&bounds,x,y))
      { ws->scale = 1<<i;
        for (j = 0; j < N_SCALES; j++)
        { sfRectangleShape_setFillColor (ws->scales[j], j==i ? sfWhite : 
                                        sfColor_fromRGB (150, 150, 150));
        }
        goto reset_running;
      }
    }
  }

  return;

reset_running:
  if (ws->running)
  { ws->start_real_time = sfTime_asSeconds(sfClock_getElapsedTime(ws->clock));
    ws->start_sim_time = ds->sim_time;
  }

  return;
}
