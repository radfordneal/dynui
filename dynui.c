/* DYNUI.C - User interface for dynamic simulations.
   Copyright 2021 by Radford M. Neal

    This program is free software: you can redistribute it and/or
    modify it under the terms of version 3 of the GNU General Public
    License as published by the Free Software Foundation.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "dynui.h"


/* ZERO VECTORS.  For convenience. */

static sfVector2f zero_vector_f = { 0, 0 };
static sfVector2i zero_vector_i = { 0, 0 };


/* CREATE AND INTERACT WITH A WINDOW ON A DYNAMIC SIMULATION. */

static void create_controls (struct window_state *ws);
static void draw_controls (struct window_state *ws);
static void destroy_controls (struct window_state *ws);

static void mouse_press (struct window_state *ws, int x, int y);
static void mouse_release (struct dynamic_state *ds, struct window_state *ws,
                           int x, int y);

static void set_start_time (struct dynamic_state *ds, struct window_state *ws);
static void append_vertex (sfVertexArray *a, int x, int y, sfColor c);


void dynui_window (struct dynamic_state *ds, struct window_state *ws)
{
  int i;

  /* Create the window, with specified size & title, or instead full screen. */

  if (ws->full_screen)
  { ws->window = sfRenderWindow_create (sfVideoMode_getDesktopMode(), ws->title,
                                        sfFullscreen, NULL);
  }
  else
  { sfVideoMode mode = { ws->width, ws->height, 32 };
    ws->window = sfRenderWindow_create (mode, ws->title, 
                                        sfClose | sfTitlebar, NULL);
  }

  if (ws->window == NULL) exit(1);

  if (ws->full_screen)
  { sfVector2u size = sfRenderWindow_getSize (ws->window);
    ws->width = size.x;
    ws->height = size.y;
  }

  if (!ws->full_screen)
  { sfVector2i window_pos = { 10, 10 };
    sfRenderWindow_setPosition (ws->window, window_pos);
  }

  sfRenderWindow_setVerticalSyncEnabled (ws->window, 1);

  /* Load a font to use for textual displays. */

  ws->font = sfFont_createFromFile("ubuntu-fonts/UbuntuMono-R.ttf");
  if (ws->font == NULL) exit(2);

  /* Create a clock used in controlling speed. */

  ws->clock = sfClock_create();

  /* Set initial values for state.  Some are assumed already set if we're
     entering full screen mode from window mode (full_screen set to -1). */

  if (ws->full_screen != -1)
  { ws->offset = zero_vector_f;
    ws->scale = 1;
    ws->sim_speed = 1;
  }

  ws->view_area_pressed = 0;
  ws->control_pressed = 0;
  ws->running_behind = 0;
  ws->exit = 0;

  /* Create user control items. */

  create_controls (ws);

  /* The main user interaction loop. */

  set_start_time (ds, ws);

  while (sfRenderWindow_isOpen(ws->window))
  {
    /* Process events */

    sfEvent event;
    while (sfRenderWindow_pollEvent (ws->window, &event))
    {
      if (event.type == sfEvtClosed)
      { ws->exit = 1;
        continue;
      }
      else if (event.type == sfEvtKeyPressed)
      { if (event.key.code == sfKeyC && event.key.control)
        { dynui_terminate();
        }
        else if (event.key.code == sfKeyEscape)
        { ws->exit = 1;
          continue;
        }
        else if (event.key.code == sfKeySpace)
        { ws->running = !ws->running;
          set_start_time (ds, ws);
          continue;
        }
      }
      else if (event.type == sfEvtMouseButtonPressed &&
               event.mouseButton.button == 0)
      { mouse_press (ws, event.mouseButton.x, event.mouseButton.y);
        continue;
      }
      else if (event.type == sfEvtMouseButtonReleased &&
               event.mouseButton.button == 0)
      { mouse_release (ds, ws, event.mouseButton.x, event.mouseButton.y);
        continue;
      }

      dynui_event (ds, ws, event);
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
        { sfCircleShape_setFillColor (ws->speeds[i], 
            ws->sim_speed*(1<<SLOW_SPEEDS) == 1<<i ? sfWhite 
             : sfColor_fromRGB(150, 150, 150));
        }
      }

      if (ws->running)
      { double current_time;
        double gap, oldgap;
        int advances;
        oldgap = HUGE_VAL;
        advances = 0;
        ws->running_behind = 0;
        for (;;) 
        { current_time 
            = 1e-6 * sfTime_asMicroseconds (sfClock_getElapsedTime (ws->clock));
          gap = ws->sim_speed * (current_time - ws->start_real_time) 
                 - (ds->sim_time - ws->start_sim_time);
          if (gap <= 0)
          { break;
          }
          if (gap > 1.1*oldgap || gap > 0.95*oldgap && advances > 1) 
          { ws->running_behind = 1;
            break;
          }
          dynui_advance (ds);
          advances += 1;
          oldgap = gap;
        }
      }

      if (ws->running_behind)
      { for (i = 0; i < N_SPEEDS; i++)
        { sfCircleShape_setFillColor (ws->speeds[i], 
            ws->sim_speed*(1<<SLOW_SPEEDS) == 1<<i ? sfRed
             : sfColor_fromRGB(150, 150, 150));
        }
      }

      char s[201];
      sprintf (s, "%10.1f", ds->sim_time);
      sfText_setString (ws->sim_time_display, s);
      strncpy (s, ds->sim_info, 200);
      s[200] = 0;
      if (ws->full_screen)
      { strncat (s, ds->sim_xinfo, 200-strlen(s));
      }
      sfText_setString (ws->sim_info_display, s);
    }
          
    /* Redraw the window. */

    sfRenderWindow_clear (ws->window, sfBlack);

    dynui_view (ds, ws);
    draw_controls (ws);

    /* Render the window onto the display. */

    sfRenderWindow_display (ws->window);

    /* Close window if user requested this. */

    if (ws->exit)
    { sfRenderWindow_close (ws->window);
    }
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
  append_vertex (ws->run_button, x,     y,     sfWhite);
  append_vertex (ws->run_button, x,     y+h,   sfWhite);
  append_vertex (ws->run_button, x+h/2, y+h/2, sfWhite);
  sfVertexArray_setPrimitiveType (ws->run_button, sfTriangles);

  ws->pause_button = sfVertexArray_create();
  append_vertex (ws->pause_button, x,       y,   sfWhite);
  append_vertex (ws->pause_button, x,       y+h, sfWhite);
  append_vertex (ws->pause_button, x+w/3,   y+h, sfWhite);
  append_vertex (ws->pause_button, x+w/3,   y,   sfWhite);
  append_vertex (ws->pause_button, x+2*w/3, y,   sfWhite);
  append_vertex (ws->pause_button, x+2*w/3, y+h, sfWhite);
  append_vertex (ws->pause_button, x+w,     y+h, sfWhite);
  append_vertex (ws->pause_button, x+w,     y,   sfWhite);
  sfVertexArray_setPrimitiveType (ws->pause_button, sfQuads);

  /* Speed controls. */

  x = 25;
  y = ws->height - c_height + 6;
  h = c_height - 12;

  for (i = 0; i < N_SPEEDS; i++)
  { ws->speeds[i] = sfCircleShape_create();
    p.x = x + (h+3)*i;
    p.y = y;
    sfCircleShape_setPosition (ws->speeds[i], p);
    sfCircleShape_setRadius (ws->speeds[i], h/2);
    sfCircleShape_setFillColor (ws->speeds[i], 
      ws->sim_speed*(1<<SLOW_SPEEDS) == 1<<i ? sfWhite 
       : sfColor_fromRGB (150, 150, 150));
    if (SLOW_SPEEDS > 0 && i == SLOW_SPEEDS)
    { ws->speed_one = sfVertexArray_create();
      append_vertex (ws->speed_one, p.x+h/2, p.y+h-1, sfBlack);
      append_vertex (ws->speed_one, p.x+h/2, p.y+1, sfBlack);
      sfVertexArray_setPrimitiveType (ws->speed_one, sfLines);
    }
  }

  /* Zoom controls. */

  sfVector2f s = { h-1, h };
  x = 35 + (h+3)*N_SPEEDS;
  y = ws->height - c_height + 6;
  h = c_height - 12;

  for (i = 0; i < N_SCALES; i++)
  { ws->scales[i] = sfRectangleShape_create();
    p.x = x + (h+3)*i;
    p.y = y;
    sfRectangleShape_setPosition (ws->scales[i], p);
    sfRectangleShape_setSize (ws->scales[i], s);
    sfRectangleShape_setFillColor (ws->scales[i], 
      ws->scale*(1<<SMALL_SCALES) == 1<<i ? sfWhite 
       : sfColor_fromRGB (150, 150, 150));
    if (SMALL_SCALES > 0 && i == SMALL_SCALES)
    { ws->scale_one = sfVertexArray_create();
      append_vertex (ws->scale_one, p.x+h/2, p.y+h-1, sfBlack);
      append_vertex (ws->scale_one, p.x+h/2, p.y+1, sfBlack);
      sfVertexArray_setPrimitiveType (ws->scale_one, sfLines);
    }
  }

  /* Text of simulation time. */

  ws->sim_time_display = sfText_create();
  sfVector2f t = { ws->width-5.2*(c_height-4)-2.7*c_height+6, 
                   ws->height-c_height-1 }; 
  sfText_setPosition (ws->sim_time_display, t);
  sfText_setFont (ws->sim_time_display, ws->font);
  sfText_setCharacterSize (ws->sim_time_display, c_height-4);
  sfText_setString (ws->sim_time_display, "");

  /* Text of simulation information. */

  ws->sim_info_display = sfText_create();
  sfVector2f ti = { 33 + (h+3)*(N_SPEEDS+N_SCALES),
                    ws->height-c_height-1 }; 
  sfText_setPosition (ws->sim_info_display, ti);
  sfText_setFont (ws->sim_info_display, ws->font);
  sfText_setCharacterSize (ws->sim_info_display, c_height-4);
  sfText_setString (ws->sim_info_display, "");

  /* Save button. */

  x = ws->width - 2.2*c_height + 3;
  y = ws->height - c_height + 3;
  h = c_height - 7;

  ws->save_button = sfVertexArray_create();
  append_vertex (ws->save_button, x+5,  y,     sfWhite);
  append_vertex (ws->save_button, x+5,  y+h,   sfWhite);
  append_vertex (ws->save_button, x+5,  y+h,   sfWhite);
  append_vertex (ws->save_button, x,    y+h-5, sfWhite);
  append_vertex (ws->save_button, x+5,  y+h,   sfWhite);
  append_vertex (ws->save_button, x+10, y+h-5, sfWhite);
  sfVertexArray_setPrimitiveType (ws->save_button, sfLines);

  /* Full screen enter/exit button. */

  x = ws->width - c_height;
  y = ws->height - c_height + 4;
  h = c_height - 8;

  ws->full_button = sfVertexArray_create();
    
  if (ws->full_screen)
  { append_vertex (ws->full_button, x,   y,   sfWhite);
    append_vertex (ws->full_button, x+h, y+h, sfWhite);
    append_vertex (ws->full_button, x,   y+h, sfWhite);
    append_vertex (ws->full_button, x+h, y,   sfWhite);
  }
  else
  { append_vertex (ws->full_button, x,       y+h,     sfWhite);
    append_vertex (ws->full_button, x,       y+h-h/2, sfWhite);
    append_vertex (ws->full_button, x,       y+h,     sfWhite);
    append_vertex (ws->full_button, x+h/2,   y+h,     sfWhite);
    append_vertex (ws->full_button, x,       y+h,     sfWhite);
    append_vertex (ws->full_button, x+h,     y,       sfWhite);
    append_vertex (ws->full_button, x+h,     y,       sfWhite);
    append_vertex (ws->full_button, x+h-h/2, y,       sfWhite);
    append_vertex (ws->full_button, x+h,     y,       sfWhite);
    append_vertex (ws->full_button, x+h,     y+h/2,   sfWhite);
  }
  sfVertexArray_setPrimitiveType (ws->full_button, sfLines);
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
  if (SLOW_SPEEDS > 0)
  { sfRenderWindow_drawVertexArray (ws->window, ws->speed_one, NULL);
  }

  for (i = 0; i < N_SCALES; i++)
  { sfRenderWindow_drawRectangleShape (ws->window, ws->scales[i], NULL);
  }
  if (SMALL_SCALES > 0)
  { sfRenderWindow_drawVertexArray (ws->window, ws->scale_one, NULL);
  }

  sfRenderWindow_drawVertexArray (ws->window, ws->save_button, NULL);

  sfRenderWindow_drawText (ws->window, ws->sim_time_display, NULL);
  sfRenderWindow_drawText (ws->window, ws->sim_info_display, NULL);

  sfRenderWindow_drawVertexArray (ws->window, ws->full_button, NULL);
}


/* DESTROY USER CONTROLS. */

static void destroy_controls (struct window_state *ws)
{ int i;
  sfRectangleShape_destroy (ws->boundary);
  sfRectangleShape_destroy (ws->controls);
  sfVertexArray_destroy (ws->pause_button);
  sfVertexArray_destroy (ws->run_button);
  for (i = 0; i < N_SPEEDS; i++)
  { sfCircleShape_destroy (ws->speeds[i]);
  }
  if (SLOW_SPEEDS > 0)
  { sfVertexArray_destroy (ws->speed_one);
  }
  for (i = 0; i < N_SCALES; i++)
  { sfRectangleShape_destroy (ws->scales[i]);
  }
  if (SMALL_SCALES > 0)
  { sfVertexArray_destroy (ws->scale_one);
  }
  sfText_destroy (ws->sim_time_display);
  sfText_destroy (ws->sim_info_display);
  sfVertexArray_destroy (ws->full_button);
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
   if the press was also in the control area.  

   Checks for whether the release was for a certain control looks slightly 
   outside the bounding box, to make missing less likely. */

static int in_bounds (sfFloatRect bounds, int x, int y)
{ bounds.left -= 3; bounds.top -= 3; 
  bounds.width += 3; bounds.height += 3;
  return sfFloatRect_contains (&bounds, x, y);
}

static void mouse_release (struct dynamic_state *ds, struct window_state *ws,
                           int x, int y)
{ 
  int in_control_area = ws->control_pressed && y > ws->height - c_height;

  if (ws->view_area_pressed)
  { ws->view_area_pressed = 0;
    set_start_time (ds, ws);
    return;
  }

  if (in_control_area)
  {
    int i, j;

    ws->control_pressed = 0;

    /* Run/pause button. */

    if (in_bounds (sfVertexArray_getBounds(ws->pause_button), x, y))
    { ws->running = !ws->running;
      set_start_time (ds, ws);
      return;
    }

    /* Save button. */

    if (in_bounds (sfVertexArray_getBounds(ws->save_button), x, y))
    { int success = dynui_save(ds);
      for (i = sfVertexArray_getVertexCount(ws->save_button)-1; i >= 0; i--)
      { sfVertexArray_getVertex(ws->save_button,i)->color 
          = success ? sfGreen : sfRed;
      }
      sfRenderWindow_drawVertexArray (ws->window, ws->save_button, NULL);
      sfRenderWindow_display (ws->window);
      sfSleep (sfSeconds(0.25));
      for (i = sfVertexArray_getVertexCount(ws->save_button)-1; i >= 0; i--)
      { sfVertexArray_getVertex(ws->save_button,i)->color = sfWhite;
      }
      set_start_time (ds, ws);
      return;
    }

    /* Speed controls. */

    for (i = 0; i < N_SPEEDS; i++)
    { if (in_bounds (sfCircleShape_getGlobalBounds(ws->speeds[i]), x, y))
      { ws->sim_speed = (double) (1<<i) / (1<<SLOW_SPEEDS);
        for (j = 0; j < N_SPEEDS; j++)
        { sfCircleShape_setFillColor (ws->speeds[j], 
            ws->sim_speed*(1<<SLOW_SPEEDS) == 1<<j ? sfWhite 
             : sfColor_fromRGB (150, 150, 150));
        }
        set_start_time (ds, ws);
        return;
      }
    }

    /* Zoom controls. */

    for (i = 0; i < N_SCALES; i++)
    { if (in_bounds (sfRectangleShape_getGlobalBounds(ws->scales[i]), x, y))
      { ws->scale = (double) (1<<i) / (1<<SMALL_SCALES);
        for (j = 0; j < N_SCALES; j++)
        { sfRectangleShape_setFillColor (ws->scales[j], 
            ws->scale*(1<<SMALL_SCALES) == 1<<j ? sfWhite 
             : sfColor_fromRGB (150, 150, 150));
        }
        set_start_time (ds, ws);
        return;
      }
    }

    /* Full screen / exit screen button. */

    if (in_bounds (sfVertexArray_getBounds(ws->full_button), x, y))
    { if (ws->full_screen)
      { ws->exit = 1;
      }
      else
      { 
        /* Create a full-screen window, and let it handle things until the
           user exits from it.  Transfer some settings from this window to 
           the full-screen one, and back. */

        struct window_state fws;
        fws.full_screen = -1;  /* full screen with some fields already set */
        fws.title = "";

        fws.running = ws->running;
        fws.offset = ws->offset;
        fws.scale = ws->scale;
        fws.sim_speed = ws->sim_speed;

        dynui_window (ds, &fws);

        ws->running = fws.running;
        ws->offset = fws.offset;
        ws->scale = fws.scale;
        ws->sim_speed = fws.sim_speed;

        for (j = 0; j < N_SPEEDS; j++)
        { sfCircleShape_setFillColor (ws->speeds[j], 
            ws->sim_speed == 1<<j ? sfWhite : sfColor_fromRGB (150, 150, 150));
        }
        for (j = 0; j < N_SCALES; j++)
        { sfRectangleShape_setFillColor (ws->scales[j], 
            ws->scale ==1<<j ? sfWhite : sfColor_fromRGB (150, 150, 150));
        }

        set_start_time (ds, ws);
      }
      return;
    }
  }
}


/* APPEND VERTEX TO VERTEX ARRAY. */

static void append_vertex (sfVertexArray *a, int x, int y, sfColor c)
{ sfVertex v;
  v.texCoords = zero_vector_f;
  v.color = c;
  v.position.x = x;
  v.position.y = y;
  sfVertexArray_append (a, v);
}


/* SET THE START TIMES FOR THE CURRENT RUNNING PERIOD. */

void set_start_time (struct dynamic_state *ds, struct window_state *ws)
{ if (ws->running)
  { ws->start_real_time = sfTime_asSeconds(sfClock_getElapsedTime(ws->clock));
    ws->start_sim_time = ds->sim_time;
  }
}
