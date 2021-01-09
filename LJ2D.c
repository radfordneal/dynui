/* LJ2D.C - Program for simulating dynamics of Lennard-Jones particles in 2D.
   Copyright 2021 by Radford M. Neal
 */

#include "dynui.h"


/* STATE OF A LENNARD-JONES SYSTEM. */

struct LJ_state
{ double W, H;		/* Width and height of area */
  int N;		/* Number of particles */
  double initial_T;	/* Initial temperature */
  double T;		/* Equilibrium temperature */
  double alpha;		/* Temperature refresh constant, 1 = no refresh */
  double *qx, *qy;	/* Particle positions */
  double *px, *py;	/* Particle momenta */
  double *gx, *gy;	/* Gradients of energy w.r.t. position */
  unsigned seed;	/* State of random number generator */
};

#define I(ds) (*(struct LJ_state *)((ds)->i))


/* GENERATE UNIFORM RANDOM NUMBER IN [0,1). Uses and updates seed pointed to. */

double unif (unsigned *seed)
{ return rand_r(seed) / (RAND_MAX+1.0);
}


/* GENERATE RANDOM N(0,1) VARIATE.  Uses Box-Muller method, but looks at
   only one of the values, to avoid needing to save more state. */

double norm (unsigned *seed)
{ double a = unif(seed);
  double b = unif(seed);
  return cos(2.0*M_PI*a) * sqrt(-2.0*log(b));
}


/* PRINT USAGE MESSAGE AND EXIT. */

void usage (void)
{ fprintf (stderr, 
   "Usage: LJ2D width height N-particles / alpha temp [ initial-temp ] [ / seed ]\n");
  exit(-1);
}


/* MAIN PROGRAM. */

int main (int argc, char **argv)
{
  struct LJ_state LJ;
  char junk;
  int i;

  /* Parse program arguments. */

  if (argc < 7)
  { usage();
  }

  LJ.initial_T = -1;
  LJ.seed = 1;

  if (sscanf (argv[1], "%lf%c", &LJ.W,     &junk) != 1 || LJ.W <= 0
   || sscanf (argv[2], "%lf%c", &LJ.H,     &junk) != 1 || LJ.H <= 0
   || sscanf (argv[3], "%d%c",  &LJ.N,     &junk) != 1 || LJ.N <= 0
   || strcmp (argv[4], "/") != 0
   || sscanf (argv[5], "%lf%c", &LJ.alpha, &junk) != 1 || LJ.alpha<0||LJ.alpha>1
   || sscanf (argv[6], "%lf%c", &LJ.T,     &junk) != 1 || LJ.T < 0)
  { usage();
  }

  if (argc > 7)
  { int p = 7;
    if (strcmp (argv[7], "/") != 0)
    { if (sscanf (argv[7], "%lf%c", &LJ.initial_T, &junk) != 1
       || LJ.initial_T < 0)
      { usage();
      }
      p = 8;
    }
    if (argc > p)
    { if (strcmp (argv[p], "/") != 0 || argc != p+2
       || sscanf (argv[p+1], "%u%c", &LJ.seed, &junk) != 1)
      { usage();
      }
    }
  }

  if (LJ.initial_T < 0) LJ.initial_T = LJ.T;

  if (1)
  { fprintf (stderr,
      "W %.3f, 1H %.3f, N %d, alpha %.3f, T %.3f, initial_T %.3f, seed %u\n",
      LJ.W, LJ.H, LJ.N, LJ.alpha, LJ.T, LJ.initial_T, LJ.seed);
  }

  /* Allocate and initialize dynamic state. */

  struct dynamic_state ds0, *ds = &ds0;
  ds->i = (void *) &LJ;
  ds->sim_time = 0;

  if ((I(ds).qx = calloc (I(ds).N, sizeof (double))) == NULL
   || (I(ds).qy = calloc (I(ds).N, sizeof (double))) == NULL
   || (I(ds).px = calloc (I(ds).N, sizeof (double))) == NULL
   || (I(ds).py = calloc (I(ds).N, sizeof (double))) == NULL
   || (I(ds).gx = calloc (I(ds).N, sizeof (double))) == NULL
   || (I(ds).gy = calloc (I(ds).N, sizeof (double))) == NULL)
  { fprintf (stderr, "Can't allocate enough memory\n");
    exit (1);
  }

  for (i = 0; i < I(ds).N; i++)
  { I(ds).qx[i] = unif(&I(ds).seed) * I(ds).W;
    I(ds).qy[i] = unif(&I(ds).seed) * I(ds).H;
    I(ds).px[i] = norm(&I(ds).seed) * sqrt(I(ds).initial_T);
    I(ds).py[i] = norm(&I(ds).seed) * sqrt(I(ds).initial_T);
  }

  /* Allocate window state and initialize application-specified fields. */

  struct window_state ws0, *ws = &ws0;

  ws->width = 600;
  ws->height = 500;
  ws->title = "LJ2D";
  ws->running = 0;
  ws->full_screen = 0;

  /* Create window and let user interact with it, with simulation running
     when requested. */

  dynui_window (ds, ws);

  dynui_terminate();
}


/* COMPUTE SQUARED DISTANCE OF THE NEAREST IMAGES OF A PAIR OF MOLECULES. */

static double squared_distance (struct dynamic_state *ds, int i, int j,
                                double *sdx, double *sdy)
{
  double dx, dy;

  dx =  I(ds).qx[i] - I(ds).qx[j];
  while (dx < 0) dx += I(ds).W;
  if (dx > I(ds).W/2) dx = I(ds).W - dx;
  *sdx = dx;

  dy =  I(ds).qy[i] - I(ds).qy[j];
  while (dy < 0) dy += I(ds).H;
  if (dy > I(ds).H/2) dy = I(ds).H - dy;
  *sdy = dy;

  return dx*dx + dy*dy;
}


/* COMPUTE CONTRIBUTION TO ENERGY FROM A PAIR, GIVEN SQUARED DISTANCE. */

static double pair_energy (double d2)
{ double t6 = 1/(d2*d2*d2);
  double t12 = t6*t6;
  return 4 * (t12 - t6);
}


/* COMPUTE DERIVATIVE OF PAIR ENERGY W.R.T. SQUARED DISTANCE. */

static double pair_energy_deriv (double d2)
{ double t6 = 1/(d2*d2*d2);
  double t12 = t6*t6;
  return (-24*t12 + 12*t6) / d2;
}


/* COMPUTE GRADIENT OF ENERGY W.R.T. MOLECULE POSITIONS. */

static void compute_gradient (struct dynamic_state *ds)
{ 
  double *gx = I(ds).gx;
  double *gy = I(ds).gy;
  int N = I(ds).N;
  double dx, dy, d2, g;
  int i, j;
  
  for (i = 0; i < N; i++) gx[i] = gy[i] = 0;

  for (i = 1; i < N; i++)
  { for (j = 0; j < i; j++)
    { d2 = squared_distance (ds, i, j, &dx, &dy);
      if (d2 < 1) continue;   /* Avoid extreme gradients */
      if (d2 > 20) continue;  /* Ignore small forces */
      g = 2 * pair_energy_deriv(d2);
      gx[i] += dx*g;
      gy[i] += dy*g;
      gx[j] -= dx*g;
      gy[j] -= dy*g;
    }
  }
}


/* ADVANCE SIMULATION BY A SMALL TIME AMOUNT. */

static double delta_t = 0.1;
static int steps = 30;

void dynui_advance (struct dynamic_state *ds)
{
  double eta = delta_t / steps;
  double half_eta = eta/2;
  double alpha = I(ds).alpha;
  double *qx = I(ds).qx;
  double *qy = I(ds).qy;
  double *px = I(ds).px;
  double *py = I(ds).py;
  double *gx = I(ds).gx;
  double *gy = I(ds).gy;
  double W = I(ds).W;
  double H = I(ds).H;
  int N = I(ds).N;
  int i, t;

  compute_gradient(ds);

  for (i = 0; i < N; i++) 
  { px[i] -= half_eta * gx[i];
    py[i] -= half_eta * gy[i];
  }

  t = 0;
  for (;;)
  { 
    for (i = 0; i < N; i++) 
    { qx[i] += eta * px[i];
      while (qx[i] < 0) qx[i] += W;
      while (qx[i] > W) qx[i] -= W;
      qy[i] += eta * py[i];
      while (qy[i] < 0) qy[i] += H;
      while (qy[i] > H) qy[i] -= H;
    }

    compute_gradient(ds);

    if (++t == steps) 
    { break;
    }

    for (i = 0; i < N; i++) 
    { px[i] -= eta * gx[i];
      py[i] -= eta * gy[i];
    }
  }

  for (i = 0; i < N; i++) 
  { px[i] -= half_eta * gx[i];
    py[i] -= half_eta * gy[i];
  }

  if (alpha < 1)
  { double s = sqrt (I(ds).T * (1-alpha*alpha));
    unsigned *seed = &I(ds).seed;
    for (i = 0; i < N; i++) 
    { px[i] = alpha * px[i] + s * norm(seed);
      py[i] = alpha * py[i] + s * norm(seed);
    }
  }

  ds->sim_time += delta_t;
}


/* DRAW VIEW OF SIMULATION. */

static int dot_radius = 2;

void dynui_view (struct dynamic_state *ds, struct window_state *ws)
{ 
  double scaled_radius = ws->scale * dot_radius;

  double *qx = I(ds).qx;
  double *qy = I(ds).qy;
  double W = I(ds).W;
  double H = I(ds).H;
  int N = I(ds).N;
  int i;

  sfVector2f dot_origin = { scaled_radius, scaled_radius };

  sfCircleShape *dot = sfCircleShape_create();
  sfCircleShape_setOrigin (dot, dot_origin);
  sfCircleShape_setRadius (dot, scaled_radius);
  sfCircleShape_setFillColor (dot, sfWhite);

  double ox = ws->width/2, oy = ws->height/2;

  for (i = 0; i < N; i++)
  {
    sfVector2f pos = { ws->scale * (ws->offset.x + qx[i] - W/2) + ox,
                       ws->scale * (ws->offset.y + qy[i] - H/2) + oy };

    while (pos.x >= -scaled_radius) pos.x -= ws->scale * W;
    while (pos.x < -scaled_radius) pos.x += ws->scale * W;
    while (pos.y >= -scaled_radius) pos.y -= ws->scale * H;
    while (pos.y < -scaled_radius) pos.y += ws->scale * H;

    sfVector2f p1, p2;
    for (p1 = pos; p1.x <= ws->width+scaled_radius; p1.x += ws->scale * W)
    { for (p2 = p1; p2.y <= ws->height+scaled_radius; p2.y += ws->scale * H)
      { sfCircleShape_setPosition (dot, p2);
        sfRenderWindow_drawCircleShape (ws->window, dot, NULL);
      }
    }
  }

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
