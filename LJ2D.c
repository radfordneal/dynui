/* LJ2D.C - Program for simulating dynamics of Lennard-Jones particles in 2D.
   Copyright 2021 by Radford M. Neal
 */

#include "dynui.h"


/* STATE OF A LENNARD-JONES SYSTEM. */

struct LJ_state
{ double W, H;		/* Width and height of area */
  int N;		/* Number of particles */
  double T;		/* Equilibrium temperature */
  double initial_T;	/* Initial temperature */
  double alpha;		/* Temperature refresh constant, 1 = no refresh */
  unsigned seed;	/* State of random number generator */
  unsigned start_seed;	/* Seed for random generator at start of run */
  double *qx, *qy;	/* Particle positions */
  double *px, *py;	/* Particle momenta */
  double *gx, *gy;	/* Gradients of energy w.r.t. position */
};

#define I(ds) (*(struct LJ_state *)((ds)->i))


/* FILE FOR READING OR SAVING. */

static char *save_file;


/* ALLOCATE MEMORY FOR STATE. */

void alloc (struct LJ_state *I)
{
  if ((I->qx = calloc (I->N, sizeof (double))) == NULL
   || (I->qy = calloc (I->N, sizeof (double))) == NULL
   || (I->px = calloc (I->N, sizeof (double))) == NULL
   || (I->py = calloc (I->N, sizeof (double))) == NULL
   || (I->gx = calloc (I->N, sizeof (double))) == NULL
   || (I->gy = calloc (I->N, sizeof (double))) == NULL)
  { fprintf (stderr, "Can't allocate enough memory\n");
    exit (1);
  }
}


/* GENERATE UNIFORM RANDOM NUMBER IN (0,1). Uses and updates seed pointed to. */

double unif (unsigned *seed)
{ return (rand_r(seed)+1.0) / (RAND_MAX+2.0);
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
   "Usage: LJ2D [ save-file ] W H N / alpha T [ initial-T ] [ / seed ]\n");
  fprintf (stderr, 
   "   or: LJ2D save-file\n");
  exit(-1);
}


/* MAIN PROGRAM. */

int main (int argc, char **argv)
{
  struct LJ_state LJ;
  int from_file;
  char junk;
  int i, p;

  /* Parse program arguments, and perhaps read saved file. */

  if (argc < 2)
  { usage();
  }

  if (sscanf (argv[1], "%lf%c", &LJ.W, &junk) != 1 || LJ.W <= 0)  /* file arg */
  { save_file = argv[1];
    p = 2;
  }
  else
  { save_file = "LJsave";
    p = 1;
  }

  if (argc == p)  /* Read saved arguments and state from file */
  { FILE *f;
    from_file = 1;
    f = fopen (save_file, "rb");
    if (f == NULL)
    { fprintf (stderr, "Can't read file\n");
      exit(-2);
    }
    if (fread (&LJ.W, sizeof(LJ.W), 1, f) != 1
     || fread (&LJ.H, sizeof(LJ.H), 1, f) != 1
     || fread (&LJ.N, sizeof(LJ.N), 1, f) != 1
     || fread (&LJ.T, sizeof(LJ.T), 1, f) != 1
     || fread (&LJ.initial_T, sizeof(LJ.initial_T), 1, f) != 1
     || fread (&LJ.alpha, sizeof(LJ.alpha), 1, f) != 1
     || fread (&LJ.seed, sizeof(LJ.seed), 1, f) != 1
     || fread (&LJ.start_seed, sizeof(LJ.start_seed), 1, f) != 1)
    { fprintf (stderr, "Error reading file header\n");
      exit(-2);
    }
    alloc (&LJ);
    if (fread (LJ.qx, sizeof(*LJ.qx), LJ.N, f) != LJ.N
     || fread (LJ.qy, sizeof(*LJ.qy), LJ.N, f) != LJ.N
     || fread (LJ.px, sizeof(*LJ.px), LJ.N, f) != LJ.N
     || fread (LJ.py, sizeof(*LJ.py), LJ.N, f) != LJ.N)
    { fprintf (stderr, "Error reading file data\n");
      exit(-2);
    }
    fclose(f);
  }

  else  /* Parse arguments on command line */
  {
    from_file = 0;

    LJ.initial_T = -1;
    LJ.seed = 1;

    if (argc < p+6)
    { usage();
    }

    if (sscanf (argv[p+0], "%lf%c", &LJ.W, &junk) != 1 || LJ.W <= 0
     || sscanf (argv[p+1], "%lf%c", &LJ.H, &junk) != 1 || LJ.H <= 0
     || sscanf (argv[p+2], "%d%c",  &LJ.N, &junk) != 1 || LJ.N <= 0
     || strcmp (argv[p+3], "/") != 0
     || sscanf (argv[p+4], "%lf%c", &LJ.alpha, &junk) != 1 
         || LJ.alpha < 0 || LJ.alpha > 1
     || sscanf (argv[p+5], "%lf%c", &LJ.T, &junk) != 1 || LJ.T < 0)
    { usage();
    }

    if (argc > p+6)
    { if (strcmp (argv[p+6], "/") != 0)
      { if (sscanf (argv[p+6], "%lf%c", &LJ.initial_T, &junk) != 1
         || LJ.initial_T < 0)
        { usage();
        }
        p += 1;
      }
      if (argc > p+6)
      { if (strcmp (argv[p+6], "/") != 0 || argc != p+8
         || sscanf (argv[p+7], "%u%c", &LJ.seed, &junk) != 1)
        { usage();
        }
      }
    }

    LJ.start_seed = LJ.seed;
    if (LJ.initial_T < 0) 
    { LJ.initial_T = LJ.T;
    }
  }

  if (1)
  { fprintf (stderr,
     "file %s, W %.1f, H %.1f, N %d, alpha %.3f, T %.2f, initial_T %.2f\n",
      save_file, LJ.W, LJ.H, LJ.N, LJ.alpha, LJ.T, LJ.initial_T);
    fprintf (stderr, "seed %u, start_seed %u\n", LJ.seed, LJ.start_seed);
  }

  /* Allocate and initialize dynamic state. */

  struct dynamic_state ds0, *ds = &ds0;
  ds->i = (void *) &LJ;
  ds->sim_time = 0;

  if (!from_file)
  {
    alloc (&I(ds));

    for (i = 0; i < I(ds).N; i++)
    { I(ds).qx[i] = unif(&I(ds).seed) * I(ds).W;
      I(ds).qy[i] = unif(&I(ds).seed) * I(ds).H;
      I(ds).px[i] = norm(&I(ds).seed) * sqrt(I(ds).initial_T);
      I(ds).py[i] = norm(&I(ds).seed) * sqrt(I(ds).initial_T);
    }
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
{ FILE *f;
  f = fopen (save_file, "wb");
  if (f == NULL)
  { return 0;
  }
  if (fwrite (&I(ds).W, sizeof(I(ds).W), 1, f) != 1
   || fwrite (&I(ds).H, sizeof(I(ds).H), 1, f) != 1
   || fwrite (&I(ds).N, sizeof(I(ds).N), 1, f) != 1
   || fwrite (&I(ds).T, sizeof(I(ds).T), 1, f) != 1
   || fwrite (&I(ds).initial_T, sizeof(I(ds).initial_T), 1, f) != 1
   || fwrite (&I(ds).alpha, sizeof(I(ds).alpha), 1, f) != 1
   || fwrite (&I(ds).seed, sizeof(I(ds).seed), 1, f) != 1
   || fwrite (&I(ds).start_seed, sizeof(I(ds).start_seed), 1, f) != 1
   || fwrite (I(ds).qx, sizeof(*I(ds).qx), I(ds).N, f) != I(ds).N
   || fwrite (I(ds).qy, sizeof(*I(ds).qy), I(ds).N, f) != I(ds).N
   || fwrite (I(ds).px, sizeof(*I(ds).px), I(ds).N, f) != I(ds).N
   || fwrite (I(ds).py, sizeof(*I(ds).py), I(ds).N, f) != I(ds).N)
  { return 0;
  }
  fclose(f);
  return 1;
}


/* TERMINATE PROGRAM. */

void dynui_terminate()
{ exit(0);
}
