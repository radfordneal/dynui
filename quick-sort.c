/* QUICK-SORT.C - Template for quick-sort procedure.
   Copyright (C) 2021 by Radford M. Neal
 */


/* Standalone test program.  Test quick_sort with

       cc -DTEST_QUICK_SORT quick-sort.c; a.out n1 n2 n3 ...
*/

#ifdef TEST_QUICK_SORT

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define quick_value int
#define quick_greater(a,b) ((a) > (b))

static void quick_sort (quick_value *data, int n);

int main (int argc, char **argv)
{
  int n = argc - 1;
  quick_value s[n];
  int i;
    
  for (i = 0; i < n; i++) s[i] = atoi(argv[i+1]);
  quick_sort (s, n);
  for (i = 0; i < n; i++) printf(" %d",s[i]);
  printf("\n");

  return 0;
}

#endif


/* Template for a quick sort procedure, which sorts 'n' items in 'data',
   storing the sorted list in the same place.

   To avoid name conflicts, 'quick_sort' should be defined to
   something appropriate before including quick-sort.c.  Also,
   <string.h> should be included before quick-sort.c.  A symbol
   'quick_value' must be defined, which is the type of an item value,
   or of an index to an item.  A macro 'quick_greater' must also be
   defined, which will be called as quick_greater(v1,v2), where 'v1'
   and 'v2' are of type 'quick_value'.

   These symbols can all be undefined, then redefined before including
   quick-sort.c again, to make another sort procedure. */

static void quick_sort (quick_value *data, int n)
{
  for (;;)
  { 
    if (n < 2)
    { return;
    }

    if (n == 2)
    { if (quick_greater (data[0], data[1]))
      { double t = data[0];
        data[0] = data[1];
        data[1] = t;
      }
      return;
    }

    quick_value m = data[n/2];
    int i, j;

    i = 0; j = n-1;
    for (;;)
    { while (i <= j && ! quick_greater (data[i], m))
      { i += 1;
      }
      while (i <= j && ! quick_greater (m, data[j]))
      { j -= 1;
      }
      if (i > j)
      { break;
      }
      double t = data[i];
      data[i] = data[j];
      data[j] = t;
    }

    if (i == 0)
    { data[n/2] = data[0];
      data[0] = m;
      data += 1;
      n -= 1;
      continue;
    }

    if (i == n)
    { data[n/2] = data[n-1];
      data[n-1] = m;
      n -= 1;
      continue;
    }

    if (i > n/2)
    { quick_sort (data+i, n-i);
      n = i;
    }
    else
    { quick_sort (data, i);
      data += i;
      n -= i;
    }
  }
}
