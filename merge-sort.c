/* MERGE-SORT.C - Template for merge sort procedure.
   Copyright (C) 2017, 2018, 2021 by Radford M. Neal
 */


/* Standalone test program.  Test merge_sort with

       cc -DTEST_MERGE_SORT merge-sort.c; a.out n1 n2 n3 ...
*/

#ifdef TEST_MERGE_SORT

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define merge_value int
#define merge_greater(a,b) ((a) > (b))

static void merge_sort (merge_value *dst, merge_value *src, int n);

int main (int argc, char **argv)
{
  int n = argc - 1;
  merge_value d[n], s[n];
  int i;
    
  for (i = 0; i < n; i++) s[i] = atoi(argv[i+1]);
  merge_sort (d, s, n);
  for (i = 0; i < n; i++) printf(" %d",d[i]);
  printf("\n");

  return 0;
}

#endif


/* Template for a merge sort procedure, which sorts 'n' items in
   'src', storing the sorted list in 'dst'.  The data in 'src' may be
   destroyed during the sort.  The sort is stable.

   To avoid name conflicts, 'merge_sort' should be defined to
   something appropriate before including merge-sort.c.  Also,
   <string.h> should be included before merge-sort.c.  A symbol
   'merge_value' must be defined, which is the type of an item value,
   or of an index to an item.  A macro 'merge_greater' must also be
   defined, which will be called as merge_greater(v1,v2), where 'v1'
   and 'v2' are of type 'merge_value'.

   These symbols can all be undefined, then redefined before including
   merge-sort.c again, to make another sort procedure. 

   The algorithm is recursive.  It can be viewed in terms of halves of
   src and dst, as follows:

           src = [ A | B ]     SORT->      dst = [ X | Y ]

   This is implemented as follows:

                B      SORT->    Y      (B destroyed)
                A      SORT->    B      (A destroyed)
               B,Y    MERGE->   dst
	       
   If src is of odd length, B is one bigger than A, and Y is one
   bigger than X.  The merge of B and Y will overwrite Y, but a write
   to an element of Y will happen only after that element is no longer
   needed.  Some writes in a merge are avoided, when the last elements
   merged are already in place. */

static void merge_sort (merge_value *dst, merge_value *src, int n)
{
  if (n == 2)
  { if (merge_greater (src[0], src[1]))
    { dst[0] = src[1];
      dst[1] = src[0];
    }
    else
    { dst[0] = src[0];
      dst[1] = src[1];
    }
  }
  else if (n > 2)
  {
    merge_sort (dst + n/2, src + n/2, n-n/2);

    if (n == 3)
    { if (!merge_greater (src[0], dst[1]))
      { dst[0] = src[0];
      }
      else
      { dst[0] = dst[1];
        if (!merge_greater (src[0], dst[2]))
        { dst[1] = src[0];
        }
        else
        { dst[1] = dst[2];
          dst[2] = src[0];
        }
      }
    }
    else
    {
      merge_sort (src + n/2, src, n/2);

      int n1 = n/2;
      int n2 = n-n1;

      merge_value *i1 = src + n1;
      merge_value *i2 = dst + n1;
      merge_value *j = dst;

      for (;;)
      { if (merge_greater (*i1, *i2))
        { *j++ = *i2++;
          n2 -= 1;
          if (n2 == 0)
          { if (n1 <= 2)
            { *j = *i1;
              if (n1 == 2)
              { *(j+1) = *(i1+1);
              }
            }
            else
            { memcpy (j, i1, n1*sizeof(merge_value));
            }
            break;
          }
        }
        else
        {
          *j++ = *i1++;
          n1 -= 1;

          /* Periodically check for a big block that can be copied.
             This helps when the data is partially sorted already. */

#         define blk 32  /* must be a power of two */
          if (n1 != 0 && (n1 & (blk-1)) == 0 
           && !merge_greater (*(i1+(blk-1)), *i2))
         {
            memcpy (j, i1, blk*sizeof(merge_value));
            j += blk;
            i1 += blk;
            n1 -= blk;
          }

          if (n1 == 0)
          { break;  /* no need to copy - data is already there */
          }
        }
      }
    }
  }
  else if (n > 0)
  { dst[0] = src[0];
  }
}
