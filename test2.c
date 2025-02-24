#include <gmp.h>
#include "flint/flint.h"

void myfun(void)
{
   /* other variable declarations */
   nn_ptr a, b;
   TMP_INIT;

   /* arbitrary code */

   TMP_START; /* we are about to do some allocation */

   /* arbitrary code */

   a = TMP_ALLOC(32*sizeof(ulong));
   b = TMP_ALLOC(64*sizeof(ulong));

   /* arbitrary code */

   TMP_END; /* cleans up a and b */

   /* arbitrary code */
}