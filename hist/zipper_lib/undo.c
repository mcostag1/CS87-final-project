/*

Undo facilities for transformation (motion) of objects.

Greg Turk, August 1993.

---------------------------------------------------------------

Copyright (c) 1994 The Board of Trustees of The Leland Stanford
Junior University.  All rights reserved.   
  
Permission to use, copy, modify and distribute this software and its   
documentation for any purpose is hereby granted without fee, provided   
that the above copyright notice and this permission notice appear in   
all copies of this software and that you do not sell the software.   
  
THE SOFTWARE IS PROVIDED "AS IS" AND WITHOUT WARRANTY OF ANY KIND,   
EXPRESS, IMPLIED OR OTHERWISE, INCLUDING WITHOUT LIMITATION, ANY   
WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.   

*/

#include <math.h>
#include <stdio.h>
#include <cyfile.h>
#include <zipper.h>
#include <matrix.h>

/* record of an object's position */

typedef struct Undo_Info {
  Scan *obj;			/* which object is being moved */
  Matrix rotmat;		/* rotation of object */
  float xtrans,ytrans,ztrans;	/* translation */
} Undo_Info;

static Undo_Info **undo_stack = NULL;
static int undo_ptr = 0;
static int undo_num = 0;
static int undo_max = 50;


/******************************************************************************
Save the current transformation of an object before it is moved.

Entry:
  scan - object to be moved
******************************************************************************/

save_for_undo(scan)
  Scan *scan;
{
  Undo_Info *item;

  /* don't add to stack if we've been redoing all the way up the stack */
  if (undo_ptr == undo_num - 1) {
    undo_ptr = undo_num;
    return;
  }

  /* allocate the stack if it hasn't been already */
  if (undo_stack == NULL)
    undo_stack = (Undo_Info **) myalloc (sizeof (Undo_Info *) * undo_max);

  /* make stack bigger if necessary */
  if (undo_num >= undo_max) {
    undo_max *= 2;
    undo_stack = (Undo_Info **)
      realloc (undo_stack, sizeof (Undo_Info *) * undo_max);
  }

  /* create new undo item and place it on the stack */

  item = (Undo_Info *) myalloc (sizeof (Undo_Info));
  item->obj = scan;
  mat_copy (item->rotmat, scan->rotmat);
  item->xtrans = scan->xtrans;
  item->ytrans = scan->ytrans;
  item->ztrans = scan->ztrans;
  undo_stack[undo_num++] = item;
  undo_ptr = undo_num;
}


/******************************************************************************
Undo the last transformation on the undo stack.
******************************************************************************/

perform_undo()
{
  Undo_Info *item;

  /* don't do anything if there is nothing on the stack */
  if (undo_num == 0)
    return;

  /* see if we have to save away current object's state */
  if (undo_ptr == undo_num) {
    save_for_undo (undo_stack[undo_num - 1]->obj);
    undo_ptr = undo_num - 1;
  }

  /* figure out where undo pointer is */
  undo_ptr--;
  if (undo_ptr < 0)
    undo_ptr = 0;

  /* undo the last transformation */

  item = undo_stack[undo_ptr];
  mat_copy (item->obj->rotmat, item->rotmat);
  item->obj->xtrans = item->xtrans;
  item->obj->ytrans = item->ytrans;
  item->obj->ztrans = item->ztrans;

  /* re-draw objects */
  draw_object();
}


/******************************************************************************
Redo the transformation.
******************************************************************************/

perform_redo()
{
  Undo_Info *item;

  /* don't do anything if there is nothing on the stack */
  if (undo_num == 0)
    return;

  /* don't redo if we're already at the top of the stack */
  if (undo_ptr == undo_num - 1) {
    return;
  }

  /* figure out where undo pointer should be */
  undo_ptr++;

  /* redo transformation */

  item = undo_stack[undo_ptr];
  mat_copy (item->obj->rotmat, item->rotmat);
  item->obj->xtrans = item->xtrans;
  item->obj->ytrans = item->ytrans;
  item->obj->ztrans = item->ztrans;

  /* re-draw objects */
  draw_object();
}

