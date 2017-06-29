/*

Graphical cursors.

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

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <fcntl.h>
#include <gl/gl.h>
#include <gl/device.h>
#include <cyfile.h>
#include <strings.h>
#include <malloc.h>
#include <zipper.h>
#include <cursors.h>
#include <view.h>

extern int pick_flag;
extern int move_num;

#define MOVE_WORLD -1
#define MOVE_LIGHT -2

static Cursor cross_cursor = {
  0x0100, 0x0100, 0x0100, 0x0100,
  0x0100, 0x0100, 0x0100, 0xffff,
  0x0100, 0x0100, 0x0100, 0x0100,
  0x0100, 0x0100, 0x0100, 0x0100,
};

static Cursor wait_cursor = {
  0x1ff0, 0x1ff0, 0x0820, 0x0820,
  0x0820, 0x0c60, 0x06c0, 0x0100,
  0x0100, 0x06c0, 0x06c0, 0x0820,
  0x0820, 0x0820, 0x1ff0, 0x1ff0,
};

static Cursor light_cursor = {
  0x0810, 0x0c30, 0x0420, 0x03c0,
  0xc7e3, 0x6c36, 0x1818, 0x1818,
  0x1818, 0x1818, 0x6c36, 0xc7e3,
  0x03c0, 0x0420, 0x0c30, 0x0810,
};

static Cursor object_cursor = {
  0x0000, 0x07e0, 0x1ff8, 0x381c,
  0x300c, 0x6006, 0x6006, 0x6006,
  0x6006, 0x6006, 0x6006, 0x300c,
  0x381c, 0x1ff8, 0x07e0, 0x0000,
};

#define CURSOR_MAX 10
static int cursor_stack[CURSOR_MAX];
static int cursor_stack_count = 0;


/******************************************************************************
Define all my cursors.
******************************************************************************/

init_cursors()
{
  drawmode (CURSORDRAW);
  mapcolor (1, 255, 0, 0);

  curstype (C16X1);

  defcursor (CROSS_CURSOR, cross_cursor);
  defcursor (WAIT_CURSOR, wait_cursor);
  defcursor (OBJECT_CURSOR, object_cursor);
  defcursor (LIGHT_CURSOR, light_cursor);

  curorigin (CROSS_CURSOR, 7, 7);
  curorigin (LIGHT_CURSOR, 7, 7);

  drawmode (NORMALDRAW);

  cursor_stack[0] = DEFAULT_CURSOR;
  cursor_stack_count = 1;
}


/******************************************************************************
Set the appropriate cursor based on the moving mode.
******************************************************************************/

set_move_cursor()
{
  if (move_num == MOVE_WORLD)
    select_cursor (DEFAULT_CURSOR);
  else if (move_num == MOVE_LIGHT)
    select_cursor (LIGHT_CURSOR);
  else if (move_num >= 0)
    select_cursor (OBJECT_CURSOR);
}


/******************************************************************************
Set the appropriate cursor based on the picking flag.
******************************************************************************/

set_pick_cursor()
{
  if (pick_flag == 0)
    set_move_cursor();
  else
    select_cursor (CROSS_CURSOR);
}


/******************************************************************************
Select from a number of cursors.

Entry:
  num - number of the desired cursor
******************************************************************************/

select_cursor(num)
  int num;
{
  winset (draw_win);

  /* pick the appropriate cursor */

  cursor_stack[cursor_stack_count - 1] = num;

  switch (num) {
    case DEFAULT_CURSOR:
      setcursor (DEFAULT_CURSOR, 0, 0);
      break;
    case CROSS_CURSOR:
      setcursor (CROSS_CURSOR, 0, 0);
      break;
    case WAIT_CURSOR:
      setcursor (WAIT_CURSOR, 0, 0);
      break;
    case OBJECT_CURSOR:
      setcursor (OBJECT_CURSOR, 0, 0);
      break;
    case LIGHT_CURSOR:
      setcursor (LIGHT_CURSOR, 0, 0);
      break;
    default:
      fprintf (stderr, "select_cursor: bad cursor number = %d\n", num);
      break;
  }
}


/******************************************************************************
Switch to a new cursor, but pop the old one on the cursor stack.

Entry:
  num - number of the desired cursor
******************************************************************************/

push_cursor(num)
  int num;
{
  if (cursor_stack_count >= CURSOR_MAX) {
    fprintf (stderr, "push_cursor: cursor stack overflow\n");
    return;
  }

  cursor_stack_count++;
  select_cursor (num);
}


/******************************************************************************
Switch to an old cursor that was pushed onto the stack earlier.
******************************************************************************/

pop_cursor()
{
  if (cursor_stack_count <= 1) {
    fprintf (stderr, "push_cursor: cursor stack underflow\n");
    return;
  }

  cursor_stack_count--;
  select_cursor (cursor_stack[cursor_stack_count - 1]);
}

