/*

External declarations for viewing routines.

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

extern Matrix global_mat; /* global transformation matrix */
extern Matrix rotmat;     /* rotations */
extern Matrix transmat;   /* translations */
extern Matrix matident;   /* identity matrix */

/* for getting cursor position */
extern Device mdev[2];
extern short mval[2];

#define Cursor_Setup  getorigin (&xorg, &yorg);  getsize (&xsize, &ysize);
#define Cursor_Pos(x,y)  getdev (2, mdev, mval); \
                         x = mval[0] - xorg; \
                         y = mval[1] - yorg;

/* drawing window */

extern long draw_win;
extern long xorg;
extern long yorg;
extern long xsize;
extern long ysize;

