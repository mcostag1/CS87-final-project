/*

Viewing routines.

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
#include <gl/gl.h>
#include <gl/device.h>
#include <fmclient.h>
#include "matrix.h"
#include "event.h"
#include "view.h"

/* viewing parameters */

float fov;              /* field of view */

float fx,fy,fz;         /* view from */
float ax,ay,az;         /* view at */
float ux,uy,uz;         /* view up */

float hither = 0.001;
float yon = 100;

Matrix global_mat; /* global transformation matrix */
Matrix rotmat;     /* rotations */
Matrix transmat;   /* translations */
Matrix matident;   /* identity matrix */

/* for getting cursor position */
Device mdev[2];
short mval[2];

/* drawing window */

long draw_win;
long xorg = 520;
long yorg = 240;

#if 1
long xsize = 750;
long ysize = 750;
#endif

#if 0
long xsize = 680;
long ysize = 486;
#endif

/* fonts */

fmfonthandle font1;
fmfonthandle font2;

int bits24;             /* are we on a 24-bit system? */

int is_orthographic = 0;


/******************************************************************************
Initialize the graphics state.
******************************************************************************/

init_graphics()
{
  int i,j;
  int nplanes;

  /* initialize for mouse picking */
  mdev[0] = MOUSEX;
  mdev[1] = MOUSEY;

  prefposition (xorg, xorg + xsize, yorg, yorg + ysize);

  foreground();
  draw_win = winopen ("Object Window");
  wintitle ("Object Window");
  winconstraints();

  RGBmode();
  doublebuffer();
  gconfig();

  mmode (MVIEWING);

  lookat (fx, fy, fz, ax, ay, az, 0);

  frontbuffer (TRUE);
  czclear (0x0, 0x0);

  zbuffer (FALSE);
  frontbuffer (FALSE);
  backbuffer (TRUE);

  /* determine if we're on a 24-bit system */

  nplanes = getplanes();
  if (nplanes > 12)
    bits24 = 1;
  else
    bits24 = 0;

  /* set the viewing parameters */
  reset_view();

  /* initialize font manager */

  fminit();
  font1 = fmfindfont ("HelveticaBold");
  font2 = fmscalefont (font1, 20.0);
  fmsetfont (font2);

  /* set up callback routines for windows */
  setup_draw_callbacks();

  /* initialize the cursors */
  init_cursors();
}


/******************************************************************************
Reset the viewing parameters.
******************************************************************************/

reset_view()
{
  int i,j;

  fov = 60;

  fx = 0;
  fy = 0;
  fz = -1;

  ax = 0;
  ay = 0;
  az = 0;

  ux = 0;
  uy = 1;
  uz = 0;

  /* initialize transformation matrices */

  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
      matident[i][j] = (i == j) ? 1.0 : 0.0;

  /* set transmat to identity */

  pushmatrix();
  loadmatrix (matident);
  getmatrix (transmat);
  popmatrix();
  transmat[3][1] = -0.15;
  transmat[3][2] = -0.7;

  /* set rotmat to 180 degree rotation around the y axis */

  pushmatrix();
  loadmatrix (matident);
  rot (180.0, 'y');
  getmatrix (rotmat);
  popmatrix();

  set_global_transformation();
}


/******************************************************************************
Set the global transformation.
******************************************************************************/

set_global_transformation()
{
  float scale, left, right, top, bottom;


  /* Why doesn't get_drawing_size() work? */
#if 0
  get_drawing_size(&xsize, &ysize);
#else
  getsize (&xsize, &ysize);
#endif

  viewport(0, xsize, 0, ysize);

  if (is_orthographic) {
      scale = 2*(1.2+transmat[3][2]);
      left = -0.15*scale;
      right = 0.15*scale;
      bottom = -0.15*scale*ysize/(float)xsize;
      top = 0.15*scale*ysize/(float)xsize;
      ortho(left, right, bottom, top, -10.0, 10.0);
  } else {
      perspective ((int) fov * 10, xsize / (float) ysize, hither, yon);
  }

  loadmatrix (matident);
  lookat (fx, fy, fz, ax, ay, az, 0);
  multmatrix (transmat);
  multmatrix (rotmat);
  getmatrix (global_mat);
}


set_orthographic(int val) {
    is_orthographic = val;
    draw_object();
}




/******************************************************************************
Set the global transformation with subpixel motion for anti-aliasing.
******************************************************************************/

set_subpixel_transformation(dx,dy)
  float dx,dy;
{
  subpixperspective((int) fov * 10, xsize / (float) ysize, hither, yon, dx, dy);
  loadmatrix (matident);
  lookat (fx, fy, fz, ax, ay, az, 0);
  multmatrix (transmat);
  multmatrix (rotmat);
  getmatrix (global_mat);
}


/******************************************************************************
Sub-pixel window position.
******************************************************************************/

subpixwindow(left,right,bottom,top,near,far,pixdx,pixdy)
  float left,right,bottom,top;
  float near,far;
  float pixdx,pixdy;
{
  short vleft,vright,vbottom,vtop;
  float xwsize,ywsize,dx,dy;
  int xpixels,ypixels;

  getviewport (&vleft, &vright, &vbottom, &vtop);
  xpixels = vright - vleft + 1;
  ypixels = vtop - vbottom + 1;
  xwsize = right - left;
  ywsize = top - bottom;
  dx = -pixdx * xwsize / xpixels;
  dy = -pixdy * ywsize / ypixels;
  window (left+dx, right+dx, bottom+dy, top+dy, near, far);
}


/******************************************************************************
Sub-pixel perspective.
******************************************************************************/

subpixperspective(fovy,aspect,near,far,pixdx,pixdy)
  Angle fovy;
  float aspect,near,far,pixdx,pixdy;
{
  float fov2,left,right,bottom,top;

  fov2 = ((fovy * M_PI) / 1800) / 2.0;
  top = near / (fcos (fov2) / fsin (fov2));
  bottom = -top;
  right = top * aspect;
  left = -right;
  subpixwindow (left, right, bottom, top, near, far, pixdx, pixdy);
}


/******************************************************************************
Return a matrix describing the transformation between world and screen
coordinates.
******************************************************************************/

world_to_screen(mat)
  Matrix mat;
{
  Matrix view_mat;

  pushmatrix();
  set_global_transformation();
  mmode (MPROJECTION);
  getmatrix (view_mat);
  mmode (MVIEWING);
  popmatrix();

  pushmatrix();
  loadmatrix (view_mat);
  multmatrix (global_mat);
  getmatrix (mat);
  popmatrix();
}


/******************************************************************************
Return the origin of the drawing window.
******************************************************************************/

get_drawing_origin(x,y)
  int *x,*y;
{
  int xx,yy;
  long last_win;

  last_win = winget();
  winset (draw_win);

  getorigin (&xorg, &yorg);
  *x = xorg;
  *y = yorg;

  winset (last_win);
}


/******************************************************************************
Return the size of the drawing window.
******************************************************************************/

get_drawing_size(x,y)
  int *x,*y;
{
  int xx,yy;
  long last_win;

  last_win = winget();
  winset (draw_win);

  getsize (&xsize, &ysize);
  *x = xsize;
  *y = ysize;

  winset (last_win);
}


/******************************************************************************
Apply a matrix to a vector, and perform the homogeneous divide.
******************************************************************************/

mat_apply_homog (mat,vec)
  Matrix mat;
  Vector vec;
{
  int i,j;
  Vector4 temp;

  /* down the columns */
  for (i = 0; i < 4; i++)
    temp[i] = vec[0] * mat[0][i] + vec[1] * mat[1][i] + vec[2] * mat[2][i]
	    + mat[3][i];

  vec[0] = temp[0] / temp[3];
  vec[1] = temp[1] / temp[3];
  vec[2] = temp[2] / temp[3];
}

 
/******************************************************************************
Print a matrix.
******************************************************************************/

print_mat (str,mat)
  char *str;
  Matrix mat;
{
  int i,j;

  printf ("\n");
  printf ("%s:\n", str);
  for (i = 0; i < 4; i++) {
    printf (" ");
    for (j = 0; j < 4; j++)
      printf ("%f ", mat[i][j]);
    printf ("\n");
  }
  printf ("\n");
}

