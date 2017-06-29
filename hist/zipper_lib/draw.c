/*

Draw various objects.

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
#include <gl/image.h>
#include <tcl.h>
#include <strings.h>
#include <malloc.h>
#ifdef VOID
#undef VOID
#endif
#include <cyfile.h>
#include <zipper.h>
#include <event.h>
#include <view.h>
#include <trackball.h>
#include "zipglobal.h"

extern Tcl_Interp *interp;
extern Vertex *find_nearest();

void left_proc();
void middle_proc();
void right_proc();
void quit_proc();
void refresh_drawing();



extern int global_dont_draw;
extern int is_orthographic;

static int poly_draw_flag = 1;	/* draw polygons? */
static int near_flag = 0;	/* show lines connecting nearby mesh points? */
static int nearest_flag = 0;	/* similar to above */
static int smooth_flag = 1;	/* use smooth shading? */
static int normal_flag = 0;	/* display surface normals? */
static int colors_flag = 1;	/* different colors for different meshes? */
static int true_color_flag = 0;	/* display the color at vertices? */
static int refine_flag = 1;	/* use progressive refinement? */
static int value_flag = 0;	/* show vertex value as color? */
static int intensity_flag = 0;	/* show vertex intensities? */
static int draw_axes = 0;	/* draw a set of axes with the object? */
static int singlebuf = 0;
int edge_flag = 0;		/* draw polygons on the edge of the mesh? */
int bound_flag = 0;		/* draw edges on boundary of mesh? */
int pick_flag = 0;		/* pick a graphical primitive? */
int back_cull = 1;		/* cull backfacing polygons? */
int back_red = 0;               /* background color */
int back_grn = 0;
int back_blu = 0;
int copied_to_back = 0;	        /* has the image been copied to backbuffer? */

float mouse_speed = 1;		/* speed of object motion based on mouse */


/* This should be a variable; a function of hardware performance */
#define MAX_TRIANGLES_FOR_INTERACTION 5000


#define MOVE_WORLD -1
#define MOVE_LIGHT -2
int move_num = MOVE_WORLD;		/* move an object? */

/* for fast drawing while moving the object */

static int pushed_mesh_level;	/* for saving the mesh level */
static int partial_flag = 0;	/* are we only drawing part of an object? */
static int drew_partial = 0;	/* did we draw partial object? */

/* lights and materials */

static float poly_material[] = {
  AMBIENT, 0.1, 0.1, 0.1,
  DIFFUSE, 0.8, 0.8, 0.8,
  SPECULAR, 0.0, 0.0, 0.0,
  SHININESS, 0,
  LMNULL
};


static float my_light1[] = {
  LCOLOR, 1, 1, 1,
  POSITION, 0, 0, 1, 0,
  LMNULL
};

static float my_light2[] = {
  LCOLOR, 1, 1, 1,
  POSITION, 0, 0, -1, 0,
  LMNULL
};

static float light_model[] = {
  AMBIENT, 0.1, 0.1, 0.1,
  LOCALVIEWER, 0,
  LMNULL
};

/* "extra" lines to be displayed each frame */
static Vector *line1 = NULL;
static Vector *line2 = NULL;
static unsigned long *line_colors = NULL;
static int line_num;
static int line_max;


/******************************************************************************
Set up the callback routines for the drawing window.
******************************************************************************/

setup_draw_callbacks()
{
  /* set up callback routines for windows */

  add_event (draw_win, LEFTMOUSE, DOWN, left_proc, 0);
  add_event (draw_win, MIDDLEMOUSE, DOWN, middle_proc, 0);
  add_event (draw_win, RIGHTMOUSE, DOWN, right_proc, 0);
  add_event (ANY, REDRAW, draw_win, refresh_drawing, 0);
  add_event (ANY, ESCKEY, draw_win, quit_proc, 0);
  qdevice (REDRAW);
  qdevice (LEFTMOUSE);
  qdevice (MIDDLEMOUSE);
  qdevice (RIGHTMOUSE);
  qdevice (ESCKEY);
}


/******************************************************************************
Draw the current object
******************************************************************************/

draw_object()
{
  long result;

  if (global_dont_draw)
      return;

  /* unset pick */
  nothing_picked();

  winset (draw_win);

  /* set the global transformation */
  pushmatrix();
  set_global_transformation();

  /* clear screen */
  drawmode (OVERDRAW);
  color (0);
  clear();
  drawmode (NORMALDRAW);
  cpack (back_red | back_grn << 8 | back_blu << 16);
  clear();

  /* draw the object */

  draw_triangles();

  if (!singlebuf)
      swapbuffers();

  popmatrix();
  copied_to_back = 0;
}


/******************************************************************************
Draw an anti-aliased version of the object.
******************************************************************************/

draw_anti_aliased()
{
  int i,j;
  float dx,dy;
  float samples = 4;

  /* unset pick */
  nothing_picked();

  winset (draw_win);

  drawmode (OVERDRAW);
  color (0);
  clear();
  drawmode (NORMALDRAW);

  subpixel (1);
  acsize (16);
  acbuf (AC_CLEAR, 0.0);

  for (i = 0; i < samples; i++)
    for (j = 0; j < samples; j++) {

      dx = (i + drand48() * 0.25) / samples;
      dy = (j + drand48() * 0.25) / samples;

      pushmatrix();
      set_subpixel_transformation (dx, dy);
      cpack (back_red | back_grn << 8 | back_blu << 16);
      clear();
      draw_triangles();
      popmatrix();
      acbuf (AC_ACCUMULATE, 1.0);
    }

  acbuf (AC_RETURN, 1.0 / (samples * samples));
  if (!singlebuf)
      swapbuffers();

  subpixel (0);
}


/******************************************************************************
Spin around an object.

Entry:
  angle   - how far to spin (in degrees)
  nsteps  - how many steps to take during spin
  prefix  - prefix of filenames to write images to, or NULL if no file writing
  samples - number of supersamples across and down (or 1 for no anti-aliasing)
******************************************************************************/

spin_object(angle,nsteps,prefix,samples)
  float angle;
  int nsteps;
  char *prefix;
  int samples;
{
  int i,j,k;
  float theta;
  float dx,dy;

  /* unset pick */
  nothing_picked();

  winset (draw_win);

  drawmode (OVERDRAW);
  color (0);
  clear();
  drawmode (NORMALDRAW);

  for (k = 0; k <= nsteps; k++) {

    theta = angle * k / (float) nsteps;

    if (samples == 1) {	/* no anti-aliasing */
      pushmatrix();
      set_global_transformation();
      rot (theta, 'y');
      cpack (back_red | back_grn << 8 | back_blu << 16);
      clear();
      draw_triangles();
      if (!singlebuf)
	  swapbuffers();
      popmatrix();
    }
    else {	/* anti-alias */

      /* get ready to use the accumulation buffer */
      subpixel (1);
      acsize (16);
      acbuf (AC_CLEAR, 0.0);

      /* jittered super-sampling */
      for (i = 0; i < samples; i++)
	for (j = 0; j < samples; j++) {

	  dx = (i + drand48() * 0.25) / samples;
	  dy = (j + drand48() * 0.25) / samples;

	  pushmatrix();
	  set_subpixel_transformation (dx, dy);
	  rot (theta, 'y');
          cpack (back_red | back_grn << 8 | back_blu << 16);
	  clear();
	  draw_triangles();
	  popmatrix();
	  acbuf (AC_ACCUMULATE, 1.0);
	}

      acbuf (AC_RETURN, 1.0 / (samples * samples));
      if (!singlebuf)
	  swapbuffers();

      subpixel (0);

    }

    if (prefix != NULL) {
      write_frame (prefix, k);
    }
  }
}


/******************************************************************************
Evaluate various flags specified in Tcl.
******************************************************************************/

eval_flags()
{
  char *str;

  /* find out various drawing modes */

  str = Tcl_GetVar (interp, "polygon_state", 0);
  if (strcmp (str, "0") == 0)
    poly_draw_flag = 0;
  else
    poly_draw_flag = 1;
}


/******************************************************************************
Set the polygon material properties.

Entry:
  ambient  - ambient coefficient
  diffuse  - diffuse coefficient
  specular - specular coefficient
  spec_pow - specular power (or 0 for no specularity)
******************************************************************************/

void
set_polygon_material(ambient,diffuse,specular,spec_pow)
  float ambient;
  float diffuse;
  float specular;
  float spec_pow;
{
  poly_material[1] = ambient;
  poly_material[2] = ambient;
  poly_material[3] = ambient;

  poly_material[5] = diffuse;
  poly_material[6] = diffuse;
  poly_material[7] = diffuse;

  poly_material[9] = specular;
  poly_material[10] = specular;
  poly_material[11] = specular;

  poly_material[13] = spec_pow;
}


new_material(ambient)
  float ambient;
{
  poly_material[1] = ambient;
  poly_material[2] = ambient;
  poly_material[3] = ambient;

}


/******************************************************************************
Draw the triangles from the scans.
******************************************************************************/

draw_triangles()
{
  int i,j;
  Scan *sc;
  Mesh *mesh;
  Matrix mat;

  /* identity matrix */
  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
      mat[i][j] = (i == j);

  /* draw white object */
  cpack (0xffffff);

  /* maybe set up for polygon drawing */

  if (poly_draw_flag) {

    lsetdepth (0, 0x7fffff);
    zbuffer (TRUE);
    zclear();

    lmdef (DEFMATERIAL, 1, 0, poly_material);

    lmdef (DEFLIGHT, 1, 0, my_light1);
    lmdef (DEFLIGHT, 2, 0, my_light2);
    lmdef (DEFLMODEL, 1, 0, light_model);

    mmode(MVIEWING);

    lmbind (MATERIAL, 1);
    lmbind (LMODEL, 1);

    pushmatrix();
    loadmatrix (mat);
    lmbind (LIGHT0, 1);
    /*
    lmbind (LIGHT1, 2);
    */
    popmatrix();
  }

  /* draw each object */

  for (i = 0; i < nscans; i++) {

    /* shorthands */
    sc = scans[i];
    mesh = sc->meshes[mesh_level];

    /* draw the mesh if the flag says to */
    if (sc->draw_flag) {
      draw_mesh (mesh, sc, i);
#if 0
      if (sc->edge_mesh)
	draw_mesh (sc->edge_mesh, sc, i);
#endif
    }
  }

  /* maybe draw lines showing nearby points between meshes */

  if (near_flag)
    draw_near_lines();

  if (nearest_flag)
    draw_nearest_lines();

  if (bound_flag)
    draw_edge_boundaries();

  if (normal_flag)
    draw_normals();

  if (draw_axes) {
    Vector vec;
    float d = 0.2;

    set_line_color (0.2, 0.99, 0.0);

    bgnline();
    vset (vec, 0.0, 0.0, -d);
    v3f (vec);
    vset (vec, 0.0, 0.0,  d);
    v3f (vec);
    endline();

    bgnline();
    vset (vec, 0.0, -d, 0.0);
    v3f (vec);
    vset (vec, 0.0,  d, 0.0);
    v3f (vec);
    endline();

    bgnline();
    vset (vec, -d, 0.0, 0.0);
    v3f (vec);
    vset (vec,  d, 0.0, 0.0);
    v3f (vec);
    endline();
  }

  /* draw the "extra" lines in a frame */
  draw_extra_lines();

  if (poly_draw_flag) {
    zbuffer (FALSE);
    lmbind (LMODEL, 0);
  }
}


choose_mesh_level()
{
  int new_level, best_level;
  best_level = best_mesh_level();
  new_level = best_level > mesh_level ? best_level : mesh_level;
  
  save_mesh_level (new_level);
}


int
best_mesh_level()
{
  int i, j;
  int num_tris;
  Scan *sc;
  Mesh *mesh;

  for (j = 0; j <= 3; j++) {
      num_tris = 0;
      for (i = 0; i < nscans; i++) {

	  sc = scans[i];
	  mesh = sc->meshes[j];
	  
	  if (sc->draw_flag && mesh != NULL) {
	      num_tris += mesh->ntris;
	  } 
	  else if (mesh == NULL) {
	      num_tris += 2*MAX_TRIANGLES_FOR_INTERACTION;
	      break;
	  }
      }
      if (num_tris < MAX_TRIANGLES_FOR_INTERACTION)
	  return j;
  }

  /* Nothing satisfied the criterion */
  return 3;
}



/******************************************************************************
Draw a mesh.

Entry:
  mesh  - mesh to draw
  sc    - scan of mesh
  index - base color index
******************************************************************************/

draw_mesh(mesh,sc,index)
  Mesh *mesh;
  Scan *sc;
  int index;
{
  int j,k;
  Vertex **v;
  Matrix mat,timat;
  int some_flag = 0;
  int some_max = 3000;
  Triangle *tri;
  Vector norm;

  /* don't draw if the appropriate flag is set */
  if (global_dont_draw)
    return;

  /* if we're in refine-mode and this is a polyfile object, we will */
  /* only draw some of the polygons */

  if (partial_flag && sc->file_type == POLYFILE && mesh->ntris > some_max) {
    some_flag = 1;
    srand48 (1234);
  }

  /* different colors? */
  if (colors_flag)
    different_colors (index);

  /* transform object appropriately */

  pushmatrix();

  translate (sc->xtrans, sc->ytrans, sc->ztrans);
  multmatrix (sc->rotmat);

  getmatrix (mat);
  mat_copy (timat, mat);
  mat_invert (timat);
  mat_transpose (timat);

  /* draw the object */
  for (j = 0; j < mesh->ntris; j++) {

    tri = mesh->tris[j];
    v = tri->verts;

    /* maybe draw a few random triangles and then bail */
    if (some_flag) {
      /* maybe break out of loop */
      if (j > some_max)
	break;
      /* pick a random triangle to draw */
      k = trunc (drand48() * mesh->ntris);
      tri = mesh->tris[k];
      v = tri->verts;
      drew_partial = 1;
    }

    /* maybe don't draw polygons on edge of mesh */
    if (edge_flag && (v[0]->on_edge || v[1]->on_edge || v[2]->on_edge))
      continue;

    /* see if we're culling lines and a polygon is backfacing */
    if (back_cull && !poly_draw_flag && backface_tri(tri, mat, timat))
      continue;

    /* draw lines or filled polygons */

    if (colors_flag && tri->mark)
      different_colors (index + tri->mark);

    if (poly_draw_flag) {
      bgnpolygon();
      if (smooth_flag) {
	if (value_flag) {
	  set_interp_color (v[0]->confidence);
	  n3f (v[0]->normal);
	  v3f (v[0]->coord);
	  set_interp_color (v[1]->confidence);
	  n3f (v[1]->normal);
	  v3f (v[1]->coord);
	  set_interp_color (v[2]->confidence);
	  n3f (v[2]->normal);
	  v3f (v[2]->coord);
	}
	else if (intensity_flag) {
	  set_vertex_intensity (v[0]->intensity);
	  n3f (v[0]->normal);
	  v3f (v[0]->coord);
	  set_vertex_intensity (v[1]->intensity);
	  n3f (v[1]->normal);
	  v3f (v[1]->coord);
	  set_vertex_intensity (v[2]->intensity);
	  n3f (v[2]->normal);
	  v3f (v[2]->coord);
	}
	else if (true_color_flag) {
	  pcol (v[0]->red / 255.0, v[0]->grn / 255.0, v[0]->blu / 255.0);
	  n3f (v[0]->normal);
	  v3f (v[0]->coord);
	  pcol (v[1]->red / 255.0, v[1]->grn / 255.0, v[1]->blu / 255.0);
	  n3f (v[1]->normal);
	  v3f (v[1]->coord);
	  pcol (v[2]->red / 255.0, v[2]->grn / 255.0, v[2]->blu / 255.0);
	  n3f (v[2]->normal);
	  v3f (v[2]->coord);
	}
	else {
	  n3f (v[0]->normal);
	  v3f (v[0]->coord);
	  n3f (v[1]->normal);
	  v3f (v[1]->coord);
	  n3f (v[2]->normal);
	  v3f (v[2]->coord);
	}
      }
      else {
	if (value_flag) {
	  float t;
	  t = 0.333 * (v[0]->confidence + v[1]->confidence + v[2]->confidence);
	  set_interp_color (t);
          /*
	  n3f (mesh->tris[j]->normal);
          */
          norm[X] = -mesh->tris[j]->aa;
          norm[Y] = -mesh->tris[j]->bb;
          norm[Z] = -mesh->tris[j]->cc;
	  n3f (norm);
	  v3f (v[0]->coord);
	  v3f (v[1]->coord);
	  v3f (v[2]->coord);
	}
	else if (intensity_flag) {
	  float t;
	  t = 0.333 * (v[0]->intensity + v[1]->intensity + v[2]->intensity);
	  set_vertex_intensity (t);
          /*
	  n3f (mesh->tris[j]->normal);
          */
          norm[X] = -mesh->tris[j]->aa;
          norm[Y] = -mesh->tris[j]->bb;
          norm[Z] = -mesh->tris[j]->cc;
	  n3f (norm);
	  v3f (v[0]->coord);
	  v3f (v[1]->coord);
	  v3f (v[2]->coord);
	}
	else {
          /*
	  n3f (mesh->tris[j]->normal);
          */
          norm[X] = -mesh->tris[j]->aa;
          norm[Y] = -mesh->tris[j]->bb;
          norm[Z] = -mesh->tris[j]->cc;
	  n3f (norm);
	  v3f (v[0]->coord);
	  v3f (v[1]->coord);
	  v3f (v[2]->coord);
	}
      }
      endpolygon();
    }
    else {
      /* line drawing */
      bgnline();
      v3f (v[0]->coord);
      v3f (v[1]->coord);
      v3f (v[2]->coord);
      v3f (v[0]->coord);
      endline();
    }

    if (colors_flag && tri->mark)
      different_colors (index);
  }

  popmatrix();
}


/******************************************************************************
Set the color based on a scalar.  To be used for values in [0,1].
******************************************************************************/

set_interp_color(t)
  float t;
{
  static float r1 = 0.0, g1 = 0.0, b1 = 0.9;
  static float r2 = 0.9, g2 = 0.0, b2 = 0.0;

  if (t < -0.5)
    pcol (0.4, 0.4, 0.4);
  else {
    pcol (r1 + t * (r2 - r1), g1 + t * (g2 - g1), b1 + t * (b2 - b1));
/*    pcol (t,t,t);*/
  }
}


/******************************************************************************
Set the intensity of a vertex.
******************************************************************************/

set_vertex_intensity(t)
  float t;
{
  static float r1 = 0.0, g1 = 0.0, b1 = 0.0;
  static float r2 = 1.0, g2 = 1.0, b2 = 1.0;

  if (t < 0.0)
    pcol (r1, g1, b1);
  else if (t > 1.0)
    pcol (r2, g2, b2);
  else
    pcol (r1 + t * (r2 - r1), g1 + t * (g2 - g1), b1 + t * (b2 - b1));
}


/******************************************************************************
Say if a triangle is front or back facing.

Entry:
  tri   - triangle to determine facing direction
  mat   - current transformation matrix
  timat - inverse transpose of current transformation

Exit:
  returns 1 if triangle is backfacing, 0 if it is frontfacing
******************************************************************************/

int backface_tri(tri,mat,timat)
  Triangle *tri;
  Matrix mat,timat;
{
  Vector norm;
  Vector pnt;
  Vertex **v = tri->verts;
  float dot;

  /* transform the normal of the triangle */
  norm[X] = -tri->aa;
  norm[Y] = -tri->bb;
  norm[Z] = -tri->cc;
  mat_apply (timat, norm);


  if (is_orthographic) {
      if (norm[Z] < 0)
	  return (1);
      else
	  return (0);
  };


  /* transform the centriod of the triangle */
  pnt[X] = (v[0]->coord[X] + v[1]->coord[X] + v[2]->coord[X]) * 0.333333;
  pnt[Y] = (v[0]->coord[Y] + v[1]->coord[Y] + v[2]->coord[Y]) * 0.333333;
  pnt[Z] = (v[0]->coord[Z] + v[1]->coord[Z] + v[2]->coord[Z]) * 0.333333;
  mat_apply (mat, pnt);

  dot = vdot (norm, pnt);

  if (dot > 0)
    return (1);
  else
    return (0);
}


/******************************************************************************
Draw the boundary edges of all the meshes.
******************************************************************************/

draw_edge_boundaries()
{
  int i,j;
  Scan *sc;
  Mesh *mesh;
  long old_width;

  /* draw the boundaries as fat green lines */
  set_line_color (0.4, 0.9, 0.3);
  old_width = getlwidth();
  linewidth (3);

  /* draw the boundary edges */

  for (i = 0; i < nscans; i++) {

    /* shorthands */
    sc = scans[i];
    mesh = sc->meshes[mesh_level];

    /* should we draw this scan? */
    if (sc->draw_flag == 0)
      continue;

    pushmatrix();
    translate (sc->xtrans, sc->ytrans, sc->ztrans);
    multmatrix (sc->rotmat);

    /* create info about the boundary edges, if necessary */
    if (!mesh->edges_valid)
      create_edge_list (mesh);

    for (j = 0; j < mesh->nedges; j++) {
      Edge *e;
      e = mesh->edges[j];
      bgnline();
      v3f (e->v1->coord);
      v3f (e->v2->coord);
      endline();
    }

    popmatrix();
  }

  /* restore old line width */
  linewidth (old_width);
}


/******************************************************************************
Draw the surface normals of the triangles.
******************************************************************************/

draw_normals()
{
  int i,j;
  Scan *sc;
  Mesh *mesh;
  Triangle *tri;
  Vector v0,v1;

  /* draw the normals of each triangle of each displayed mesh */

  for (i = 0; i < nscans; i++) {

    /* shorthands */
    sc = scans[i];
    mesh = sc->meshes[mesh_level];

    /* should we draw this scan? */
    if (sc->draw_flag == 0)
      continue;

    pushmatrix();
    translate (sc->xtrans, sc->ytrans, sc->ztrans);
    multmatrix (sc->rotmat);

    set_line_color (0.2, 0.99, 0.0);

    for (j = 0; j < mesh->ntris; j++) {

      tri = mesh->tris[j];

      /* find center of triangle */
      vadd (tri->verts[0], tri->verts[1], v0);
      vadd (tri->verts[2], v0, v0);
      vscale (v0, 0.3333333);

      /* add surface normal to center */
      v1[X] = -tri->aa;
      v1[Y] = -tri->bb;
      v1[Z] = -tri->cc;
      vscale (v1, 0.005);
      vadd (v0, v1, v1);

      /* draw surface normal */
      bgnline();
      v3f (v0);
      v3f (v1);
      endline();
    }

    popmatrix();
  }
}


/******************************************************************************
Set color for a line.
******************************************************************************/

set_line_color(r,g,b)
  float r,g,b;
{
  int red,grn,blu;

  red = (int) (255 * r);
  grn = (int) (255 * g);
  blu = (int) (255 * b);
  cpack (red | (grn << 8) | (blu << 16));
}


/******************************************************************************
Set the color.
******************************************************************************/

pcol(r,g,b)
  float r,g,b;
{
  unsigned int col;
  int red,grn,blu;
  static float my_material[] = {
    AMBIENT, 0.1, 0.1, 0.1,
    DIFFUSE, 0.8, 0.8, 0.8,
    SPECULAR, 0.0, 0.0, 0.0,
    SHININESS, 0,
  };

  if (poly_draw_flag) {
    my_material[1] = poly_material[1];
    my_material[2] = poly_material[2];
    my_material[3] = poly_material[3];
    my_material[5] = r*poly_material[5];
    my_material[6] = g*poly_material[6];
    my_material[7] = b*poly_material[7];
    my_material[9] = poly_material[9];
    my_material[10] = poly_material[10];
    my_material[11] = poly_material[11];
    my_material[13] = poly_material[13];
    lmdef (DEFMATERIAL, 1, 0, my_material);
    lmbind (MATERIAL, 1);
  }
  else {
    red = (int) (255 * r);
    grn = (int) (255 * g);
    blu = (int) (255 * b);
    cpack (red | (grn << 8) | (blu << 16));
  }
}


/******************************************************************************
Use different colors, based on an index.
******************************************************************************/

different_colors(index)
  int index;
{
  switch (index) {
    case 0:
      pcol (0.9*1.1, 0.9*1.1, 0.9*1.1);
      break;
    case 1:
      pcol (0.9*1.1, 0.15*1.1, 0.35*1.1);
      break;
    case 2:
      pcol (0.0*1.1, 0.9*1.1, 0.0*1.1);
      break;
    case 3:
      pcol (0.0*1.1, 0.9*1.1, 0.9*1.1);
      break;
    case 4:
      pcol (0.9*1.1, 0.0*1.1, 0.9*1.1);
      break;
    case 5:
      pcol (0.9*1.1, 0.9*1.1, 0.0*1.1);
      break;
    case 6:
      pcol (0.0*1.1, 0.0*1.1, 0.9*1.1);
      break;
    case 7:
      pcol (0.9*1.1, 0.5*1.1, 0.2*1.1);
      break;
    case 8:
      pcol (0.5*1.1, 0.9*1.1, 0.2*1.1);
      break;
    case 9:
      pcol (0.2*1.1, 0.5*1.1, 0.9*1.1);
      break;
    case 10:
      pcol (0.5*1.1, 0.2*1.1, 0.9*1.1);
      break;
    case 11:
      pcol (0.9*1.1, 0.2*1.1, 0.5*1.1);
      break;
    case 12:
      pcol (0.2*1.1, 0.9*1.1, 0.5*1.1);
      break;
    case 13:
      pcol (0.9*1.1, 0.7*1.1, 0.7*1.1);
      break;
    case 14:
      pcol (0.7*1.1, 0.9*1.1, 0.7*1.1);
      break;
    case 15:
      pcol (0.7*1.1, 0.7*1.1, 0.9*1.1);
      break;
    case 16:
      pcol (0.7*1.1, 0.9*1.1, 0.9*1.1);
      break;
    case 17:
      pcol (0.9*1.1, 0.7*1.1, 0.9*1.1);
      break;
    case 18:
      pcol (0.9*1.1, 0.9*1.1, 0.7*1.1);
      break;
    case 19:
      pcol (0.6*1.1, 0.6*1.1, 0.6*1.1);
      break;
    case 20:
      pcol (0.6*1.1, 0.3*1.1, 0.1*1.1);
      break;
    case 21:
      pcol (0.6*1.1, 0.1*1.1, 0.3*1.1);
      break;
    case 22:
      pcol (0.3*1.1, 0.1*1.1, 0.6*1.1);
      break;
    case 23:
      pcol (0.3*1.1, 0.6*1.1, 0.1*1.1);
      break;
    case 24:
      pcol (0.1*1.1, 0.3*1.1, 0.6*1.1);
      break;
    case 25:
      pcol (0.1*1.1, 0.6*1.1, 0.3*1.1);
      break;
    case 99:
      pcol (0.4*1.1, 0.9*1.1, 0.3*1.1);
      break;
    case 100:
      pcol (0.2*1.1, 0.99*1.1, 0.0*1.1);
      break;
    default:
      pcol (0.9*1.1, 0.9*1.1, 0.9*1.1);
      break;
  }
}


/******************************************************************************
Draw lines connecting nearby points between meshes.
******************************************************************************/

draw_near_lines()
{
  int i;
  Mesh *m1,*m2;
  Vertex *near;
  Vector vec1,vec1w,vec2,vec3;
  Vector norm1,norm1w;
  float dist;

  pushmatrix();

  /* red lines */
  cpack (0xff);

  /* shorthand for two meshes */
  m1 = scans[0]->meshes[mesh_level];
  m2 = scans[1]->meshes[mesh_level];

  /* examine all vertices in first mesh for proximity to second mesh */

  for (i = 0; i < m1->nverts; i++) {

    /* find nearest point on other mesh */

    if (m1->verts[i]->ntris == 0) /* forget if vert not part of any triangle */
      continue;
    mesh_to_world (scans[0], m1->verts[i]->coord, vec1w);
    world_to_mesh (scans[1], vec1w, vec1);

    mesh_to_world_normal (scans[0], m1->verts[i]->normal, norm1w);
    world_to_mesh_normal (scans[1], norm1w, norm1);

    near = find_nearest (m2, NULL, vec1, norm1, FIND_COS);
    if (near == NULL)
      continue;

    /* if it is near enough, draw a line */
    if (near->ntris == 0) /* forget if vert not part of any triangle */
      continue;
    mesh_to_world (scans[1], near->coord, vec2);
    vsub (vec1w, vec2, vec3);
    dist = vlen (vec3);

    if (dist < 0.01) {
      bgnline();
      v3f (vec1w);
      v3f (vec2);
      endline();
    }
  }

  popmatrix();
}


/******************************************************************************
Draw lines connecting nearby points between meshes.
******************************************************************************/

draw_nearest_lines()
{
  int i;
  Mesh *m1,*m2;
  Vector vec1;
  Vector norm1;
  int result;
  NearPosition near_info;
  int inc;

  inc = level_to_inc (mesh_level);

  pushmatrix();

  /* green lines */
  cpack (0xff00);

  /* shorthand for two meshes */
  m1 = scans[0]->meshes[mesh_level];
  m2 = scans[1]->meshes[mesh_level];

  /* examine all vertices in first mesh for proximity to second mesh */

  for (i = 0; i < m1->nverts; i++) {

    /* find nearest point on other mesh */

    if (m1->verts[i]->ntris == 0) /* forget if vert not part of any triangle */
      continue;
    mesh_to_world (scans[0], m1->verts[i]->coord, vec1);
    mesh_to_world_normal (scans[0], m1->verts[i]->normal, norm1);
    result = nearest_on_mesh (scans[1], scans[1]->meshes[mesh_level], NULL,
			      vec1, norm1, inc*ZIPPER_RESOLUTION, FIND_COS, &near_info);

    /* draw the line between these position if we got a nearby point */
    if (result && !near_info.on_edge) {
      bgnline();
      v3f (vec1);
      v3f (near_info.pos);
      endline();
    }
  }

  popmatrix();
}


/******************************************************************************
Transform point from meshes space to world space.

Entry:
  sc    - scan of mesh
  invec - point to transform

Exit:
  outvec - transformed point
******************************************************************************/

mesh_to_world(sc,invec,outvec)
  Scan *sc;
  Vector invec,outvec;
{
  int i;
  Vector v;

  for (i = 0; i < 3; i++)
    v[i] = invec[X] * sc->rotmat[X][i] + invec[Y] * sc->rotmat[Y][i] +
	   invec[Z] * sc->rotmat[Z][i];

  outvec[X] = v[X] + sc->xtrans;
  outvec[Y] = v[Y] + sc->ytrans;
  outvec[Z] = v[Z] + sc->ztrans;
}


/******************************************************************************
Transform point from world space to a meshes space.

Entry:
  sc    - scan of mesh
  invec - point to transform

Exit:
  outvec - transformed point
******************************************************************************/

world_to_mesh(sc,invec,outvec)
  Scan *sc;
  Vector invec,outvec;
{
  int i;
  Vector v,w;

  v[X] = invec[X] - sc->xtrans;
  v[Y] = invec[Y] - sc->ytrans;
  v[Z] = invec[Z] - sc->ztrans;

  for (i = 0; i < 3; i++)
    outvec[i] = v[X] * sc->rotmat[i][X] + v[Y] * sc->rotmat[i][Y] +
	        v[Z] * sc->rotmat[i][Z];
}


/******************************************************************************
Transform normal from world space to a meshes space.

Entry:
  sc      - scan of mesh
  in_norm - surface normal to transform

Exit:
  out_norm - transformed point
******************************************************************************/

world_to_mesh_normal(sc,in_norm,out_norm)
  Scan *sc;
  Vector in_norm,out_norm;
{
  int i;
  Vector v;

  v[X] = in_norm[X];
  v[Y] = in_norm[Y];
  v[Z] = in_norm[Z];

  for (i = 0; i < 3; i++)
    out_norm[i] = v[X] * sc->rotmat[i][X] + v[Y] * sc->rotmat[i][Y] +
		  v[Z] * sc->rotmat[i][Z];
}


/******************************************************************************
Transform normal from meshes space to world space.

Entry:
  sc      - scan of mesh
  in_norm - point to transform

Exit:
  out_norm - transformed point
******************************************************************************/

mesh_to_world_normal(sc,in_norm,out_norm)
  Scan *sc;
  Vector in_norm,out_norm;
{
  int i;
  Vector v;

  for (i = 0; i < 3; i++)
    v[i] = in_norm[X] * sc->rotmat[X][i] + in_norm[Y] * sc->rotmat[Y][i] +
	   in_norm[Z] * sc->rotmat[Z][i];

  out_norm[X] = v[X];
  out_norm[Y] = v[Y];
  out_norm[Z] = v[Z];
}


/******************************************************************************
Toggle flag saying whether to draw polygons or wireframe.
******************************************************************************/

toggle_polygons()
{
  poly_draw_flag = 1 - poly_draw_flag;
  draw_object();
}


/******************************************************************************
Toggle flag saying whether to draw smooth-shaded polygons.
******************************************************************************/

toggle_smooth()
{
  smooth_flag = 1 - smooth_flag;
  draw_object();
}


/******************************************************************************
Toggle flag saying whether to use different colors for different meshes.
******************************************************************************/

toggle_colors()
{
  colors_flag = 1 - colors_flag;
  draw_object();
}


/******************************************************************************
Toggle flag saying whether to draw lines showing nearness of two meshes.
******************************************************************************/

toggle_near()
{
  near_flag = 1 - near_flag;
  draw_object();
}


/******************************************************************************
Toggle another flag similar to above.
******************************************************************************/

toggle_nearest()
{
  nearest_flag = 1 - nearest_flag;
  draw_object();
}


/******************************************************************************
Toggle if polygons on the edge of a mesh should be drawn.
******************************************************************************/

toggle_edge()
{
  edge_flag = 1 - edge_flag;
  draw_object();
}


/******************************************************************************
Set one of various drawing flags.
******************************************************************************/

int tcl_set_drawing(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  char *str;
  int flag;

  if (argc != 3) {
    return (TCL_ERROR);
  }

  /* set one of the drawing flags */

  Tcl_GetBoolean (interp, argv[2], &flag);

  if (strcmp (argv[1], "polygons") == 0)
    poly_draw_flag = flag;
  else if (strcmp (argv[1], "normal") == 0)
    normal_flag = flag;
  else if (strcmp (argv[1], "smooth") == 0)
    smooth_flag = flag;
  else if (strcmp (argv[1], "colors") == 0)
    colors_flag = flag;
  else if (strcmp (argv[1], "true_color") == 0)
    true_color_flag = flag;
  else if (strcmp (argv[1], "near") == 0)
    near_flag = flag;
  else if (strcmp (argv[1], "nearest") == 0)
    nearest_flag = flag;
  else if (strcmp (argv[1], "edge") == 0)
    edge_flag = flag;
  else if (strcmp (argv[1], "bounds") == 0)
    bound_flag = flag;
  else if (strcmp (argv[1], "refine") == 0)
    refine_flag = flag;
  else if (strcmp (argv[1], "value") == 0)
    value_flag = flag;
  else if (strcmp (argv[1], "intensity") == 0)
    intensity_flag = flag;
  else {
    return (TCL_ERROR);
  }

  draw_object();

  return (TCL_OK);
}


/******************************************************************************
Switch to another level-of-detail for a mesh, but save the old value.

Entry:
  level - new value of mesh level-of-detail
******************************************************************************/

save_mesh_level(level)
  int level;
{
  pushed_mesh_level = mesh_level;

  /* only switch levels if we are doing progressive refinement */
  if (refine_flag) {
    mesh_level = level;
    partial_flag = 1;
  }

  drew_partial = 0;
}


/******************************************************************************
Restore the mesh level-of-detail to the old, saved value.  Re-draw objects
if we are actually changing the value.
******************************************************************************/

restore_mesh_level()
{
  partial_flag = 0;

  if (mesh_level != pushed_mesh_level) {
    mesh_level = pushed_mesh_level;
    draw_object();
  }
  else if (drew_partial) {
    draw_object();
  }

  drew_partial = 0;
}


/******************************************************************************
Middle mouse button down in drawing window.
******************************************************************************/

void middle_proc()
{
  int i,j;
  int x,y;
  int xold,yold;
/*
  float delta = 0.002;
  float tdelta = 0.0002;
*/
  float delta = 0.0002;
  float tdelta = 0.0002;
  Matrix tmat;
  Matrix rotinv;

  /* see if we're doing picking */
  if (pick_flag) {
    char str[80];
    sprintf (str, "execute_edit");
    Tcl_Eval (interp, str, 0, (char **) NULL);
    return;
  }

  delta *= mouse_speed;
  tdelta *= mouse_speed;

  /* exit if we're in light-moving mode */
  if (move_num == MOVE_LIGHT)
    return;

  /* switch to coarse mesh level */
  choose_mesh_level();

  /* make inverse of rotmat */
  mat_copy (rotinv, rotmat);
  mat_transpose (rotinv);

  /* store position in undo stack if we're moving one mesh */
  if (move_num >= 0)
    save_for_undo (scans[move_num]);

  /* find cursor location */
  winset (draw_win);
  Cursor_Setup;
  Cursor_Pos (x, y);
  xold = x;
  yold = y;

  while (getbutton (MIDDLEMOUSE)) {

    Cursor_Pos (x, y);
    if ((xold == x) && (yold == y))
      continue;

    if (move_num >= 0) {

      /* move one mesh */

      pushmatrix();
	loadmatrix (rotinv);
	translate ((xold - x) * tdelta, (y - yold) * tdelta, 0.0);
	multmatrix (rotmat);
	translate (scans[move_num]->xtrans, scans[move_num]->ytrans,
		   scans[move_num]->ztrans);
	getmatrix (tmat);
      popmatrix();

      scans[move_num]->xtrans = tmat[3][0];
      scans[move_num]->ytrans = tmat[3][1];
      scans[move_num]->ztrans = tmat[3][2];
    }
    else {

      /* translate all meshes together */

      pushmatrix();
      loadmatrix (matident);
      translate ((xold - x) * delta, (y - yold) * delta, 0.0);
      multmatrix (transmat);
      getmatrix (transmat);
      popmatrix();
    }

    /* draw the object */
    draw_object();

    xold = x;
    yold = y;
  }

  /* restore earlier mesh level */
  restore_mesh_level();
}


/******************************************************************************
Left mouse button down in drawing window.
******************************************************************************/

void left_proc()
{
  int x,y;
  int x_old,y_old;
  float s,t;
  float s_old,t_old;
  float r[4];
  Matrix tmat;
  Matrix rotinv;
  Quaternion rvec;

  /* see if we're doing picking */
  if (pick_flag) {
    pick_polygon();
    return;
  }

  /* switch to coarse mesh level */
  choose_mesh_level();

  /* make inverse of rotmat */
  mat_copy (rotinv, rotmat);
  mat_transpose (rotinv);

  /* store position in undo stack if we're moving one mesh */
  if (move_num >= 0)
    save_for_undo (scans[move_num]);

  /* find cursor location */
  winset (draw_win);
  Cursor_Setup;
  Cursor_Pos (x, y);
  s = (2 * x - xsize) / (float) xsize;
  t = (2 * y - ysize) / (float) ysize;
  t = -t;
  s *= mouse_speed;
  t *= mouse_speed;
  s_old = s;
  t_old = t;
  x_old = x;
  y_old = y;

  while (getbutton (LEFTMOUSE)) {

    Cursor_Pos (x, y);
    s = (2 * x - xsize) / (float) xsize;
    t = (2 * y - ysize) / (float) ysize;
    t = -t;
    s *= mouse_speed;
    t *= mouse_speed;
    if ((x_old == x) && (y_old == y))
      continue;

    if (move_num >= 0) {

      /*** rotate one mesh ***/

      trackball (r, s_old, t_old, s, t);
      quat_to_mat (r, tmat);

      pushmatrix();
	loadmatrix (rotinv);
	multmatrix (tmat);
	multmatrix (rotmat);
	multmatrix (scans[move_num]->rotmat);
	getmatrix (scans[move_num]->rotmat);
      popmatrix();

    }
    else if (move_num == MOVE_WORLD) {

      /*** below is for global viewer position ***/

      /* build new rotation matrix from accumulated quaternion rotation */
      mat_to_quat (rotmat, rvec);
      trackball (r, s_old, t_old, s, t);
      add_quats (r, rvec, rvec);

      /* the following stomps old rotation matrix!!! */
      quat_to_mat (rvec, rotmat);
    }
    else if (move_num == MOVE_LIGHT) {
      my_light1[5] = 2 * x / (float) xsize - 1.0;
      my_light1[6] = 2 * y / (float) ysize - 1.0;
    }

    /* draw the object */
    draw_object();

    s_old = s;
    t_old = t;
    x_old = x;
    y_old = y;
  }

  /* restore earlier mesh level */
  restore_mesh_level();
}



/******************************************************************************
Right mouse button down in drawing window.
******************************************************************************/

void right_proc()
{
  int x,y;
  int xold,yold;
/*
  float delta = 0.005;
  float tdelta = 0.0002;
*/
  float delta = 0.00025;
  float tdelta = 0.00025;
  Matrix rotinv;
  Matrix tmat;

  delta *= mouse_speed;
  tdelta *= mouse_speed;

  /* switch to coarse mesh level */
  choose_mesh_level();

  /* make inverse of rotmat */
  mat_copy (rotinv, rotmat);
  mat_transpose (rotinv);

  /* exit if we're in light-moving mode */
  if (move_num == MOVE_LIGHT)
    return;

  /* store position in undo stack if we're moving one mesh */
  if (move_num >= 0)
    save_for_undo (scans[move_num]);

  /* find cursor location */
  winset (draw_win);
  Cursor_Setup;
  Cursor_Pos (x, y);
  xold = x;
  yold = y;

  while (getbutton (RIGHTMOUSE)) {

    Cursor_Pos (x, y);
    if ((xold == x) && (yold == y))
      continue;

    /* translate the object */

    if (move_num >= 0) {

      /* move one mesh forward and back */

      pushmatrix();
	loadmatrix (rotinv);
	translate (0.0, 0.0, (yold - y) * tdelta);
	multmatrix (rotmat);
	translate (scans[move_num]->xtrans, scans[move_num]->ytrans,
		   scans[move_num]->ztrans);
	getmatrix (tmat);
      popmatrix();

      scans[move_num]->xtrans = tmat[3][0];
      scans[move_num]->ytrans = tmat[3][1];
      scans[move_num]->ztrans = tmat[3][2];
    }
    else {

      /* move all meshes forward and back */

      pushmatrix();
      loadmatrix (matident);
      translate (0.0, 0.0, (yold - y) * delta);
      multmatrix (transmat);
      getmatrix (transmat);
      popmatrix();
    }

    /* draw the object */
    draw_object();

    xold = x;
    yold = y;
  }

  /* restore earlier mesh level */
  restore_mesh_level();
}


/******************************************************************************
Redraw event for drawing window.
******************************************************************************/

void refresh_drawing()
{
  draw_object();
}


/******************************************************************************
Quit the program.
******************************************************************************/

void quit_proc()
{
  exit (0);
}


/******************************************************************************
Initialize the data structures for "extra" lines to be displayed each frame.
******************************************************************************/

init_extra_lines()
{
  /* maybe free up old space */
  if (line1 != NULL) {
    free (line1);
    free (line2);
    free (line_colors);
    line1 = NULL;
    line2 = NULL;
    line_colors = NULL;
  }

  line_num = 0;
}


/******************************************************************************
Add a line to the "extra" lines to be displayed each frame.

Entry:
  p1,p2 - endpoints of line segment to be displayed
  color - line color
******************************************************************************/

add_extra_line(p1,p2,color)
  Vector p1,p2;
  unsigned int color;
{
  /* create new space for lines if necessary */
  if (line1 == NULL) {
    line_max = 200;
    line1 = (Vector *) myalloc (sizeof (Vector) * line_max);
    line2 = (Vector *) myalloc (sizeof (Vector) * line_max);
    line_colors = (unsigned int *) myalloc (sizeof (unsigned int) * line_max);
  }

  /* reallocate more space if necessary */

  if (line_num == line_max) {
    line_max *= 2;
    line1 = (Vector *) realloc (line1, sizeof (Vector) * line_max);
    line2 = (Vector *) realloc (line2, sizeof (Vector) * line_max);
    line_colors = (unsigned int *)
		  realloc (line_colors, sizeof (unsigned int) * line_max);
  }

  /* add the new line to the list */
  vcopy (p1, line1[line_num]);
  vcopy (p2, line2[line_num]);
  line_colors[line_num] = color;
  line_num++;
}


/******************************************************************************
Draw the "extra" lines for a frame.
******************************************************************************/

draw_extra_lines()
{
  int i;
  unsigned int old_color = 0xffffffff;

  cpack (old_color);

  for (i = 0; i < line_num; i++) {
    if (line_colors[i] != old_color) {
      old_color = line_colors[i];
      cpack (old_color);
    }
    bgnline();
    v3f (line1[i]);
    v3f (line2[i]);
    endline();
  }
}


/******************************************************************************
Print how many extra lines were created.
******************************************************************************/

print_extra_lines()
{
  printf ("%d extra lines\n", line_num);
}


/******************************************************************************
Set the flag that says whether or not to draw axes.
******************************************************************************/

set_axes(val)
  int val;
{
  draw_axes = val;
}


/******************************************************************************
Write the current image to a file.

Entry:
  initial - initial letters of file name
  num     - frame number
******************************************************************************/

old_write_frame(initial,num)
  char *initial;
  int num;
{
  int i,j;
  IMAGE *oimage;
  char outfile[80];
  static short sbuf[4096];
  unsigned long *scrn;
  int base;

  /* allocate room to save screen */

  scrn = (unsigned long *)
         malloc (sizeof (long) * xsize * ysize);

  /* read pixels from window */

  winset (draw_win);
  lrectread (0, 0, xsize - 1, ysize - 1, scrn);

  /* open image file */

  printf ("writing to %s%03d\n", initial, num);
  sprintf (outfile, "%s%03d", initial, num);
  oimage = iopen (outfile, "w", RLE(1), 3, xsize, ysize, 3);

  /* write out image line-by-line */

  base = 0;
  for (j = 0; j < ysize; j++) {

    /* red portion of line */
    for (i = 0; i < xsize; i++)
      sbuf[i] = (scrn[i+base] & 0x000000ff);
    putrow (oimage, sbuf, j, 0);

    /* green portion of line */
    for (i = 0; i < xsize; i++)
      sbuf[i] = (scrn[i+base] & 0x0000ff00) >> 8;
    putrow (oimage, sbuf, j, 1);

    /* blue portion of line */
    for (i = 0; i < xsize; i++)
      sbuf[i] = (scrn[i+base] & 0x00ff0000) >> 16;
    putrow (oimage, sbuf, j, 2);

    base += xsize;
  }

  free (scrn);
  iclose (oimage);
}


/******************************************************************************
Write the current image to a file.

Entry:
  initial - initial letters of file name
  num     - frame number
******************************************************************************/

write_frame(initial,num)
  char *initial;
  int num;
{
  int i,j;
  IMAGE *oimage;
  char outfile[80];
  static short sbuf[4096];
  unsigned long *scrn;

  /* allocate room to save screen */

  scrn = (unsigned long *)
         malloc (sizeof (long) * xsize);

  /* read pixels from window */

  winset (draw_win);

  /* open image file */

  printf ("writing to %s%03d\n", initial, num);
  sprintf (outfile, "%s%03d", initial, num);
  oimage = iopen (outfile, "w", RLE(1), 3, xsize, ysize, 3);

  /* write out image line-by-line */

  for (j = 0; j < ysize; j++) {

    lrectread (0, j, xsize - 1, j, scrn);

    /* red portion of line */
    for (i = 0; i < xsize; i++)
      sbuf[i] = (scrn[i] & 0x000000ff);
    putrow (oimage, sbuf, j, 0);

    /* green portion of line */
    for (i = 0; i < xsize; i++)
      sbuf[i] = (scrn[i] & 0x0000ff00) >> 8;
    putrow (oimage, sbuf, j, 1);

    /* blue portion of line */
    for (i = 0; i < xsize; i++)
      sbuf[i] = (scrn[i] & 0x00ff0000) >> 16;
    putrow (oimage, sbuf, j, 2);
  }

  free (scrn);
  iclose (oimage);
}




set_singlebuffer(int val) {

    singlebuf = val;
    if (singlebuf) {
	singlebuffer();
	gconfig();
	draw_object();
    } else {
	doublebuffer();
	gconfig();
	draw_object();
    }
}


