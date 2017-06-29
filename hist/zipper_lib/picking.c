/*

Pick parts of an object.

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
#include <GL/gl.h>
//#include <device.h>
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

extern int edge_flag;
extern int pick_flag;
extern int copied_to_back;
extern int bits24;

/* picking memory */
static Scan *pick_scan = NULL;		/* scan of the picked item */
static Vertex *pick_vert1 = NULL;	/* first picked vertex */
static Vertex *pick_vert2 = NULL;	/* second picked vertex */
static Vertex *pick_vert3 = NULL;	/* third picked vertex */
static Triangle *pick_tri = NULL;	/* the picked triangle, if any */
static int pick_loop = -1;		/* the picked boundary, if any */
static int item_picked = 0;		/* is an item currently picked? */
static int box_x0,box_y0;		/* one corner of "delete" box */
static int box_x1,box_y1;		/* other corner of "delete" box */

#define PICK_TRI      1
#define PICK_VERT     2
#define PICK_BOUNDARY 3
#define PICK_BOUND2   4
#define PICK_3_VERTS  5
#define PICK_BOX      6

/******************************************************************************
Pick a polygon or vertex with the cursor.
******************************************************************************/

pick_polygon()
{
  int i,j;
  Scan *sc;
  Mesh *mesh;
  int sc_index;
  float dist,max;
  int index,max_index;
  Triangle *tri;
  Vertex *v;

  winset (draw_win);

  /* box picking needs to be handled specially */
  if (pick_flag == PICK_BOX) {
    pick_box();
    return;
  }

  /* wait until user lifts up on mouse button */
  while (getbutton (LEFTMOUSE)) ;

  /* examine each scan for picked items */
  sc_index = -1;
  max = -1e20;
  for (i = 0; i < nscans; i++) {

    /* shorthands */
    sc = scans[i];
    mesh = sc->meshes[mesh_level];

    /* pick from the mesh if the flag says the mesh is visible */
    if (sc->draw_flag) {

      if (pick_flag == PICK_TRI)
	pick_tri_from_mesh (sc, mesh, &index, &dist);
      else if (pick_flag == PICK_VERT || pick_flag == PICK_3_VERTS)
	pick_vert_from_mesh (sc, mesh, &index, &dist);
      else if (pick_flag == PICK_BOUNDARY || pick_flag == PICK_BOUND2) {
	pick_boundary_from_mesh (sc, mesh, &index);
	dist = 0;
      }

      /* see if this item is nearer than others */
      if (index > -1 && dist > max) {
	sc_index = i;
	max_index = index;
	max = dist;
      }
    }
  }

  /* highlight the selected triangle or vertex */

  if (sc_index > -1) {
    mesh = scans[sc_index]->meshes[mesh_level];
    if (pick_flag == PICK_TRI) {
      printf ("triangle %d picked on mesh '%s'\n", max_index,
	       scans[sc_index]->name);
      tri = mesh->tris[max_index];
      highlight_tri (tri, scans[sc_index]);
      pick_tri = tri;
      pick_scan = scans[sc_index];
      item_picked = 1;
    }
    else if (pick_flag == PICK_VERT) {
      v = mesh->verts[max_index];
      vertex_edge_test (v);
      printf ("vertex %d picked on mesh '%s'\n", max_index,
	       scans[sc_index]->name);
      highlight_vert (v, scans[sc_index], 1);
      pick_vert1 = v;
      pick_scan = scans[sc_index];
      item_picked = 1;
    }
    else if (pick_flag == PICK_BOUNDARY || pick_flag == PICK_BOUND2) {
      printf ("boundary %d picked on mesh '%s'\n", max_index,
	      scans[sc_index]->name);
      highlight_bounary (max_index, scans[sc_index]);
      pick_loop = max_index;
      pick_scan = scans[sc_index];
      item_picked = 1;
    }
    else if (pick_flag == PICK_3_VERTS) {
      v = mesh->verts[max_index];
      vertex_edge_test (v);
      printf ("vertex %d picked on mesh '%s'\n", max_index,
	       scans[sc_index]->name);
      pick_vert3 = pick_vert2;
      pick_vert2 = pick_vert1;
      pick_vert1 = v;
      highlight_vert (pick_vert1, scans[sc_index], 1);
      highlight_vert (pick_vert2, scans[sc_index], 0);
      highlight_vert (pick_vert3, scans[sc_index], 0);
      pick_scan = scans[sc_index];
      item_picked = 1;
    }
  }
  else {
    printf ("nothing picked\n");
  }
}


/******************************************************************************
Delete the currently picked item, if there is one.
******************************************************************************/

delete_picked_item()
{
  if (item_picked && pick_scan) {
    if (pick_vert1 && pick_flag == PICK_VERT) {
      remove_a_vertex (pick_scan, pick_vert1);
      pick_scan->meshes[mesh_level]->edges_valid = 0;
    }
    else if (pick_tri) {
      delete_triangle (pick_tri, pick_scan->meshes[mesh_level], 1);
      pick_scan->meshes[mesh_level]->edges_valid = 0;
    }
    else
      fprintf (stderr, "delete_picked_item: can't determine what to delete\n");
  }
  else
    printf ("no picked item to delete\n");
}


/******************************************************************************
Fill in the currently picked loop, if there is one.
******************************************************************************/

fill_picked_loop()
{
  int result;

  if (item_picked && pick_scan && pick_loop >= 0) {
    result = fill_loop (pick_loop, pick_scan);
    if (result)
      fprintf (stderr, "can't fill this loop\n");
  }
  else
    printf ("no picked loop to fill\n");
}


/******************************************************************************
Use more elaborate method to fill in the currently picked loop.
******************************************************************************/

better_fill_picked_loop()
{
  if (item_picked && pick_scan && pick_loop >= 0)
    better_fill_loop (pick_loop, pick_scan);
  else
    printf ("no picked loop to fill\n");
}


/******************************************************************************
Pick the nearest triangle from a given mesh.

Entry:
  sc   - scan of given mesh
  mesh - mesh in question

Exit:
  index - index of nearest polygon picked (or -1 if nothing picked)
  dist  - distance of nearest polygon
******************************************************************************/

pick_tri_from_mesh(sc,mesh,index,dist)
  Scan *sc;
  Mesh *mesh;
  int *index;
  float *dist;
{
  int i,j;
  Vertex **v;
  int npicks;
  static short int pick_buf[99];
  float max;
  int max_index;
  Matrix mat;
  Vector vec;
  int ind;

  /* pass all triangles through the picker */

  picksize (1, 1);
  pick (pick_buf, 99);

  /* set up transformation */

  pushmatrix();
  pushviewport();

  set_global_transformation();
  translate (sc->xtrans, sc->ytrans, sc->ztrans);
  multmatrix (sc->rotmat);
  getmatrix (mat);

  /* draw all triangles */
  for (j = 0; j < mesh->ntris; j++) {

    v = mesh->tris[j]->verts;

    /* maybe don't draw polygons on edge of mesh */
    if (edge_flag && (v[0]->on_edge || v[1]->on_edge || v[2]->on_edge))
      continue;

    ind = mesh->tris[j]->index;
    loadname (ind % 1000);
    pushname (ind / 1000);
    bgnpolygon();
    v3f (v[0]->coord);
    v3f (v[1]->coord);
    v3f (v[2]->coord);
    v3f (v[0]->coord);
    endpolygon();
    popname();
  }

  popmatrix();
  popviewport();

  npicks = endpick (pick_buf);

  /* look for nearest triangle in pick buffer */
  max = -1e20;
  max_index = -1;
  for (i = 0; i < npicks; i++) {
    ind = pick_buf[3*i+1] + 1000 * pick_buf[3*i+2];
    for (j = 0; j < 3; j++) {
      vcopy (mesh->tris[ind]->verts[j]->coord, vec);
      mat_apply (mat, vec);
      if (vec[Z] > max) {
	max = vec[Z];
	max_index = ind;
      }
    }
  }

  /* return pick values */
  *index = max_index;
  *dist = max;
}


/******************************************************************************
Pick the nearest vertex from a given mesh.

Entry:
  sc   - scan of given mesh
  mesh - mesh in question

Exit:
  index - index of nearest vertex picked (or -1 if nothing picked)
  dist  - distance of nearest vertex
******************************************************************************/

pick_vert_from_mesh(sc,mesh,index,dist)
  Scan *sc;
  Mesh *mesh;
  int *index;
  float *dist;
{
  int i,j;
  Vertex *v;
  int npicks;
  static short int pick_buf[99];
  float max;
  int max_index;
  Matrix mat;
  Vector vec;
  int ind;

  /* pass all vertices through the picker */

  picksize (6, 6);
  pick (pick_buf, 99);

  /* set up transformation */

  pushmatrix();
  pushviewport();

  set_global_transformation();
  translate (sc->xtrans, sc->ytrans, sc->ztrans);
  multmatrix (sc->rotmat);
  getmatrix (mat);

  /* draw all vertices */
  for (j = 0; j < mesh->nverts; j++) {

    v = mesh->verts[j];
    ind = v->index;
    loadname (ind % 1000);
    pushname (ind / 1000);
    bgnpoint();
    v3f (v->coord);
    endpoint();
    popname();
  }

  popmatrix();
  popviewport();

  npicks = endpick (pick_buf);

  /* look for nearest vertex in pick buffer */
  max = -1e20;
  max_index = -1;
  for (i = 0; i < npicks; i++) {
    ind = pick_buf[3*i+1] + 1000 * pick_buf[3*i+2];
    for (j = 0; j < 3; j++) {
      vcopy (mesh->verts[ind]->coord, vec);
      mat_apply (mat, vec);
      if (vec[Z] > max) {
	max = vec[Z];
	max_index = ind;
      }
    }
  }

  /* return pick values */
  *index = max_index;
  *dist = max;
}


/******************************************************************************
Pick the boundary of edges from a given mesh.

Entry:
  sc   - scan of given mesh
  mesh - mesh in question

Exit:
  index - index of nearest boundary loop picked (or -1 if nothing picked)
******************************************************************************/

pick_boundary_from_mesh(sc,mesh,index)
  Scan *sc;
  Mesh *mesh;
  int *index;
{
  int i,j;
  Vertex *v;
  int npicks;
  static short int pick_buf[99];
  float max;
  int max_index;
  Matrix mat;
  Vector vec;
  int ind;
  Edge *e,*fedge;
  int been_around;

  /* pass all boundaries through the picker */

  picksize (6, 6);
  pick (pick_buf, 99);

  /* set up transformation */

  pushmatrix();
  pushviewport();

  set_global_transformation();
  translate (sc->xtrans, sc->ytrans, sc->ztrans);
  multmatrix (sc->rotmat);
  getmatrix (mat);

  /* make sure the boundaries are up-to-date */
  if (!mesh->edges_valid)
    create_edge_list (mesh);

  /* draw all boundaries */
  for (i = 0; i < mesh->looplist.nloops; i++) {

    loadname (i);

    fedge = mesh->looplist.loops[i];
    been_around = 0;

    for (e = fedge; e != fedge || !been_around; e = e->next) {
      been_around = 1;
      bgnline();
      v3f (e->v1->coord);
      v3f (e->v2->coord);
      endline();
    }
  }

  popmatrix();
  popviewport();

  npicks = endpick (pick_buf);

  /* return pick value */
  if (npicks)
    *index = pick_buf[1];
  else
    *index = -1;
}


/******************************************************************************
Highlight a triangle in the overlay plane.
******************************************************************************/

highlight_tri(tri,sc)
  Triangle *tri;
  Scan *sc;
{
  begin_overdraw (1);

  /* set up the appropriate transformations */
  pushmatrix();
  set_global_transformation();
  translate (sc->xtrans, sc->ytrans, sc->ztrans);
  multmatrix (sc->rotmat);

  /* draw the triangle */

  bgnpolygon();
  v3f (tri->verts[0]->coord);
  v3f (tri->verts[1]->coord);
  v3f (tri->verts[2]->coord);
  v3f (tri->verts[0]->coord);
  endpolygon();

  /* reset the matrices */
  popmatrix();

  end_overdraw();
}


/******************************************************************************
Highlight a vertex in the overlay plane.

Entry:
  vert - pointer to vertex
  sc   - scan containing vertex
  clr  - flag saying whether to clear the overlay planes, 1 = clear
******************************************************************************/

highlight_vert(vert,sc,clr)
  Vertex *vert;
  Scan *sc;
  int clr;
{
  Matrix mat;
  long int xorg,yorg;
  long int xsize,ysize;
  Vector vec;

  /* don't do anything if "vert" is NULL */
  if (!vert)
    return;

  object_to_screen (sc, mat);

  begin_overdraw (clr);

  /* set up the appropriate transformations */
  getsize (&xsize, &ysize);
  getorigin (&xorg, &yorg);
  pushmatrix();
  loadmatrix (matident);
  ortho2 (xorg - 0.5, xorg + xsize + 0.5,
          yorg - 0.5, yorg + ysize + 0.5);

  /* draw a circle at the vertex */
  vcopy (vert->coord, vec);
  mat_apply_homog (mat, vec);
  vec[X] = (vec[X] + 1) * xsize / 2 + xorg;
  vec[Y] = (vec[Y] + 1) * ysize / 2 + yorg;
  circf (vec[X], vec[Y], 3.0);

  /* reset the matrices */
  popmatrix();

  end_overdraw();
}


/******************************************************************************
Return a matrix describing the transformation between object and screen
coordinates.

Entry:
  sc - scan that defines the object's space

Exit:
  mat - the matrix describing the transformation
******************************************************************************/

object_to_screen(sc,mat)
  Scan *sc;
  Matrix mat;
{
  Matrix obj_mat;
  Matrix view_mat;

  pushmatrix();
  set_global_transformation();
  translate (sc->xtrans, sc->ytrans, sc->ztrans);
  multmatrix (sc->rotmat);
  getmatrix (obj_mat);
  popmatrix();

  pushmatrix();
  mmode (MPROJECTION);
  getmatrix (view_mat);
  mmode (MVIEWING);
  popmatrix();

  pushmatrix();
  loadmatrix (view_mat);
  multmatrix (obj_mat);
  getmatrix (mat);
  popmatrix();
}


/******************************************************************************
Highlight a boundary using the overlay plane.

Entry:
  index - index of boundary loop to highlight
  sc    - scan containing bounary
******************************************************************************/

highlight_bounary(index,sc)
  int index;
  Scan *sc;
{
  int i;
  Mesh *mesh = sc->meshes[mesh_level];
  Edge *e,*fedge;
  int been_around;
  long old_width;

  begin_overdraw (1);

  /* set up the appropriate transformations */
  pushmatrix();
  set_global_transformation();
  translate (sc->xtrans, sc->ytrans, sc->ztrans);
  multmatrix (sc->rotmat);

  /* draw the boundary */
  old_width = getlwidth();
  linewidth (3);

  fedge = mesh->looplist.loops[index];
  been_around = 0;
  for (e = fedge; e != fedge || !been_around; e = e->next) {
    been_around = 1;
    bgnline();
    v3f (e->v1->coord);
    v3f (e->v2->coord);
    endline();
  }

  /* reset the matrices and line width */
  popmatrix();
  linewidth (old_width);

  end_overdraw();
}


/******************************************************************************
Specify that nothing is picked.
******************************************************************************/

nothing_picked()
{
  item_picked = 0;
  pick_scan = NULL;
  pick_vert1 = NULL;
  pick_vert2 = NULL;
  pick_vert3 = NULL;
  pick_tri = NULL;
  pick_loop = -1;

  /* clear the overlay plane */
  
  if (bits24) {
    drawmode (OVERDRAW);
    singlebuffer();
    gconfig();
    color (0);
    clear();
    drawmode (NORMALDRAW);
  }
}


/******************************************************************************
Make a triangle out of the three currently picked vertices, if there are such.
******************************************************************************/

make_picked_tri()
{
  int fore,back;
  Mesh *mesh;
  Vertex *v1,*v2,*v3;
  Scan *scan;

  if (!pick_vert1 || !pick_vert2 || !pick_vert3) {
    fprintf (stderr, "three vertices must be picked to make a triangle\n");
    return;
  }

  mesh = pick_scan->meshes[mesh_level];

  /* check if edges are okay in forward or reverse direction */
  fore = check_proposed_tri (pick_vert1, pick_vert2, pick_vert3);
  back = check_proposed_tri (pick_vert3, pick_vert2, pick_vert1);

  /* maybe make triangle */
  if (fore && !back) {
    make_triangle (mesh, pick_vert1, pick_vert2, pick_vert3, 1e20);
  }
  else if (!fore && back) {
    make_triangle (mesh, pick_vert3, pick_vert2, pick_vert1, 1e20);
  }
  else {
    fprintf (stderr, "triangle can't be oriented\n");
    return;
  }

  vertex_edge_test (pick_vert1);
  vertex_edge_test (pick_vert2);
  vertex_edge_test (pick_vert3);
  find_vertex_normal (pick_vert1);
  find_vertex_normal (pick_vert2);
  find_vertex_normal (pick_vert3);

  v1 = pick_vert1;
  v2 = pick_vert2;
  v3 = pick_vert3;
  scan = pick_scan;

  draw_object();

  pick_vert1 = v1;
  pick_vert2 = v2;
  pick_vert3 = v3;
  pick_scan = scan;
  item_picked = 1;

  highlight_vert (pick_vert1, pick_scan, 1);
  highlight_vert (pick_vert2, pick_scan, 0);
  highlight_vert (pick_vert3, pick_scan, 0);

  /* mark that the edges are not valid in this mesh */
  mesh->edges_valid = 0;
}


/******************************************************************************
Pick a bunch of triangles to delete with a box.
******************************************************************************/

pick_box()
{
  int x0,y0;
  int x1,y1;
  int xold,yold;
  Matrix mat;
  long xorg,yorg;
  long xsize,ysize;
  Vector vec;

  winset (draw_win);
  Cursor_Setup;
  Cursor_Pos (x0, y0);
  x0 += xorg;
  y0 += yorg;

  begin_overdraw (1);

  /* set up the appropriate transformations */
  getsize (&xsize, &ysize);
  getorigin (&xorg, &yorg);
  pushmatrix();
  loadmatrix (matident);
  ortho2 (xorg - 0.5, xorg + xsize + 0.5,
          yorg - 0.5, yorg + ysize + 0.5);

  /* draw a rubber-band box */
  while (getbutton (LEFTMOUSE)) {
    Cursor_Pos (x1, y1);
    x1 += xorg;
    y1 += yorg;
    if (x1 != xold || y1 != yold) {
#if 0
      color (0);
      clear();
      color (3);
#endif
      bgnline();
      vset (vec, (float) x0, (float) y0, 0.0);
      v3f (vec);
      vset (vec, (float) x1, (float) y0, 0.0);
      v3f (vec);
      vset (vec, (float) x1, (float) y1, 0.0);
      v3f (vec);
      vset (vec, (float) x0, (float) y1, 0.0);
      v3f (vec);
      vset (vec, (float) x0, (float) y0, 0.0);
      v3f (vec);
      endline();
    }
    xold = x1;
    yold = y1;
  }

  /* remember where this delete box is */
  item_picked = 1;
  box_x0 = x0;
  box_y0 = y0;
  box_x1 = x1;
  box_y1 = y1;

  /* reset stuff */
  popmatrix();
  end_overdraw();
}


/******************************************************************************
Delete those triangles that fall within the delete box.
******************************************************************************/

delete_in_box()
{
  int i,j,k;
  Matrix mat;
  long xorg,yorg;
  long xsize,ysize;
  Vector vec;
  Scan *sc;
  Mesh *mesh;
  Vertex *vert;
  int temp;
  float x,y;

  winset (draw_win);

  /* order the box boundaries */

  if (box_x0 > box_x1) {
    temp = box_x0;
    box_x0 = box_x1;
    box_x1 = temp;
  }
  if (box_y0 > box_y1) {
    temp = box_y0;
    box_y0 = box_y1;
    box_y1 = temp;
  }

  /* get screen size and origin */
  getsize (&xsize, &ysize);
  getorigin (&xorg, &yorg);

  /* examine each scan to see which vertices fall into box */

  for (i = 0; i < nscans; i++) {

    /* shorthands */
    sc = scans[i];
    mesh = sc->meshes[mesh_level];
    if (sc->draw_flag == 0)
      continue;

    /* get transformation from the object to the screen */
    object_to_screen (sc, mat);

    for (j = 0; j < mesh->nverts; j++) {

      /* figure out where vertex lies on the screen */
      vert = mesh->verts[j];
      vcopy (vert->coord, vec);
      mat_apply_homog (mat, vec);
      x = (vec[X] + 1) * xsize / 2 + xorg;
      y = (vec[Y] + 1) * ysize / 2 + yorg;

      /* if it lies within the delete box, delete its triangles */
      if (x >= box_x0 && x <= box_x1 && y >= box_y0 && y <= box_y1) {
	/* delete the triangles */
	for (k = vert->ntris - 1; k >= 0; k--)
	  delete_triangle (vert->tris[k], mesh, 0);
	/* mark that the edges are not valid in this mesh */
	mesh->edges_valid = 0;
      }
    }

    /* delete any unused vertices */
    remove_unused_verts (mesh);
  }

  draw_object();
}


/******************************************************************************
Get ready for simulating overdraw planes.  This is done by making a copy
of the current image in the backbuffer.

Entry:
  clear_flag - whether to clear the overdraw planes
******************************************************************************/

begin_overdraw(clear_flag)
  int clear_flag;
{
  static unsigned long *parray = NULL;

  winset (draw_win);

  if (bits24) {
    /* set drawing to overlay, single buffer */
    drawmode (OVERDRAW);
    singlebuffer();
    gconfig();

    mapcolor (1, 255, 255, 255);
    mapcolor (2, 255,   0,   0);
    mapcolor (3,   0, 255,   0);

    if (clear_flag) {
      color (0);
      clear();
    }
    color (3);
    return;
  }

  if (parray == NULL)
    parray = (unsigned long *) myalloc (sizeof (unsigned long) * xsize * ysize);

  /* see if we've already copied the current image to the backbuffer */
  if (!copied_to_back) {
    copied_to_back = 1;
    readsource (SRC_FRONT);
    frontbuffer (FALSE);
    backbuffer (TRUE);
    lrectread (0, 0, xsize-1, ysize-1, parray);
    lrectwrite (0, 0, xsize-1, ysize-1, parray);
  }

  /* copy the current image to the front buffer */
  frontbuffer (TRUE);
  backbuffer (FALSE);
  if (clear_flag) {
    readsource (SRC_BACK);
    lrectread (0, 0, xsize-1, ysize-1, parray);
    lrectwrite (0, 0, xsize-1, ysize-1, parray);
    readsource (SRC_AUTO);
  }

  /* set the drawing color to green */
  cpack (0x00ff00);
}


/******************************************************************************
End overdraw simulation.
******************************************************************************/

end_overdraw()
{
  winset (draw_win);

  if (bits24) {
    drawmode (NORMALDRAW);
    return;
  }

  /* re-set drawing into the backbuffer */
  frontbuffer (FALSE);
  backbuffer (TRUE);
}

