#include <cyfile.h>
#include <tcl.h>
#include <zipper.h>
#include <gl/gl.h>
#include <limits.h>
#include <gl/gl.h>
#include <gl/device.h>
#include <view.h>
#include <matrix.h>
#include "cursors.h"

extern int mesh_level;
extern int move_num;
extern int isatty();
extern Scan *find_scan();
extern float measure_error();
extern create_match_list();
extern Matrix transmat,rotmat;	/* viewing parameters */
extern int zipper_old;		/* use old zipper routines? */
extern int pick_flag;		/* enable picking? */
extern int back_cull;		/* cull backfacing outlines of triangles? */
extern float mouse_speed;	/* speed of object motion based on mouse */
extern int back_red,back_grn,back_blu;  /* background color */
extern int write_intensity_flag;        /* whether to write out intensities */
extern int align_on;            /* are we aligning? */
extern int draw_during_ops;     /* draw during slow operations? */
extern int write_normals_flag;  /* write out surface normals? */
extern int write_colors_flag;   /* write out colors? */
extern int parallel_procs_max;  /* maximum number of parallel processes */
extern int global_dont_draw;


Scan *target_scan = NULL; /* default target for zippering */
Scan *next_merge_scan = NULL; /* default source for zippering */


/******************************************************************************
Eat edges of redundant geometry.
******************************************************************************/

int tcl_eat_action(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  eat_edge_proc();
  return (TCL_OK);
}


/******************************************************************************
Do final stage of zipping.
******************************************************************************/

int tcl_only_merge(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  int i;
  Mesh *mesh;
  Scan *sc1,*sc2;

  if (argc != 1 && argc != 2 && argc != 3) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }

  if (argc == 1) {
    sc1 = target_scan;
    sc2 = next_merge_scan;
  }
  else if (argc == 2) {
    sc1 = target_scan;
    sc2 = find_scan (argv[1]);
  }
  else {
    sc1 = find_scan (argv[1]);
    sc2 = find_scan (argv[2]);
  }

  if (sc1 == NULL || sc1->meshes[mesh_level] == NULL) {
      interp->result = "target mesh is invalid";
      return (TCL_ERROR);
  }
  if (sc2 == NULL || sc2->meshes[mesh_level] == NULL) {
      interp->result = "merging mesh is invalid";
      return (TCL_ERROR);
  }

  if (!global_dont_draw)
      push_cursor (WAIT_CURSOR);

  if (zipper_old) {
    /* old zippering steps */
    zipper_meshes (sc1, sc2);
    fill_in_holes (sc1, sc2);
    move_vertices (sc2, sc1);
  }
  else {
    /* new zippering */
    fix_bows (sc1);
    fix_bows (sc2);
    create_edge_list (sc1->meshes[mesh_level]);
    gather_triangles (sc1, sc2);
    clip_triangles (sc1, sc2);

#if 1

#ifdef DEBUG_CLIP
  printf ("filling small holes...\n");
#endif

  fill_small_holes (sc1);

#ifdef DEBUG_CLIP
  printf ("removing cut vertices...\n");
#endif
  remove_cut_vertices (sc1);

#ifdef DEBUG_CLIP
  printf ("removing short edges...\n");
#endif
  remove_short_edges (sc1, 0.333);

  mesh = sc1->meshes[mesh_level];
  for (i = 0; i < mesh->ntris; i++)
    mesh->tris[i]->mark = 0;

#ifdef DEBUG_CLIP
  printf ("cleaning up mesh...\n");
#endif
  clean_up_mesh (sc1);

#endif

    /* Is this really necessary?  Appears so... */
  fill_small_holes (sc1);

  }

  if (argc == 1)
      next_merge_scan = NULL;

  if (!global_dont_draw)
      pop_cursor();

  draw_object();

  return (TCL_OK);
}


/******************************************************************************
Do final stage of zipping.
******************************************************************************/

int tcl_merge(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  int result;

  if (argc > 2) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }

/*
fprintf (stderr, "This new version of 'merge' also performs edge-eating.\n");
fprintf (stderr,
  "Use 'only_merge' for the version of merging that doesn't eat edges.\n");
*/

  /* eat edges */

  result = tcl_eat_edges (dummy, interp, argc, argv);
  if (result == TCL_ERROR)
    return (TCL_ERROR);

  /* perform zippering */

  result = tcl_only_merge (dummy, interp, argc, argv);

  return (result);
}




/******************************************************************************
Perform intersection (by mutual clipping) of two meshes.
******************************************************************************/

int tcl_intersect(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  int i;
  Mesh *mesh;
  Scan *sc1,*sc2;

  if (argc != 3) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }

  /* find the meshes */
  sc1 = find_scan (argv[1]);
  sc2 = find_scan (argv[2]);

  if (sc1 == NULL || sc2 == NULL) {
    interp->result = "can't find mesh name";
    return (TCL_ERROR);
  }

  intersect_meshes (sc1, sc2);

  draw_object();

  return (TCL_OK);
}


/******************************************************************************
Finish intersection (by mutual clipping) of two meshes.
******************************************************************************/

int tcl_ifinish(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  int i;
  Mesh *mesh;
  Scan *sc1,*sc2;

  if (argc != 3) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }

  /* find the meshes */
  sc1 = find_scan (argv[1]);
  sc2 = find_scan (argv[2]);

  if (sc1 == NULL || sc2 == NULL) {
    interp->result = "can't find mesh name";
    return (TCL_ERROR);
  }

  finish_intersect_meshes (sc1, sc2);

  draw_object();

  return (TCL_OK);
}


/******************************************************************************
Do all stages of zippering together: eat edges, share and merge.
******************************************************************************/

int tcl_zip(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  do_it_all();
  return (TCL_OK);
}



/******************************************************************************
Join two meshes into one.

Entry:
  argv[1] - first mesh
  argv[2] - second mesh
******************************************************************************/

int tcl_join(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  Scan *sc1,*sc2;

  if (argc != 3) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }

  /* find the meshes */
  sc1 = find_scan (argv[1]);
  sc2 = find_scan (argv[2]);

  if (sc1 == NULL || sc2 == NULL) {
    interp->result = "can't find mesh name";
    return (TCL_ERROR);
  }

  zipper_meshes (sc1, sc2);
  fill_in_holes (sc1, sc2);
  move_vertices (sc2, sc1);

  draw_object();

  return (TCL_OK);
}


/******************************************************************************
Repeatedly eat away at a pair of edges.

Entry:
  argv[1] - first mesh to eat away at
  argv[2] - second scan
******************************************************************************/

int tcl_eat_edges(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  Scan *sc1,*sc2;
  extern find_mesh_edges();
  extern lower_edge_confidence();

  if (argc > 3) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }

  if (argc == 1) {
    sc1 = target_scan;
    sc2 = next_merge_scan;
  }
  else if (argc == 2) {
    sc1 = target_scan;
    sc2 = find_scan (argv[1]);
  }
  else {
    /* find the meshes */
    sc1 = find_scan (argv[1]);
    sc2 = find_scan (argv[2]);
  }

  if (sc1 == NULL || sc1->meshes[mesh_level] == NULL) {
      interp->result = "target mesh is invalid";
      return (TCL_ERROR);
  }
  if (sc2 == NULL|| sc2->meshes[mesh_level] == NULL ) {
      interp->result = "merging mesh is invalid";
      return (TCL_ERROR);
  }

  /* An experiment!  It shouldn't relower already lowered confidences,
     but this is only a test */
/*
  find_mesh_edges(sc1->meshes[mesh_level]);
  lower_edge_confidence(sc1->meshes[mesh_level], mesh_level);
  find_mesh_edges(sc2->meshes[mesh_level]);
  lower_edge_confidence(sc2->meshes[mesh_level], mesh_level);
*/

  eat_edge_pair (sc1, sc2);

  return (TCL_OK);
}


/******************************************************************************
Repeatedly eat away at one mesh where it overlaps with another

Entry:
  argv[1] - mesh to eat away at
  argv[2] - other mesh, which remains uneaten
******************************************************************************/

int tcl_eat_one(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  int num;
  Scan *sc1,*sc2;

  if (argc != 3) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }

  /* find the meshes */
  sc1 = find_scan (argv[1]);
  sc2 = find_scan (argv[2]);

  if (sc1 == NULL || sc2 == NULL) {
    interp->result = "can't find mesh name";
    return (TCL_ERROR);
  }

  /* eat away at one mesh */

  init_eating (sc1);
  num = 1;
  while (num)
    num = eat_mesh_edges (sc2, sc1, 1, 0, 0, NEAR_DIST);
  done_eating (sc1);

  return (TCL_OK);
}


/******************************************************************************
Perform one eating step upon one mesh's edge where it overlaps with another.

Entry:
  argv[1] - mesh to eat edges of once
  argv[2] - other mesh, which remains uneaten

Exit:
  returns (in interp->result) the number of triangles eaten, so zero value
  means there is nothing more to eat away
******************************************************************************/

int tcl_eat_step(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  int num_eaten;
  Scan *sc1,*sc2;

  if (argc != 3) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }

  /* find the meshes */
  sc1 = find_scan (argv[1]);
  sc2 = find_scan (argv[2]);

  if (sc1 == NULL || sc2 == NULL) {
    interp->result = "can't find mesh name";
    return (TCL_ERROR);
  }

  /* eat away at one mesh */
  init_eating (sc1);
  num_eaten = eat_mesh_edges (sc2, sc1, 1, 0, 0, NEAR_DIST);
  done_eating (sc1);

  /* return number of triangles eaten */
  sprintf (interp->result, "%d", num_eaten);

  return (TCL_OK);
}



/******************************************************************************
Set the "Target" for zippering.  This then becomes the default scan to use.
******************************************************************************/

int tcl_target(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  Scan *sc;
  char str[PATH_MAX];

  if (argc != 2) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }

  sc = find_scan (argv[1]);
  if (sc == NULL) {
    fprintf (stderr, "Can't find scan called %s\n", argv[1]);
    return (TCL_ERROR);
  }

  target_scan = sc;

  sprintf (str, "select_target %d", target_scan->button_index);
  Tcl_Eval (interp, str, 0, (char **) NULL);

  return (TCL_OK);
}



/******************************************************************************
Set the "Next merge" for zippering.  This then becomes the default scan to use.
******************************************************************************/

int tcl_next_merge(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  Scan *sc;
  char str[PATH_MAX];

  if (argc != 2) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }

  sc = find_scan (argv[1]);
  if (sc == NULL) {
    fprintf (stderr, "Can't find scan called %s\n", argv[1]);
    return (TCL_ERROR);
  }

  next_merge_scan = sc;

  sprintf (str, "select_merge %d", next_merge_scan->button_index);
  Tcl_Eval (interp, str, 0, (char **) NULL);

  return (TCL_OK);
}

