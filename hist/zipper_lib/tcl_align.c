#include <cyfile.h>
#include <tcl.h>
#include <zipper.h>
#include <limits.h>
#include <gl/gl.h>
#include <gl/device.h>
#include <view.h>
#include <matrix.h>


extern Tcl_Interp *interp;

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


Scan *anchor_scan = NULL; /* default target for alignment */
Scan *next_align_scan = NULL; /* default source for alignment */


/******************************************************************************
Trigger the matching algorithm.
******************************************************************************/

int tcl_match(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  match_proc();
  return (TCL_OK);
}


/******************************************************************************
Invoke Besl's matching algorithm.
******************************************************************************/

int tcl_close(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  closest_proc();
  return (TCL_OK);
}



/******************************************************************************
Match a list of meshes to another list, possibly dragging along other meshes.

Entry:
  argv[1] - meshes to move
  argv[2] - meshes to match to
  argv[3] - meshes to drag along
******************************************************************************/

int tcl_match_list(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  int i;
  char **names;
  Scan **move_list;
  Scan **match_list;
  Scan **drag_list;
  int nmove,nmatch,ndrag;
  
  if (argc < 3 || argc > 4) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }

  /* get list of meshes to move */
  if (Tcl_SplitList (interp, argv[1], &nmove, &names) != TCL_OK)
    return (TCL_ERROR);
  move_list = make_scan_list (nmove, names);
  free (names);

  /* get list of meshes to match */
  if (Tcl_SplitList (interp, argv[2], &nmatch, &names) != TCL_OK)
    return (TCL_ERROR);
  match_list = make_scan_list (nmatch, names);
  free (names);

  /* error check */
  if (move_list == NULL || match_list == NULL) {
    interp->result = "couldn't find mesh name";
    return (TCL_ERROR);
  }

  /* maybe get list of meshes to drag along with the ones to move */
  if (argc == 4) {
    if (Tcl_SplitList (interp, argv[1], &ndrag, &names) != TCL_OK)
      return (TCL_ERROR);
    drag_list = make_scan_list (ndrag, names);
    free (names);
  }

  /* call the match routine */
  if (argc == 3)
    new_iterative_match(move_list, match_list, NULL, nmove, nmatch, 0);
  else
    new_iterative_match(move_list, match_list, drag_list, nmove, nmatch, ndrag);

  /* free the lists */
  free (move_list);
  free (match_list);
  if (argc == 4)
    free (drag_list);

  return (TCL_OK);
}


/******************************************************************************
******************************************************************************/

int tcl_align_error(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  int i;
  Mesh *mesh;
  Scan *sc1,*sc2;
  char str[80];
  float error;
  Matrix ident;

  if (argc != 1 && argc != 2 && argc != 3) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }

  if (argc == 1) {
    sc1 = anchor_scan;
    sc2 = next_align_scan;
  }
  else if (argc == 2) {
    sc1 = anchor_scan;
    sc2 = find_scan (argv[1]);
  }
  else {
    sc1 = find_scan (argv[1]);
    sc2 = find_scan (argv[2]);
  }

  if (sc1 == NULL || sc1->meshes[mesh_level] == NULL) {
      interp->result = "anchor mesh is invalid";
      return (TCL_ERROR);
  }
  if (sc2 == NULL || sc2->meshes[mesh_level] == NULL) {
      interp->result = "aligning mesh is invalid";
      return (TCL_ERROR);
  }

  create_match_list (sc1, sc2, 0, 0);
  mat_ident(ident);
  error = measure_error(ident,sc2);
  sprintf(str, "%g mm", error*1000);

  Tcl_SetResult(interp, str, TCL_VOLATILE);

  return (TCL_OK);
}


/******************************************************************************
Match a list of meshes to another list, possibly dragging along other meshes.
This verison uses Besl's matching algorithm.

Entry:
  argv[1] - meshes to move
  argv[2] - meshes to match to
  argv[3] - meshes to drag along
******************************************************************************/

int tcl_align(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  int i;
  char str[80];
  char **names;
  Scan **move_list;
  Scan **match_list;
  Scan **drag_list;
  int nmove,nmatch,ndrag;
  
  if (argc > 4) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }


  /* get list of meshes to move */
  if (argc > 1) {
      if (Tcl_SplitList (interp, argv[1], &nmove, &names) != TCL_OK)
	  return (TCL_ERROR);
  } 
  else {
      if (next_align_scan == NULL) {
	  interp->result = "couldn't find mesh to align";
	  return (TCL_ERROR);
      }

      if (next_align_scan->meshes[mesh_level] == NULL) {
	  interp->result = "couldn't find mesh to align";
	  return (TCL_ERROR);
      }

      sprintf (str, "%s", next_align_scan->name);
      if (Tcl_SplitList (interp, str, &nmove, &names) != TCL_OK)
	  return (TCL_ERROR);
  }
  move_list = make_scan_list (nmove, names);
  free (names);

  if (argc <= 2) {
    /* create list containing just the anchor mesh */
      if (anchor_scan == NULL) {
	  interp->result = "couldn't find anchor mesh";
	  return (TCL_ERROR);
      }

      if (anchor_scan->meshes[mesh_level] == NULL) {
	  interp->result = "couldn't find anchor mesh";
	  return (TCL_ERROR);
      }

    sprintf (str, "%s", anchor_scan->name);
    if (Tcl_SplitList (interp, str, &nmatch, &names) != TCL_OK)
      return (TCL_ERROR);
  }
  else {
    /* get list of meshes to match */
    if (Tcl_SplitList (interp, argv[2], &nmatch, &names) != TCL_OK)
      return (TCL_ERROR);
  }
  match_list = make_scan_list (nmatch, names);
  free (names);

  /* error check */
  if (move_list == NULL || match_list == NULL) {
    interp->result = "couldn't find mesh name";
    return (TCL_ERROR);
  }

  /* maybe get list of meshes to drag along with the ones to move */
  if (argc == 4) {
    if (Tcl_SplitList (interp, argv[1], &ndrag, &names) != TCL_OK)
      return (TCL_ERROR);
    drag_list = make_scan_list (ndrag, names);
    free (names);
  }

  /* call the match routine */
  /* (note the reverse ordering in the call!!) */

  Tcl_Eval (interp, "disable_all $button_list", 0, (char **) NULL);
  check_tk();
  six_degree_match (&move_list[0], &match_list[0], NULL, 1, 1, 0, 10);
  Tcl_Eval (interp, "enable_all $button_list", 0, (char **) NULL);

  /* free the lists */
  free (move_list);
  free (match_list);
  if (argc == 4)
    free (drag_list);

  return (TCL_OK);
}


/******************************************************************************
Match a list of meshes to another list, possibly dragging along other meshes.
This verison uses Besl's matching algorithm.  Only performs one iteration.

Entry:
  argv[1] - meshes to move
  argv[2] - meshes to match to
  argv[3] - meshes to drag along
******************************************************************************/

int tcl_align_step(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  int i;
  char str[80];
  char **names;
  Scan **move_list;
  Scan **match_list;
  Scan **drag_list;
  int nmove,nmatch,ndrag;
  
  if (argc > 4) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }


  /* get list of meshes to move */
  if (argc > 1) {
      if (Tcl_SplitList (interp, argv[1], &nmove, &names) != TCL_OK)
	  return (TCL_ERROR);
  } 
  else {
      if (next_align_scan == NULL) {
	  interp->result = "couldn't find mesh to align";
	  return (TCL_ERROR);
      }

      if (next_align_scan->meshes[mesh_level] == NULL) {
	  interp->result = "couldn't find mesh to align";
	  return (TCL_ERROR);
      }

      sprintf (str, "%s", next_align_scan->name);
      if (Tcl_SplitList (interp, str, &nmove, &names) != TCL_OK)
	  return (TCL_ERROR);
  }
  move_list = make_scan_list (nmove, names);
  free (names);

  if (argc <= 2) {
    /* create list containing just the anchor mesh */
      if (anchor_scan == NULL) {
	  interp->result = "couldn't find anchor mesh";
	  return (TCL_ERROR);
      }

      if (anchor_scan->meshes[mesh_level] == NULL) {
	  interp->result = "couldn't find anchor mesh";
	  return (TCL_ERROR);
      }

    sprintf (str, "%s", anchor_scan->name);
    if (Tcl_SplitList (interp, str, &nmatch, &names) != TCL_OK)
      return (TCL_ERROR);
  }
  else {
    /* get list of meshes to match */
    if (Tcl_SplitList (interp, argv[2], &nmatch, &names) != TCL_OK)
      return (TCL_ERROR);
  }
  match_list = make_scan_list (nmatch, names);
  free (names);

  /* error check */
  if (move_list == NULL || match_list == NULL) {
    interp->result = "couldn't find mesh name";
    return (TCL_ERROR);
  }

  /* maybe get list of meshes to drag along with the ones to move */
  if (argc == 4) {
    if (Tcl_SplitList (interp, argv[1], &ndrag, &names) != TCL_OK)
      return (TCL_ERROR);
    drag_list = make_scan_list (ndrag, names);
    free (names);
  }

  /* call the match routine */
  /* (note the reverse ordering in the call!!) */

  Tcl_Eval (interp, "disable_all $button_list", 0, (char **) NULL);
  check_tk();
  /* The only difference between this routine and the previous is */
  /* the last argument */
  six_degree_match (&move_list[0], &match_list[0], NULL, 1, 1, 0, 1);
  Tcl_Eval (interp, "enable_all $button_list", 0, (char **) NULL);

  /* free the lists */
  free (move_list);
  free (match_list);
  if (argc == 4)
    free (drag_list);

  return (TCL_OK);
}


/******************************************************************************
Match a list of meshes to another list, possibly dragging along other meshes.
This verison uses Besl's matching algorithm.

Entry:
  argv[1] - meshes to move
  argv[2] - meshes to match to
  argv[3] - meshes to drag along
******************************************************************************/

int tcl_besl_old(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  fprintf (stderr, "The 'besl' command is now called 'align'.\n");
  return (TCL_OK);
}


/******************************************************************************
Set the "Anchor" for alignment.  This then becomes the default scan to use.
******************************************************************************/

int tcl_anchor(dummy,interp,argc,argv)
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

  anchor_scan = sc;

  sprintf (str, "select_anchor %d", anchor_scan->button_index);
  Tcl_Eval (interp, str, 0, (char **) NULL);

  return (TCL_OK);
}


/******************************************************************************
Set the "Next align" for alignment.  This then becomes the default scan to use.
******************************************************************************/

int tcl_next_align(dummy,interp,argc,argv)
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

  next_align_scan = sc;

  sprintf (str, "select_align %d", next_align_scan->button_index);
  Tcl_Eval (interp, str, 0, (char **) NULL);

  return (TCL_OK);
}


/******************************************************************************
Stop the alignment process (assuming it is going on).
******************************************************************************/

int tcl_stop_align(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  align_on = 0;

  return (TCL_OK);
}

