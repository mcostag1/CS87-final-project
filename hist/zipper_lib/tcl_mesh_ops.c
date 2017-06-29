#include <cyfile.h>
#include <tcl.h>
#include <zipper.h>
#include <gl/gl.h>
#include <gl/device.h>
#include <view.h>
#include <matrix.h>
#include <stdlib.h>

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

extern Scan *target_scan; /* default target for zippering */
extern Scan *next_merge_scan; /* default source for zippering */
extern Scan *anchor_scan; /* default target for alignment */
extern Scan *next_align_scan; /* default source for alignment */


static int button_count = 0;     /* button indices for scans */


/******************************************************************************
Remove any short edges between vertices of two recently zippered meshes.

Entry:
  argv[1] - name of scan to remove short edges from
  argv[2] - (optional) fraction of spacing that will be min allowed distance
******************************************************************************/

int tcl_edge_remove(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  Scan *scan;
  float fract;

  if (argc != 2 && argc != 3) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }

  scan = find_scan (argv[1]);
  if (scan == NULL) {
    interp->result = "can't find mesh name";
    return (TCL_ERROR);
  }

  if (argc == 2)
    remove_short_edges (scan, 0.333);
  else {
    fract = atof (argv[2]);
    remove_short_edges (scan, fract);
  }

  draw_object();

  return (TCL_OK);
}


/******************************************************************************
Remove "sliver" triangles from a mesh.

Entry:
  argv[1] - name of scan to remove slivers from
  argv[2] - (optional) fraction of spacing that will be min allowed distance
******************************************************************************/

int tcl_sliver_remove(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  Scan *scan;
  float fract;

  if (argc != 2 && argc != 3) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }

  scan = find_scan (argv[1]);
  if (scan == NULL) {
    interp->result = "can't find mesh name";
    return (TCL_ERROR);
  }

  if (argc == 2)
    remove_sliver_tris (scan, 0.333);
  else {
    fract = atof (argv[2]);
    remove_sliver_tris (scan, fract);
  }

  draw_object();

  return (TCL_OK);
}

int tcl_bad_aspect_remove(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  Scan *scan;
  float fract;
  float max_aspect, min_cos;
  int diff;

  if (argc != 5) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }

  scan = find_scan (argv[1]);
  if (scan == NULL) {
    interp->result = "can't find mesh name";
    return (TCL_ERROR);
  }

  max_aspect = atof(argv[2]);
  min_cos = atof(argv[3]);
  diff = atoi(argv[4]);
/*
  printf("Max aspect = %f.\nMin cos = %f\n", max_aspect, min_cos);
*/
  remove_bad_aspect_tris(scan, max_aspect, min_cos, diff);

  draw_object();

  return (TCL_OK);
}


int tcl_flat_verts_remove(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  Scan *scan;
  float max_cos;

  if (argc != 3) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }

  scan = find_scan (argv[1]);
  if (scan == NULL) {
    interp->result = "can't find mesh name";
    return (TCL_ERROR);
  }

  max_cos = atof(argv[2]);

  remove_flat_verts(scan, max_cos);

  draw_object();

  return (TCL_OK);
}


/******************************************************************************
Remove the "extra" vertices that were introduced during zippering.
******************************************************************************/

int tcl_vertex_remove(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  Scan *scan;

  if (argc != 2) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }

  scan = find_scan (argv[1]);
  if (scan == NULL) {
    interp->result = "can't find mesh name";
    return (TCL_ERROR);
  }

  remove_cut_vertices (scan);
  draw_object();

  return (TCL_OK);
}




/******************************************************************************
Fix vertices that have too many unshared edges.

Entry:
  argv[1] - name of mesh to fix
******************************************************************************/

int tcl_fix_bows(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  Scan *sc;

  if (argc != 2) {
    return (TCL_ERROR);
  }

  sc = find_scan (argv[1]);
  if (sc == NULL) {
    fprintf (stderr, "Can't find scan called %s\n", argv[1]);
    return (TCL_ERROR);
  }

  fix_bows (sc);
  draw_object();

  return (TCL_OK);
}


/******************************************************************************
Clean up a mesh by doing all sorts of consistancy checks and removing
offending triangles and vertices.

Entry:
  argv[1] - name of mesh to clean up
******************************************************************************/

int tcl_clean_mesh(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  Scan *sc;

  if (argc != 2) {
    return (TCL_ERROR);
  }

  sc = find_scan (argv[1]);
  if (sc == NULL) {
    fprintf (stderr, "Can't find scan called %s\n", argv[1]);
    return (TCL_ERROR);
  }

  clean_up_mesh (sc);
  draw_object();

  return (TCL_OK);
}


/******************************************************************************
Fill in small holes in a mesh

Entry:
  argv[1] - name of mesh to fix
******************************************************************************/

int tcl_fill_holes(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  Scan *sc;

  if (argc != 2) {
    return (TCL_ERROR);
  }

  sc = find_scan (argv[1]);
  if (sc == NULL) {
    fprintf (stderr, "Can't find scan called %s\n", argv[1]);
    return (TCL_ERROR);
  }

  fill_small_holes (sc);
  draw_object();

  return (TCL_OK);
}


/******************************************************************************
Test out the triangle split routine.

Entry:
  argv[1] - name of scan to split triangles in
******************************************************************************/

int tcl_split_test(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  Scan *sc;

  if (argc != 2) {
    return (TCL_ERROR);
  }

  sc = find_scan (argv[1]);
  if (sc == NULL) {
    fprintf (stderr, "Can't find scan called %s\n", argv[1]);
    return (TCL_ERROR);
  }

  split_test (sc);
  draw_object();

  return (TCL_OK);
}


/******************************************************************************
Select the resolution of the meshes being displayed.
******************************************************************************/

int tcl_mesh_resolution(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  char *str;
  char str2[80];
  int resolution;

  if (argc != 1 && argc != 2) {
    return (TCL_ERROR);
  }

  if (argc == 1) {
    str = Tcl_GetVar (interp, "mesh_resolution", 0);
    resolution = atoi (str);
    if ((resolution < 0) || (resolution > 3))
      return (TCL_ERROR);
  }
  else {
    resolution = atoi (argv[1]);
    resolution -= 1;
    if ((resolution < 0) || (resolution > 3)) {
      fprintf (stderr, "bad mesh resolution\n");
      return (TCL_ERROR);
    }
    sprintf (str2, "set mesh_resolution %d", resolution);
    Tcl_Eval (interp, str2, 0, (char **) NULL);
  }

  mesh_level = resolution;
  create_current_level();
  draw_object();

  return (TCL_OK);
}





/******************************************************************************
Cut each of the triangles in a mesh into four smaller triangles.
******************************************************************************/

int tcl_quarter_mesh(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  Scan *sc;

  if (argc != 2) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }

  sc = find_scan (argv[1]);
  if (sc == NULL) {
    fprintf (stderr, "Can't find scan called %s\n", argv[1]);
    return (TCL_ERROR);
  }

  quarter_mesh (sc);
  draw_object();

  return (TCL_OK);
}



/******************************************************************************
Looks for pairs of polygons that share an edge and see if the surface
should be changed to use the other diagonal.
******************************************************************************/

int tcl_swap_edges(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  Scan *sc;

  if (argc != 2) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }

  sc = find_scan (argv[1]);
  if (sc == NULL) {
    fprintf (stderr, "Can't find scan called %s\n", argv[1]);
    return (TCL_ERROR);
  }

  swap_edges (sc);
  draw_object();

  return (TCL_OK);
}


/******************************************************************************
Smooth the vertex positions by taking local averages.
******************************************************************************/

int tcl_smooth_vertices(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  Scan *sc;

  if (argc != 2) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }

  sc = find_scan (argv[1]);
  if (sc == NULL) {
    fprintf (stderr, "Can't find scan called %s\n", argv[1]);
    return (TCL_ERROR);
  }

  smooth_vertices (sc);
  draw_object();

  return (TCL_OK);
}




/******************************************************************************
Remove a mesh from the screen, from memory and from the list of buttons.
******************************************************************************/

int tcl_remove_mesh(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  int i;
  char str[80];
  int index;
  Scan *scan;

  if (argc != 2) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }

  /* find the mesh */

  index = -1;
  for (i = 0; i < nscans; i++)
    if (strcmp (argv[1], scans[i]->name) == 0) {
      index = i;
      break;
    }
  if (index == -1) {
    interp->result = "can't find mesh name";
    return (TCL_ERROR);
  }


  /* Clear some pointers that may reference this scan */

  if (anchor_scan == scan)
      anchor_scan = NULL;
  if (next_align_scan == scan)
      next_align_scan = NULL;
  if (target_scan == scan)
      target_scan = NULL;
  if (next_merge_scan == scan)
      next_merge_scan = NULL;




  /* remove the scan from the list of scans */

  scan = scans[index];
  for (i = index; i < nscans; i++)
    scans[i] = scans[i+1];
  nscans--;

  /* free up memory from meshes in scan */

  for (i = 0; i < MAX_MESH_LEVELS; i++)
    if (scan->meshes[i] != NULL)
      clear_mesh (scan->meshes[i]);

  /* free up memory of raw data (if any) */

  switch (scan->file_type) {
    case CYFILE:
      free (scan->gs);
      break;
    case RAWFILE:
      free (scan->raw_geom);
      break;
    case PLYRANGEFILE:
      delete_ply_geom(scan->ply_geom);
      break;
    case POLYFILE:
      break;
    default:
      printf ("unknown scan type in tcl_remove_mesh: %d\n", scan->file_type);
      break;
  }

  /* remove the button */

  sprintf (str, "destroy .meshbuttons.obj_list.frame%d", scan->button_index);
  Tcl_Eval (interp, str, 0, (char **) NULL);

  /* free up scan itself */

  free (scan);

  /* re-display everything */

  draw_object();

  return (TCL_OK);
}


/******************************************************************************
Rename a mesh.
******************************************************************************/

int tcl_rename_mesh(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  int i;
  char str[80];
  int index;
  Scan *scan;

  if (argc != 3) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }

  scan = rename_scan(argv[1], argv[2]);
  if (scan != NULL) {
      sprintf (str, "rename_mesh_buttons .meshbuttons %s %d",
	       scan->name, scan->button_index);
      Tcl_Eval (interp, str, 0, (char **) NULL);
  } else {
    interp->result = "can't rename";
    return (TCL_ERROR);
  }
}



/******************************************************************************
Help create a scan by parsing the transformation information and making
the buttons for the scan.

Entry:
  sc - scan into which to place transformation
  argc - number of arguments in command line
  argv - arguments in command line
******************************************************************************/

make_scan_help(sc,argc,argv)
  Scan *sc;
  int argc;
  char *argv[];
{
  char str[200];
  Quaternion quat;
  double xt,yt,zt;
  double rot;
  double q0,q1,q2,q3;

  sc->draw_flag = 1;
  sc->edge_mesh = NULL;

  /* parse the transformation info */

  /* argc == 2, no transformation */
  /* argc == 3, <rot-about-y-axis> */
  /* argc == 6, <tx> <ty> <tz> <rot-about-y-axis> */
  /* argc == 9, <tx> <ty> <tz> <q0> <q1> <q2> <q3> */

  /* translation */

  if (argc == 2 || argc == 3) {
    sc->xtrans = 0;
    sc->ytrans = 0;
    sc->ztrans = 0;
  }
  else {
    Tcl_GetDouble (interp, argv[2], &xt);
    Tcl_GetDouble (interp, argv[3], &yt);
    Tcl_GetDouble (interp, argv[4], &zt);
    sc->xtrans = xt;
    sc->ytrans = yt;
    sc->ztrans = zt;
  }

  /* rotation */

  if (argc == 9) {
    Tcl_GetDouble (interp, argv[5], &q0);
    Tcl_GetDouble (interp, argv[6], &q1);
    Tcl_GetDouble (interp, argv[7], &q2);
    Tcl_GetDouble (interp, argv[8], &q3);
    quat[0] = q0;
    quat[1] = q1;
    quat[2] = q2;
    quat[3] = q3;
    quat_to_mat (quat, sc->rotmat);
  }
  else if (argc == 6) {
    Tcl_GetDouble (interp, argv[5], &rot);
    sc->rotate = rot;
    build_rotmat (sc);
  }
  else if (argc == 3) {
    Tcl_GetDouble (interp, argv[2], &rot);
    sc->rotate = rot;
    build_rotmat (sc);
  }
  else {
    sc->rotate = 0;
    build_rotmat (sc);
  }

  /* pick a button index */
  sc->button_index = button_count;
  button_count++;

  /* make the buttons for this scan */
  sprintf (str, "make_mesh_buttons .meshbuttons %s %d",
           sc->name, sc->button_index);
  Tcl_Eval (interp, str, 0, (char **) NULL);
  sprintf (str, "set mesh_draw_%d 1", sc->button_index);
  Tcl_Eval (interp, str, 0, (char **) NULL);

  return (TCL_OK);
}

