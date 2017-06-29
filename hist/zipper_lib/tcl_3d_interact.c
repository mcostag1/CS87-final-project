#include <cyfile.h>
#include <tcl.h>
#include <zipper.h>
#include <gl/gl.h>
#include <gl/device.h>
#include <view.h>
#include <matrix.h>
#include <stdlib.h>
#include <stdio.h>

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




/******************************************************************************
Set the speed of translation and rotation based on mouse motion.

Entry:
  argv[1] - speed of motion (1 = default speed)
******************************************************************************/

int tcl_mouse_speed(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  if (argc != 2) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }

  mouse_speed = atof (argv[1]);
  printf ("mouse speed = %g\n", mouse_speed);

  return (TCL_OK);
}


/******************************************************************************
Turn on/off backface culling of line-drawn triangles.
******************************************************************************/

int tcl_cull(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  if (argc != 2) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }

  Tcl_GetBoolean (interp, argv[1], &back_cull);
  printf ("backface culling = %d\n", back_cull);
  draw_object();

  return (TCL_OK);
}



/******************************************************************************
Undo transformation (motion) of an object.
******************************************************************************/

int tcl_undo(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  perform_undo();
  return (TCL_OK);
}


/******************************************************************************
Redo (undo the undo) of a transformation (motion) of an object.
******************************************************************************/

int tcl_redo(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  perform_redo();
  return (TCL_OK);
}




/******************************************************************************
Specify whether to draw a given mesh.

Entry:
  argv[1] - name of mesh
  argv[2] - boolean saying whether to draw mesh
******************************************************************************/

int tcl_mesh_draw(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  char *str;
  int resolution;
  Scan *sc;

  if (argc != 3) {
    return (TCL_ERROR);
  }

  sc = find_scan (argv[1]);
  if (sc == NULL) {
    fprintf (stderr, "Can't find scan called %s\n", argv[1]);
    return (TCL_ERROR);
  }

  /* set drawing flag to the appropriate value */
  Tcl_GetBoolean (interp, argv[2], &sc->draw_flag);

  /* re-draw things */
  draw_object();

  return (TCL_OK);
}


/******************************************************************************
Select which object is being moved by user.

Entry:
  argv[1] - name of scan to move, "world" for moving everything, "light"
	    for moving the light
******************************************************************************/

int tcl_move_object(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  int i;
  int index;
  char *str;
  int resolution;
  Scan *sc;

  if (argc != 2) {
    return (TCL_ERROR);
  }

  if (strcmp ("world", argv[1]) == 0) {
    move_num = -1;
    set_move_cursor();
    return (TCL_OK);
  }
  else if (strcmp ("light", argv[1]) == 0) {
    move_num = -2;
    set_move_cursor();
    return (TCL_OK);
  }

  /* look for the scan to move */
  index = -1;
  for (i = 0; i < nscans; i++)
    if (strcmp (scans[i]->name, argv[1]) == 0) {
      index = i;
      break;
    }

  if (index == -1) {
    fprintf (stderr, "tcl_move_object: Can't find scan called %s\n", argv[1]);
    return (TCL_ERROR);
  }

  move_num = index;
  set_move_cursor();

  return (TCL_OK);
}


/******************************************************************************
Draw an anti-aliased version of the object.
******************************************************************************/

int tcl_anti_alias(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{

  draw_anti_aliased();

  return (TCL_OK);
}



/******************************************************************************
Spin the object around.

Entry:
  argv[1] - angle to spin (optional)
  argv[2] - number of steps in spin (optional)
  argv[3] - prefix of filenames to write frames to (optional)
******************************************************************************/

int tcl_spin(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  float theta;
  int nsteps;

  if (argc != 1 && argc != 3 && argc != 4) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }

  if (argc == 1) {
    spin_object (360.0, 90, NULL, 1);
  }
  else {
    theta = atof (argv[1]);
    nsteps = atoi (argv[2]);
    if (argc == 3)
      spin_object (theta, nsteps, NULL, 1);
    else
      spin_object (theta, nsteps, argv[3], 3);
  }

  return (TCL_OK);
}


/******************************************************************************
Whether or not to draw a set of axes.
******************************************************************************/

int tcl_axes(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  int val;

  if (argc != 2) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }

  /* set drawing flag to the appropriate value */
  Tcl_GetBoolean (interp, argv[1], &val);

  set_axes (val);
  draw_object();

  return (TCL_OK);
}


/******************************************************************************
Set the background color.
******************************************************************************/

int tcl_set_background(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  int val;

  if (argc != 4) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }

  back_red = atoi (argv[1]);
  back_grn = atoi (argv[2]);
  back_blu = atoi (argv[3]);

  draw_object();

  return (TCL_OK);
}


/******************************************************************************
Reset the viewing parameters.
******************************************************************************/

int tcl_reset_view(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  reset_view();
  draw_object();

  return (TCL_OK);
}



/******************************************************************************
Set flag saying whether to draw the objects during various operations.
******************************************************************************/

int tcl_set_fast_ops(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  int val;

  if (argc != 2) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }

  /* set flag */
  Tcl_GetBoolean (interp, argv[1], &val);
  printf ("fast_ops = %d\n", val);
  draw_during_ops = 1 - val;

  return (TCL_OK);
}


/******************************************************************************
Set the maximum number of parallel processes.
******************************************************************************/

int tcl_parallel_max(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  if (argc != 2) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }

  /* set maximum value */
  parallel_procs_max = atoi (argv[1]);
  printf ("maximum parallel processes = %d\n", parallel_procs_max);

  return (TCL_OK);
}


/******************************************************************************
Set the polygon surface material.

argv[1] - ambient coefficient
argv[2] - diffuse coefficient
argv[3] - specular coefficient
argv[4] - specular power
******************************************************************************/

int tcl_material(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  float amb,diff,spec,spec_pow;

  if (argc != 5) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }

  amb = atof (argv[1]);
  diff = atof (argv[2]);
  spec = atof (argv[3]);
  spec_pow = atof (argv[4]);
  printf ("amb diff spec spec_pow: %f %f %f %f\n",
          amb, diff, spec, spec_pow);
  set_polygon_material (amb, diff, spec, spec_pow);

  draw_object();

  return (TCL_OK);
}



tcl_set_singlebuffer(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  int val;
  extern set_singlebuffer(int);

  if (argc != 2) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }

  /* set flag */
  Tcl_GetBoolean (interp, argv[1], &val);
  printf ("singlebuffer = %d\n", val);

  set_singlebuffer(val);

  return (TCL_OK);
}

tcl_set_orthographic(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  int val;
  extern set_orthographic(int);

  if (argc != 2) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }

  /* set flag */
  Tcl_GetBoolean (interp, argv[1], &val);
  printf ("orthographic = %d\n", val);

  set_orthographic(val);

  return (TCL_OK);
}

