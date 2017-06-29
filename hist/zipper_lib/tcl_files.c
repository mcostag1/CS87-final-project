#include <cyfile.h>
#include <tcl.h>
#include <zipper.h>
#include <limits.h>
#include <gl/gl.h>
#include <gl/device.h>
#include <view.h>
#include <matrix.h>
#include "raw.h"
#include "ply.h"

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


/******************************************************************************
Write out polygon file.

Entry:
  argv[1] - (optional) scan to write out
  argv[2] - (optional) filename to write to
******************************************************************************/

int tcl_fileout(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  Scan *scan;
  char filename[80];

  if (argc > 3) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }

  /* default filename is out.poly */
  strcpy (filename, "out");

  if (argc == 3) {  /* write out a named scan to a named file */
    scan = find_scan (argv[1]);
    if (scan == NULL) {
      interp->result = "can't find mesh name";
      return (TCL_ERROR);
    }
    strcpy (filename, argv[2]);
    scan_to_file (scan, filename);
  }
  else if (argc == 2) {  /* write out a named scan */
    scan = find_scan (argv[1]);
    if (scan == NULL) {
      interp->result = "can't find mesh name";
      return (TCL_ERROR);
    }
    scan_to_file (scan, filename);
  }
  else {	/* write out first scan, by default */
    scan_to_file (scans[0], filename);
  }

  return (TCL_OK);
}


/******************************************************************************
Write out polygon file in binary format.

Entry:
  argv[1] - (optional) scan to write out
  argv[2] - (optional) filename to write to
******************************************************************************/

int tcl_bin_poly_write(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  Scan *scan;
  char filename[80];

  if (argc > 3) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }

  /* default filename is out.poly */
  strcpy (filename, "out");

  if (argc == 3) {  /* write out a named scan to a named file */
    scan = find_scan (argv[1]);
    if (scan == NULL) {
      interp->result = "can't find mesh name";
      return (TCL_ERROR);
    }
    strcpy (filename, argv[2]);
    write_bin_polyfile (scan, filename);
  }
  else if (argc == 2) {  /* write out a named scan */
    scan = find_scan (argv[1]);
    if (scan == NULL) {
      interp->result = "can't find mesh name";
      return (TCL_ERROR);
    }
    write_bin_polyfile (scan, filename);
  }
  else {	/* write out first scan, by default */
    write_bin_polyfile (scans[0], filename);
  }

  return (TCL_OK);
}


/******************************************************************************
Write out polygon file in PLY format.

Entry:
  argv[1] - scan to write out
  argv[2] - filename to write to
******************************************************************************/

int tcl_ply_write(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  Scan *scan;
  char filename[80];
  int writeInfo;

  if (argc > 4) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }

  writeInfo = FALSE;
  if (argc == 4)
      writeInfo = atoi(argv[3]);

  scan = find_scan (argv[1]);
  if (scan == NULL) {
    interp->result = "can't find mesh name";
    return (TCL_ERROR);
  }
  strcpy (filename, argv[2]);
  write_ply (scan, filename, writeInfo);

  return (TCL_OK);
}


/******************************************************************************
Read in a PLY polygon file.

Entry:
  argv[1]       - name of new scan to read
  argv[2,3,4]   - (optional) translation of object
  argv[5] or
  argv[5,6,7,8] - (optional) rotation in degrees or quaternions
******************************************************************************/

int tcl_ply_read(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  Scan *sc;
  int result;

  if (argc != 2 && argc != 3 && argc != 6 && argc != 9) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }

  result = read_ply (argv[1]);
  if (result) {
    interp->result = "couldn't read file";
    return (TCL_ERROR);
  }
  sc = scans[nscans-1];

  /* parse the transformation info and make the scan's buttons */
  make_scan_help (sc, argc, argv);

  /* draw all objects */
  draw_object();

  interp->result = "";

  return (TCL_OK);
}



/******************************************************************************
Read in a new scan and create the appropriate buttons for it.

Entry:
  argv[1]       - name of new scan to read
  argv[2,3,4]   - (optional) translation of object
  argv[5] or
  argv[5,6,7,8] - (optional) rotation in degrees or quaternions
******************************************************************************/

int tcl_new_scan(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  int i,j;
  Scan *sc;
  int result;
  GSPEC *read_gs_file();
  Quaternion quat;
  double xt,yt,zt;
  double rot;
  double q0,q1,q2,q3;
  char str[PATH_MAX];

  if (argc != 2 && argc != 3 && argc != 6 && argc != 9) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }

  /* allocate new scan */
  scans[nscans] = (Scan *) myalloc (sizeof (Scan));
  sc = scans[nscans];
  nscans++;

  /* read in the geometry */

  strcpy (sc->name, argv[1]);
  sc->gs = read_gs_file (sc->name);
  if (sc->gs == NULL) {
    nscans--;
    interp->result = "can't read file";
    return (TCL_ERROR);
  }

  sc->file_type = CYFILE;

  /* parse the transformation info and make the scan's buttons */
  make_scan_help (sc, argc, argv);

  if (anchor_scan == NULL) {
      sprintf (str, "anchor %s", sc->name);
      Tcl_Eval (interp, str, 0, (char **) NULL);
  } 
  else if (next_align_scan == NULL) {
      sprintf (str, "next_align %s", sc->name);
      Tcl_Eval (interp, str, 0, (char **) NULL);
  }
  
  if (target_scan == NULL) {
      sprintf (str, "target %s", sc->name);
      Tcl_Eval (interp, str, 0, (char **) NULL);
  } 
  else if (next_merge_scan == NULL) {
      sprintf (str, "next_merge %s", sc->name);
      Tcl_Eval (interp, str, 0, (char **) NULL);
  }

  /* pre-compute values for quick 3-space position computation */
  setup_geometry_info (sc);

  /* zero out mesh entries */
  for (j = 0; j < MAX_MESH_LEVELS; j++)
    sc->meshes[j] = NULL;

  /* create meshes from scan info */
  for (i = 3; i >= mesh_level; i--)
    create_scan_mesh (sc, i);

  /* draw all objects */
  draw_object();

  return (TCL_OK);
}



int tcl_read_ply_mesh(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  Scan *sc;
  int result;
  char str[PATH_MAX];
  
  if (argc < 2) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }
  if (is_range_grid_file(argv[1])) {
      result = tcl_ply_range_mesh(dummy,interp,argc,argv);
  }
  else {
      result = tcl_ply_read(dummy,interp,argc,argv);
  }

  if (result == TCL_OK) {
      sc = scans[nscans-1];
      if (anchor_scan == NULL) {
	  sprintf (str, "anchor %s", sc->name);
	  Tcl_Eval (interp, str, 0, (char **) NULL);
      } 
      else if (next_align_scan == NULL) {
	  sprintf (str, "next_align %s", sc->name);
	  Tcl_Eval (interp, str, 0, (char **) NULL);
      }

      if (target_scan == NULL) {
	  sprintf (str, "target %s", sc->name);
	  Tcl_Eval (interp, str, 0, (char **) NULL);
      } 
      else if (next_merge_scan == NULL) {
	  sprintf (str, "next_merge %s", sc->name);
	  Tcl_Eval (interp, str, 0, (char **) NULL);
      }
  }

  return result;
}



/******************************************************************************
Read in a PLY file containing range data.

Entry:
  argv[1]       - name of new scan to read
  argv[2,3,4]   - (optional) translation of object
  argv[5] or
  argv[5,6,7,8] - (optional) rotation in degrees or quaternions
******************************************************************************/

int tcl_ply_range_mesh(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  int i,j;
  Scan *sc;
  int result;
  Quaternion quat;
  GSPEC *read_gs_file();
  RangeData *read_ply_geom();

  if (argc != 2 && argc != 3 && argc != 6 && argc != 9) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }

  /* allocate new scan */
  scans[nscans] = (Scan *) myalloc (sizeof (Scan));
  sc = scans[nscans];
  nscans++;

  /* read in the geometry */

  strcpy (sc->name, argv[1]);
  sc->ply_geom = read_ply_geom (sc->name);
  if (sc->ply_geom == NULL) {
    nscans--;
    interp->result = "can't read file";
    return (TCL_ERROR);
  }

  sc->num_obj_info = sc->ply_geom->num_obj_info;
  for (i = 0; i < sc->ply_geom->num_obj_info; i++) {
    strcpy(sc->obj_info[i], sc->ply_geom->obj_info[i]);
  }
  

  sc->file_type = PLYRANGEFILE;

  /* parse the transformation info and make the scan's buttons */
  make_scan_help (sc, argc, argv);

  /* zero out mesh entries */
  for (j = 0; j < MAX_MESH_LEVELS; j++)
    sc->meshes[j] = NULL;

  /* create mesh from scan info */
  for (i = 3; i >= mesh_level; i--)
    create_scan_mesh (sc, i);

  /* draw all objects */
  draw_object();

  return (TCL_OK);
}


/******************************************************************************
Read in a raw data file and create the appropriate buttons for it.

Entry:
  argv[1]       - name of new scan to read
  argv[2,3,4]   - (optional) translation of object
  argv[5] or
  argv[5,6,7,8] - (optional) rotation in degrees or quaternions
******************************************************************************/

int tcl_raw_mesh(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  int i,j;
  Scan *sc;
  int result;
  Quaternion quat;
  GSPEC *read_gs_file();
  RawData *read_raw_geom();

  if (argc != 2 && argc != 3 && argc != 6 && argc != 9) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }

  /* allocate new scan */
  scans[nscans] = (Scan *) myalloc (sizeof (Scan));
  sc = scans[nscans];
  nscans++;

  /* read in the geometry */

  strcpy (sc->name, argv[1]);
  sc->raw_geom = read_raw_geom (sc->name);
  if (sc->raw_geom == NULL) {
    nscans--;
    interp->result = "can't read file";
    return (TCL_ERROR);
  }

  sc->file_type = RAWFILE;

  /* parse the transformation info and make the scan's buttons */
  make_scan_help (sc, argc, argv);

  /* zero out mesh entries */
  for (j = 0; j < MAX_MESH_LEVELS; j++)
    sc->meshes[j] = NULL;

  /* create mesh from scan info */
  for (i = 3; i >= mesh_level; i--)
    create_scan_mesh (sc, i);

  /* draw all objects */
  draw_object();

  return (TCL_OK);
}


/******************************************************************************
Read in a polygon file.
******************************************************************************/

int tcl_polygon_file(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  Scan *sc;
  int result;

  if (argc != 2 && argc != 3 && argc != 6 && argc != 9) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }

  /* read in the polygons */

  result = file_to_scan (argv[1]);
  if (result) {
    interp->result = "couldn't read file";
    return (TCL_ERROR);
  }
  sc = scans[nscans-1];

  /* parse the transformation info and make the scan's buttons */
  make_scan_help (sc, argc, argv);

  /* draw all objects */
  draw_object();

  interp->result = "";

  return (TCL_OK);
}


/******************************************************************************
Read in a binary polygon file.
******************************************************************************/

int tcl_bin_poly_read(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  Scan *sc;
  int result;

  if (argc != 2 && argc != 3 && argc != 6 && argc != 9) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }

  /* read in the polygons */

  result = read_bin_polyfile (argv[1]);
  if (result) {
    interp->result = "couldn't read file";
    return (TCL_ERROR);
  }
  sc = scans[nscans-1];

  /* parse the transformation info and make the scan's buttons */
  make_scan_help (sc, argc, argv);

  /* draw all objects */
  draw_object();

  interp->result = "";

  return (TCL_OK);
}



/******************************************************************************
Specify whether to write out the intensity at each vertex.

Entry:
  argv[1] - boolean saying whether to write out intensity
******************************************************************************/

int tcl_intensity_write(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{

  if (argc != 2) {
    return (TCL_ERROR);
  }

  /* set the flag */
  Tcl_GetBoolean (interp, argv[1], &write_intensity_flag);

  return (TCL_OK);
}


/******************************************************************************
Specify whether to write out the colors of each vertex.

Entry:
  argv[1] - boolean saying whether to write out colors
******************************************************************************/

int tcl_colors_write(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{

  if (argc != 2) {
    return (TCL_ERROR);
  }

  /* set the flag */
  Tcl_GetBoolean (interp, argv[1], &write_colors_flag);

  return (TCL_OK);
}


/******************************************************************************
Specify whether to write out the surface normals at each vertex.

Entry:
  argv[1] - boolean saying whether to write out intensity
******************************************************************************/

int tcl_normals_write(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{

  if (argc != 2) {
    return (TCL_ERROR);
  }

  /* set the flag */
  Tcl_GetBoolean (interp, argv[1], &write_normals_flag);

  return (TCL_OK);
}

