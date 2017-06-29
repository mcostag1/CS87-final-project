#include <cyfile.h>
#include <tcl.h>
#include <zipper.h>
#include <stdio.h>
#include <fcntl.h>
#include <gl/gl.h>
#include <gl/device.h>
#include <view.h>
#include <matrix.h>

#define PERMS 0644

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
Set the viewing parameters.

Entry:
  argv[1]-argv[3] - translation
  argv[4]-argv[7] - rotation
******************************************************************************/

int tcl_camera(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  Quaternion quat;
  double xt,yt,zt;
  double q0,q1,q2,q3;

  if (argc != 8) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }

  /* get translation */
  Tcl_GetDouble (interp, argv[1], &xt);
  Tcl_GetDouble (interp, argv[2], &yt);
  Tcl_GetDouble (interp, argv[3], &zt);
  transmat[3][0] = xt;
  transmat[3][1] = yt;
  transmat[3][2] = zt;

  /* get rotation */
  Tcl_GetDouble (interp, argv[4], &q0);
  Tcl_GetDouble (interp, argv[5], &q1);
  Tcl_GetDouble (interp, argv[6], &q2);
  Tcl_GetDouble (interp, argv[7], &q3);
  quat[0] = q0;
  quat[1] = q1;
  quat[2] = q2;
  quat[3] = q3;
  quat_to_mat (quat, rotmat);

  /* re-draw the objects */
  draw_object();

  return (TCL_OK);
}


/******************************************************************************
Print out the positions of the various scans.
******************************************************************************/

int tcl_print_positions(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  FILE *fp;
  char filename[80];
  char filename2[80];
  int fileid;

  if (argc == 1) {
    print_positions (stdout);
  }
  else if (argc == 2) {

    /* append .conf if necessary */
    strcpy (filename, argv[1]);
    if (strlen (filename) < 5 ||
        strcmp (filename + strlen (filename) - 5, ".conf") != 0)
        strcat (filename, ".conf");

    /* see if file already exists, and rename the old file if it does */
    fileid = open (filename, O_RDONLY, PERMS);
    if (fileid != -1) {
      close (fileid);
      sprintf (filename2, "%s.bak", filename);
      fprintf (stderr, "Renaming old '%s' to '%s'\n", filename, filename2);
      rename (filename, filename2);
    }

    /* write out configuration */
    fp = fopen (filename, "w");
    if (fp == NULL) {
      fprintf (stderr, "Can't open file '%s'\n", filename);
      return (TCL_ERROR);
    }
    fprintf (stderr, "Writing to '%s'\n", filename);
    print_positions (fp);
    fclose (fp);

  }
  else {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }

  return (TCL_OK);
}


/******************************************************************************
Print out the positions of the various scans in matrix form.
******************************************************************************/

int tcl_mprint_positions(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  FILE *fp;
  char filename[80];
  char filename2[80];
  int fileid;

  if (argc == 1) {
    print_mat_positions (stdout);
  }
  else if (argc == 2) {

    /* append .conf if necessary */
    strcpy (filename, argv[1]);
    if (strlen (filename) < 4 ||
        strcmp (filename + strlen (filename) - 4, ".mat") != 0)
        strcat (filename, ".mat");

    /* see if file already exists, and rename the old file if it does */
    fileid = open (filename, O_RDONLY, PERMS);
    if (fileid != -1) {
      close (fileid);
      sprintf (filename2, "%s.bak", filename);
      fprintf (stderr, "Renaming old '%s' to '%s'\n", filename, filename2);
      rename (filename, filename2);
    }

    /* write out matrices */
    fp = fopen (filename, "w");
    if (fp == NULL) {
      fprintf (stderr, "Can't open file '%s'\n", filename);
      return (TCL_ERROR);
    }
    fprintf (stderr, "Writing to '%s'\n", filename);
    print_mat_positions (fp);
    fclose (fp);

  }
  else {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }

  return (TCL_OK);
}


/******************************************************************************
Translate a named scan.
******************************************************************************/

int tcl_translate(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  Scan *sc;

  if (argc != 5) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }

  sc = find_scan (argv[1]);
  if (sc == NULL) {
    fprintf (stderr, "Can't find scan called %s\n", argv[1]);
    return (TCL_ERROR);
  }

  sc->xtrans += atof (argv[2]);
  sc->ytrans += atof (argv[3]);
  sc->ztrans += atof (argv[4]);

  draw_object();

  return (TCL_OK);
}


/******************************************************************************
Rotate a named scan around an axis.
******************************************************************************/

int tcl_rotate(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  Scan *sc;
  float angle;
  Matrix mat;
  Quaternion quat;
  double q0, q1, q2, q3;

  if (argc != 3 && argc != 6) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }

  sc = find_scan (argv[1]);
  if (sc == NULL) {
    fprintf (stderr, "Can't find scan called %s\n", argv[1]);
    return (TCL_ERROR);
  }

  if (argc == 3) {
      angle = atof (argv[2]);
      switch (argv[0][0]) {
      case 'x':
	  mat_rotate (mat, angle, 'x');
	  break;
      case 'y':
	  mat_rotate (mat, angle, 'y');
	  break;
      case 'z':
	  mat_rotate (mat, angle, 'z');
	  break;
      }
  } 
  else {
    Tcl_GetDouble (interp, argv[2], &q0);
    Tcl_GetDouble (interp, argv[3], &q1);
    Tcl_GetDouble (interp, argv[4], &q2);
    Tcl_GetDouble (interp, argv[5], &q3);
    quat[0] = q0;
    quat[1] = q1;
    quat[2] = q2;
    quat[3] = q3;
    quat_to_mat (quat, mat);
  }

  mat_mult (sc->rotmat, mat, sc->rotmat);

  draw_object();

  return (TCL_OK);
}


/******************************************************************************
Absorb the transformation into the current meshes vertices, resetting
the object transformation to the identity when finished
******************************************************************************/

int tcl_absorb_transform(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  Scan *sc;
  float angle;
  Matrix mat;
  Quaternion quat;
  double q0, q1, q2, q3;

  if (argc != 2) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }

  sc = find_scan (argv[1]);
  if (sc == NULL) {
    fprintf (stderr, "Can't find scan called %s\n", argv[1]);
    return (TCL_ERROR);
  }

  absorb_transform(sc);

  return (TCL_OK);
}


/* tcl_tranform
 * 
 * this function composes two transforms and returns the resulting
 * one 
 *
 * Originally by Brian Curless
 * Modified:  Afra Zomorodian 7/5/95
 * Essentially, given two transforms of (translation and rotation in form
 * of a quaternion), we have:  
 * Final Rotation = R(q1)*R(q2)
 * Final Translation = R(q1)*T2 + T1
 */
tcl_transform(dummy, interp, argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char *argv[];
{
  char str[200];
  Quaternion quat1, quat2, quat3;
  Matrix rotmat1, rotmat2, rotmatFinal;
  Matrix transmat1, transmat2;
  Matrix mat1, mat2;
  Vector transVector1, transVector2, transVectorFinal;
  double xt, yt, zt, q0, q1, q2, q3;

  if (argc != 15) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }
  /* parse the transformation info */
  Tcl_GetDouble (interp, argv[1], &xt);
  Tcl_GetDouble (interp, argv[2], &yt);
  Tcl_GetDouble (interp, argv[3], &zt);
  Tcl_GetDouble (interp, argv[4], &q0);
  Tcl_GetDouble (interp, argv[5], &q1);
  Tcl_GetDouble (interp, argv[6], &q2);
  Tcl_GetDouble (interp, argv[7], &q3);
  
  /* Convert to float */
  transVector1[0] = xt;
  transVector1[1] = yt;
  transVector1[2] = zt;
  quat1[0] = q0;
  quat1[1] = q1;
  quat1[2] = q2;
  quat1[3] = q3;

  /* printf("The first transl. vector is %f, %f, %f\n",
	 transVector1[0], transVector1[1], transVector1[2]);
  printf("The first quaternion is %f, %f, %f, %f\n", quat1[0], quat1[1], 
	 quat1[2], quat1[3]);*/

  quat_to_mat (quat1, rotmat1);
  /* printf("The first rot matrix is this:\n");
  mat_print(rotmat1); */
  
  /* Parse Info */
  Tcl_GetDouble (interp, argv[8], &xt);
  Tcl_GetDouble (interp, argv[9], &yt);
  Tcl_GetDouble (interp, argv[10], &zt);
  Tcl_GetDouble (interp, argv[11], &q0);
  Tcl_GetDouble (interp, argv[12], &q1);
  Tcl_GetDouble (interp, argv[13], &q2);
  Tcl_GetDouble (interp, argv[14], &q3);

  /* Convert to float */
  transVector2[0] = xt;
  transVector2[1] = yt;
  transVector2[2] = zt;
  quat2[0] = q0;
  quat2[1] = q1;
  quat2[2] = q2;
  quat2[3] = q3;

  /* printf("The second transl. vector is %f, %f, %f\n",
	 transVector2[0], transVector2[1], transVector2[2]);
  printf("The first quaternion is %f, %f, %f, %f\n", quat2[0], quat2[1], 
	 quat2[2], quat2[3]); */
  quat_to_mat (quat2, rotmat2);

  /*printf("The second rot matrix is this:\n");
  mat_print(rotmat2);*/

  /* apply R(q1)*T2 */
  mat_apply(rotmat1, transVector2);

  /*printf("R(q1)*T2 is %f, %f, %f\n",
	 transVector2[0], transVector2[1], transVector2[2]);*/

  /* Final = R(q1)*T2 + T1 */
  transVectorFinal[0] = transVector1[0] + transVector2[0];
  transVectorFinal[1] = transVector1[1] + transVector2[1];
  transVectorFinal[2] = transVector1[2] + transVector2[2];

  /*printf("Final vector is T1+T2:  %f, %f, %f\n", transVectorFinal[0], 
	 transVectorFinal[1], transVectorFinal[2]);*/

  mat_mult(rotmatFinal, rotmat1, rotmat2);
  /*printf("Final matrix is this:\n");
  mat_print(rotmatFinal);
  */

  mat_to_quat(rotmatFinal, quat3);
  /*printf("Final quaternion is:  %f, %f, %f, %f\n", quat3[0], quat3[1],
	 quat3[2], quat3[3]);*/
  sprintf(str, "%f %f %f %f %f %f %f",
	  transVectorFinal[0], transVectorFinal[1],
	  transVectorFinal[2], quat3[0], quat3[1], quat3[2], quat3[3]);

  Tcl_SetResult(interp, str, TCL_VOLATILE);
  return (TCL_OK);
}


/* tcl_inv_tranform
 * 
 * this function composes the inverse of one transform with
 * another transform and returns the result 
 *
 * Originally by Brian Curless
 * Modified:  Afra Zomorodian 7/5/95
 * Essentially, given two transforms of (translation and rotation in form
 * of a quaternion), we have:  
 * Final Rotation = R1^transpose*R2
 * Final Translation = R1^transpose*(T2-T1)
 * NOTE:  Get_Double returns a double (we need to convert it to float)
 */
tcl_inv_transform(dummy, interp, argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char *argv[];
{
  char str[200];
  Quaternion quat1, quat2, quat3;
  Matrix rotmat1, rotmat2, rotmatFinal;
  Matrix transmat1, transmat2;
  Matrix mat1, mat2;
  Vector transVector1, transVector2, transVectorFinal;
  double xt, yt, zt, q0, q1, q2, q3;

  if (argc != 15) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }
  /* parse the transformation info */
  Tcl_GetDouble (interp, argv[1], &xt);
  Tcl_GetDouble (interp, argv[2], &yt);
  Tcl_GetDouble (interp, argv[3], &zt);
  Tcl_GetDouble (interp, argv[4], &q0);
  Tcl_GetDouble (interp, argv[5], &q1);
  Tcl_GetDouble (interp, argv[6], &q2);
  Tcl_GetDouble (interp, argv[7], &q3);
  
  /* Convert to float */
  transVector1[0] = xt;
  transVector1[1] = yt;
  transVector1[2] = zt;
  quat1[0] = q0;
  quat1[1] = q1;
  quat1[2] = q2;
  quat1[3] = q3;

  /*
  printf("The first transl. vector is %f, %f, %f\n",
	 transVector1[0], transVector1[1], transVector1[2]);
  printf("The first quaternion is %f, %f, %f, %f\n", quat1[0], quat1[1], 
	 quat1[2], quat1[3]);
  */
  quat_to_mat (quat1, rotmat1);
  /* printf("The first rot matrix is this:\n");
  mat_print(rotmat1);
  */

  /* take the transpose, cause that's what we need */
  mat_transpose(rotmat1);
  /* printf("Transposed, it's this:\n");
  mat_print(rotmat1);
   */
  /* Parse Info */
  Tcl_GetDouble (interp, argv[8], &xt);
  Tcl_GetDouble (interp, argv[9], &yt);
  Tcl_GetDouble (interp, argv[10], &zt);
  Tcl_GetDouble (interp, argv[11], &q0);
  Tcl_GetDouble (interp, argv[12], &q1);
  Tcl_GetDouble (interp, argv[13], &q2);
  Tcl_GetDouble (interp, argv[14], &q3);

  /* Convert to float */
  transVector2[0] = xt;
  transVector2[1] = yt;
  transVector2[2] = zt;
  quat2[0] = q0;
  quat2[1] = q1;
  quat2[2] = q2;
  quat2[3] = q3;

  /*
  printf("The second transl. vector is %f, %f, %f\n",
	 transVector2[0], transVector2[1], transVector2[2]);
  printf("The first quaternion is %f, %f, %f, %f\n", quat2[0], quat2[1], 
	 quat2[2], quat2[3]);
  */
  quat_to_mat (quat2, rotmat2);

  /* printf("The second rot matrix is this:\n");
  mat_print(rotmat2); */

  /* Find T2-T1 */
  transVectorFinal[0] = -transVector1[0] + transVector2[0];
  transVectorFinal[1] = -transVector1[1] + transVector2[1];
  transVectorFinal[2] = -transVector1[2] + transVector2[2];
  /* printf("T2-T1 is %f, %f, %f\n",
	 transVector2[0], transVector2[1], transVector2[2]); */

  /* finally, we have Final Translation = R1^transpose(T2-T1) */
  mat_apply(rotmat1, transVectorFinal);

  /* printf("Final vector is: %f, %f, %f\n", transVectorFinal[0], 
	 transVectorFinal[1], transVectorFinal[2]); */

  mat_mult(rotmatFinal, rotmat1, rotmat2);
  /* printf("Final matrix is this:\n");
  mat_print(rotmatFinal);
   */
  mat_to_quat(rotmatFinal, quat3);
  /* printf("Final quaternion is:  %f, %f, %f, %f\n", quat3[0], quat3[1],
	 quat3[2], quat3[3]); */
  sprintf(str, "%f %f %f %f %f %f %f",
	  transVectorFinal[0], transVectorFinal[1],
	  transVectorFinal[2], quat3[0], quat3[1], quat3[2], quat3[3]);

  Tcl_SetResult(interp, str, TCL_VOLATILE);
  return (TCL_OK);
}
