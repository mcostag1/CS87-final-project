#include <cyfile.h>
#include <tcl.h>
#include <zipper.h>
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


/******************************************************************************
Set the editing flag.
******************************************************************************/

int tcl_set_edit(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  char *str;
  int result;

  str = Tcl_GetVar (interp, "edit_state", 0);
  result = atoi (str);

  if (result == 0) {
    pick_flag = 0;
    nothing_picked();
  }
  else if (result == 1) {
    pick_flag = 1;
    nothing_picked();
  }
  else if (result == 2) {
    pick_flag = 2;
    nothing_picked();
  }
  else if (result == 3) {
    pick_flag = 3;
    nothing_picked();
  }
  else if (result == 4) {
    pick_flag = 4;
    nothing_picked();
  }
  else if (result == 5) {
    pick_flag = 5;
    nothing_picked();
  }
  else if (result == 6) {
    pick_flag = 6;
    nothing_picked();
  }
  else
    return (TCL_ERROR);

  /* set the appropriate cursor */
  set_pick_cursor();

  return (TCL_OK);
}


/******************************************************************************
Execute the correct edit command.
******************************************************************************/

int tcl_execute_edit(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  if (argc != 1) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }

  switch (pick_flag) {
    case 0:
      printf ("no editing feature selected\n");
      break;
    case 1:
      delete_picked_item();
      draw_object();
      break;
    case 2:
      delete_picked_item();
      draw_object();
      break;
    case 3:
      fill_picked_loop();
      draw_object();
      break;
    case 4:
      better_fill_picked_loop();
      draw_object();
      break;
    case 5:
      make_picked_tri();
      break;
    case 6:
      delete_in_box();
      break;
    default:
      printf ("tcl_execute_edit: invalid switch value = %d\n", pick_flag);
      break;
  }

  return (TCL_OK);
}


/******************************************************************************
Toggle whether the edit window is being displayed.
******************************************************************************/

int tcl_edit_window(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  char str[80];
  static int edit_toggle = 0;

  edit_toggle = 1 - edit_toggle;

  if (edit_toggle) {
    /* set the pick flag to the appropriate value, based on radio buttons */
    sprintf (str, "set_edit");
    Tcl_Eval (interp, str, 0, (char **) NULL);
    /* display the edit controls */
    sprintf (str, "wm deiconify .editbuttons");
    Tcl_Eval (interp, str, 0, (char **) NULL);
    set_pick_cursor();
  }
  else {
    /* clear the overlay plane */
    drawmode (OVERDRAW);
    color (0);
    clear();
    drawmode (NORMALDRAW);
    /* un-display the edit controls */
    sprintf (str, "wm withdraw .editbuttons");
    Tcl_Eval (interp, str, 0, (char **) NULL);
    pick_flag = 0;
    set_move_cursor();
  }

  return (TCL_OK);
}
