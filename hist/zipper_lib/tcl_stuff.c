/*

Set up TCL and TK stuff.

Greg Turk, July 1993

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
#include <tcl.h>
#include <tkConfig.h>
#include <tk.h>
#include <gl/gl.h>
#include <gl/device.h>
#include <matrix.h>
#include <cyfile.h>
#include <zipper.h>
#include <cursors.h>
#include <view.h>
#include <raw.h>
#include <limits.h>
#include "ply.h"

#ifndef EQSTR
#define EQSTR(x, y)  (strcmp((x),(y)) == 0)
#endif


/* externals */

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

Tcl_Interp *interp;

/*
#define MAX_ALIGN_DRAG 100
static Scan *drag_list[MAX_ALIGN_DRAG];
*/

static Tk_Window w;
static Tcl_DString command;	/* Used to assemble lines of terminal input
				 * into Tcl commands. */
static char errorExitCmd[] = "exit 1";
char *tcl_RcFileName = NULL;	/* Name of a user-specific startup script
				 * to source if the application is being run
				 * interactively (e.g. "~/.wishrc").  Set
				 * by Tcl_AppInit.  NULL means don't source
				 * anything ever. */


static int tty;			/* 0 = stdin from file, not 0 means terminal */

static void Prompt _ANSI_ARGS_((Tcl_Interp *interp, int partial));
static void StdinProc _ANSI_ARGS_((ClientData clientData, int mask));
static void DelayedMap _ANSI_ARGS_((ClientData clientData));
static void StructureProc _ANSI_ARGS_((ClientData clientData,
                                      XEvent *eventPtr));

/* 3D interaction */
extern int tcl_set_drawing _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_mesh_draw _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_move_object _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_undo _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_redo _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_cull _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_mouse_speed _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_anti_alias _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_set_fast_ops _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_set_singlebuffer _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_set_orthographic _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_parallel_max _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_material _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_spin _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_axes _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_set_background _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_reset_view _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));


/* Editing meshes */
extern int tcl_set_edit _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_edit_window _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_execute_edit _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));


/* Zipper = eat + merge */
extern int tcl_eat_action _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_merge _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_only_merge _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_zip _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_eat_one _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_eat_step _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_join _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_eat_edges _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_target _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_next_merge _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_intersect _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_ifinish _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));


/* File I/O */
extern int tcl_raw_mesh _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_bin_poly_write _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_bin_poly_read _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_ply_read _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_ply_write _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_ply_range_mesh _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_read_ply_mesh _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_polygon_file _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_intensity_write _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_normals_write _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_colors_write _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_fileout _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_new_scan _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));


/* Alignment */
extern int tcl_align _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_align_error _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_align_step _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_next_align _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_stop_align _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_anchor _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_besl_old _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_match_list _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_match _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_close _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));


/* "Misc." */
extern int tcl_zipparam _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_consensus _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_pair_consensus _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_quit _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_beep _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_help _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));


/* Mesh operations */
extern int tcl_mesh_resolution _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_fix_bows _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_fill_holes _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_split_test _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_quarter_mesh _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_swap_edges _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_smooth_vertices _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_sliver_remove _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_bad_aspect_remove _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_flat_verts_remove _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_clean_mesh _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_remove_mesh _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_rename_mesh _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_edge_remove _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_vertex_remove _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));


/* Transformations */
extern int tcl_camera _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_translate _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_rotate _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_absorb_transform _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_transform _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_inv_transform _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_print_positions _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));
extern int tcl_mprint_positions _ANSI_ARGS_((ClientData clientData,
  Tcl_Interp *interp, int argc, char **argv));


extern float get_zipper_resolution();
extern float get_max_edge_length_factor();
extern float get_fill_edge_length_factor();
extern float get_conf_angle();
extern float get_conf_exponent();
extern float get_conf_edge_count_factor();
extern int get_conf_edge_zero();
extern float get_align_near_dist_factor();
extern float get_align_near_cos();
extern float get_eat_near_dist_factor();
extern float get_eat_near_cos();
extern int get_eat_start_iters();
extern float get_eat_start_factor();    
extern float get_clip_near_dist_factor();
extern float get_clip_near_cos();
extern float get_clip_boundary_dist_factor();
extern float get_clip_boundary_cos();
extern float get_consensus_position_dist_factor();
extern float get_consensus_normal_dist_factor();
extern float get_consensus_jitter_dist_factor();
extern float get_range_data_sigma_factor();
extern float get_range_data_min_intensity();
extern int get_range_data_horizontal_erode();

extern set_zipper_resolution(float);



/******************************************************************************
Set up widgets to control movie.
******************************************************************************/

setup_tcl(int argc, char **argv)
{
  int result;
  int code;
  char progname[80];
  long xorg,yorg;
  long xsize,ysize;
  char str[200];
  char *args;
  char buf[20];
  char *zipper_dir;
  static char *name = NULL;
  static char *display = NULL;

  /* create interpreter */
  interp = Tcl_CreateInterp();

  /*
   * Make command-line arguments available in the Tcl variables "argc"
   * and "argv".  Also set the "geometry" variable from the geometry
   * specified on the command line.
   */
  
  args = Tcl_Merge(argc-1, argv+1);
  Tcl_SetVar(interp, "argv", args, TCL_GLOBAL_ONLY);
  ckfree(args);
  sprintf(buf, "%d", argc-1);
  Tcl_SetVar(interp, "argc", buf, TCL_GLOBAL_ONLY);

/* Should be something like

  Tcl_SetVar(interp, "argv0", (fileName != NULL) ? fileName : argv[0],
 	     TCL_GLOBAL_ONLY);

   but instead: */
  Tcl_SetVar(interp, "argv0", argv[0], TCL_GLOBAL_ONLY);


  /* add a new Tcl commands */

  Tcl_CreateCommand (interp, "zipparam", tcl_zipparam, (ClientData) w,
	      (void (*)()) NULL);


  Tcl_CreateCommand (interp, "transform", tcl_transform, (ClientData) w,
	      (void (*)()) NULL);
  Tcl_CreateCommand (interp, "inv_transform", tcl_inv_transform, (ClientData) w,
	      (void (*)()) NULL);

  Tcl_CreateCommand (interp, "set_drawing", tcl_set_drawing, (ClientData) w,
	      (void (*)()) NULL);
  Tcl_CreateCommand (interp, "set_edit", tcl_set_edit, (ClientData) w,
	      (void (*)()) NULL);

  Tcl_CreateCommand (interp, "match", tcl_match, (ClientData) w,
	      (void (*)()) NULL);
  Tcl_CreateCommand (interp, "close_proc", tcl_close, (ClientData) w,
	      (void (*)()) NULL);
  Tcl_CreateCommand (interp, "print", tcl_print_positions,
	      (ClientData) w, (void (*)()) NULL);
  Tcl_CreateCommand (interp, "mprint", tcl_mprint_positions,
	      (ClientData) w, (void (*)()) NULL);
  Tcl_CreateCommand (interp, "pwrite", tcl_fileout, (ClientData) w,
	      (void (*)()) NULL);

  Tcl_CreateCommand (interp, "eat_action", tcl_eat_action, (ClientData) w,
	      (void (*)()) NULL);
  Tcl_CreateCommand (interp, "only_merge", tcl_only_merge, (ClientData) w,
	      (void (*)()) NULL);
  Tcl_CreateCommand (interp, "merge", tcl_merge, (ClientData) w,
	      (void (*)()) NULL);
  Tcl_CreateCommand (interp, "zip", tcl_zip, (ClientData) w,
	      (void (*)()) NULL);

  Tcl_CreateCommand (interp, "mesh_resolution", tcl_mesh_resolution,
              (ClientData) w, (void (*)()) NULL);
  Tcl_CreateCommand (interp, "mesh_draw", tcl_mesh_draw, (ClientData) w,
	      (void (*)()) NULL);
  Tcl_CreateCommand (interp, "move_object", tcl_move_object, (ClientData) w,
	      (void (*)()) NULL);

  Tcl_CreateCommand (interp, "mesh", tcl_new_scan, (ClientData) w,
	      (void (*)()) NULL);
  Tcl_CreateCommand (interp, "rmesh", tcl_raw_mesh, (ClientData) w,
	      (void (*)()) NULL);
  Tcl_CreateCommand (interp, "bmesh", tcl_read_ply_mesh, (ClientData) w,
	      (void (*)()) NULL);
/*
  Tcl_CreateCommand (interp, "bmesh", tcl_ply_range_mesh, (ClientData) w,
	      (void (*)()) NULL);
*/
  Tcl_CreateCommand (interp, "match_list", tcl_match_list, (ClientData) w,
	      (void (*)()) NULL);
  Tcl_CreateCommand (interp, "eat_edges", tcl_eat_edges, (ClientData) w,
	      (void (*)()) NULL);
  Tcl_CreateCommand (interp, "eat_one", tcl_eat_one, (ClientData) w,
	      (void (*)()) NULL);
  Tcl_CreateCommand (interp, "eat_step", tcl_eat_step, (ClientData) w,
	      (void (*)()) NULL);
  Tcl_CreateCommand (interp, "join", tcl_join, (ClientData) w,
	      (void (*)()) NULL);
  Tcl_CreateCommand (interp, "besl", tcl_besl_old, (ClientData) w,
	      (void (*)()) NULL);
  Tcl_CreateCommand (interp, "align_error", tcl_align_error, (ClientData) w,
	      (void (*)()) NULL);
  Tcl_CreateCommand (interp, "align", tcl_align, (ClientData) w,
	      (void (*)()) NULL);
  Tcl_CreateCommand (interp, "align_step", tcl_align_step, (ClientData) w,
	      (void (*)()) NULL);
  Tcl_CreateCommand (interp, "undo", tcl_undo, (ClientData) w,
	      (void (*)()) NULL);
  Tcl_CreateCommand (interp, "redo", tcl_redo, (ClientData) w,
	      (void (*)()) NULL);
  Tcl_CreateCommand (interp, "fix_bows", tcl_fix_bows, (ClientData) w,
	      (void (*)()) NULL);
  Tcl_CreateCommand (interp, "fill_holes", tcl_fill_holes, (ClientData) w,
	      (void (*)()) NULL);
  Tcl_CreateCommand (interp, "split_test", tcl_split_test, (ClientData) w,
	      (void (*)()) NULL);
  Tcl_CreateCommand (interp, "camera", tcl_camera, (ClientData) w,
	      (void (*)()) NULL);
  Tcl_CreateCommand (interp, "cull", tcl_cull, (ClientData) w,
	      (void (*)()) NULL);
  Tcl_CreateCommand (interp, "vremove", tcl_vertex_remove, (ClientData) w,
	      (void (*)()) NULL);
  Tcl_CreateCommand (interp, "eremove", tcl_edge_remove, (ClientData) w,
	      (void (*)()) NULL);
  Tcl_CreateCommand (interp, "sremove", tcl_sliver_remove, (ClientData) w,
	      (void (*)()) NULL);
  Tcl_CreateCommand (interp, "flatremove", tcl_flat_verts_remove, (ClientData) w,
	      (void (*)()) NULL);
  Tcl_CreateCommand (interp, "badasremove", tcl_bad_aspect_remove, (ClientData) w,
	      (void (*)()) NULL);
  Tcl_CreateCommand (interp, "speed", tcl_mouse_speed, (ClientData) w,
	      (void (*)()) NULL);
  Tcl_CreateCommand (interp, "edit_window", tcl_edit_window, (ClientData) w,
	      (void (*)()) NULL);
  Tcl_CreateCommand (interp, "quarter", tcl_quarter_mesh, (ClientData) w,
	      (void (*)()) NULL);
  Tcl_CreateCommand (interp, "consensus", tcl_consensus, (ClientData) w,
	      (void (*)()) NULL);
  Tcl_CreateCommand (interp, "pair_consensus", tcl_pair_consensus, (ClientData) w,
	      (void (*)()) NULL);
  Tcl_CreateCommand (interp, "quit", tcl_quit, (ClientData) w,
	      (void (*)()) NULL);
  Tcl_CreateCommand (interp, "polygon_file", tcl_polygon_file, (ClientData) w,
	      (void (*)()) NULL);
  Tcl_CreateCommand (interp, "swap_edges", tcl_swap_edges, (ClientData) w,
	      (void (*)()) NULL);
  Tcl_CreateCommand (interp, "smooth_vertices", tcl_smooth_vertices,
	      (ClientData) w, (void (*)()) NULL);
  Tcl_CreateCommand (interp, "execute_edit", tcl_execute_edit,
	      (ClientData) w, (void (*)()) NULL);
  Tcl_CreateCommand (interp, "intersect", tcl_intersect,
	      (ClientData) w, (void (*)()) NULL);
  Tcl_CreateCommand (interp, "ifinish", tcl_ifinish,
	      (ClientData) w, (void (*)()) NULL);
  Tcl_CreateCommand (interp, "old_bpolygon", tcl_bin_poly_read,
	      (ClientData) w, (void (*)()) NULL);
  Tcl_CreateCommand (interp, "old_bpwrite", tcl_bin_poly_write,
	      (ClientData) w, (void (*)()) NULL);
  Tcl_CreateCommand (interp, "clean_mesh", tcl_clean_mesh,
	      (ClientData) w, (void (*)()) NULL);
  Tcl_CreateCommand (interp, "anti_alias", tcl_anti_alias,
	      (ClientData) w, (void (*)()) NULL);
  Tcl_CreateCommand (interp, "beep", tcl_beep,
	      (ClientData) w, (void (*)()) NULL);
  Tcl_CreateCommand (interp, "target", tcl_target,
	      (ClientData) w, (void (*)()) NULL);
  Tcl_CreateCommand (interp, "next_merge", tcl_next_merge,
	      (ClientData) w, (void (*)()) NULL);
  Tcl_CreateCommand (interp, "anchor", tcl_anchor,
	      (ClientData) w, (void (*)()) NULL);
  Tcl_CreateCommand (interp, "next_align", tcl_next_align,
	      (ClientData) w, (void (*)()) NULL);
  Tcl_CreateCommand (interp, "spin", tcl_spin,
	      (ClientData) w, (void (*)()) NULL);
  Tcl_CreateCommand (interp, "axes", tcl_axes,
	      (ClientData) w, (void (*)()) NULL);
  Tcl_CreateCommand (interp, "bpolygon", tcl_read_ply_mesh,
	      (ClientData) w, (void (*)()) NULL);
/*
  Tcl_CreateCommand (interp, "bpolygon", tcl_ply_read,
	      (ClientData) w, (void (*)()) NULL);
*/
  Tcl_CreateCommand (interp, "bpwrite", tcl_ply_write,
	      (ClientData) w, (void (*)()) NULL);
  Tcl_CreateCommand (interp, "background", tcl_set_background,
	      (ClientData) w, (void (*)()) NULL);
  Tcl_CreateCommand (interp, "reset_view", tcl_reset_view,
	      (ClientData) w, (void (*)()) NULL);
  Tcl_CreateCommand (interp, "remove_mesh", tcl_remove_mesh,
	      (ClientData) w, (void (*)()) NULL);
  Tcl_CreateCommand (interp, "rename_mesh", tcl_rename_mesh,
	      (ClientData) w, (void (*)()) NULL);
  Tcl_CreateCommand (interp, "intensity_write", tcl_intensity_write,
	      (ClientData) w, (void (*)()) NULL);
  Tcl_CreateCommand (interp, "normals_write", tcl_normals_write,
	      (ClientData) w, (void (*)()) NULL);
  Tcl_CreateCommand (interp, "colors_write", tcl_colors_write,
	      (ClientData) w, (void (*)()) NULL);
  Tcl_CreateCommand (interp, "translate", tcl_translate,
	      (ClientData) w, (void (*)()) NULL);
  Tcl_CreateCommand (interp, "xrotate", tcl_rotate,
	      (ClientData) w, (void (*)()) NULL);
  Tcl_CreateCommand (interp, "yrotate", tcl_rotate,
	      (ClientData) w, (void (*)()) NULL);
  Tcl_CreateCommand (interp, "zrotate", tcl_rotate,
	      (ClientData) w, (void (*)()) NULL);
  Tcl_CreateCommand (interp, "qrotate", tcl_rotate,
	      (ClientData) w, (void (*)()) NULL);
  Tcl_CreateCommand (interp, "absorb_transform", tcl_absorb_transform,
	      (ClientData) w, (void (*)()) NULL);
  Tcl_CreateCommand (interp, "stop_align", tcl_stop_align,
	      (ClientData) w, (void (*)()) NULL);
  Tcl_CreateCommand (interp, "fast_ops", tcl_set_fast_ops,
	      (ClientData) w, (void (*)()) NULL);
  Tcl_CreateCommand (interp, "singlebuffer", tcl_set_singlebuffer,
	      (ClientData) w, (void (*)()) NULL);
  Tcl_CreateCommand (interp, "orthographic", tcl_set_orthographic,
	      (ClientData) w, (void (*)()) NULL);
  Tcl_CreateCommand (interp, "parallel_max", tcl_parallel_max,
	      (ClientData) w, (void (*)()) NULL);
  Tcl_CreateCommand (interp, "material", tcl_material,
	      (ClientData) w, (void (*)()) NULL);
  Tcl_CreateCommand (interp, "help", tcl_help,
	      (ClientData) w, (void (*)()) NULL);
  Tcl_CreateCommand (interp, "?", tcl_help,
	      (ClientData) w, (void (*)()) NULL);

  /* create main window */
  strcpy (progname, "Button Window");
  name = progname;
  w = Tk_CreateMainWindow (interp, display, name, "Tk");
  if (w == NULL) {
    fprintf (stderr, "%s\n", interp->result);
    exit (-1);
  }


  Tk_SetClass(w, "Tk");
  Tk_CreateEventHandler(w, StructureNotifyMask, StructureProc,
      (ClientData) NULL);
  Tk_DoWhenIdle(DelayedMap, (ClientData) NULL);

  /* read commands from file */
  zipper_dir = getenv ("ZIPPER_DIR");
  if (zipper_dir == NULL)
    /* look for file in current directory */
    code = Tcl_EvalFile (interp, "zipper.tcl");
  else {
    /* look for file in ZIPPER_DIR directory */
    strcpy (str, zipper_dir);
    strcat (str, "/zipper.tcl");
    code = Tcl_EvalFile (interp, str);
  }

  if (*interp->result != 0) {
    printf ("%s\n", interp->result);
  }
  if (code != TCL_OK) {
    fprintf (stderr, "error\n");
    exit (-1);
  }

  /* find location and size of the main drawing window */
  if (!global_dont_draw) {
      winset (draw_win);
      getorigin (&xorg, &yorg);
      getsize (&xsize, &ysize);
  } else {
      /* Dummy values; hacks until a no GUI version is in place */
      xorg = 500;
      yorg = 500;
      xsize = 200;
      ysize = 200;
  }

  /* position the widget window below the main window */
  sprintf (str, "wm geometry .topbuttons +%d+%d\n", xorg + 2, 1024 - yorg + 39);
  Tcl_Eval (interp, str, 0, (char **) NULL);

  /* position the scans window to the left of the main window */
  sprintf (str, "wm geometry .meshbuttons +%d+%d\n",
	   xorg - 350, 1024 - (yorg + 550));
  Tcl_Eval (interp, str, 0, (char **) NULL);

  /* position the scans window to the left of the main window */
  sprintf (str, "wm geometry .editbuttons +%d+%d\n",
	   xorg - 500, 1024 - (yorg + 550));
  Tcl_Eval (interp, str, 0, (char **) NULL);

  /* set up things to read commands from keyboard */
  tty = isatty(0);
  Tk_CreateFileHandler(0, TK_READABLE, StdinProc, (ClientData) 0);
  if (tty)
    printf("zipper> ");


  Tcl_SetVar(interp, "tcl_interactive",
	     (tty) ? "1" : "0", TCL_GLOBAL_ONLY);

  fflush(stdout);
  Tcl_DStringInit(&command);

  if (global_dont_draw) {
      Tcl_Eval(interp, "wm withdraw .meshbuttons");
      Tcl_Eval(interp, "wm withdraw .topbuttons");
      Tcl_Eval(interp, "wm withdraw .editbuttons");
      Tcl_Eval(interp, "wm withdraw .");
  }

#if 0
  /* event handler */
  Tk_MainLoop();
#endif
}


/*
 *----------------------------------------------------------------------
 *
 * StdinProc --
 *
 *	This procedure is invoked by the event dispatcher whenever
 *	standard input becomes readable.  It grabs the next line of
 *	input characters, adds them to a command being assembled, and
 *	executes the command if it's complete.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	Could be almost arbitrary, depending on the command that's
 *	typed.
 *
 *----------------------------------------------------------------------
 */

    /* ARGSUSED */
static void
StdinProc(clientData, mask)
    ClientData clientData;		/* Not used. */
    int mask;				/* Not used. */
{
#define BUFFER_SIZE 4000
    char input[BUFFER_SIZE+1];
    static int gotPartial = 0;
    char *cmd;
    int code, count;

    count = read(fileno(stdin), input, BUFFER_SIZE);
    if (count <= 0) {
	if (!gotPartial) {
	    if (tty) {
		Tcl_Eval(interp, "exit");
		exit(1);
	    } else {
		Tk_DeleteFileHandler(0);
	    }
	    return;
	} else {
	    count = 0;
	}
    }
    cmd = Tcl_DStringAppend(&command, input, count);
    if (count != 0) {
	if ((input[count-1] != '\n') && (input[count-1] != ';')) {
	    gotPartial = 1;
	    goto prompt;
	}
	if (!Tcl_CommandComplete(cmd)) {
	    gotPartial = 1;
	    goto prompt;
	}
    }
    gotPartial = 0;

    /*
     * Disable the stdin file handler while evaluating the command;
     * otherwise if the command re-enters the event loop we might
     * process commands from stdin before the current command is
     * finished.  Among other things, this will trash the text of the
     * command being evaluated.
     */

    Tk_CreateFileHandler(0, 0, StdinProc, (ClientData) 0);
    code = Tcl_RecordAndEval(interp, cmd, 0);
    Tk_CreateFileHandler(0, TK_READABLE, StdinProc, (ClientData) 0);
    Tcl_DStringFree(&command);
    if (*interp->result != 0) {
	if ((code != TCL_OK) || (tty)) {
	    printf("%s\n", interp->result);
	}
    }

    /*
     * Output a prompt.
     */

    prompt:
    if (tty) {
	Prompt(interp, gotPartial);
    }
}


/*
 *----------------------------------------------------------------------
 *
 * Prompt --
 *
 *	Issue a prompt on standard output, or invoke a script
 *	to issue the prompt.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	A prompt gets output, and a Tcl script may be evaluated
 *	in interp.
 *
 *----------------------------------------------------------------------
 */

static void
Prompt(interp, partial)
    Tcl_Interp *interp;			/* Interpreter to use for prompting. */
    int partial;			/* Non-zero means there already
					 * exists a partial command, so use
					 * the secondary prompt. */
{
    char *promptCmd;
    int code;

    promptCmd = Tcl_GetVar(interp,
	partial ? "tcl_prompt2" : "tcl_prompt1", TCL_GLOBAL_ONLY);
    if (promptCmd == NULL) {
	defaultPrompt:
	if (!partial) {
	    fputs("zipper> ", stdout);
	}
    } else {
	code = Tcl_Eval(interp, promptCmd);
	if (code != TCL_OK) {
	    Tcl_AddErrorInfo(interp,
		    "\n    (script that generates prompt)");
	    fprintf(stderr, "%s\n", interp->result);
	    goto defaultPrompt;
	}
    }
    fflush(stdout);
}



/******************************************************************************
This procedure is invoked by the event dispatcher whenever standard input
becomes readable.  It grabs the next line of input characters, adds them
to a command being assembled, and executes the command if it's complete.
******************************************************************************/
 
/*
static void StdinProc(clientData, mask)
  ClientData clientData;
  int mask;
{
#define BUFFER_SIZE 4000
    char input[BUFFER_SIZE+1];
    static int gotPartial = 0;
    char *cmd;
    int result, count;

    count = read(fileno(stdin), input, BUFFER_SIZE);
    if (count <= 0) {
	if (!gotPartial) {
	    if (tty) {
		Tcl_Eval(interp, "destroy .", 0, (char **) NULL);
		exit(0);
	    } else {
		Tk_DeleteFileHandler(0);
	    }
	    return;
	} else {
	    input[0] = 0;
	}
    } else {
	input[count] = 0;
    }
    cmd = Tcl_AssembleCmd(buffer, input);
    if (cmd == NULL) {
	gotPartial = 1;
	return;
    }
    gotPartial = 0;
    result = Tcl_RecordAndEval(interp, cmd, 0);
    if (*interp->result != 0) {
	if ((result != TCL_OK) || (tty)) {
	    printf("%s\n", interp->result);
	}
    }
    if (tty) {
	printf("zipper> ");
	fflush(stdout);
    }
}
*/


/******************************************************************************
This procedure is invoked whenever a structure-related event
occurs on the main window.  If the window is deleted, the
procedure modifies "w" to record that fact.
******************************************************************************/

static void StructureProc(clientData, eventPtr)
  ClientData clientData;      /* Information about window. */
  XEvent *eventPtr;           /* Information about event. */
{
  if (eventPtr->type == DestroyNotify) {
    w = NULL;
    exit (0);
  }
}


/******************************************************************************
This procedure is invoked by the event dispatcher once the
startup script has been processed.  It waits for all other
pending idle handlers to be processed (so that all the
geometry information will be correct), then maps the
application's main window.
******************************************************************************/

static void DelayedMap(clientData)
  ClientData clientData;      /* Not used. */
{

  while (Tk_DoOneEvent(TK_IDLE_EVENTS) != 0) {
    /* Empty loop body. */
  }
  if (w == NULL) {
    return;
  }
  Tk_MapWindow(w);
}


int tcl_zipparam(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
    float fval;
    int ival, count;

    if (argc == 1) {
	printf("Command: zipparam [-option value]\n");
	printf("  -zipper_resolution <float> (%g)\n", get_zipper_resolution());
	printf("  -max_edge_length_factor <float> (%g)\n", 
	       get_max_edge_length_factor());
	printf("  -fill_edge_length_factor <float> (%g)\n", 
	       get_fill_edge_length_factor());
	printf("  -conf_angle <float> (%g)\n", 
	       get_conf_angle());
	printf("  -conf_exponent <float> (%g)\n", 
	       get_conf_exponent());
	printf("  -conf_edge_count_factor <float> (%g)\n", 
	       get_conf_edge_count_factor());
	printf("  -conf_edge_zero <boolean> (%d)\n", 
	       get_conf_edge_zero());
	printf("  -align_near_dist_factor <float> (%g)\n", 
	       get_align_near_dist_factor());
	printf("  -align_near_cos <float> (%g)\n", 
	       get_align_near_cos());
	printf("  -eat_near_dist_factor <float> (%g)\n", 
	       get_eat_near_dist_factor());
	printf("  -eat_near_cos <float> (%g)\n", get_eat_near_cos());
	printf("  -eat_start_iters <int> (%d)\n", get_eat_start_iters());
	printf("  -eat_start_factor <float> (%g)\n", get_eat_start_factor());
	printf("  -clip_near_dist_factor <float> (%g)\n", 
	       get_clip_near_dist_factor());
	printf("  -clip_near_cos <float> (%g)\n", get_clip_near_cos());
	printf("  -clip_boundary_dist_factor <float> (%g)\n", 
	       get_clip_boundary_dist_factor());
	printf("  -clip_boundary_cos <float> (%g)\n", get_clip_boundary_cos());
	printf("  -consensus_position_dist_factor <float> (%g)\n", 
	       get_consensus_position_dist_factor());
	printf("  -consensus_normal_dist_factor <float> (%g)\n", 
	       get_consensus_normal_dist_factor());
	printf("  -consensus_jitter_dist_factor <float> (%g)\n", 
	       get_consensus_jitter_dist_factor());
	printf("  -range_data_sigma_factor <float> (%g)\n", 
	       get_range_data_sigma_factor());
	printf("  -range_data_min_intensity <float> (%g)\n", 
	       get_range_data_min_intensity());
	printf("  -range_data_horizontal_erode <int> (%d)\n", 
	       get_range_data_horizontal_erode());
	return TCL_OK;
    }

    count = 1;
    while (count < argc) {
	if (EQSTR(argv[count], "-zipper_resolution")) {
	    count++;
	    fval = atof(argv[count]);
		set_zipper_resolution(fval);
	}
	else if (EQSTR(argv[count], "-max_edge_length_factor")) {
	    count++;
	    fval = atof(argv[count]);
	    set_max_edge_length_factor(fval);
	}
	else if (EQSTR(argv[count], "-fill_edge_length_factor")) {
	    count++;
	    fval = atof(argv[count]);
	    set_fill_edge_length_factor(fval);
	}
	else if (EQSTR(argv[count], "-conf_angle")) {
	    count++;
	    fval = atof(argv[count]);
	    set_conf_angle(fval);
	}
	else if (EQSTR(argv[count], "-conf_exponent")) {
	    count++;
	    fval = atof(argv[count]);
	    set_conf_exponent(fval);
	}
	else if (EQSTR(argv[count], "-conf_edge_count_factor")) {
	    count++;
	    fval = atof(argv[count]);
	    set_conf_edge_count_factor(fval);
	}
	else if (EQSTR(argv[count], "-conf_edge_zero")) {
	    count++;
	    ival = atoi(argv[count]);
	    set_conf_edge_zero(ival);
	}
	else if (EQSTR(argv[count], "-align_near_dist_factor")) {
	    count++;
	    fval = atof(argv[count]);
	    set_align_near_dist_factor(fval);
	}
	else if (EQSTR(argv[count], "-align_near_cos")) {
	    count++;
	    fval = atof(argv[count]);
	    set_align_near_cos(fval);
	}
	else if (EQSTR(argv[count], "-eat_near_dist_factor")) {
	    count++;
	    fval = atof(argv[count]);
	    set_eat_near_dist_factor(fval);
	}
	else if (EQSTR(argv[count], "-eat_near_cos")) {
	    count++;
	    fval = atof(argv[count]);
	    set_eat_near_cos(fval);
	}
	else if (EQSTR(argv[count], "-eat_start_iters")) {
	    count++;
	    ival = atoi(argv[count]);
	    set_eat_start_iters(ival);
	}
	else if (EQSTR(argv[count], "-eat_start_factor")) {
	    count++;
	    fval = atof(argv[count]);
	    set_eat_start_factor(fval);    
	}
	else if (EQSTR(argv[count], "-clip_near_dist_factor")) {
	    count++;
	    fval = atof(argv[count]);
	    set_clip_near_dist_factor(fval);
	}
	else if (EQSTR(argv[count], "-clip_near_cos")) {
	    count++;
	    fval = atof(argv[count]);
	    set_clip_near_cos(fval);
	}
	else if (EQSTR(argv[count], "-clip_boundary_dist_factor")) {
	    count++;
	    fval = atof(argv[count]);
	    set_clip_boundary_dist_factor(fval);
	}
	else if (EQSTR(argv[count], "-clip_boundary_cos")) {
	    count++;
	    fval = atof(argv[count]);
	    set_clip_boundary_cos(fval);
	}
	else if (EQSTR(argv[count], "-consensus_position_dist_factor")) {
	    count++;
	    fval = atof(argv[count]);
	    set_consensus_position_dist_factor(fval);
	}
	else if (EQSTR(argv[count], "-consensus_normal_dist_factor")) {
	    count++;
	    fval = atof(argv[count]);
	    set_consensus_normal_dist_factor(fval);
	}
	else if (EQSTR(argv[count], "-consensus_jitter_dist_factor")) {
	    count++;
	    fval = atof(argv[count]);
	    set_consensus_jitter_dist_factor(fval);
	}
	else if (EQSTR(argv[count], "-range_data_sigma_factor")) {
	    count++;
	    fval = atof(argv[count]);
	    set_range_data_sigma_factor(fval);
	}
	else if (EQSTR(argv[count], "-range_data_min_intensity")) {
	    count++;
	    fval = atof(argv[count]);
	    set_range_data_min_intensity(fval);
	}
	else if (EQSTR(argv[count], "-range_data_horizontal_erode")) {
	    count++;
	    ival = atoi(argv[count]);
	    set_range_data_horizontal_erode(ival);
	}
	else {
	    printf("Do not recognize option '%s'.\n", argv[count]);
	    count++;
	}
	count++;
    }
    return TCL_OK;
}



/******************************************************************************
Create a list of pointers to scans, given a list of scan names.

Entry:
  num_scans  - number of names in list of names
  scan_names - the list of names

Exit:
  returns a list of pointers to scans, or NULL if there was an error
******************************************************************************/

Scan **make_scan_list(num_scans,scan_names)
  int num_scans;
  char **scan_names;
{
  int i;
  Scan **list;

  /* allocate the list to return */
  list = (Scan **) myalloc (sizeof (Scan *) * num_scans);

  /* look for each scan name and get pointer to each */

  for (i = 0; i < num_scans; i++) {
    list[i] = find_scan (scan_names[i]);
    if (list[i] == NULL) {
      free (list);
      return (NULL);
    }
  }

  return (list);
}





/******************************************************************************
Create consensus surface position by weighted average of all scans.

Entry:
  argv[1] - name of scan containing mesh to average
  argv[2] - level of detail to use
  argv[3] - scaling factor for neighborhood search distance
******************************************************************************/

int tcl_consensus(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  Scan *sc;
  int level, num_scans, i;
  float k_scale;
  Scan **scan_list;
  char *scan_name;
  int *use_old_mesh;

  if (argc < 3) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }

  sc = find_scan (argv[1]);
  if (sc == NULL) {
      fprintf (stderr, "Can't find scan called %s\n", argv[1]);
      return (TCL_ERROR);
  }
  
  level = atoi (argv[2]);

  if (argc > 3) {
      num_scans = argc-3;
      scan_list = (Scan **)malloc(num_scans*sizeof(Scan*));
      use_old_mesh = (int *)malloc(num_scans*sizeof(int));
      for (i = 0; i < num_scans; i++) {
	  scan_name = argv[i+3];
	  scan_list[i] = find_scan (scan_name);
	  if (scan_list[i] == NULL) {
	      fprintf (stderr, "Can't find scan called %s\n", scan_name);
	      free(scan_list);
	      return (TCL_ERROR);
	  }
	  use_old_mesh[i] = 0;
      }
      new_consensus_surface (sc, scan_list, use_old_mesh, num_scans, level - 1);
      free(scan_list);
      free(use_old_mesh);
  } else {
      consensus_surface (sc, level - 1, 1.0);
  }

  draw_object();

  return (TCL_OK);
}



/******************************************************************************
Create consensus surface position by weighted average of all scans.

Entry:
  argv[1] - name of scan containing mesh to average
  argv[2] - level of detail to use
  argv[3] - scaling factor for neighborhood search distance
******************************************************************************/

int tcl_pair_consensus(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  Scan *sc1, *sc2;
  int level, num_scans, i;
  float k_scale;
  Scan *scan_list[2];
  int use_old_mesh[2];

  if (argc < 4) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }

  sc1 = find_scan (argv[1]);
  if (sc1 == NULL) {
      fprintf (stderr, "Can't find scan called %s\n", argv[1]);
      return (TCL_ERROR);
  }
  
  sc2 = find_scan (argv[2]);
  if (sc2 == NULL) {
      fprintf (stderr, "Can't find scan called %s\n", argv[1]);
      return (TCL_ERROR);
  }
  
  level = atoi (argv[3]);

  scan_list[0] = sc1;
  scan_list[1] = sc2;
  use_old_mesh[0] = 1;
  use_old_mesh[1] = 0;

  new_consensus_surface (sc2, scan_list, use_old_mesh, 2, level - 1);
  new_consensus_surface (sc1, scan_list, use_old_mesh, 2, level - 1);

  draw_object();

  return (TCL_OK);
}


/******************************************************************************
Quit the program.
******************************************************************************/

int tcl_quit(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
    extern char *configuration_filename;
    char str[PATH_MAX];

  if (configuration_filename != NULL) {
      printf("\nSave configuration file: %s? [y/n] ", configuration_filename);
      scanf("%s", str);
      if (!strcmp(str, "y")) {
	  sprintf(str, "print %s", configuration_filename);
	  Tcl_Eval (interp, str, 0, (char **) NULL);	  
      }
  }

  exit (0);
}



/******************************************************************************
Make a beep.
******************************************************************************/

int tcl_beep(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  printf ("beep!\007\n");

  return (TCL_OK);
}



/******************************************************************************
Print out help information.
******************************************************************************/

int tcl_help(dummy,interp,argc,argv)
  ClientData dummy;
  Tcl_Interp *interp;
  int argc;
  char **argv;
{
  if (argc > 2) {
    interp->result = "wrong number of arguments";
    return (TCL_ERROR);
  }

  if (argc == 1)
    help (NULL);

  if (argc == 2)
    help (argv[1]);
  
  return (TCL_OK);
}

