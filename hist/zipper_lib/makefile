# Zipper together multiple meshes that were created from depth maps

#TCL_LIB = /usr/common/lib
#TCL_INCLUDE = /usr/common/include
TCL_LIB = /usr/lib/x86_64-linux-gnu/
TCL_INCLUDE = /usr/include/tcl8.6

#COPTIONS = -g -cckr
COPTIONS = -g 
#COPTIONS = -O2 -cckr
#CFLAGS = -I. -I$(TCL_INCLUDE) $(COPTIONS) -float -DIRIS -DHAVE_UNISTD_H
#CFLAGS = -I. -I$(TCL_INCLUDE) -I/usr/lib/perl5/Tk/pTk $(COPTIONS) 
CFLAGS = -I. -I$(TCL_INCLUDE)  $(COPTIONS) 
LIBS = libply.a libmeschach.a $(TCL_LIB)/libtk.a $(TCL_LIB)/libtcl.a -lmpc -lmalloc \
	-lc_s -lfm_s -lgl_s -lm -limage

current: zipper

	 #triangulate.c invert.c picking.c consensus.c intersect.c \

ZIPPER_C = zipper.c myalloc.c cyfile.c strings.c event.c \
	 invert.c consensus.c intersect.c \
	 near.c match_parallel.c remove.c edges.c polyfile.c quaternion.c \
	 rotations.c close_match.c symmeig.c closest.o \
	 undo.c mesh.c meshops.c clip.c trackball.c vect.c vector.c matrix.c \
	 cursors.c fill.c raw.c ply_wrapper.c help.c zipglobal.c triangulate.c \
         tcl_3d_interact.c tcl_transform.c tcl_align.c \
	 tcl_edit.c tcl_mesh_ops.c tcl_files.c tcl_eat_merge.c tcl_stuff.c \
     picking.c view.c draw.c

ZIPPER_OBJS = $(ZIPPER_C:.c=.o)

zipper: $(ZIPPER_OBJS) makefile
	rm -f zipper
	cc $(CFLAGS) -o zipper $(ZIPPER_OBJS) $(LIBS)

.c.o:
	cc $(CFLAGS) -c $<

# DO NOT DELETE THIS LINE -- make depend depends on it.
