# wish.tcl --
#
# This script is invoked by the "wish" program whenever it starts up.
# It invokes initialization scripts for Tcl and Tk, then does a few
# wish-specific things like setting the window geometry, if one was
# specified.
#
# $Header: /usr/graphics/project/scanner/zipper/tcl_and_tk/wish.tcl,v 1.1 1995/02/10 23:58:04 curless Exp $ SPRITE (Berkeley)
#
# Copyright 1992 Regents of the University of California
# Permission to use, copy, modify, and distribute this
# software and its documentation for any purpose and without
# fee is hereby granted, provided that this copyright
# notice appears in all copies.  The University of California
# makes no representations about the suitability of this
# software for any purpose.  It is provided "as is" without
# express or implied warranty.
#

source [info library]/init.tcl
source $tk_library/tk.tcl

if [info exists geometry] {
    wm geometry . $geometry
}
