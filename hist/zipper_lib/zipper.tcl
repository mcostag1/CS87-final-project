# Zippering polygonal meshes

source [info library]/init.tcl
source $tk_library/tk.tcl

# don't display the top-level window
wm withdraw .

# create a separate window for the buttons
set w .topbuttons
toplevel $w
wm title $w "Button Window"

# pick a nice, big font
set fnt "-*-new century schoolbook-medium-r-*-*-*-140-*-*-*-*-*-*"
option add *Button.font $fnt
option add *Radiobutton.font $fnt
option add *Checkbutton.font $fnt
option add *Label.font $fnt

set smallfont "-*-new century schoolbook-medium-r-*-*-*-100-*-*-*-*-*-*"

# make radiobuttons have different color than checkbuttons
option add *Radiobutton.selector limegreen

set bset {top pady 4 frame w}

frame $w.frame1 -borderwidth 10
pack append $w.frame1 \
    [checkbutton $w.frame1.b1 -relief flat -text "Polygons" \
     -command { set_drawing polygons $polygon_state } \
     -variable polygon_state ] $bset \
    [checkbutton $w.frame1.b2 -relief flat -text "Smooth" \
     -command { set_drawing smooth $smooth_state } \
     -variable smooth_state ] $bset \
    [checkbutton $w.frame1.b3 -relief flat -text "FalseColor" \
     -command { set_drawing colors $colors_state } \
     -variable colors_state ] $bset \
    [checkbutton $w.frame1.b4 -relief flat -text "TrueColor" \
     -command { set_drawing true_color $tcolor_state } \
     -variable tcolor_state ] $bset

frame $w.frame4 -borderwidth 10
pack append $w.frame4 \
    [checkbutton $w.frame4.b1 -relief flat -text "Confidence" \
     -command { set_drawing value $value_state } \
     -variable value_state ] $bset \
    [checkbutton $w.frame4.b2 -relief flat -text "Normals" \
     -command { set_drawing normal $normal_state } \
     -variable normal_state ] $bset \
    [checkbutton $w.frame4.b3 -relief flat -text "Bounds" \
     -command { set_drawing bounds $bounds_state } \
     -variable bounds_state ] $bset \
    [checkbutton $w.frame4.b4 -relief flat -text "Refine" \
     -command { set_drawing refine $refine_state } \
     -variable refine_state ] $bset

frame $w.frame6 -borderwidth 10
pack append $w.frame6 \
    [checkbutton $w.frame6.b1 -relief flat -text "Update" \
     -command { fast_ops [expr !$update_state] } \
     -variable update_state ] $bset \
    [checkbutton $w.frame6.b2 -relief flat -text "Orthographic" \
     -command { orthographic $ortho_state } \
     -variable ortho_state ] $bset \
    [checkbutton $w.frame6.b3 -relief flat -text "Singlebuffer" \
     -command { singlebuffer $singlebuffer_state } \
     -variable singlebuffer_state ] $bset \
    [checkbutton $w.frame6.b4 -relief flat -text "Shiny" \
     -command { set_shiny $shiny_state } \
     -variable shiny_state ] $bset


#    [checkbutton $w.frame4.b2 -relief flat -text "Intensity" \
#     -command { set_drawing intensity $intensity_state } \
#     -variable intensity_state ] $bset \

#    [checkbutton $w.frame4.b3 -relief flat -text "Nearest" \
#     -command { set_drawing nearest $nearest_state } \
#     -variable nearest_state ] $bset \

#    [checkbutton $w.frame4.b2 -relief flat -text "Edges" \
#     -command { set_drawing edge $edges_state } \
#     -variable edges_state ] $bset \

set polygon_state 1
set smooth_state 1
set colors_state 1
set tcolor_state 0
set refine_state 1
set value_state 0
set intensity_state 0
set shiny_state 0 
set singlebuffer_state 0 
set orthographic_state 0
set update_state 1

set bset {top pady 4 fill}

frame $w.frame2 -borderwidth 10
pack append $w.frame2 \
    [button $w.frame2.b1 -text "Edit" -command edit_window] $bset \
    [button $w.frame2.b2 -text "Undo" -command undo ] $bset \
    [button $w.frame2.b3 -text "Redo" -command redo ] $bset

#   [button $w.frame2.b4 -text "Print" -command print_positions] $bset

frame $w.frame3 -borderwidth 10
pack append $w.frame3 \
    [button $w.frame3.b1 -text "Align" -command align] $bset \
    [button $w.frame3.b2 -text "StopAlign" -command stop_align] $bset \
    [button $w.frame3.b3 -text "Eat Edge" -command eat_edges] $bset \
    [button $w.frame3.b4 -text "Merge" -command only_merge] $bset

.topbuttons.frame3.b2 configure -state disabled

#   [button $w.frame3.b4 -text "Fileout" -command pwrite] $bset
#   [button $w.frame3.b3 -text "Zip" -command zip] $bset

frame $w.frame5 -borderwidth 10
pack append $w.frame5 \
    [radiobutton $w.frame5.b1 -relief flat -text "Mesh1" \
     -variable mesh_resolution -value 0 \
     -command mesh_res] $bset \
    [radiobutton $w.frame5.b2 -relief flat -text "Mesh2" \
     -variable mesh_resolution -value 1 \
     -command mesh_res] $bset \
    [radiobutton $w.frame5.b3 -relief flat -text "Mesh3" \
     -variable mesh_resolution -value 2 \
     -command mesh_res] $bset \
    [radiobutton $w.frame5.b4 -relief flat -text "Mesh4" \
     -variable mesh_resolution -value 3 \
     -command mesh_res] $bset

set mesh_resolution 3

pack append $w \
  $w.frame6 {left expand fill} \
  $w.frame1 {left expand fill} \
  $w.frame4 {left expand fill} \
  $w.frame2 {left expand fill} \
  $w.frame3 {left expand fill} \
  $w.frame5 {left expand fill}

# create a window for mesh buttons
set w .meshbuttons
toplevel $w
wm title $w "Selection Window"

set bset {right}

frame $w.frameobjlbl -borderwidth 10
label $w.frameobjlbl.label -font $smallfont -text "Vis   Mv   Anc   Aln   Tgt  Mrg  "
#label $w.frameobjlbl.label -font $smallfont -text "Vis    Mv  Anc    Al    Tgt    Mg  "

pack append $w.frameobjlbl $w.frameobjlbl.label $bset

frame $w.obj_list

frame $w.frameobj -borderwidth 10
label $w.frameobj.label -text "Viewer       "
radiobutton $w.frameobj.radio -relief flat -text "                             " -variable mesh_move \
  -value -1 -command {move_object world}

pack append $w.frameobj \
	$w.frameobj.radio {right}

pack append $w.frameobj \
	$w.frameobj.label {right}

frame $w.frameobj2 -borderwidth 10
label $w.frameobj2.label -text "Light       "
radiobutton $w.frameobj2.radio -relief flat -text "                             " -variable mesh_move \
  -value -2 -command {move_object light}
pack append $w.frameobj2 \
  $w.frameobj2.radio $bset \
  $w.frameobj2.label $bset 

frame $w.frameobjcmds -borderwidth 10
button $w.frameobjcmds.align -text "Align" -command align
button $w.frameobjcmds.stopalign -text "StopAlign" -command stop_align
button $w.frameobjcmds.eat -text "Eat" -command eat_edges
button $w.frameobjcmds.merge -text "Merge" -command only_merge
button $w.frameobjcmds.eat_merge -text "Eat/Merge" -command merge
pack append $w.frameobjcmds \
	$w.frameobjcmds.align {left padx 10 fill} \
	$w.frameobjcmds.stopalign {left padx 10  fill} \
	$w.frameobjcmds.eat {left padx 10  fill} \
	$w.frameobjcmds.merge {left padx 10  fill} \
	$w.frameobjcmds.eat_merge {left padx 10  fill}

pack append $w \
  $w.frameobj {bottom expand fill} \
  $w.frameobj2 {bottom expand fill} \
  $w.obj_list {bottom expand fill} \
  $w.frameobjlbl {bottom expand fill} 

#pack append $w \
#  $w.frameobjcmds {bottom expand fill} \
#  $w.frameobj {bottom expand fill} \
#  $w.frameobj2 {bottom expand fill} \
#  $w.obj_list {bottom expand fill} \
#  $w.frameobjlbl {bottom expand fill} 

set mesh_move -1
set anchor_id -1
set align_id -1
set target_id -1
set merge_id -1

proc make_mesh_buttons {w name index} \
{
  global smallfont

  frame $w.obj_list.frame$index -borderwidth 10
  label $w.obj_list.frame$index.label -text $name

  checkbutton $w.obj_list.frame$index.vis -text "" \
    -relief flat -variable mesh_draw_$index \
    -command "mesh_draw $name \$mesh_draw_$index"
  radiobutton $w.obj_list.frame$index.radio -text "" \
     -relief flat -variable mesh_move  -value $name -command "move_object $name"
  radiobutton $w.obj_list.frame$index.anchor -text "" \
    -relief flat -variable anchor_id  -value $name -command "anchor $name"
  radiobutton $w.obj_list.frame$index.align -text "" \
    -relief flat -variable align_id  -value $name -command "next_align $name"
  radiobutton $w.obj_list.frame$index.target -text "" \
    -relief flat -variable target_id  -value $name -command "target $name"
  radiobutton $w.obj_list.frame$index.merge -text "" \
    -relief flat -variable merge_id  -value $name -command "next_merge $name"

  set bset {right}
  pack append $w.obj_list.frame$index \
    $w.obj_list.frame$index.merge $bset \
    $w.obj_list.frame$index.target $bset \
    $w.obj_list.frame$index.align $bset \
    $w.obj_list.frame$index.anchor $bset \
    $w.obj_list.frame$index.radio $bset \
    $w.obj_list.frame$index.vis $bset \
    $w.obj_list.frame$index.label $bset
  pack append $w.obj_list $w.obj_list.frame$index {top expand fill } 
}


proc select_anchor {index} {
    .meshbuttons.obj_list.frame$index.anchor select
}


proc select_align {index} {
    .meshbuttons.obj_list.frame$index.align select
}


proc select_target {index} {
    .meshbuttons.obj_list.frame$index.target select
}


proc select_merge {index} {
    .meshbuttons.obj_list.frame$index.merge select
}



proc rename_mesh_buttons {w name index} \
{
  $w.obj_list.frame$index.label configure -text $name
  $w.obj_list.frame$index.vis configure \
    -command "mesh_draw $name \$mesh_draw_$index"
  $w.obj_list.frame$index.radio configure \
    -value $name -command "move_object $name"
  $w.obj_list.frame$index.anchor configure \
    -value $name -command "anchor $name"
  $w.obj_list.frame$index.align configure \
    -value $name -command "next_align $name"
  $w.obj_list.frame$index.target configure \
    -value $name -command "target $name"
  $w.obj_list.frame$index.merge configure \
    -value $name -command "next_merge $name"
}

bind . <Escape> { destroy . }

proc eat {mesh1 mesh2} {
  set num 1
  while {$num != 0} {
    set num [eat_step $mesh1 $mesh2]
  }
}

# create a window for editing the triangle mesh
set w .editbuttons
toplevel $w
wm title $w "edit"
wm withdraw .editbuttons

set bset {top pady 4 frame w}

frame $w.frame1 -borderwidth 10
pack append $w.frame1 \
    [radiobutton $w.frame1.b1 -relief flat -text "Null Mode" \
     -variable edit_state -value 0 \
     -command set_edit ] $bset \
    [radiobutton $w.frame1.b2 -relief flat -text "Delete Triangle" \
     -variable edit_state -value 1 \
     -command set_edit ] $bset \
    [radiobutton $w.frame1.b3 -relief flat -text "Delete Vertex" \
     -variable edit_state -value 2 \
     -command set_edit ] $bset \
    [radiobutton $w.frame1.b4 -relief flat -text "Fill Hole" \
     -variable edit_state -value 3 \
     -command set_edit ] $bset \
    [radiobutton $w.frame1.b5 -relief flat -text "Fill Better" \
     -variable edit_state -value 4 \
     -command set_edit ] $bset \
    [radiobutton $w.frame1.b6 -relief flat -text "Make Triangle" \
     -variable edit_state -value 5 \
     -command set_edit ] $bset \
    [radiobutton $w.frame1.b7 -relief flat -text "Delete Box" \
     -variable edit_state -value 6 \
     -command set_edit ] $bset \
    [button $w.frame1.b8 -text "Execute Edit" -command execute_edit] $bset \
    [button $w.frame1.b9 -text "Dismiss" -command edit_window] $bset

set edit_state 0

pack append $w \
  $w.frame1 {left expand fill}

set button_list {
  .topbuttons.frame1.b1
  .topbuttons.frame1.b2
  .topbuttons.frame1.b3
  .topbuttons.frame1.b4

  .topbuttons.frame2.b1
  .topbuttons.frame2.b2
  .topbuttons.frame2.b3

  .topbuttons.frame3.b1
  .topbuttons.frame3.b2
  .topbuttons.frame3.b3
  .topbuttons.frame3.b4

  .topbuttons.frame4.b1
  .topbuttons.frame4.b2
  .topbuttons.frame4.b3
  .topbuttons.frame4.b4

  .topbuttons.frame5.b1
  .topbuttons.frame5.b2
  .topbuttons.frame5.b3
  .topbuttons.frame5.b4
}

proc disable_all {button_list} {
  foreach i $button_list { $i configure -state disabled }
  .topbuttons.frame3.b2 configure -state normal
}

proc enable_all {button_list} {
  foreach i $button_list { $i configure -state normal }
  .topbuttons.frame3.b2 configure -state disabled
}


proc set_shiny {shiny_state} {
    if {$shiny_state} {
	material 0.1 0.6 0.4 30
    } else {
	material 0.1 0.8 0 0
    }
}


proc align_pair {mesh1 mesh2 numiters start_level target_level} {
    anchor $mesh1
    for {set level $start_level} {$level >= $target_level} \
	    {set level [expr $level-1]} {
	mesh_resolution $level
	puts ""
	puts "Aligning $mesh2 against $mesh1 at level $level in $numiters passes"
	for {set i 0} {$i < $numiters} {incr i} {
	    align $mesh2
	}
    }
}


proc align_step_pair {mesh1 mesh2 numiters start_level target_level} {
    anchor $mesh1
    for {set level $start_level} {$level >= $target_level} \
	    {set level [expr $level-1]} {
	mesh_resolution $level
	puts ""
	puts "Aligning $mesh2 against $mesh1 at level $level in $numiters passes"
	for {set i 0} {$i < $numiters} {incr i} {
	    align_step $mesh2
	}
    }
}


proc align_conf {confile iters start_level stop_level} {
    set fileid [open $confile "r"]
    set count 0
    set numchars [gets $fileid line($count)]
    incr count
    while {$numchars > 0} {
	set numchars [gets $fileid line($count)]
	incr count
    }
    close $fileid
    
    set count [expr $count -1]
    
    for {set index 0} {$index < $count} {incr index} {
	set curline $line($index)
	if {"bmesh" == [lindex $curline 0]} {
	    set anchor_mesh [lindex $curline 1]
	    break
	}
    }
    
    set first $index
    set second [expr $first + 1]
    
    for {set index $second} {$index < $count} {incr index} {
	set curline $line($index)
	if {"bmesh" == [lindex $curline 0]} {
	    set align_mesh [lindex $curline 1]
	    align_pair $anchor_mesh $align_mesh $iters $start_level $stop_level
	    set anchor_mesh $align_mesh
	    print $confile
	}
    }
}


proc align_conf_to_first {confile iters start_level stop_level} {
    set fileid [open $confile "r"]
    set count 0
    set numchars [gets $fileid line($count)]
    incr count
    while {$numchars > 0} {
	set numchars [gets $fileid line($count)]
	incr count
    }
    close $fileid
    
    set count [expr $count -1]
    
    for {set index 0} {$index < $count} {incr index} {
	set curline $line($index)
	if {"bmesh" == [lindex $curline 0]} {
	    set anchor_mesh [lindex $curline 1]
	    break
	}
    }
    
    set first $index
    set second [expr $first + 1]
    
    for {set index $second} {$index < $count} {incr index} {
	set curline $line($index)
	if {"bmesh" == [lindex $curline 0]} {
	    set align_mesh [lindex $curline 1]
	    align_pair $anchor_mesh $align_mesh $iters $start_level $stop_level
	    print $confile
	}
    }
}


proc merge_conf {confile zipper_file use_first} {
    set fileid [open $confile "r"]
    set count 0
    set numchars [gets $fileid line($count)]
    incr count
    while {$numchars > 0} {
	set numchars [gets $fileid line($count)]
	incr count
    }
    close $fileid
    
    set count [expr $count -1]
    
    if {$use_first == 1} {
	for {set index 0} {$index < $count} {incr index} {
	    set curline $line($index)
	    if {"bmesh" == [lindex $curline 0]} {
		set target_mesh [lindex $curline 1]
		break
	    }
	}
    } else {
	set target_mesh $zipper_file
    }
    
    target $target_mesh
    puts "Using $target_mesh as target mesh."
    
    for {set index 0} {$index < $count} {incr index} {
	set curline $line($index)
	if {"bmesh" == [lindex $curline 0]} {
	    set mesh_name [lindex $curline 1]
	    if {$mesh_name != $target_mesh} {
		puts "Merging $mesh_name into $target_mesh..."
		merge $mesh_name
		puts "Writing $target_mesh to $zipper_file..."
		bpwrite $target_mesh $zipper_file
	    }
	}
    }
}


proc precon_merge_conf {confile zipper_file use_first resolution} {
    mesh_resolution $resolution
    set fileid [open $confile "r"]
    set count 0
    set numchars [gets $fileid line($count)]
    incr count
    while {$numchars > 0} {
	set numchars [gets $fileid line($count)]
	incr count
    }
    close $fileid
    
    set count [expr $count -1]
    
    if {$use_first == 1} {
	for {set index 0} {$index < $count} {incr index} {
	    set curline $line($index)
	    if {"bmesh" == [lindex $curline 0]} {
		set target_mesh [lindex $curline 1]
		break
	    }
	}
    } else {
	set target_mesh $zipper_file
    }
    
    for {set index 0} {$index < $count} {incr index} {
	set curline $line($index)
	if {"bmesh" == [lindex $curline 0]} {
	    set mesh_name [lindex $curline 1]
	    puts "Performing pre-consensus pass on $mesh_name"	    
	    consensus $mesh_name $resolution
	}
    }

    target $target_mesh
    puts "Using $target_mesh as target mesh."
    
    for {set index 0} {$index < $count} {incr index} {
	set curline $line($index)
	if {"bmesh" == [lindex $curline 0]} {
	    set mesh_name [lindex $curline 1]
	    if {$mesh_name != $target_mesh} {
#		puts "Mutual consensus: $mesh_name and $target_mesh..."
#		consensus $target_mesh $resolution $target_mesh $mesh_name
#		consensus $mesh_name $resolution $target_mesh $mesh_name
		puts "Merging $mesh_name into $target_mesh..."
		merge $mesh_name
		puts "Writing $target_mesh to $zipper_file..."
		bpwrite $target_mesh $zipper_file
	    }
	}
    }
}


proc pairwise_precon_merge_conf {confile zipper_file use_first resolution} {
    mesh_resolution $resolution
    set fileid [open $confile "r"]
    set count 0
    set numchars [gets $fileid line($count)]
    incr count
    while {$numchars > 0} {
	set numchars [gets $fileid line($count)]
	incr count
    }
    close $fileid
    
    set count [expr $count -1]
    
    if {$use_first == 1} {
	for {set index 0} {$index < $count} {incr index} {
	    set curline $line($index)
	    if {"bmesh" == [lindex $curline 0]} {
		set target_mesh [lindex $curline 1]
		break
	    }
	}
    } else {
	set target_mesh $zipper_file
    }
    
    target $target_mesh
    puts "Using $target_mesh as target mesh."
    
    for {set index 0} {$index < $count} {incr index} {
	set curline $line($index)
	if {"bmesh" == [lindex $curline 0]} {
	    set mesh_name [lindex $curline 1]
	    if {$mesh_name != $target_mesh} {
		puts "Mutual consensus: $mesh_name and $target_mesh..."
		pair_consensus $target_mesh $mesh_name $resolution
#		consensus $target_mesh $resolution $target_mesh $mesh_name 
		puts "Merging $mesh_name into $target_mesh..."
		merge $mesh_name
		puts "Writing $target_mesh to $zipper_file..."
		bpwrite $target_mesh $zipper_file
	    }
	}
    }
}


proc cons_conf {zipper_mesh conf_files cons_mesh} {
    bmesh $zipper_mesh
    set num_files [llength $conf_files]
    for {set i 0} {$i < $num_files} {incr i} {
	set confile [lindex $conf_files $i]
	set fileid [open $confile "r"]
	set count 0
	set numchars [gets $fileid line($count)]
	incr count
	while {$numchars > 0} {
	    set numchars [gets $fileid line($count)]
	    incr count
	}
	close $fileid
	
	set count [expr $count -1]
    
	for {set index 0} {$index < $count} {incr index} {
	    set curline $line($index)
	    if {"bmesh" == [lindex $curline 0]} {
		set mesh_name [lindex $curline 1]
		if {$mesh_name != $zipper_mesh} {
		    eval $curline
		}
	    }
	}
    }
    consensus $zipper_mesh 1
    bpwrite $zipper_mesh $cons_mesh
}

# countProcessors
#
# Brian Curless
# Added by Afra Zomorodian 7/18/95
# counts the processors on the machine
proc countProcessors {} {
    catch {exec hinv -c processor} msg
    scan $msg "%s" numProcessors
    return $numProcessors
}
