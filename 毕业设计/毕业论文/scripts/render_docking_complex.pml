# Figure 3-3: P450-substrate docking complex (overview + zoom-in)
# Source: PDB 2HI4 (CYP1A2 + alpha-naphthoflavone)
# Outputs: docking_overview.png, docking_zoom.png

reinitialize
bg_color white

load 2HI4.pdb, complex

remove resn HOH
remove (hydro)
hide everything

select protein, polymer.protein
select heme, resn HEM
select ligand, resn BHF
select fe_atom, resn HEM and name FE
select pocket, byres (protein within 5.0 of ligand)
select pocket_label, pocket and name CA

# ---------- Common rendering settings ----------
set ray_shadows, 0
set specular, 0.3
set ambient, 0.4
set ray_trace_mode, 1
set ray_trace_color, grey40
set antialias, 2
set ray_opaque_background, on
set ray_trace_fog, 0
set depth_cue, 0
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1

# ---------- Heme styling (used in both views) ----------
show sticks, heme
color grey60, heme and elem C
color blue, heme and elem N
color red, heme and elem O
util.cnc heme

show spheres, fe_atom
color orange, fe_atom
set sphere_scale, 0.5, fe_atom

# ---------- Ligand styling ----------
show sticks, ligand
color forest, ligand and elem C
color blue, ligand and elem N
color red, ligand and elem O
set stick_radius, 0.22, ligand

# ============== VIEW 1: OVERVIEW ==============
show cartoon, protein
color lightblue, protein
set cartoon_transparency, 0.0, protein
hide sticks, pocket and not heme and not ligand

orient protein
zoom protein, 2

ray 1800, 1800
png docking_overview.png, dpi=300

# ============== VIEW 2: ZOOM-IN ==============
# Pocket residues as thin sticks
show sticks, pocket and (sidechain or name CA)
color cyan, pocket and elem C and (not resn HEM) and (not resn BHF)
util.cnc (pocket and not resn HEM and not resn BHF)
set stick_radius, 0.16, pocket and not ligand and not heme

# More transparent cartoon for context
set cartoon_transparency, 0.7, protein

# Distance labels Fe to ligand
distance fe_dist, fe_atom, ligand, 5.5
color black, fe_dist
set dash_width, 1.8
set dash_radius, 0.04
set dash_gap, 0.25
set label_size, 16
set label_color, black
set label_font_id, 7
set label_outline_color, white

# Label key pocket residues (representative subset)
label (pocket_label and resi 226), "Phe226"
label (pocket_label and resi 305), "Asp320"
label (pocket_label and resi 322), "Phe260"
label (pocket_label and resi 124), "Phe125"

set label_position, (0, 1.8, 0)

# Camera: zoom on ligand-heme region
orient ligand or fe_atom
zoom (ligand or fe_atom), 5
turn x, -5
turn y, 10

ray 2000, 1800
png docking_zoom.png, dpi=300

quit
