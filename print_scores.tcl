# Carregar a trajetória
mol addfile "traj.pdb" type pdb waitfor all

set num_frame [molinfo top get numframes]

# Função para obter S (score)
proc score { d1 d2 phi psi } {

   ## Evita valores de ângulo negativos (mas não é preciso nesse caso)
   #if {$phi < 0} {
   #  set phi_graus [expr {$phi + 360}]
   #} else {
   #  set phi_graus $phi
   #}
  
   set pi [expr {acos(-1)}]

   set phi_rad [expr {$phi * $pi / 180}]
   set psi_rad [expr {$psi * $pi / 180}]
 
   # puts "phi: $phi"
   # puts "phi_graus: $phi_graus"
   # puts "phi_rad: $phi_rad"

   if { $phi_rad >= 2.0 } {
     if { $psi_rad >= 1.5 && $psi_rad <= 2.0 } {  
       set score [expr { 10 * ( (1/$d1) + (1/$d2) + (5.0 * $psi_rad) + (5.0 * $phi_rad) ) }]
     } else {
       set score [expr { 10 * ( (1/$d1) + (1/$d2) + $psi_rad + (5.0 * $phi_rad) ) }]
     }  	     
   } else {
     set score [expr { 10 * ( (1/$d1) + (1/$d2) + $psi_rad + $phi_rad ) }]
   }

   return $score
}

# Função para calcular a distância entre dois pontos
proc distance {c1 c2} {
   set x_diff [expr {[lindex $c2 0] - [lindex $c1 0]}]
   set y_diff [expr {[lindex $c2 1] - [lindex $c1 1]}]
   set z_diff [expr {[lindex $c2 2] - [lindex $c1 2]}]

   set distance [expr {sqrt($x_diff**2 + $y_diff**2 + $z_diff**2)}]
   return $distance
}

# Função para subtrair vetores
proc vecsub {vec1 vec2} {
    return [list [expr {[lindex $vec1 0] - [lindex $vec2 0]}] [expr {[lindex $vec1 1] - [lindex $vec2 1]}] [expr {[lindex $vec1 2] - [lindex $vec2 2]}]]
}

# Função para calcular o vetor entre dois átomos
proc calculate_vector {coord1 coord2} {
    return [vecsub $coord1 $coord2]
}

# Função para calcular o produto cruzado de dois vetores
proc veccrossprod {u v} {
    set x [expr {[lindex $u 1]*[lindex $v 2] - [lindex $u 2]*[lindex $v 1]}]
    set y [expr {[lindex $u 2]*[lindex $v 0] - [lindex $u 0]*[lindex $v 2]}]
    set z [expr {[lindex $u 0]*[lindex $v 1] - [lindex $u 1]*[lindex $v 0]}]
    return [list $x $y $z]
}

# Função para calcular o vetor normal a partir de três pontos
proc normal_vector {p1 p2 p3} {
    set u [vecsub $p2 $p1]
    set v [vecsub $p3 $p1]
    return [veccrossprod $u $v]
}

# Função para calcular o produto escalar de dois vetores
proc vecdot {vec1 vec2} {
    set dot 0
    for {set i 0} {$i < [llength $vec1]} {incr i} {
        set dot [expr {$dot + [lindex $vec1 $i] * [lindex $vec2 $i]}]
    }
    return $dot
}

# Função para calcular a norma de um vetor
proc vecnorm {vec} {
    return [expr {sqrt([vecdot $vec $vec])}]
}

# Função para calcular o ângulo entre dois vetores
proc angle_between_vectors {vec1 vec2} {
    set dot [vecdot $vec1 $vec2]
    set norm1 [vecnorm $vec1]
    set norm2 [vecnorm $vec2]
    
    if {$norm1 == 0 || $norm2 == 0} {
        return -1 ; # Indica que não foi possível calcular o ângulo
    }

    set cos_angle [expr {$dot / ($norm1 * $norm2)}]
    
    # Garantir que o valor está no intervalo [-1, 1]
    if {$cos_angle > 1} { set cos_angle 1 }
    if {$cos_angle < -1} { set cos_angle -1 }
    
    return [expr {acos($cos_angle) * 180 / [expr {acos(-1)}]}]
}


# Função para calcular o ângulo de rotação entre dois vetores
proc rotation_angle_between_vectors {vec1 vec2} {
    set cross_prod [veccrossprod $vec1 $vec2]
    set norm_cross [vecnorm $cross_prod]
    set angle [angle_between_vectors $vec1 $vec2]

    # O ângulo de rotação é calculado considerando o ângulo e a magnitude do vetor perpendicular
    # Se a magnitude do vetor perpendicular for zero, a orientação não pode ser determinada
    if {$norm_cross == 0} {
        return -1 ;# Indica que não foi possível calcular a orientação
    }

    # return [list $angle $norm_cross]
    return $norm_cross
}

set out_file1 [open "scores.csv" w]
set out_file2 [open "scores_cutoff.csv" w]


for {set frame 0} { $frame < $num_frame } {incr frame} {

  animate goto $frame

  puts $frame
  
  set bgl1_C7 [atomselect top "resname BGL and resid 1 and name C7"]
  set bgl2_C7 [atomselect top "resname BGL and resid 2 and name C7"]

  set bgl1_C11 [atomselect top "resname BGL and resid 1 and name C11"]
  set bgl2_C11 [atomselect top "resname BGL and resid 2 and name C11"]

  set bgl1_C7_index [$bgl1_C7 get index]
  set bgl2_C7_index [$bgl2_C7 get index]

  set bgl1_C11_index [$bgl1_C11 get index]
  set bgl2_C11_index [$bgl2_C11 get index]
  
  set d1_1 [measure bond "$bgl1_C7_index $bgl2_C11_index"]
  set d1_2 [measure bond "$bgl2_C7_index $bgl1_C11_index"]

  if { $d1_1 <= $d1_2 } {
    set d1 $d1_1
  } else {
    set d1 $d1_2
  }

  set sel1 [atomselect 0 "resname BGL"] 
  set sel2 [atomselect 0 "protein and resid 166 377"] 
  set c1 [measure center $sel1 weight [$sel1 get mass]]
  set c2 [measure center $sel2 weight [$sel2 get mass]]
  set d2 [distance $c1 $c2]

  set ring1 [atomselect top "resname BGL and resid 1 and name C7 C8 C9 C10 C11 O1"]
  set ring2 [atomselect top "resname BGL and resid 2 and name C7 C8 C9 C10 C11 O1"]

  # Obter as coordenadas dos átomos dos anéis
  set coords1 [$ring1 get {x y z}]
  set coords2 [$ring2 get {x y z}]

  # Calcular os vetores normais aos planos
  set normal1 [normal_vector [lindex $coords1 0] [lindex $coords1 1] [lindex $coords1 2]]
  set normal2 [normal_vector [lindex $coords2 0] [lindex $coords2 1] [lindex $coords2 2]]

  # Calcular o ângulo de inclinação
  set phi [angle_between_vectors $normal1 $normal2]

  puts "The angle of inclination between the planes is: $phi degrees"

  # Selecionar os átomos para obter os vetores C8 -> O1
  set bgl1_C8 [atomselect top "resname BGL and resid 1 and name C8"]
  set bgl1_O1 [atomselect top "resname BGL and resid 1 and name O1"]
  set bgl2_C8 [atomselect top "resname BGL and resid 2 and name C8"]
  set bgl2_O1 [atomselect top "resname BGL and resid 2 and name O1"]

  # Obter as coordenadas dos átomos
  set coords1_C8 [$bgl1_C8 get {x y z}]
  set coords1_O1 [$bgl1_O1 get {x y z}]
  set coords2_C8 [$bgl2_C8 get {x y z}]
  set coords2_O1 [$bgl2_O1 get {x y z}]

  # Calcular vetores
  set vector1 [calculate_vector [lindex $coords1_C8 0] [lindex $coords1_O1 0]]
  set vector2 [calculate_vector [lindex $coords2_C8 0] [lindex $coords2_O1 0]]

  # Calcular o ângulo entre os vetores
  set psi [angle_between_vectors $vector1 $vector2]

  puts "The angle between the vectors is: $psi degrees"

  set S [score $d1 $d2 $phi $psi]
  puts $S
 
  # Frame, distances (d1 and d2), phi, and S: 
  puts $out_file1 "$frame,$d1,$d2,$phi,$psi,$S"

  # set cutoff 100

  # if { $S >= $cutoff } {
  #  puts $out_file2 "$frame,$d1,$d2,$phi,$psi,$S"
  # } else {
  #  continue
  # }

}

close $out_file1
close $out_file2

# Clear memory
mol delete 0

exit

