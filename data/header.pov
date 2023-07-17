global_settings {
  max_trace_level 10
  assumed_gamma 2.2
}

camera {
  orthographic
  location <0.0000, 0.0000, -2.0000>
  look_at <-0.0000, -0.0000, 2.0000>
  right x*image_width/image_height
}
light_source { 
  <-0.1000, 0.1000, -1.0000> 
  color rgb<1.000, 1.000, 1.000> 
  parallel 
  point_at <0.0, 0.0, 0.0> 
}
light_source { 
  <1.0000, 2.0000, -0.5000> 
  color rgb<1.000, 1.000, 1.000> 
  parallel 
  point_at <0.0, 0.0, 0.0> 
}

#default { texture {
 finish { ambient 0.000 diffuse 0.620 phong 0.1 phong_size 38.905 specular 0.000 }
} }

