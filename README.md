# RayTracing
simple RayTracing

support Sphere,Plane,Cuboid 

all the classes:

vec : simple 3D vector

    vec.x,vec.y,vec.z:3D coordinate
  
    vec.len: vec's length
  
    vec.slen: vec.len*vec.len
  
color : color=vec

Ray : the light ray

    Ray.st: the start coordinate of the ray
  
    Ray.dir: the direction
  
shape : the fathe class of all props in the environment

    kA : ambient Reflection Coefficient [0,1]
    
    kD : diffuse Reflection Coefficient [0,1]
     
    kS : specular Reflection Coefficient [0,1]
    
    kT : refraction Reflection Coeeficient [0,1]
    
    n : refractive index
    
circle : the sphere class

    center: the center of the sphere
    
    dim : diameter
    
plane : the plane

     n : the normal vector of the plane
   
     dis : the distance from the plane to the (0,0,0)
    
rectangle: the rectangle

    a,b,c,d : the four vertex of the rectangle
   
    n : the normal vector 
   
    dis : same as plane
    
cuboid : the cuboid

    r[6] : the six rectangle of the cuboid
    
camera : the camera

    eye : the position of the camera
   
    forward : the forward vector(normalized)  
   
    up : the up vector(normalized)
   
    right : the right vector(normalized)
   
    stdUp : used for the first building to find the up and
    right vector
   
    fov : the field of angle 
    
    scale : the length of the field
    
    light : the point light in the environment

    pos : the position of the light
    
    color : the color of the light
