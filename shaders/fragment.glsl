////بِسْمِ ٱللَّٰهِ ٱلرَّحْمَٰنِ ٱلرَّحِيمِ//


 ///This is a simple ray tracing


          #ifdef GL_ES
        precision mediump float;
     #endif


   uniform vec2 u_resolution;
   uniform vec2 u_mouse;
   uniform float u_time;
   
   float c = .1;

   //bluenoise///

   float random (vec2 st) {
    return fract(sin(1.*dot(st.xy,
                         vec2(1222.9898,72.233))+2.*u_time)*
    2222.5453123);
}


   float random1 (vec2 st) {
    return fract(sin(1.2*dot(st.xy,
                         vec2(22222.9898,72.233))+2.*u_time)*
    2222.5453123);
}



   float random2 (vec2 st) {
    return fract(sin(1.333*dot(st.xy,
                         vec2(122.9898,72.233))+2.*u_time)*
    2222.5453123);
}


   /// ACES tone mapping 

   
   vec3 aces(vec3 x) {
  const float a = 2.51;
  const float b = 0.03;
  const float c = 2.43;
  const float d = 0.59;
  const float e = 0.14;
  return clamp((x * (a * x + b)) / (x * (c * x + d) + e), 0.0, 1.0);
}


   float H = .046778;

    float r1 = .6; vec3 pos1 = vec3(0.0, 0.0, 0.0);
    ////////sphere///////////////////
    float sphere(vec2 st, float r, vec3 pos){

        /// there UV needs to be from -1....1 
          float h = pow(length(vec3(st,0.)+pos),2.);  
                 h = pow(r,2.) - h;
                 h = sqrt(h);
                     
      return h;
    }

 ////WaveFunction////////////





 #define MAX_OCTAVES 8








 /// LightRay ////
vec3 lightray1 = normalize(vec3(0.2275, 0.2275, 0.2275));
   vec3 pointlight2 = vec3(-0.9118, -.0, -0.508);
     vec3 pointlight3 = vec3(-0.9118, -.3,1.508);

///Rotational Function////

mat2 rotate(float th) {

  return mat2(
    cos(th), -sin(th),
    sin(th), cos(th)
  );
}




void main () {
    
           vec2 st = gl_FragCoord.xy/u_resolution.xy;
                vec2 uv = 2.*(st-.5);
                 uv.x *= u_resolution.x/u_resolution.y;
               st.x *= u_resolution.x/u_resolution.y;


              ///Pre Dethering effect////  Optional
                  st *= 1234.0; // Scale the coordinate system by 10
                   vec2 ipos = floor(st);  // get the integer coords
                    vec2 fpos = fract(st);  // get the fractional coords
                     float noise = random(ipos);
                     float noise11 = random1(ipos);
                     float noise2 = random2(ipos);
                      vec3 vecnoise = vec3(noise,noise11,noise2);
                         vecnoise = normalize(vecnoise);
                          // uv = uv- H*(vecnoise.xy);



               ///camera Direction
               float f = .6;
               vec3 ray1 = (normalize(vec3(uv.x,uv.y,f)));
                vec2 ray = ray1.xy;
                 vec2 mouse = u_mouse.xy;



                 float theta = -clamp((mouse.x/mouse.y),0.,3.)*c;
                 ray = rotate(theta)*ray;


     ////Scene 
        //sphere 01/////
         float sph1 = sphere(ray , r1,pos1);
         bool m1 = (sph1 >0.);
         vec3 n1 = vec3(ray.x+pos1.x,sph1+pos1.z,ray.y+pos1.y) * vec3(m1);
           ///light ray    
            float lightmap1 = clamp((dot(n1 , lightray1)*1.),0.,1.);
              //point light
               vec3 HitPoint = vec3(ray1)*vec3(1.-(sph1+0.0000001));
                pointlight2 = normalize(pointlight2-HitPoint);  
                 float lightmap2 = clamp((dot(n1,pointlight2)*1.),.0,1.);
                  pointlight3 = normalize(pointlight3-HitPoint);
                    float lightmap3 = clamp((dot(n1,pointlight3)*1.),.0,1.);



                    vec3 light1 = vec3(1.3)*vec3(lightmap1);
                    vec3 light2 = vec3(0.2745, 0.1137, 0.1137)*vec3(lightmap2);
                    vec3 light3 = vec3(0.0235, 0.0, 0.2471)*vec3(lightmap3);
                    vec3 lightmap = light2+light1+light3;
                       

                 ///Wave function 


                
 
       




                 lightmap = aces(lightmap);

    gl_FragColor = vec4( vec3(lightmap),1.);
}

