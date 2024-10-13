////بِسْمِ ٱللَّٰهِ ٱلرَّحْمَٰنِ ٱلرَّحِيمِ//


 ///This is a simple ray tracing


          #ifdef GL_ES
        precision mediump float;
     #endif


   uniform vec2 u_resolution;
   uniform vec2 u_mouse;
   uniform float u_time;
   #define PI 3.14159265359
   #define TWO_PI 6.28318530718
   float c = .1;

   //bluenoise///

   float random (vec2 st) {
    return fract(sin(1.222*dot(st.xy,
                         vec2(1222.9898,72.233))+u_time)*
    111.5453123);
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


  




 #define MAX_OCTAVES 8
   







 /// LightRay ////
vec3 lightray1 = normalize(vec3(0.2745, 02.2275, 03.2275));
   vec3 pointlight2 = vec3(-0.9118, -.0, -0.508);
     vec3 pointlight3 = vec3(-0.9118, -.3,1.508);

///2D Rotational Matrix////

mat2 rotate(float th) {

  return mat2(
    cos(th), -sin(th),
    sin(th), cos(th)
  );
}

///3D Rotation Matrix
mat3 rotate3D(float theta) {
  mat3 Angelx = mat3(1.,0.,0.,
              0., cos(theta), -sin(theta),
              0., sin(theta), cos(theta));
  mat3 Angely = mat3(cos(theta),0.,sin(theta),
              0., 1., 0.,
              -sin(theta), 0., cos(theta));
  mat3 Angelz = mat3(cos(theta),-sin(theta),0.,
              sin(theta), cos(theta), 0.,
              0., 0., 1.);
  
       return Angelz;
  
}
 ////WaveFunction////////////
   float wave(vec2 inpu) {
        vec3 inp = vec3(inpu,0.);
        float x = (inp.x+inp.y+inp.z)/3.;
        vec3 amplituds = vec3(-.4,1.,2.);
        vec3 frequency = vec3(2.8,.3,1.43);
         float phz = 7.4;   
        float sin1 = amplituds.x*sin((x*frequency.x)+u_time*phz);
        float sin2 = amplituds.y*sin((x*frequency.y)+u_time*phz);
        float sin3 = amplituds.z*sin((x*frequency.z)+u_time*5.);

        return (sin1-sin2+sin3)*.3-.3;
    }
        float sumwave(vec2 inpv) {
          vec3 positon = vec3(11.56,12.1,0.);
           vec3 position2 = vec3(1.16,0.,0.);
       
            float wave1 = wave(inpv*positon.xy);
             float wave2 = wave((inpv+position2.xy)*positon.xy);

               return wave1*wave2;
           }
          
         float finalwave(vec2 inpt) {
           vec3 timeposition = vec3(.0,u_time,0.);
              inpt = inpt + timeposition.xy;
              vec2 inpt2 = rotate(2.79253)*inpt;

              float wave1f = sumwave(inpt);
              float wave2f = sumwave(inpt2);


              return wave1f*wave2f;

         } 

void main () {
    
           vec2 st = gl_FragCoord.xy/u_resolution.xy;
                vec2 uv = 2.*(st-.5);
                
                 uv.x *= u_resolution.x/u_resolution.y;
                
               st.x *= u_resolution.x/u_resolution.y;
                vec2 uv1 = st;

              ///Pre Dethering effect////  Optional
                  st *= 111.0; // Scale the coordinate system by 10
                   vec2 ipos = floor(st);  // get the integer coords
                    vec2 fpos = fract(st);  // get the fractional coords
                     float noise = random(ipos);
                     float noise11 = random1(fpos);
                     float noise2 = random2(fpos);
                      vec3 vecnoise = vec3(noise,noise11,noise2);
                         vecnoise = normalize(vecnoise);
                          //uv = uv- H*(vecnoise.xy);



               ///camera Direction
               float f = .6;
               vec3 ray1 = (normalize(vec3(uv.x,uv.y,f)));
                vec2 ray = ray1.xy;
                 vec2 mouse = u_mouse.xy;



                 float theta = -clamp((mouse.x/mouse.y),0.,3.)*c;
                 ray = rotate(theta)*ray;


     ////Scene //////////
        //sphere 01/////
         float sph1 = sphere(ray , r1,pos1);
         bool m1 = (sph1 >0.);
         vec3 n1 = vec3(ray.x+pos1.x,sph1+pos1.z,ray.y+pos1.y) * vec3(m1);

                     //mouse interaction 3D
                      vec2 mouse3D = u_mouse;
                      float phi = -(mouse3D.x/mouse3D.y)+1.;
                       //float phi = u_time;
                        n1 = rotate3D(phi)*n1;
                           n1 = rotate3D(phi)*n1;


                       

                 ///Wave scene
                  
                 float d = 0.0;
                   // Remap the space to -1. to 1.
                   //uv1 = uv1- H*(vecnoise.xy);
                  uv1 = 2.*(uv1-.5);  float th = u_time;
                  vec2 uv2 = (uv1 * vec2(1.,4.));
                     uv2 = uv2*.6;
                  ///postion !
                   uv2 += vec2(-.4,-2.);



                   uv2 = rotate(-phi+1.) * uv2;
                       uv2 = rotate(3.14/2.)*uv2;
                     // Number of sides of your shape
                    int N = 4;

                    // Angle and radius from the current pixel
                   float a = atan(uv2.x,uv2.y)+PI;
                   float r = TWO_PI/float(N);

                  // Shaping function that modulate the distance
                  d = cos(floor(.5+a/r)*r-a)*length(uv2);

                  vec3 color = vec3(1.0-smoothstep(.4,.41,d));
     
                  vec2 squv = uv2*color.x;
                        
                      vec3 sqn = vec3(0.,.0,1.)*color.x;

                       vec3 scenenormal = n1; 
                 
                    
                    ///Wave noise implement:

                      float wvf = finalwave(squv)*color.x;
                      vec3 wvn = vec3(squv.x*wvf,squv.y*wvf,1.-abs(wvf))*color;
                    //////////////////     

                        if (color.x >0.) {
                        //scenenormal = wvn.xzy;

                       } 

                    ///Lets calculate the refraction's 

                                  /// uv sphere : 
                          
                                   vec2 spuv = vec2(-n1.xz);
                                     squv = rotate(2.)*squv;
                                   ///sphere caustic
                                     float wvsq = finalwave(spuv);
                                     vec3 wvsn = vec3(spuv.x*wvsq,spuv.y*wvsq,1.-abs(wvsq));
                   
                                      vec3 spref = refract(normalize(vec3(1, 0., 1.9)),normalize(wvsn),1.3333);
                                        vec3 sprr = vec3(2.*spref.y);
                     ////lets Make a Metric Equation : 
                                    /// first metric equation:
                                    float wvsq2 = finalwave(spuv);
                                     vec3 wvsn2 = vec3(spuv.x*wvsq2,spuv.y*wvsq2,1.-abs(wvsq2));
                   
                                         vec3 spref2 = refract(normalize(vec3(0.,0.,1.)),normalize(wvsn2),1.3333);
                                          spref2 = rotate3D(phi)*spref2;
                                       // vec3 sprr2 = vec3(2.*spref2.z);
                                                  float sprr2 = dot(spref2, lightray1);  
 
                                             float wvsq3 = finalwave(spuv);
                                             vec3 wvsn3 = vec3(spuv.x*wvsq3,spuv.y*wvsq3,1.-abs(wvsq3));
                                      ///Second row matrix
                                                  vec3 spref3 = refract(normalize(vec3(0.4,0.8,1.)),normalize(wvsn3),1.3333);
                                                   // vec3 sprr2 = vec3(2.*spref2.z);
                                                  float sprr3 = dot(spref3, lightray1); 
                                                  // MAT 2  
                                                  vec3 ispref3 = refract(normalize(vec3(-0.6,0.2,1.)),normalize(wvsn3),1.3333);
                                                 // vec3 sprr2 = vec3(2.*spref2.z);
                                                  float isprr3 = dot(ispref3, lightray1);

                                      ///First row matrix
                                                  vec3 spref4 = refract(normalize(vec3(0.2,0.1,1.)),normalize(wvsn3),1.3333);
                                                   // vec3 sprr2 = vec3(2.*spref2.z);
                                                  float sprr4 = dot(spref4, lightray1); 
                                                  // MAT 2  
                                                  vec3 ispref4 = refract(normalize(vec3(-0.1,-0.3,1.)),normalize(wvsn3),1.3333);
                                                 // vec3 sprr2 = vec3(2.*spref2.z);
                                                  float isprr4 = dot(ispref4, lightray1); 


                                                  mat2 caustic = mat2(sprr3,isprr3,
                                                                    sprr4 ,isprr4);         

                                                                    ///evaluate the matrix to the metric equation : 
                                                              float invsum = (sprr3+isprr3)*(sprr4+isprr4);
                                                            float sumcaustic = ((sprr3+isprr3)+(sprr4+isprr4))/4.;                                                 
                                                               // sprr = sprr+sumcaustic;      

           ///light ray    
            float lightmap1 = clamp((dot(scenenormal , lightray1)*1.),0.,1.);
              //point light
               vec3 HitPoint = vec3(ray1)*vec3(1.-(sph1+0.0000001));
                pointlight2 = normalize(pointlight2-HitPoint);  
                 float lightmap2 = clamp((dot(scenenormal,pointlight2)*1.),.0,1.);
                  pointlight3 = normalize(pointlight3-HitPoint);
                    float lightmap3 = clamp((dot(scenenormal,pointlight3)*1.),.0,1.);



                    vec3 light1 = vec3(0.7647, 0.6549, 0.4863)*vec3(lightmap1);
                    vec3 light2 = vec3(0.1176, 0.1137, 0.1176)*vec3(lightmap2);
                    
                    vec3 light3 = vec3(0.0902, 0.1725, 0.2157)*vec3(lightmap3);

                    light1 = light1+clamp(((sprr*light1)-color),.0,1.);
                    vec3 lightmap = light2+light1+light3;                
                  


                  ///Post Proccessing////
                 
                
                 float flimgain = noise-.5;
                   //lightmap = lightmap+flimgain;




                   lightmap = aces(lightmap);

    gl_FragColor = vec4( vec3(lightmap),1.);
}

