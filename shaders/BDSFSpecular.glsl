
////بِسْمِ ٱللَّٰهِ ٱلرَّحْمَٰنِ ٱلرَّحِيمِ//



   ////Arttibutes////
  

            #ifdef GL_ES
        precision mediump float;
     #endif
   
    


   #define PI 3.14159265359
   #define TWO_PI 6.28318530718
   #define EYEPATHLENGTH 1
   #define FOV 66
    #define SAMPLES 64 

   //Global Varriables
   float c = .1;
   const float limit = float(EYEPATHLENGTH)*1000.;
   ////Uniform////

   uniform vec2 u_resolution;
   uniform vec2 u_mouse;
   uniform float u_time;
        
// Filmic Tonemapping Operators http://filmicworlds.com/blog/filmic-tonemapping-operators/
vec3 filmic(vec3 x) {
  vec3 X = max(vec3(0.0), x - 0.004);
  vec3 result = (X * (6.2 * X + 0.5)) / (X * (6.2 * X + 1.7) + 0.06);
  return pow(result, vec3(2.2));
}

vec3 aces(vec3 x) {
  const float a = 2.51;
  const float b = 0.03;
  const float c = 2.43;
  const float d = 0.59;
  const float e = 0.14;
  return clamp((x * (a * x + b)) / (x * (c * x + d) + e), 0.0, 1.0);
}



      ///uniformic Function
        ///hitPoint Function////
           vec3 hitpoint(vec3 ray, vec3 rayo, float dist) {
              return rayo+(ray*vec3(dist));
           }
           
           ///vector to light function////

            vec4 vectorlight(vec3 pos, vec3 hp) {
            float len = length(pos-hp);
               return vec4((normalize(pos-hp)),len);
            }

         ///Seed Function;
       float hash1(inout float seed) {
    return fract(sin(seed += 0.1)*43758.5453123);
         }

    vec2 hash2(in float seed) {
    return fract(sin(vec2(seed+=0.1,seed+=0.1))*vec2(43758.5453123,22578.1459123));
}

            ///White Noise///
         float rand(vec2 st) {
    return fract(sin(dot(2.*st.xy,
                         vec2(12.9898,78.233))+u_time*33.)*
        4372238.5453123);
         }
      

       ////Unifrom 3D RGB Noise////
      vec3 Noise(in vec2 st){
          float x = rand(st.xy+u_time);
          float y = rand(st.yx);
          float z = rand(vec2(-st.x,st.y)-u_time);

          return vec3(x,y,z);
      }

//3D Rotational metric Equations


mat3 rotatex(float theta) {
  return mat3(1.,0.,0.,
              0., cos(theta), -sin(theta),
              0., sin(theta), cos(theta));
}
  mat3 rotatey(float theta) {
  return mat3(cos(theta),0.,sin(theta),
              0., 1., 0.,
              -sin(theta), 0., cos(theta));  }



mat2 rotate(float th) {

  return mat2(
    cos(th), -sin(th),
    sin(th), cos(th)
  );
}


      ////Global Function////
       ///LightFunction/// 
                   ///type Function////
         vec4 Pointlight(in vec3 rd,in vec3 ro,in vec3 normal, float d,in vec3 lightposition,in vec4 Lightcolor) {
                vec3 hp = hitpoint(rd,ro,d);  
                 vec3 glvectorlight = normalize(lightposition-hp);
                   float len = length(lightposition-hp);
                    float lm = max(dot(normal,glvectorlight),0.);
                    float glshade = Lightcolor.w*(lm/(pow(len,2.)));
                     


                     vec3 col = Lightcolor.xyz*vec3(glshade);
                     
            return vec4(vec3(col),lm);           ////the lightmap ---> .w // and LightColor ----> .rgb is Stored in a vec4 information: 
                              }  

                      ////General Direct LightSource////
                         vec4 DirLight(in vec3 normal, in vec3 lightpostion, in vec4 Lightcolor) {

                      vec3 gldirVector = normalize(lightpostion);
                    
                       float lm = Lightcolor.w*(dot(normal,gldirVector));
                      
                      vec3 col = Lightcolor.xyz*vec3(lm);
                      return vec4(vec3(col),lm);
                             }



    ///Sphere Functions////

     vec4 intsphere(in vec3 ray,in vec3 rayo,const vec3 pos,const float radius,in vec3 Normal,in float WSdist) {
         ///length intersection//
        float l = max((dot((pos-rayo),ray)),0.);
        float b = distance(pos,hitpoint(ray,rayo,l));  
         float f = float(bool(distance(rayo,pos)>radius));
         float i = float(bool(b < radius));    
             i = i*f;

          float a = l-(sqrt((radius*radius)-(b*b))); 
           
             float dist = a;
                if (float(i) < 1.){dist = limit;}
                  bool n = (dist<WSdist);
               vec3 normal = normalize(hitpoint(ray,rayo,dist)-pos);
                  if (float(n)<1.) {normal=Normal;}
                    
                      dist = min(dist,WSdist);  


       return vec4(vec3(normal),dist);
     }
       ///Global Generate Point on Hemisphere////
         vec3 coshemfun(in vec3 n,in vec2 st) {

           float r1 = hash1(st.x);
           float r2 = hash1(st.y); 
              ///unifrom distribution///
                float ud = sqrt(r2);  
              ///Shoot 2D point from r1 
              vec2 net1 = vec2(cos(TWO_PI*r1),sin(TWO_PI*r1));
                net1.x *= ud;
                net1.y *= ud;
              //pick hem///
               float net2 = sqrt(1.-r2);
                  vec3 hnet = vec3(net1.x,net1.y,net2);

            const vec3 in1 = vec3(-1.,-1.,1.);
            const vec3 in2 = vec3(0.,0.,1.);
             vec3 v1 = hnet*in1;
             vec3 v2 = n*in2;
             float dp = dot(v1,v2)/v2.z;
             vec3 normal = dp*n;

            return normal-v1;
         }
  
       ///Intersection Function
     ///Visual matrix
       vec3 pos1 = vec3(0.+sin(u_time), -0.6706, 3.2+cos(u_time)); float radi1 = .8;
       const vec3 pos2 = vec3(-0.3294, 0.60, 2.20745); float radi2 = 0.7;
        const vec3 pos3 = vec3(0.8,1.0,1.8); float radi3 = 0.8;
         const vec3 GlobalNormal = vec3(0.);



          
         
         vec4 intscene(in vec3 rd, in vec3 ro){
            vec4 isphere1 = intsphere(rd,ro,pos1,radi1,GlobalNormal,limit);
            vec4 isphere2 = intsphere(rd,ro,pos2,radi2,isphere1.xyz,isphere1.w);
            vec4 isphere3 = intsphere(rd,ro,pos3,radi3,isphere2.xyz,isphere2.w);
            vec3 Normal = isphere3.xyz;   
            float dist = isphere3.w; 



          return vec4(Normal,dist); 
         }     

         mat2 Mintid(in vec3 rd, in vec3 ro){
            vec4 isphere1 = intsphere(rd,ro,pos1,radi1,GlobalNormal,limit);
               bool i1 = (isphere1.w<limit);         
            vec4 isphere2 = intsphere(rd,ro,pos2,radi2,isphere1.xyz,isphere1.w);
               bool i2 = (isphere2.w<isphere1.w);
            vec4 isphere3 = intsphere(rd,ro,pos3,radi3,isphere2.xyz,isphere2.w);
               bool i3 = (isphere3.w<isphere2.w);
            vec3 Normal = isphere3.xyz;   
            float dist = isphere3.w; 



          return mat2(float(i1),float(i2),float(i3),0.); 
         }     


        
             ////Material Index////
        vec4 MatData(vec3 rd,vec3 ro) {
            vec4 index0 = vec4(vec3(0.1176, 0.1882, 0.5294),float (.0)); 
            vec4 index1 = vec4(vec3(1.0, 1.0, 1.0),float (.3)); 
            vec4 index2 = vec4(vec3(0.4196, 0.8667, 0.0549),float (.6)); 
             
               mat2 id = Mintid(rd,ro);
              
              vec3 c1 = index0.rgb * vec3(id[0][0]);
                vec3 c2 = index1.rgb * vec3(id[1][0]);
                 vec3 c3 = index2.rgb * vec3(id[0][1]);
                     
               vec3 col = c1+c2+c3;
                 if ( float(id[0][0]) < float(id[1][0])) {col = c2;}
                 if ( float(id[1][0]) < float(id[0][1])) {col = c3;}
               



            return vec4(col,1.);
        }

                   ///Light Source Matrix////
              
         const vec3 plightpos1 = vec3(0.8706, 02.6667, -0.6667); const vec4 icol1 = vec4(vec3(1.0, 1.0, 1.0),13.1);  
         const vec3 plightpos2 = vec3(0.2627, 0.0784, 1.22849); const vec4 icol2 = vec4(vec3(0.5294, 0.1333, 0.0549),0.1);   
         const vec3 plightpos3 = vec3(-1.2627, 01.0784, 0.849); const vec4 icol3 = vec4(vec3(0.2275, 0.1373, 0.2824),1.12);
         const vec3 plightpos4 = vec3(-2.2627, .0784, 0.5849); const vec4 icol4 = vec4(vec3(0.1608, 0.1608, 0.1608),.612);
      ///gldiffuse_function_Direct+indirect///
        vec4 gldiffuse(in vec3 rd, in vec3 ro,in vec3 normal,const float dist,in vec3 glhp,in vec2 st1){
             

            //Hitpoint
            vec3 hp = hitpoint(rd,ro,dist);
            vec3 hemdir1 = coshemfun(normal,st1);
               
                 vec4 gi = intscene(hemdir1,hp);
                 vec3 hp2 = hitpoint(hemdir1,hp,gi.w);



             /// Light Jitter Offsets 
               vec3 jit = 1.*Noise(st1*102.);
             ///vector to light
                   //Direct Light 
                vec4 ilight1 = Pointlight(rd,ro,normal,dist,plightpos1+jit,icol1);
                vec4 ilight2 = Pointlight(rd,ro,normal,dist,plightpos2+jit,icol2);
                vec4 ilight3 = Pointlight(rd,ro,normal,dist,plightpos3+jit,icol3); 
                vec4 ilight4 = Pointlight(rd,ro,normal,dist,plightpos4+jit,icol4); 
                      
                      ///indirect Light
                vec4 gilight1 = Pointlight(hemdir1,hp,gi.xyz,gi.w,plightpos1+jit,icol1);
                vec4 gilight2 = Pointlight(hemdir1,hp,gi.xyz,gi.w,plightpos2+jit,icol2);
                vec4 gilight3 = Pointlight(hemdir1,hp,gi.xyz,gi.w,plightpos3+jit,icol3); 
                vec4 gilight4 = Pointlight(hemdir1,hp,gi.xyz,gi.w,plightpos4+jit,icol4);

                     



                  vec4 vectorData1 = vectorlight(plightpos1+jit,glhp);
                  vec4 vectorData2 = vectorlight(plightpos2+jit,glhp);
                  vec4 vectorData3 = vectorlight(plightpos3+jit,glhp);
                  vec4 vectorData4 = vectorlight(plightpos4+jit,glhp);
                  //indirect data
                  vec4 giData1 = vectorlight(plightpos1+jit,hp2);
                  vec4 giData2 = vectorlight(plightpos2+jit,hp2);
                  vec4 giData3 = vectorlight(plightpos3+jit,hp2);
                  vec4 giData4 = vectorlight(plightpos4+jit,hp2);


                     vec4 iShadow1 = intscene(vectorData1.xyz,glhp);
                       bool cs1 = (iShadow1.w > vectorData1.w);
                        ilight1 = ilight1*vec4(cs1);
                      vec4 iShadow2 = intscene(vectorData2.xyz,glhp);
                        bool cs2 = (iShadow2.w > vectorData2.w);
                         ilight2 = ilight2*vec4(cs2);
                      vec4 iShadow3 = intscene(vectorData3.xyz,glhp);
                        bool cs3 = (iShadow3.w > vectorData3.w);
                         ilight3 = ilight3*vec4(cs3);
                      vec4 iShadow4 = intscene(vectorData4.xyz,glhp);
                        bool cs4 = (iShadow4.w > vectorData4.w);
                         ilight4 = ilight4*vec4(cs4);



                     vec4 giShadow1 = intscene(vectorData1.xyz,glhp);
                       float gcs1 = float(giShadow1.w > 3000.);
                       // gilight1 = gilight1*vec4(gcs1);
                      vec4 giShadow2 = intscene(vectorData2.xyz,glhp);
                        bool gcs2 = (giShadow2.w > 3000.);
                        // gilight2 = gilight2*vec4(gcs2);
                      vec4 giShadow3 = intscene(vectorData3.xyz,glhp);
                        bool gcs3 = (giShadow3.w > 3000.);
                        // gilight3 = gilight3*vec4(gcs3);
                      vec4 giShadow4 = intscene(vectorData4.xyz,glhp);
                        bool gcs4 = (giShadow4.w > 3000.);
                         //gilight4 = gilight4*vec4(gcs4);


              vec4 sumgili = gilight1+gilight2+gilight3+gilight4;
              vec4 SumLight = ilight1+ilight2+ilight3+ilight4;

              vec4 final = sumgili+SumLight;
            return final;
        }



    ///Reflections

    vec4 RTReflection(vec3 ro, vec3 hp,vec3 normal,vec3 rd, vec2 st, vec3 diff,in float r){


    float roughness = r;
        vec3 jit = 1.*Noise(st*102.);   
      
        vec3 ref = reflect(rd,normal+roughness*jit);
        vec4 refhit = intscene(ref,hp);
         vec3 hhp = hitpoint(rd,ro,refhit.w); 
             
           vec4 ColBuffer = MatData(rd,ro);

           vec4 refdiff = gldiffuse(rd,ro,refhit.xyz,refhit.w,hhp,st)*ColBuffer;

            ///Frenel Effect
            float a =.11;
            float f =  0.;
             if ( a > .001) {f = pow((dot(normal,rd)+1.),4.54);}
                 refdiff.xyz = refdiff.xyz*7.*diff;

            return vec4(refdiff.xyz,f);
    }




 ////Render Buffer Function/////


 vec3 rBuffer(vec3 RayD, vec3 RayO, vec2 st,vec2 seed) {



          vec4 tod; 
          vec3 col;
          vec3 glNormal;

           ///Scene 
          tod = intscene(RayD,RayO);
           glNormal = tod.xyz;    vec3 glHitPoint = hitpoint(RayD,RayO,tod.w);
        vec4 ColBuffer = MatData(RayD,RayO);  
         vec4 DiffBuffer = gldiffuse(RayD,RayO,glNormal,tod.w,glHitPoint,st+seed)*ColBuffer;
          vec4 RefBuffer = RTReflection(RayO,glHitPoint,glNormal,RayD,st+seed,DiffBuffer.xyz,float(0.2));  
            vec3 glrender = mix(DiffBuffer.rgb,RefBuffer.rgb,RefBuffer.w);

  return glrender;
 }




              
    void main() {
           
         ///THIS IS NOW 3D space Screen the idea is no longer using the SCREENuv we are using the object space to shootrays;
        


              vec2 st = gl_FragCoord.xy/u_resolution.xy;
              st = 2.*(st-.5);
              st.x *= u_resolution.x/u_resolution.y;    


         

          //RayConst 
          vec3 RayO = vec3(0.0, 0.0, -02.0);  

         
           //Fov to e 
         float e = 1./(tan((radians(float(FOV)))/2.));
         float ef = e; 
          
                   vec3 RayD = normalize(vec3(st.x,st.y,ef));
                     
                              RayD = rotatey(.003*u_mouse.x-1.)*RayD;
                              RayD = rotatex(.003*u_mouse.y-1.)*RayD;


             float seed = st.x + st.y * 3.43121412313;        
               





             ///Final Pipline Render/////   

              vec3 bufferA = vec3(0.0);
               for( int i = 1; i <= SAMPLES; i++ )
               
               
                bufferA = mix(bufferA, rBuffer(RayD,RayO,st,vec2(seed+float(i))).xyz, 1.0/float(i));
                   

                   
         

             vec3 color = aces(bufferA);

             color = pow(color,vec3(1./2.2));
                                      
          gl_FragColor = vec4(vec3(color),1.);
    }