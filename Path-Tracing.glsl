
//بِسْمِ ٱللَّٰهِ ٱلرَّحْمَٰنِ ٱلرَّحِيمِ


 precision highp float;

   #define PI 3.141592653589
   #define motionblur true
   
                uniform mat4 projectionMatrix;

   #define isBouncing true
    
   #define f 0.0000001
  uniform vec2 u_resolution; 
  uniform float u_time;
  uniform int u_frame;
  uniform vec3 u_camera;

   //Given' 
     ///Ray Origin'
       ///Ray Direction' 

   #define LIMIT 1000
   #define FOV 50

    ///Pre-Struct///
      vec3 vRayDirection;
      vec3 vRayOrigin;
       vec3 vNormal;
        vec3 vPosition;
        vec3 vColor;
         float vDistance;
  
    /////////////////////////////
    /////Uniformic Function//////
    /////////////////////////////
    
    // Narkowicz 2015, "ACES Filmic Tone Mapping Curve"
vec3 aces(vec3 x) {
  const float a = 2.51;
  const float b = 0.03;
  const float c = 2.43;
  const float d = 0.59;
  const float e = 0.14;
  return clamp((x * (a * x + b)) / (x * (c * x + d) + e), 0.0, 1.0);
}
//// rotational Matrix 

vec2 rotational(vec2 v, float a) {
	float s = sin(a);
	float c = cos(a);
	mat2 m = mat2(c, s, -s, c);
	return m * v;
}

    ///hitpoint Function /// 

     vec3 gHitpoint(vec3 rayDirection, vec3 rayOrigin, float gDistance ){
         return (rayDirection*gDistance)+rayOrigin;
     }     
        /// Vector to light equation   /// For point Lights
       vec3 vecToLight(vec3 hitpoint, vec3 lightPosition){
         vec3 lightVectorConst = lightPosition-hitpoint;
         return normalize(lightVectorConst);
      } 
       
        ///Random Intersection Test Function/// 
        bool intersected(vec3 rayDirection, vec3 rayOrigin, float gDistance) {
            vec3 inHitpoint = gHitpoint(rayDirection,rayOrigin,gDistance); 
             float setPointinSpace = inHitpoint.z+gDistance;
             float dotTestIntersection = dot((normalize(rayOrigin)),inHitpoint);
            if (dotTestIntersection > gDistance) {
                return true;
            } 
            else {
                return false;
            }
        }
    ///////////////////////////////       
    /////GeneralRandom Function//// 
    //////////////////////////////
//---------------------------------------------------------------------------------------------------------------
vec4 hash( vec2 gridcell )
{
//    gridcell is assumed to be an integer coordinate
   vec2 OFFSET = vec2( 26.0, 161.0222222222 );
const float DOMAIN = 71.0;
const float SOMELARGEFLOAT = 9222.135664;
vec4 P = vec4( gridcell.xy, gridcell.xy + vec2(45.0) );
P = P - floor(P * ( 1.0 / DOMAIN )) * DOMAIN;    //    truncate the domain
P += OFFSET.xyxy;                                //    offset to interesting part of the noise
P *= P;                                          //    calculate and return the hash
return fract( P.xzxz * P.yyww * ( .2220 / SOMELARGEFLOAT ) );
}


    ////Global///////

    int FresnelEffect_float(vec3 Normal, vec3 ViewDir, float Power,float bias, float scale, out float Out)
{
    Out = max(0.,min(1.,bias+scale*(1.0+dot(ViewDir,Normal)*Power)));

    return 1;
}

     vec3 inUnitSphere(vec3 Normal, vec2 seed) {

         seed = rotational(seed, 0221.*u_time);

       vec2 r = hash(seed*333.-123.*u_time).xy*1.0;
	vec3  uu = normalize( cross( Normal, vec3(0.0,1.0,1.0) ) );
	vec3  vv = cross( uu, Normal );
	float ra = sqrt(r.y);
	float rx = ra*cos(6.2831*r.x); 
	float ry = ra*sin(6.2831*r.x);
	float rz = sqrt( 1.0-r.y );
	vec3  rr = vec3( rx*uu + ry*vv + rz*Normal );
    return normalize( rr );

     } 

     int iSphere(vec3 rayDirection, vec3 rayOrigin, vec3 vaPosition, vec3 vaColor,
      float vaRadius, vec3 wNormal, float wDistance, vec3 wPosition, vec3 wColor
      ,out vec3 gNormal,out vec3 gPosition,out vec3 gColor,out float gDistance) {

   
        float dotPosition = max(dot((vaPosition-rayOrigin),rayDirection),0.0);  /// length ////   
           float bDistance = distance(vaPosition,gHitpoint(rayDirection,rayOrigin,dotPosition));
           float aIntersection = sqrt(pow(vaRadius,2.)-pow(bDistance,2.));
        float Distance = (dotPosition-aIntersection); 
          float rays_distance = distance(rayOrigin,vaPosition);

         bool rayHit = (rays_distance > (vaRadius+f)); 

        int relativeDistance = int(bDistance < vaRadius)*int(rayHit);

        vaColor = vaColor*vec3(relativeDistance);
          ///intersection bool/// 
        if(relativeDistance != 0) {
            Distance = Distance;
        }
        else {
          Distance = float(LIMIT);
        }
           gDistance = min(Distance,wDistance);
        ///Normal Struct/// 
        vec3 Normal = gHitpoint(rayDirection,rayOrigin,gDistance)-vaPosition;
         Normal = normalize(Normal);

           bool cDistance = (gDistance > wDistance);
          if (gDistance >= wDistance) {gNormal = wNormal;
                                      gColor = wColor;
                                      gPosition = wPosition;       
                            }
          else {gNormal = Normal;
                gColor = vaColor;
                gPosition = vaPosition;
    }
    return 1;

     }  
       ///Draw Function /// 
      int drawIntersection(vec3 rayDirection, vec3 rayOrigin, out vec3 gNormal, out vec3 gPosition, out vec3 gColor, out float gDistance){
          //varaiable//

         float bounce = 01.0*(abs(cos((4./1.0*u_time)))*abs(cos(u_time/6.0))); 

          ///S1
            vec3 position_1 = vec3(0.2020, bounce-.22+f, 0.984);
            vec3 color_1 = vec3(0.16, 0.42, 0.42);
            float radius_1 = 0.2222;
            float roughness_1 = 1.0;
            float reflectivity_1 = 0.0;
            //S2
            vec3 position_3 = vec3(-0.5320, -0.25, 0.9);
            vec3 color_3 = vec3(0.65);
            float roughness_2 = 1.0;
            float reflectivity_2 = 0.0;
            float radius_3 = 0.2145;
            //S3
            vec3 position_2 = vec3(0.310, -12.6620259, 1.1);
            vec3 color_2 = vec3(0.35, 0.19, 0.57);
            float radius_2 = 012.22;
            float roughness_3 = 1.0;
            float reflectivity_3 = 0.0;
            //S4
           vec3 position_4 = vec3(-0.2020, -0.250259, 0.984+0.5);
            vec3 color_4 = vec3(1.0);
            float radius_4 = 0.2;
            float roughness_4 = 1.0;
            float reflectivity_4 = 0.0;        

          //Arttibutes//
           vec3 globalNormal = vec3(0.0); //Apply at First functions///
           float infDistance = float(LIMIT); 
         ///Draw Scene//// 
        ///Call Function !
        
        int sphere_1 = iSphere(rayDirection,rayOrigin,position_1,color_1,radius_1,vec3(0.),infDistance,vec3(0.0),vec3(0.0),gNormal,gPosition,gColor,gDistance);
        int sphere_2 = iSphere(rayDirection,rayOrigin,position_2,color_2,radius_2,gNormal,gDistance,gPosition,gColor,gNormal,gPosition,gColor,gDistance); 
        int sphere_3 = iSphere(rayDirection,rayOrigin,position_3,color_3,radius_3,gNormal,gDistance,gPosition,gColor,gNormal,gPosition,gColor,gDistance);
        int sphere_4 = iSphere(rayDirection,rayOrigin,position_4,color_4,radius_4,gNormal,gDistance,gPosition,gColor,gNormal,gPosition,gColor,gDistance);  
        return 1;
      }
           ///Draw indirect Loop///    
      vec3 inShadeData(in vec3 hitpoint, vec3 lightPosition, vec3 Normal, vec3 inColor, float Intensity) {
                         vec3 inVectorLight = vecToLight(hitpoint,lightPosition);   
            vec3 inShadeMap = vec3(float(max(dot(Normal,inVectorLight),0.0)));
            float inRelativePoint = length(lightPosition-hitpoint);
            inShadeMap = inShadeMap/pow(inRelativePoint,2.);

            ///Render Shadow Map///
            float inRayDistance;
            vec3 inRayBreak;
            int rayGen = drawIntersection(inVectorLight,hitpoint,
            inRayBreak,inRayBreak,inRayBreak,inRayDistance);  
            bool inShadowRelative = (inRayDistance >= inRelativePoint);
                 inShadeMap = inShadeMap*inColor;
                inShadeMap = (inShadeMap*float(inShadowRelative))*Intensity; 
              vec3 inDirect_Lighting = 0.554321*vec3(inShadeMap);   

              return inDirect_Lighting;
      }


     ///Indirect Test Loop Function////
      
          int inDirectLightIntersectionLoop(inout vec3 Normal, inout vec3 hitpoint,out vec3 color, in vec2 seed,vec3 lightPosition,float Intensity , out vec3 Final) {
            
                          vec3 rayOrigin = hitpoint; 
                          vec3 rayDirection = inUnitSphere(Normal,seed);
                           vec4 inbreak;

                          int indirectIntersection = drawIntersection(rayDirection,rayOrigin,Normal,inbreak.zxy,color,inbreak.x);
                          hitpoint = gHitpoint(rayDirection,rayOrigin,inbreak.x);
                        
                        Normal = Normal;

                      //Draw//
                     Final = inShadeData(hitpoint,lightPosition,Normal,color,Intensity);
                      
            return 1;
          }



      //Generate Light Function/// 
      vec3 pixelShade(in vec3 Normal,in vec3 lightPosition, vec3 hitpoint, vec3 lightColor, float Intensity, inout vec3 wShadeData, in vec2 seed)
      {

           seed = rotational(seed, 0221.*u_time);
           vec3 jitter = ((hash(floor(seed*343.-123.*u_time)).xyz)-0.5)*1.;
           lightPosition = jitter+lightPosition;
           vec3 vectorLight = vecToLight(hitpoint,lightPosition);
           float shadeMap = max(dot(Normal,vectorLight),0.0);
           float relativePoint = length(lightPosition-hitpoint);
            shadeMap = shadeMap/pow(relativePoint,2.);

            ///Render Shadow Map///
            float rayDistance;
            vec3 rayBreak;
            int rayGen = drawIntersection(vectorLight,hitpoint,
            rayBreak,rayBreak,rayBreak,rayDistance);  
            bool shadowRelative = (rayDistance >= relativePoint);
            shadeMap = max((shadeMap*float(shadowRelative)),0.0)*Intensity;


             ///Calculate Indirect lighting /// 
              vec3 inDirect_Lighting = vec3(0.0, 0.0, 0.0); 

              const int Global_illumination_bounce = 3;


            if (Global_illumination_bounce >= 1) { 
                          vec3 inNormal;
                          vec3 inColor;
                          float inDistance;
                          vec3 inBreak;
                          vec3 rayOrigin = hitpoint; 
                          vec3 rayDirection = inUnitSphere(Normal,seed);

                          int perPixel = drawIntersection(rayDirection,rayOrigin,inNormal,inBreak,inColor,inDistance);
                         vec3 vectorLight = vecToLight(hitpoint,lightPosition);

                 vec3 inVectorLight = vecToLight(gHitpoint(rayDirection,rayOrigin,inDistance),lightPosition);   
            vec3 inShadeMap = vec3(float(max(dot(inNormal,inVectorLight),0.0)));
            float inRelativePoint = length(lightPosition-gHitpoint(rayDirection,rayOrigin,inDistance));
            inShadeMap = inShadeMap/pow(inRelativePoint,2.);

            ///Render Shadow Map///
            float inRayDistance;
            vec3 inRayBreak;
            int rayGen = drawIntersection(inVectorLight,gHitpoint(rayDirection,rayOrigin,inDistance),
            inRayBreak,inRayBreak,inRayBreak,inRayDistance);  
            bool inShadowRelative = (inRayDistance >= inRelativePoint);
                 inShadeMap = inShadeMap*inColor;
                inShadeMap = (inShadeMap*float(inShadowRelative))*Intensity; 
                inDirect_Lighting = 0.654321*vec3(inShadeMap);  

                      if (Global_illumination_bounce > 1) {
                         vec3 Final;

                       hitpoint = gHitpoint(rayDirection,rayOrigin,inDistance);
                         
                        for( int i = 1; i <= Global_illumination_bounce ; i++ ) {
                          int j = inDirectLightIntersectionLoop(inNormal,hitpoint,inColor,seed+float(i),lightPosition,Intensity,Final)
                          ;
                                 Final *= inColor;
                          
                          inDirect_Lighting += max(Final,0.00000000000001);
                            
                           
                        }
                      } 
                      
                            } 

     return (lightColor*(shadeMap+inDirect_Lighting))+wShadeData;
       // return inDirect_Lighting;
      } 
  
      ///Draw Light data///
      vec3 pixelData(in vec3 hitpoint, vec3 Normal, vec3 gColor, in vec2 seed){
        vec3 spaceShade ;
        vec3 inShade = vec3(1);
          const vec3 light_position01 = vec3(1.0, 1.0, -0.01);
          const vec3 light_color01 = vec3(0.73);
          const float light_intansity01 = 4.;
           
           ///Draw Scene///
          vec3 shade_01 = pixelShade(Normal,light_position01,hitpoint,light_color01,light_intansity01,spaceShade,seed);
          return shade_01*gColor;
      }  

      ///Generate Reflection///
      int inReflections(inout vec3 rayDirection, inout vec3 rayOrigin, inout vec3 hitpoint, inout vec3 Normal,inout vec2 seed, out vec3 Final) {
           vec3 inbreak;  vec3 normal; vec3 position; vec3 color; float distance;

           seed = rotational(seed, 0221.*u_time);
           vec3 jitter = ((hash(floor(seed*343.-123.*u_time)).xyz)-0.5)*.1;

           const float ior = 1.490;
               
               float roughness = 0.61;
           
              rayDirection = reflect(rayDirection,Normal+(jitter*roughness));
              rayOrigin = hitpoint;
           int intersectionReflections = drawIntersection(rayDirection,rayOrigin,Normal,position,color,distance);
               hitpoint = gHitpoint(rayDirection,rayOrigin,distance);

                seed += (521.+f);
               Final = pixelData(hitpoint,Normal,color,seed);

               Final = max(Final,0.00000000000001);

           return 1;
      } 



       ////Reflection Arry//// 
        int reflectionArray( vec3 rayDirection,  vec3 rayOrigin, in float rayDistance, in vec2 seed, in vec3 Normal, inout vec3 pixelFragment) {

            vec3 normal; vec3 position; vec3 color; float distance;
               vec3 hitpoint = gHitpoint(rayDirection,rayOrigin,rayDistance);
                vec3 Final = vec3(0.0);
                

                float fresnal = pow((dot(Normal,rayDirection)+1.0),4.5);
   
                
             if (true) {

               for(int r = 1; r <= 1; r++) {
                 int j = inReflections(rayDirection,rayOrigin,hitpoint,Normal,seed,Final);
                
               }

             }


                  float reflectivity = 0.1;

                        reflectivity = mix(reflectivity,1.0,fresnal);
   
                pixelFragment = mix(pixelFragment,Final,reflectivity  );
          

          
          return 1;
        }
        
      ////Renderer Piple line//// 
      vec3 rBuffer(inout vec3 rayDirection, vec3 rayOrigin, vec2 seed) {

              
                    int render = drawIntersection(vRayDirection,vRayOrigin,
                        vNormal,vPosition,vColor,vDistance);
                           vec3 hitPoint = gHitpoint(rayDirection,rayOrigin,vDistance);  
                             vec3 pixelFragment = pixelData(hitPoint,vNormal,vColor,seed);  
                          
                          int vReflection = reflectionArray(rayDirection,rayOrigin,vDistance,seed,vNormal,pixelFragment);

                          if (motionblur != false) {

                            float h = 0.001;
                               vec2 momentum ;

                          }
            return pixelFragment;
      } 

      vec3 getRayDirection(vec3 cameraDir, vec3 up, vec2 pixelCoord, vec2 resolution, float fov) {
    float aspect_ratio = resolution.x / resolution.y;
    float x = (2.0 * (pixelCoord.x + 0.5) / resolution.x - 1.0) * aspect_ratio * tan(fov / 2.0);
    float y = (1.0 - 2.0 * (pixelCoord.y + 0.5) / resolution.y) * tan(fov / 2.0);

    // Basis vectors for the camera
    vec3 w = normalize(cameraDir);
    vec3 u = normalize(cross(up, w));
    vec3 v = cross(w, u);

    // Calculate the ray direction
    return normalize(x * u + y * v + w);
}


       void main(){ 
              vec2 coord = gl_FragCoord.xy/u_resolution.xy;
                    coord = (2.0*coord)-1.0;   
                     float a = u_resolution.x/u_resolution.y;
                      coord.x *= a;
                
                float rFocal = ((1.0/tan(radians(float(FOV))/2.)));
                   rFocal *= 1.0;

               vRayOrigin = vec3(0.0,0.0,-0.60);
            

            vec3 viewPoint = vec4(u_camera,1.0).xyz-vRayOrigin-2.5;

             vRayDirection =  normalize(vec3(coord.x,coord.y,rFocal));


               vec3 up = vec3(0.0, 1.0, 0.0);

              //vRayDirection = getRayDirection(u_camera,up,coord,u_resolution,radians(float(FOV)));


               

  
                  
      
            ///  vec3 render = rBuffer(vRayDirection,vRayOrigin,coord);
                vec3 render = vec3(0.0);


                ///Pipline Renderer////
                const int sample = 8 ;

                      for( int i = 1; i <= sample ; i++ )
               
               
                render = mix(render, rBuffer(vRayDirection,vRayOrigin,coord+vec2(i)).xyz, 1.0/float(i));
                  render = aces(render);
                  render = pow(render,vec3(1./2.2));

          gl_FragColor = vec4(vec3(render),1.);  
       
       }

 


      
       
       