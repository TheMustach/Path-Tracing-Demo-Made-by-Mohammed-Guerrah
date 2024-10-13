



////بِسْمِ ٱللَّٰهِ ٱلرَّحْمَٰنِ ٱلرَّحِيمِ//
///Hello 

          #ifdef GL_ES
        precision mediump float;
     #endif


   uniform vec2 u_resolution;
   uniform vec2 u_mouse;
   uniform float u_time;

 










#define eps 0.0002
#define EYEPATHLENGTH 4
#define SAMPLES 25

#define SHOWSPLITLINE
#define FULLBOX

#define DOF
#define ANIMATENOISE
#define MOTIONBLUR

#define MOTIONBLURFPS 12.

#define denoiser .0

#define LIGHTCOLOR vec3(16.86, 10.76, 8.2)*3.
#define WHITECOLOR vec3(1.0, 1.0, 1.0)*0.7
#define GREENCOLOR vec3(1.0, 0.0, 0.0)*.12
#define REDCOLOR vec3(1.0, 0.2667, 0.0235)*0.7





//ToneMap 
// 0 == sRGB
// 1 == Filmic
// 2 == ACES.


#define TONEMAP 2

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


 /// ACES tone mapping 

   
   vec3 aces(vec3 x) {
  const float a = 2.51;
  const float b = 0.03;
  const float c = 2.43;
  const float d = 0.59;
  const float e = 0.14;
  return clamp((x * (a * x + b)) / (x * (c * x + d) + e), 0.0, 1.0);
}



///Seed Function;
float hash1(inout float seed) {
    return fract(sin(seed += 0.1)*43758.5453123);
}

vec2 hash2(inout float seed) {
    return fract(sin(vec2(seed+=0.1,seed+=0.1))*vec2(43758.5453123,22578.1459123));
}

vec3 hash3(inout float seed) {
    return fract(sin(vec3(seed+=0.1,seed+=0.1,seed+=0.1))*vec3(43758.5453123,22578.1459123,19642.3490423));
}

//-----------------------------------------------------
// Intersection functions (by iq)
//-----------------------------------------------------

vec3 nSphere( in vec3 pos, in vec4 sph ) {
    return (pos-sph.xyz)/sph.w;
}

float iSphere( in vec3 ro, in vec3 rd, in vec4 sph ) {
    vec3 oc = ro - sph.xyz;
    float b = dot(oc, rd);
    float c = dot(oc, oc) - sph.w * sph.w;
    float h = b * b - c;
    if (h < 0.0) return -1.0;

	float s = sqrt(h);
	float t1 = -b - s;
	float t2 = -b + s;
	
	return t1 < 0.0 ? t2 : t1;
}

vec3 nPlane( in vec3 ro, in vec4 obj ) {
    return obj.xyz;
}

float iPlane( in vec3 ro, in vec3 rd, in vec4 pla ) {
    return (-pla.w - dot(pla.xyz,ro)) / dot( pla.xyz, rd );
}

//-----------------------------------------------------
// scene
//-----------------------------------------------------

vec3 cosWeightedRandomHemisphereDirection( const vec3 n, inout float seed ) {
  	vec2 r = hash2(seed);
    
	vec3  uu = normalize( cross( n, vec3(0.0392, 0.0627, 0.4941) ) );
	vec3  vv = cross( uu, n );
	
	float ra = sqrt(r.y);
	float rx = ra*cos(6.2831*r.x); 
	float ry = ra*sin(6.2831*r.x);
	float rz = sqrt( 1.0-r.y );
	vec3  rr = vec3( rx*uu + ry*vv + rz*n );
    
    return normalize( rr );
}

vec3 randomSphereDirection(inout float seed) {
    vec2 h = hash2(seed) * vec2(2.,6.28318530718)-vec2(1,0);
    float phi = h.y;
	return vec3(sqrt(1.-h.x*h.x)*vec2(sin(phi),cos(phi)),h.x);
}

vec3 randomHemisphereDirection( const vec3 n, inout float seed ) {
	vec3 dr = randomSphereDirection(seed);
	return dot(dr,n) * dr;
}

//-----------------------------------------------------
// light
//-----------------------------------------------------

vec4 lightSphere;
vec4 lightSphere02;
void initLightSphere( float time ) {
              
     
	lightSphere = vec4( -3.0,2.8+2.*sin(2.*0.9),-cos(u_time), .5 );
     
        lightSphere.xyz = rotatey(u_time)*lightSphere.xyz;
        lightSphere02 = vec4(-lightSphere.x,lightSphere.y,-lightSphere.z,.7);
}

vec3 sampleLight( const in vec3 ro, inout float seed ) {
    vec3 n = randomSphereDirection( seed ) * lightSphere.w;
    vec3 n2 = randomSphereDirection( seed ) * lightSphere02.w;
    return (lightSphere.xyz + n)+(lightSphere02.xyz + n2);
}

//-----------------------------------------------------
// scene
//-----------------------------------------------------

vec2 intersect( in vec3 ro, in vec3 rd, inout vec3 normal ) {
	vec2 res = vec2( 1e20, -1.0 );
    float t;
	
	t = iPlane( ro, rd, vec4(0.0, 0.3255, 0.0706, -0.20) ); if( t>eps && t<res.x ) { res = vec2( t, 1. ); normal = vec3( 0., 1., 0.); }
	t = iPlane( ro, rd, vec4( 0.0, 0.0,-1.0,8.0 ) ); if( t>eps && t<res.x ) { res = vec2( t, 1.3 ); normal = vec3( 0., 0.,-1.); }
    //t = iPlane( ro, rd, vec4( 1.0, 0.0, 0.0,0.0 ) ); if( t>eps && t<res.x ) { res = vec2( t, 2. ); normal = vec3( 1., 0., 0.); }
#ifdef FULLBOX
   //t = iPlane( ro, rd, vec4( 0.0,-1.0, 0.0,5.49) ); if( t>eps && t<res.x ) { res = vec2( t, 1. ); normal = vec3( 0., -1., 0.); }
    t = iPlane( ro, rd, vec4(-1.0, 0.0, 0.0,5.59) ); if( t>eps && t<res.x ) { res = vec2( t, 3. ); normal = vec3(-1., 0., 0.); }
#endif

	t = iSphere( ro, rd, vec4( 1.5,1.0, 2.7, 1.0) ); if( t>eps && t<res.x ) { res = vec2( t, 1. ); normal = nSphere( ro+t*rd, vec4( 1.5,1.0, 2.7,1.0) ); }
    t = iSphere( ro, rd, vec4( 4.0,1.0, 4.0, 1.0) ); if( t>eps && t<res.x ) { res = vec2( t, 6. ); normal = nSphere( ro+t*rd, vec4( 4.0,1.0, 4.0,1.0) ); }
    t = iSphere( ro, rd, vec4( -4.0,2.0, 4.0, 1.0) ); if( t>eps && t<res.x ) { res = vec2( t, 6.1 ); normal = nSphere( ro+t*rd, vec4( -4.0,2.0, 4.0,1.0) ); }
     t = iSphere( ro, rd, vec4( -2.0,1.5, .0, 1.0) ); if( t>eps && t<res.x ) { res = vec2( t, 6.2 ); normal = nSphere( ro+t*rd, vec4( -2.0,1.5, .0,1.0) ); }
    t = iSphere( ro, rd, lightSphere ); if( t>eps && t<res.x ) { res = vec2( t, 0.0 );  normal = nSphere( ro+t*rd, lightSphere ); }
	t = iSphere( ro, rd, lightSphere02 ); if( t>eps && t<res.x ) { res = vec2( t, 0.0 );  normal = nSphere( ro+t*rd, lightSphere ); }				  
    return res;					  
}

bool intersectShadow( in vec3 ro, in vec3 rd, in float dist ) {
    float t;
	
	t = iSphere( ro, rd, vec4( 1.5,1.0, 2.7,1.0) );  if( t>eps && t<dist ) { return true; }
    t = iSphere( ro, rd, vec4( 4.0,1.0, 4.0,1.0) );  if( t>eps && t<dist ) { return true; }
    t = iSphere( ro, rd, vec4( -4.0,2.0, 4.0,1.0) );  if( t>eps && t<dist ) { return true; }
    t = iSphere( ro, rd, vec4( -2.0,1.5, .0,1.0) );  if( t>eps && t<dist ) { return true; }
    return false; // optimisation: planes don't cast shadows in this scene
}

//-----------------------------------------------------
// materials
//-----------------------------------------------------

vec3 matColor( const in float mat ) {
	vec3 nor = vec3(1.0, 1.0, 1.0);
	
    if (mat ==6.1) nor = vec3(0.3137, 0.8667, 0.2431);
    if (mat ==6.2) nor = vec3(0.0431, 0.0, 0.5216);
	if( mat<3.5 ) nor = REDCOLOR;
    if( mat<2.5 ) nor = GREENCOLOR;
	if( mat<1.5 ) nor = WHITECOLOR;
    if(mat == 1.3) nor = vec3(1.0, 0.2353, 0.0);
    if(mat == 1.2 ) nor = vec3(0.2);
	if( mat<0.5 ) nor = LIGHTCOLOR;
					  
    return nor;					  
}

bool matIsSpecular( const in float mat ) {
    return mat > 4.5;
    return mat == 1.2;
}

bool matIsLight( const in float mat ) {
    return mat < .5;
}

//-----------------------------------------------------
// brdf
//-----------------------------------------------------

vec3 getBRDFRay( in vec3 n, const in vec3 rd, const in float m, inout bool specularBounce, inout float seed ) {
    specularBounce = false;
    
    vec3 r = cosWeightedRandomHemisphereDirection( n, seed );
    if(  !matIsSpecular( m ) ) {
        return r;
    } else {
        specularBounce = true;
        
        float n1, n2, ndotr = dot(rd,n);
        
        if( ndotr > 0. ) {
            n1 = 1.1; 
            n2 = 1.5;
            n = -n;
        } else {
            n1 = 1.5;
            n2 = 1.0; 
        }
                
        float r0 = (n1-n2)/(n1+n2); r0 *= r0;
		float fresnel = r0 + (1.-r0) * pow(1.0-abs(ndotr),5.);
        
        vec3 ref =  refract( rd, n, n2/n1 );        
        if( ref == vec3(0) || hash1(seed) < fresnel ) {
            ref = reflect( rd, n );
        } 
        
        return ref;
	}
}

//-----------------------------------------------------
// eyepath
//-----------------------------------------------------

vec3 traceEyePath( in vec3 ro, in vec3 rd, const in bool directLightSampling, inout float seed ) {
    vec3 tcol = vec3(0.0, 0.0, 0.0);
    vec3 fcol  = vec3(1.0, 1.0, 1.0);
    
    bool specularBounce = true;
    
    for( int j=0; j<EYEPATHLENGTH; ++j ) {
        vec3 normal;
        
        vec2 res = intersect( ro, rd, normal );
        if( res.y < -0.5 ) {
            return tcol;
        }
        
        if( matIsLight( res.y ) ) {
            if( directLightSampling ) {
            	if( specularBounce ) tcol += fcol*LIGHTCOLOR;
            } else {
                tcol += fcol*LIGHTCOLOR;
            }
            return tcol;
        }
        
        ro = ro + res.x * rd;
        rd = getBRDFRay( normal, rd, res.y, specularBounce, seed );
            
        if(!specularBounce || dot(rd,normal) < 0.) {  
        	fcol *= matColor( res.y );
        }
        
        if( directLightSampling ) {
        	vec3 ld = sampleLight( ro, seed ) - ro;
			vec3 nld = normalize(ld);
            if( !specularBounce && j < EYEPATHLENGTH-1 && !intersectShadow( ro, nld, length(ld)) ) {

                float cos_a_max = sqrt(1. - clamp(lightSphere.w * lightSphere.w / dot(lightSphere.xyz-ro, lightSphere.xyz-ro), 0., 1.));
                float weight = 2. * (1. - cos_a_max);

                tcol += (fcol * LIGHTCOLOR) * (weight * clamp(dot( nld, normal ), 0., 1.));
            }
        }
    }    
    return tcol;
}

//-----------------------------------------------------
// main
//-----------------------------------------------------



























void main() {

    
    vec2 iResolution = u_resolution.xy; float iTime = u_time;
      vec2 fragCoord = gl_FragCoord.xy;
	vec2 q = fragCoord.xy / iResolution.xy;
    
	float splitCoord = denoiser*3000.;
    bool directLightSampling = fragCoord.x < splitCoord;
    
    //-----------------------------------------------------
    // camera
    //-----------------------------------------------------

    vec2 p = -1.0 + 3.0 * (fragCoord.xy) / iResolution.xy;
    p.x *= iResolution.x/iResolution.y;

#ifdef ANIMATENOISE
    float seed = p.x + p.y * 3.43121412313 + fract(1.12345314312*iTime);
#else
    float seed = p.x + p.y * 3.43121412313;
#endif
    

    vec2 phi = .003*vec2(u_mouse.x,u_mouse.y);


    vec3 ro = vec3(2.78, 6.73, -8.00);
    vec3 ti = vec3(2.78, 2.3,  0.00);
    vec3 ta = rotatey(phi.x-1.)*ti*rotatex(phi.y-1.);
    vec3 ww = normalize( ta - ro );
    vec3 uu = normalize( cross(ww,vec3(0.0,1.0,0.0) ) );
    vec3 vv = normalize( cross(uu,ww));

    //-----------------------------------------------------
    // render
    //-----------------------------------------------------

    vec3 col = vec3(0.0);
    vec3 tot = vec3(0.0);
    vec3 uvw = vec3(0.0);
    
    for( int a=0; a<SAMPLES; a++ ) {

        vec2 rpof = 2.*(hash2(seed)-vec2(0.5)) / iResolution.y;
	    vec3 rd = normalize( (p.x+rpof.x)*uu + (p.y+rpof.y)*vv + 3.0*ww );
        
#ifdef DOF
	    vec3 fp = ro + rd * 12.0;
   		vec3 rof = ro + (uu*(hash1(seed)-0.5) + vv*(hash1(seed)-0.5))*0.125;
    	rd = normalize( fp - rof );
#else
        vec3 rof = ro;
#endif        
        
#ifdef MOTIONBLUR
        initLightSphere( iTime + hash1(seed) / MOTIONBLURFPS );
#else
        initLightSphere( iTime );        
#endif
        
        col = traceEyePath( rof, rd, directLightSampling, seed );

        tot += col;
        
        seed = mod( seed*1.1234567893490423, 13. );
    }
    
    tot /= float(SAMPLES);
    
#ifdef SHOWSPLITLINE
	if (abs(fragCoord.x - splitCoord) < 1.0) {
		tot.x = 1.0;
	}
#endif
    
	tot =  clamp(tot,0.0,1.0);



    int deftone = TONEMAP;

     if ( deftone == 2 ) 
     { tot = aces(tot);
     tot = pow(tot,vec3(1./2.2));
      }
        else {
           if (deftone ==1 ) {} else { 
            tot = pow(tot,vec3(1./2.2));
           }
        }
       
     
          


    gl_FragColor = vec4( tot, 1.0 );
}