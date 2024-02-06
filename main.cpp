#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <deque>
#include <stack>
#include "Transform.h"
#include <FreeImage.h>
#include <algorithm>
#include <vector>

using namespace std;
 

#define MAINPROGRAM
#include "variables.h"
#include "readfile.h"
 
const float PI = 3.14159265;
class Ray
{
public:
    vec3 p0;
    vec3 p1;
    Ray(vec3 p1, vec3 p0);
    Ray();
    virtual ~Ray();
    
};
 
class Intersection
{
public:
    bool hit;
    float dist;
    int obji;
    vec3 norm;
    vec3 P;
    vec3 dir;
    Intersection();
    Intersection(bool hit, float dist, int obji, vec3 norm, vec3 P, vec3 dir);
    virtual ~Intersection();
    
};

Intersection::Intersection(bool hit, float dist, int obji, vec3 norm , vec3 P, vec3 dir){
        this->hit = hit;
        this->dist = dist;
        this-> obji = obji;
        this->norm = norm;
        this-> P = P;
        this->dir = dir;
}
 
Intersection::~Intersection(){
    
}

Intersection::Intersection(){
    
}
 
Ray::Ray(vec3 p1, vec3 p0){
    this->p1 = p1;
    this->p0 = p0;
}
Ray::~Ray(){
    
}
 
Ray::Ray(){
    
}
 
 
void saveScreenshot(string fname,BYTE * pixels)
{
    FIBITMAP *img = FreeImage_ConvertFromRawBits(pixels, w, h, w * 3, 24, 0xFF0000, 0x00FF00, 0x0000FF, true);
    std::cout << "Saving screenshot: " << fname << "\n";
    FreeImage_Save(FIF_PNG, img, fname.c_str(), 0);
}

struct UGrid {
    vector<object *> grid[5][5][5];
    vec3 t_min, t_max;
    vec3 gridDim;
    vec3 cellsize;
};
 
UGrid ObjtoGrid() {
    UGrid gridObj;
    gridObj.t_max = vec3(-INFINITY);
    gridObj.t_min = vec3(INFINITY);
    for (int i=0; i<numobjects; i++) {
        object * obj = &(objects[i]);
        vec3 objMin, objMax;
        if(obj->type == sphere){
           vec3 position =  vec3(obj->transform * vec4(obj->centerSphere, 1));
           objMin[0] = position[0] - obj->size;
           objMin[1]  = position[1] - obj->size;
           objMin[2]  = position[2] - obj->size;
           objMax[0] = position[0] + obj->size;
           objMax[1] = position[1] + obj->size;
           objMax[2] = position[2] + obj->size;
        }
        else if(obj->type == triangle) {
            vec3 A = vec3(obj->transform * vec4(vertexpoint[(obj->point)[0]],1));
            vec3 B = vec3(obj->transform * vec4(vertexpoint[(obj->point)[1]],1));
            vec3 C = vec3(obj->transform * vec4(vertexpoint[(obj->point)[2]],1));
            objMin[0] = std::min(A[0], std::min(B[0], C[0]));
            objMin[1] = std::min(A[1], std::min(B[1], C[1]));
            objMin[2] = std::min(A[2], std::min(B[2], C[2]));
            objMax[0] = std::max(A[0], std::max(B[0], C[0]));
            objMax[1] = std::max(A[1], std::max(B[1], C[1]));
            objMax[2] = std::max(A[2], std::max(B[2], C[2]));
        }
        if(objMin[0] < gridObj.t_min[0]){
            gridObj.t_min[0] = objMin[0];
        }
        if(objMin[1] < gridObj.t_min[1]){
            gridObj.t_min[1] = objMin[1];
        }
        if(objMin[2] < gridObj.t_min[2]){
            gridObj.t_min[2] = objMin[2];
        }
        if(objMax[0] > gridObj.t_max[0]){
            gridObj.t_max[0] = objMax[0];
        }
        if(objMax[1] > gridObj.t_max[1]){
            gridObj.t_max[1] = objMax[1];
        }
        if(objMax[2] > gridObj.t_max[2]){
            gridObj.t_max[2] = objMax[2];
        }
    }
    gridObj.cellsize = (gridObj.t_max - gridObj.t_min)/(float)5;
    if (gridObj.cellsize.x == 0) {
        gridObj.cellsize.x = 0.2;
    }
    else if (gridObj.cellsize.y == 0) {
        gridObj.cellsize.y = 0.2;
    }
    else if (gridObj.cellsize.z == 0) {
        gridObj.cellsize.z = 0.2;
    }
    for (int i=0; i<numobjects; i++) {
        object * obj = &(objects[i]);
        vec3 objMin, objMax;
        if(obj->type ==  sphere){
            vec3 position =  vec3(obj->transform * vec4(obj->centerSphere, 1));
            objMin.x = position[0] - obj->size;
            objMin.y  = position[1] - obj->size;
            objMin.z  = position[2] - obj->size;
            objMax.x = position[0] + obj->size;
            objMax.y = position[1] + obj->size;
            objMax.z = position[2] + obj->size;
        }
        else if(obj->type == triangle) {
            vec3 A = vec3(obj->transform * vec4(vertexpoint[(obj->point)[0]],1));
            vec3 B = vec3(obj->transform * vec4(vertexpoint[(obj->point)[1]],1));
            vec3 C = vec3(obj->transform * vec4(vertexpoint[(obj->point)[2]],1));
            objMin[0] = std::min(A[0], std::min(B[0], C[0]));
            objMin[1] = std::min(A[1], std::min(B[1], C[1]));
            objMin[2] = std::min(A[2], std::min(B[2], C[2]));
            objMax[0] = std::max(A[0], std::max(B[0], C[0]));
            objMax[1] = std::max(A[1], std::max(B[1], C[1]));
            objMax[2] = std::max(A[2], std::max(B[2], C[2]));
        }
        
        vec3 min = (objMin - gridObj.t_min)/gridObj.cellsize;
        //cout << min.x << " " << min.y << " " << min.z <<'\n';
        vec3 max = (objMax - gridObj.t_min)/gridObj.cellsize;
        //cout << max.x << " " << max.y << " " << max.z <<'\n';
        min = glm::max(min, vec3(0));
        max = glm::min(max, vec3(4));
        for (int z = min.z; z <= max.z; ++z){
            for (int y = min.y; y <= max.y; ++y){
                for (int x = min.x; x <= max.x; ++x){
                    cout << x << " " << y << " " << z <<'\n';
                    gridObj.grid[x][y][z].push_back(obj);
                }
            }
        }
    }

    return gridObj;
}

Ray RaythruPixel(int j, int i){
    double fovyR= fovy * PI/180;
    float alpha = tan(fovyR/2)*((float)w/h) * (((i + 0.5) - (float)w/2)/(w/2));
    float beta = tan(fovyR/2)* (((float)h/2 - (j + 0.5 ))/(h/2));
    vec3 p0 = eye;
    vec3 x = normalize(eye - center);
    vec3 u = normalize(cross(up, x));
    vec3 v = cross(x, u);
    vec3 dir = normalize((alpha* u) + (beta *v) - x);
    return Ray(dir, p0);
}
 
Intersection intersect(Ray ray, vector<object *> list ){
    float minimum = 0.00001;
    float min = INFINITY;
    int indexHit = -1;
    vec3 normal;
    vec3 newP;
    bool hit;
    vec3 dir;
    float t;
    vec3 P;
    vec3 N;
    vec3 d;
    object * hitObject = NULL;
    for (int i = 0; i< list.size(); i++) {
        object * obj = list[i];
        if(obj->type == triangle){
            vec3 p0 = ray.p0;
            d = ray.p1;
            vec3 A = vec3(obj->transform * vec4(vertexpoint[(obj->point)[0]],1));
            vec3 B = vec3(obj->transform * vec4(vertexpoint[(obj->point)[1]],1));
            vec3 C = vec3(obj->transform * vec4(vertexpoint[(obj->point)[2]],1));
            N = glm::normalize(cross(B-A, C-A));
            t =  (dot(A,N) - dot(p0,N))/(dot(d,N));
            P = p0 + t* d;
            vec3 alpha = cross(N, C-B)/dot(cross(N, C-B), A-C);
            float a = dot(alpha,P) + -dot(alpha, C);
            vec3 beta = cross(N, A-C)/dot(cross(N, A-C), B-A);
            float b = dot(beta,P) + -dot(beta, A);
            vec3 gamma = cross(N, B-A)/dot(cross(N, B-A), C-B);
            float g = dot(gamma, P) + -dot(gamma,B);
         
            if(a < 0 || a > 1 || b < 0 || b > 1||g < 0 || g > 1 ){
                hit = false;
                t = INFINITY;
            }
            else{
                hit = true;
            }
            
        }else if(obj->type == sphere){
            mat4 M_inv = inverse(obj->transform);
            vec3 ctr = obj->centerSphere;
            float r = obj->size;
            vec3 p0 = ray.p0;
            d = ray.p1;
            vec3 origin = vec3(M_inv * vec4(p0, 1));
            vec3 direction = vec3(M_inv * vec4(d, 0));
            float a = dot(direction,direction);
            float b = 2.0 * dot(direction,origin - ctr);
            float c = dot(origin - ctr,origin - ctr) - (float)(r * r);
            float descriminate = (b * b)-(4.0*a*c);
            if( descriminate < 0){
                hit = false;
                t = INFINITY;
            }
            else {
                hit = true;
                float s1 = (-b + sqrt(descriminate))/(2.0*a);
                float s2 = (-b - sqrt(descriminate))/(2.0*a);
                if(s1 >= 0 && s2 >= 0) {
                    t = std::min(s1, s2);
                }
                vec3 p = origin + t*direction;
                P = vec3(obj->transform * vec4(p,1));
                N =  normalize(vec3(transpose(M_inv)* vec4(p - obj->centerSphere,0)));
            }
        }
        if(hit == true){
            if(min > t && t > minimum){
                min = t;
                hitObject = obj;
                indexHit = i;
                dir = d;
                newP = P;
                normal = N;
            }
        }
    }
    return Intersection(hit,min,indexHit,normal,newP,dir);
}

Intersection DDA(Ray ray, UGrid grid){
    int t_x, t_y, t_z, t ;
    vec3 deltaT;
    Intersection rayIntersect;
    
   // cout << "eye" <<  ray.p0[0] << " " << ray.p0[1] << " " << ray.p0[2] << '\n';
    //cout << "t_min" << grid.t_min[0]  << " " << grid.t_min[1] << " " << grid.t_min[2] << '\n';
    //cout << "t_max" << grid.t_max[0]  << " " << grid.t_max[1] << " " << grid.t_max[2] << '\n';
    //cout << "cellsize" << grid.cellsize[0]  << " " << grid.cellsize[1] << " " << grid.cellsize[2] << '\n';
    vec3 dir = ray.p1;
    vec3 tnear = vec3(0,0,0);
    if( ray.p0[0] > grid.t_max[0] ){
        tnear[0] = (grid.t_max[0]- ray.p0[0]);
    }
    if( ray.p0[1] > grid.t_max[1] ){
        tnear[1] = (grid.t_max[1]- ray.p0[1]);
    }
    if( ray.p0[2] > grid.t_max[2] ){
        tnear[2] = (grid.t_max[2]- ray.p0[2]);
    }
    if(ray.p0[0] < grid.t_min[0] ){
        tnear[0] = grid.t_min[0] - ray.p0[0];
    }
    if(ray.p0[1] < grid.t_min[1] ){
        tnear[1] = grid.t_min[1] - ray.p0[1];
    }
    if(ray.p0[2] < grid.t_min[2] ){
        tnear[2] = grid.t_min[2] - ray.p0[2];
    }
   // cout << "tnear[2]" << tnear[2] << '\n';
    vec3 index =
    floor((ray.p0 + tnear - grid.t_min)/ grid.cellsize);
    if (index[0] == 5) {
        index[0] = 4;
    }
    if (index[1] == 5) {
        index[1] = 4;
    }
    if (index[2] == 5) {
        index[2] = 4;
    }
    
    //cout <<(ray.p0 + tnear - grid.t_min)[2] << '/n';
    vec3 cellDim = grid.cellsize;
    vec3 OGrid = ray.p0- grid.t_min; // rayOrigCell
    
    if(dir[0] < 0){
        deltaT[0] = - cellDim[0]/dir[0];
        t_x = (floor(OGrid[0]/cellDim[0]) * cellDim[0] - OGrid[0])/dir[0];
    }
    else{
        deltaT[0] = cellDim[0]/dir[0];
        t_x = ((floor(OGrid[0]/cellDim[0]) + 1) * cellDim[0] - OGrid[0])/dir[0];
    }
    if(dir[1] < 0){
        deltaT[1] = - cellDim[1]/dir[1];
        t_y = (floor(OGrid[1]/cellDim[1]) * cellDim[1] - OGrid[1])/dir[1];
    }
    else{
        deltaT[1] = cellDim[1]/dir[1];
        t_y = ((floor(OGrid[0]/cellDim[0]) + 1) * cellDim[1] - OGrid[1])/dir[1];
    }
    if(dir[2] < 0){
        deltaT[2] = - cellDim[2]/dir[2];
        t_z = (floor(OGrid[2]/cellDim[2]) * cellDim[2] - OGrid[2])/dir[2];
    }
    else{
        deltaT[2] = cellDim[2]/dir[2];
        t_z = ((floor(OGrid[0]/cellDim[0]) + 1) * cellDim[2] - OGrid[2])/dir[2];
    }
    //cout << " index" << index[0] << " " << index[1] << " " << index[2] << '\n';
    while(1) {
        //cout << grid.grid[3][2][4].size();
        //cout <<" index" <<  index[0] << " " << index[1] << " " << index[2] << '\n';
        //cout << "size " << grid.grid[(int)index.x][(int)index.y][(int)index.z].size() << '\n';
        rayIntersect = intersect(ray, grid.grid[(int)index.x][(int)index.y][(int)index.z]);
        if(std::min(t_x, std::min(t_y, t_z)) == t_x){
            t = t_x;
            t_x += deltaT[0];
            if(dir[0] < 0){
                index[0] = index[0] - 1;
                
            }
            else {
                index[0] += 1;
            }
        }
       else if(std::min(t_x, std::min(t_y, t_z)) == t_y){
            t = t_y;
            t_y += deltaT[1];
            if(dir[1] < 0){
                index[1] -= 1;
            }
            else {
                index[1] += 1;
            }
        }
        else if(std::min(t_x, std::min(t_y, t_z)) == t_z){
            t = t_z;
            t_z += deltaT[2];
            if(dir[2] < 0){
                index[2] -= 1;
            }
            else {
                index[2] += 1;
            }
        }
        
        if(rayIntersect.hit == true && rayIntersect.dist < t){
          //  cout << "hit " << '\n';
            break;
        }
        if(index[0] < 0 || index[1] < 0 || index[2] < 0 || index[0] > 4 || index[1] > 4 || index[2] > 4 ){
          //  cout << "end " << '\n';
            break;
        }
    }
    return rayIntersect;
}

// used in HW2, just copied
vec4 ComputeLight(vec3 direction, vec4 lightcolor, vec3 normal, vec3 halfvec, vec4 mydiffuse, vec4 myspecular, float myshininess) {
    
    float nDotL = dot(normal, direction)  ;
    vec4 lambert = mydiffuse * lightcolor * std::max(nDotL, 0.0f);
    float nDotH = dot(normal, halfvec);
    vec4 phong = myspecular * lightcolor * pow(std::max(nDotH, 0.0f), myshininess);
    vec4 retval = lambert + phong;
    return retval;
}

struct OptionalRay {
    bool valid;
    Ray ray;
};
 
OptionalRay findColor(Intersection * hit, BYTE * finalcolor ,float * ks, int depth, UGrid grid){
    if(hit->dist != INFINITY){
        object * obj = &(objects[hit->obji]);
        ks[(depth)*3] = (obj->specular)[0];
        ks[(depth)*3+1] = (obj->specular)[1];
        ks[(depth)*3+2] = (obj->specular)[2];
        vec4 col = vec4(0,0,0,0);
        vec3 color = vec3(0,0,0);
        vec3 rayDir = normalize(hit->dir);
        vec3 recDir= normalize(rayDir - 2 * dot(hit->norm,rayDir) * hit->norm);
        Ray recursiveRay = Ray(recDir,hit->P);
        for (int i = 0 ; i < numused ; i++) {
            vec3 position = vec3(lightposn[i * 4 + 0],lightposn[i * 4 + 1],lightposn[i * 4 + 2]) ;
            vec3 direction;
            float L = 1.0;
            float V = 1.0;
            if(lightposn[i * 4 + 3] == 0){
                direction = normalize(position);
                Ray rTL = Ray( direction, hit->P + (direction *(float) 0.0003));
                Intersection inshadow = DDA(rTL, grid);
                if(inshadow.dist != INFINITY){
                    V=0.0;
                }
            }else{
                direction = normalize(position - hit->P);
                float r = length(position - hit->P);
                L = 1/(attenuation.x + (attenuation.y*r) + (attenuation.z*r*r));
                Ray rTL = Ray(position - hit->P, hit->P + ((position - hit->P)*(float)0.0003));
                Intersection inshadow = DDA(rTL, grid);
                if(inshadow.dist < 1){
                    V=0.0;
                }
                
            }
            vec4 lightcolorv = vec4(lightcolor[i*4],lightcolor[i*4+1],lightcolor[i*4+2],1);
            vec3 half = normalize (direction - rayDir);
            vec4 diffusev = vec4((obj->diffuse)[0],(obj->diffuse)[1],(obj->diffuse)[2],1.0);
            vec4 specularv = vec4((obj->specular)[0],(obj->specular)[1],(obj->specular)[2],1.0);
            col += ComputeLight(direction,lightcolorv , hit->norm, half, diffusev, specularv, obj->shininess)*V*L;
        }
        
        color[0] = ((obj->ambient)[2] + (obj->emission)[2] + col[2]) *255;
        color[1] = ((obj->ambient)[1] + (obj->emission)[1] + col[1]) *255;
        color[2] = ((obj->ambient)[0] + (obj->emission)[0] + col[0]) *255;
       if(depth>0){
            for (int i=0; i < depth; i++) {
                color[0] = color[0] * ks[(i)*3 + 2];
                color[1] = color[1] * ks[(i)*3 + 1];
                color[2] = color[2] * ks[(i)*3];
            }
            finalcolor[0] += color[0];
            finalcolor[1] += color[1];
            finalcolor[2] += color[2];
            
        }else{
            finalcolor[0] = color[0];
            finalcolor[1] = color[1];
            finalcolor[2] = color[2];
        }
        OptionalRay ray;
        ray.valid = true;
        ray.ray = recursiveRay;
        return ray;
    }else{
       if(depth == 0){
            finalcolor[0] = 0;
            finalcolor[1] = 0;
            finalcolor[2] = 0;
        }
        
        OptionalRay ray;
        ray.valid = false;
        return ray;
    }
    
}
 
BYTE* rayTrace(){
    int doneVal=0;
    float ks[maxdepth*3];
    int pix = w * h;
    BYTE *image = new BYTE[3*pix];
    UGrid grid = ObjtoGrid();
    for (int i=0; i<h; i++) {
        for (int j=0; j<w; j++) {
            //depth = 0;
            Ray ray = RaythruPixel(i,j);
            Intersection intersection = DDA(ray, grid);
            for (int depth = 0; depth <= maxdepth; depth++){
                OptionalRay recursiveRay = findColor(&intersection, &image[(i*w+j)*3], ks, depth, grid);
                if(!recursiveRay.valid){
                    break;
                }
                intersection = DDA(recursiveRay.ray, grid);
            }
            //**print progress**
            int percent= (i*j/(double)pix)*100;
            if(percent%100 >= doneVal){
                printf("Done %i %%\n" , percent);
                doneVal+=2;
            }
            //******
        }
    }
    return image;
}
 
    
int main(int argc, char* argv[])
{
    maxdepth = 5;
    attenuation = vec3(1,0,0);
    FreeImage_Initialise();
    string fn = readfile(argv[1]);
    BYTE * image = rayTrace();
    saveScreenshot(fn, image);
    delete[] image;
    FreeImage_DeInitialise();
    return 0;
}

