// Transform.cpp: implementation of the Transform class.

// Note: when you construct a matrix using mat4() or mat3(), it will be COLUMN-MAJOR
// Keep this in mind in readfile.cpp and display.cpp
// See FAQ for more details or if you're having problems.

#include "Transform.h"

// Helper rotation function.  Please implement this.  
mat3 Transform::rotate(const float degrees, const vec3& axis) 
{
    float radians =  glm::radians(degrees);
      vec3 a = glm::normalize(axis);
      // make this a lot shorter if possible
      float cosine = cos(radians);
      float sine = sin(radians);
    vec3 row1 = vec3(cosine + (a.x*a.x) * (1 - cosine),(a.x*a.y) * (1-cosine) - (a.z*sine),(a.x*a.z) * (1-cosine) + (a.y*sine));
    vec3 row2 = vec3((a.y*a.x) * (1-cosine) + (a.z*sine),cosine + (a.y*a.y) * (1 - cosine),(a.y*a.z) * (1-cosine) - (a.x*sine));
    vec3 row3 = vec3((a.z*a.x) * (1-cosine) - (a.y*sine),(a.z*a.y) * (1-cosine) + (a.x*sine),cosine + (a.z*a.z) * (1 - cosine));
        
      return glm::transpose(mat3(row1,row2,row3));
}

void Transform::left(float degrees, vec3& eye, vec3& up) 
{
  // YOUR CODE FOR HW2 HERE
  // Likely the same as in HW 1.
    mat3 rotateMatrix = rotate(degrees, up);
        eye = rotateMatrix * eye;
}

void Transform::up(float degrees, vec3& eye, vec3& up) 
{
    vec3 cross = glm::cross(up, eye);
    mat3 rotateMatrix = rotate(-degrees, cross);
    eye = rotateMatrix * eye;
    up = rotateMatrix * up;
}

mat4 Transform::lookAt(const vec3 &eye, const vec3 &center, const vec3 &up) 
{
    vec3 c = glm::normalize(eye);
    vec3 b = glm::normalize(glm::cross(up,eye));
    vec3 a = glm::normalize(glm::cross(c,b));
    
    vec4 row2 = vec4(a.x,a.y,a.z,glm::dot(a, -eye));
    vec4 row1 = vec4(b.x,b.y,b.z,glm::dot(b, -eye));
    vec4 row3 = vec4(c.x,c.y,c.z,glm::dot(c, -eye));
    vec4 row4 = vec4(0,0,0,1);
    
    return glm::transpose(mat4(row1, row2, row3, row4));
}

mat4 Transform::perspective(float fovy, float aspect, float zNear, float zFar)
{
  mat4 ret;
    float d = 1/(tan(glm::radians(fovy)/2));
    float A = -(zFar + zNear)/(zFar - zNear);
    float B = -(2 * zFar * zNear)/(zFar - zNear);
    ret = mat4(d/aspect, 0, 0, 0,
               0, d, 0, 0,
               0, 0, A , -1,
               0, 0, B, 0);
  // YOUR CODE FOR HW2 HERE
  // New, to implement the perspective transform as well.  
  return ret;
}

mat4 Transform::scale(const float &sx, const float &sy, const float &sz) 
{
  mat4 ret;
  ret = mat4(sx, 0, 0, 0,
             0, sy, 0, 0,
             0, 0, sz, 0,
             0, 0, 0, 1);
  // Implement scaling
  return ret;
}

mat4 Transform::translate(const float &tx, const float &ty, const float &tz) 
{
  mat4 ret;
  ret = mat4(1, 0, 0, 0,
             0, 1, 0, 0,
             0, 0, 1, 0,
             tx, ty, tz, 1);
  // YOUR CODE FOR HW2 HERE
  // Implement translation 
  return ret;
}

// To normalize the up direction and construct a coordinate frame.  
// As discussed in the lecture.  May be relevant to create a properly 
// orthogonal and normalized up. 
// This function is provided as a helper, in case you want to use it. 
// Using this function (in readfile.cpp or display.cpp) is optional.  

vec3 Transform::upvector(const vec3 &up, const vec3 & zvec) 
{
  vec3 x = glm::cross(up,zvec); 
  vec3 y = glm::cross(zvec,x); 
  vec3 ret = glm::normalize(y); 
  return ret; 
}


Transform::Transform()
{

}

Transform::~Transform()
{

}
