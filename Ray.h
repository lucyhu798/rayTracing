#include <glm/glm.hpp>
#include "Camera.h"
class Ray {
  public:
    Ray();
    glm::vec3 p0;
    glm::vec3 p1;
    Ray(glm::vec3 center, glm::vec3 direction);
    Ray rayThruPixel(Camera cam, int i, int j);
    virtual ~Ray();
};