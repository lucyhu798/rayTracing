#include <glm/glm.hpp>
class Camera {
    public:
        glm::vec3 lookFrom; 
        glm::vec3 lookAt; 
        glm::vec3 upp; 
        float fovyy; 
        Camera();
        Camera(glm::vec3 eyew, glm::vec3 centerw, glm::vec3 upw, float fovyw);
        virtual ~Camera();
};