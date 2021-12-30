#pragma once

#include <glad/glad.h>


// Abstraction of opengl vertex attribute array object
class VAO {
    private:
        unsigned int ID;

    public:
        VAO();
        ~VAO();

        void Bind();
        void Unbind();
        void VertexAttribPointer(GLuint index, GLint size, GLenum type, GLboolean normalized, GLsizei stride, const void* pointer);
};
