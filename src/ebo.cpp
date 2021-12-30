#include <glad/glad.h>

#include "ebo.h"


// Abstraction of opengl vertex attribute array object
EBO::EBO(const void* data, unsigned int size, GLenum usage)
{
    glGenBuffers(1, &ID);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ID);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, size, data, usage);
}

EBO::~EBO()
{
    glDeleteBuffers(1, &ID);
}

void EBO::Bind()
{
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ID);
}

void EBO::Buffer(const void* data, unsigned int size, GLenum usage)
{
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, size, data, usage);
}

void EBO::Unbind()
{
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}
