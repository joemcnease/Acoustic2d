#include <glad/glad.h>

#include "vbo.h"


// Abstraction of opengl vertex attribute array object
VBO::VBO(const void* data, unsigned int size, GLenum draw_style)
{
    glGenBuffers(1, &ID);
    glBindBuffer(GL_ARRAY_BUFFER, ID);
    glBufferData(GL_ARRAY_BUFFER, size, data, draw_style);
}

VBO::~VBO()
{
    glDeleteBuffers(1, &ID);
}

void VBO::Bind()
{
    glBindBuffer(GL_ARRAY_BUFFER, ID);
};

void VBO::Buffer(const void* data, unsigned int size, GLenum draw_style)
{
    glBufferData(GL_ARRAY_BUFFER, size, data, draw_style);
}

void VBO::SubBuffer(const void* data, unsigned int size, GLintptr offset)
{
    glBindBuffer(GL_ARRAY_BUFFER, ID);
    glBufferSubData(GL_ARRAY_BUFFER, offset, size, data);
}

void VBO::Unbind()
{
    glBindBuffer(GL_ARRAY_BUFFER, 0);
};
