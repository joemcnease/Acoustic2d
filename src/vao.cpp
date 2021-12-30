#include <glad/glad.h>

#include "vao.h"


// Abstraction of opengl vertex attribute array object
VAO::VAO()
{
    glGenVertexArrays(1, &ID);
    glBindVertexArray(ID);
}

VAO::~VAO()
{
    glDeleteVertexArrays(1, &ID);
}

void VAO::Bind()
{
    glBindVertexArray(ID);
}

void VAO::Unbind()
{
    glBindVertexArray(0);
};

void VAO::VertexAttribPointer(GLuint index, GLint size, GLenum type, GLboolean normalized, GLsizei stride, const void* pointer)
{
    glVertexAttribPointer(index, size, type, normalized, stride, pointer);
    glEnableVertexAttribArray(index);
};
