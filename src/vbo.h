#pragma once

#include <glad/glad.h>


// Abstraction of opengl vertex attribute array object
class VBO {
    private:
        // private
        unsigned int ID;

    public:
        // public
        VBO(const void* data, const unsigned int size, GLenum draw_style);
        ~VBO();

        void Bind();
        void Buffer(const void* data, unsigned int size, GLenum draw_style);
        void SubBuffer(const void* data, unsigned int size, GLintptr offset);
        void Unbind();
};
