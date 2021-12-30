#pragma once

#include <glad/glad.h>


// Abstraction of opengl vertex attribute array object
class EBO {
    private:
        // private
        unsigned int ID;

    public:
        // public
        EBO(const void* data, unsigned int size, GLenum usage);
        ~EBO();

        void Bind();
        void Buffer(const void* data, unsigned int size, GLenum usage);
        void Unbind();
};
