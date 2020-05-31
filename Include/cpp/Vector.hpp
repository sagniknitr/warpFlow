
#ifndef VECTOR_H
#define VECTOR_H
 
#include <vector>
 
using namespace std;
 
template <class T>
class Vector
{
protected:
    int dimension;
    T* p_data;
public:
    Vector(void);
    Vector(int dim,const T *data=NULL);
    Vector(const Vector<T>& vect);
    ~Vector(void);
 
    void init(int dim);
    void reset();
    int dim() const {return dimension;};
 
    inline T& operator[](int index){return *(p_data + index);};
};
 
template <class T>
Vector<T>::Vector(void)
{
    dimension = 0;
    p_data = NULL;
}
 
template <class T>
Vector<T>::Vector(int dim, const T *data)
{
    dimension = dim;
    p_data = new T[dimension];
    if(data != NULL)
        memcpy(p_data, data, sizeof(T)*dimension);
    else
        memset(p_data, 0, sizeof(T)*dimension);
}
 
template <class T>
Vector<T>::~Vector(void)
{
    if(p_data != NULL)
        delete[] p_data;
 
    p_data = NULL;
    dimension = 0;
}
 
template <class T>
void Vector<T>::init(int dim)
{
    if(p_data!=NULL)
            delete[] p_data;
 
    p_data = NULL;
 
    dimension = dim;
    if(dimension > 0)
    {
        p_data = new T[dimension];
        reset();
    }
}
 
template <class T>
void Vector<T>::reset()
{
    if(p_data!=NULL)
        memset(p_data,0,sizeof(T)*dimension);
}
#endif //VECTOR_H
  