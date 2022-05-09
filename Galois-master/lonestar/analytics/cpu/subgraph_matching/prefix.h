#pragma once
#include <cstring>
#include <vector>

class Prefix {
public:
    Prefix();
    ~Prefix();
    void init(int input_size, const int* input_data);
    bool operator ==(const Prefix& pre) const;
    bool operator !=(const Prefix& pre) const;
    bool equal(int input_size, const int* input_data) const;
    inline int get_size() const { return size;}
    inline const int* get_data_ptr() const { return data;}
    inline int get_data(int index) const { return data[index];}

private:
    int size;
    int* data;
};

Prefix::Prefix() {
    data = nullptr;
}

Prefix::~Prefix() {
    // if (data != nullptr)
    //     delete[] data;
}

void Prefix::init(int input_size, const int* input_data) {
    size = input_size;
    data = new int[input_size];
    memcpy(data, input_data, size * sizeof(int));
}

bool Prefix::operator==(const Prefix& pre) const {
    if (size != pre.get_size())
        return false;
    for (int i = 0; i < size; ++i)
        if (data[i] != pre.get_data(i))
            return false;
    return true;
}

bool Prefix::operator!=(const Prefix& pre) const {
    if (size != pre.get_size())
        return true;
    for (int i = 0; i < size; ++i)
        if (data[i] != pre.get_data(i))
            return true;
    return false;
}

bool Prefix::equal(int input_size, const int* input_data) const {
    if (size != input_size)
        return false;
    for (int i = 0; i < size; ++i)
        if (data[i] != input_data[i])
            return false;
    return true;
}