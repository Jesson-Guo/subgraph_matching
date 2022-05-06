#pragma once
#include <cstring>
#include <vector>

class Prefix {
public:
    Prefix() {};
    ~Prefix() {};
    void init(int input_size, std::vector<int> input_data);
    bool operator==(const Prefix &pre) const;
    bool operator!=(const Prefix &pre) const;
    bool equal(int input_size, std::vector<int> input_data) const;
    inline int get_size() const { return size; }
    inline std::vector<int> get_data_ptr() const { return data; }
    inline int get_data(int index) const { return data[index]; }

private:
    int size;
    std::vector<int> data;
};

void Prefix::init(int input_size, std::vector<int> input_data) {
    size = input_size;
    data = input_data;
}

bool Prefix::operator==(const Prefix &pre) const {
    if (size != pre.get_size())
        return false;
    for (int i = 0; i < size; ++i)
        if (data[i] != pre.get_data(i))
            return false;
    return true;
}

bool Prefix::operator!=(const Prefix &pre) const {
    if (size != pre.get_size())
        return true;
    for (int i = 0; i < size; ++i)
        if (data[i] != pre.get_data(i))
            return true;
    return false;
}

bool Prefix::equal(int input_size, std::vector<int> input_data) const {
    if (size != input_size)
        return false;
    for (int i = 0; i < size; ++i)
        if (data[i] != input_data[i])
            return false;
    return true;
}