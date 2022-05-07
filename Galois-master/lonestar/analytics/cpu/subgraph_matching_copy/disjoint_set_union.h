#pragma once

class DisjointSetUnion {
public:
    DisjointSetUnion(int n);
    ~DisjointSetUnion();
    void init();
    void merge(int id1, int id2);
    inline int get_set_size() const { return set_size;}
    inline int get_size() const { return size;}

private:
    int get_father(int a);
    
    int size;
    int set_size;
    int* father;
};
DisjointSetUnion::DisjointSetUnion(int n) {
    size = n;
    set_size = n;
    father = new int[n];
}

DisjointSetUnion::~DisjointSetUnion() {
    if( size > 0) {
        delete[] father;
    }
}

void DisjointSetUnion::init() {
    for(int i = 0; i < size; ++i)
        father[i] = i;
    set_size = size;
}

void DisjointSetUnion::merge(int id1, int id2) {
    int fa1 = get_father(id1);
    int fa2 = get_father(id2);
    if(fa1 != fa2) {
        --set_size;
        father[fa1] = fa2;
    }
}

int DisjointSetUnion::get_father(int a) {
    if( father[a] == a) return a;
    else return father[a] = get_father(father[a]);
}
