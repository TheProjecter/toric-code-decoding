#ifndef OPERATOR_H
#define OPERATOR_H

#include <stdlib.h>
#include <string>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <assert.h>
#include <vector>


using namespace std;

//typedef
enum pauli {I, X, Y, Z};
pauli operator*(pauli p, pauli q);


class Operator
{
public:
    int size;
    vector<pauli> ops;

    Operator (int newSize);
    Operator (const Operator & p);
    void pushBack (pauli & op);
    bool commute (Operator & op);
    void generateRandom (double & p);
    void generate_random_evolution(double & p);
    pauli generateRandomPauli (double & p);
    friend Operator operator*(const Operator & p, const Operator & q);
    friend bool operator==(const Operator & p, const Operator & q);
    friend bool operator<(const Operator & p, const Operator & q);
    bool identity (void);
    void printState (void);
};


struct operatorCompare
{
    bool operator()(const Operator* p, const Operator* q) const
    {
        return *p < *q;
    }
};


#endif
