// Copyright 2023 Kulikov Artem
#ifndef TASKS_TASK_2_LESNIKOV_HYPERCUBE_HYPERCUBE_H_
#define TASKS_TASK_2_LESNIKOV_HYPERCUBE_HYPERCUBE_H_

#include <utility>
#include <algorithm>
#include <vector>
#include <queue>
#include <unordered_set>
#include <exception>
#include <boost/serialization/vector.hpp>

void visualize(const std::vector<std::vector<bool>>& m);
void visualizeQ(std::queue<int> q);
bool is2Degree(int n);
void writeIdentity(std::vector<std::vector<bool>>* m, int start_i, int end_i, int start_j);
void fillMatrix(std::vector<std::vector<bool>>* m, int start, int end);
std::vector<std::vector<bool>> createHypercube(int n);
std::queue<int> getPath(const std::vector<std::vector<bool>>& m, int start, int end);
std::pair<std::vector<int>, std::vector<int>> getTransitionsAndExpectations(std::queue<int> path, int numProcess);

template<class T>
void sendData(int source, int dest, int tag, const T& value) {
    boost::mpi::communicator world;

    std::vector<std::vector<bool>> hypercube = createHypercube(world.size());
    std::pair<std::vector<int>, std::vector<int>> transitionsAndExpectations =
    getTransitionsAndExpectations(getPath(hypercube, source, dest), world.size());

    T data;

    world.recv(transitionsAndExpectations.second[world.rank()], tag, data);
    world.send(transitionsAndExpectations.first[world.rank()], tag, data);
}

#endif  // TASKS_TASK_2_LESNIKOV_HYPERCUBE_HYPERCUBE_H_
