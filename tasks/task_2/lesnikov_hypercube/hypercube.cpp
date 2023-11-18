// Copyright 2023 Kulikov Artem
#include <iostream>
#include <utility>
#include <vector>
#include <string>
#include <random>
#include <algorithm>
#include <functional>
#include <boost/mpi/collectives.hpp>
#include <boost/mpi/communicator.hpp>
#include "task_2/lesnikov_hypercube/hypercube.h"



void visualize(const std::vector<std::vector<bool>>& m) {
    for (const auto& row : m) {
        for (const int e : row) {
            std::cout << (e ? "1" : " ") << " ";
        }
        std::cout << std::endl;
    }
}

void visualizeQ(std::queue<int> q) {
    while (!q.empty()) {
        std::cout << q.front() << " ";
        q.pop();
    }
    std::cout << std::endl;
}

void visualizeV(const std::vector<int>& v) {
    for (int e : v) {
        std::cout << e << " ";
    }

    std::cout << std::endl;
}

bool is2Degree(int n) {
    int copy_n = n;

    while (copy_n % 2 == 0 && copy_n) {
        copy_n /= 2;
    }

    return static_cast<bool>(copy_n);
}

void writeIdentity(std::vector<std::vector<bool>>* m, int start_i, int end_i, int start_j) {
    for (int i = start_i, j = start_j; i < end_i; i++, j++) {
        (*m)[i][j] = true;
    }
}

void fillMatrix(std::vector<std::vector<bool>>* m, int start, int end) {
    writeIdentity(m, start, start + (end - start) / 2, start + (end - start) / 2);
    writeIdentity(m, start + (end - start) / 2, end, start);

    if (end - start != 2) {
        fillMatrix(m, start, start + (end - start) / 2);
        fillMatrix(m, start + (end - start) / 2, end);
    }
}

std::vector<std::vector<bool>> createHypercube(int n) {
    if (n == 0) {
        throw std::invalid_argument("hypercube size cannot be 0");
    } else if (!is2Degree(n)) {
        throw std::invalid_argument("hypercube size must be 2 power");
    }

    std::vector<std::vector<bool>> m(n);

    for (auto& vec : m) {
        vec = std::vector<bool>(n);
    }

    fillMatrix(&m, 0, n);

    return m;
}

std::queue<int> getPath(const std::vector<std::vector<bool>>& m, int start, int end) {
    std::unordered_set<int> labelled;

    std::queue<std::pair<int, std::queue<int>>> verticiesAndPath;

    verticiesAndPath.push(std::make_pair(start, std::queue<int>()));
    labelled.insert(start);

    while (!verticiesAndPath.empty()) {
        int size = verticiesAndPath.size();

        for (int i = 0; i < size; i++) {
            for (int j = 0; j < m[0].size(); j++) {
                if (m[verticiesAndPath.front().first][j] && labelled.find(j) == labelled.end()) {
                    labelled.insert(j);
                    std::queue<int> newPath(verticiesAndPath.front().second);
                    newPath.push(verticiesAndPath.front().first);

                    if (j == end) {
                        newPath.push(end);
                        return newPath;
                    }

                    verticiesAndPath.push(std::make_pair(j, newPath));
                }
            }

            verticiesAndPath.pop();
        }
    }

    return std::queue<int>();
}

std::pair<std::vector<int>, std::vector<int>>
getTransitionsAndExpectations(std::queue<int> path, int numProcess) {
    std::vector<int> transitions(numProcess, -1);
    std::vector<int> expectations(numProcess, -1);

    int prev;
    int cur;

    if (!path.empty()) {
        prev = path.front();
        path.pop();
    }

    while (!path.empty()) {
        cur = path.front();
        path.pop();
        transitions[prev] = cur;
        expectations[cur] = prev;
        prev = cur;
    }

    return std::make_pair(transitions, expectations);
}















std::vector<int> getRandomVector(int n) {
    std::random_device dev;
    std::mt19937 gen(dev());
    std::vector<int> random_vector;

    for (size_t i = 0; i < n; i++) {
        random_vector.push_back(gen() % 100);
    }

    return random_vector;
}


std::vector<std::pair<size_t, size_t>> getSequentialMostDifferentElements(std::vector<int> v) {
    std::vector<std::pair<size_t, size_t>> most_different_elements;
    int max_different_elements_value = 0;

    for (size_t i = 0; i < v.size() - 1; i++) {
        max_different_elements_value = std::max(max_different_elements_value, std::abs(v[i] - v[i+1]));
    }

    for (size_t i = 0; i < v.size() - 1; i++) {
        if (std::abs(v[i] - v[i+1]) == max_different_elements_value) {
            most_different_elements.push_back(std::make_pair(i, i + 1));
        }
    }

    return most_different_elements;
}

std::vector<std::pair<size_t, size_t>> getParallelMostDifferentElements(std::vector<int> v, int n) {
    boost::mpi::communicator world;
    const int useful_world_size = std::min(n / 2 + n % 2, world.size());

    if (n < 2 || useful_world_size <= world.rank()) {
        return std::vector<std::pair<size_t, size_t>>();
    }

    std::vector<size_t> temp2;
    std::vector<int> temp;
    std::vector<int> loc_v;
    std::vector<size_t> glob_most_different_elements;
    std::vector<std::pair<size_t, size_t>> glob_most_different_elements_pairs;
    std::vector<size_t> loc_most_different_elements;
    const int part_size = n / useful_world_size;
    const int remainder = n % useful_world_size;

    if (world.rank() == 0) {
        loc_v = v;
        int loc_part_size_root = useful_world_size > 1 ? part_size + remainder + 1 : part_size + remainder;
        loc_v.resize(loc_part_size_root);

        for (size_t i = 1; i < useful_world_size; i++) {
            int loc_part_size = i == useful_world_size - 1 ? part_size : part_size + 1;
            int start = remainder + i * part_size;
            int end = remainder + i * part_size + loc_part_size;
            temp = std::vector<int>(v.begin() + start, v.begin() + end);

            world.send(i, 0, temp);
        }

    } else {
        world.recv(0, 0, loc_v);
    }

    int max_different_elements_value = 0;

    for (int i = 0; i < loc_v.size() - 1; i++) {
        max_different_elements_value = std::max(max_different_elements_value, std::abs(loc_v[i] - loc_v[i+1]));
    }

    for (int i = 0; i < loc_v.size() - 1; i++) {
        if (std::abs(loc_v[i] - loc_v[i+1]) == max_different_elements_value) {
            loc_most_different_elements.push_back(i + world.rank() * part_size
            + remainder * static_cast<int>(world.rank() != 0));
        }
    }

    glob_most_different_elements = std::move(loc_most_different_elements);

    if (world.rank() == 0) {
        for (int i = 1; i < useful_world_size; i++) {
            world.recv(i, 0, temp2);

            if (!glob_most_different_elements.size()) {
                glob_most_different_elements = std::move(temp2);
            } else if (temp2.size()) {
                int glob_diff = std::abs(v[glob_most_different_elements[0]] - v[glob_most_different_elements[0] + 1]);
                int temp_diff = std::abs(v[temp2[0]] - v[temp2[0] + 1]);
                if (glob_diff < temp_diff) {
                    glob_most_different_elements = std::move(temp2);
                } else if (glob_diff == temp_diff) {
                    glob_most_different_elements.insert(
                    glob_most_different_elements.end(), temp2.begin(), temp2.end());
                }
            }
        }
    } else {
        world.send(0, 0, glob_most_different_elements);
    }

    for (const size_t& i : glob_most_different_elements) {
        glob_most_different_elements_pairs.push_back(std::make_pair(i, i + 1));
    }

    return glob_most_different_elements_pairs;
}
