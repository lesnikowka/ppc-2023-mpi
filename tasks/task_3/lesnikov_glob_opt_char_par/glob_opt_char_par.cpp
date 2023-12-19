
#include "task_3/lesnikov_glob_opt_char_par/glob_opt_char_par.h"

#include <cmath>
#include <functional>
#include <iostream>
#include <utility>
#include <vector>
#include <string>
#include <random>
#include <algorithm>
#include <boost/mpi/collectives.hpp>
#include <boost/mpi/communicator.hpp>



double findM(const std::vector<double>& X, const std::function<double(double)>& f) {
	double M = 0;
	for (int i = 0; i < X.size() - 1; i++) {
		M = std::max(M, std::abs(f(X[i]) - f(X[i + 1])) / std::abs(X[i + 1] - X[i]));
	}
	return M;
}

double getm(double M, double r) {
	if (M) {
		return r * M;
	}
	return 1;
}

std::vector<double> getR(double m, const std::vector<double>& X, std::function<double(double)> f) {
	std::vector<double> R(X.size() - 1);
	for (int i = 0; i < X.size() - 1; i++) {
		R[i] = m * abs(X[i + 1] - X[i]) + (f(X[i]) - f(X[i + 1])) * (f(X[i]) - f(X[i + 1]))
			/ m / abs(X[i + 1] - X[i]) - 2 * (f(X[i]) + f(X[i + 1]));
	}
	return R;
}

double getXk_1(double xt_1, double xt, double r,
	double m, std::function<double(double)> f) {

	return (xt_1 + xt) / 2 - ((f(xt) - f(xt_1)) > 0 ? 1 : -1) * r / m * std::abs(f(xt) - f(xt_1)) / (2*r);
}

int getMax(const std::vector<double>& v) {
	int maxIndex = 0;
	for (int i = 0; i < v.size(); i++) {
		if (v[maxIndex] < v[i]) {
			maxIndex = i;
		}
	}
	return maxIndex;
}

double getMinSequential(std::function<double(double)> f, double leftBound,
	double rightBound, double eps, int maxIterations, double r) {

	std::vector<double> X = { leftBound, rightBound };

	double lastX = 0;

	for (int i = 0; i < maxIterations; i++) {
		double M = findM(X, f);
		double m = getm(M, r);
		std::vector<double> R = getR(m, X, f);
		int maxRindex = getMax(R);

		double xk_1 = getXk_1(X[maxRindex], X[maxRindex + 1], r, m, f);

		if (xk_1 < leftBound) {
			xk_1 = leftBound;
		}
		else if (xk_1 > rightBound) {
			xk_1 = rightBound;
		}

		lastX = xk_1;


		if (std::abs(xk_1 - X[maxRindex]) < eps) {
			return xk_1;
		}

		auto xk_1Place = std::upper_bound(X.begin(), X.end(), xk_1);

		X.insert(xk_1Place, xk_1);
	}

	return lastX;
}

double getMinParallel(std::function<double(double)> f, double leftBound,
	double rightBound, double eps, int maxIterations, double r) {
    
    boost::mpi::communicator world;

	std::vector<double> X;
    std::vector<double> loc_X;
    if (world.size() == 0) {
        X = { leftBound, rightBound };
    }

	double lastX = 0;

	for (int i = 0; i < maxIterations; i++) {
        int usefulWorldSize;
        if (world.rank == 0) {
            usefulWorldSize = std::min(world.size(), X.size());
        }
        boost::mpi::broadcast(world, usefulWorldSize, 0);

        if (world.rank() >= usefulWorldSize) {
            continue;
        }

        const int part_size; 
        const int remainder;

        if (world.rank() == 0) { 
            remainder = X.size() % usefulWorldSize; 
            part_size = X.size() / usefulWorldSize;

            loc_X.resize(part_size + remainder);

            for (int i = 1; i < usefulWorldSize; i++) {
                int not_end = static_cast<int>(i != world.size() - 1);
                std::vector<double> temp(X.begin() + remainder + part_size * i, X.begin() + remainder + part_size * (i + 1) + not_end);
                world.send(i, 0, temp);
            }
        }
        else {
            world.recv(0, 0, loc_X);
        }

		double loc_M = findM(loc_X, f);
        double M;
        double m;

        boost::mpi::reduce(world, loc_M, boost::mpi::maximum<double>(), M, 0);

		if (world.rank() == 0) {
            m = getm(M, r);
        }
        
        boost::mpi::broadcast(world, m, 0);
        
		std::vector<double> R = getR(m, loc_X, f);
		int maxRindex = getMax(R);
        double maxRValue = R[maxRindex];

        if (world.rank() == 0) {
            for (int i = 1; i < usefulWorldSize; i++) {
                int temp_index;
                world.recv(i, 0, temp_index);
                int realIndex = temp_index + remainder + i * part_size;
                double locR;
                world.recv(i, 0, locR);

                if (maxRValue < locR) {
                    locR = maxRValue;
                    maxRindex = real_index;
                }
            }

            double xk_1 = getXk_1(X[maxRindex], X[maxRindex + 1], r, m, f);

            if (xk_1 < leftBound) {
                xk_1 = leftBound;
            }
            else if (xk_1 > rightBound) {
                xk_1 = rightBound;
            }

            lastX = xk_1;

            if (std::abs(xk_1 - X[maxRindex]) < eps) {
                return xk_1;
            }

            auto xk_1Place = std::upper_bound(X.begin(), X.end(), xk_1);

            X.insert(xk_1Place, xk_1);
        }
        else {
            world.send(0, 0, maxRindex);
            world.send(0, 0, maxRValue);
        }
	}

	return lastX;
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