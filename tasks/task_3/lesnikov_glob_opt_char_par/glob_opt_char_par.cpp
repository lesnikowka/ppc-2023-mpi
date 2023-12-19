
#include "task_3/lesnikov_glob_opt_char_par/glob_opt_char_par.h"

#include <cmath>
#include <functional>
#include <iostream>
#include <utility>
#include <cmath>
#include <vector>
#include <string>
#include <random>
#include <algorithm>
#include <boost/mpi/collectives.hpp>
#include <boost/mpi/communicator.hpp>




double findM(const std::vector<double>& X, const std::function<double(double)>& f) {
	double M = 0;
	for (int i = 0; i < X.size() - 1; i++) {
	    double abs_ = std::abs(f(X[i]) - f(X[i + 1]));
	    double abs_x = std::abs(X[i + 1] - X[i]);
		M = std::max(M, abs_ / abs_x);
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
	    double abs_x = std::abs(X[i + 1] - X[i]);
	    double sq = (f(X[i]) - f(X[i + 1])) * (f(X[i]) - f(X[i + 1]));
	    double fsum = (f(X[i]) + f(X[i + 1]));
		R[i] = m * abs_x + sq / m / abs_x - 2 * fsum;
	}
	return R;
}

double getXk_1(double xt_1, double xt, double r,
	double m, std::function<double(double)> f) {

	int sign;
	if (f(xt) - f(xt_1) > 0) {
		sign = 1;
	}
	else {
		sign = -1;
	}

	double mean = ((xt_1 + xt) / 2);

	double abs_ = std::abs(f(xt) - f(xt_1));

	return mean - sign * abs_ / m / 2;
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

		//std::cout << "i=" << i << " xk=" << xk_1 << " M=" << M << " m=" << m << " maxRindex=" << maxRindex
		//<< "x: " << X[maxRindex] << "\n";

		double abs_ = std::abs(xk_1 - X[maxRindex]);

		if (abs_ < eps) {
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

    printf("world_size: %d\n\n", static_cast<int>(world.size()));

	std::vector<double> X;
    std::vector<double> loc_X;
    bool end = false;

    if (world.rank() == 0) {
        X = { leftBound, rightBound };
    }

	double lastX = 0;

	for (int i = 0; i < maxIterations; i++) {

        if (world.size() != 1) {
            boost::mpi::broadcast(world, end, 0);
        }

        if (end) {
            break;
        }

        int usefulWorldSize = 0;
        if (world.rank() == 0) {
            usefulWorldSize = std::min(static_cast<int>(world.size()), static_cast<int>(X.size()));
        }
        if (world.size() != 1) {
            boost::mpi::broadcast(world, usefulWorldSize, 0);
        }
        
        if (world.rank() >= usefulWorldSize) {
            if (usefulWorldSize != 1) {
                double dummy = 0;
                double dummy2 = 0;
                dummy = std::numeric_limits<double>::lowest();
                boost::mpi::reduce(world, dummy, dummy2, boost::mpi::maximum<double>(), 0);
                boost::mpi::broadcast(world, dummy, 0);
            }
            continue;
        }

        int part_size = 0; 
        int remainder = 0;

        if (world.rank() == 0) { 
            remainder = X.size() % usefulWorldSize; 
            part_size = X.size() / usefulWorldSize;

            printf("remain %d part_s %d x_size %d\n", remainder, part_size, (int)X.size());

            int not_end = static_cast<int>(usefulWorldSize != 1);
            loc_X = std::vector<double>(X.begin(), X.begin() + part_size + remainder + not_end);

            printf("locX root size: %d USEFUL WORLD SIZE %d\n",static_cast<int>(loc_X.size()), (int)usefulWorldSize);

            for (int j = 1; j < usefulWorldSize; j++) {
                not_end = static_cast<int>(j != usefulWorldSize - 1);
                printf("i=%d start=%d end=%d\n", j, remainder + part_size * j, remainder + part_size * (j + 1) + not_end);
                std::vector<double> temp(X.begin() + remainder + part_size * j, X.begin() + remainder + part_size * (j + 1) + not_end);
                size_t old_temp_size = temp.size();
                world.send(j, 0, temp);
                if (old_temp_size < 2) {
                    usefulWorldSize--;
                    break;
                }
            }
            
        }
        else {
            world.recv(0, 0, loc_X);
        }

        if (loc_X.size() < 2) {
            printf("PROCESS I=%d usefulsize=%d locxsize=%d\n", (int)world.rank(), (int)usefulWorldSize, (int)loc_X.size());
            if (loc_X.size() == 1) {
                usefulWorldSize--;
            }
            if (usefulWorldSize != 1) {
                double dummy = 0;
                double dummy2 = 0;
                dummy = std::numeric_limits<double>::lowest();
                boost::mpi::reduce(world, dummy, dummy2, boost::mpi::maximum<double>(), 0);
                boost::mpi::broadcast(world, dummy, 0);
            }
            continue;
        }

        printf("STEP OVER POINT 1 PROCESS I=%d\n", (int)world.rank());

		double loc_M = findM(loc_X, f);
        double M = loc_M;
        double m = 0;

        printf("STEP OVER POINT 2 PROCESS I=%d\n", (int)world.rank());

        if (usefulWorldSize != 1) {
            boost::mpi::reduce(world, loc_M, M, boost::mpi::maximum<double>(), 0);
        }
        
        printf("STEP OVER POINT 3 PROCESS I=%d\n", (int)world.rank());

		if (world.rank() == 0) {
            m = getm(M, r);
        }
        
        if (usefulWorldSize != 1) {
            boost::mpi::broadcast(world, m, 0);
        }

        printf("STEP OVER POINT 4 PROCESS I=%d\n", (int)world.rank());
        
		std::vector<double> R = getR(m, loc_X, f);
		int maxRindex = getMax(R);
        printf("maxRindex=%d Rsize=%d\n", maxRindex, (int)R.size());
        double maxRValue = R[maxRindex];

        printf("STEP OVER POINT 5 PROCESS I=%d\n", (int)world.rank());

        if (world.rank() == 0) {
            for (int j = 1; j < usefulWorldSize; j++) {
                int temp_index = 0;
                world.recv(j, 0, temp_index);
                int realIndex = temp_index + remainder + j * part_size;
                double locR = 0;
                world.recv(j, 0, locR);

                if (maxRValue < locR) {
                    locR = maxRValue;
                    maxRindex = realIndex;
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
                end = true;
                if (world.size() != 1) {
                    boost::mpi::broadcast(world, end, 0);
                }
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