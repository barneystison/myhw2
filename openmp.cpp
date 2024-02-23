#include "common.h"
#include <cmath>
#include <list>
#include <vector>
#include <stdio.h>
#include <omp.h>
#include <mutex>
typedef struct Bin {
    std::list<particle_t*> plist;
}Bin;

Bin** bin_array = NULL;
double bin_size = 0.0;
int bin_num = 0;

// Apply the force from neighbor to particle
void apply_force(particle_t& particle, particle_t& neighbor) {
    // Calculate Distance
    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;

    // Check if the two particles should interact
    if (r2 > cutoff * cutoff)
        return;

    r2 = fmax(r2, min_r * min_r);
    double r = sqrt(r2);

    // Very simple short-range repulsive force
    double coef = (1 - cutoff / r) / r2 / mass;
    particle.ax += coef * dx;
    particle.ay += coef * dy;
}

// Integrate the ODE
void move(particle_t& p, double size) {
    // Slightly simplified Velocity Verlet integration
    // Conserves energy better than explicit Euler method
    p.vx += p.ax * dt;
    p.vy += p.ay * dt;
    p.x += p.vx * dt;
    p.y += p.vy * dt;

    // Bounce from walls
    while (p.x < 0 || p.x > size) {
        p.x = p.x < 0 ? -p.x : 2 * size - p.x;
        p.vx = -p.vx;
    }

    while (p.y < 0 || p.y > size) {
        p.y = p.y < 0 ? -p.y : 2 * size - p.y;
        p.vy = -p.vy;
    }
}

void test() {
    printf("test\n");
    //printf("bin_num=%d\n", bin_num);
    for (int i = 0; i < bin_num; ++i) {
        for (int j = 0; j < bin_num; ++j) {
            printf("Bin[%d][%d] has %d particles\n", i, j, (int)bin_array[i][j].plist.size());
        }
    }
    exit(EXIT_FAILURE);
}

void init_simulation(particle_t* parts, int num_parts, double size) {
    // You can use this space to initialize static, global data objects
    // that you may need. This function will be called once before the
    // algorithm begins. Do not do any particle simulation here
    bin_size = cutoff;
    bin_num = (int)ceil(size / bin_size);
    //printf("bin_num=%d\n", bin_num);
    
    bin_array = new Bin * [bin_num];
    for (int i = 0; i < bin_num; ++i) {
        bin_array[i] = new Bin[bin_num];
    }
    //put part into Bin
    int index_row = 0, index_col = 0;

    for (int i = 0; i < num_parts; i++)
    {
        index_col = int(parts[i].x / bin_size);
        index_row = int(parts[i].y / bin_size);
        bin_array[index_row][index_col].plist.push_back(parts + i);
    }
    //test

}


void simulate_one_step(particle_t* parts, int num_parts, double size) {

    // Compute Forces
  // Compute Forces
#pragma omp for collapse(2)
    for (int i = 0; i < bin_num; ++i) {
        for (int j = 0; j < bin_num; ++j) {
            for (std::list<particle_t*>::iterator it = bin_array[i][j].plist.begin(); it != bin_array[i][j].plist.end(); ++it) {
                particle_t* p = *it;
                p->ax = 0;
                p->ay = 0;

                ////Traverse neighboring bines
                // 
                //center
                for (std::list<particle_t*>::iterator it2 = bin_array[i][j].plist.begin(); it2 != bin_array[i][j].plist.end(); ++it2) {
                    if (p != *it2) {
                        apply_force(*p, **it2);
                    }

                }
                //left up
                if (i - 1 >= 0 && j - 1 >= 0) {
                    for (std::list<particle_t*>::iterator it2 = bin_array[i - 1][j - 1].plist.begin(); it2 != bin_array[i - 1][j - 1].plist.end(); ++it2) {
                        apply_force(*p, **it2);
                    }
                }
                //left
                if (i - 1 >= 0) {
                    for (std::list<particle_t*>::iterator it2 = bin_array[i - 1][j].plist.begin(); it2 != bin_array[i - 1][j].plist.end(); ++it2) {
                        apply_force(*p, **it2);
                    }
                }
                //left down
                if (i - 1 >= 0 && j + 1 < bin_num) {
                    for (std::list<particle_t*>::iterator it2 = bin_array[i - 1][j + 1].plist.begin(); it2 != bin_array[i - 1][j + 1].plist.end(); ++it2) {
                        apply_force(*p, **it2);
                    }
                }
                //up
                if (j - 1 >= 0) {
                    for (std::list<particle_t*>::iterator it2 = bin_array[i][j - 1].plist.begin(); it2 != bin_array[i][j - 1].plist.end(); ++it2) {
                        apply_force(*p, **it2);
                    }
                }
                //down
                if (j + 1 < bin_num) {
                    for (std::list<particle_t*>::iterator it2 = bin_array[i][j + 1].plist.begin(); it2 != bin_array[i][j + 1].plist.end(); ++it2) {
                        apply_force(*p, **it2);
                    }
                }
                //right up
                if (i + 1 < bin_num && j - 1 >= 0) {
                    for (std::list<particle_t*>::iterator it2 = bin_array[i + 1][j - 1].plist.begin(); it2 != bin_array[i + 1][j - 1].plist.end(); ++it2) {
                        apply_force(*p, **it2);
                    }
                }
                //right
                if (i + 1 < bin_num) {
                    for (std::list<particle_t*>::iterator it2 = bin_array[i + 1][j].plist.begin(); it2 != bin_array[i + 1][j].plist.end(); ++it2) {
                        apply_force(*p, **it2);
                    }
                }
                //right down
                if (i + 1 < bin_num && j + 1 < bin_num) {
                    for (std::list<particle_t*>::iterator it2 = bin_array[i + 1][j + 1].plist.begin(); it2 != bin_array[i + 1][j + 1].plist.end(); ++it2) {
                        apply_force(*p, **it2);
                    }
                }
            }
        }
    }
#pragma omp for
    // Move Particles
    for (int i = 0; i < num_parts; ++i) {
        move(parts[i], size);

    }
#pragma omp for collapse(2)
    for (int i = 0; i < bin_num; ++i) {
        for (int j = 0; j < bin_num; ++j) {
            for (std::list<particle_t*>::iterator it = bin_array[i][j].plist.begin(); it != bin_array[i][j].plist.end();) {
                particle_t* p = *it;
                // move(*p, size);
                int index_row = int(p->y / bin_size);
                int index_col = int(p->x / bin_size);
                if (index_row != i || index_col != j) {
                    std::list<particle_t*>::iterator nextIt = std::next(it);
                    bin_array[i][j].plist.erase(it);
                    bin_array[index_row][index_col].plist.push_back(p);
                    it = nextIt;
                }
                else
                {
                    ++it;
                }
            }
        }
    }

}
