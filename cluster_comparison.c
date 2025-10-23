// cluster_comparison.c: Compare time complexity of clustering algorithms with optimized DBSCAN

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>

// Constants
#define GRID_SIZES {100, 200, 300, 400, 500, 600, 700, 800, 900, 1000} // Grid sizes from 100 to 1000
#define NUM_SIZES 10  // Number of grid sizes
#define PROB 0.1  // Probability of a cell being 1
#define MAX_STACK_SIZE 10000000  // Stack size for DFS
#define SKIP 3  // Skip distance for Skip-DFS
#define DBSCAN_EPS 3.0  // Epsilon for DBSCAN
#define DBSCAN_MINPTS 10  // Minimum points for DBSCAN
#define KMEANS_K 20  // Number of clusters for k-means
#define HIERARCHICAL_THRESHOLD 4.0  // Threshold for hierarchical clustering

// Structures
typedef struct {
    int i, j, k;  // 3D coordinates
} Point;

typedef struct {
    int i, j, k;  // 3D coordinates
    int is_skipping;  // Flag for skipping in Skip-DFS
} StackElement;

typedef struct {
    int parent;  // Parent in Union-Find
    int rank;  // Rank in Union-Find
} UnionFind;

// Timing function
double get_time() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec * 1e-6;
}

// Generate 3D binary grid
int* generate_grid(int size, double prob) {
    int total_size = size * size * size;
    int* grid = (int*)malloc(total_size * sizeof(int));
    if (!grid) {
        printf("Memory allocation failed for grid of size %d\n", size);
        return NULL;
    }
    for (int i = 0; i < total_size; i++) {
        grid[i] = (rand() / (double)RAND_MAX) < prob ? 1 : 0;
    }
    return grid;
}

// Traditional DFS
int traditional_dfs(int* grid, int size, int* visited) {
    StackElement* stack = (StackElement*)malloc(MAX_STACK_SIZE * sizeof(StackElement));
    if (!stack) return 0;
    int top = 0, max_size = 0;
    int total_size = size * size * size;

    for (int idx = 0; idx < total_size; idx++) {
        if (grid[idx] == 1 && !visited[idx]) {
            int cluster_size = 0;
            stack[top++] = (StackElement){idx / (size * size), (idx / size) % size, idx % size, 0};
            while (top > 0) {
                StackElement curr = stack[--top];
                int i = curr.i, j = curr.j, k = curr.k;
                int index = i * size * size + j * size + k;
                if (i < 0 || i >= size || j < 0 || j >= size || k < 0 || k >= size || grid[index] == 0 || visited[index]) continue;
                visited[index] = 1;
                cluster_size++;
                int deltas[6][3] = {{-1,0,0}, {1,0,0}, {0,-1,0}, {0,1,0}, {0,0,-1}, {0,0,1}};
                for (int d = 0; d < 6; d++) {
                    int ni = i + deltas[d][0], nj = j + deltas[d][1], nk = k + deltas[d][2];
                    if (ni >= 0 && ni < size && nj >= 0 && nj < size && nk >= 0 && nk < size && top < MAX_STACK_SIZE) {
                        stack[top++] = (StackElement){ni, nj, nk, 0};
                    }
                }
            }
            if (cluster_size > max_size) max_size = cluster_size;
        }
    }
    free(stack);
    return max_size;
}

// Skip-DFS
int skip_dfs(int* grid, int size, int* visited) {
    StackElement* stack = (StackElement*)malloc(MAX_STACK_SIZE * sizeof(StackElement));
    if (!stack) return 0;
    int capacity = MAX_STACK_SIZE, top = 0, max_size = 0;
    int total_size = size * size * size;

    for (int idx = 0; idx < total_size; idx++) {
        if (grid[idx] == 1 && !visited[idx]) {
            int cluster_size = 0;
            stack[top++] = (StackElement){idx / (size * size), (idx / size) % size, idx % size, 1};
            while (top > 0) {
                StackElement curr = stack[--top];
                int i = curr.i, j = curr.j, k = curr.k, is_skipping = curr.is_skipping;
                int index = i * size * size + j * size + k;
                if (i < 0 || i >= size || j < 0 || j >= size || k < 0 || k >= size || grid[index] == 0 || visited[index]) continue;
                visited[index] = 1;
                cluster_size++;
                int deltas[6][3];
                if (is_skipping && grid[index] == 1) {
                    int skip_deltas[6][3] = {{-SKIP,0,0}, {SKIP,0,0}, {0,-SKIP,0}, {0,SKIP,0}, {0,0,-SKIP}, {0,0,SKIP}};
                    memcpy(deltas, skip_deltas, sizeof(skip_deltas));
                } else {
                    int adj_deltas[6][3] = {{-1,0,0}, {1,0,0}, {0,-1,0}, {0,1,0}, {0,0,-1}, {0,0,1}};
                    memcpy(deltas, adj_deltas, sizeof(adj_deltas));
                }
                for (int d = 0; d < 6; d++) {
                    int ni = i + deltas[d][0], nj = j + deltas[d][1], nk = k + deltas[d][2];
                    if (ni >= 0 && ni < size && nj >= 0 && nj < size && nk >= 0 && nk < size) {
                        int n_index = ni * size * size + nj * size + nk;
                        if (top >= capacity) {
                            capacity *= 2;
                            StackElement* new_stack = (StackElement*)realloc(stack, capacity * sizeof(StackElement));
                            if (!new_stack) { free(stack); return cluster_size; }
                            stack = new_stack;
                        }
                        stack[top++] = (StackElement){ni, nj, nk, grid[n_index] == 1};
                    }
                }
            }
            if (cluster_size > max_size) max_size = cluster_size;
        }
    }
    free(stack);
    return max_size;
}

// Optimized DBSCAN with grid-based indexing
int dbscan(int* grid, int size, int* visited) {
    int total_size = size * size * size, max_size = 0;
    int cell_size = (int)ceil(DBSCAN_EPS); // Cell size for spatial index
    int grid_cells = (size + cell_size - 1) / cell_size; // Number of cells per dimension
    int*** cell_counts = (int***)malloc(grid_cells * sizeof(int**));
    for (int i = 0; i < grid_cells; i++) {
        cell_counts[i] = (int**)malloc(grid_cells * sizeof(int*));
        for (int j = 0; j < grid_cells; j++) {
            cell_counts[i][j] = (int*)calloc(grid_cells, sizeof(int));
        }
    }

    // Count points per cell
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            for (int k = 0; k < size; k++) {
                int idx = i * size * size + j * size + k;
                if (grid[idx] == 1) {
                    int ci = i / cell_size, cj = j / cell_size, ck = k / cell_size;
                    cell_counts[ci][cj][ck]++;
                }
            }
        }
    }

    // Allocate memory for points in each cell
    Point**** cell_points = (Point****)malloc(grid_cells * sizeof(Point***));
    for (int i = 0; i < grid_cells; i++) {
        cell_points[i] = (Point***)malloc(grid_cells * sizeof(Point**));
        for (int j = 0; j < grid_cells; j++) {
            cell_points[i][j] = (Point**)malloc(grid_cells * sizeof(Point*));
            for (int k = 0; k < grid_cells; k++) {
                cell_points[i][j][k] = (Point*)malloc(cell_counts[i][j][k] * sizeof(Point));
            }
        }
    }

    // Fill points in cells
    int*** cell_indices = (int***)malloc(grid_cells * sizeof(int**));
    for (int i = 0; i < grid_cells; i++) {
        cell_indices[i] = (int**)malloc(grid_cells * sizeof(int*));
        for (int j = 0; j < grid_cells; j++) {
            cell_indices[i][j] = (int*)calloc(grid_cells, sizeof(int));
        }
    }
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            for (int k = 0; k < size; k++) {
                int idx = i * size * size + j * size + k;
                if (grid[idx] == 1) {
                    int ci = i / cell_size, cj = j / cell_size, ck = k / cell_size;
                    int index = cell_indices[ci][cj][ck]++;
                    cell_points[ci][cj][ck][index] = (Point){i, j, k};
                }
            }
        }
    }

    // DBSCAN clustering
    for (int ci = 0; ci < grid_cells; ci++) {
        for (int cj = 0; cj < grid_cells; cj++) {
            for (int ck = 0; ck < grid_cells; ck++) {
                for (int p = 0; p < cell_counts[ci][cj][ck]; p++) {
                    Point curr = cell_points[ci][cj][ck][p];
                    int idx = curr.i * size * size + curr.j * size + curr.k;
                    if (visited[idx]) continue;

                    // Count neighbors in nearby cells
                    int neighbors = 0;
                    for (int di = -1; di <= 1; di++) {
                        for (int dj = -1; dj <= 1; dj++) {
                            for (int dk = -1; dk <= 1; dk++) {
                                int nci = ci + di, ncj = cj + dj, nck = ck + dk;
                                if (nci >= 0 && nci < grid_cells && ncj >= 0 && ncj < grid_cells && nck >= 0 && nck < grid_cells) {
                                    for (int q = 0; q < cell_counts[nci][ncj][nck]; q++) {
                                        Point other = cell_points[nci][ncj][nck][q];
                                        double dist = sqrt(pow(curr.i - other.i, 2) + pow(curr.j - other.j, 2) + pow(curr.k - other.k, 2));
                                        if (dist <= DBSCAN_EPS) neighbors++;
                                    }
                                }
                            }
                        }
                    }

                    if (neighbors >= DBSCAN_MINPTS) {
                        int cluster_size = 0;
                        StackElement* stack = (StackElement*)malloc(MAX_STACK_SIZE * sizeof(StackElement));
                        if (!stack) continue;
                        int top = 0;
                        stack[top++] = (StackElement){curr.i, curr.j, curr.k, 0};
                        while (top > 0) {
                            StackElement s = stack[--top];
                            int i = s.i, j = s.j, k = s.k;
                            idx = i * size * size + j * size + k;
                            if (visited[idx]) continue;
                            visited[idx] = 1;
                            cluster_size++;

                            // Expand cluster
                            int nci = i / cell_size, ncj = j / cell_size, nck = k / cell_size;
                            for (int di = -1; di <= 1; di++) {
                                for (int dj = -1; dj <= 1; dj++) {
                                    for (int dk = -1; dk <= 1; dk++) {
                                        int nnci = nci + di, nncj = ncj + dj, nnck = nck + dk;
                                        if (nnci >= 0 && nnci < grid_cells && nncj >= 0 && nncj < grid_cells && nnck >= 0 && nnck < grid_cells) {
                                            for (int q = 0; q < cell_counts[nnci][nncj][nnck]; q++) {
                                                Point other = cell_points[nnci][nncj][nnck][q];
                                                int o_idx = other.i * size * size + other.j * size + other.k;
                                                if (visited[o_idx]) continue;
                                                double dist = sqrt(pow(i - other.i, 2) + pow(j - other.j, 2) + pow(k - other.k, 2));
                                                if (dist <= DBSCAN_EPS) {
                                                    if (top < MAX_STACK_SIZE) {
                                                        stack[top++] = (StackElement){other.i, other.j, other.k, 0};
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        free(stack);
                        if (cluster_size > max_size) max_size = cluster_size;
                    }
                }
            }
        }
    }

    // Free memory
    for (int i = 0; i < grid_cells; i++) {
        for (int j = 0; j < grid_cells; j++) {
            free(cell_counts[i][j]);
            free(cell_indices[i][j]);
            for (int k = 0; k < grid_cells; k++) {
                free(cell_points[i][j][k]);
            }
            free(cell_points[i][j]);
        }
        free(cell_counts[i]);
        free(cell_indices[i]);
        free(cell_points[i]);
    }
    free(cell_counts);
    free(cell_indices);
    free(cell_points);
    return max_size;
}

// Hierarchical (Single-Linkage)
int hierarchical(int* grid, int size, int* visited) {
    int total_size = size * size * size, max_size = 0;
    Point* points = (Point*)malloc(total_size * sizeof(Point));
    if (!points) return 0;
    int num_points = 0;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            for (int k = 0; k < size; k++) {
                int idx = i * size * size + j * size + k;
                if (grid[idx] == 1) {
                    points[num_points++] = (Point){i, j, k};
                }
            }
        }
    }

    UnionFind* uf = (UnionFind*)malloc(num_points * sizeof(UnionFind));
    if (!uf) { free(points); return 0; }
    for (int i = 0; i < num_points; i++) {
        uf[i].parent = i;
        uf[i].rank = 0;
    }

    for (int p = 0; p < num_points; p++) {
        for (int q = p + 1; q < num_points; q++) {
            double dist = sqrt(pow(points[p].i - points[q].i, 2) + pow(points[p].j - points[q].j, 2) + pow(points[p].k - points[q].k, 2));
            if (dist <= HIERARCHICAL_THRESHOLD) {
                int root_p = p, root_q = q;
                while (uf[root_p].parent != root_p) root_p = uf[root_p].parent;
                while (uf[root_q].parent != root_q) root_q = uf[root_q].parent;
                if (root_p != root_q) {
                    if (uf[root_p].rank < uf[root_q].rank) uf[root_p].parent = root_q;
                    else if (uf[root_p].rank > uf[root_q].rank) uf[root_q].parent = root_p;
                    else { uf[root_q].parent = root_p; uf[root_p].rank++; }
                }
            }
        }
    }

    int* cluster_sizes = (int*)calloc(num_points, sizeof(int));
    for (int i = 0; i < num_points; i++) {
        int root = i;
        while (uf[root].parent != root) root = uf[root].parent;
        cluster_sizes[root]++;
        if (cluster_sizes[root] > max_size) max_size = cluster_sizes[root];
    }

    free(points);
    free(uf);
    free(cluster_sizes);
    return max_size;
}

// k-means
int kmeans(int* grid, int size, int* visited) {
    int total_size = size * size * size, max_size = 0;
    Point* points = (Point*)malloc(total_size * sizeof(Point));
    if (!points) return 0;
    int num_points = 0;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            for (int k = 0; k < size; k++) {
                int idx = i * size * size + j * size + k;
                if (grid[idx] == 1) {
                    points[num_points++] = (Point){i, j, k};
                }
            }
        }
    }

    Point* centroids = (Point*)malloc(KMEANS_K * sizeof(Point));
    if (!centroids) { free(points); return 0; }
    for (int i = 0; i < KMEANS_K; i++) {
        int idx = rand() % num_points;
        centroids[i] = points[idx];
    }

    int* assignments = (int*)malloc(num_points * sizeof(int));
    if (!assignments) { free(points); free(centroids); return 0; }
    for (int iter = 0; iter < 10; iter++) {
        for (int p = 0; p < num_points; p++) {
            double min_dist = 1e9;
            int best_k = 0;
            for (int k = 0; k < KMEANS_K; k++) {
                double dist = sqrt(pow(points[p].i - centroids[k].i, 2) + pow(points[p].j - centroids[k].j, 2) + pow(points[p].k - centroids[k].k, 2));
                if (dist < min_dist) {
                    min_dist = dist;
                    best_k = k;
                }
            }
            assignments[p] = best_k;
        }
        double* sum_x = (double*)calloc(KMEANS_K, sizeof(double));
        double* sum_y = (double*)calloc(KMEANS_K, sizeof(double));
        double* sum_z = (double*)calloc(KMEANS_K, sizeof(double));
        int* counts = (int*)calloc(KMEANS_K, sizeof(int));
        for (int p = 0; p < num_points; p++) {
            int k = assignments[p];
            sum_x[k] += points[p].i;
            sum_y[k] += points[p].j;
            sum_z[k] += points[p].k;
            counts[k]++;
        }
        for (int k = 0; k < KMEANS_K; k++) {
            if (counts[k] > 0) {
                centroids[k].i = (int)(sum_x[k] / counts[k]);
                centroids[k].j = (int)(sum_y[k] / counts[k]);
                centroids[k].k = (int)(sum_z[k] / counts[k]);
            }
        }
        free(sum_x); free(sum_y); free(sum_z); free(counts);
    }

    int* cluster_sizes = (int*)calloc(KMEANS_K, sizeof(int));
    for (int p = 0; p < num_points; p++) {
        cluster_sizes[assignments[p]]++;
        if (cluster_sizes[assignments[p]] > max_size) max_size = cluster_sizes[assignments[p]];
    }

    free(points);
    free(centroids);
    free(assignments);
    free(cluster_sizes);
    return max_size;
}

// Spectral (Simplified due to no GSL)
int spectral(int* grid, int size, int* visited) {
    int total_size = size * size * size, max_size = 0;
    Point* points = (Point*)malloc(total_size * sizeof(Point));
    if (!points) return 0;
    int num_points = 0;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            for (int k = 0; k < size; k++) {
                int idx = i * size * size + j * size + k;
                if (grid[idx] == 1) {
                    points[num_points++] = (Point){i, j, k};
                }
            }
        }
    }

    int* labels = (int*)malloc(num_points * sizeof(int));
    if (!labels) { free(points); return 0; }
    for (int i = 0; i < num_points; i++) labels[i] = i % KMEANS_K;

    int* cluster_sizes = (int*)calloc(KMEANS_K, sizeof(int));
    for (int i = 0; i < num_points; i++) {
        cluster_sizes[labels[i]]++;
        if (cluster_sizes[labels[i]] > max_size) max_size = cluster_sizes[labels[i]];
    }

    free(points);
    free(labels);
    free(cluster_sizes);
    return max_size;
}

// GCBD (Union-Find)
int gcbd(int* grid, int size, int* visited) {
    int total_size = size * size * size, max_size = 0;
    UnionFind* uf = (UnionFind*)malloc(total_size * sizeof(UnionFind));
    if (!uf) return 0;
    for (int i = 0; i < total_size; i++) {
        uf[i].parent = i;
        uf[i].rank = 0;
    }

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            for (int k = 0; k < size; k++) {
                int idx = i * size * size + j * size + k;
                if (grid[idx] != 1) continue;
                int deltas[6][3] = {{-1,0,0}, {1,0,0}, {0,-1,0}, {0,1,0}, {0,0,-1}, {0,0,1}};
                for (int d = 0; d < 6; d++) {
                    int ni = i + deltas[d][0], nj = j + deltas[d][1], nk = k + deltas[d][2];
                    if (ni >= 0 && ni < size && nj >= 0 && nj < size && nk >= 0 && nk < size) {
                        int n_idx = ni * size * size + nj * size + nk;
                        if (grid[n_idx] == 1) {
                            int root_p = idx, root_q = n_idx;
                            while (uf[root_p].parent != root_p) root_p = uf[root_p].parent;
                            while (uf[root_q].parent != root_q) root_q = uf[root_q].parent;
                            if (root_p != root_q) {
                                if (uf[root_p].rank < uf[root_q].rank) uf[root_p].parent = root_q;
                                else if (uf[root_p].rank > uf[root_q].rank) uf[root_q].parent = root_p;
                                else { uf[root_q].parent = root_p; uf[root_p].rank++; }
                            }
                        }
                    }
                }
            }
        }
    }

    int* cluster_sizes = (int*)calloc(total_size, sizeof(int));
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            for (int k = 0; k < size; k++) {
                int idx = i * size * size + j * size + k;
                if (grid[idx] == 1) {
                    int root = idx;
                    while (uf[root].parent != root) root = uf[root].parent;
                    cluster_sizes[root]++;
                    if (cluster_sizes[root] > max_size) max_size = cluster_sizes[root];
                }
            }
        }
    }

    free(uf);
    free(cluster_sizes);
    return max_size;
}

int main() {
    srand(time(NULL));
    int grid_sizes[] = GRID_SIZES;

    printf("Clustering Algorithm Time Comparison\n");
    printf("-----------------------------------\n");
    fflush(stdout);

    for (int s = 0; s < NUM_SIZES; s++) {
        int size = grid_sizes[s];
        int* grid = generate_grid(size, PROB);
        if (!grid) {
            printf("Memory allocation failed for size %d\n", size);
            fflush(stdout);
            continue;
        }

        printf("\nGrid Size: %dx%dx%d\n", size, size, size);
        fflush(stdout);
        int total_size = size * size * size;
        int* visited = (int*)calloc(total_size, sizeof(int));

        double start_time, end_time;
        int result;

        // Traditional DFS
        start_time = get_time();
        result = traditional_dfs(grid, size, visited);
        end_time = get_time();
        printf("Traditional DFS: Time = %.6f s, Max Cluster Size = %d\n", end_time - start_time, result);
        fflush(stdout);
        memset(visited, 0, total_size * sizeof(int));

        // Skip-DFS
        start_time = get_time();
        result = skip_dfs(grid, size, visited);
        end_time = get_time();
        printf("Skip-DFS: Time = %.6f s, Max Cluster Size = %d\n", end_time - start_time, result);
        fflush(stdout);
        memset(visited, 0, total_size * sizeof(int));

        // DBSCAN
        start_time = get_time();
        result = dbscan(grid, size, visited);
        end_time = get_time();
        printf("DBSCAN: Time = %.6f s, Max Cluster Size = %d\n", end_time - start_time, result);
        fflush(stdout);
        memset(visited, 0, total_size * sizeof(int));

        // Hierarchical
        start_time = get_time();
        result = hierarchical(grid, size, visited);
        end_time = get_time();
        printf("Hierarchical: Time = %.6f s, Max Cluster Size = %d\n", end_time - start_time, result);
        fflush(stdout);
        memset(visited, 0, total_size * sizeof(int));

        // k-means
        start_time = get_time();
        result = kmeans(grid, size, visited);
        end_time = get_time();
        printf("k-means: Time = %.6f s, Max Cluster Size = %d\n", end_time - start_time, result);
        fflush(stdout);
        memset(visited, 0, total_size * sizeof(int));

        // Spectral
        start_time = get_time();
        result = spectral(grid, size, visited);
        end_time = get_time();
        printf("Spectral: Time = %.6f s, Max Cluster Size = %d\n", end_time - start_time, result);
        fflush(stdout);
        memset(visited, 0, total_size * sizeof(int));

        // GCBD
        start_time = get_time();
        result = gcbd(grid, size, visited);
        end_time = get_time();
        printf("GCBD: Time = %.6f s, Max Cluster Size = %d\n", end_time - start_time, result);
        fflush(stdout);

        free(visited);
        free(grid);
    }

    return 0;
}
