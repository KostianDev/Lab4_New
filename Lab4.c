#include <stdio.h>
#include <stdlib.h>
#include <SDL2/SDL.h>
#include <math.h>
#include <SDL2/SDL_ttf.h>

#define n1 3
#define n2 1
#define n3 0
#define n4 3
#define N (10 + n3)
#define max_rand 199
#define min_rand 0
#define k1 (1.0 - n3 * 0.01 - n4 * 0.01 - 0.3)
#define k2 (1.0 - n3 * 0.005 - n4 * 0.005 - 0.27)
#define WIDTH 800
#define HEIGHT 600
#define seed (1000 * n1 + 100 * n2 + 10 * n3 + n4)
#define arrowAngle  (M_PI / 5.0)
#define sizeMult 3
#define shiftAngle (M_PI / 18.0)

typedef struct passMatrix {
    double matrix[N][N];
} matrix;

typedef struct node {
    int key;
    int x;
    int y;
    int pos;
    int degree;
    int outDegree;
    int inDegree;
    struct node *next_node;
} l_list;

l_list *l_list_init(int key, int x, int y, int pos, int degree, int outDegree, int inDegree);

l_list *addto_list(l_list *l_pointer, int key, int x, int y, int pos, int degree, int outDegree, int inDegree);

l_list *delfrom_start(l_list *l_pointer);

l_list *find_num(l_list *l_pointer, int key);

matrix generateDirectedMatrix(int isNew);

matrix generateUndirectedMatrix(matrix passMatrix);

matrix generateIdentityMatrix();

matrix generateReachabilityMatrix(matrix passMatrix);

matrix matrixMult(matrix matrixA, matrix matrixB);

matrix powMatrix(matrix passMatrix, int degree);

matrix boolTransform(matrix passMatrix);

matrix transposeMatrix(matrix passMatrix);

matrix multByElements(matrix matrixA, matrix matrixB);

matrix matrixUnion(matrix matrixA, matrix matrixB);

matrix makeIrreflexive(matrix passMatrix);

matrix switchLines(matrix passMatrix, int lines[N], int *number1, int *number2);

matrix generateCondensedMatrix(matrix graphMatrix, matrix KMatrix, int *n);

int findMidVertex(matrix matrixA, matrix matrixB, int i, int j, int shift, int isSquare);

void printMatrix(matrix passMatrix, FILE *fptr, const char[]);

void printDegrees(FILE *fptr, l_list *list_ptr, const char graphName[]);

void printHalfDegrees(FILE *fptr, l_list *list_ptr, const char graphName[]);

void printVertexStates(FILE *fptr, l_list *list_ptr, const char graphName[]);

void printDirTrailsL2(matrix squaredMatrix, matrix passMatrix, FILE *fptr);

void printDirTrailsL3(matrix cubedMatrix, matrix squaredMatrixTransformed, matrix passMatrix, FILE *fptr);

matrix printStrConnComponents(int vertices[N], FILE *fptr);

void clearScreen(SDL_Renderer *Renderer);

void drawCircle(SDL_Renderer *renderer, int32_t centreX, int32_t centreY, int32_t radius);

void drawBezierCurve(SDL_Renderer *renderer, int x1, int x2, int x3, int x4, int y1, int y2, int y3, int y4);

void drawArrowHead(SDL_Renderer *renderer, int endX, int endY, int gap, double angle);

void drawVertexNumber(SDL_Renderer *renderer, int number, int x, int y, int gap);

l_list *drawGraph(SDL_Renderer *renderer, SDL_Window *window, matrix matrix, l_list *list_ptr, int size);

void drawConnections(SDL_Renderer *renderer, l_list *node1, l_list *node2,
                     int r, int dir, int width, int height, int gap2, int size);

int main(int argc, char *argv[]) {
    srand(seed);

    FILE *fptr;

    TTF_Init();

    const char dirG[] = "directed graph";
    const char undirG[] = "undirected graph";
    const char newDirG[] = "new directed graph";
    const char sqrNewDirG[] = "squared new directed graph";
    const char cbNewDirG[] = "cubed new directed graph";
    int flag = 0, vertices[N], lineToSwitch1, lineToSwitch2, n = 0;

    matrix directedMatrix = generateDirectedMatrix(0);
    matrix undirectedMatrix = generateUndirectedMatrix(directedMatrix);
    matrix newDirectedMatrix = generateDirectedMatrix(1);
    matrix squaredMatrix = powMatrix(newDirectedMatrix, 2);
    matrix cubedMatrix = powMatrix(newDirectedMatrix, 3);
    matrix reachabilityMatrix = generateReachabilityMatrix(newDirectedMatrix);
    matrix stronglyConnectedMatrix = multByElements(reachabilityMatrix,
                                                    transposeMatrix(reachabilityMatrix));
    matrix switchedMatrix = switchLines(stronglyConnectedMatrix, vertices,
                                        &lineToSwitch1, &lineToSwitch2);


    printMatrix(directedMatrix, fptr, dirG);
    printMatrix(undirectedMatrix, fptr, undirG);
    printMatrix(newDirectedMatrix, fptr, newDirG);
    printMatrix(squaredMatrix, fptr, sqrNewDirG);
    printMatrix(cubedMatrix, fptr, cbNewDirG);
    printMatrix(reachabilityMatrix, fptr, "reachability");
    printMatrix(stronglyConnectedMatrix, fptr, "strongly connected");
    printMatrix(switchedMatrix, fptr, "switched strongly connected");
    matrix KMatrix = printStrConnComponents(vertices, fptr);
    matrix condensedMatrix = generateCondensedMatrix(newDirectedMatrix, KMatrix, &n);
    printMatrix(condensedMatrix, fptr, "condensed graph");

    l_list *list1_ptr, *list2_ptr, *list3_ptr, *list4_ptr;

    SDL_Init(SDL_INIT_EVERYTHING);

    SDL_Window *directedWindow = SDL_CreateWindow("Directed Graph", SDL_WINDOWPOS_UNDEFINED,
                                                  SDL_WINDOWPOS_UNDEFINED, WIDTH,
                                                  HEIGHT, SDL_WINDOW_SHOWN);
    SDL_Renderer *directedRenderer = SDL_CreateRenderer(directedWindow, -1,
                                                        SDL_RENDERER_ACCELERATED);

    SDL_Window *undirectedWindow = SDL_CreateWindow("Undirected Graph", SDL_WINDOWPOS_UNDEFINED,
                                                    SDL_WINDOWPOS_UNDEFINED, WIDTH,
                                                    HEIGHT, SDL_WINDOW_SHOWN);
    SDL_Renderer *undirectedRenderer = SDL_CreateRenderer(undirectedWindow, -1,
                                                          SDL_RENDERER_ACCELERATED);

    SDL_Window *newDirectedWindow = SDL_CreateWindow("New Directed Graph", SDL_WINDOWPOS_UNDEFINED,
                                                     SDL_WINDOWPOS_UNDEFINED, WIDTH,
                                                     HEIGHT, SDL_WINDOW_SHOWN);
    SDL_Renderer *newDirectedRenderer = SDL_CreateRenderer(newDirectedWindow, -1,
                                                           SDL_RENDERER_ACCELERATED);

    SDL_Window *condensedGraphWindow = SDL_CreateWindow("Graph Condensation", SDL_WINDOWPOS_UNDEFINED,
                                                        SDL_WINDOWPOS_UNDEFINED, WIDTH,
                                                        HEIGHT, SDL_WINDOW_SHOWN);
    SDL_Renderer *condensedGraphRenderer = SDL_CreateRenderer(condensedGraphWindow, -1,
                                                              SDL_RENDERER_ACCELERATED);

    SDL_SetWindowResizable(directedWindow, SDL_TRUE);
    SDL_SetWindowResizable(undirectedWindow, SDL_TRUE);
    SDL_SetWindowResizable(newDirectedWindow, SDL_TRUE);
    SDL_SetWindowResizable(condensedGraphWindow, SDL_TRUE);

    SDL_Event event;
    int quit = 0;

    while (!quit) {
        while (SDL_PollEvent(&event)) {
            if ((event.type == SDL_WINDOWEVENT) && (event.window.event == SDL_WINDOWEVENT_CLOSE)) {
                quit = 1;
            }
        }

        clearScreen(directedRenderer);
        clearScreen(undirectedRenderer);
        clearScreen(newDirectedRenderer);
        clearScreen(condensedGraphRenderer);

        SDL_SetRenderDrawColor(undirectedRenderer, 255, 255, 255, 0);
        SDL_SetRenderDrawColor(directedRenderer, 255, 255, 255, 0);
        SDL_SetRenderDrawColor(newDirectedRenderer, 255, 255, 255, 0);
        SDL_SetRenderDrawColor(condensedGraphRenderer, 255, 255, 255, 0);

        list1_ptr = drawGraph(directedRenderer, directedWindow,
                              directedMatrix, list1_ptr, N);
        list2_ptr = drawGraph(undirectedRenderer, undirectedWindow,
                              undirectedMatrix, list2_ptr, N);
        list3_ptr = drawGraph(newDirectedRenderer, newDirectedWindow,
                              newDirectedMatrix, list3_ptr, N);
        list4_ptr = drawGraph(condensedGraphRenderer, condensedGraphWindow,
                              stronglyConnectedMatrix, list4_ptr, n);

        SDL_RenderPresent(undirectedRenderer);
        SDL_RenderPresent(directedRenderer);
        SDL_RenderPresent(newDirectedRenderer);
        SDL_RenderPresent(condensedGraphRenderer);

        if (!flag) {
            printDegrees(fptr, list1_ptr, dirG);
            printDegrees(fptr, list2_ptr, undirG);
            printHalfDegrees(fptr, list1_ptr, dirG);
            printVertexStates(fptr, list1_ptr, dirG);
            printVertexStates(fptr, list2_ptr, undirG);
            printHalfDegrees(fptr, list3_ptr, newDirG);
            printDirTrailsL2(squaredMatrix, newDirectedMatrix, fptr);
            printDirTrailsL3(cubedMatrix, boolTransform(squaredMatrix),
                             newDirectedMatrix, fptr);
            flag = 1;
        }
    }

    SDL_DestroyWindow(directedWindow);
    SDL_DestroyRenderer(directedRenderer);

    SDL_DestroyWindow(undirectedWindow);
    SDL_DestroyRenderer(undirectedRenderer);

    SDL_DestroyWindow(newDirectedWindow);
    SDL_DestroyRenderer(newDirectedRenderer);

    SDL_DestroyWindow(condensedGraphWindow);
    SDL_DestroyRenderer(condensedGraphRenderer);

    SDL_Quit();

    while (list1_ptr != NULL)
        list1_ptr = delfrom_start(list1_ptr);
    while (list2_ptr != NULL)
        list2_ptr = delfrom_start(list2_ptr);
    while (list3_ptr != NULL)
        list3_ptr = delfrom_start(list3_ptr);
    while (list4_ptr != NULL)
        list4_ptr = delfrom_start(list4_ptr);

    unlink("Output.txt");
    return 0;
}

l_list *l_list_init(int key, int x, int y, int pos, int degree, int outDegree, int inDegree) {
    l_list *l_pointer;
    l_pointer = malloc(sizeof(struct node));
    *l_pointer = (l_list) {
            .key = key,
            .x = x,
            .y = y,
            .pos = pos,
            .degree = degree,
            .outDegree = outDegree,
            .inDegree = inDegree,
            .next_node = NULL
    };
    return l_pointer;
}

l_list *addto_list(l_list *l_pointer, int key, int x, int y, int pos, int degree, int outDegree, int inDegree) {
    l_list *new_node;
    new_node = malloc(sizeof(struct node));
    new_node->key = key;
    new_node->x = x;
    new_node->y = y;
    new_node->pos = pos;
    new_node->degree = degree;
    new_node->next_node = l_pointer;
    new_node->outDegree = outDegree;
    new_node->inDegree = inDegree;
    return new_node;
}

l_list *delfrom_start(l_list *l_pointer) {
    l_list *node_ptr;
    node_ptr = l_pointer->next_node;
    free(l_pointer);
    return node_ptr;
}

l_list *find_num(l_list *l_pointer, int key) {
    l_list *this_node = l_pointer;

    while (this_node != NULL) {
        if (this_node->key == key) return this_node;
        else this_node = this_node->next_node;
    }
    return NULL;
}

matrix generateDirectedMatrix(int isNew) {
    matrix passMatrix;
    double a;
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) {
            a = (rand() % (max_rand + 1 - min_rand) + min_rand) / 100.0;
            a *= (isNew) ? k2 : k1;
            passMatrix.matrix[i][j] = a < 1.0 ? 0 : 1;
        }
    return passMatrix;
}

matrix generateUndirectedMatrix(matrix passMatrix) {
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            if (passMatrix.matrix[i][j] == 1) passMatrix.matrix[j][i] = 1;
    return passMatrix;
}

matrix generateIdentityMatrix() {
    matrix I;

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) {
            if (i == j) I.matrix[i][j] = 1;
            else I.matrix[i][j] = 0;
        }
    return I;
}

matrix generateReachabilityMatrix(matrix passMatrix) {
    matrix expMatrix;
    matrix midMatrix = generateIdentityMatrix();

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            for (int k = 1; k < N; ++k) {
                expMatrix = powMatrix(passMatrix, k);
                midMatrix.matrix[i][j] += expMatrix.matrix[i][j];
            }
    return boolTransform(midMatrix);
}

matrix matrixMult(matrix matrixA, matrix matrixB) {
    matrix midMatrix = {{0}};

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            for (int l = 0; l < N; ++l)
                midMatrix.matrix[i][j] += (int) (matrixA.matrix[i][l] * matrixB.matrix[l][j]);
    return midMatrix;
}

matrix powMatrix(matrix passMatrix, int degree) {
    if (degree > 1) {
        if (degree % 2 == 0)
            passMatrix = powMatrix(matrixMult(passMatrix, passMatrix), degree / 2);
        else
            passMatrix = matrixMult(passMatrix,
                                    powMatrix(matrixMult(passMatrix, passMatrix),
                                              (degree - 1) / 2));
    }
    return passMatrix;
}

matrix boolTransform(matrix passMatrix) {
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) {
            if (passMatrix.matrix[i][j] >= 1) passMatrix.matrix[i][j] = 1;
            else passMatrix.matrix[i][j] = 0;
        }
    return passMatrix;
}

matrix transposeMatrix(matrix passMatrix) {
    matrix midMatrix;

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) {
            midMatrix.matrix[i][j] = passMatrix.matrix[j][i];
        }
    return midMatrix;
}

matrix multByElements(matrix matrixA, matrix matrixB) {
    matrix midMatrix;

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) {
            midMatrix.matrix[i][j] = matrixA.matrix[i][j] * matrixB.matrix[i][j];
        }
    return midMatrix;
}

matrix matrixUnion(matrix matrixA, matrix matrixB) {
    matrix midMatrix;

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) {
            if (matrixA.matrix[i][j] >= 1 || matrixB.matrix[i][j] >= 1)
                midMatrix.matrix[i][j] = (matrixA.matrix[i][j] >= matrixB.matrix[i][j]) ?
                                         matrixA.matrix[i][j] : matrixB.matrix[i][j];
        }
    return midMatrix;
}

matrix makeIrreflexive(matrix passMatrix) {
    for (int i = 0; i < N; ++i)
        passMatrix.matrix[i][i] = 0;
    return passMatrix;
}

matrix switchLines(matrix passMatrix, int lines[N], int *number1, int *number2) {
    matrix midMatrix = passMatrix;
    int a, line = 0, n = 0;

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            line += (int) midMatrix.matrix[j][i] * (j + 1);
        }
        lines[i] = line;
        line = 0;
    }

    for (int i = 1; i < N - 1; ++i)
        if (lines[i - 1] == lines[i + 1] && lines[i - 1] != lines[i]) {
            *number1 = i;
            line = lines[i - 1];
            break;
        } else if (lines[i] == lines[i + 1] && lines[i - 1] != lines[i]) {
            *number1 = i - 1;
            line = lines[i];
            break;
        } else if (lines[i - 1] == lines[i] && lines[i + 1] != lines[i]) {
            *number1 = i + 1;
            line = lines[i];
            break;
        }
    for (int i = N - 1; i >= 0; --i)
        if (lines[i] == line) {
            *number2 = i;
            break;
        }

    for (int i = 0; i < N; ++i) {
        a = (int) midMatrix.matrix[i][*number2];
        midMatrix.matrix[i][*number2] = midMatrix.matrix[i][*number1];
        midMatrix.matrix[i][*number1] = a;
    }
    for (int j = 0; j < N; ++j) {
        a = (int) midMatrix.matrix[*number2][j];
        midMatrix.matrix[*number2][j] = midMatrix.matrix[*number1][j];
        midMatrix.matrix[*number1][j] = a;
    }
    return midMatrix;
}

matrix printStrConnComponents(int vertices[N], FILE *fptr) {
    int n = 1, vertex;
    matrix KMatrix;

    fptr = fopen("Output.txt", "a");
    fprintf(fptr, "\nHere are your strongly connected components:\n\n");

    while (vertices[N - 1] != -1) {
        fprintf(fptr, "K%d = {", n);
        for (int i = 0; i < N; ++i) {
            if (vertices[i] != -1) {
                vertex = vertices[i];
                for (int j = 0; j < N; ++j) {
                    if (vertices[j] == vertex && vertices[j] != 0) {
                        fprintf(fptr, " %d ", j + 1);
                        KMatrix.matrix[n - 1][j] = 1;
                        vertices[j] = -1;
                    } else
                        KMatrix.matrix[n - 1][j] = 0;
                }
                break;
            }
        }
        fprintf(fptr, "}\n");
        n++;
    }
    fprintf(fptr, "\n");

    fclose(fptr);
    return KMatrix;
}

matrix generateCondensedMatrix(matrix graphMatrix, matrix KMatrix, int *n) {

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            if (KMatrix.matrix[i][j] != 0) {
                *n = *n + 1;
                i++;
                break;
            }
    *n = *n + 1;

    matrix condensedMatrix;

    for (int k = 0; k < *n; ++k)
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j) {
                if (KMatrix.matrix[k][i] == 1 && KMatrix.matrix[k + 1][j] == 1)
                    if (graphMatrix.matrix[i][j] == 1) {
                        condensedMatrix.matrix[k][k + 1] = 1;
                        break;
                    } else
                        condensedMatrix.matrix[i][j] = 0;
            }

    return condensedMatrix;
}

int findMidVertex(matrix matrixA, matrix matrixB, int i, int j, int shift, int isSquare) {
    int vertexKey = 0;

    for (int k = 0 + shift; k < N; ++k) {
        if (matrixA.matrix[i][k] == 1 && matrixB.matrix[k][j] == 1) {
            vertexKey = (isSquare && (k == i && k == j)) ? 0 : k + 1;
            break;
        }
    }
    return vertexKey;
}

void printMatrix(matrix passMatrix, FILE *fptr, const char graphName[]) {
    fptr = fopen("Output.txt", "a");

    fprintf(fptr, "\nHere is your %s matrix:\n\n", graphName);

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            fprintf(fptr, " %2.0lf", passMatrix.matrix[i][j]);
        }
        fprintf(fptr, "\n");
    }
    fprintf(fptr, "\n");

    fclose(fptr);
}

void printDegrees(FILE *fptr, l_list *list_ptr, const char graphName[]) {
    l_list *prev_node = list_ptr;
    l_list *this_node = prev_node->next_node;
    int isRegular = 1;

    fptr = fopen("Output.txt", "a");

    fprintf(fptr, "\nHere are your %s vertices' degrees:\n\n", graphName);
    fprintf(fptr, "vertex number:%d vertex degree:%d;\n", prev_node->key, prev_node->degree);

    while (this_node != NULL) {
        fprintf(fptr, "vertex number:%d vertex degree:%d;\n", this_node->key, this_node->degree);
        if (prev_node->degree != this_node->degree) isRegular = 0;
        prev_node = this_node;
        this_node = this_node->next_node;
    }

    if (isRegular) fprintf(fptr, "\nThis graph is regular with degree %d\n\n", prev_node->degree);
    else fprintf(fptr, "\nThis graph is not regular\n\n");

    fclose(fptr);
}

void printHalfDegrees(FILE *fptr, l_list *list_ptr, const char graphName[]) {
    l_list *this_node = list_ptr;

    fptr = fopen("Output.txt", "a");

    fprintf(fptr, "\nHere are your %s vertices' half-degrees:\n\n", graphName);

    while (this_node != NULL) {
        fprintf(fptr, "vertex number:%d vertex outdegree:%d vertex indegree:%d;\n",
                this_node->key, this_node->outDegree, this_node->inDegree);
        this_node = this_node->next_node;
    }
    fprintf(fptr, "\n");

    fclose(fptr);
}

void printVertexStates(FILE *fptr, l_list *list_ptr, const char graphName[]) {
    l_list *this_node = list_ptr;

    fptr = fopen("Output.txt", "a");

    fprintf(fptr, "\nHere is a list of %s's isolated and leaf vertices:\n\n", graphName);

    while (this_node != NULL) {
        if (this_node->degree == 0)
            fprintf(fptr, "vertex number:%d, vertex is isolated;\n", this_node->key);
        if (this_node->degree == 1)
            fprintf(fptr, "vertex number:%d, this is a leaf vertex;\n", this_node->key);
        this_node = this_node->next_node;
    }
    fprintf(fptr, "\n");

    fclose(fptr);
}

void printDirTrailsL2(matrix squaredMatrix, matrix passMatrix, FILE *fptr) {
    int vertexKey, shift = 0;

    fptr = fopen("Output.txt", "a");
    fprintf(fptr, "\nHere are your graph's directed trails with length 2:\n\n");

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) {
            if (squaredMatrix.matrix[i][j] >= 1) {
                vertexKey = N;
                while (vertexKey != 0) {
                    vertexKey = findMidVertex(passMatrix, passMatrix, i, j, shift, 1);
                    if (vertexKey != 0) {
                        fprintf(fptr, "Trail: %d-%d-%d\n", i + 1, vertexKey, j + 1);
                        shift = vertexKey;
                    }
                }
                shift = 0;
            }
        }
    fprintf(fptr, "\n");

    fclose(fptr);
}

void printDirTrailsL3(matrix cubedMatrix, matrix squaredMatrixTransformed, matrix passMatrix, FILE *fptr) {
    int vertexKey1, vertexKey2, shift1 = 0, shift2 = 0;

    fptr = fopen("Output.txt", "a");
    fprintf(fptr, "\nHere are your graph's directed trails with length 3:\n\n");

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) {
            if (cubedMatrix.matrix[i][j] >= 1) {
                vertexKey2 = N;
                while (vertexKey2 != 0) {
                    vertexKey2 = findMidVertex(squaredMatrixTransformed,
                                               passMatrix, i, j, shift2, 0);
                    if (vertexKey2 != 0) {
                        vertexKey1 = N;
                        while (vertexKey1 != 0) {
                            vertexKey1 = findMidVertex(passMatrix, passMatrix,
                                                       i, vertexKey2 - 1, shift1, 1);
                            if (vertexKey1 != 0) {
                                if (vertexKey1 == vertexKey2 && vertexKey2 == j + 1)
                                    shift1 = vertexKey1;
                                else {
                                    fprintf(fptr, "Trail: %d-%d-%d-%d\n", i + 1, vertexKey1,
                                            vertexKey2, j + 1);
                                    shift1 = vertexKey1;
                                }
                            }
                        }
                        shift1 = 0;
                        shift2 = vertexKey2;
                    }
                }
                shift2 = 0;
            }
        }
    fprintf(fptr, "\n");

    fclose(fptr);
}

void clearScreen(SDL_Renderer *Renderer) {
    SDL_SetRenderDrawColor(Renderer, 0, 0, 0, 255);
    SDL_RenderClear(Renderer);
}

void drawCircle(SDL_Renderer *renderer, int32_t centreX, int32_t centreY, int32_t radius) {
    const int32_t diameter = (radius * 2);

    int32_t x = (radius - 1);
    int32_t y = 0;
    int32_t tx = 1;
    int32_t ty = 1;
    int32_t error = (tx - diameter);

    while (x >= y) {
        SDL_RenderDrawPoint(renderer, centreX + x, centreY - y);
        SDL_RenderDrawPoint(renderer, centreX + x, centreY + y);
        SDL_RenderDrawPoint(renderer, centreX - x, centreY - y);
        SDL_RenderDrawPoint(renderer, centreX - x, centreY + y);
        SDL_RenderDrawPoint(renderer, centreX + y, centreY - x);
        SDL_RenderDrawPoint(renderer, centreX + y, centreY + x);
        SDL_RenderDrawPoint(renderer, centreX - y, centreY - x);
        SDL_RenderDrawPoint(renderer, centreX - y, centreY + x);

        if (error <= 0) {
            ++y;
            error += ty;
            ty += 2;
        }

        if (error > 0) {
            --x;
            tx += 2;
            error += (tx - diameter);
        }
    }
}

void drawBezierCurve(SDL_Renderer *renderer, int x1, int x2, int x3, int x4, int y1, int y2, int y3, int y4) {
    double xu, yu, u;
    for (u = 0.0; u <= 1.0; u += 0.0005) {
        xu = pow(1 - u, 3) * x1 + 3 * u * pow(1 - u, 2) * x2 + 3 * pow(u, 2) * (1 - u) * x3
             + pow(u, 3) * x4;
        yu = pow(1 - u, 3) * y1 + 3 * u * pow(1 - u, 2) * y2 + 3 * pow(u, 2) * (1 - u) * y3
             + pow(u, 3) * y4;
        SDL_RenderDrawPoint(renderer, (int) xu, (int) yu);
    }
}

void drawArrowHead(SDL_Renderer *renderer, int endX, int endY, int gap, double angle) {
    double arrowSize = (double) gap / 3 / sizeMult;

    double x1 = endX - arrowSize * cos(angle + arrowAngle);
    double y1 = endY - arrowSize * sin(angle + arrowAngle);
    double x2 = endX - arrowSize * cos(angle - arrowAngle);
    double y2 = endY - arrowSize * sin(angle - arrowAngle);

    SDL_Vertex vertex_1 = {{(float) endX, (float) endY},
                           {255,          255, 255, 255},
                           {1,            1}};
    SDL_Vertex vertex_2 = {{(float) x1, (float) y1},
                           {255,        255, 255, 255},
                           {1,          1}};
    SDL_Vertex vertex_3 = {{(float) x2, (float) y2},
                           {255,        255, 255, 255},
                           {1,          1}};

    SDL_Vertex vertices[] = {
            vertex_1,
            vertex_2,
            vertex_3
    };

    SDL_RenderGeometry(renderer, NULL, vertices, 3, NULL, 0);
}

void drawVertexNumber(SDL_Renderer *renderer, int number, int x, int y, int gap) {
    SDL_Color color = {255, 255, 255, 255};
    TTF_Font *font = TTF_OpenFont("arial.ttf", gap / 3);
    char numberString[4];

    snprintf(numberString, sizeof(numberString), "%d", number);

    SDL_Surface *surface = TTF_RenderText_Solid(font, numberString, color);
    SDL_Texture *texture = SDL_CreateTextureFromSurface(renderer, surface);

    int texW, texH;
    SDL_QueryTexture(texture, NULL, NULL, &texW, &texH);

    SDL_Rect rect = {x - texW / 2, y - texH / 2, texW, texH};
    SDL_RenderCopy(renderer, texture, NULL, &rect);

    SDL_FreeSurface(surface);
    SDL_DestroyTexture(texture);
    TTF_CloseFont(font);
}

l_list *drawGraph(SDL_Renderer *renderer, SDL_Window *window, matrix matrix, l_list *list_ptr, int size) {
    int width, height, xTopCircles, xLowCircles, yCircles, turn = 1, mid;
    int gap, gap1, gap2, gap3, r, key = 1, pos = 1, degree = 0, outDegree = 0, inDegree = 0, flag = 1;

    mid = (abs(size - 5)) / 4;
    yCircles = 2 + mid * 2;
    xLowCircles = (size - yCircles) / 2;
    xTopCircles = size - xLowCircles - yCircles;
    SDL_GetWindowSize(window, &width, &height);

    r = (height > width || (size >= 7 && size < 9)) ? width / (int) ((sizeMult + 1) * (xTopCircles + 1) + sizeMult) :
        height / ((sizeMult + 1) * (yCircles / 2 + 1) + sizeMult);

    gap = r * sizeMult;
    gap1 = (width - 2 * gap - (xTopCircles + 1) * 2 * r) / xTopCircles;
    gap2 = (height - (yCircles / 2 + 1) * 2 * r) / (yCircles / 2 + 2);
    gap3 = (width - 2 * gap - (xLowCircles + 1) * 2 * r) / xLowCircles;

    drawCircle(renderer, gap + r, gap2 + r, r);
    drawVertexNumber(renderer, key, gap + r, gap2 + r, gap);
    list_ptr = l_list_init(key, gap + r, gap2 + r, pos, degree, outDegree, inDegree);
    key++;

    while (key <= size) {
        switch (turn % 4) {
            case 1:
                for (int j = 1; j < xTopCircles; ++j) {
                    drawCircle(renderer, gap + r + j * (2 * r + gap1),
                               gap2 + r, r);
                    drawVertexNumber(renderer, key, gap + r + j * (2 * r + gap1),
                                     gap2 + r, gap);
                    list_ptr = addto_list(list_ptr, key,
                                          gap + r + j * (2 * r + gap1), gap2 + r, pos,
                                          degree, outDegree, inDegree);
                    key++;
                }
                break;
            case 2:
                for (int j = 0; j < yCircles / 2; ++j) {
                    drawCircle(renderer, width - gap - r,
                               (gap2 + r) + j * (gap2 + 2 * r), r);
                    drawVertexNumber(renderer, key, width - gap - r,
                                     (gap2 + r) + j * (gap2 + 2 * r), gap);
                    list_ptr = addto_list(list_ptr, key,
                                          width - gap - r, (gap2 + r) + j * (gap2 + 2 * r), pos,
                                          degree, outDegree, inDegree);
                    key++;
                }
                break;
            case 3:
                for (int j = 0; j < xLowCircles; ++j) {
                    drawCircle(renderer, width - gap - r - j * (2 * r + gap3),
                               height - gap2 - r, r);
                    drawVertexNumber(renderer, key, width - gap - r - j * (2 * r + gap3),
                                     height - gap2 - r, gap);
                    list_ptr = addto_list(list_ptr, key,
                                          width - gap - r - j * (2 * r + gap3), height - gap2 - r, pos,
                                          degree, outDegree, inDegree);
                    key++;
                }
                break;
            case 0:
                for (int j = 0; j < yCircles / 2; ++j) {
                    drawCircle(renderer, gap + r,
                               height - (r + gap2) - j * (2 * r + gap2), r);
                    drawVertexNumber(renderer, key, gap + r,
                                     height - (r + gap2) - j * (2 * r + gap2), gap);
                    list_ptr = addto_list(list_ptr, key,
                                          gap + r, height - (r + gap2) - j * (2 * r + gap2), pos,
                                          degree, outDegree, inDegree);
                    key++;
                }
                break;
        }
        turn++;
        pos++;
    }
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
            if (matrix.matrix[i][j] == 1) {
                l_list *node1 = find_num(list_ptr, i + 1);
                l_list *node2 = find_num(list_ptr, j + 1);
                if (i == j) {
                    if (node1->pos >= 3) flag = -1;
                    drawCircle(renderer,
                               (int) (node1->x - flag * (r + (double) r / 2 - flag)
                                                 * cos(M_PI / 4)),
                               (int) (node1->y - flag * (r + (double) r / 2 - flag)
                                                 * sin(M_PI / 4)),
                               r / 2);
                    drawArrowHead(renderer, (int) (node1->x - flag * r * cos(M_PI / 4)),
                                  (int) (node1->y - flag * r * sin(M_PI / 4)), gap,
                                  3 * M_PI / 4 - flag * M_PI / 15);
                    node1->degree += 2;
                    node1->outDegree++;
                    node1->inDegree++;
                } else {
                    int dir = 1;
                    drawConnections(renderer, node1, node2, r, dir, width, height, gap2, size);
                    if (matrix.matrix[j][i] == 1) {
                        matrix.matrix[j][i] = 0;
                        dir = -1;
                        drawConnections(renderer, node2, node1, r, dir, width, height, gap2, size);
                    }
                }
            }
    return list_ptr;
}

void drawConnections(SDL_Renderer *renderer, l_list *node1, l_list *node2,
                     int r, int dir, int width, int height, int gap2, int size) {
    int startX = node1->x, startY = node1->y, endX = node2->x, endY = node2->y;
    int mid1X, mid1Y, mid2X, mid2Y;
    int midX = (node1->x + node2->x) / 2;
    int midY = (node1->y + node2->y) / 2;
    int dirX = (node1->x - width / 2) <= 0 ? -1 : 1;
    int dirY = (node1->y - height / 2) <= 0 ? -1 : 1;
    int gap = sizeMult * r;
    double angle1, angle2;
    if ((abs(node1->key - node2->key) == 1) ||
        ((node1->key == 1 || node2->key == 1) && (node1->key == size || node2->key == size))) {
        angle1 = atan2(endY - startY, endX - startX);
        SDL_RenderDrawLine(renderer, startX + (int) ((double) r * cos(angle1 + shiftAngle)),
                           startY + (int) ((double) r * sin(angle1 + shiftAngle)),
                           endX - (int) ((double) r * cos(angle1 - shiftAngle)),
                           endY - (int) ((double) r * sin(angle1 - shiftAngle)));
        drawArrowHead(renderer, endX - (int) ((double) r * cos(angle1 - shiftAngle)),
                      endY - (int) ((double) r * sin(angle1 - shiftAngle)), gap, angle1);
    } else {
        if (startX == endX && (startX == gap + r || startX == width - gap - r)) {
            mid1X = (int) ((double) (startX + midX) / 2 + gap * dir * dirX);
            mid1Y = (startY + midY) / 2;
            mid2X = (int) ((double) (midX + endX) / 2 + gap * dir * dirX);
            mid2Y = (midY + endY) / 2;
            angle1 = atan2(mid1Y - startY, mid1X - startX);
            angle2 = atan2(endY - mid2Y, endX - mid2X);
            drawBezierCurve(renderer, startX + (int) ((double) r * cos(angle1)), mid1X,
                            mid2X, endX - (int) ((double) r * cos(angle2)),
                            startY + (int) ((double) r * sin(angle1)), mid1Y,
                            mid2Y, endY - (int) ((double) r * sin(angle2)));
            drawArrowHead(renderer, endX - (int) ((double) r * cos(angle2)),
                          endY - (int) ((double) r * sin(angle2)), gap, angle2);
        } else if (startY == endY && (startY == gap2 + r || startY == height - gap2 - r)) {
            mid1X = (startX + midX) / 2;
            mid1Y = (int) ((double) (startY + midY) / 2 + gap * dir * dirY);
            mid2X = (midX + endX) / 2;
            mid2Y = (int) ((double) (midY + endY) / 2 + gap * dir * dirY);
            angle1 = atan2(mid1Y - startY, mid1X - startX);
            angle2 = atan2(endY - mid2Y, endX - mid2X);
            drawBezierCurve(renderer, startX + (int) ((double) r * cos(angle1)), mid1X,
                            mid2X, endX - (int) ((double) r * cos(angle2)),
                            startY + (int) ((double) r * sin(angle1)), mid1Y,
                            mid2Y, endY - (int) ((double) r * sin(angle2)));
            drawArrowHead(renderer, endX - (int) ((double) r * cos(angle2)),
                          endY - (int) ((double) r * sin(angle2)), gap, angle2);
        } else {
            angle1 = atan2(endY - startY, endX - startX);
            SDL_RenderDrawLine(renderer, startX + (int) ((double) r * cos(angle1 + shiftAngle)),
                               startY + (int) ((double) r * sin(angle1 + shiftAngle)),
                               endX - (int) ((double) r * cos(angle1 - shiftAngle)),
                               endY - (int) ((double) r * sin(angle1 - shiftAngle)));
            drawArrowHead(renderer, endX - (int) ((double) r * cos(angle1 - shiftAngle)),
                          endY - (int) ((double) r * sin(angle1 - shiftAngle)), gap, angle1);
        }
    }
    node1->degree++;
    node2->degree++;
    node1->outDegree++;
    node2->inDegree++;
}