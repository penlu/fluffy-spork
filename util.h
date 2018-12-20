struct graph;

int urand(int max); // generate a uniform random integer [0, max)
double absd(double x); // absolute value

// return number of steps walked
int walk(int max_steps, int p_param, int print_freq, struct graph *graph, int **vr, int **cr);
int walk2(int max_steps, int p_param, int print_freq, struct graph *graph, int **vr, int **cr);
