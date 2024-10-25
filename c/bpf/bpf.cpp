#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <array>
#include <string>
#include <cmath>
#include <armadillo>

using namespace arma;
#define TOL 1E-8

/*==============================================================================
                            STRUCTURE DEFINITIONS
==============================================================================*/
class Meas{
public:
    int dim; 
    vec mean;
    mat cov;
    double T;

    Meas(int DIM, const std::string& M_DIR, const std::string& M_FILE){
        std::string M_PATH = M_DIR + "/" + M_FILE;
        std::ifstream m_file(M_PATH);
        if (!m_file.is_open()) {
            std::cerr << "Error: could not open file " << M_PATH << std::endl;
            exit(EXIT_FAILURE);
        }
        dim = DIM; 
        mean.zeros(DIM);
        cov.zeros(DIM, DIM); 

        std::string line;
        std::getline(m_file, line); // skip label line
        std::getline(m_file, line); // mean vector
        std::istringstream meanStream(line);
        for (int i = 0; i < dim; ++i) {
            meanStream >> mean(i);
        }
        std::getline(m_file, line); // skip blank line
        std::getline(m_file, line); // skip label line
        for (int i = 0; i < dim; ++i) { // read covariance matrix
            std::getline(m_file, line);
            std::istringstream covStream(line);
            for (int j = 0; j < dim; ++j) {
                covStream >> cov(i, j);
            }
        }
        std::getline(m_file, line); // skip blank line
        std::getline(m_file, line); // skip label line
        std::getline(m_file, line); // read T value
        T = std::stod(line);
        m_file.close();
    }
};

class Traj {
public:
    vec coef;

    Traj(vec COEF){
        coef = COEF;
    }
};

class BPF{
public:
    vec P; 
    mat x; 
    int n; 
    int dim; 

    BPF(int N, Meas M){
        x = arma::mvnrnd(M.mean, M.cov, N);
        P.zeros(N);
        double sum = 0; 
        for(int i = 0; i < N; i++){
            P(i) = arma::as_scalar(exp(-0.5*trans(x.col(i)-M.mean)*arma::inv(M.cov)*(x.col(i)-M.mean)));
        }
        P = P/arma::sum(P);
        
        n = N; 
        dim = M.dim; 
    }
};
/*==============================================================================
                            RUNGE-KUTTA 8(7) DEFINITIONS
==============================================================================*/
const vec c1 = {1.0/18, 1.0/12, 1.0/8, 5.0/16, 3.0/8, 59.0/400, 93.0/200, 5490023248.0/9719169821, 13.0/20, 1201146811.0/1299019798, 1.0, 1.0};
mat aT = {{1.0/18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
          {1.0/48, 1.0/16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
          {1.0/32, 0, 3.0/32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
          {5.0/16, 0, -75.0/64, 75.0/64, 0, 0, 0, 0, 0, 0, 0, 0, 0},
          {3.0/80, 0, 0, 3.0/16, 3.0/20, 0, 0, 0, 0, 0, 0, 0, 0},
          {29443841.0/614563906, 0, 0, 77736538.0/692538347, -28693883.0/1125000000, 23124283.0/1800000000, 0, 0, 0, 0, 0, 0, 0},
          {16016141.0/946692911, 0, 0, 61564180.0/158732637, 22789713.0/633445777, 545815736.0/2771057229, -180193667.0/1043307555, 0, 0, 0, 0, 0, 0},
          {39632708.0/573591083, 0, 0, -433636366.0/683701615, -421739975.0/2616292301, 100302831.0/723423059, 790204164.0/839813087, 800635310.0/3783071287, 0, 0, 0, 0, 0},
          {246121993.0/1340847787, 0, 0, -37695042795.0/15268766246, -309121744.0/1061227803, -12992083.0/490766935, 6005943493.0/2108947869, 393006217.0/1396673457, 123872331.0/1001029789, 0, 0, 0, 0},
          {-1028468189.0/846180014, 0, 0, 8478235783.0/508512852, 1311729495.0/1432422823, -10304129995.0/1701304382, -48777925059.0/3047939560, 15336726248.0/1032824649, -45442868181.0/3398467696, 3065993473.0/597172653, 0, 0, 0},
          {185892177.0/718116043, 0, 0, -3185094517.0/667107341, -477755414.0/1098053517, -703635378.0/230739211, 5731566787.0/1027545527, 5232866602.0/850066563, -4093664535.0/808688257, 3962137247.0/1805957418, 65686358.0/487910083, 0, 0},
          {403863854.0/491063109, 0, 0, -5068492393.0/434740067, -411421997.0/543043805, 652783627.0/914296604, 11173962825.0/925320556, -13158990841.0/6184727034, 3936647629.0/1978049680, -160528059.0/685178525, 248638103.0/1413531060, 0, 0}};
const mat a = trans(aT);
const vec b7 = {14005451.0/335480064, 0.0, 0.0, 0.0, 0.0, -59238493.0/1068277825, 181606767.0/758867731, 561292985.0/797845732, -1041891430.0/1371343529, 760417239.0/1151165299, 118820643.0/751138087, -528747749.0/2220607170, 1.0/4};
const vec b8 = {13451932.0/455176623, 0.0, 0.0, 0.0, 0.0, -808719846.0/976000145, 1757004468.0/5645159321, 656045339.0/265891186, -3867574721.0/1518517206, 465885868.0/322736535, 53011238.0/667516719, 2.0/45, 0.0};
const double EPS = 2.220446049250313e-16;
const double POW = 1.0/8;

/*==============================================================================
                            BPF FUNCTION DEFINITIONS
==============================================================================*/
void record_data(const std::string& filename, const BPF p, double t){
    std::ofstream myfile; myfile.open(filename);
    myfile << std::fixed << std::setprecision(10);
    myfile << t << std::endl;

    for(int i=0; i < p.n; i++){
        myfile << p.P(i) << " ";
        for(int j=0; j < p.dim; j++){
            myfile << p.x(j, i) << " ";
        }
        myfile << std::endl; 
    }
}

void RK87(BPF& p, vec (*f)(vec, double, vec), double t0, double t1, double dt, Traj T, int DIM){

    double tol = 1e-6;
    double dtmax = 1e-3;
    int n_reject;
    int reject;
    vec y; 
    double t; 
    double tfinal; 
    mat F; 
    vec sol2;
    vec sol1;
    double error_1; 
    double error_step; 
    double tau; 

    for(int i = 0; i < p.n; i++){
        y = p.x.col(i); 
        n_reject = 0;
        reject = 0; 
        F.zeros(DIM,13);
        t = t0;
        tfinal = t1;
        while (t < tfinal) {
            if ((t + dt) > tfinal) {
                dt = tfinal - t;
            }

            F.col(0) = f(y, t, T.coef);
            for (int j = 0; j <= 11; j++) {
                vec yh = y + dt * F * a.col(j);
                F.col(j + 1) = f(yh, t + c1(j) * dt, T.coef);
            }

            sol2 = y + dt * F * b8;
            sol1 = y + dt * F * b7;

            error_1 = norm(sol1 - sol2);
            error_step = abs(error_1);
            tau = tol * std::max(norm(y, "inf"), 1.0);

            if (error_step <= tau) {
                t = t + dt;
                y = sol2;
                reject--;
            } else {
                n_reject++;
                reject = 1;
            }

            if (error_step == 0.0) {
                error_step = EPS * 10.0;
            }
            dt = std::min(dtmax, 0.9 * dt * pow((tau / error_step), POW));
            if (abs(dt) <= EPS) {
                if (reject == 0) {
                    std::cout << "WARNING!!! ode87. Step is very small!!!" << std::endl;
                    dt = EPS * 100;
                } else {
                    std::cout << "Error in ode87. Step is too small." << std::endl;
                }
            }
        }
        p.x.col(i) = y; 
    }

    return;
}

double effectiveSampleSize(const vec& weights) {
    vec weights_squared = pow(weights, 2);
    double ess = 1.0 / sum(weights_squared);
    return ess;
}

double computeLikelihood(const vec& x, vec (*h)(vec, vec), Meas M, Traj T) {
    vec diff = h(x, T.coef) - M.mean;
    return  exp(-0.5 * as_scalar(trans(diff) * inv(M.cov) * diff));
}

void measurement_update(BPF& p, vec (*h)(vec, vec), Meas M, Traj T, double N_THR){

    // Correction Step
    vec normalized_weights(p.n);
    vec likelihoods(p.n);
    
    for (int i = 0; i < p.n; ++i) {
        likelihoods(i) = computeLikelihood(p.x.col(i), h, M, T);
    }

    // Normalize weights
    double sum_likelihoods = sum(likelihoods);
    normalized_weights = likelihoods / sum_likelihoods;
    p.P = normalized_weights;

    // Compute ESS
    double ess = effectiveSampleSize(p.P);
    
    if (ess < N_THR) {
        // Resampling Step
        vec cdf(p.n);
        cdf(0) = p.P(0);
        for (int i = 1; i < p.n; ++i) {
            cdf(i) = cdf(i - 1) + p.P(i);
        }

        vec resampled_indices(p.n);
        std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);
        std::default_random_engine generator;
        for (int i = 0; i < p.n; ++i) {
            double rand_val = uniform_dist(generator);
            auto it = std::lower_bound(cdf.begin(), cdf.end(), rand_val);
            int index = std::distance(cdf.begin(), it);
            resampled_indices(i) = index;
        }

        // Resample particles
        mat resampled_particles(p.dim, p.n);
        vec resampled_particles_p(p.n);
        for (int i = 0; i < p.n; ++i) {
            resampled_particles.col(i) = p.x.col(resampled_indices(i));
            resampled_particles_p(i) = p.P(resampled_indices(i));
        }
        p.x = resampled_particles;
        p.P.fill(1.0 / p.n);  // Equal weights after resampling
    }

    return;
}

void run_bpf(vec (*f)(vec, double, vec), vec (*h)(vec, vec), int N, double N_THR, int DIM_f, int DIM_h, mat Q, Meas M, Traj T, std::string P_DIR, std::string M_DIR, int NUM_DIST, int NUM_MEAS, bool OUTPUT, bool RECORD, bool MEASURE){
    arma::arma_rng::set_seed_random();
    std::string P_PATH; 

    printf("Initializing BPF...\n\n");

    BPF p = BPF(N, M); 

    printf("Entering time marching...\n\n");

    clock_t start = clock(); 
    double tt = 0;
    for(int nm = 0; nm < NUM_MEAS; nm++){
        printf("Timestep: %d-0, Program time: %f s, Sim. time: %f TU, Size: %d\n", nm, ((double)(clock()-start))/CLOCKS_PER_SEC, tt, p.n); 
        if(RECORD){P_PATH = P_DIR + "/P" + std::to_string(nm) + "/pdf_0.txt"; record_data(P_PATH, p, tt);};

        vec tspan = linspace<vec>(0, M.T, NUM_DIST);

        int record_count = 1; int step_count = 1; double dt = 1e-3; 
        for(int i=1; i < tspan.n_elem; i++){
            
            RK87(p, f, tspan(i-1), tspan(i), dt, T, DIM_f);

            if (OUTPUT) { // print size to terminal
                printf("Timestep: %d-%d, Program time: %f s, Sim. time: %f TU, Size: %d\n", nm, step_count, ((double)(clock()-start))/CLOCKS_PER_SEC, tt + tspan(i), p.n); 
            }

            if (RECORD) { // record PDF
                P_PATH = P_DIR + "/P" + std::to_string(nm) + "/pdf_" + std::to_string(record_count) + ".txt"; 
                record_data(P_PATH, p, tt + tspan(i));
                record_count += 1;
            }

            step_count += 1; 
        }
        tt += M.T;

        if (!OUTPUT){ // print size to terminal
            printf("Timestep: %d-%d, Program time: %f s, Sim. time: %f TU, Size: %d\n", nm, step_count - 1, ((double)(clock()-start))/CLOCKS_PER_SEC, tt, p.n); 
        }

        if ((MEASURE) && (nm < NUM_MEAS - 1)) { // peform discrete measurement update

            printf("\nPERFORMING BAYESIAN UPDATE AT: %f TU...\n\n", tt);                      

            // reading new measurement
            std::string M_FILE = "measurement"+ std::to_string(nm + 1) + ".txt";
            M = Meas(DIM_h, M_DIR, M_FILE);                                       

            // performing discrete update
            measurement_update(p, h, M, T, N_THR);                                
        }
    }

    printf("Time marching complete.\n");

    return;
}