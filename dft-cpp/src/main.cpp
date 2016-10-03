#include <armadillo>
#include <array>
#include <boost/numeric/odeint.hpp>
#include <fstream>
#include <functional>
#include <iostream>
#include <vector>
#include "numerov.h"
#include "quad/quad_simpson.h"
#include "verlet.h"

#include <boost/algorithm/string.hpp>


using namespace std;
using namespace boost::numeric::odeint;


double rmax = 10;
double h;


void plot(const arma::vec &col, const arma::vec &r);

typedef array<double, 3> state_type;


class radial_funct {

    arma::vec _V;
    int _z;
public:
    radial_funct(arma::vec V, int z) : _V(V), _z(z) {}

    double operator()(double r, double E) const {

        int index = r / h + 0.5;
        if (index >= _V.n_rows) {

            cout << "HEREEEE!!\n";
        }
        assert(index < _V.n_rows);
        double val = -2 * (E + _z / r - _V[index]);
        return val;
    }

};


void write_lorentz(const state_type &x, double t) {

    cout << t << '\t' << x[0] << '\t' << x[1] << '\t' << x[2] << endl;
}


void harm_oscilator(const double &x, const double &v, double &a, double t) {

    a = -2 * x;

}


void
poisson_radial(const double &x, const double &v, double &a, double r, double tstart, double h, const arma::vec &n)
{

    //compute index
    int index = r / h + 0.5;
    a = -4 * M_PI * r * n(index);
//    cout << "a = " << t << '\t' << a << endl;
}

void testodeint() {
    FILE *pipe = popen("/usr/local/bin/gnuplot", "w");
    assert(pipe);
    fprintf(pipe, "plot '-' w l \n");
    state_type x = {10.0, 1.0, 1.0};

    velocity_verlet<double> stepper;
    double x0 = 0;
    double v0 = 1;
    double t = 0;

    stepper.initialize(harm_oscilator, x0, v0, t);

    //or
    std::cout << stepper.order() << endl;
//    integrate_n_steps(stepper, harm_oscilator, std::make_pair(std::ref(x0), std::ref(v0)), t, 0.1, 1000, std::bind(write_oscilator,std::placeholders::_1,std::placeholders::_2, pipe));

    fprintf(pipe, "e\n");
}

void solve_radial_v1(function<double(double)> radial, arma::vec &u) {

    int nsteps = u.n_rows;

//    numerov_stepper stepper;
//    velocity_verlet<double> stepper2;

    auto ex = [](double r) { return r * exp(-r); };
//    arma::vec w = arma::zeros(u.n_rows);
    arma::vec w(u.n_rows);
    arma::vec r = arma::linspace(0, rmax, u.n_rows);
    arma::vec f = r.transform(radial);

    f(0) = 0;
//    cout << endl << endl;
//    cout << f << endl;


    u[nsteps - 1] = ex(rmax);
    u[nsteps - 2] = ex(rmax - h);
    w(nsteps - 1) = (1 - h * h * f(nsteps - 1) / 12.) * u(nsteps - 1);
    w(nsteps - 2) = (1 - h * h * f(nsteps - 2) / 12.) * u(nsteps - 2);

    for (int i = 3; i < nsteps; i++) {
        double r = rmax - (i - 1) * h;
        w[nsteps - i] = 2. * w[nsteps + 1 - i] - w[nsteps + 2 - i] + h * h * f(nsteps + 1 - i) * u[nsteps + 1 - i];
        u[nsteps - i] = w[nsteps - i] / (1. - h * h * (f(nsteps - i)) / 12.);

    }



    u(0) = 2 * u(1) - u(2) + f(1) * u(1) * h * h;


}


void solve_radial_v3(function<double(double)> radial, arma::vec &u) {

    int nsteps = u.n_rows;

    numerov_stepper stepper;


    auto ex = [](double r) { return r * exp(-r); };


    arma::vec r = arma::linspace(0, rmax, u.n_rows);
    arma::vec f = r.transform(radial);

    f(0) = 0;

    u[nsteps - 1] = ex(rmax);
    u[nsteps - 2] = ex(rmax - h);
    stepper.init({u(nsteps -1), u(nsteps -2)}, { radial(rmax),radial(rmax - h)},-h);
    for (int i = 3; i < nsteps; i++) {
        double r = rmax - (i - 2) * h;
        stepper.do_step(f, u(nsteps + 1 - i), 0, u(nsteps -i ),r);

    }



    u(0) = 2 * u(1) - u(2) + f(1) * u(1) * h * h;


}



void solve_radial_v2(function<double(double)> radial, arma::vec &u) {


    int nsteps = u.n_rows;

    velocity_verlet<double> stepper2;

    auto ex = [](double r) { return r * exp(-r); };


    double ucurrent = ex(rmax);
    double dudrcurrent = exp(-rmax) - rmax * exp(-rmax);
    stepper2.initialize(radial(rmax) * u(nsteps - 1));

    integrate_n_steps(stepper2, [&radial](const double &u, const double &dudt, double &acc, double r) {
                          acc = radial(r) * u;
                      }, std::make_pair(std::ref(ucurrent), std::ref(dudrcurrent)), rmax, -h, nsteps,
                      [&u](const std::pair<double, double> &state, double r) {
                          int index = r / h + 0.5;
                          u[index] = state.first;
                          if (index == 0) {
                              cout << "HEREEE!\n";
                          }
                      });

}


///
/// \param [in] radial the radial function
/// \param [out] u the value of u is placed here
void solve_radial(function<double(double)> radial, arma::vec &u) {

    solve_radial_v1(radial, u);
}

void plot(const arma::vec &x, const arma::vec &y) {

//    return;
    auto pip = popen("/usr/local/bin/gnuplot -p", "w");
    assert(pip);
//
    fprintf(pip, "set terminal x11\n");
    fprintf(pip, "plot  0, '-' with l\n");
    for (int j = 0; j < x.n_rows; ++j) {
        fprintf(pip, "%lf\t%lf\n", x(j), y(j));

    }
    fprintf(pip, "e\n");
    fflush(pip);
    cout << "Press a key to continue \n";
    getchar();


}

///
/// \param func a vector of the func sampled in equal spaced steps
/// \return a vector with the pair of values where function contains zeros
vector<pair<int, int>> bracket_zeros(const arma::vec &func) {


    typedef pair<int, int> range;
    vector<range> intervals;
    double f = func(0);
    range r = {0, 0};


    for (int i = 1; i < func.n_rows; i++) {
//        fprintf(pip,"%lf\n",func[i]);
        if ((func(i) * f) < 0) {
            intervals.push_back({i - 1, i});
            f = func(i);
        }
    }

//    fprintf(pip,"e\n");
//    fflush(pip);
//    getchar();
    assert(intervals.size());

    return intervals;
};


/**
 * @brief Solve the poisson equation using verlet and place the value in V
 * @param r  the radial
 * @param n the charge density
 * @param V the value computed for the potential
 */


void solve_poisson(const arma::vec &r, const arma::vec &n, arma::vec &V) {


    velocity_verlet<double> stepper;


    double dv0dt = 0.5;
    assert(V.n_rows == n.n_rows);


    int nsteps = n.n_rows - 1;


    auto f = [&n](const double &x, const double &v, double &a, double t) {
        poisson_radial(x, v, a, t, 0.0, h, n);

    };

    VerletSolver<double, double> stepper2([f](double x, double t) {
        double a;
        double v = 0;
        f(x, v, a, t);
        return a;

    }, 0, 1, h);

    verlet_stepper stepper3;


    V(0) = 0;
    stepper.initialize(f, V[0], dv0dt, 0.);
    stepper3.initialize(f, 0, 0, 0);

    auto observer = [&V](const pair<double, double> &v, double r) {
        int index = r / h + 0.5;

        V[index] = v.first;
//        cout << r << '\t' << v.first << endl;
    };

    double rr = 0;

//    integrate_n_steps(stepper, f, std::make_pair(std::ref(V(0)), std::ref(dv0dt)), rr, h, nsteps ,observer);
    double vv = h;
    auto xinout = std::make_pair(vv, dv0dt);
    for (int i = 1; i <= nsteps; i++) {
        V(i) = xinout.first;
        stepper3.do_step(f, xinout, i * h, h);
//        stepper2.Iterate();

    }
    V(0) = 0;
//    V(0) = 0;
    //compute alpha such that V(r) =  V(r) + alpha *r
    double qmax = 1;
    double alpha = (qmax - V(nsteps)) / rmax;
    cout << "alpha " << alpha << endl;

//    cout << "Vmax = " << V(V.n_rows - 1) << endl;
//    for (int j = 0; j < V.n_rows; ++j) {
//        cout << r[j] << '\t' << V[j] << endl;
//    }
//    terminate();
//    cout << "V= " << V<< endl;

    V += alpha * r;

}

/// Found a zero of a function inside the given interval and returns it when a given tolerance is
/// achieved
/// \param f the function to
/// \param start the start of the interval to find the root
/// \param end  the end of the interval to find the u
/// \param tol the tolerace used as stopping criteria
/// \param max_iter the maximun number of iterations to be performed
/// \return the zero found
/// \throw  invalid_argument when the ny
double bisec(function<double(double)> f, double start, double end, double tol = 1e-8, int max_iter = 100) {


    assert(start < end);
    double fstart = f(start);
    double fend = f(end);
    double mid = (start + end) * 0.5;
    double fmid = f(mid);
    if ((fstart * fend) > 0) {
        throw invalid_argument("bisec: invalid interval given for root findind");

    }

    for (int i = 0; i < max_iter; i++) {
        if ((end - start) < tol) {
            cout << "f(mid) = " << f(mid) << endl;
            return mid;
        }
        if ((fmid * fstart) > 0) {
            fstart = fmid;
            start = mid;
            mid = (mid + end) * 0.5;
            fmid = f(mid);
        } else {
            fend = fmid;
            end = mid;
            mid = (start + mid) * .5;

            fmid = f(mid);
        }

    }
    throw runtime_error("bisec: failed to converge after " + to_string(max_iter) + " iterations");

}


//retunrs  u(r), E0
//
std::pair<arma::vec, double> compute_ground_state(const arma::vec &r, const arma::vec &V) {

    int nsteps = V.n_rows - 1;


    double Emin = -3;
    double Emax = -0.000001;


    arma::vec E = arma::linspace<arma::vec>(Emin, Emax, nsteps + 1);
    arma::vec u0(nsteps + 1);
    arma::vec u(nsteps + 1);

    assert(E.n_rows > 1);


    auto radial = radial_funct(V, 2);

//    double ee = -0.5;
//    solve_radial([ee,radial](double r){return radial(r,ee);}, u);
//    plot(r,u);



    for (int i = 0; i < E.n_rows; i++) {
//        cout << "computing u(0) for E = " << E(i) << endl;
        solve_radial([&E, i, radial](double r) { return radial(r, E[i]); }, u);
        u0[i] = u[0];
    }





//    std::cout << "Finding zeros in the interval from Emin = " << Emin << " to Emax = " << Emax << endl;

    auto interval = bracket_zeros(u0).front();

    auto start = interval.first;
    auto end = interval.second;
    int n = end - start + 1;
    auto Estart = E(start);
    auto Eend = E(end);

    auto f = [nsteps, &radial](double EE) {
        arma::vec uu(nsteps + 1);

        solve_radial([EE, radial](double r) { return radial(r, EE); }, uu);
//        cout << u << endl;
        return uu[0];

    };

//    f(-2);
//    terminate();



    auto Ebind = bisec(f, Estart, Eend, 1e-8);


    solve_radial([&Ebind, radial](double r) { return radial(r, Ebind); }, u);

//    plot(arma::linspace(0,rmax,u.n_rows),u);
    cout << "u0 " << u(0) << endl;
    return {u, Ebind};


}



void compute_exchange_potential( arma::vec n, arma::vec & Vx){
    double C = - pow((6./M_PI),1./3.);
   n.transform([C](double x){  return C*pow(x,1./3.); });
   Vx = n;

}

void compute_correlation_potential( arma::vec n, arma::vec & Vc, arma::vec & Ec){




    n.transform([](double x){ return 3./(4. * M_PI * pow(2 * x, 1./3.));});
    auto & rs = n;


    //parameters for calculation of correlation in LDA (Ceperley-Alder y Perdew-Zunger)

    double A = 0.0311;
    double B = -0.048;
    double C =  0.002 ;
    double D =  - 0.0116;
    double gamma = -0.1423;
    double beta1  =  1.0529;
    double beta2 = 0.3334;
    

    auto ecf= [=](double rs){
        if(rs < 1){
            return A*log(rs) + B + C * rs * log(rs) + D * rs;

        }else{

            return gamma/(1 + beta1* sqrt(rs) + beta2* rs);

        }

    };

    Ec = n;
    Ec.transform(ecf);
    Ec(0)= 0;


    rs.transform([=](double rsx){

        if(rsx >= 1.0){
            double num  = 1 + 7./6 * beta1 * sqrt(rsx) + 4. /3. * beta2 * rsx;
            double den = 1 + beta1 * sqrt(rsx) + beta2 * rsx;
            return ecf(rsx)*num/den;
        } else{
            return A *log(rsx) + B - A/3. +2./3 * C * rsx * log(rsx) + ( 2 * D - C) * rsx / 3.;

        }

    });
    Vc = rs;
    Vc(0)  = 0;
}


void partA(){


    ofstream fout("U0.E.dat");
    int n = 1001;
    arma::vec u = arma::zeros(n,1);
    auto E = arma::linspace(-0.6, -0.01,200);

    for(int i=0; i < E.size(); i++){

        double Ei = E[i];
        cout << "Solving for E = " << Ei << endl;
        solve_radial([Ei](double r){


                return -2 *(Ei + 1/r);

        },u);
        fout << E[i] << '\t' << u[0] << endl;
    }
    fout.close();
}

int main() {

    std::cout << "Hello world\n";
    

    partA();



    // int n = 1001;

    // arma::vec V = arma::zeros(n);
    // arma::vec Vx = arma::zeros(n);
    // arma::vec Vc = arma::zeros(n);
    // arma::vec Ec = arma::zeros(n);

    // auto r = arma::linspace(0, rmax, n);
    // h = r[1] - r[0];


    // double Etotal = 0., EtotalPrev = 0;
    // double tol = 1e-3;
    // int max_iter = 20;
    // int iter = 0;
    // arma::vec u;
    // double EE;


    // while (true) {


    //     std::tie(u,EE) = compute_ground_state(r, V);

    //     //compute density

    //     double t = simpson(arma::vec(arma::square(u))) * h;

    //     u = u / sqrt(t * 4 * M_PI);
    //     arma::vec rho = arma::square(u) / arma::square(r);

    //     rho(0) = 0;
    //     solve_poisson(r, rho, V);
    //     compute_exchange_potential(rho,Vx);
    //     compute_correlation_potential(rho, Vc, Ec);

    //     // compute from V to Vh
    //     V = V / r;

    //     V(0) = 0;


    //     arma::vec u2 = arma::square(u);

    //     double Eh = 4 * M_PI * simpson(arma::vec(V % u2)) * h;
    //     double Ex =  4* M_PI * simpson(arma::vec(Vx % u2)) * h;

    //     arma::vec tmp =(Ec - Vc) % u2;
    //     cout << tmp << endl;
    //     double EEc = 4* M_PI * simpson(tmp) * h;


    //     Etotal = 2 * EE - 2*Eh  - Ex/2. ;// - Ec;

    //     V =  2*V + Vx + Vc;


    //     cout << "E = " << EE << endl;
    //     cout << "Ehartee = " << 2*Eh << endl;
    //     cout << "Ex = " << Ex << endl;
    //     cout << "Ec = " << EEc << endl;

    //     cout << "Etotal = " << Etotal << endl;
    //     cout << endl;
    //     iter ++;
    //     if (abs(Etotal - EtotalPrev) < tol) {
    //         cout << "Convergence achieved\n";
    //         cout << "Etotal =  " << Etotal << endl;
    //         cout << "Number of iterations: " << iter << endl;
    //         break;
    //     }
    //     EtotalPrev = Etotal;
    //     if(iter > max_iter){
    //         throw runtime_error("max number of iterations exceeded");
    //     }

    // }


}