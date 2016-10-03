#ifndef NUMEROV_H
#define NUMEROV_H

#include <functional>



class numerov_stepper{

    double w[2];
    double acc; // equal to fx
    double h;

public:

    typedef std::function<double(double)> system; //t => f(t)

    numerov_stepper(){

    }
    void do_step(system f, double pos, double vel, double &posout,double t){

        double tmp = -w[0] + 2 * w[1] + acc * pos *h*h;
        this->acc = f(t + h);
        w[0] = w[1];
        w[1] = tmp;
        posout = w_to_x(w[1],this->acc, h);
        nsteps ++;
    }

    void init(std::pair<double,double> x, std::pair<double,double> fx, double h){

        w[0] = x_to_w(x.first, fx.first, h);
        w[1] = x_to_w(x.second, fx.second, h);
        acc = fx.second;
        this->h = h;
        nsteps = 0;

    }

    //public property
    int nsteps = 0;

   static double x_to_w(double x, double fx, double h){
           return (1 - h*h*fx/12)*x;

    }
    static double w_to_x(double ww, double fx, double h){
        double  xout = ww/ (1.0 - h*h*fx/12);
        return xout;
    }

};

#endif //NUMEROV_H