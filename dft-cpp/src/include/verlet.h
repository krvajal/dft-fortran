//
// Created by Miguel Angel Carvajal on 8/30/16.
//

#ifndef DFT_VERLET_H
#define DFT_VERLET_H

#include <functional>

template<typename Real,typename RealN>
struct VerletSolver                                         //Para resolver la ecuacion x''(t)=A(t,x)
{
    std::function<RealN(Real,RealN)> A;
    RealN x, p, a;
    Real dt;
    long it=0;
    VerletSolver(){}
    VerletSolver(std::function<RealN(Real,RealN)> Acel, RealN x0, RealN p0, Real dt0)
    {
        A=Acel;
        x=x0;
        p=p0;
        a=Acel(0.0,x0);
        dt=dt0;
    }

    void Iterate()
    {
        p+=0.5*dt*a;
        x+=dt*p;
        a=A(Time(),x);
        p+=0.5*dt*a;
        it++;
    }
    inline Real Time() const {return it*dt;}
};

class verlet_stepper{

    double pos,vel, acc;

public:
    verlet_stepper()= default;
    typedef std::function<void(double,double,double&, double)> system_type;

    void initialize(system_type sys, double pos, double vel, double t){
        this->pos = pos;
        this->vel  = vel;
        sys(pos,vel, this->acc, t);
    }
    void do_step(system_type sys, std::pair<double,double> &xinout, double t, double h){
            sys(xinout.first, xinout.second, this->acc, t);
            auto aux = xinout.first;
            xinout.first =  -this->pos + 2 * xinout.first + this->acc * h *h;
            pos =  aux;
            //
    }

};
#endif //DFT_VERLET_H
